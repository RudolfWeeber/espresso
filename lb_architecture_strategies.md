# Color-Gradient LB — Architecture Strategies

Audit and refactor proposals for the two-component (color-gradient, CG) lattice-Boltzmann implementation in `src/walberla_bridge`, alongside the existing single-component (SC) LB.

---

## 1. Current state (summary of audit)

The CG and SC models share a single class, `LBWalberlaImpl<FloatType, Architecture>`. Mode is selected at construction time via a runtime flag `m_two_components` (initialised from `viscosity.size() == 2`). The class carries the **union** of fields, kernels, and communicators for both models and branches on the flag at ~68 sites in the non-generated bridge code.

### Key files

| Area | File | Notes |
|---|---|---|
| Interface | `include/walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp` | 333 LoC, 62 virtuals; CG-only methods mixed in |
| Implementation | `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImpl.hpp` | 1249 LoC, 98 CG-tagged |
| Kernels | `src/walberla_bridge/src/lattice_boltzmann/generated_kernels/` | parallel `ColorGradient*` family alongside `StreamCollide*` family |
| Core wrapper | `src/core/lb/{Solver,LBWalberla,LBNone}.{hpp,cpp}` | `has_two_components()` flows through Solver |
| Particle coupling | `src/core/lb/particle_coupling.cpp` | solvation force gated on `has_two_components()` |
| Script interface | `src/script_interface/walberla/LBFluid.{hpp,cpp}` | decides SC vs CG by `viscosity.size()` |
| Python | `src/python/espressomd/lb.py` | `kinematic_viscosity` tuple, `sigma`, `beta`, `init_two_component()` |

### CG-only public surface today
- `set_collision_model_two_component(sigma, beta)`
- `init_two_component()`
- `has_two_components()`
- `get_color_gradients_at_pos(positions)`
- `add_solvation_forces_at_pos(positions, delta_mus)`

### Gaps vs single-component LB

| Feature | SC | CG |
|---|---|---|
| Lees–Edwards kernel | yes | **missing** |
| Thermalized collision | yes | declared but no thermalized CG kernel |
| `set_node_velocity` | works | **throws at runtime** |
| VTK writers | density/velocity/pressure | only component A; no φ, ρ_a/ρ_b, ∇φ |
| Pressure tensor | from PDF | from PDF[0] only |
| Checkpoint format | PDF + boundaries | unaudited; both PDFs + φ + ρ_{a,b}? |
| Particle drag | barycentric | barycentric (same) |
| Solvation force (∇φ) | n/a | yes, CG-only |
| EK coupling | yes | unclear |
| Automated tests | moderate | essentially absent |

The core complaint: the two models are co-housed under one class with ~68 branches, the public interface has CG-only methods that throw on SC and vice versa, and several signatures (`set_node_velocity`, `get_node_density`) violate LSP. CG has effectively no test coverage.

---

## 2. Strategy A — Shared slim base + sibling implementations + CG extension interface

**Goal:** replace the single `LBWalberlaImpl<F,A>` with two sibling implementations sharing infrastructure via a CRTP helper, behind a slimmed common interface plus a CG-specific extension interface.

```
LatticeModel
  └── LBWalberlaBase                          (common fluid interface, slimmed)
        └── LBWalberlaCGBase                  (adds CG-only virtuals)

implementations (templates):
  detail::LBWalberlaCommon<Derived, F, A>     (CRTP — shared infrastructure)
    ├── LBWalberlaImplSC<F, A>  : LBWalberlaBase
    └── LBWalberlaImplCG<F, A>  : LBWalberlaCGBase
```

`LBWalberlaCommon` is a compile-time mixin, not a polymorphic base — no diamond, no vtable cost.

### Plan

**Phase 0 — preparation**
1. Smoke test for CG locked into CI before any refactor.
2. Tag every `has_two_components()` site by category: integration loop, field allocation, ghost comm, accessor with different return shape, error branch.

**Phase 1 — extract `LBWalberlaCommon` (mechanical)**
3. Create `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaCommon.hpp` as a CRTP template carrying: block storage, lattice geometry, halo, boundary flag + UBB fields, `clear_boundaries`, `reallocate_ubb_field`, `update_boundary_from_shape`, `get_boundary_force*`, B-spline interpolation helpers, slice iteration helpers, `KernelTrait` consumption, `is_double_precision`, `is_gpu`, `stencil_size`, lattice-position checker.
4. Have the existing `LBWalberlaImpl` derive from it. No semantic change. Run tests.

**Phase 2 — split the leaf class**
5. Create `LBWalberlaImplSC<F,A>` by copying `LBWalberlaImpl` and deleting the CG arm of every branch. Keep thermalized + Lees-Edwards collision paths.
6. Create `LBWalberlaImplCG<F,A>` by copying `LBWalberlaImpl` and deleting the SC arm. Keep CG fields, kernels, comms, and the `stream → color_gradient → collide → reset_force` integration loop.
7. Delete the original `LBWalberlaImpl`.

**Phase 3 — split the public interface**
8. In `LBWalberlaBase.hpp`, delete `has_two_components`, `set_collision_model_two_component`, `init_two_component`, `get_color_gradients_at_pos`, `add_solvation_forces_at_pos`. Narrow return types of `get_node_density`/`set_node_density`/`get_viscosity`/`set_viscosity` from `vector<double>` to scalars (see §2.3 below).
9. Add `LBWalberlaCGBase.hpp` with the CG-only virtuals and the component-aware accessors.
10. `LBWalberlaImplSC` overrides `LBWalberlaBase`. `LBWalberlaImplCG` overrides `LBWalberlaCGBase`.

**Phase 4 — fix the dispatch sites**
11. `script_interface/walberla/LBFluid.cpp::make_instance()`: branch on `viscosity.size()` and instantiate the concrete type. Store as `shared_ptr<LBWalberlaBase>` for SC, `shared_ptr<LBWalberlaCGBase>` for CG.
12. `core/lb/Solver`: cache a `LBWalberlaCGBase*` (non-owning, null for SC). Replace `has_two_components()` queries with `solver.cg() != nullptr`. Optional follow-up: promote to `std::variant<LBWalberla, LBWalberlaCG, LBNone>`.
13. `core/lb/particle_coupling.cpp`: replace `has_two_components()` with `solver.cg()`; solvation path becomes unconditional within that branch.
14. `LBNone`: drop `get_color_gradients_at_pos`/`has_two_components`.

**Phase 5 — clean up newly-typeable asymmetries**
15. In `LBWalberlaImplCG`, either implement `set_node_velocity` correctly or move it to a follow-up `LBWalberlaSCBase` interface so it disappears at compile time.
16. Split VTK writers (`LBVTKSC.impl.hpp`, `LBVTKCG.impl.hpp`). CG writes φ, ρ_a, ρ_b, ∇φ.
17. Split checkpoint format with explicit version tag. SC stays byte-identical; CG bumps version.
18. Split node-access translation units.

**Phase 6 — tests, docs**
19. Add `LBWalberlaImplCG_unit_tests.cpp` mirroring the SC tests.
20. Add Python testsuite cases for CG.
21. Update `lb.py` docstrings.

### Diff size estimate

- Phase 1: ~600 LoC moved, ~0 deleted.
- Phase 2: ~1200 LoC duplicated then deduped via Phase 1; ~400 LoC of branching deleted.
- Phase 3: ~10 LoC interface shuffle.
- Phase 4: ~50 LoC in Solver / particle_coupling.
- Phases 5–6: incremental, each landable as its own PR.

Only Phases 1 and 2 must land together.

### 2.1 `LBWalberlaBase` skeleton (slimmed)

```cpp
// include/walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp
#pragma once

#include <walberla_bridge/LatticeModel.hpp>
#include <walberla_bridge/lattice_boltzmann/LeesEdwardsPack.hpp>
#include <utils/Vector.hpp>

#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

class LBWalberlaBase : public LatticeModel {
public:
  ~LBWalberlaBase() override = default;

  // --- integration & ghost layers ---
  virtual void integrate() = 0;
  virtual void ghost_communication() = 0;
  virtual void ghost_communication_pdf() = 0;
  virtual void ghost_communication_vel() = 0;
  virtual void ghost_communication_laf() = 0;

  // --- lattice / precision / device ---
  virtual std::size_t stencil_size() const noexcept = 0;
  [[nodiscard]] virtual bool is_double_precision() const noexcept = 0;
  [[nodiscard]] virtual bool is_gpu() const noexcept = 0;
  virtual std::function<bool(Utils::Vector3d const &)>
  make_lattice_position_checker(bool consider_points_in_halo) const = 0;

  // --- interpolated fields at off-lattice positions ---
  virtual std::optional<Utils::Vector3d>
  get_velocity_at_pos(Utils::Vector3d const &, bool consider_points_in_halo = false) const = 0;
  virtual std::vector<Utils::Vector3d>
  get_velocities_at_pos(std::vector<Utils::Vector3d> const &) = 0;
  virtual std::optional<double>
  get_density_at_pos(Utils::Vector3d const &, bool consider_points_in_halo = false) const = 0;
  virtual std::vector<double>
  get_densities_at_pos(std::vector<Utils::Vector3d> const &) = 0;

  // --- forces at off-lattice positions ---
  virtual bool add_force_at_pos(Utils::Vector3d const &, Utils::Vector3d const &) = 0;
  virtual void add_forces_at_pos(std::vector<Utils::Vector3d> const &positions,
                                 std::vector<Utils::Vector3d> const &forces) = 0;
  virtual void add_density_weighted_forces_at_pos(
      std::vector<Utils::Vector3d> const &positions,
      std::vector<Utils::Vector3d> const &forces) = 0;

  // --- nodal force accessors ---
  virtual std::optional<Utils::Vector3d>
  get_node_force_to_be_applied(Utils::Vector3i const &) const = 0;
  virtual std::optional<Utils::Vector3d>
  get_node_last_applied_force(Utils::Vector3i const &, bool consider_ghosts = false) const = 0;
  virtual bool set_node_last_applied_force(Utils::Vector3i const &,
                                           Utils::Vector3d const &) = 0;
  virtual std::vector<double>
  get_slice_last_applied_force(Utils::Vector3i const &, Utils::Vector3i const &) const = 0;
  virtual void set_slice_last_applied_force(Utils::Vector3i const &,
                                            Utils::Vector3i const &,
                                            std::vector<double> const &) = 0;

  // --- populations (variable length is intrinsic; vector is appropriate) ---
  virtual std::optional<std::vector<double>>
  get_node_population(Utils::Vector3i const &, bool consider_ghosts = false) const = 0;
  virtual bool set_node_population(Utils::Vector3i const &,
                                   std::vector<double> const &) = 0;
  virtual std::vector<double>
  get_slice_population(Utils::Vector3i const &, Utils::Vector3i const &) const = 0;
  virtual void set_slice_population(Utils::Vector3i const &, Utils::Vector3i const &,
                                    std::vector<double> const &) = 0;

  // --- velocity (barycentric in both models) ---
  virtual std::optional<Utils::Vector3d>
  get_node_velocity(Utils::Vector3i const &, bool consider_ghosts = false) const = 0;
  virtual bool set_node_velocity(Utils::Vector3i const &, Utils::Vector3d const &) = 0;
  virtual std::vector<double>
  get_slice_velocity(Utils::Vector3i const &, Utils::Vector3i const &) const = 0;
  virtual void set_slice_velocity(Utils::Vector3i const &, Utils::Vector3i const &,
                                  std::vector<double> const &) = 0;

  // --- density (scalar = barycentric / total; CG component split lives on CGBase) ---
  virtual std::optional<double>
  get_node_density(Utils::Vector3i const &, bool consider_ghosts = false) const = 0;
  virtual bool set_node_density(Utils::Vector3i const &, double) = 0;
  virtual std::vector<double>
  get_slice_density(Utils::Vector3i const &, Utils::Vector3i const &) const = 0;
  virtual void set_slice_density(Utils::Vector3i const &, Utils::Vector3i const &,
                                 std::vector<double> const &) = 0;

  // --- boundaries (UBB) ---
  virtual std::optional<Utils::Vector3d>
  get_node_velocity_at_boundary(Utils::Vector3i const &,
                                bool consider_ghosts = false) const = 0;
  virtual bool set_node_velocity_at_boundary(Utils::Vector3i const &,
                                             Utils::Vector3d const &) = 0;
  virtual std::vector<std::optional<Utils::Vector3d>>
  get_slice_velocity_at_boundary(Utils::Vector3i const &,
                                 Utils::Vector3i const &) const = 0;
  virtual void set_slice_velocity_at_boundary(
      Utils::Vector3i const &, Utils::Vector3i const &,
      std::vector<std::optional<Utils::Vector3d>> const &) = 0;
  virtual std::optional<Utils::Vector3d>
  get_node_boundary_force(Utils::Vector3i const &) const = 0;
  virtual bool remove_node_from_boundary(Utils::Vector3i const &) = 0;
  virtual std::optional<bool>
  get_node_is_boundary(Utils::Vector3i const &, bool consider_ghosts = false) const = 0;
  virtual std::vector<bool>
  get_slice_is_boundary(Utils::Vector3i const &, Utils::Vector3i const &) const = 0;
  virtual void reallocate_ubb_field() = 0;
  virtual void clear_boundaries() = 0;
  virtual void update_boundary_from_shape(std::vector<int> const &,
                                          std::vector<double> const &) = 0;

  // --- collision model (currently throws on CG; candidate for SC-only follow-up) ---
  virtual void set_collision_model(double kT, unsigned int seed) = 0;
  virtual void set_collision_model(std::unique_ptr<LeesEdwardsPack> &&) = 0;
  virtual void check_lebc(unsigned int shear_direction,
                          unsigned int shear_plane_normal) const = 0;

  // --- pressure tensor / momentum / external force ---
  virtual std::optional<Utils::VectorXd<9>>
  get_node_pressure_tensor(Utils::Vector3i const &) const = 0;
  virtual std::vector<double>
  get_slice_pressure_tensor(Utils::Vector3i const &, Utils::Vector3i const &) const = 0;
  virtual Utils::Vector3d
  get_boundary_force_from_shape(std::vector<int> const &) const = 0;
  virtual Utils::Vector3d get_boundary_force() const = 0;
  virtual Utils::VectorXd<9> get_pressure_tensor() const = 0;
  virtual Utils::Vector3d get_momentum() const = 0;
  virtual void set_external_force(Utils::Vector3d const &) = 0;
  virtual Utils::Vector3d get_external_force() const noexcept = 0;

  // --- bulk parameters (scalar viscosity; CG component-wise viscosity on CGBase) ---
  virtual void set_viscosity(double) = 0;
  [[nodiscard]] virtual double get_viscosity() const = 0;
  virtual double get_density() const noexcept = 0;
  virtual double get_kT() const noexcept = 0;
  virtual unsigned int get_seed() const noexcept = 0;
  [[nodiscard]] virtual std::optional<uint64_t> get_rng_state() const = 0;
  virtual void set_rng_state(uint64_t) = 0;

  // --- field IDs (e.g. for EK coupling) ---
  [[nodiscard]] virtual std::size_t get_velocity_field_id() const noexcept = 0;
  [[nodiscard]] virtual std::size_t get_force_field_id() const noexcept = 0;
};
```

### 2.2 `LBWalberlaCGBase` skeleton (CG extension)

```cpp
// include/walberla_bridge/lattice_boltzmann/LBWalberlaCGBase.hpp
#pragma once

#include "LBWalberlaBase.hpp"
#include <utils/Vector.hpp>

#include <array>
#include <vector>

/**
 * Public interface of the two-component (color-gradient) LB solver.
 * Inherits the generic fluid interface and adds model-specific
 * configuration, initialization, and per-component accessors.
 */
class LBWalberlaCGBase : public LBWalberlaBase {
public:
  ~LBWalberlaCGBase() override = default;

  // --- CG-specific configuration ---
  virtual void set_collision_model_cg(double sigma, double beta) = 0;
  virtual void init_pdfs_from_components() = 0;

  // --- CG-specific observables / forces ---
  virtual std::vector<Utils::Vector3d>
  get_color_gradients_at_pos(std::vector<Utils::Vector3d> const &positions) = 0;
  virtual void
  add_solvation_forces_at_pos(std::vector<Utils::Vector3d> const &positions,
                              std::vector<double> const &delta_mus) = 0;

  // --- per-component accessors (the inherited get_node_density returns the
  //     total ρ = ρ_a + ρ_b; this returns the split) ---
  virtual std::optional<std::array<double, 2>>
  get_node_component_densities(Utils::Vector3i const &,
                               bool consider_ghosts = false) const = 0;
  virtual bool
  set_node_component_densities(Utils::Vector3i const &,
                               std::array<double, 2> const &) = 0;
  virtual std::vector<double>
  get_slice_component_densities(Utils::Vector3i const &,
                                Utils::Vector3i const &) const = 0;

  // --- component-wise viscosity ---
  virtual void set_component_viscosities(std::array<double, 2> const &) = 0;
  [[nodiscard]] virtual std::array<double, 2>
  get_component_viscosities() const = 0;
};
```

### 2.3 Note on `vector<double>` returns

The original draft kept `optional<vector<double>>` for `get_node_density` etc. Two problems:

- **Allocation cost.** Each per-node call is a ~32-byte heap allocation. Negligible per-call but visible in a tight Python loop over many nodes. Slice APIs amortise this and stay as `vector<double>` — they are bulk and correct.
- **API surprise.** Single-component Python users see `[1.0]` instead of `1.0` from `lbf[node].density`.

Fix in the skeletons above:
- `get_node_density` / `set_node_density` / `get_viscosity` / `set_viscosity` on `LBWalberlaBase` use `double`.
- Per-component CG accessors live on `LBWalberlaCGBase` with `array<double, 2>` returns.
- `get_node_population` keeps `vector<double>` — variable Q is intrinsic.

This is a Python-visible behaviour change; deprecation cycle recommended.

### 2.4 Detecting CG from the core

Two equally clean options:

```cpp
// Option 1: dynamic_cast at construction, cache.
struct Solver {
  std::shared_ptr<LBWalberlaBase> base;
  LBWalberlaCGBase *cg = nullptr;  // non-owning, == base.get() iff CG
};

// Option 2 (follow-up): variant in the actor.
using HydrodynamicsActor =
    std::variant<std::monostate, LBWalberla, LBWalberlaCG, LBNone>;
```

`particle_coupling.cpp` becomes:

```cpp
if (auto *cg = m_lb.cg()) {
  auto grads = cg->get_color_gradients_at_pos(positions);
  cg->add_solvation_forces_at_pos(positions, delta_mus);
}
```

### 2.5 Duplication budget

| Site | Near-duplicate LoC | Notes |
|---|---:|---|
| `LBWalberlaBase.hpp` | 0 | single source for all common virtuals |
| `LBWalberlaCGBase.hpp` | ~10 | header preamble, class shell |
| `LBWalberlaImplSC` / `LBWalberlaImplCG` override decls | ~60 | ~30 common methods × 2 `override` lines |
| Trivial CRTP-delegating one-liners | ~20 | repeated in both leaves |
| `core/lb/LBWalberla` | 0 | unchanged; adds one `cg()` accessor |
| `core/lb/Solver` | 0 | unchanged |
| `script_interface/walberla/LBFluid` | 0 | unchanged; existing internal branches |
| `script_interface/walberla/LBFluidNode` / `LBFluidSlice` | 0 | unchanged |
| `python/espressomd/lb.py` | 0 | unchanged |
| **Total** | **~90 LoC** | |

The override-declaration cost is irreducible: each leaf must declare its overrides.

---

## 3. Strategy C — Independent top-level classes with shared infrastructure

**Goal:** treat SC and CG as **two unrelated top-level features**. No shared virtual base. Common waLBerla scaffolding lives in a non-polymorphic helper that both classes inherit from (CRTP). Core dispatches via `std::variant`.

```
LatticeModel                                  (kept for the Lattice abstraction;
  └─ (nothing else — no shared LB base)        not an LB interface)

src/walberla_bridge/src/common/
  WalberlaFluidInfrastructure<F, A>           (CRTP / mixin)
    - BlockForest / StructuredBlockStorage
    - boundary flag field, UBB field, mask helpers
    - B-spline interpolation
    - slice / cuboid iteration
    - PackInfo + RegularCommunicator scaffolding
    - precision / arch traits
    - VTK & checkpoint primitives (no field bindings)

src/walberla_bridge/src/lattice_boltzmann/
  LBWalberlaBase                              (SC's own interface)
    └── LBWalberlaImpl<F, A>                  (SC implementation)

src/walberla_bridge/src/multicomponent/
  LBWalberlaCGBase                            (independent interface)
    └── LBWalberlaCGImpl<F, A>                (CG implementation)

src/core/lb/
  Solver { std::variant<std::monostate, LBWalberla, LBWalberlaCG, LBNone> actor; }
```

The two `*Base` classes do **not** share a parent. They look similar where the physics overlaps but that is type-level coincidence, not a Liskov contract.

### Contents of `WalberlaFluidInfrastructure`

Strict criterion: nothing here may know that a "PDF" exists. It is the waLBerla glue layer.

- block storage, lattice geometry, halo width
- creation/destruction of `GhostLayerField<T,N>` by tag (callers register their own field IDs)
- generic `RegularFullCommunicator` factory + `PackInfo` for registered fields
- boundary flag field, UBB-field allocation, raster-to-mask, `clear_boundaries`, `update_boundary_from_shape`, `get_boundary_force*` (parameterised on the force-field id)
- B-spline interpolation kernels (templated on field type)
- slice iteration helpers
- precision/arch traits
- GPU block-data allocation, host/device sync helpers
- checkpoint I/O primitives (stream open/close, version tag) — **not** the schema
- VTK file-writer primitives — **not** the field bindings

What stays in the concrete LB classes: `integrate()` loop, ghost-comm sequence, collision/stream kernel ownership, `get_node_density`, `get_node_velocity`, pressure tensor extraction, particle-coupling hooks.

### Plan

**Phase 0 — preparation**
1. CG smoke tests locked into CI first.
2. Bucket every method in `LBWalberlaImpl.hpp` as (i) infrastructure, (ii) SC-specific, (iii) CG-specific.

**Phase 1 — extract `WalberlaFluidInfrastructure`**
3. Create `src/walberla_bridge/src/common/WalberlaFluidInfrastructure.hpp` as a CRTP mixin.
4. Move bucket-(i) into it. Existing `LBWalberlaImpl` inherits from it. No semantic change; run tests.
5. Move `LatticeModel` to live behind this helper.

**Phase 2 — fork the implementations**
6. Create `LBWalberlaCGImpl<F,A>` by copying `LBWalberlaImpl`, deleting the SC branches; inherit `WalberlaFluidInfrastructure` and `LBWalberlaCGBase`.
7. Trim original `LBWalberlaImpl` to SC only; delete `m_two_components`, CG fields, CG comms, CG kernels, CG integration sequence, `init_two_component`.
8. Move CG kernel files to `src/walberla_bridge/src/multicomponent/generated_kernels/` (optional cosmetic step).

**Phase 3 — split the public interfaces**
9. Slim `LBWalberlaBase.hpp` (SC): remove `has_two_components`, `set_collision_model_two_component`, `init_two_component`, `get_color_gradients_at_pos`, `add_solvation_forces_at_pos`.
10. New `include/walberla_bridge/multicomponent/LBWalberlaCGBase.hpp` — its own abstract class, no `LBWalberlaBase` inheritance. Re-declares overlapping common methods with the natural CG signatures (size 2 returns where applicable).

**Phase 4 — split the core boundary**
11. `core/lb/`: split into `LBWalberla.{hpp,cpp}` and new `LBWalberlaCG.{hpp,cpp}`. Both wrap their respective bridge type.
12. `core/lb/Solver`:
    ```cpp
    using HydrodynamicsActor =
        std::variant<std::monostate, LBWalberla, LBWalberlaCG, LBNone>;
    ```
    Each public Solver method either visits all arms (common methods) or visits only the CG arm and throws elsewhere (CG-only methods).

**Phase 5 — split the model-specific feature code**
13. `particle_coupling.cpp` — `std::get_if<LBWalberlaCG>(&solver.actor)`.
14. VTK: `LBVTK_SC.cpp`, `LBVTK_CG.cpp`. CG binds φ, ρ_a, ρ_b, ∇φ.
15. Checkpoint: two schemas, version-tagged.
16. Node access: `LBNodeAccessSC.impl.hpp`, `LBNodeAccessCG.impl.hpp`.

**Phase 6 — script interface & Python**
17. Split `LBFluid` script-interface class into `LBFluid` (SC) and `LBFluidCG` (CG). Move `sigma`, `beta`, `init_two_component` to CG.
18. Likely split `LBFluidNode`/`LBFluidSlice` too, or template them on the bridge type.
19. In Python, add `class LBFluidCG(HydrodynamicInteraction)`. `LBFluid` stops accepting tuple viscosity / `sigma` / `beta`. **User-visible API break** — deprecation cycle required.

**Phase 7 — tests, docs**
20. Independent test directories.
21. Separate user-guide chapters.

### `LBWalberlaCGBase` skeleton under Strategy C

No shared base, so the CG interface re-declares the overlapping common surface. Trimmed deliberately — only methods CG actually supports today.

```cpp
// include/walberla_bridge/multicomponent/LBWalberlaCGBase.hpp
#pragma once

#include <walberla_bridge/LatticeModel.hpp>
#include <utils/Vector.hpp>

#include <array>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

class LBWalberlaCGBase : public LatticeModel {
public:
  ~LBWalberlaCGBase() override = default;

  // --- integration & ghost layers ---
  virtual void integrate() = 0;
  virtual void ghost_communication() = 0;

  // --- lattice / device ---
  virtual std::size_t stencil_size() const noexcept = 0;
  [[nodiscard]] virtual bool is_double_precision() const noexcept = 0;
  [[nodiscard]] virtual bool is_gpu() const noexcept = 0;
  virtual std::function<bool(Utils::Vector3d const &)>
  make_lattice_position_checker(bool consider_points_in_halo) const = 0;

  // --- interpolated fields ---
  virtual std::optional<Utils::Vector3d>
  get_velocity_at_pos(Utils::Vector3d const &, bool = false) const = 0;
  virtual std::vector<Utils::Vector3d>
  get_velocities_at_pos(std::vector<Utils::Vector3d> const &) = 0;
  virtual std::vector<std::array<double, 2>>
  get_densities_at_pos(std::vector<Utils::Vector3d> const &) = 0;

  // --- CG-only: color gradient & solvation ---
  virtual std::vector<Utils::Vector3d>
  get_color_gradients_at_pos(std::vector<Utils::Vector3d> const &) = 0;
  virtual void
  add_solvation_forces_at_pos(std::vector<Utils::Vector3d> const &positions,
                              std::vector<double> const &delta_mus) = 0;

  // --- per-node accessors (no set_node_velocity on purpose) ---
  virtual std::optional<std::array<double, 2>>
  get_node_density(Utils::Vector3i const &, bool = false) const = 0;
  virtual bool set_node_density(Utils::Vector3i const &,
                                std::array<double, 2> const &) = 0;
  virtual std::optional<Utils::Vector3d>
  get_node_velocity(Utils::Vector3i const &, bool = false) const = 0;
  virtual std::optional<std::vector<double>>
  get_node_population(Utils::Vector3i const &, bool = false) const = 0;
  virtual bool set_node_population(Utils::Vector3i const &,
                                   std::vector<double> const &) = 0;

  // --- boundaries ---
  virtual bool set_node_velocity_at_boundary(Utils::Vector3i const &,
                                             Utils::Vector3d const &) = 0;
  virtual bool remove_node_from_boundary(Utils::Vector3i const &) = 0;
  virtual void reallocate_ubb_field() = 0;
  virtual void clear_boundaries() = 0;
  virtual void update_boundary_from_shape(std::vector<int> const &,
                                          std::vector<double> const &) = 0;
  virtual Utils::Vector3d get_boundary_force() const = 0;

  // --- forces at off-lattice positions ---
  virtual void add_forces_at_pos(std::vector<Utils::Vector3d> const &,
                                 std::vector<Utils::Vector3d> const &) = 0;
  virtual void add_density_weighted_forces_at_pos(
      std::vector<Utils::Vector3d> const &,
      std::vector<Utils::Vector3d> const &) = 0;

  // --- CG-specific configuration ---
  virtual void set_collision_model_cg(double sigma, double beta) = 0;
  virtual void init_pdfs_from_components() = 0;

  // --- bulk parameters ---
  virtual void set_viscosity(std::array<double, 2> const &) = 0;
  [[nodiscard]] virtual std::array<double, 2> get_viscosity() const = 0;

  // --- diagnostics ---
  virtual Utils::Vector3d get_momentum() const = 0;
  virtual Utils::VectorXd<9> get_pressure_tensor() const = 0;
};
```

Differences from `LBWalberlaBase` are encoded in the **type system**:

- `get_node_density` returns `optional<array<double,2>>`, not `optional<double>`.
- `set_viscosity` takes `array<double,2>`, not `double`.
- `set_node_velocity` is **deliberately absent** — operation ill-defined for current CG implementation; the compiler enforces this.
- Thermalization (`set_collision_model(kT,seed)`, `check_lebc`, `get_kT`, RNG state) is absent — gets added back explicitly when/if CG gains it.
- Boundary surface is trimmed to what is actually wired up today.

### Duplication budget

| Site | Near-duplicate LoC | Notes |
|---|---:|---|
| `LBWalberlaBase.hpp` (slimmed SC) | 0 | single source for SC |
| `LBWalberlaCGBase.hpp` overlap | ~140 | ~25 base virtuals duplicated with identical signatures + docs; ~10 with similar-but-not-identical signatures |
| Override decls in both Impl classes | ~60 | same as A |
| Trivial CRTP-delegating one-liners | ~20 | same as A |
| `core/lb/LBWalberla` vs new `LBWalberlaCG` | ~200 | new wrapper mirrors ~⅔ of the SC wrapper |
| `core/lb/Solver` visitor lambdas | ~80 | visit arms for ~30 common Solver methods |
| `script_interface/walberla/LBFluidCG` (new) | ~350 | mirror of `LBFluid` (623 LoC) |
| `script_interface/walberla/LBFluidNodeCG` (new) | ~180 | mirror of `LBFluidNode` (269 LoC) |
| `script_interface/walberla/LBFluidSliceCG` (new) | ~200 | mirror of `LBFluidSlice` (303 LoC) |
| `python/espressomd/lb.py` — new `LBFluidCG` | ~130 | mirror of `LBFluid` |
| **Total** | **~1,360 LoC** | |

Possible mitigations inside Strategy C (templated CRTP wrappers in script-interface node/slice, common SI base for `LBFluid`/`LBFluidCG`, Python mixin) bring the floor down to roughly **~900 LoC**.

---

## 4. Side-by-side

| Dimension | Strategy A | Strategy C |
|---|---|---|
| Shared virtual base | yes (`LBWalberlaBase`) | no — two unrelated trees |
| Polymorphic dispatch in core | through `LBWalberlaBase` + `dynamic_cast` for CG | `std::variant` + `std::visit` |
| Python API stability | preserved | break (`LBFluidCG` added; tuple-viscosity removed from `LBFluid`) |
| `set_node_velocity` on CG | still has to exist (throws or implements) | doesn't exist (compile-time enforced) |
| `set_collision_model(kT, seed)` on CG | exists on base; behaviour TBD | doesn't exist |
| Adding a third multicomponent model | another sibling under `LBWalberlaCGBase` or peer; Liskov negotiation | clean drop-in; pays off around N=3 models |
| Near-duplicate boilerplate | **~90 LoC** | **~1,360 LoC** (or ~900 with mitigations) |
| Diff size of the refactor itself | smaller (bridge-only) | larger (reaches Python) |
| Footprint of follow-ups (VTK, checkpoint, tests split) | same | same |
| Risk to existing SC behaviour | low | low (SC bridge stays put; only Solver dispatch changes) |
| Compile-time enforcement of model differences | partial (CG-only methods on extension; common base still has cross-cutting methods) | full |

### When to pick which

**Pick A if all of:**
1. CG is the only multicomponent model on the roadmap.
2. Python API stability matters more than C++ API ergonomics.
3. You want the minimum diff that removes the 68 internal branches.

**Pick C if any of:**
1. A third or fourth multicomponent model (Shan–Chen, ternary, Allen–Cahn) is on the ~12-month roadmap.
2. The `set_node_velocity`-throws-at-runtime kind of behaviour is unacceptable; compile-time enforcement is required.
3. You are willing to break (with deprecation) the Python `LBFluid` constructor to expose `LBFluidCG`.
4. The team prefers explicit `std::visit` dispatch over polymorphic virtual calls at the core/bridge boundary.

### Recommendation

Land **Strategy A** first. Phases 1 (extract `LBWalberlaCommon`) and 2 (split leaves) are identical between A and C — they are the cheap, mechanical part with the biggest payoff (deleting the ~68 internal branches). Re-evaluate before Phase 3.

If you later decide to upgrade to C, the additional work is bounded: delete the shared `LBWalberlaBase`, duplicate the signatures into `LBWalberlaCGBase`, switch core to `std::variant`, split Python. That is a tractable second PR. The reverse direction (C → A) is not as clean, so A is the safer first commitment.

---

## 5. Open questions / follow-ups (orthogonal to A vs C)

These are real gaps that should land regardless of which strategy is chosen:

- Thermalized collision kernel for CG (does it exist?).
- Lees–Edwards kernel for CG (almost certainly missing).
- VTK output for φ, ρ_a, ρ_b, ∇φ.
- Component-wise pressure tensor — is the current `get_pressure_tensor` meaningful for CG?
- Checkpoint format for CG — audit and version-tag.
- Python testsuite coverage for CG (currently absent).
- EK + CG interaction — is it functional?
