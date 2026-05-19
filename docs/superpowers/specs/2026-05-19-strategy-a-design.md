# Strategy A — Single-component / Color-gradient LB split (Phases 1–4)

**Date:** 2026-05-19
**Branch:** `strategy-a` (worktree at `/tikhome/weeber/es-strategy-a`), branched from `color_gradient`
**Source design:** `lb_architecture_strategies.md` in repo root, Strategy A
**Scope:** Phases 1–4 of Strategy A. Phase 5 (VTK / checkpoint / node-access split, fixing `set_node_velocity` on color-gradient) and Phase 6 (C++ unit tests for color-gradient, additional Python tests) are explicitly deferred to follow-up PRs.

---

## 1. Goal

Replace the single class `LBWalberlaImpl<FloatType, Architecture>`, which carries the union of single-component (SC) and two-component color-gradient (CG) fields and branches on a runtime `m_two_components` flag at ~68 sites, with:

- a CRTP infrastructure mixin `LBWalberlaCommon<Derived, F, A>` holding shared waLBerla scaffolding,
- two sibling leaf templates `LBWalberlaImplSingleComponent<F, A>` and `LBWalberlaImplColorGradient<F, A>`,
- a slimmed common interface `LBWalberlaBase`,
- a color-gradient extension interface `LBWalberlaColorGradientBase : public LBWalberlaBase`.

Core (`src/core/lb/Solver`, `src/core/lb/particle_coupling.cpp`, `LBNone`) and the script-interface (`LBFluid`) dispatch via the new typed accessor `solver.color_gradient()` instead of the runtime flag.

## 2. Out of scope

- Splitting VTK writers (φ, ρ_a, ρ_b, ∇φ output for CG).
- Versioned checkpoint format for CG.
- Implementing `set_node_velocity` correctly on color-gradient — it continues to throw.
- Adding a thermalized color-gradient kernel.
- Adding Lees–Edwards collision for color-gradient.
- C++ unit tests for color-gradient.
- A `std::variant`-based Solver upgrade (Strategy C territory).
- Renaming the generated `ColorGradient*` kernel files.

## 3. Naming convention

Full-word names throughout. No `SC` / `CG` abbreviations in new identifiers.

| Identifier | Value |
|---|---|
| Leaf class (single-component) | `LBWalberlaImplSingleComponent<FloatType, Architecture>` |
| Leaf class (color-gradient) | `LBWalberlaImplColorGradient<FloatType, Architecture>` |
| Common interface (slimmed) | `LBWalberlaBase` (existing, slimmed) |
| Color-gradient extension interface | `LBWalberlaColorGradientBase : public LBWalberlaBase` |
| Shared CRTP mixin | `LBWalberlaCommon<Derived, FloatType, Architecture>` |
| Solver accessor | `solver.color_gradient()` returning `LBWalberlaColorGradientBase *` (nullptr in SC) |
| Cached pointer member | `m_color_gradient` |
| CG configuration method | `set_collision_model_color_gradient(double sigma, double beta)` |
| Per-component density accessor | `get_node_component_densities` / `set_node_component_densities` |

File names introduced:
- `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaCommon.hpp` (new).
- `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImplSingleComponent.hpp` (new).
- `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImplColorGradient.hpp` (new).
- `include/walberla_bridge/lattice_boltzmann/LBWalberlaColorGradientBase.hpp` (new).
- `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImpl.hpp` is deleted in Phase 2c.

Generated kernel files keep their existing names.

## 4. `LBWalberlaCommon` — contents of the CRTP mixin

Inclusion criterion: used by both leaves today and does not require the leaf's PDF layout to compile. CRTP `static_cast<Derived&>(*this)` is used wherever the helper needs the leaf's typed field IDs (boundary handling and force/velocity accessors).

**Storage and lattice geometry**
- `m_lattice` (`shared_ptr<LatticeWalberla>`), `get_lattice()`, `get_blocks()`.
- Halo width, ghost-layer accessors.
- Traits: `is_double_precision()`, `is_gpu()`, `stencil_size()`.
- `make_lattice_position_checker(bool consider_points_in_halo)`.

**Boundary / UBB scaffolding**
- Boundary flag field id + UBB-velocity field id (registered by the leaf, owned via the leaf's `BlockDataID`s; the mixin holds typedefs and helpers).
- `reallocate_ubb_field()`, `clear_boundaries()`, `update_boundary_from_shape(raster, values)`.
- `get_boundary_force_from_shape(raster)`, `get_boundary_force()`.
- `get_node_velocity_at_boundary` / `set_node_velocity_at_boundary` and slice equivalents.
- `get_node_is_boundary`, `remove_node_from_boundary`, `get_node_boundary_force`.

**Off-lattice interpolation (read paths and weight math)**
- `BSplineWeights::compute(Vector3d pos)` returning the 8 lattice nodes and weights.
- `get_velocity_at_pos` / `get_velocities_at_pos` (one velocity field exists in both leaves; this is barycentric in CG).
- `get_density_at_pos` / `get_densities_at_pos` (returns total ρ; per-component positional accessors live on the CG base if/when added).

**Ghost-comm scaffolding (not the schemas)**
- `RegularFullCommunicator` factory.
- `PackInfo` registration helper.
- A virtual `ghost_communication()` hook that leaves override; the mixin holds no comm-sequence ordering (SC and CG diverge here).

**RNG + external force**
- `m_rng_counter`, `get_rng_state()`, `set_rng_state(uint64_t)`, `get_seed()`.
- `m_external_force`, `set_external_force`, `get_external_force`.

**Not in the mixin — stays in each leaf**
- All `pystencils`-generated kernel ownership and the `integrate()` loop body.
- `get_node_density` / `set_node_density` / `get_node_velocity` / `set_node_velocity` (different field shapes; `set_node_velocity` throws in `LBWalberlaImplColorGradient` — fixing it is deferred).
- `get_node_population` / `set_node_population` and slice equivalents (Q differs; color-gradient has two PDF fields).
- `get_node_pressure_tensor`, `get_pressure_tensor`, `get_momentum`.
- `set_collision_model(kT, seed)` and `set_collision_model(LeesEdwardsPack)` — only `LBWalberlaImplSingleComponent` supports them; `LBWalberlaImplColorGradient` overrides both to throw.
- **Force-at-position setters** (`add_force_at_pos`, `add_forces_at_pos`, `add_density_weighted_forces_at_pos`): in the leaves, not the mixin. Single-component has one force field; color-gradient has two per-component force fields, and the barycentric splitting rule is leaf-local. Each leaf uses `BSplineWeights::compute` from the mixin and writes into its own field(s).
- For `LBWalberlaImplColorGradient`: additionally `set_collision_model_color_gradient(sigma, beta)`, `init_pdfs_from_components`, `get_color_gradients_at_pos`, `add_solvation_forces_at_pos`, `get_node_component_densities` and friends.

**Estimated sizes:** mixin ~550 LoC, single-component leaf ~400 LoC, color-gradient leaf ~500 LoC. The current `LBWalberlaImpl.hpp` is 1249 LoC.

## 5. Interface split (Phase 3)

### 5.1 `LBWalberlaBase` (slimmed)

Deletes: `has_two_components`, `set_collision_model_two_component`, `init_two_component`, `get_color_gradients_at_pos`, `add_solvation_forces_at_pos`.

Narrows (§2.3 of `lb_architecture_strategies.md`):
- `get_node_density` returns `optional<double>` (was `optional<vector<double>>`).
- `set_node_density(Vector3i, double)` (was `vector<double> const &`).
- `get_viscosity()` returns `double`.
- `set_viscosity(double)`.

Slice accessors (`get_slice_density`, `get_slice_population`, etc.) stay `vector<double>` — bulk APIs are naturally flat.

Per-node viscosity is *not* added — node-level viscosity is not implemented in espresso and is not planned short-term. Viscosity remains a fluid-wide scalar on `LBWalberlaBase`.

### 5.2 `LBWalberlaColorGradientBase`

Inherits `LBWalberlaBase`. Adds:

```cpp
virtual void set_collision_model_color_gradient(double sigma, double beta) = 0;
virtual void init_pdfs_from_components() = 0;

virtual std::vector<Utils::Vector3d>
get_color_gradients_at_pos(std::vector<Utils::Vector3d> const &positions) = 0;
virtual void
add_solvation_forces_at_pos(std::vector<Utils::Vector3d> const &positions,
                            std::vector<double> const &delta_mus) = 0;

// Per-component density (inherited get_node_density returns total ρ)
virtual std::optional<std::array<double, 2>>
get_node_component_densities(Utils::Vector3i const &,
                             bool consider_ghosts = false) const = 0;
virtual bool
set_node_component_densities(Utils::Vector3i const &,
                             std::array<double, 2> const &) = 0;
virtual std::vector<double>
get_slice_component_densities(Utils::Vector3i const &,
                              Utils::Vector3i const &) const = 0;

// Per-component viscosity
virtual void set_component_viscosities(std::array<double, 2> const &) = 0;
[[nodiscard]] virtual std::array<double, 2>
get_component_viscosities() const = 0;
```

### 5.3 Python deprecation shim

Only `density` needs a shim (viscosity is fluid-scope, not per-node).

Approach:
- `LBFluidNode.density` returns a `float` directly for the natural Python path.
- If the user does `lbf[x,y,z].density[0]` (legacy index-into-scalar pattern), the property returns a `_LegacyScalarWrapper` that emits `DeprecationWarning("indexing into a scalar single-component LB density is deprecated; in a future release this will be a plain float")` on `__getitem__`, but still returns the value.
- Removal of the wrapper is **out of scope for this PR**; it is documented as a planned follow-up after one release.
- For color-gradient: `LBFluidNode.density` returns the barycentric total (scalar). Per-component access uses the new `component_densities` property returning a length-2 array.
- The wrapper class lives in `src/python/espressomd/lb.py` (no separate `_deprecations.py` for one wrapper).

New test: `testsuite/python/lb_density_api_deprecation.py` asserting (a) scalar read works, (b) `[0]` indexing emits `DeprecationWarning` and still returns the value, (c) the wrapper does not measurably slow down bare scalar reads (basic ratio assertion).

## 6. Dispatch sites (Phase 4)

### 6.1 `core/lb/Solver`

```cpp
struct Implementation {
  std::shared_ptr<LBWalberlaBase> base;       // owns the bridge
  LBWalberlaColorGradientBase *color_gradient = nullptr;  // non-owning,
                                                          // == base.get() iff color-gradient
};
```

Set at construction via `dynamic_cast<LBWalberlaColorGradientBase*>(base.get())`. Never mutated.
Replace `has_two_components()` with `bool has_color_gradient() const noexcept { return m_impl->color_gradient != nullptr; }` for callers that want a bool, and with direct `m_impl->color_gradient` use elsewhere.

### 6.2 `core/lb/particle_coupling.cpp`

```cpp
if (auto *cg = m_solver.color_gradient()) {
  auto grads = cg->get_color_gradients_at_pos(positions);
  cg->add_solvation_forces_at_pos(positions, delta_mus);
}
```

The previous `has_two_components()` guard is replaced; the conditional doubles as the typed entry point.

### 6.3 `core/lb/LBWalberla` (wrapper)

Keep one wrapper class `LBWalberla` holding `shared_ptr<LBWalberlaBase>`. Add `LBWalberlaColorGradientBase *color_gradient()` accessor that returns the cached down-cast. **Do not** split the wrapper — that is the Strategy C path.

### 6.4 `core/lb/LBNone`

Delete `get_color_gradients_at_pos` and `has_two_components` (both currently throw `runtime_error`). The fact that `LBNone` does not inherit `LBWalberlaColorGradientBase` is what expresses "no color-gradient when LB is disabled" — at the type level rather than at runtime.

### 6.5 `script_interface/walberla/LBFluid.cpp::make_instance()`

```cpp
if (viscosity.size() == 2) {
  m_lb_fluid = std::make_shared<LBWalberlaImplColorGradient<F, A>>(...);
} else {
  m_lb_fluid = std::make_shared<LBWalberlaImplSingleComponent<F, A>>(...);
}
```

This is the only runtime branch in script-interface. Everywhere else, script-interface holds `shared_ptr<LBWalberlaBase>` and dispatches via virtual calls.

### 6.6 `script_interface/walberla/LBFluid.cpp` — color-gradient-only methods

`sigma`, `beta`, `init_two_component` script-interface methods stay on the unified script-interface class. They `dynamic_cast` (or use a cached `LBWalberlaColorGradientBase *m_color_gradient`) to dispatch. They throw a clear error if invoked on a single-component instance:
`"this LB instance is single-component; sigma is only valid for two-component (color-gradient) LB"`.

### 6.7 `LBFluidNode` and `LBFluidSlice`

No structural change. Forward the now-scalar density value to Python; the deprecation wrapper lives in `lb.py`, not in the C++ script-interface layer.

## 7. Commit sequence

Each commit builds independently and the verification battery passes after each.

1. **Baseline snapshot.** Run the verification battery (§8) on `color_gradient` HEAD. Record timings. No code change.
2. **Phase 1.** Add `LBWalberlaCommon.hpp`. Make existing `LBWalberlaImpl` inherit from it. Move the ~550 LoC of shared infrastructure into the mixin. Pure motion.
3. **Phase 2a.** Add `LBWalberlaImplSingleComponent.hpp` (copy of `LBWalberlaImpl` with color-gradient arms deleted). Wire `make_instance` to instantiate it when `viscosity.size() == 1`.
4. **Phase 2b.** Add `LBWalberlaImplColorGradient.hpp` (copy of `LBWalberlaImpl` with single-component arms deleted). Wire `make_instance` to instantiate it when `viscosity.size() == 2`. Original `LBWalberlaImpl` no longer referenced.
5. **Phase 2c.** Delete `LBWalberlaImpl.hpp` and the `m_two_components` member. This commit removes the ~68 internal branches.
6. **Phase 3a.** Add `LBWalberlaColorGradientBase.hpp`. Make `LBWalberlaImplColorGradient` inherit from it. Move the 5 CG-only virtuals (plus the per-component accessors) from `LBWalberlaBase` to `LBWalberlaColorGradientBase`. Keep `has_two_components()` *temporarily* on `LBWalberlaBase` as a thin shim (false in SC override, true in CG override) so the rest of the code keeps compiling.
7. **Phase 3b.** Narrow density/viscosity accessors on `LBWalberlaBase` to scalar. Add `get_node_component_densities` to `LBWalberlaColorGradientBase`. Update both leaves. Add the Python `_LegacyScalarWrapper` deprecation shim. Add `testsuite/python/lb_density_api_deprecation.py`.
8. **Phase 4a.** Add `Solver::color_gradient()` accessor + cached pointer. Replace `has_two_components()` callers in `core/lb/` with the new accessor. `LBNone` loses `has_two_components` and `get_color_gradients_at_pos`. `particle_coupling.cpp` uses the typed pointer.
9. **Phase 4b.** Delete `has_two_components()` from `LBWalberlaBase`. The script-interface `LBFluid` calls the typed CG methods via cached `m_color_gradient`. Improve throw messages.
10. **Final.** Update or remove `lb_architecture_strategies.md` to reflect what landed.

## 8. Verification battery

After each commit:

- `make -j8` build clean.
- `testsuite/python/lb.py`
- `testsuite/python/lb_stats.py`
- `testsuite/python/lb_color_gradient.py`
- `testsuite/python/lb_color_gradient_particle_coupling.py`
- `testsuite/python/lb_planar_couette.py`
- `testsuite/python/lb_interpolation.py`
- `testsuite/python/lb_lees_edwards.py`

After commits 7 and 8 (user-visible Python surface changes), also:

- `testsuite/python/lb_mass_conservation.py`
- `testsuite/python/engine_lb.py`

## 9. Risks and mitigations

- **Subtle regression in CG integration loop ordering during the Phase 2 split.** Commits 4 and 5 explicitly only diff one model's arm at a time; visual diff review per commit.
- **`dynamic_cast` in `Solver` construction returns nullptr under `-fno-rtti`.** Espresso already relies on RTTI elsewhere; verify with `grep` for `-fno-rtti` in the CMake setup. If present, fall back to a virtual `is_color_gradient()` returning a typed pointer.
- **Deprecation shim adds measurable per-call overhead.** The wrapper is created only on `__getitem__` access; bare numerical use takes the scalar path with no overhead. The deprecation test includes a basic ratio assertion.
- **A dispatch site (besides Solver / particle_coupling / script-interface) is missed.** Before deleting `has_two_components()` in commit 9, `git grep has_two_components` across the full tree. Compile failure after deletion acts as a second enforcement layer.

## 10. Rollback

Each commit is independently buildable, so any single commit can be `git revert`-ed without cascading. The branch is in a worktree (`/tikhome/weeber/es-strategy-a`) and can be deleted wholesale without affecting other workspaces.

## 11. Open questions deferred to the implementation plan

- Whether to keep a cached `m_color_gradient` pointer in the script-interface `LBFluid` class itself, or always `dynamic_cast` at each color-gradient method call. (Construction-time cache is preferred but small; this is a writing-plans-level decision.)
- The exact shape of the `LBWalberlaColorGradientBase` per-component slice accessor return type — `vector<double>` interleaved vs. `vector<array<double,2>>`. Today's `get_slice_density` is interleaved; matching that is the default.
- Whether the `_LegacyScalarWrapper` raises `DeprecationWarning` once per access site (using `warnings.warn(..., stacklevel=3)`) or once globally. Default: per access site, stacklevel matched to the user's caller.
