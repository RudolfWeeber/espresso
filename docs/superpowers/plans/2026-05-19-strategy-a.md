# Strategy A — LB Single-Component / Color-Gradient Split Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor the single `LBWalberlaImpl<FloatType, Architecture>` class into two sibling leaves (`LBWalberlaImplSingleComponent`, `LBWalberlaImplColorGradient`) sharing infrastructure via a CRTP mixin (`LBWalberlaCommon`), behind a slimmed `LBWalberlaBase` and a new `LBWalberlaColorGradientBase` extension interface. Replace all `has_two_components()` dispatch with a typed `solver.color_gradient()` accessor. Land Phases 1–4 of Strategy A.

**Architecture:** CRTP mixin + sibling leaves + slim base + extension base. Core layer dispatches via a `LBWalberlaColorGradientBase*` cached at construction (`nullptr` ⇔ single-component). One narrow Python-visible behaviour change (`density` accessor narrowed from `vector<double>` to `double`), guarded with a `DeprecationWarning` shim.

**Tech Stack:** C++20, waLBerla, espresso script-interface, pybind / espresso's call_method bridge, Python testsuite, CMake, make.

**Source spec:** `docs/superpowers/specs/2026-05-19-strategy-a-design.md` (commit `dd1014b10d`).

**Working tree:** `/tikhome/weeber/es-strategy-a/` (a sibling git worktree of `/tikhome/weeber/es`). All paths below are relative to that worktree root. Branch: `strategy-a`, branched from `color_gradient`.

**Verification battery** (run after every commit; defined once here, referenced from each task):

```bash
# Build (always from build/)
make -j8

# Python testsuite — minimum gate
ctest -R '^lb$|^lb_stats$|^lb_color_gradient$|^lb_color_gradient_particle_coupling$|^lb_planar_couette$|^lb_interpolation$|^lb_lees_edwards$' --output-on-failure

# Extended gate (after Task 7 and Task 8 specifically)
ctest -R '^lb_mass_conservation$|^engine_lb$' --output-on-failure
```

Expected: all tests pass. Build clean (no new warnings).

---

## Task 0: Set up the strategy-a worktree and capture baseline

**Files:**
- Create worktree: `/tikhome/weeber/es-strategy-a/`
- No source files modified.

- [ ] **Step 1: Create the worktree and branch**

From `/tikhome/weeber/es`:

```bash
git -C /tikhome/weeber/es worktree add /tikhome/weeber/es-strategy-a -b strategy-a color_gradient
```

Expected output: `Preparing worktree (new branch 'strategy-a' from 'color_gradient')` followed by `HEAD is now at <sha>`. Verify with:

```bash
git -C /tikhome/weeber/es worktree list
```

You should see both `/tikhome/weeber/es` (on `color_gradient`) and `/tikhome/weeber/es-strategy-a` (on `strategy-a`). **All subsequent steps run inside `/tikhome/weeber/es-strategy-a`.**

- [ ] **Step 2: Create a fresh build directory inside the worktree**

```bash
mkdir -p /tikhome/weeber/es-strategy-a/build
cd /tikhome/weeber/es-strategy-a/build
cmake .. <same flags as the main build — copy from /tikhome/weeber/es/build/CMakeCache.txt if unsure>
```

If you don't know the cmake invocation, look at `/tikhome/weeber/es/build/CMakeCache.txt` and extract the non-default cache variables. **Do not** copy or share build artifacts between worktrees.

- [ ] **Step 3: Run the baseline verification battery**

From `/tikhome/weeber/es-strategy-a/build`:

```bash
make -j8
ctest -R '^lb$|^lb_stats$|^lb_color_gradient$|^lb_color_gradient_particle_coupling$|^lb_planar_couette$|^lb_interpolation$|^lb_lees_edwards$|^lb_mass_conservation$|^engine_lb$' --output-on-failure
```

Expected: all tests pass. Record this list as the green baseline. If anything fails, **stop and fix the baseline before proceeding** — every later task's verification compares against this baseline.

- [ ] **Step 4: No commit yet — Task 0 is setup only.**

The branch is on `color_gradient`'s commit. Subsequent tasks introduce commits.

---

## Task 1: Phase 1 — extract `LBWalberlaCommon` CRTP mixin

**Files:**
- Create: `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaCommon.hpp`
- Modify: `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImpl.hpp` (inherit from mixin; move members and methods out)
- Modify: `src/walberla_bridge/src/lattice_boltzmann/LBBoundaryAccess.impl.hpp` (likely qualify methods on the mixin)
- Modify: `src/walberla_bridge/src/lattice_boltzmann/LBInterpolation.impl.hpp` (likely partial move into mixin for read paths)
- No file deletions.

**Goal of this commit:** pure code motion. `LBWalberlaImpl` continues to satisfy the same interface, instantiated identically, and `m_two_components` and all CG branches stay where they are. Only the **shared infrastructure** identified in spec §4 moves into the CRTP mixin.

- [ ] **Step 1: Read the current state**

Read these files top-to-bottom before changing anything:

```
src/walberla_bridge/include/walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp
src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImpl.hpp
src/walberla_bridge/src/lattice_boltzmann/LBBoundaryAccess.impl.hpp
src/walberla_bridge/src/lattice_boltzmann/LBInterpolation.impl.hpp
src/walberla_bridge/src/lattice_boltzmann/LBNodeAccess.impl.hpp
src/walberla_bridge/src/lattice_boltzmann/LBSliceAccess.impl.hpp
src/walberla_bridge/src/lattice_boltzmann/LBCollisionSetup.impl.hpp
src/walberla_bridge/src/lattice_boltzmann/LBVTK.impl.hpp
src/walberla_bridge/src/lattice_boltzmann/lb_fields.hpp
src/walberla_bridge/src/lattice_boltzmann/lb_kernels.hpp
src/walberla_bridge/src/lattice_boltzmann/ResetForce.hpp
src/walberla_bridge/src/lattice_boltzmann/InterpolateAndShiftAtBoundary.hpp
```

This is mandatory — the spec gives categories, but the precise method-level placement decisions live in these files.

- [ ] **Step 2: Create `LBWalberlaCommon.hpp` skeleton**

Create `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaCommon.hpp` with this skeleton:

```cpp
/*
 * Copyright (C) 2010-2026 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file
 * Shared waLBerla scaffolding for both single-component
 * (LBWalberlaImplSingleComponent) and color-gradient
 * (LBWalberlaImplColorGradient) lattice-Boltzmann leaves.
 * Used as a CRTP base; the leaf's typed field IDs are reached via
 * static_cast<Derived &>(*this).
 */

#pragma once

#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>

// add includes for: BlockForest, StructuredBlockStorage, GhostLayerField,
// RegularFullCommunicator, FlagField, UBB field types, Utils::Vector, etc.
// Copy-paste from the corresponding includes at the top of LBWalberlaImpl.hpp.

#include <memory>

namespace walberla {
namespace lbm {

template <class Derived, typename FloatType, lbmpy::Arch Architecture>
class LBWalberlaCommon : public LBWalberlaBase {
protected:
  LBWalberlaCommon(std::shared_ptr<LatticeWalberla> lattice)
      : LBWalberlaBase(lattice) {}

  Derived &derived() noexcept { return static_cast<Derived &>(*this); }
  Derived const &derived() const noexcept {
    return static_cast<Derived const &>(*this);
  }

  // === Storage / lattice geometry (members moved from LBWalberlaImpl) ===
  // m_lattice and accessors — note LBWalberlaBase already owns m_lattice;
  // do not duplicate the member, just expose helpers if needed.

  // === RNG + external force ===
  uint64_t m_rng_counter{0};
  Utils::Vector3d m_external_force{};

public:
  // === Traits (moved verbatim from LBWalberlaImpl) ===
  [[nodiscard]] bool is_double_precision() const noexcept override {
    return std::is_same_v<FloatType, double>;
  }
  [[nodiscard]] bool is_gpu() const noexcept override {
    return Architecture == lbmpy::Arch::GPU;
  }
  // stencil_size(): moved verbatim.
  // make_lattice_position_checker: moved verbatim.

  // === Off-lattice interpolation (READ paths) ===
  // get_velocity_at_pos / get_velocities_at_pos: moved verbatim. They access
  // the velocity field via derived().m_velocity_field_id (CRTP).
  // get_density_at_pos / get_densities_at_pos: moved verbatim, returning the
  // total ρ. Per-component positional accessors stay on LBWalberlaColorGradientBase.

  // === Boundary / UBB scaffolding ===
  // reallocate_ubb_field, clear_boundaries, update_boundary_from_shape,
  // get_boundary_force_from_shape, get_boundary_force,
  // get_node_velocity_at_boundary, set_node_velocity_at_boundary,
  // get_slice_velocity_at_boundary, set_slice_velocity_at_boundary,
  // get_node_boundary_force, remove_node_from_boundary,
  // get_node_is_boundary, get_slice_is_boundary
  // — moved verbatim; they refer to derived().m_boundary_handling_id,
  // derived().m_force_field_id, etc.

  // === External force getter/setter ===
  void set_external_force(Utils::Vector3d const &f) override {
    m_external_force = f;
  }
  Utils::Vector3d get_external_force() const noexcept override {
    return m_external_force;
  }

  // === RNG state ===
  [[nodiscard]] std::optional<uint64_t> get_rng_state() const override {
    // moved verbatim
  }
  void set_rng_state(uint64_t v) override { m_rng_counter = v; }
  unsigned int get_seed() const noexcept override {
    // moved verbatim — returns derived().m_seed if seed lives in leaf,
    // or the cached value here if hoisted.
  }
};

} // namespace lbm
} // namespace walberla
```

This is intentionally a skeleton — fill in the actual method bodies by **moving** them out of `LBWalberlaImpl.hpp` (cut, not copy). The spec §4 "stays in each leaf" list tells you what NOT to move.

**Important:** force-at-position setters (`add_force_at_pos`, `add_forces_at_pos`, `add_density_weighted_forces_at_pos`) **stay in `LBWalberlaImpl`**, not in the mixin (spec §4). Only `BSplineWeights::compute` and equivalent weight-math helpers go into the mixin if they exist as a separable function — if they're inlined into the force methods today, leave them inlined for now and revisit in Task 2.

- [ ] **Step 3: Make `LBWalberlaImpl` inherit from `LBWalberlaCommon`**

In `LBWalberlaImpl.hpp`, change the class declaration from:

```cpp
template <typename FloatType, lbmpy::Arch Architecture>
class LBWalberlaImpl : public LBWalberlaBase {
```

to:

```cpp
template <typename FloatType, lbmpy::Arch Architecture>
class LBWalberlaImpl
    : public LBWalberlaCommon<LBWalberlaImpl<FloatType, Architecture>,
                              FloatType, Architecture> {
  using Base = LBWalberlaCommon<LBWalberlaImpl<FloatType, Architecture>,
                                FloatType, Architecture>;
  friend Base; // grants Base access to private field IDs via CRTP
```

Add `#include "LBWalberlaCommon.hpp"` near the top.

Delete the moved methods from `LBWalberlaImpl`. Verify nothing is duplicated.

- [ ] **Step 4: Adjust `.impl.hpp` files**

For each method that you moved into the mixin, its out-of-line definition lives in one of the `.impl.hpp` files. Move that definition into the mixin's body as well (or, if you prefer, leave the body in the `.impl.hpp` but change the qualifier from `LBWalberlaImpl<F,A>::name` to `LBWalberlaCommon<D,F,A>::name`). The cleaner choice is to inline into the mixin header — these methods are typically small.

If you keep the `.impl.hpp` indirection for any of these, ensure it's included from `LBWalberlaCommon.hpp` rather than from `LBWalberlaImpl.hpp`.

For methods that stay in `LBWalberlaImpl`, their `.impl.hpp` definitions stay unchanged.

- [ ] **Step 5: Build**

```bash
cd build && make -j8
```

Expected: clean build. Likely error categories to watch for:
- Missing field-ID accessor on `Derived` — fix by adding the field as a public member in `LBWalberlaImpl` or by friending `Base`.
- Method ambiguity from leftover definitions — delete the duplicate in `LBWalberlaImpl`.
- Missing `override` keyword now that base is the mixin — add `override` (the mixin re-declares the virtuals; the leaf still overrides them only if it overrides further).

- [ ] **Step 6: Run the verification battery**

```bash
ctest -R '^lb$|^lb_stats$|^lb_color_gradient$|^lb_color_gradient_particle_coupling$|^lb_planar_couette$|^lb_interpolation$|^lb_lees_edwards$' --output-on-failure
```

Expected: all green.

- [ ] **Step 7: Commit**

```bash
git -C /tikhome/weeber/es-strategy-a add \
  src/walberla_bridge/src/lattice_boltzmann/LBWalberlaCommon.hpp \
  src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImpl.hpp \
  src/walberla_bridge/src/lattice_boltzmann/LBBoundaryAccess.impl.hpp \
  src/walberla_bridge/src/lattice_boltzmann/LBInterpolation.impl.hpp

git -C /tikhome/weeber/es-strategy-a commit -m "$(cat <<'EOF'
LB: extract LBWalberlaCommon CRTP mixin from LBWalberlaImpl

Moves shared waLBerla scaffolding (lattice/storage, UBB boundary
handling, off-lattice read-side interpolation, RNG state, external
force) into a new CRTP mixin LBWalberlaCommon<Derived, F, A>.
LBWalberlaImpl now inherits from it. Pure code motion: no semantic
change, no interface change.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Phase 2a — add `LBWalberlaImplSingleComponent`

**Files:**
- Create: `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImplSingleComponent.hpp`
- Modify: `src/script_interface/walberla/LBFluid.cpp` (dispatch in `make_instance`)
- Modify: `src/walberla_bridge/src/lattice_boltzmann/lb_walberla_init.cpp` (if this file does the construction)
- Modify: `src/walberla_bridge/src/lattice_boltzmann/lb_walberla_init.cu` (CUDA arm, if applicable)
- Modify: `src/walberla_bridge/src/lattice_boltzmann/CMakeLists.txt` (no — header-only)

**Goal of this commit:** add the single-component leaf in parallel to the existing class. The original `LBWalberlaImpl` is **still present** and still handles color-gradient. After this commit, single-component instances flow through the new leaf; color-gradient still flows through `LBWalberlaImpl`.

- [ ] **Step 1: Create `LBWalberlaImplSingleComponent.hpp`**

Copy `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImpl.hpp` to `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImplSingleComponent.hpp`. Then:

1. Rename the class from `LBWalberlaImpl` to `LBWalberlaImplSingleComponent`.
2. Update the `Base` typedef: `using Base = LBWalberlaCommon<LBWalberlaImplSingleComponent<FloatType, Architecture>, FloatType, Architecture>;`
3. **Delete every color-gradient arm.** Specifically:
   - Delete `m_two_components` member.
   - Delete the `has_two_components()` override or change it to `return false;`.
   - For every `if (has_two_components()) { ... } else { ... }`: keep the `else` (single-component) arm only and unindent.
   - For every `if (has_two_components()) { ... }` with no `else`: delete the whole block.
   - Delete CG-only fields: `m_pdf_field_id_a`, `m_pdf_field_id_b`, `m_phase_field_id`, color-gradient field IDs, two-component density fields, two-component force fields, the second velocity field if present, the `ColorGradient*` kernel members, the color-gradient communicator(s).
   - Delete `init_two_component`, `set_collision_model_two_component`, `get_color_gradients_at_pos`, `add_solvation_forces_at_pos`.
   - Delete CG-only branches inside `integrate()`, `set_collision_model`, `set_viscosity`, `set_external_force` etc.
4. For methods that previously threw on single-component (e.g. `set_node_velocity` was *not* one of those — it's CG that throws — but check the codebase): they stay implemented here.

Use the audit comments (`// CG-only` / `// SC-only` in your head while reading) to drive the deletions. A useful guard: every time you see `m_two_components` or `has_two_components`, you should either delete code or simplify.

- [ ] **Step 2: Wire up the construction site**

Find the place that constructs `LBWalberlaImpl<FloatType, Architecture>` from the script-interface. It's in one of:
- `src/walberla_bridge/src/lattice_boltzmann/lb_walberla_init.cpp` (most likely)
- `src/script_interface/walberla/LBFluid.cpp`

Search:

```bash
git -C /tikhome/weeber/es-strategy-a grep -n "LBWalberlaImpl<" src/
```

Wherever you find a construction (e.g. `std::make_shared<LBWalberlaImpl<F, A>>(...)`), wrap it in a runtime dispatch on the viscosity tuple size. The exact form depends on the constructor signature — if the constructor takes `std::vector<double> const &viscosity`:

```cpp
std::shared_ptr<LBWalberlaBase> make_lb_walberla(/* ... */,
                                                 std::vector<double> const &viscosity,
                                                 /* ... */) {
  if (viscosity.size() == 1) {
    return std::make_shared<LBWalberlaImplSingleComponent<FloatType, Architecture>>(
        /* arguments minus what only-CG-needs (sigma, beta) */);
  } else {
    return std::make_shared<LBWalberlaImpl<FloatType, Architecture>>(/* unchanged */);
  }
}
```

If the construction is via a factory template, parameterise the factory or duplicate it.

- [ ] **Step 3: Build**

```bash
cd build && make -j8
```

Likely errors:
- `LBWalberlaImplSingleComponent` references a CG-only method you forgot to delete — delete it.
- The single-component constructor signature now has unused parameters (`sigma`, `beta`) — accept this for now; clean up in Task 4.

- [ ] **Step 4: Run the verification battery**

```bash
ctest -R '^lb$|^lb_stats$|^lb_color_gradient$|^lb_color_gradient_particle_coupling$|^lb_planar_couette$|^lb_interpolation$|^lb_lees_edwards$' --output-on-failure
```

Expected: all green. The color-gradient tests still pass because `LBWalberlaImpl` is still wired up for them. The single-component tests now pass through `LBWalberlaImplSingleComponent`.

If `lb.py` or `lb_planar_couette.py` fails, your deletion of CG arms was probably too aggressive in a shared method body — diff `LBWalberlaImplSingleComponent.hpp` against `LBWalberlaImpl.hpp` and verify.

- [ ] **Step 5: Commit**

```bash
git -C /tikhome/weeber/es-strategy-a add \
  src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImplSingleComponent.hpp \
  src/walberla_bridge/src/lattice_boltzmann/lb_walberla_init.cpp \
  src/walberla_bridge/src/lattice_boltzmann/lb_walberla_init.cu \
  src/script_interface/walberla/LBFluid.cpp

# Drop any .cu file if it doesn't exist on disk

git -C /tikhome/weeber/es-strategy-a commit -m "$(cat <<'EOF'
LB: add LBWalberlaImplSingleComponent leaf

Forks the single-component path out of LBWalberlaImpl into a new
sibling class LBWalberlaImplSingleComponent. The script-interface
construction site now dispatches on viscosity.size() and instantiates
the single-component leaf when size == 1. Color-gradient still
flows through the original LBWalberlaImpl until the next commit.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Phase 2b — add `LBWalberlaImplColorGradient`

**Files:**
- Create: `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImplColorGradient.hpp`
- Modify: same construction-site file(s) touched in Task 2
- The original `LBWalberlaImpl.hpp` stays for now (deleted in Task 4).

**Goal of this commit:** add the color-gradient leaf in parallel. After this commit, both leaves exist and both are wired up; `LBWalberlaImpl` is still present but no longer referenced.

- [ ] **Step 1: Create `LBWalberlaImplColorGradient.hpp`**

Copy `LBWalberlaImpl.hpp` → `LBWalberlaImplColorGradient.hpp`. Then:

1. Rename the class from `LBWalberlaImpl` to `LBWalberlaImplColorGradient`.
2. Update `using Base = LBWalberlaCommon<LBWalberlaImplColorGradient<FloatType, Architecture>, FloatType, Architecture>;`
3. **Delete every single-component arm.** Specifically:
   - Delete `m_two_components` member (it is always-true here).
   - Delete the `has_two_components()` override or change it to `return true;`.
   - For every `if (has_two_components()) { ... } else { ... }`: keep the `if` arm, unindent it, delete the `else`.
   - For every `if (!has_two_components()) { ... }` or equivalent: delete.
   - Delete SC-only kernel members if they exist (most don't — but check `m_stream_collide_kernel` etc.).
4. **Override SC-only methods to throw**:
   - `set_collision_model(double kT, unsigned int seed)`: throw `std::runtime_error("thermalized collision is not supported on the color-gradient LB model")`.
   - `set_collision_model(std::unique_ptr<LeesEdwardsPack> &&)`: throw `std::runtime_error("Lees-Edwards is not supported on the color-gradient LB model")`.
   - `check_lebc(...)`: throw the same.
   - `set_node_velocity(...)`: keep the existing throw (it already throws in the current CG-arm code).
5. CG-only public methods (`set_collision_model_two_component`, `init_two_component`, `get_color_gradients_at_pos`, `add_solvation_forces_at_pos`) become first-class methods of this class. **Rename** `set_collision_model_two_component` → `set_collision_model_color_gradient`. **Rename** `init_two_component` → `init_pdfs_from_components`. **Search the codebase** for callers of the old names (script-interface, possibly tests) and rename them too:

```bash
git -C /tikhome/weeber/es-strategy-a grep -n "set_collision_model_two_component\|init_two_component" src/
```

Update each caller.

- [ ] **Step 2: Wire up the color-gradient construction**

In the dispatch from Task 2, change the `else` branch to construct `LBWalberlaImplColorGradient` instead of `LBWalberlaImpl`:

```cpp
if (viscosity.size() == 1) {
  return std::make_shared<LBWalberlaImplSingleComponent<FloatType, Architecture>>(/* ... */);
} else {
  return std::make_shared<LBWalberlaImplColorGradient<FloatType, Architecture>>(/* ... */);
}
```

The original `LBWalberlaImpl` should no longer appear in any call site.

- [ ] **Step 3: Build**

```bash
cd build && make -j8
```

Watch for: callers of the old method names you renamed but missed in Step 1. Fix each.

- [ ] **Step 4: Verification battery**

```bash
ctest -R '^lb$|^lb_stats$|^lb_color_gradient$|^lb_color_gradient_particle_coupling$|^lb_planar_couette$|^lb_interpolation$|^lb_lees_edwards$' --output-on-failure
```

Expected: all green.

If a `lb_color_gradient*` test fails, your SC-arm-deletion was too aggressive in a shared method. Diff `LBWalberlaImplColorGradient.hpp` against `LBWalberlaImpl.hpp` and verify.

- [ ] **Step 5: Commit**

```bash
git -C /tikhome/weeber/es-strategy-a add \
  src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImplColorGradient.hpp \
  src/walberla_bridge/src/lattice_boltzmann/lb_walberla_init.cpp \
  src/walberla_bridge/src/lattice_boltzmann/lb_walberla_init.cu \
  src/script_interface/walberla/LBFluid.cpp

git -C /tikhome/weeber/es-strategy-a commit -m "$(cat <<'EOF'
LB: add LBWalberlaImplColorGradient leaf

Forks the color-gradient path out of LBWalberlaImpl into a new
sibling class LBWalberlaImplColorGradient. The script-interface
construction site now dispatches to this leaf when viscosity.size() == 2.
The old LBWalberlaImpl is no longer referenced from any call site
and will be deleted in the next commit. Renames
set_collision_model_two_component -> set_collision_model_color_gradient
and init_two_component -> init_pdfs_from_components on the bridge
interface; callers updated in the same commit.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Phase 2c — delete `LBWalberlaImpl`

**Files:**
- Delete: `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImpl.hpp`
- Modify: any `.impl.hpp` file that still references `LBWalberlaImpl<F,A>::`
- Modify: clean up unused `sigma`/`beta` parameters from `LBWalberlaImplSingleComponent` constructor if you let them sit in Task 2.

**Goal of this commit:** the original is gone. `git grep` of `LBWalberlaImpl[^A-Za-z]` (i.e., not followed by an identifier character) should return zero hits in non-deleted code. This commit is the one that physically removes the ~68 internal branches.

- [ ] **Step 1: Verify `LBWalberlaImpl` is unreferenced**

```bash
git -C /tikhome/weeber/es-strategy-a grep -n "LBWalberlaImpl<" -- 'src/**'
git -C /tikhome/weeber/es-strategy-a grep -n "LBWalberlaImpl::" -- 'src/**'
```

Expected: no hits (or only hits inside `LBWalberlaImpl.hpp` itself, which is what we're deleting). If anything else turns up, that call site needs to be updated to one of the new leaves first.

- [ ] **Step 2: Migrate `.impl.hpp` files**

The `.impl.hpp` files define methods like:

```cpp
template <typename FloatType, lbmpy::Arch Architecture>
auto LBWalberlaImpl<FloatType, Architecture>::get_node_density(...) {
```

For methods that now live in **both leaves**, fork the `.impl.hpp` file:
- `LBNodeAccess.impl.hpp` → `LBNodeAccessSingleComponent.impl.hpp` + `LBNodeAccessColorGradient.impl.hpp`
- `LBSliceAccess.impl.hpp` → `LBSliceAccessSingleComponent.impl.hpp` + `LBSliceAccessColorGradient.impl.hpp`

In each forked file, change `LBWalberlaImpl` to the appropriate leaf class and delete the irrelevant arm of every `has_two_components()` branch. The fork mirrors what you did in Tasks 2 and 3, but for the out-of-line definitions.

For files that hold methods now in `LBWalberlaCommon`, change the qualifier to `LBWalberlaCommon<Derived, F, A>::` and add the appropriate `template <class Derived, ...>` parameters. Or, easier: move those method bodies inline into `LBWalberlaCommon.hpp` and delete them from the `.impl.hpp` entirely.

For each `.impl.hpp`, ensure it is now `#include`d only by the leaf that needs it (single-component leaf includes `LBNodeAccessSingleComponent.impl.hpp` at the bottom of its class definition, etc.).

- [ ] **Step 3: Delete `LBWalberlaImpl.hpp`**

```bash
git -C /tikhome/weeber/es-strategy-a rm src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImpl.hpp
```

- [ ] **Step 4: Clean up unused parameters**

If `LBWalberlaImplSingleComponent`'s constructor still accepts `sigma`/`beta` parameters that it ignores, drop them. Update the construction site to not pass them.

If the dispatch site has a way to express "I built a single-component bridge" cleanly without dragging CG parameters along, prefer that.

- [ ] **Step 5: Build**

```bash
cd build && make -j8
```

The link/build should be clean.

- [ ] **Step 6: Verification battery**

```bash
ctest -R '^lb$|^lb_stats$|^lb_color_gradient$|^lb_color_gradient_particle_coupling$|^lb_planar_couette$|^lb_interpolation$|^lb_lees_edwards$' --output-on-failure
```

Expected: all green. This is the commit that materializes the win (~68 branches gone) — verify with:

```bash
git -C /tikhome/weeber/es-strategy-a diff color_gradient..HEAD --stat
```

You should see a substantial net deletion of `has_two_components()` branches (each leaf has fewer than the original).

- [ ] **Step 7: Commit**

```bash
git -C /tikhome/weeber/es-strategy-a add -A
git -C /tikhome/weeber/es-strategy-a commit -m "$(cat <<'EOF'
LB: delete LBWalberlaImpl, completing the leaf split

Removes the original LBWalberlaImpl class and its m_two_components
runtime flag. Forks LBNodeAccess.impl.hpp and LBSliceAccess.impl.hpp
into per-leaf variants. Net effect: the ~68 has_two_components()
branches in the bridge layer are gone; each leaf carries only its
own logic. Bridge interface (LBWalberlaBase) is unchanged in this
commit; that comes next.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Phase 3a — add `LBWalberlaColorGradientBase`, move CG-only virtuals

**Files:**
- Create: `src/walberla_bridge/include/walberla_bridge/lattice_boltzmann/LBWalberlaColorGradientBase.hpp`
- Modify: `src/walberla_bridge/include/walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp`
- Modify: `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImplColorGradient.hpp`
- Modify: `src/core/lb/LBNone.hpp` (keep `has_two_components()` returning false — it stays as a shim this commit)
- Modify: `src/core/lb/LBWalberla.hpp` / `LBWalberla.cpp` (keep `has_two_components()` shim; will be deleted in Task 8)
- Modify: `src/core/lb/Solver.cpp` (keep the `std::visit` call to `has_two_components()` — will be deleted in Task 7)

**Goal of this commit:** introduce the extension interface and move the 5 CG-only virtuals from `LBWalberlaBase` to `LBWalberlaColorGradientBase`. Keep `has_two_components()` on `LBWalberlaBase` *temporarily* as a shim with a default implementation so existing call sites compile. Add the per-component density accessors to `LBWalberlaColorGradientBase`.

- [ ] **Step 1: Create `LBWalberlaColorGradientBase.hpp`**

Create `src/walberla_bridge/include/walberla_bridge/lattice_boltzmann/LBWalberlaColorGradientBase.hpp`:

```cpp
/*
 * Copyright (C) 2010-2026 The ESPResSo project
 *
 * (license boilerplate matching LBWalberlaBase.hpp)
 */

/**
 * @file
 * Color-gradient (two-component) LB extension interface.
 * Inherits the common fluid interface (LBWalberlaBase) and adds
 * model-specific configuration, initialization, and per-component
 * accessors.
 */

#pragma once

#include "LBWalberlaBase.hpp"
#include <utils/Vector.hpp>

#include <array>
#include <optional>
#include <vector>

class LBWalberlaColorGradientBase : public LBWalberlaBase {
public:
  ~LBWalberlaColorGradientBase() override = default;

  // --- CG-specific configuration ---
  virtual void set_collision_model_color_gradient(double sigma, double beta) = 0;
  virtual void init_pdfs_from_components() = 0;

  // --- CG-specific observables / forces ---
  virtual std::vector<Utils::Vector3d>
  get_color_gradients_at_pos(std::vector<Utils::Vector3d> const &positions) = 0;
  virtual void
  add_solvation_forces_at_pos(std::vector<Utils::Vector3d> const &positions,
                              std::vector<double> const &delta_mus) = 0;

  // --- per-component density accessors ---
  virtual std::optional<std::array<double, 2>>
  get_node_component_densities(Utils::Vector3i const &node,
                               bool consider_ghosts = false) const = 0;
  virtual bool
  set_node_component_densities(Utils::Vector3i const &node,
                               std::array<double, 2> const &rho) = 0;
  virtual std::vector<double>
  get_slice_component_densities(Utils::Vector3i const &lower,
                                Utils::Vector3i const &upper) const = 0;

  // --- per-component viscosity ---
  virtual void set_component_viscosities(std::array<double, 2> const &nu) = 0;
  [[nodiscard]] virtual std::array<double, 2>
  get_component_viscosities() const = 0;
};
```

Use the same namespace / formatting as the existing `LBWalberlaBase.hpp` (no namespace at file scope in espresso, typically).

- [ ] **Step 2: Move CG-only virtuals out of `LBWalberlaBase`**

In `LBWalberlaBase.hpp`, delete:
- `set_collision_model_two_component`
- `init_two_component`
- `get_color_gradients_at_pos`
- `add_solvation_forces_at_pos`

Keep `has_two_components()` for one more commit — add a default implementation:

```cpp
[[nodiscard]] virtual bool has_two_components() const noexcept { return false; }
```

The CG leaf overrides this; SC leaf inherits the default. `LBNone` already returns false. Mark the method as deprecated in the doxygen comment:

```cpp
/**
 * @deprecated Use a dynamic_cast to LBWalberlaColorGradientBase instead.
 *             Removed after Task 8.
 */
```

- [ ] **Step 3: Update `LBWalberlaImplColorGradient`**

Change the class declaration from inheriting `LBWalberlaCommon<Self, F, A>` (which inherits `LBWalberlaBase`) to:

```cpp
template <typename FloatType, lbmpy::Arch Architecture>
class LBWalberlaImplColorGradient
    : public LBWalberlaCommon<LBWalberlaImplColorGradient<FloatType, Architecture>,
                              FloatType, Architecture>,
      public virtual LBWalberlaColorGradientBase {
```

This introduces multiple inheritance: the mixin and the extension interface. The mixin already inherits `LBWalberlaBase`, and `LBWalberlaColorGradientBase` also inherits `LBWalberlaBase`, so use `virtual` inheritance on the path that matters. Adjust as the compiler tells you — likely:
- `LBWalberlaCommon` inherits `LBWalberlaBase` virtually.
- `LBWalberlaColorGradientBase` inherits `LBWalberlaBase` virtually.

Modify `LBWalberlaCommon.hpp` to `class LBWalberlaCommon : public virtual LBWalberlaBase {`.

Add the new per-component virtuals (`get_node_component_densities`, `set_node_component_densities`, `get_slice_component_densities`, `set_component_viscosities`, `get_component_viscosities`) to `LBWalberlaImplColorGradient` as overrides. Implement them by **adapting** existing methods that today return `vector<double>` with a 2-element layout: copy the body, replace return-type / arg-type with `array<double, 2>` and call the inherited `vector`-returning method internally if it's faster to compose than to duplicate. (Don't delete the `vector`-returning methods yet — Phase 3b will deal with the type narrowing.)

Make `has_two_components()` on the CG leaf `return true;` (override).

- [ ] **Step 4: Update LBNone, LBWalberla, Solver to compile**

Since `has_two_components()` still exists on `LBWalberlaBase`, no changes needed to LBNone, LBWalberla.hpp/cpp, Solver.hpp/cpp for this commit. Verify with:

```bash
git -C /tikhome/weeber/es-strategy-a grep -n "has_two_components" src/core/ src/script_interface/
```

These callers continue working.

- [ ] **Step 5: Build**

```bash
cd build && make -j8
```

Watch for:
- Diamond inheritance ambiguity — fix with `virtual` inheritance per Step 3.
- The script-interface `LBFluid.cpp` calls `set_collision_model_two_component` (old name) — replace with `set_collision_model_color_gradient` and add a `dynamic_cast<LBWalberlaColorGradientBase*>(...)` to access it. Pattern:
```cpp
auto *cg = dynamic_cast<LBWalberlaColorGradientBase*>(m_instance.get());
if (cg == nullptr) {
  throw std::runtime_error(
      "this LB instance is single-component; sigma is only valid for "
      "two-component (color-gradient) LB");
}
cg->set_collision_model_color_gradient(sigma, beta);
```

Apply the same pattern to `init_two_component`, `get_color_gradients_at_pos`, `add_solvation_forces_at_pos` callers.

- [ ] **Step 6: Verification battery**

```bash
ctest -R '^lb$|^lb_stats$|^lb_color_gradient$|^lb_color_gradient_particle_coupling$|^lb_planar_couette$|^lb_interpolation$|^lb_lees_edwards$' --output-on-failure
```

Expected: all green.

- [ ] **Step 7: Commit**

```bash
git -C /tikhome/weeber/es-strategy-a add -A
git -C /tikhome/weeber/es-strategy-a commit -m "$(cat <<'EOF'
LB: add LBWalberlaColorGradientBase extension interface

Introduces a new abstract base LBWalberlaColorGradientBase inheriting
LBWalberlaBase, holding the 5 CG-only virtuals plus per-component
density/viscosity accessors. LBWalberlaImplColorGradient now also
inherits LBWalberlaColorGradientBase (virtual inheritance on
LBWalberlaBase). The five color-gradient virtuals are removed from
the common base; has_two_components() is kept on LBWalberlaBase as a
deprecated shim for one more commit. Script-interface code that
called the removed methods now dynamic_casts to
LBWalberlaColorGradientBase and throws a clear error on
single-component instances.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Phase 3b — narrow density/viscosity types, add Python deprecation shim

**Files:**
- Modify: `src/walberla_bridge/include/walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp`
- Modify: `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImplSingleComponent.hpp`
- Modify: `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImplColorGradient.hpp`
- Modify: `src/walberla_bridge/src/lattice_boltzmann/LBNodeAccessSingleComponent.impl.hpp` (and CG variant)
- Modify: `src/core/lb/LBWalberla.hpp` / `LBWalberla.cpp`
- Modify: `src/core/lb/LBNone.hpp`
- Modify: `src/script_interface/walberla/LBFluidNode.hpp` / `.cpp`
- Modify: `src/script_interface/walberla/LBFluidSlice.hpp`
- Modify: `src/python/espressomd/lb.py` (add `_LegacyScalarWrapper`)
- Create: `testsuite/python/lb_density_api_deprecation.py`
- Modify: `testsuite/python/CMakeLists.txt` (register new test)

**Goal of this commit:** narrow `get_node_density` / `set_node_density` / `get_viscosity` / `set_viscosity` on `LBWalberlaBase` from vector to scalar; add per-component density accessors (already declared on `LBWalberlaColorGradientBase` in Task 5 — now implement); add the Python `_LegacyScalarWrapper` shim that emits `DeprecationWarning`; add a new test that pins both the new behavior and the deprecation path.

- [ ] **Step 1: Write the failing deprecation test first (TDD)**

Create `testsuite/python/lb_density_api_deprecation.py`:

```python
#
# Copyright (C) 2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# (license boilerplate)
#

import unittest as ut
import unittest_decorators as utx
import warnings
import espressomd
import espressomd.lb


@utx.skipIfMissingFeatures("WALBERLA")
class TestLBDensityScalarDeprecation(ut.TestCase):

    def setUp(self):
        self.system = espressomd.System(box_l=[8.0, 8.0, 8.0])
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.1
        self.lbf = espressomd.lb.LBFluidWalberla(
            agrid=1.0, density=1.0, kinematic_viscosity=1.0, tau=0.01)
        self.system.lb = self.lbf

    def tearDown(self):
        self.system.lb = None

    def test_single_component_density_returns_float(self):
        """Reading lbf[node].density returns a plain float."""
        rho = self.lbf[0, 0, 0].density
        # Bare float comparison must work without indexing.
        self.assertAlmostEqual(float(rho), 1.0, places=10)

    def test_legacy_indexing_emits_deprecation_warning(self):
        """lbf[node].density[0] still returns the value, but warns."""
        rho = self.lbf[0, 0, 0].density
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            value = rho[0]
            self.assertAlmostEqual(float(value), 1.0, places=10)
            deprecation_warnings = [
                w for w in caught if issubclass(w.category, DeprecationWarning)
            ]
            self.assertGreater(len(deprecation_warnings), 0,
                               "expected a DeprecationWarning on density[0]")
            self.assertIn("density", str(deprecation_warnings[0].message))


@utx.skipIfMissingFeatures("WALBERLA")
class TestLBColorGradientDensityAccessor(ut.TestCase):

    def setUp(self):
        self.system = espressomd.System(box_l=[8.0, 8.0, 8.0])
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.1
        self.lbf = espressomd.lb.LBFluidWalberla(
            agrid=1.0, density=1.0,
            kinematic_viscosity=[1.0, 1.0],  # two-component
            tau=0.01, sigma=0.0, beta=0.7)
        self.system.lb = self.lbf

    def tearDown(self):
        self.system.lb = None

    def test_total_density_is_scalar(self):
        """On color-gradient, lbf[node].density is the barycentric total."""
        rho = self.lbf[0, 0, 0].density
        self.assertAlmostEqual(float(rho), 1.0, places=10)

    def test_component_densities_is_length_two(self):
        """lbf[node].component_densities returns the per-component split."""
        rho_ab = self.lbf[0, 0, 0].component_densities
        self.assertEqual(len(rho_ab), 2)


if __name__ == "__main__":
    ut.main()
```

Register the test by adding `lb_density_api_deprecation.py` to `testsuite/python/CMakeLists.txt` next to the other `lb_*.py` entries.

- [ ] **Step 2: Run the new test, confirm it fails**

```bash
cd build && make -j8
ctest -R '^lb_density_api_deprecation$' --output-on-failure
```

Expected: FAIL. Likely failure: `AttributeError` on `.component_densities` (doesn't exist yet), or the density is returned as `[1.0]` not `1.0`.

If the test fails for the **wrong reason** (e.g. `LBFluidWalberla` constructor signature changed), fix that first.

- [ ] **Step 3: Narrow the C++ interface**

In `src/walberla_bridge/include/walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp`, change:

```cpp
virtual std::optional<std::vector<double>>
get_node_density(Utils::Vector3i const &node, bool consider_ghosts = false) const = 0;
virtual bool set_node_density(Utils::Vector3i const &node,
                              std::vector<double> const &rho) = 0;
virtual void set_viscosity(std::vector<double> const &nu) = 0;
[[nodiscard]] virtual std::vector<double> get_viscosity() const = 0;
```

to:

```cpp
virtual std::optional<double>
get_node_density(Utils::Vector3i const &node, bool consider_ghosts = false) const = 0;
virtual bool set_node_density(Utils::Vector3i const &node, double rho) = 0;
virtual void set_viscosity(double nu) = 0;
[[nodiscard]] virtual double get_viscosity() const = 0;
```

(Keep slice accessors as `std::vector<double>` — they are intrinsically flat.)

In each leaf:
- `LBWalberlaImplSingleComponent`: change signatures to match; bodies likely already returned single-element vectors, so unwrap to scalar.
- `LBWalberlaImplColorGradient`: `get_node_density` returns the **barycentric total** (sum of component densities); `set_node_density(double)` is implemented by splitting the new value equally? **Decision point:** the spec doesn't fully pin this. The safest behaviour: `set_node_density(double rho)` throws (`"on color-gradient LB, use set_node_component_densities to set per-component densities"`), and `get_node_density` returns the total. Implement both this way unless reading the existing two-component path suggests a different default.

Implement `get_node_component_densities` / `set_node_component_densities` / `get_slice_component_densities` / `set_component_viscosities` / `get_component_viscosities` on `LBWalberlaImplColorGradient` (declared in Task 5). Internally, these wrap the existing vector-of-2 internals.

- [ ] **Step 4: Update bridge consumers**

Now that the type changed, update everywhere that calls `get_node_density` / `set_node_density` / `get_viscosity` / `set_viscosity` on the bridge:

```bash
git -C /tikhome/weeber/es-strategy-a grep -n "->get_node_density\|->set_node_density\|->get_viscosity\|->set_viscosity" src/
```

Likely files: `src/core/lb/LBWalberla.cpp`, `src/script_interface/walberla/LBFluidNode.cpp`, `src/script_interface/walberla/LBFluidSlice.hpp`. Update each:
- `LBWalberla::get_node_density` (if it exists) — change return type and unwrap.
- `LBFluidNode::density` (script-interface): the variant-conversion now wraps a `double` instead of a 1-element vector. If LBFluidNode needs to deliver a length-2 list for the color-gradient case, route through `dynamic_cast<LBWalberlaColorGradientBase*>(m_lb_fluid.get())->get_node_component_densities(...)`. Otherwise it delivers a scalar.

In `src/script_interface/walberla/LBFluidSlice.hpp:129`, the lines:

```cpp
m_shape_val["density"] = {m_lb_fluid->has_two_components() ? 2 : 1};
m_shape_val["pop"] = {m_lb_fluid->has_two_components() ? ... : ...};
```

Change the density branch to:

```cpp
auto *cg = dynamic_cast<LBWalberlaColorGradientBase*>(m_lb_fluid.get());
m_shape_val["density"] = {cg ? 2 : 1};
```

(`pop` will be handled later if necessary; for density specifically, this is the right shape.)

- [ ] **Step 5: Add the Python deprecation shim in `lb.py`**

In `src/python/espressomd/lb.py`, add (near the top, after imports):

```python
import warnings as _warnings


class _LegacyScalarWrapper(float):
    """Compatibility shim for the SC single-component density API.

    Returned by the legacy code path where users may write
    ``lbf[i, j, k].density[0]``. Calling :py:meth:`__getitem__` emits
    a :class:`DeprecationWarning` but still returns the underlying
    value, so existing scripts keep working until the wrapper is
    removed in a future release.
    """

    def __getitem__(self, idx):
        _warnings.warn(
            "indexing into a scalar single-component LB density is "
            "deprecated; in a future release this will be a plain float",
            DeprecationWarning,
            stacklevel=2,
        )
        if idx == 0 or idx == -1:
            return float(self)
        raise IndexError(
            "scalar single-component LB density has only one element")

    def __len__(self):
        _warnings.warn(
            "len() on a scalar single-component LB density is deprecated; "
            "in a future release this will be a plain float",
            DeprecationWarning,
            stacklevel=2,
        )
        return 1
```

Find the `density` property on `LBFluidNode` in `lb.py`. Wrap the C++ return value:

```python
@property
def density(self):
    raw = self.call_method("get_density")
    if isinstance(raw, (int, float)):
        return _LegacyScalarWrapper(raw)
    # Color-gradient bridge returned a 2-element list — that path no
    # longer exists at the script-interface level. If we get here,
    # we'll get a list and the test will fail loudly. Return as-is
    # for forward compatibility with a future per-component scalar.
    return raw
```

Add a `component_densities` property to `LBFluidNode`:

```python
@property
def component_densities(self):
    """Per-component density for color-gradient LB. Returns a 2-array."""
    return self.call_method("get_component_densities")
```

The exact `call_method` keys must match the C++ script-interface side. If `LBFluidNode` doesn't yet expose `get_component_densities`, add it on the C++ side (parameters-less call routing to the CG leaf).

- [ ] **Step 6: Wire `get_component_densities` through the script interface**

In `src/script_interface/walberla/LBFluidNode.cpp`, add a `get_component_densities` method handler:

```cpp
} else if (name == "get_component_densities") {
  auto *cg = dynamic_cast<LBWalberlaColorGradientBase*>(m_lb_fluid.get());
  if (cg == nullptr) {
    throw std::runtime_error(
        "component_densities is only available on two-component (color-gradient) LB");
  }
  auto const rho_ab = cg->get_node_component_densities(/* coordinates */);
  // pack into a Variant of length 2; mirror how other length-2 returns are packed
}
```

Adapt to the existing dispatch style in `LBFluidNode.cpp`.

- [ ] **Step 7: Build**

```bash
cd build && make -j8
```

Likely errors:
- Slice access code paths that returned vectors no longer compile against scalar `set_node_density` — keep slice paths returning vectors (slices are intrinsically flat); only node-level signatures changed.
- A few callers in `LBWalberla.cpp` need adaptation.

- [ ] **Step 8: Run the new deprecation test, confirm it passes**

```bash
ctest -R '^lb_density_api_deprecation$' --output-on-failure
```

Expected: PASS.

- [ ] **Step 9: Run the verification battery**

```bash
ctest -R '^lb$|^lb_stats$|^lb_color_gradient$|^lb_color_gradient_particle_coupling$|^lb_planar_couette$|^lb_interpolation$|^lb_lees_edwards$|^lb_mass_conservation$|^engine_lb$' --output-on-failure
```

Expected: all green. (This is an extended-gate task because Python surface changed.)

If `lb_mass_conservation.py` or `engine_lb.py` fails, the most likely cause is a missed `vector → scalar` adaptation in the script-interface. Check `LBFluidNode` and `LBFluidSlice` adaptations again.

- [ ] **Step 10: Commit**

```bash
git -C /tikhome/weeber/es-strategy-a add -A
git -C /tikhome/weeber/es-strategy-a commit -m "$(cat <<'EOF'
LB: narrow node density/viscosity to scalar; add deprecation shim

Narrows get_node_density / set_node_density / get_viscosity /
set_viscosity on LBWalberlaBase from optional<vector<double>> /
vector<double> to optional<double> / double. Slice accessors keep
their vector return type. Implements get_node_component_densities,
set_node_component_densities, get_slice_component_densities,
set_component_viscosities, get_component_viscosities on
LBWalberlaImplColorGradient (interface declared in Task 5).

Python-visible behaviour change: lbf[i, j, k].density now returns
a plain float. To preserve backwards compatibility for scripts that
write lbf[i, j, k].density[0], wraps the scalar in a
_LegacyScalarWrapper(float) subclass that emits DeprecationWarning
on __getitem__ / __len__. New espressomd.lb.LBFluidNode property
component_densities exposes the per-component split for
color-gradient LB.

Adds testsuite/python/lb_density_api_deprecation.py asserting both
the new scalar behaviour and the deprecation path.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Phase 4a — `Solver::color_gradient()` accessor, replace `has_two_components` callers

**Files:**
- Modify: `src/core/lb/Solver.hpp` / `Solver.cpp`
- Modify: `src/core/lb/LBWalberla.hpp` / `LBWalberla.cpp`
- Modify: `src/core/lb/LBNone.hpp`
- Modify: `src/core/lb/particle_coupling.cpp`
- Modify: `src/walberla_bridge/src/lattice_boltzmann/LBInterpolation.impl.hpp` (still has 3 `has_two_components()` callsites)
- Modify: any other internal bridge callers still referencing `has_two_components()`

**Goal of this commit:** add the typed accessor `Solver::color_gradient()` returning `LBWalberlaColorGradientBase*` (cached at LB attach time). Replace every external use of `has_two_components()` with the typed accessor. Drop `get_color_gradients_at_pos` and `has_two_components` from `LBNone`. After this commit, only `has_two_components()` *declarations* remain (on `LBWalberlaBase` and `LBWalberla`); Task 8 deletes them.

- [ ] **Step 1: Add a cached pointer to `LBWalberla`**

In `src/core/lb/LBWalberla.hpp`:

```cpp
class LBWalberla : public LBImplementation {
public:
  // ...existing...
  [[nodiscard]] LBWalberlaColorGradientBase *color_gradient() noexcept {
    return m_color_gradient;
  }
  [[nodiscard]] LBWalberlaColorGradientBase const *color_gradient() const noexcept {
    return m_color_gradient;
  }

private:
  std::shared_ptr<LBWalberlaBase> lb_fluid;
  LBWalberlaColorGradientBase *m_color_gradient = nullptr;  // non-owning
};
```

In `LBWalberla.cpp`, set `m_color_gradient` in the constructor (or wherever `lb_fluid` is assigned):

```cpp
LBWalberla::LBWalberla(std::shared_ptr<LBWalberlaBase> bridge, ...)
    : lb_fluid(std::move(bridge)) {
  m_color_gradient = dynamic_cast<LBWalberlaColorGradientBase*>(lb_fluid.get());
}
```

Include `walberla_bridge/lattice_boltzmann/LBWalberlaColorGradientBase.hpp` at the top.

- [ ] **Step 2: Expose the accessor through `Solver`**

In `src/core/lb/Solver.hpp`, add:

```cpp
[[nodiscard]] LBWalberlaColorGradientBase *color_gradient() noexcept;
[[nodiscard]] LBWalberlaColorGradientBase const *color_gradient() const noexcept;
```

In `Solver.cpp`:

```cpp
LBWalberlaColorGradientBase *Solver::color_gradient() noexcept {
  return std::visit(
      [](auto &ptr) -> LBWalberlaColorGradientBase * {
        using T = std::decay_t<decltype(*ptr)>;
        if constexpr (std::is_same_v<T, LBWalberla>) {
          return ptr->color_gradient();
        } else {
          return nullptr;
        }
      },
      *impl->solver);
}
// repeat for const-overload
```

Adjust to the exact variant types used in `Solver::impl->solver` (check the existing `has_two_components` implementation for the visitor pattern).

- [ ] **Step 3: Replace external callers**

For each file listed at the top of this task, replace `m_lb.has_two_components()` / `m_lb_fluid->has_two_components()` / `solver.has_two_components()` calls with the appropriate typed pointer check. Pattern:

```cpp
// before
if (m_lb.has_two_components()) {
  m_lb.add_solvation_forces_at_pos(positions, delta_mus);
}
// after
if (auto *cg = m_lb.color_gradient()) {
  cg->add_solvation_forces_at_pos(positions, delta_mus);
}
```

For `particle_coupling.cpp` (5 callsites), the pattern simplifies the code substantially because the typed pointer carries the CG-only methods. Re-read the file after edits to verify nothing was lost.

For `LBInterpolation.impl.hpp` (3 callsites at lines 498, 580, 595): those are *inside* a leaf, not at an external boundary. The callsites should already have been removed during the Task 2/3 forks; if any remain, that's a leftover branch in the leaf — remove the now-dead conditional.

For `LBFluid.cpp:228`: replace with a `dynamic_cast` or, if the script-interface already cached a CG pointer (Task 5 added the dynamic-cast pattern), use the cached pointer.

For `LBFluidSlice.hpp:129-131`: already changed in Task 6 if you followed Step 4 there; if not, change now.

- [ ] **Step 4: Drop CG methods from `LBNone`**

In `src/core/lb/LBNone.hpp`:
- Delete `bool has_two_components() const { return false; }`
- Delete `get_color_gradients_at_pos` (it currently throws).
- `LBNone` already does **not** inherit from `LBWalberlaColorGradientBase`, so `Solver::color_gradient()` correctly returns nullptr for it.

- [ ] **Step 5: Build**

```bash
cd build && make -j8
```

Watch for remaining `has_two_components()` callers — should be zero after this commit (except the declaration on `LBWalberlaBase` and `LBWalberla::has_two_components()`).

- [ ] **Step 6: Verification battery**

```bash
ctest -R '^lb$|^lb_stats$|^lb_color_gradient$|^lb_color_gradient_particle_coupling$|^lb_planar_couette$|^lb_interpolation$|^lb_lees_edwards$|^lb_mass_conservation$|^engine_lb$' --output-on-failure
```

Expected: all green. The particle_coupling changes in `lb_color_gradient_particle_coupling.py` and the solvation path are the most likely regression points — review the diff against the previous task closely if anything fails.

- [ ] **Step 7: Commit**

```bash
git -C /tikhome/weeber/es-strategy-a add -A
git -C /tikhome/weeber/es-strategy-a commit -m "$(cat <<'EOF'
LB: add Solver::color_gradient() typed accessor; rewire callers

Caches a non-owning LBWalberlaColorGradientBase* in core/lb/LBWalberla
(set via dynamic_cast at LB attach time) and exposes it through
Solver::color_gradient(). Replaces all external has_two_components()
callers in core/lb/ and script_interface/walberla/ with the typed
accessor. Color-gradient-only branches in particle_coupling.cpp now
use the typed pointer directly. LBNone loses get_color_gradients_at_pos
and has_two_components — its non-inheritance from
LBWalberlaColorGradientBase expresses "no color-gradient when LB is
disabled" at the type level.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Phase 4b — delete `has_two_components()` from interface and wrapper

**Files:**
- Modify: `src/walberla_bridge/include/walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp`
- Modify: `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImplSingleComponent.hpp`
- Modify: `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImplColorGradient.hpp`
- Modify: `src/core/lb/LBWalberla.hpp` / `LBWalberla.cpp`
- Modify: `src/core/lb/Solver.hpp` / `Solver.cpp`
- Modify: `src/script_interface/walberla/LBFluid.cpp` — improve throw messages

**Goal of this commit:** delete `has_two_components()` everywhere. The compiler is the second enforcement layer — any remaining caller turns into a compile error. Improve throw messages on CG-only operations attempted against single-component instances.

- [ ] **Step 1: Verify no remaining callers**

```bash
git -C /tikhome/weeber/es-strategy-a grep -n "has_two_components" src/
```

Expected output: only declarations on `LBWalberlaBase.hpp`, `LBWalberla.hpp`, `LBWalberla.cpp`, `Solver.hpp`, `Solver.cpp`, the override in `LBWalberlaImplColorGradient.hpp`. If there are external callers, you missed them in Task 7 — go back.

- [ ] **Step 2: Delete the declarations**

In `LBWalberlaBase.hpp`: delete `[[nodiscard]] virtual bool has_two_components() const noexcept;` (the default-implementation line from Task 5).

In `LBWalberlaImplColorGradient.hpp`: delete the `bool has_two_components() const noexcept override { return true; }` override.

In `LBWalberlaImplSingleComponent.hpp`: if it has an explicit override (should not after Task 5; it would have inherited the default), delete it.

In `LBWalberla.hpp` / `.cpp`: delete `has_two_components()` declaration and definition.

In `Solver.hpp` / `Solver.cpp`: delete `has_two_components()` declaration and `std::visit` body.

- [ ] **Step 3: Tighten throw messages**

In `src/script_interface/walberla/LBFluid.cpp`, wherever a CG-only operation is attempted against a single-component instance (the `dynamic_cast` returning `nullptr` from Task 5), make sure the throw message is informative:

```cpp
auto *cg = m_color_gradient;  // cached at construction
if (cg == nullptr) {
  throw std::runtime_error(
      "this LB instance is single-component; '<operation>' is only valid "
      "for two-component (color-gradient) LB");
}
```

Fill in `<operation>` with the actual name (`sigma`, `beta`, `init_two_component`, `get_color_gradients`, `add_solvation_forces`).

Consider adding a cached `LBWalberlaColorGradientBase *m_color_gradient = nullptr;` field on the script-interface `LBFluid` class itself (set via dynamic_cast in `make_instance`), so each operation doesn't do a fresh `dynamic_cast`. Open question §11.1 of the spec — implement this here.

- [ ] **Step 4: Build**

```bash
cd build && make -j8
```

The compiler should pass cleanly. Any leftover caller surfaces as a compile error.

- [ ] **Step 5: Verification battery (extended)**

```bash
ctest -R '^lb$|^lb_stats$|^lb_color_gradient$|^lb_color_gradient_particle_coupling$|^lb_planar_couette$|^lb_interpolation$|^lb_lees_edwards$|^lb_mass_conservation$|^engine_lb$' --output-on-failure
```

Expected: all green.

- [ ] **Step 6: Commit**

```bash
git -C /tikhome/weeber/es-strategy-a add -A
git -C /tikhome/weeber/es-strategy-a commit -m "$(cat <<'EOF'
LB: delete has_two_components() — typed dispatch only

Removes has_two_components() from LBWalberlaBase, the
LBWalberla wrapper, the LBWalberlaImplColorGradient override, and
Solver. The compiler now enforces that color-gradient operations go
through Solver::color_gradient() / the typed extension interface;
any forgotten caller turns into a compile error rather than a
runtime branch. Hardens the script-interface error messages on
CG-only operations attempted against single-component instances.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Finalize — update / archive the architecture-strategies doc

**Files:**
- Modify or delete: `lb_architecture_strategies.md` (in repo root)

**Goal of this commit:** the architecture-strategies doc fulfilled its purpose. Either delete it or archive it under `docs/`.

- [ ] **Step 1: Decide**

If the team uses the doc for future onboarding / context on why Strategy A was chosen, move it: `git mv lb_architecture_strategies.md docs/lb-architecture-strategies.md` and add a header noting "implemented in PR <number>".

If the doc is purely the deliberation artefact, delete it: `git rm lb_architecture_strategies.md`.

Either is defensible. **Default: archive** (`git mv` to `docs/`) — the history of the decision is useful for future readers.

- [ ] **Step 2: Commit**

```bash
git -C /tikhome/weeber/es-strategy-a add -A
git -C /tikhome/weeber/es-strategy-a commit -m "$(cat <<'EOF'
docs: archive LB architecture strategies decision document

The deliberation between Strategy A and Strategy C is preserved in
docs/lb-architecture-strategies.md for future readers. Strategy A
was implemented in this branch (Phases 1-4); Phases 5-6 and the §5
follow-up gaps (thermalized CG, Lees-Edwards CG, VTK splits,
checkpoint versioning, C++ unit tests for CG) remain as separate
work.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 3: Final verification battery**

```bash
make -j8
ctest -R '^lb$|^lb_stats$|^lb_color_gradient$|^lb_color_gradient_particle_coupling$|^lb_planar_couette$|^lb_interpolation$|^lb_lees_edwards$|^lb_mass_conservation$|^engine_lb$|^lb_density_api_deprecation$' --output-on-failure
```

Expected: all green. This is the final gate on the branch.

- [ ] **Step 4: Summarize the branch**

```bash
git -C /tikhome/weeber/es-strategy-a log --oneline color_gradient..HEAD
git -C /tikhome/weeber/es-strategy-a diff color_gradient..HEAD --stat
```

Read the summary. You should see ~9 focused commits, a net deletion of internal branches, the new mixin and leaf headers, and the new deprecation test.

---

## Notes for the executor

- **Build parallelism:** always `make -j8` (per the user's standing preference). Do not use `-j$(nproc)`.
- **Worktree hygiene:** all work happens in `/tikhome/weeber/es-strategy-a`. Do not inspect or modify `/tikhome/weeber/es-pip-install` or other sibling worktrees.
- **Naming:** full words throughout. No `SC` / `CG` abbreviations in new identifiers (use `SingleComponent` / `ColorGradient` for class names, `single_component` / `color_gradient` for variable names).
- **Commit discipline:** every task ends in exactly one commit. Don't fold tasks together; the bisect-friendliness of one-task-per-commit is the whole point.
- **If verification fails:** stop and fix before moving on. Each task assumes the previous task's verification passed.
- **If a step can't be completed as written:** check the spec (`docs/superpowers/specs/2026-05-19-strategy-a-design.md`) for intent, then adapt. The spec is the contract; this plan is the execution path.
- **Open questions** (§11 of the spec) are decided inline in this plan: (1) cache the CG pointer in script-interface — done in Task 8; (2) `vector<double>` interleaved for slice — keep current behaviour; (3) `_LegacyScalarWrapper` warns per access site with `stacklevel=2` — done in Task 6.
