# Centralized Active-Features State Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `System::active_features` (issue #5410): one component that answers "is feature X in use right now?", with particle-derived state cached in one MPI-reduced bitmask and system-derived state queried live.

**Architecture:** A new `ActiveFeatures : System::Leaf<ActiveFeatures>` holds one `unsigned` bitmask filled by a single Kokkos-parallel sweep over local particles plus one MPI all-reduce (bitwise OR). The sweep absorbs `System::update_used_propagations()`, so one dirty flag and one sweep serve both the propagation mask and the feature mask. System-derived predicates (`orientation_ghosts_needed()`, `has_gay_berne()`) are live methods reading sibling state through `get_system()` — no stored copies.

**Tech Stack:** C++20, Kokkos (host execution space) via `reduce_over_local_particles`, Boost.MPI, Boost.Test unit tests, ESPResSo python testsuite.

**Spec:** `docs/superpowers/specs/2026-08-31-active-features-design.md`

**Deviations from the spec (approved improvements, flagged during plan review):**

1. The rotation-propagation arm of `orientation_ghosts_needed()` gains `and particles_can_rotate()`. Reason: `Propagation::update_default_propagation()` puts `ROT_EULER`/`ROT_LANGEVIN` into `default_propagation` for every NVT/BD run when `ESPRESSO_ROTATION` is compiled in, and every default particle carries `SYSTEM_DEFAULT`. Today that means QUAT push and TORQUE reduce on *every* ghost exchange of a plain LJ simulation. All rotational integrators skip particles with `can_rotate() == false` (`Particle` default is `rotation = 0b000`), so the tightened arm is exact.
2. `dipolar_pair_kernel_active` in `forces.cpp` gains `and particles_have_dipole_moment()`, mirroring the tightened dipole arm of `orientation_ghosts_needed()`. This keeps the invariant asserted at `forces.cpp:541` (`torque-scatter marking implies GHOSTTRANS_TORQUE`) intact: both sides tighten together.
3. `System::on_thermostat_param_change()` also sets the dirty flag. The spec listed the invalidation sites as unchanged, but the integrate-start update becomes flag-guarded where the old one was unconditional, and the thermostat switch feeds `default_propagation` — without this invalidation a thermostat change could leave a stale `SYSTEM_DEFAULT` expansion.

## Global Constraints

- Working directory for every command: `/tikhome/weeber/es/.claude/worktrees/active_feautes` (a git worktree on branch `worktree-active_feautes`). Never touch `/tikhome/weeber/es` itself or its other worktrees.
- Use `git -C /tikhome/weeber/es/.claude/worktrees/active_feautes ...`; do not chain `cd <path> && git ...`.
- Build with `make -j8` (never `-j$(nproc)`; the machine is shared).
- Feature macros in C++ code use the `ESPRESSO_` prefix (`ESPRESSO_ROTATION`); myconfig files use unprefixed names (`ROTATION`).
- Identifiers use full words, not abbreviations.
- A pre-commit hook runs clang-format/autopep8 on `git commit`. If the commit fails because the hook reformatted files, `git add` the reformatted files and run the same commit command again.
- End every commit message with: `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`
- This is a shared machine: run tests, but no performance measurements.

---

### Task 1: `ActiveFeatures` component, System wiring, unit tests

Creates the class with the feature bitmask, the sweep, and the MPI reduction. No production call sites yet (those come with Task 2), so this task is self-contained and safe. The class carries its own dirty flag `m_recalc` in this task; Task 2 replaces it with the shared flag on `Propagation`.

**Files:**
- Create: `src/core/system/ActiveFeatures.hpp`
- Create: `src/core/system/ActiveFeatures.cpp`
- Modify: `src/core/CMakeLists.txt` (source list, next to `system/System.cpp` at line ~51)
- Modify: `src/core/system/System.hpp` (forward declaration ~line 79, member ~line 328)
- Modify: `src/core/system/System.cpp` (ctor ~line 83, `initialize()` ~line 120)
- Modify: `src/core/system/System.impl.hpp` (include list)
- Create: `src/core/unit_tests/ActiveFeatures_test.cpp`
- Modify: `src/core/unit_tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `System::Leaf` (`src/core/system/Leaf.hpp`), `reduce_over_local_particles` (`src/core/particle_reduction.hpp`), `Particle` accessors, `::comm_cart`.
- Produces (used by every later task):
  - `class System::ActiveFeatures` with member `std::shared_ptr<ActiveFeatures> System::System::active_features`.
  - `void update()` — unconditional collective recompute.
  - `void update_if_needed()` — recompute only when dirty; collective.
  - `void invalidate()` / `bool needs_update() const` — transitional (removed in Task 2).
  - Queries (all `const`, all read the cached mask):
    `particles_are_virtual()`, `particles_can_rotate()`, `particles_are_fixed()`, `particles_have_charge()` (unconditional);
    `particles_have_exclusions()` [`ESPRESSO_EXCLUSIONS`], `particles_have_stoner_wohlfarth()` [`ESPRESSO_THERMAL_STONER_WOHLFARTH`], `particles_have_dipole_moment()` [`ESPRESSO_DIPOLES`], `particles_have_ext_force()` [`ESPRESSO_EXTERNAL_FORCES`], `particles_have_ext_torque()` [`ESPRESSO_EXTERNAL_FORCES` and `ESPRESSO_ROTATION`], `particles_are_swimmers()` [`ESPRESSO_ENGINE`]. Each returns `bool`.

- [ ] **Step 1: Configure the build (once for the whole plan)**

```bash
mkdir -p build
cp maintainer/configs/maxset.hpp build/myconfig.hpp
cmake -B build -D CMAKE_BUILD_TYPE=RelWithAssert -D ESPRESSO_BUILD_WITH_CUDA=OFF
```

The maxset config enables every feature this plan touches (`THERMAL_STONER_WOHLFARTH`, `EXCLUSIONS`, `ROTATION`, `DIPOLES`, `ENGINE`, `EXTERNAL_FORCES`, `GAY_BERNE`, `VIRTUAL_SITES_RELATIVE`). If cmake fails because a myconfig feature requires an external dependency that is not installed (e.g. SCAFACOS, GSL, walBerla-only features), delete that `#define` from `build/myconfig.hpp` and re-run the cmake command — but keep the features listed above.

- [ ] **Step 2: Write the failing unit test**

Create `src/core/unit_tests/ActiveFeatures_test.cpp`:

```cpp
/*
 * Copyright (C) 2026 The ESPResSo project
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

#define BOOST_TEST_MODULE ActiveFeatures test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "EspressoCoreGlobalConfig.hpp"
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "exclusions.hpp"
#include "particle_node.hpp"
#include "system/ActiveFeatures.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi.hpp>

struct GlobalConfig : public EspressoCoreGlobalConfig {
  GlobalConfig() {
    auto system = System::System::create();
    system->set_box_l(Utils::Vector3d{10., 10., 10.});
    system->set_cell_structure_topology(CellStructureType::REGULAR);
    ::System::set_system(system);
  }
  ~GlobalConfig() { ::System::reset_system(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

/** Apply @p mutation on the rank that owns the particle, then mark the
 *  cached state stale on all ranks (collective). */
static void mutate_particle(int p_id, auto &&mutation) {
  auto &system = System::get_system();
  auto *p = system.cell_structure->get_local_particle(p_id);
  if (p != nullptr and not p->is_ghost()) {
    mutation(*p);
  }
  system.active_features->invalidate();
}

/** Remove all particles between test cases. */
struct ParticleCleanup {
  ~ParticleCleanup() { ::remove_all_particles(); }
};

BOOST_AUTO_TEST_SUITE(suite)

BOOST_FIXTURE_TEST_CASE(default_particle_sets_no_bits, ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  active_features.invalidate();
  active_features.update_if_needed();
  BOOST_CHECK(not active_features.needs_update());
  BOOST_CHECK(not active_features.particles_are_virtual());
  BOOST_CHECK(not active_features.particles_can_rotate());
  BOOST_CHECK(not active_features.particles_are_fixed());
  BOOST_CHECK(not active_features.particles_have_charge());
#ifdef ESPRESSO_EXCLUSIONS
  BOOST_CHECK(not active_features.particles_have_exclusions());
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  BOOST_CHECK(not active_features.particles_have_stoner_wohlfarth());
#endif
#ifdef ESPRESSO_DIPOLES
  BOOST_CHECK(not active_features.particles_have_dipole_moment());
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  BOOST_CHECK(not active_features.particles_have_ext_force());
#endif
#ifdef ESPRESSO_ENGINE
  BOOST_CHECK(not active_features.particles_are_swimmers());
#endif
}

#ifdef ESPRESSO_EXCLUSIONS
BOOST_FIXTURE_TEST_CASE(exclusion_bit_is_reduced_across_ranks,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  // The particle lives on exactly one rank; the bit must be visible on all.
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  mutate_particle(0, [](Particle &p) { add_exclusion(p, 42); });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.particles_have_exclusions());
  // Removing the last exclusion clears the bit again.
  mutate_particle(0, [](Particle &p) { delete_exclusion(p, 42); });
  active_features.update_if_needed();
  BOOST_CHECK(not active_features.particles_have_exclusions());
}
#endif

#ifdef ESPRESSO_ROTATION
BOOST_FIXTURE_TEST_CASE(update_if_needed_honors_dirty_flag, ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  active_features.update();
  BOOST_CHECK(not active_features.particles_can_rotate());
  // Mutate without invalidating: the stale value must survive.
  auto *p = system.cell_structure->get_local_particle(0);
  if (p != nullptr and not p->is_ghost()) {
    p->set_can_rotate_all_axes();
  }
  active_features.update_if_needed();
  BOOST_CHECK(not active_features.particles_can_rotate());
  // After invalidation the update must pick the change up.
  active_features.invalidate();
  active_features.update_if_needed();
  BOOST_CHECK(active_features.particles_can_rotate());
}
#endif

BOOST_FIXTURE_TEST_CASE(property_bits_follow_particle_state, ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  mutate_particle(0, [](Particle &p) {
    p.stoner_wohlfarth_is_enabled() = true;
  });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.particles_have_stoner_wohlfarth());
#endif
#ifdef ESPRESSO_DIPOLES
  mutate_particle(0, [](Particle &p) { p.dipm() = 1.5; });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.particles_have_dipole_moment());
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  mutate_particle(0, [](Particle &p) {
    p.ext_force() = Utils::Vector3d{1., 0., 0.};
    p.set_fixed_along(0, true);
  });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.particles_have_ext_force());
  BOOST_CHECK(active_features.particles_are_fixed());
#endif
#ifdef ESPRESSO_ENGINE
  mutate_particle(0, [](Particle &p) { p.swimming().swimming = true; });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.particles_are_swimmers());
#endif
  // Removing the carrier particle clears every bit.
  ::remove_particle(0);
  active_features.invalidate();
  active_features.update_if_needed();
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  BOOST_CHECK(not active_features.particles_have_stoner_wohlfarth());
#endif
#ifdef ESPRESSO_DIPOLES
  BOOST_CHECK(not active_features.particles_have_dipole_moment());
#endif
  BOOST_CHECK(not active_features.particles_are_fixed());
}

BOOST_AUTO_TEST_SUITE_END()
```

Register it in `src/core/unit_tests/CMakeLists.txt`, next to the `HaloExchange_test` block (line ~104), using the same 1-and-N-rank pattern:

```cmake
espresso_unit_test_executable(
  NAME ActiveFeatures_test SRC ActiveFeatures_test.cpp DEPENDS espresso::core
  Kokkos::kokkos Boost::mpi MPI::MPI_CXX)
foreach(TEST_NUM_PROC 1 4)
  if(${TEST_NUM_PROC} LESS_EQUAL ${ESPRESSO_TEST_NP})
    espresso_unit_test_register(
      NAME ActiveFeatures_test_${TEST_NUM_PROC}_mpi_ranks TARGET
      ActiveFeatures_test NUM_PROC ${TEST_NUM_PROC})
  endif()
endforeach()
```

- [ ] **Step 3: Run the test to verify it fails**

```bash
make -C build -j8 ActiveFeatures_test
```

Expected: compilation FAILS with `system/ActiveFeatures.hpp: No such file or directory`.

- [ ] **Step 4: Create the `ActiveFeatures` class**

Create `src/core/system/ActiveFeatures.hpp` (use the standard GPL header block from `src/core/system/Leaf.hpp`, copyright `2026 The ESPResSo project`):

```cpp
#pragma once

#include <config/config.hpp>

#include "system/Leaf.hpp"

struct Particle;

namespace System {

/**
 * @brief Centralized runtime feature state.
 *
 * Answers "is feature X in use right now?".
 * Particle-derived predicates read a cached bitmask that is filled by a
 * single sweep over local particles and reduced across all MPI ranks with
 * a bitwise OR. System-derived predicates read live system state through
 * @ref Leaf::get_system and store no copies.
 *
 * Queries return the state as of the last update. They never trigger the
 * reduction lazily: the reduction is collective, and query sites are not
 * guaranteed to be collective.
 */
class ActiveFeatures : public Leaf<ActiveFeatures> {
public:
  /** @brief Recompute the particle-derived state if invalidated.
   *  Collective call: all ranks must enter together. */
  void update_if_needed();
  /** @brief Recompute the particle-derived state unconditionally.
   *  Collective call: all ranks must enter together. */
  void update();
  /** @brief Mark the particle-derived state as stale. */
  void invalidate() { m_recalc = true; }
  bool needs_update() const { return m_recalc; }

  bool particles_are_virtual() const {
    return (m_particle_features & IS_VIRTUAL) != 0u;
  }
  bool particles_can_rotate() const {
    return (m_particle_features & CAN_ROTATE) != 0u;
  }
  bool particles_are_fixed() const {
    return (m_particle_features & IS_FIXED) != 0u;
  }
  bool particles_have_charge() const {
    return (m_particle_features & HAS_CHARGE) != 0u;
  }
#ifdef ESPRESSO_EXCLUSIONS
  bool particles_have_exclusions() const {
    return (m_particle_features & HAS_EXCLUSIONS) != 0u;
  }
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  bool particles_have_stoner_wohlfarth() const {
    return (m_particle_features & HAS_STONER_WOHLFARTH) != 0u;
  }
#endif
#ifdef ESPRESSO_DIPOLES
  bool particles_have_dipole_moment() const {
    return (m_particle_features & HAS_DIPOLE_MOMENT) != 0u;
  }
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  bool particles_have_ext_force() const {
    return (m_particle_features & HAS_EXT_FORCE) != 0u;
  }
#ifdef ESPRESSO_ROTATION
  bool particles_have_ext_torque() const {
    return (m_particle_features & HAS_EXT_TORQUE) != 0u;
  }
#endif
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_ENGINE
  bool particles_are_swimmers() const {
    return (m_particle_features & IS_SWIMMER) != 0u;
  }
#endif

private:
  enum ParticleFeature : unsigned {
    IS_VIRTUAL = 1u << 0u,
    CAN_ROTATE = 1u << 1u,
    IS_FIXED = 1u << 2u,
    HAS_CHARGE = 1u << 3u,
#ifdef ESPRESSO_EXCLUSIONS
    HAS_EXCLUSIONS = 1u << 4u,
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
    HAS_STONER_WOHLFARTH = 1u << 5u,
#endif
#ifdef ESPRESSO_DIPOLES
    HAS_DIPOLE_MOMENT = 1u << 6u,
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    HAS_EXT_FORCE = 1u << 7u,
#ifdef ESPRESSO_ROTATION
    HAS_EXT_TORQUE = 1u << 8u,
#endif
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_ENGINE
    IS_SWIMMER = 1u << 9u,
#endif
  };

  static unsigned particle_bits(Particle const &p);

  /** Bitwise OR of @ref ParticleFeature over all particles on all ranks. */
  unsigned m_particle_features = 0u;
  bool m_recalc = true;
};

} // namespace System
```

Note the forward declaration `struct Particle;` between the includes and `namespace System` in the block above: `Particle` is a global type that appears only in the private static declaration, and the tag must be `struct` (clang warns on mismatched tags).

Create `src/core/system/ActiveFeatures.cpp` (same GPL header):

```cpp
#include <config/config.hpp>

#include "ActiveFeatures.hpp"

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "particle_reduction.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <functional>

namespace System {

unsigned ActiveFeatures::particle_bits(Particle const &p) {
  auto features = 0u;
  if (p.is_virtual()) {
    features |= IS_VIRTUAL;
  }
  if (p.can_rotate()) {
    features |= CAN_ROTATE;
  }
  if (p.has_fixed_coordinates()) {
    features |= IS_FIXED;
  }
  if (p.q() != 0.) {
    features |= HAS_CHARGE;
  }
#ifdef ESPRESSO_EXCLUSIONS
  if (not p.exclusions().empty()) {
    features |= HAS_EXCLUSIONS;
  }
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  if (p.stoner_wohlfarth_is_enabled()) {
    features |= HAS_STONER_WOHLFARTH;
  }
#endif
#ifdef ESPRESSO_DIPOLES
  if (p.dipm() != 0.) {
    features |= HAS_DIPOLE_MOMENT;
  }
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  if (p.ext_force() != Utils::Vector3d{0., 0., 0.}) {
    features |= HAS_EXT_FORCE;
  }
#ifdef ESPRESSO_ROTATION
  if (p.ext_torque() != Utils::Vector3d{0., 0., 0.}) {
    features |= HAS_EXT_TORQUE;
  }
#endif
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_ENGINE
  if (p.swimming().swimming) {
    features |= IS_SWIMMER;
  }
#endif
  return features;
}

void ActiveFeatures::update() {
  auto const &system = get_system();
  auto const local_features = reduce_over_local_particles<unsigned>(
      *system.cell_structure,
      [](unsigned &acc, Particle const &p) { acc |= particle_bits(p); },
      [](unsigned &acc, unsigned const &other) { acc |= other; });
  m_particle_features = boost::mpi::all_reduce(::comm_cart, local_features,
                                               std::bit_or<unsigned>());
  m_recalc = false;
}

void ActiveFeatures::update_if_needed() {
  if (m_recalc) {
    update();
  }
}

} // namespace System
```

- [ ] **Step 5: Wire the component into `System`**

In `src/core/system/System.hpp`, inside `namespace System` directly before `class System` (line ~81):

```cpp
class ActiveFeatures;
```

Add the member next to `propagation` (line ~328):

```cpp
  std::shared_ptr<ActiveFeatures> active_features;
```

In `src/core/system/System.impl.hpp`, add to the include list (alphabetical, before `"BoxGeometry.hpp"` at the bottom group or with the other `system/` includes — match surrounding style):

```cpp
#include "system/ActiveFeatures.hpp"
```

In `src/core/system/System.cpp`, in the ctor after `propagation = std::make_shared<Propagation>();` (line ~83):

```cpp
  active_features = std::make_shared<ActiveFeatures>();
```

In `initialize()` after `cell_structure->bind_system(handle);` (line ~117):

```cpp
  active_features->bind_system(handle);
```

In `src/core/CMakeLists.txt`, add to the source list next to `system/System.cpp` (line ~51):

```cmake
  system/ActiveFeatures.cpp
```

- [ ] **Step 6: Build and run the test**

```bash
make -C build -j8 ActiveFeatures_test
ctest --test-dir build -R ActiveFeatures --output-on-failure
```

Expected: both registered tests (1 and 4 ranks) PASS.

- [ ] **Step 7: Build the rest of the core to catch integration breakage**

```bash
make -C build -j8 espresso_core
```

Expected: clean build (the new header is included by `System.impl.hpp`, so `System.cpp` must compile).

- [ ] **Step 8: Commit**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/active_feautes add src/core/system/ActiveFeatures.hpp src/core/system/ActiveFeatures.cpp src/core/CMakeLists.txt src/core/system/System.hpp src/core/system/System.cpp src/core/system/System.impl.hpp src/core/unit_tests/ActiveFeatures_test.cpp src/core/unit_tests/CMakeLists.txt
git -C /tikhome/weeber/es/.claude/worktrees/active_feautes commit -m "core: add ActiveFeatures component (#5410)

One MPI-reduced bitmask answers 'does any particle have X'.
No production consumers yet.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: Absorb `update_used_propagations()` into the ActiveFeatures sweep

One sweep now fills both the feature mask and `Propagation::used_propagations`; one shared dirty flag (renamed on `Propagation`) serves both. The old update function and the transitional `m_recalc` disappear.

**Files:**
- Modify: `src/core/system/ActiveFeatures.hpp` (drop `m_recalc`/`invalidate`/`needs_update`)
- Modify: `src/core/system/ActiveFeatures.cpp` (sweep both masks, one all-reduce)
- Modify: `src/core/integrators/Propagation.hpp:35,49` (rename flag)
- Modify: `src/core/system/System.hpp:312-315` (remove `update_used_propagations` declaration)
- Modify: `src/core/system/System.cpp:288` (`on_thermostat_param_change`), `:347`, `:363` (renamed flag), `:376-380` (`update_dependent_particles`), `:410-413` (`on_observable_calc`)
- Modify: `src/core/integrate.cpp:189-201` (delete `update_used_propagations`), `:633-634` (integrate start)
- Modify: `src/core/collision_detection/BindAtPointOfCollision.cpp:191`, `src/core/collision_detection/GlueToSurface.cpp:207`
- Modify: `src/core/unit_tests/ActiveFeatures_test.cpp` (adapt Task-1 cases, add propagation cases)

**Interfaces:**
- Consumes: `Propagation::update_default_propagation(int thermo_switch)`, `Propagation::used_propagations`, `Propagation::default_propagation`, `PropagationMode::SYSTEM_DEFAULT`, `Thermostat::thermo_switch`.
- Produces:
  - `Propagation::recalc_active_features` (bool, replaces `recalc_used_propagations`; set by `Propagation::set_integ_switch()`, `System::on_particle_change()`, `System::on_particle_local_change()`, `System::on_thermostat_param_change()`; cleared only by `ActiveFeatures::update()`).
  - `ActiveFeatures::update()` now also calls `propagation.update_default_propagation(thermo_switch)` and stores `propagation.used_propagations`.
  - `System::System::update_used_propagations()` no longer exists; `ActiveFeatures::invalidate()`/`needs_update()` no longer exist.

- [ ] **Step 1: Extend the unit test (failing first)**

In `src/core/unit_tests/ActiveFeatures_test.cpp`:

Add includes:

```cpp
#include "PropagationMode.hpp"
#include "integrators/Propagation.hpp"
#include "thermostat.hpp"
```

Replace every `invalidate()` call on the active-features object with a direct write to the flag: in `mutate_particle`, the line `system.active_features->invalidate();` becomes `system.propagation->recalc_active_features = true;`; in test bodies, `active_features.invalidate();` becomes `system.propagation->recalc_active_features = true;`. Replace `active_features.needs_update()` with `system.propagation->recalc_active_features`.

Add these cases before `BOOST_AUTO_TEST_SUITE_END()`:

```cpp
BOOST_FIXTURE_TEST_CASE(used_propagations_from_same_sweep, ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  auto &propagation = *system.propagation;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  mutate_particle(0, [](Particle &p) {
    p.propagation() = PropagationMode::TRANS_LANGEVIN;
  });
  active_features.update_if_needed();
  BOOST_CHECK(propagation.used_propagations &
              PropagationMode::TRANS_LANGEVIN);
  BOOST_CHECK((propagation.used_propagations &
               PropagationMode::SYSTEM_DEFAULT) == 0);
}

BOOST_FIXTURE_TEST_CASE(system_default_expansion_and_thermostat_invalidation,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  auto &propagation = *system.propagation;
  // A default particle carries SYSTEM_DEFAULT, which expands to
  // default_propagation, which depends on the thermostat switch.
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  system.propagation->recalc_active_features = true;
  active_features.update_if_needed();
  BOOST_CHECK(propagation.used_propagations & PropagationMode::TRANS_NEWTON);
  BOOST_CHECK(not propagation.recalc_active_features);
  // Changing the thermostat must invalidate and change the expansion.
  system.thermostat->thermo_switch = THERMO_LANGEVIN;
  system.on_thermostat_param_change();
  BOOST_CHECK(propagation.recalc_active_features);
  active_features.update_if_needed();
  BOOST_CHECK(propagation.used_propagations &
              PropagationMode::TRANS_LANGEVIN);
#ifdef ESPRESSO_ROTATION
  BOOST_CHECK(propagation.used_propagations & PropagationMode::ROT_LANGEVIN);
#endif
  // Restore for other test cases.
  system.thermostat->thermo_switch = THERMO_OFF;
  system.on_thermostat_param_change();
}

BOOST_FIXTURE_TEST_CASE(set_integ_switch_invalidates, ParticleCleanup) {
  auto &system = System::get_system();
  system.active_features->update();
  BOOST_CHECK(not system.propagation->recalc_active_features);
  system.propagation->set_integ_switch(INTEG_METHOD_NVT);
  BOOST_CHECK(system.propagation->recalc_active_features);
}

BOOST_FIXTURE_TEST_CASE(particle_creation_invalidates, ParticleCleanup) {
  auto &system = System::get_system();
  system.active_features->update();
  BOOST_CHECK(not system.propagation->recalc_active_features);
  // make_new_particle fires on_particle_change, which must set the flag.
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  BOOST_CHECK(system.propagation->recalc_active_features);
}
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
make -C build -j8 ActiveFeatures_test
```

Expected: compilation FAILS (`recalc_active_features` does not exist yet).

- [ ] **Step 3: Rename the flag on `Propagation`**

In `src/core/integrators/Propagation.hpp`, replace line 35:

```cpp
  bool recalc_used_propagations = true;
```

with:

```cpp
  /** If true, the particle-derived active-feature state, including
   *  @ref used_propagations, will be recomputed at the next collective
   *  update point. Cleared only by System::ActiveFeatures::update(). */
  bool recalc_active_features = true;
```

and in `set_integ_switch()` (line 49), replace `recalc_used_propagations = true;` with `recalc_active_features = true;`.

- [ ] **Step 4: Move the propagation update into the sweep**

In `src/core/system/ActiveFeatures.hpp`: delete `invalidate()`, `needs_update()`, and the `bool m_recalc = true;` member. Update the class doc comment's staleness paragraph to name the invalidation flag:

```
 * The cached state is invalidated via Propagation::recalc_active_features
 * (set on particle changes, integrator changes, and thermostat changes)
 * and refreshed at collective update points: the start of integrate(),
 * on_observable_calc(), update_dependent_particles(), and after
 * collision-detection topology changes.
```

In `src/core/system/ActiveFeatures.cpp`, add includes:

```cpp
#include "PropagationMode.hpp"
#include "integrators/Propagation.hpp"
#include "thermostat.hpp"
```

Replace `update()` and `update_if_needed()` with:

```cpp
void ActiveFeatures::update() {
  auto &system = get_system();
  auto &propagation = *system.propagation;
  propagation.update_default_propagation(system.thermostat->thermo_switch);
  struct SweepResult {
    int propagations = PropagationMode::NONE;
    unsigned features = 0u;
  };
  auto const local = reduce_over_local_particles<SweepResult>(
      *system.cell_structure,
      [](SweepResult &acc, Particle const &p) {
        acc.propagations |= p.propagation();
        acc.features |= particle_bits(p);
      },
      [](SweepResult &acc, SweepResult const &other) {
        acc.propagations |= other.propagations;
        acc.features |= other.features;
      });
  int const local_masks[2] = {local.propagations,
                              static_cast<int>(local.features)};
  int global_masks[2];
  boost::mpi::all_reduce(::comm_cart, local_masks, 2, global_masks,
                         std::bit_or<int>());
  auto used_propagations = global_masks[0];
  if (used_propagations & PropagationMode::SYSTEM_DEFAULT) {
    used_propagations |= propagation.default_propagation;
  }
  propagation.used_propagations = used_propagations;
  m_particle_features = static_cast<unsigned>(global_masks[1]);
  propagation.recalc_active_features = false;
}

void ActiveFeatures::update_if_needed() {
  if (get_system().propagation->recalc_active_features) {
    update();
  }
}
```

(`update()` now takes the non-const `get_system()`; both methods stay non-const members.)

- [ ] **Step 5: Delete the old update function and redirect all call sites**

In `src/core/integrate.cpp`, delete the whole function `System::System::update_used_propagations()` (lines 189-201).

In `src/core/integrate.cpp` in `System::System::integrate()`, replace (lines ~633-634):

```cpp
  propagation.update_default_propagation(thermostat->thermo_switch);
  update_used_propagations();
```

with:

```cpp
  active_features->update_if_needed();
```

(`update_default_propagation` now runs inside the update; the call is flag-guarded where the old one was unconditional — the flag is set by every particle, integrator, and thermostat change.)

`integrate.cpp` needs the include `#include "system/ActiveFeatures.hpp"` — add it next to the other `system/` includes at the top of the file.

In `src/core/system/System.hpp`, delete the declaration and its doc comment (lines ~312-315):

```cpp
  /**
   * @brief Update the global propagation bitmask.
   */
  void update_used_propagations();
```

In `src/core/system/System.cpp`:

1. Lines 347 and 363: replace `propagation->recalc_used_propagations = true;` with `propagation->recalc_active_features = true;`.
2. Line 288, extend `on_thermostat_param_change` — the thermostat switch feeds `default_propagation`, so a thermostat change must invalidate the cached expansion:

```cpp
void System::on_thermostat_param_change() {
  reinit_thermo = true;
  propagation->recalc_active_features = true;
}
```

3. In `update_dependent_particles()` (lines ~376-380), replace:

```cpp
#ifdef ESPRESSO_VIRTUAL_SITES
  if (propagation->recalc_used_propagations) {
    update_used_propagations();
  }
```

with an unconditional flag-guarded update as the function's first statement, before the `#ifdef ESPRESSO_VIRTUAL_SITES` block:

```cpp
  active_features->update_if_needed();
#ifdef ESPRESSO_VIRTUAL_SITES
```

4. In `on_observable_calc()` (line ~410), add as the first statement, before the ghost update — this is a new update point; today the first ghost exchange there can read a stale propagation mask:

```cpp
  active_features->update_if_needed();
```

In `src/core/collision_detection/BindAtPointOfCollision.cpp` (line 191) and `src/core/collision_detection/GlueToSurface.cpp` (line 207), replace:

```cpp
    system.update_used_propagations();
```

with an unconditional update — collisions just created virtual particles mid-step without going through `on_particle_change()`:

```cpp
    system.active_features->update();
```

Both files need `#include "system/ActiveFeatures.hpp"` next to their other includes.

- [ ] **Step 6: Verify no stale references remain**

```bash
grep -rn "update_used_propagations\|recalc_used_propagations" src/
```

Expected: no matches.

- [ ] **Step 7: Build and run the tests**

```bash
make -C build -j8 ActiveFeatures_test espresso_core
ctest --test-dir build -R ActiveFeatures --output-on-failure
```

Expected: PASS at 1 and 4 ranks.

- [ ] **Step 8: Run neighboring unit tests that exercise propagation**

```bash
make -C build -j8 check_unit_tests
```

Expected: all core unit tests PASS.

- [ ] **Step 9: Commit**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/active_feautes add -u src/core
git -C /tikhome/weeber/es/.claude/worktrees/active_feautes commit -m "core: absorb used_propagations into ActiveFeatures sweep (#5410)

One sweep and one all_reduce fill both masks; one dirty flag
(Propagation::recalc_active_features) serves both. Thermostat changes
now invalidate the cached SYSTEM_DEFAULT expansion.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Migrate `orientation_ghosts_needed()` and the Gay-Berne aggregate

The file-static function moves onto `ActiveFeatures` with exact arms instead of proxies. The two duplicate system-wide Gay-Berne spellings unify into `has_gay_berne()`.

**Files:**
- Modify: `src/core/system/ActiveFeatures.hpp` (two new queries)
- Modify: `src/core/system/ActiveFeatures.cpp` (implementations)
- Modify: `src/core/system/System.cpp:540-672` (delete static function, redirect ghost-flag functions)
- Modify: `src/core/forces.cpp:522-542` (`gay_berne_active`, `dipolar_pair_kernel_active`)
- Modify: `src/core/unit_tests/ActiveFeatures_test.cpp`

**Interfaces:**
- Consumes: Task 1/2 queries; `Propagation::used_propagations`; `Dipoles::Solver` (`system.dipoles.impl->solver`); `InteractionsNonBonded::pair_potential_active(PairPotential)`; `LB::Solver::is_solver_set()`.
- Produces:
  - `bool ActiveFeatures::orientation_ghosts_needed() const` [`ESPRESSO_ROTATION`]
  - `bool ActiveFeatures::has_gay_berne() const` (unconditional; false when no type pair has Gay-Berne configured)

- [ ] **Step 1: Extend the unit test (failing first)**

Add to `src/core/unit_tests/ActiveFeatures_test.cpp`. New includes:

```cpp
#include "cell_system/CellStructure.hpp" // already present
#include "ghosts.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
```

New cases:

```cpp
#ifdef ESPRESSO_ROTATION
BOOST_FIXTURE_TEST_CASE(orientation_ghosts_off_for_default_particles,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  // A default particle cannot rotate, so the ROT_* bits in the SYSTEM_DEFAULT
  // expansion must not trigger orientation ghost exchange.
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  system.propagation->recalc_active_features = true;
  active_features.update_if_needed();
  BOOST_CHECK(not active_features.orientation_ghosts_needed());
  BOOST_CHECK((system.get_global_ghost_flags() & Cells::DATA_PART_QUAT) == 0u);
  BOOST_CHECK_EQUAL(system.get_force_reduce_ghost_flags(), GHOSTTRANS_FORCE);
}

BOOST_FIXTURE_TEST_CASE(orientation_ghosts_on_for_rotating_particles,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  mutate_particle(0, [](Particle &p) { p.set_can_rotate_all_axes(); });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.orientation_ghosts_needed());
  BOOST_CHECK((system.get_global_ghost_flags() & Cells::DATA_PART_QUAT) != 0u);
  BOOST_CHECK((system.get_force_reduce_ghost_flags() & GHOSTTRANS_TORQUE) !=
              0u);
}

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
BOOST_FIXTURE_TEST_CASE(orientation_ghosts_on_for_stoner_wohlfarth,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  // The Stoner-Wohlfarth arm must fire from the particle bit alone,
  // without any rotating particle and without the Langevin thermostat.
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  mutate_particle(0, [](Particle &p) {
    p.stoner_wohlfarth_is_enabled() = true;
  });
  active_features.update_if_needed();
  BOOST_CHECK(active_features.orientation_ghosts_needed());
}
#endif

#ifdef ESPRESSO_DIPOLES
BOOST_FIXTURE_TEST_CASE(dipole_moment_alone_needs_no_orientation_ghosts,
                        ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  // The dipole arm requires a dipolar solver AND a nonzero moment;
  // a moment without a solver must not fire it.
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  mutate_particle(0, [](Particle &p) { p.dipm() = 1.; });
  active_features.update_if_needed();
  BOOST_CHECK(not active_features.orientation_ghosts_needed());
}
#endif

#ifdef ESPRESSO_GAY_BERNE
BOOST_FIXTURE_TEST_CASE(gay_berne_aggregate, ParticleCleanup) {
  auto &system = System::get_system();
  auto &active_features = *system.active_features;
  ::make_new_particle(0, Utils::Vector3d{1., 1., 1.});
  system.propagation->recalc_active_features = true;
  active_features.update_if_needed();
  BOOST_CHECK(not active_features.has_gay_berne());
  system.nonbonded_ias->make_particle_type_exist(0);
  auto &ia_params = system.nonbonded_ias->get_ia_param(0, 0);
  ia_params.gay_berne = GayBerne_Parameters(1., 1., 2., 1., 1., 1., 1.);
  system.on_non_bonded_ia_change();
  BOOST_CHECK(active_features.has_gay_berne());
  BOOST_CHECK(active_features.orientation_ghosts_needed());
  // Reset so later test cases see a clean interaction table.
  ia_params.gay_berne = GayBerne_Parameters();
  system.on_non_bonded_ia_change();
}
#endif
#endif // ESPRESSO_ROTATION
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
make -C build -j8 ActiveFeatures_test
```

Expected: compilation FAILS (`orientation_ghosts_needed` is not a member).

- [ ] **Step 3: Add the queries to `ActiveFeatures`**

In `src/core/system/ActiveFeatures.hpp`, in the public section after `update()`:

```cpp
  /** @brief True when any type pair has the Gay-Berne potential configured.
   *  Live system-derived query; no cached copy. */
  bool has_gay_berne() const;

#ifdef ESPRESSO_ROTATION
  /** @brief True when any active physics requires orientation of ghost
   *  particles. Used by both @c System::get_global_ghost_flags (QUAT push)
   *  and @c System::get_force_reduce_ghost_flags (TORQUE reduce).
   *  Live system-derived query combined with cached particle bits. */
  bool orientation_ghosts_needed() const;
#endif
```

In `src/core/system/ActiveFeatures.cpp`, add includes:

```cpp
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
```

Add the implementations. Move the entire rationale comment block from `src/core/system/System.cpp:540-570` ("Return true when any active physics requires orientation of ghost particles ... harmless (bytes only, value is zero)") onto `orientation_ghosts_needed` and keep it up to date: the Stoner-Wohlfarth paragraph no longer talks about the Langevin proxy but about the particle bit, and add a sentence to the header comment that the rotational-propagation arm additionally requires a particle that can rotate.

```cpp
bool ActiveFeatures::has_gay_berne() const {
  return get_system().nonbonded_ias->pair_potential_active(
      PairPotential::GayBerne);
}

#ifdef ESPRESSO_ROTATION
bool ActiveFeatures::orientation_ghosts_needed() const {
  auto const &system = get_system();

  // Rotational propagation active on a particle that can actually rotate:
  // ghost quat needed for any kernel that uses the director of a ghost
  // particle. The can-rotate check matters because default_propagation
  // carries ROT_* bits for every SYSTEM_DEFAULT particle, while all
  // rotational integrators skip particles with rotation flags 0b000.
  if ((system.propagation->used_propagations &
       (PropagationMode::ROT_EULER | PropagationMode::ROT_LANGEVIN |
        PropagationMode::ROT_BROWNIAN | PropagationMode::ROT_STOKESIAN |
        PropagationMode::ROT_VS_RELATIVE |
        PropagationMode::ROT_VS_INDEPENDENT)) and
      particles_can_rotate()) {
    return true;
  }

  // Virtual sites relative uses p_ref.quat() where p_ref may be a ghost.
  if (system.propagation->used_propagations &
      (PropagationMode::TRANS_VS_RELATIVE | PropagationMode::ROT_VS_RELATIVE |
       PropagationMode::ROT_VS_INDEPENDENT)) {
    return true;
  }

#ifdef ESPRESSO_DIPOLES
  // Dipolar solver set and at least one particle carries a dipole moment:
  // short_range_cabana reads quat -> director for dipole pair kernels on
  // ghosts, and the dipolar solvers read p.calc_dip() on unique_particles.
  if (system.dipoles.impl and system.dipoles.impl->solver and
      particles_have_dipole_moment()) {
    return true;
  }
#endif

#ifdef ESPRESSO_GAY_BERNE
  // Gay-Berne anisotropic nonbonded interaction configured: short_range_cabana
  // commits ghost quat -> director for the Cabana pair kernel.
  if (has_gay_berne()) {
    return true;
  }
#endif

#ifdef ESPRESSO_ENGINE
  // LB active with a swimmer present: the coupling loop includes ghost
  // particles and calls p.calc_director() for swimmer_force_on_fluid
  // particles.
  if (system.lb.is_solver_set() and particles_are_swimmers()) {
    return true;
  }
#endif

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  // Stoner-Wohlfarth integrate_magnetodynamics calls
  // get_reference_particle / get_local_particle which can return a ghost,
  // then reads p_ref->calc_director() on it.
  if (particles_have_stoner_wohlfarth()) {
    return true;
  }
#endif

  return false;
}
#endif // ESPRESSO_ROTATION
```

- [ ] **Step 4: Redirect the consumers in `System.cpp` and `forces.cpp`**

In `src/core/system/System.cpp`:

1. Delete the whole `static bool orientation_ghosts_needed(System const &sys)` function together with its comment block and its `#ifdef ESPRESSO_ROTATION` wrapper (lines ~540-627).
2. In `get_global_ghost_flags()` (line ~656), replace `orientation_ghosts_needed(*this)` with `active_features->orientation_ghosts_needed()`.
3. In `get_force_reduce_ghost_flags()` (line ~667), same replacement.

In `src/core/forces.cpp` (inside `System::System::calculate_forces`, lines ~522-532), replace:

```cpp
  auto const gay_berne_active =
      nonbonded_ias->pair_potential_active(PairPotential::GayBerne);
  auto const dipolar_pair_kernel_active =
#ifdef ESPRESSO_DIPOLES
      get_ptr(dipoles_kernel) != nullptr;
#else
      false;
#endif
```

with:

```cpp
  auto const gay_berne_active = active_features->has_gay_berne();
  // Mirrors the dipole arm of ActiveFeatures::orientation_ghosts_needed():
  // both must tighten together or the GHOSTTRANS_TORQUE assert below fires.
  auto const dipolar_pair_kernel_active =
#ifdef ESPRESSO_DIPOLES
      get_ptr(dipoles_kernel) != nullptr and
      active_features->particles_have_dipole_moment();
#else
      false;
#endif
```

`forces.cpp` needs `#include "system/ActiveFeatures.hpp"` next to its other `system/` includes.

- [ ] **Step 5: Check for other callers of the deleted static function**

```bash
grep -rn "orientation_ghosts_needed" src/
```

Expected: matches only in `ActiveFeatures.hpp`, `ActiveFeatures.cpp`, `System.cpp` (the two redirected calls), the comment in `forces.cpp`, and the unit test.

- [ ] **Step 6: Build and run the tests**

```bash
make -C build -j8 ActiveFeatures_test espresso_core
ctest --test-dir build -R ActiveFeatures --output-on-failure
make -C build -j8 check_unit_tests
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/active_feautes add -u src/core
git -C /tikhome/weeber/es/.claude/worktrees/active_feautes commit -m "core: move orientation_ghosts_needed to ActiveFeatures (#5410)

Exact arms replace proxies: Stoner-Wohlfarth uses the particle bit
instead of THERMO_LANGEVIN, swimmers and dipole moments come from the
particle sweep. The rotation arm now also requires a particle that can
rotate, so plain simulations on ROTATION builds stop exchanging
QUAT/TORQUE ghosts. has_gay_berne() unifies the duplicate spellings.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: Stoner-Wohlfarth guard, sanity check, python runtime-activation test

**Files:**
- Modify: `src/core/magnetostatics/stoner_wohlfarth_thermal.cpp:258` (early return)
- Modify: `src/core/integrate.cpp:324-334` (sanity check)
- Modify: `testsuite/python/thermal_stoner_wohlfarth.py` (new test case)

**Interfaces:**
- Consumes: `ActiveFeatures::particles_have_stoner_wohlfarth()` (fresh at both sites: `integrator_sanity_checks()` runs inside `on_integration_start()` after the integrate-start update; `integrate_magnetodynamics()` runs inside the integration loop, and the collision-detection sites re-update after mid-step topology changes).
- Produces: no new interfaces.

- [ ] **Step 1: Guard `integrate_magnetodynamics()`**

In `src/core/magnetostatics/stoner_wohlfarth_thermal.cpp`, at the top of `System::System::integrate_magnetodynamics()` (line ~258), before `auto const ext_fld = ...`:

```cpp
  if (not active_features->particles_have_stoner_wohlfarth()) {
    return;
  }
```

Add `#include "system/ActiveFeatures.hpp"` to the file's includes. This removes a full per-step particle sweep from every simulation without Stoner-Wohlfarth particles on `THERMAL_STONER_WOHLFARTH` builds.

- [ ] **Step 2: Replace the sanity-check sweep**

In `src/core/integrate.cpp`, replace (lines ~324-334):

```cpp
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  if ((thermo_switch & THERMO_LANGEVIN) == 0) {
    for (auto const &p : cell_structure->local_particles()) {
      if (p.stoner_wohlfarth_is_enabled()) {
        runtimeErrorMsg() << "The thermal Stoner-Wohlfarth model requires the "
                             "Langevin thermostat";
        break;
      }
    }
  }
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH
```

with the O(1) query on the reduced bit, which also makes the error MPI-consistent (today each rank errors independently):

```cpp
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  if (active_features->particles_have_stoner_wohlfarth() and
      (thermo_switch & THERMO_LANGEVIN) == 0) {
    runtimeErrorMsg() << "The thermal Stoner-Wohlfarth model requires the "
                         "Langevin thermostat";
  }
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH
```

- [ ] **Step 3: Add the python runtime-activation test**

In `testsuite/python/thermal_stoner_wohlfarth.py`, add after `test_minimal_no_field` (line ~127):

```python
    def test_runtime_activation(self):
        # Exercise the active-features cache: integrate first with no
        # Stoner-Wohlfarth particle in the system, then add the first one
        # mid-run. The ghost-flag and magnetodynamics guards must pick the
        # change up (issue #5410).
        self.system.part.clear()
        self.system.part.add(pos=[1, 1, 1])
        self.system.integrator.run(10)
        p1, p2 = self._init_virtual_site_pair()
        self.system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(p1.director), np.copy(p2.director), atol=1e-06)
        test_flags = self._check_zero_field_flips(p2, 0.)
        self.assertEqual(all(test_flags), True)
```

- [ ] **Step 4: Build, then run the C++ and python tests**

```bash
make -C build -j8
ctest --test-dir build -R "thermal_stoner_wohlfarth" --output-on-failure
ctest --test-dir build -R ActiveFeatures --output-on-failure
```

Expected: PASS (both SW python tests plus the new case; the fluid variant is skipped if walBerla is off — that is fine).

- [ ] **Step 5: Commit**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/active_feautes add -u src/core testsuite/python/thermal_stoner_wohlfarth.py
git -C /tikhome/weeber/es/.claude/worktrees/active_feautes commit -m "core: guard Stoner-Wohlfarth paths with the particle bit (#5410)

integrate_magnetodynamics() skips its per-step sweep when no particle
uses the model; the sanity check reads the reduced bit and errors
consistently on all ranks.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: Exclusion gate migration, remove the AoSoA aggregate

The specialized-kernel gate reads the globally reduced bit instead of the rank-local AoSoA atomic. Behavior change (spec-approved): one particle with exclusions anywhere disables the specialized kernel on all ranks; the decision becomes globally consistent instead of per-rank.

**Files:**
- Modify: `src/core/forces.cpp:296-304`
- Modify: `src/core/aosoa_pack.hpp:154-183`
- Modify: `src/core/short_range_cabana.hpp` (three sites: the `commit_particle` marking ~line 73-80, the two `reset_any_exclusion()` blocks at ~lines 197-202 and ~299-305)

**Interfaces:**
- Consumes: `ActiveFeatures::particles_have_exclusions()` (fresh at the gate: `calculate_forces()` runs after the integrate-start / `on_observable_calc()` update points).
- Produces: `AoSoA_pack::any_exclusion`, `reset_any_exclusion()`, `has_any_exclusion()`, `mark_any_exclusion()` no longer exist. `set_has_exclusion()`/`has_exclusion()` (per-index flags) remain — the generic kernel still needs them.

- [ ] **Step 1: Redirect the gate**

In `src/core/forces.cpp`, in `create_specialized_verlet_pair_loop` (lines ~296-304), replace:

```cpp
#ifdef ESPRESSO_EXCLUSIONS
  // The specialized kernel has no exclusion handling. The commit sweep (run by
  // update_verlet_state earlier in this same force call) accumulates whether
  // any packed particle carries an exclusion, so this is an O(1) read of the
  // same population the old per-particle sweep covered (local + ghosts).
  if (aosoa.has_any_exclusion())
    return {};
#endif
```

with:

```cpp
#ifdef ESPRESSO_EXCLUSIONS
  // The specialized kernel has no exclusion handling. The bit is globally
  // reduced: one particle with exclusions on any rank disables the
  // specialized kernel on all ranks, so every rank takes the same path.
  if (system.active_features->particles_have_exclusions())
    return {};
#endif
```

The local `auto const &aosoa = ...` above stays — it is still passed to `make_specialized_verlet_pair_loop`.

- [ ] **Step 2: Remove the AoSoA aggregate machinery**

In `src/core/aosoa_pack.hpp`, delete lines ~154-183: the whole comment block starting "Aggregate of the per-particle exclusion flag over the whole pack", the `std::atomic<bool> any_exclusion{false};` member, and the three methods `reset_any_exclusion()`, `has_any_exclusion()`, `mark_any_exclusion()` with their comments. Keep `set_has_exclusion()` and `has_exclusion()`. If `<atomic>` is now unused in the file, remove the include.

In `src/core/short_range_cabana.hpp`:

1. In `commit_particle` (lines ~66-80), reduce the exclusions block to:

```cpp
#ifdef ESPRESSO_EXCLUSIONS
  aosoa.set_has_exclusion(index, not p.exclusions().empty());
#else
  aosoa.flags(index) = 0;
#endif
```

2. Delete both `reset_any_exclusion()` blocks with their comments (the full-rebuild sweep at lines ~197-202 and the partial-update sweep at lines ~299-305), including the now-empty `#ifdef ESPRESSO_EXCLUSIONS` / `#endif` wrappers.

- [ ] **Step 3: Verify nothing else references the aggregate**

```bash
grep -rn "any_exclusion" src/
```

Expected: no matches.

- [ ] **Step 4: Build and run the exclusion tests**

```bash
make -C build -j8
ctest --test-dir build -R "exclusions|auto_exclusions" --output-on-failure
ctest --test-dir build -R ActiveFeatures --output-on-failure
```

Expected: PASS. The exclusions python tests exercise the generic-kernel fallback with the new gate.

- [ ] **Step 5: Commit**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/active_feautes add -u src/core
git -C /tikhome/weeber/es/.claude/worktrees/active_feautes commit -m "core: exclusion gate reads the reduced ActiveFeatures bit (#5410)

The specialized-kernel dispatch becomes globally consistent; the
rank-local AoSoA any-exclusion aggregate is removed.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: Full verification pass

**Files:** none (verification only; fix regressions where they surface).

**Interfaces:** none.

- [ ] **Step 1: Full build**

```bash
make -C build -j8
```

Expected: clean.

- [ ] **Step 2: Full core unit test suite**

```bash
make -C build -j8 check_unit_tests
```

Expected: all PASS.

- [ ] **Step 3: Targeted python regression tests**

List what exists, then run the migrated consumers' coverage — rotation, virtual sites, dipoles, thermostats, exclusions, Stoner-Wohlfarth, integrators:

```bash
ctest --test-dir build -N -R "thermal_stoner_wohlfarth|exclusion|rotat|virtual_sites|dipol|langevin|brownian|integrator|engine|collision"
ctest --test-dir build --output-on-failure -R "thermal_stoner_wohlfarth|exclusion|rotat|virtual_sites|dipol|langevin|brownian|integrator|engine|collision"
```

Expected: all listed tests PASS (walBerla-dependent tests may be skipped). Investigate every failure with the systematic-debugging skill before touching code; a failure here most likely means a stale-cache path (a mutation site that does not set `recalc_active_features`) or an over-tightened ghost-flag arm.

- [ ] **Step 4: Verify the diff is complete and self-consistent**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/active_feautes log --oneline python..HEAD
git -C /tikhome/weeber/es/.claude/worktrees/active_feautes diff python...HEAD --stat
grep -rn "update_used_propagations\|recalc_used_propagations\|any_exclusion" src/
```

Expected: 5 feature commits plus the spec/plan commits; the grep returns nothing.

- [ ] **Step 5: Report**

Summarize: what landed, test evidence (exact ctest output lines), the two behavior changes (globally consistent exclusion gate; QUAT/TORQUE ghost exchange now off for non-rotating systems on ROTATION builds), and remaining follow-ups from the spec. Do not push and do not open a PR; the user decides integration (their workflow pushes to the `me` remote when they say so).
