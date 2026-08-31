# Design: centralized runtime feature state (`System::active_features`)

Implements issue #5410 ("Unify runtime feature config/state").

## Problem

Decision logic for whether a compiled-in feature is actually used in the
current simulation is scattered across the code base. Examples: the
specialized short-range kernel gate (`src/core/forces.cpp:262-326`), the
ghost-flag arms in `orientation_ghosts_needed()`
(`src/core/system/System.cpp:572-626`), and per-step particle sweeps such as
`integrate_magnetodynamics()` (`src/core/integrate.cpp:788`), which today runs
unguarded even with zero Stoner-Wohlfarth particles. The scattered checks
drift out of sync. Three arms of `orientation_ghosts_needed()` are proxies
(Stoner-Wohlfarth via `THERMO_LANGEVIN`, swimmers via "LB is set", dipole
moments via "a dipolar solver is set") because no exact particle-derived
predicate exists.

## Goal

Create one component, `system.active_features`, that answers "is feature X in
use right now?". It distinguishes two kinds of predicates:

* System-derived predicates read live `System` state on every call. They
  store no copies.
* Particle-derived predicates come from one cached bitmask, filled by a
  single sweep over local particles and reduced across MPI ranks with a
  logical OR.

This deviates from the issue's suggestion deliberately: there is no
`update_from_system()`, because system-derived state is queried directly
instead of copied.

## Scope

This is the "component plus focused rollout" scope. The component and its
full particle-derived query interface land now. Only consumers where the
migration fixes something real move now (section "Consumers migrated in
v1"). The remaining catalog migrates in follow-ups.

## Component

New class `ActiveFeatures : public System::Leaf<ActiveFeatures>` in
`src/core/system/ActiveFeatures.hpp` and `.cpp`.

Wiring follows `src/core/system/System.dox`: forward declaration in
`System.hpp`, include in `System.impl.hpp`, member
`std::shared_ptr<ActiveFeatures> active_features`, constructed in the
`System` ctor, bound in `System::initialize()`. No script-interface
counterpart in v1.

## Particle-derived state

### Storage

One `unsigned` bitmask member. Bits, each guarded by its feature ifdef where
the particle field is conditional:

* `HAS_EXCLUSIONS`
* `HAS_STONER_WOHLFARTH`
* `IS_VIRTUAL`
* `CAN_ROTATE`
* `HAS_DIPOLE_MOMENT` (`p.dipm() != 0`)
* `HAS_CHARGE` (`p.q() != 0`)
* `IS_FIXED`
* `HAS_EXT_FORCE`
* `HAS_EXT_TORQUE`
* `IS_SWIMMER`

### Update

One collective method, `update_if_needed()`. When the dirty flag is set, it:

1. Sweeps local particles once via `reduce_over_local_particles`
   (`src/core/particle_reduction.hpp`), ORing the feature bits and
   `p.propagation()` in the same pass.
2. Folds `Propagation::default_propagation` into the propagation mask when
   any particle carries `PropagationMode::SYSTEM_DEFAULT`.
3. Reduces both masks in one `boost::mpi::all_reduce` call over
   `::comm_cart` with a bitwise OR.
4. Stores the propagation mask into `Propagation::used_propagations` and the
   feature mask into the member, then clears the dirty flag.

This absorbs `System::update_used_propagations()`
(`src/core/integrate.cpp:189-201`), which is deleted. It removes one full
particle sweep and one all-reduce per recompute, and it keeps the two masks
in lockstep by construction.

### Queries

Named methods reading the cached mask: `particles_have_exclusions()`,
`particles_have_stoner_wohlfarth()`, `particles_are_virtual()`,
`particles_can_rotate()`, `particles_have_dipole_moment()`,
`particles_have_charge()`, `particles_are_fixed()`,
`particles_have_ext_force()`, `particles_have_ext_torque()`,
`particles_are_swimmers()`. Methods for bits whose feature is not compiled
in are compiled out with the bit.

## Invalidation, update points, staleness

### Dirty flag

One dirty flag serves both masks; one sweep fills both, so a single flag is
correct by construction. It remains on `Propagation` (today's
`recalc_used_propagations`, renamed to reflect the broader meaning), because
`Propagation::set_integ_switch()` must set it without a `System` handle.
`ActiveFeatures::update_if_needed()` reads and clears it. Invalidation sites
are unchanged:
`System::on_particle_change()`, `System::on_particle_local_change()`,
`Propagation::set_integ_switch()`.

### Update points

`update_if_needed()` runs at:

* start of `System::integrate()` (replacing the
  `update_used_propagations()` call at `src/core/integrate.cpp:633-634`;
  note the old call ran unconditionally, the new one honors the dirty
  flag),
* top of `System::on_observable_calc()` (new; see below),
* inside `System::update_dependent_particles()` (replacing the guarded call
  at `src/core/system/System.cpp:381-397`),
* the two collision-detection sites that force a recompute after creating
  virtual particles (`src/core/collision_detection/BindAtPointOfCollision.cpp:191`,
  `GlueToSurface.cpp:207`).

### Staleness semantics

Queries return the state as of the last update. They never trigger the
reduction lazily, because the reduction is collective and query sites are not
guaranteed to be collective. This matches how `used_propagations` is consumed
today: `on_particle_local_change()` reads it stale, and the next collective
point refreshes it.

The new `on_observable_calc()` update point tightens freshness relative to
the status quo: today the first ghost exchange there can see a stale
propagation mask. The update sweeps only local particles, so it is safe to
run before the ghost exchange.

## System-derived queries

Live methods on `ActiveFeatures`, reading sibling state through
`get_system()`. No stored copies. In v1, only queries with migrated
consumers move in:

* `orientation_ghosts_needed()` moves from its file-static home in
  `src/core/system/System.cpp:572-626`. Three proxy arms become exact:
  * Stoner-Wohlfarth: `particles_have_stoner_wohlfarth()` replaces the
    `THERMO_LANGEVIN` proxy.
  * Swimmers: `particles_are_swimmers() and lb.is_solver_set()` replaces
    "LB is set".
  * Dipole moments: "dipolar solver set `and`
    `particles_have_dipole_moment()`" replaces "dipolar solver set".
  The rotation-propagation, virtual-sites-relative, and Gay-Berne arms carry
  over unchanged, as does the rationale comment.
* `has_gay_berne()` unifies the two duplicate system-wide spellings
  (`src/core/forces.cpp:525-526` and `src/core/system/System.cpp:599-603`).
  The name deliberately differs from the per-pair device function
  `gay_berne_configured(IA_parameters const &)` in
  `src/core/short_range_cabana_helpers.hpp`, which stays where the kernels
  need it.

Existing cheap queries stay where they live today; `ActiveFeatures` calls
them instead of duplicating them: `System::has_npt_enabled()`,
`Thermostat::thermo_switch` bits,
`InteractionsNonBonded::combined_active_pair_mask()`.

## Consumers migrated in v1

1. `System::get_global_ghost_flags()` and
   `System::get_force_reduce_ghost_flags()`
   (`src/core/system/System.cpp:634-672`) call
   `active_features->orientation_ghosts_needed()`.
2. `integrate_magnetodynamics()` (`src/core/integrate.cpp:788`) gets a guard:
   skip the sweep when `not particles_have_stoner_wohlfarth()`.
3. The Stoner-Wohlfarth sanity check (`src/core/integrate.cpp:324-333`)
   becomes `particles_have_stoner_wohlfarth() and not (thermo_switch &
   THERMO_LANGEVIN)`. This replaces a per-rank particle sweep and makes the
   error MPI-consistent; today each rank errors independently.
4. The specialized-kernel exclusion gate (`src/core/forces.cpp:302`) reads
   `particles_have_exclusions()` instead of the rank-local AoSoA atomic.
   Behavior change: the gate becomes globally consistent. One particle with
   exclusions anywhere disables the specialized kernel on all ranks; today
   the decision is rank-local (locals plus ghosts). The `any_exclusion`
   machinery in `src/core/aosoa_pack.hpp:164-183` and its marking in
   `src/core/short_range_cabana.hpp` are removed; no consumers remain.

Not migrated in v1:

* The rotation and virtual-sites-COM sanity sweeps
  (`src/core/integrate.cpp:280-320`). They test per-particle conjunctions,
  which global per-property bits cannot express.
* The virtual-particle scan in `src/core/analysis/statistics_chain.cpp:147-156`.
  It runs outside the paths with guaranteed freshness.
* All `used_propagations` mask consumers. Storage and reads are unchanged;
  only the update path moves.

## Error handling

Runtime error reporting is unchanged. The migrated Stoner-Wohlfarth sanity
check gains MPI consistency as described above. Queries perform no MPI
communication, so no new failure modes appear on the query path.

## Testing

* New `src/core/unit_tests/ActiveFeatures_test.cpp`, registered for 1 and 4
  ranks, using the `EspressoCoreGlobalConfig` fixture pattern
  (`src/core/unit_tests/particle_reduction_test.cpp:39-56` is the template).
  Cases:
  * A feature present on only one rank is visible on all ranks after
    `update_if_needed()`.
  * Particle add, remove, and property change set the dirty flag, and the
    next update reflects the change.
  * Bits clear again when the last carrier particle is removed.
  * The absorbed propagation mask matches what
    `update_used_propagations()` produced, including `SYSTEM_DEFAULT`
    expansion.
* Python level: existing rotation, Stoner-Wohlfarth, and exclusion tests
  cover the migrated consumers. Extend one Stoner-Wohlfarth test to add the
  first Stoner-Wohlfarth particle mid-simulation; this exercises runtime
  feature activation through the ghost-flag path.
* Run `maintainer/format` on changed files before any push.

## Follow-ups (out of scope)

* Migrate the remaining catalog: the duplicate NPT spellings
  (`has_npt_enabled()` versus `used_propagations & TRANS_LANGEVIN_NPT`), the
  full specialized-kernel gate composite, the `short_range_only` Verlet
  criterion choice (`src/core/short_range_verlet.cpp:54-71`), the analysis
  scans.
* Script-interface exposure for python-level querying.
* Incremental single-particle OR updates (`update_from_particle`) as an
  optimization; monotone additions could skip the full resweep.
