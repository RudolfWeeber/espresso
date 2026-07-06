# Array-Based Particle Storage in the ESPResSo Core — Design

Date: 2026-07-03
Status: approved (all sections reviewed during brainstorming)

## Goals

- Replace the `Particle` struct (AoS) as the core's particle storage with
  array-based (SoA) storage that is both SIMD- and GPU-friendly, using Kokkos.
- Keep `Particle` objects in the Python interface; the user-facing API does
  not change.
- Clear separation between **state** (position, quaternion, velocity, mass,
  ...) and **observables** (forces, torques, ...): they live in separate
  containers.

## Decisions made during brainstorming

| Topic | Decision |
|---|---|
| Migration strategy | Incremental via proxy: SoA becomes source of truth; a `Particle`-compatible view/proxy keeps all subsystems compiling; hot paths migrate to direct array access subsystem-by-subsystem. |
| Vector layout | Component-major: one `Kokkos::View<double*[3], LayoutLeft>` per quantity (memory: all x, all y, all z — "x/y/z as separate arrays" with a single handle). |
| Residency | Per-field dual residency (DualView-style host+device pair with per-field modified/sync tracking). `Kokkos::ScatterView` for force/torque accumulation. |
| Ragged data (bonds, exclusions) | Mutation-friendly host sidecar containers as source of truth; flattened CSR-style Kokkos views rebuilt for kernels (existing `LocalBondState` pattern). |
| Store shape | Cell-sorted flat store: locals contiguous sorted by cell, ghosts appended; `Cell` = (offset, count) range; resort = stable permutation of all columns. |
| Dependencies | Storage layer is pure Kokkos; Cabana kept only for neighbor lists (`Cabana::VerletList`, `neighbor_parallel_for`). |
| Success bar per checkpoint | Test suite always green. Perf: at most 5% cumulative regression relative to the phase-0 baseline on LJ and P3M benchmarks at: (1) `--particles_per_core 1000` on 1 and 4 MPI ranks; (2) `--particles_per_core 4000` with 4 OMP threads. |

**Amendment (2026-07-05, user-approved):** for the phase-3 and phase-4
checkpoints only, the `--particles_per_core 1000` configurations (all rank
counts) are allowed up to 8% cumulative regression under the sequential gate
protocol (the `--particles_per_core 4000` configurations stay at 5%).
Same-state interleaved A/B measurements place lj-1rank at ~+3%; the gate's
sequential protocol adds cross-run machine drift on this P/E-core host.
Rationale: measured, bounded residual cost of columnar storage on the
ghost-heavy small-per-rank regime (multi-stream ghost packing, store rebuild
machinery) after seven verified optimizations; it amortizes with per-rank
particle count and is scheduled to shrink structurally in phase 4
(commit_particle removal) and phase 5 (columnar PROPERTIES ghost exchange).
The 5% budget for all configurations is re-tightened at the phase-5
checkpoint. The phase-0 baseline was re-recorded 2026-07-05 from a rebuilt
phase-0 binary (machine-state drift ~3-4%; see baselines/machine-info).

**Amendment (2026-07-06, user-approved):** the phase-5 checkpoint re-tightening
did not hold for the `--particles_per_core 1000` LJ configurations. Final gate
(3 reps, quiet machine): lj-1rank 1.113, lj-4rank 1.103, lj-omp 1.044,
p3m-1rank 0.991, p3m-4rank 1.051, p3m-omp 1.027. A perf-recovery round
(pack-contiguous type cache, solver-gated packed charge/dipm, integrator
accessor hoists) recovered p3m-4rank fully (interleaved A/B: 4.4% FASTER than
phase 4; its 1.051 gate reading is cross-run baseline drift) but left an
lj residual that sub-slot wall-timer instrumentation shows is diffuse SoA
multi-stream cost (per-column PROPRTS ghost serialization, integrator
half-step, pair/pack cache streams) with no single recoverable hot spot —
re-confirming the phase-3.5 finding of an inherent ~6-11% columnar penalty at
1000 particles/rank on this host. Amended bar from phase 5 onward: the
`--particles_per_core 1000` configurations are allowed up to 12% cumulative
regression; the `--particles_per_core 4000` configurations stay at 5%. The
1000-ppc budget is re-evaluated at the phase-7 checkpoint, where the Particle
struct, migration carriers, and accessor attached/detached branches are
removed.

## Non-goals

- No change to the Python user interface or its semantics.
- No algorithmic changes to integrators, thermostats, or solvers — this is a
  storage migration; any numerical difference beyond reduction-order effects
  (phases 7–8) is a bug.
- Double precision remains canonical for all stored quantities; solvers that
  internally use single precision keep doing their own casting.
- MPI domain decomposition and the ghost-exchange communication topology are
  unchanged — only the pack/unpack representation changes.

## Section 1: Architecture & data model

**ParticleStore.** A new class (new directory `src/core/particle_store/`) owns
every per-particle quantity in a single index space: local particles first,
contiguously sorted by cell, ghost particles appended after `n_local`.
`Cell` degrades to an `(offset, count)` range into the store;
`Utils::Bag<Particle>` ownership disappears. A consequence of the layout:
`is_ghost(i)` becomes `i >= n_local` — the per-particle ghost flag is deleted.

```
index:   0 ................ n_local-1 | n_local .... n_total-1
order:   [cell 0][cell 1][cell 2]...  | [ghost cell 0][...]

Cell #k  = { offset, count }   // range into store
resort   = stable permutation of all columns
id2index = Kokkos::View<int*> (rebuilt on resort)
```

**Three field categories** (state/observables separation, refined one step):

1. **Parameters** — user-set, constant during integration: id, mol_id, type,
   propagation, rotation/fix bitfields, mass, charge, dipole moment,
   rotational inertia, gamma, gamma_rot, external force/torque, mu_E.
   Cold, feature-gated ones (swimming, Stoner-Wohlfarth, vs_relative) become
   host-only sidecar columns.
2. **State** — evolves during integration: position, image box, quaternion,
   velocity, omega, Lees-Edwards offset/flag, position at last Verlet update,
   RATTLE previous position.
3. **Observables** — recomputed or accumulated every step, never part of
   state proper: force, torque, RATTLE correction, dipole field. Written
   through `Kokkos::ScatterView` in parallel kernels, then contributed into
   the column.

Each hot field is a `FieldColumn<T, N>`: a DualView-style host/device pair of
`Kokkos::View<T*[N], LayoutLeft>` (component-major) plus per-field
modified/sync flags. In CPU-only builds host and device alias the same
memory — zero overhead. Scalars are `View<T*>`. Optional fields stay behind
the existing `ESPRESSO_*` ifdefs. Ragged data (bonds, exclusions) stays in
mutation-friendly host sidecars, flattened to Kokkos views on rebuild.

**One honest API cost:** with component-major layout there is no contiguous
`Vector3d` in memory, so today's reference-returning accessors (`p.pos()`)
cannot survive as-is. The proxy returns a small `VectorRef` supporting read,
`=`, `+=`, and conversion to `Utils::Vector3d` — most call sites compile
unchanged; the rest are mechanical fixes.

**id→index map:** `Kokkos::View<int*>` indexed by particle id (replacing
today's `std::vector<Particle*>` `m_particle_index`), rebuilt on resort;
serves Python access, bond partner resolution, and ghost packing.

## Section 2: Data flow

**Timestep loop** (all kernels are `Kokkos::parallel_for` over index ranges):

1. **Integrate** — velocity-Verlet kernel over `[0, n_local)`: reads
   parameters + observables (force), writes state (pos, v). Propagation-mode
   filtering via the propagation column, as today.
2. **Resort** (when triggered) — the decomposition computes each particle's
   target cell; a stable permutation is built and applied to every column;
   cell `(offset, count)` ranges, the id→index map, ghost rows, and flattened
   bond views are rebuilt. Resort already implies Verlet-rebuild today, so no
   new invalidation cadence is introduced.
3. **Ghost update** — per-field: gather send-rows into `CommBuf` by index,
   MPI exchange, scatter into ghost rows (indices ≥ `n_local`). The
   `GHOSTTRANS_*` flags map onto field groups: POSITION → pos+image+quat,
   MOMENTUM → v+omega, PROPRTS → parameter columns, FORCE → reverse path with
   reduction into local rows. Boost-serialization of `Particle` sub-structs is
   replaced by per-column memcpy-style gathers — bitwise, no archive overhead.
   Bonds-in-ghosts uses the sidecar as today.
4. **Neighbor lists** — Cabana `VerletList` builds directly from the position
   column. The entire `commit_particle` copy-in layer and the `AoSoA_pack`
   mirror are deleted — this is the payoff that funds the ≤5% regression
   budget.
5. **Forces** — zero/init observables columns (external forces, thermostat),
   short-range pair kernel via `Cabana::neighbor_parallel_for`, bonded kernels
   via the flattened bond views, long-range solvers (P3M) read charge/pos
   columns. All contributions go through ScatterView, then one contribute +
   ghost-force reduction.

**Python path** — unchanged UI. `ParticleHandle` resolves id→index on the
owning rank and reads/writes columns through the same MPI machinery as today.
Writes mark the affected field host-modified; the next device-side kernel
syncs only stale fields.

**GPU builds** — kernels launch in the default execution space; per-field
sync tracking means host-only subsystems (Python, ghosts before GPU-aware
MPI) pull only the fields they touch.

## Section 3: Migration phasing

The struct is not replaced in one step; it is **emptied field-by-field**.
Each phase moves one field group's source of truth from `Particle` into
`ParticleStore` columns; the `Particle` accessors for that group redirect
into the columns via a transitional `store_index` member, so all call sites
keep compiling and behavior is identical. When the struct is empty, it is
replaced by the view class. Every phase ends at a mergeable checkpoint:
full test suite green, LJ/P3M benchmarks within the 5% budget.

**Field-eviction order** (simplest ghost-communication direction first):

| Phase | Content | Checkpoint risk retired |
|---|---|---|
| 0 | Record LJ/P3M benchmark baselines for the agreed configs; add a comparison script; C++ unit-test scaffolding for `src/core/particle_store/`. | Perf gate exists before anything moves. |
| 1 | Preparatory refactor, no storage change: centralize every particle insert/remove/move behind a small set of `CellStructure` primitives (today `RegularDecomposition::resort` manipulates `ParticleList`s directly). | Single hook points for mirroring rows; the riskiest transitional machinery gets one home. |
| 2 | `ParticleStore` skeleton: row bookkeeping mirrors the phase-1 primitives; `store_index` member added to `Particle`; store order rebuilt (cell-sorted) at resort; topology changes mark the store dirty for rebuild. Evict **force + torque** (observables): columns become source of truth, `ScatterView` contribution, per-field ghost-force reduction. | Store/row-index consistency machinery proven on the field with the weakest persistence requirement (forces are recomputed after every resort, so resort may simply reallocate the column). |
| 3 | Evict **position, image box, quaternion** (+ `p_old`, `p_last_timestep`, Lees-Edwards offset/flag). Verlet/link-cell and short-range kernels read the position column directly; `commit_particle` position copies deleted; per-field POSITION ghost exchange. | The hot read path and the biggest copy-in cost. First phase where resort must permute state columns. |
| 4 | Evict **velocity + omega**. Thermostats and integrators via proxy accessors on columns; per-field MOMENTUM ghost exchange; remaining `commit_particle`/`AoSoA_pack` mirror deleted. | Dual-storage copy layer fully gone. |
| 5 | Evict **parameters** (id, mol_id, type, propagation, rotation/fix bitfields, mass, charge, dipm, rinertia, gamma, gamma_rot, ext_force/torque, mu_E) and move cold groups (swim, magnetodynamics, vs_relative) to host sidecar columns. Per-field PROPRTS ghost exchange; Python setters write columns. | id→index map replaces `m_particle_index`. |
| 6 | Evict **ragged data**: `BondList` and exclusions become sidecar columns indexed by store row; flattened kernel views unchanged; ghost-bond transfer from sidecar. | `Particle` is now empty except `store_index`. |
| 7 | **Kill the struct**: `Particle` becomes the view class (store pointer + index, same accessor API); `Utils::Bag<Particle>` in cells replaced by (offset, count) ranges; resort rewritten as pure column permutation; inter-rank particle exchange (global resort) switches from boost serialization to per-field pack; `GpuParticleData` replaced by device views of columns. | AoS storage no longer exists anywhere in the core. |
| 8 | **De-proxy hot paths & GPU enablement**: rewrite integrators, thermostats, and remaining per-particle loops as direct column kernels; validate CUDA build; enable device-resident execution with per-field sync. | End-state performance goals. |

**Transitional invariants** (checked in debug builds from phase 2 on):
`store.id(p.store_index()) == p.id()` for every particle after any resort or
topology change; column extents equal `n_local + n_ghost`; a field group is
owned by exactly one side (struct or column) at any commit — accessors of an
evicted field never touch struct memory.

Each phase is a separate sub-project with its own implementation plan
(phases 0 and 1 may share one); phase boundaries are the only merge points.

## Section 4: Error handling

**User-facing error behavior does not change.** Python-level exceptions
(accessing a nonexistent particle, invalid attribute values) and the runtime
error collection machinery (`errorhandling`) keep their exact semantics; bond
resolution/breakage errors already work from Kokkos kernels
(`bond_broken_error`) and are unaffected.

**Internal invariants become debug-build checks** rather than silent
assumptions:

- **Row consistency:** `store.id(p.store_index()) == p.id()` after every
  resort and topology change; column extents equal `n_local + n_ghost`.
  Adjudicated deviation (phase 2): the id-based cross-check
  (`store.id(p.store_index()) == p.id()`) is deferred until the id column
  exists (migration phase 5); phase 2 asserts only the row-index bounds
  (`0 <= row < n_local + n_ghost`) in `ParticleStore::assign_row`.
- **Single ownership:** an evicted field's struct accessor must never touch
  struct memory — enforced by removing the struct member in the same commit
  that evicts the field (compile-time, not runtime).
  Adjudicated deviation (phase 2, observables): migration-carrier members
  ferry force/torque values through whole-particle serialization (inter-rank
  migration, fetch cache); accessors never read them; `assign_row` seeds a
  detached arrival's row from its carriers.
  Adjudicated deviation (phase 3, state fields): for DETACHED particles
  (freshly constructed, deserialized, not yet assigned a row), state-field
  accessors read/write the migration carriers — detached state lives in
  carriers, attach seeds columns from carriers. Single ownership holds
  unconditionally for ATTACHED particles: their accessors touch only
  columns. New invariant this creates: never `assign_row` a particle into a
  foreign store without refreshing its carriers (serialize) or re-writing
  its state through accessors afterwards — stale carriers would seed the row.
  Adjudicated deviation (phase 3, kernels): `commit_particle`'s per-step
  position/image/director copies into the kernel pack were retained at the
  phase-3 flip (kernels still read the pack, not the store columns; the
  pack-index→store-row translation view was not built). Re-adjudication with
  benchmark data is owed at the phase-3 checkpoint; the copy layer is
  scheduled to die with the velocity eviction (phase 4) or a dedicated
  kernel-rewiring step.
- **Sync-state misuse:** in debug builds, reading a field in a memory space
  where it is stale aborts with the field name and the modifying site.
  In release builds sync is automatic (per-field, copy only when stale).
- **id→index lookups:** a missing id yields the sentinel `-1`, mapped to the
  same "particle not here" behavior as today's `nullptr` from
  `get_local_particle()`.
- **Ghost pack/unpack:** per-field byte counts are computed from the
  `GHOSTTRANS_*` flags on both sides; debug builds assert the received
  buffer size matches the expectation (replacing the implicit trust the
  boost-archive layer provided).
- **Allocation:** Kokkos host-side allocation failures propagate as
  exceptions before any column is partially resized (resize order:
  allocate-new, copy, swap).

## Section 5: Testing

- **Primary gate:** the full existing Python test suite (ctest) at every
  phase checkpoint, on 1 and 4 MPI ranks. No test may be weakened or skipped
  to make a checkpoint.
- **New C++ unit tests** for `src/core/particle_store/`: column permutation
  under resort, id→index rebuild, per-field DualView sync semantics,
  `VectorRef` proxy semantics (read/`=`/`+=`/conversion), per-field-group
  ghost pack/unpack roundtrips, and randomized add/remove/resort sequences
  validating the Section 3 invariants.
- **Sanitizers:** ASan/UBSan builds of the unit tests plus a representative
  subset of the Python suite once per phase (row-index bugs manifest as
  out-of-bounds column access — exactly what ASan catches).
- **Benchmarks:** LJ and P3M at `--particles_per_core 1000` (1 and 4 MPI
  ranks) and `--particles_per_core 4000` (4 OMP threads), compared against
  the phase-0 recorded baseline; budget ≤5% cumulative regression.
- **Shared-machine protocol:** the machine is shared. The benchmark
  comparison script checks load (`uptime`, foreign heavy processes) before
  measuring and refuses to record results on a loaded machine; runs are
  repeated several times and the minimum is used. No baseline or gate
  decision is taken from a busy machine.
- **Numerical identity check:** for phases 2–6 (pure storage moves, no
  algorithm change), a short LJ+P3M trajectory must be bitwise-identical
  (positions/forces) before and after the phase on a fixed seed, serial run.
  Phases 7–8 may relax this to statistical agreement where reduction order
  legitimately changes.
