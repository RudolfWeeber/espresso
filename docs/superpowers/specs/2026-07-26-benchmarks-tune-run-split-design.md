# Benchmarks: tune/run separation, state files, and flexible parametrization

Date: 2026-07-26
Scope: `maintainer/benchmarks/{lj,p3m,lb,benchmarks}.py` and a small
ESPResSo script-interface addition (OMP thread count accessor).

## Motivation

The `lj`, `p3m`, and `lb` benchmarks currently do tuning/equilibration and the
timed integration loop in a single process invocation. This makes it impossible
to profile *only* the timed portion without the tuning noise. We also want more
control over particle count, skin, mesh, and skin-retuning from the command
line, and consistent particle-count semantics across MPI and OpenMP
parallelism.

## Goals

1. Separate **tuning** from **running**:
   - `lj.py tune --particles_per_core ... --state_file lj.npz` tunes and saves
     parameters + particle/field state.
   - `lj.py run --state_file lj.npz` runs the timed loop directly from the
     state, with no warmup/tuning (so a profiler sees only the hot loop).
   - `lj.py` with no subcommand does both in one go (current behavior).
2. Flexibility flags:
   - all: optional `--skin`.
   - lj: configurable skin-retuning interval.
   - p3m: explicit `--mesh`, or `--lowest_mesh` / `--highest_mesh`.
3. Consistent particle-count handling: `--n_particles` (fixed total) or
   `--particles_per_core` (counting cores as `mpi_ranks * omp_threads`).
4. State captures parallel topology (OMP threads + MPI node grid); `run`
   refuses on mismatch and restores the exact node grid.

## Non-goals

- No change to `runner.sh` / `suite.sh` / `CMakeLists.txt` behavior. The bare
  (no-subcommand) invocation stays 100% backward compatible, including
  `--output` and `--particles_per_core`, so the existing CMake tests keep
  working unchanged.
- No new benchmark scenarios; `ferrofluid.py` and `mc_acid_base_reservoir.py`
  are out of scope.

## Architecture

### New ESPResSo script-interface API: OMP thread count

Expose the OpenMP thread count in `cell_system.get_state()`, alongside the
existing `n_nodes` / `node_grid` parallel info. In
`CellSystem::do_call_method` (`src/script_interface/cell_system/CellSystem.cpp`,
the `get_state` branch), add `state["omp_num_threads"] =
omp_get_max_threads();` right next to `state["n_nodes"] = ...`, and add
`#include <omp.h>`. Python access is then simply
`system.cell_system.get_state()["omp_num_threads"]` — no Python-side change
needed, and it matches how the benchmarks already read `n_nodes`/`node_grid`.

Rationale: the benchmarks must record and verify the OpenMP thread count
without scraping the `OMP_NUM_THREADS` environment variable (which may be unset
even when threads are active). `benchmarks.py` uses this for `resolve_n_part`
and for state save/verify; `write_report` is updated to use it too for a
consistent thread count in the CSV.

### Shared helpers in `benchmarks.py`

- `get_omp_num_threads(system)` — reads
  `system.cell_system.get_state()["omp_num_threads"]`.
- `n_cores(system)` — `n_nodes * omp_threads`.
- `resolve_n_part(system, args)` — returns `args.n_particles` if set, else
  `args.particles_per_core * n_cores(system)`. `--n_particles` and
  `--particles_per_core` are a mutually-exclusive argparse group;
  `particles_per_core` keeps each script's current default so the bare
  invocation is unchanged.
- `add_common_args(parser, default_ppc)` — registers the optional positional
  subcommand (`tune` | `run`, absent = both), `--skin`, `--state_file`,
  `--output`, and the particle-count group. Returns the subparsers object so
  each script can attach script-specific options.
- `save_state(path, meta, **arrays)` — `np.savez(path, meta=np.array(meta,
  dtype=object), **arrays)`. `meta` is a plain dict of scalars / short lists;
  large data (positions, velocities, charges, LB fields) are passed as named
  arrays and stored in numpy binary form (suitable for large systems).
- `load_state(path)` — returns `(meta, handle)` where
  `meta = np.load(path, allow_pickle=True)["meta"].item()` and `handle` gives
  lazy access to the arrays.
- `verify_topology(system, meta)` — raise (non-zero exit) if the current OMP
  thread count differs from `meta["omp_num_threads"]` or the current MPI rank
  count differs from `meta["n_nodes"]`; otherwise set
  `system.cell_system.node_grid = meta["node_grid"]` to restore the exact grid.

### Mode dispatch (each of lj/p3m/lb)

- **`tune`**: build system → warmup/minimize → equilibrate → tune skin (+ P3M
  for p3m) → `save_state`. Requires `--state_file`. No timing, no CSV.
- **`run`**: `load_state` → `verify_topology` → rebuild the system from state
  (box, interactions, particles, solver) with **no warmup and no tuning** →
  `get_timings` → `write_report`. Requires `--state_file`.
- **bare (both)**: current end-to-end path. If `--state_file` is given, also
  `save_state` after tuning (so the state can be reused by a later `run`).

Common structure: a `build_and_tune(system, args)` function used by `tune` and
`both`, and a `rebuild_from_state(system, meta, handle)` function used by `run`.
`both` calls `build_and_tune` then times directly (no save/reload round-trip
required).

### State file contents (single `.npz`)

Common meta: `skin`, `box_l`, `n_part`, `time_step`, `measurement_steps`,
`n_iterations`, `n_nodes`, `node_grid`, `omp_num_threads`, LJ params
(`epsilon`, `sigma`, `cutoff`), thermostat (`kT`, `gamma`, `seed`).
Common arrays: `pos`, `vel`.

- **lj**: meta adds `bonds` (bool) and, when bonds are used, `harmonic`
  (`r_0`, `k`); array adds `bond_pairs` (list of id pairs).
- **p3m**: meta adds tuned P3M params (`prefactor`, `accuracy`, `mesh`, `cao`,
  `alpha`, `r_cut`, `gpu`) and particle `types`; array adds `charge`. `run`
  re-instantiates `P3M(..., tune=False)` from the saved mesh/cao/alpha/r_cut —
  this is precisely what lets `run` skip P3M tuning.
- **lb**: meta adds LB params (`agrid`, `tau`, `kinematic_viscosity`,
  `density`, `single_precision`, `blocks_per_mpi_rank`, `gpu`, `has_particles`);
  arrays add `lb_velocity` (full field) and `lb_last_applied_force`. `run`
  restores the field via `lbf[:, :, :].velocity = ...` and
  `lbf[:, :, :].last_applied_force = ...`. When `has_particles`, particle
  pos/vel are restored as well.

(Exact ESPResSo accessor names — particle `.v`, LB slice attributes — verified
against the installed API during implementation.)

### Flexibility flags

- **all** — `--skin FLOAT` (optional): when provided, set this skin and
  **disable** skin auto-tuning in `tune`/`both`. Stored in state.
- **lj** — `--retune_skin_after N` (default 5, `0` disables): retune the skin
  every N timing *iterations* inside `get_timings` (current semantics; the
  parameter name is clarified from the misleading `retune_skin_after_steps`).
  Passed through to `get_timings`.
- **p3m** — `--mesh M` (fixed P3M mesh, disables mesh tuning) **XOR**
  `--lowest_mesh` / `--highest_mesh` (override the current
  `tune_limits=[12, 160]`). `--mesh` and the lowest/highest pair are mutually
  exclusive.

### Particle-count semantics

`cores = mpi_ranks * omp_threads`, where `mpi_ranks =
cell_system.get_state()["n_nodes"]` and `omp_threads` comes from the new API.
`--particles_per_core N` → `n_part = N * cores` (default, backward compatible).
`--n_particles N` → `n_part = N` regardless of parallelism. Mutually exclusive.

## Error handling

- `run`/`tune` without `--state_file` → argparse error.
- `run` on a state file whose OMP thread count or MPI rank count differs from
  the current process → explicit error message and non-zero exit
  (`verify_topology`).
- `--n_particles` together with `--particles_per_core` → argparse mutually
  exclusive error.
- `--mesh` together with `--lowest_mesh`/`--highest_mesh` → argparse error.
- Existing assertions (volume fraction bounds, minimum step counts, feature
  checks) are preserved.

## Testing

- Manual smoke test per script on a small system: `tune` writes a state file;
  `run` loads it and produces timings; bare run reproduces current behavior.
- Round-trip check: values in the state file (skin, box, n_part, p3m params)
  match what `tune` computed; `run` restores node_grid and refuses on a
  mismatched rank/thread count.
- Verify the bare invocations used by `CMakeLists.txt` still succeed
  (a couple of the existing lj/p3m/lb argument combinations).
- Build check for the C++/script-interface addition; confirm the new accessor
  returns the expected value under `OMP_NUM_THREADS=1` and `>1`.

## Backward compatibility

Bare invocation is unchanged in behavior and arguments. `benchmarks.py` gains
new functions but keeps `minimize`, `get_timings`, `get_average_time`,
`write_report` signatures working (`get_timings` already accepts
`retune_skin_after_steps`; the lj CLI maps `--retune_skin_after` onto it).
`CMakeLists.txt`, `runner.sh`, `suite.sh` are untouched.
