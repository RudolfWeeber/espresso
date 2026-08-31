# Lees-Edwards Dynamic Ghost Neighborship Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the physical-cell -> ghost-cell wiring for Lees-Edwards (LE) shear a function of the current position offset, rebuilt as the offset evolves, so the short-range loop sees every in-range pair across the shear boundary without `fully_connected_boundary`. The MPI grid keeps the existing restriction `node_grid[shear_dir] == 1` (see Scope revision).

**Scope revision (user decision 2026-08-31):** Splitting the MPI domain along the shear direction is OUT of scope — `node_grid[shear_dir]` must stay 1, the same limit `fully_connected` already enforces. Rationale: the offset's fractional part forces a reach of up to ~1.5 cells across the shear boundary, which a single ghost layer cannot supply when the shear axis is split; keeping that axis unsplit makes the full shear-direction row local, so the single normal-ghost-layer already holds every column. This removes the `fully_connected` *attribute* (the real goal) while keeping its parallelization limit — no regression. Grids that split along the shear-plane normal (`[1,N,1]`) or the neutral axis (`[1,1,N]`), plus OpenMP threading, are all still in scope and proven.

**Architecture:** Keep the LE offset in `LeesEdwardsBC::distance()` (ghost positions carry only periodic-image shifts). All dynamic behavior lives in `RegularDecomposition::init_cell_interactions()` plus a rebuild hook — `make_halo_plan()` is NOT changed. Two changes: (1) `init_cell_interactions()` shifts the shear-direction neighbor window to the sheared columns (`cur + shift ± reach`) at the shear-plane-normal boundary and folds the shear index into the full local row (valid because `node_grid[shear_dir] == 1`); the ghost fill stays geometric, so the single normal-ghost-layer already holds the needed columns. (2) `resort()` rebuilds neighborships when the integer offset shift changes. Correctness is proven by bitwise-identical trajectories (per fixed decomposition, vs the retained `fully_connected` reference), by `non_bonded_loop_trace` pair-set agreement, and by average shear stress — across serial, split-normal, split-neutral, and OpenMP thread counts.

**Tech Stack:** C++20 (ESPResSo core, Kokkos, boost::mpi), Python 3 (espressomd, numpy), CMake, ctest, mpiexec.

**Spec:** `docs/superpowers/specs/2026-08-19-lees-edwards-dynamic-ghost-neighborship-design.md`

## Global Constraints

- Every shell that builds, runs `pypresso`, `ctest`, or `mpiexec` must first `source /tikhome/weeber/es-env/bin/activate` (cmake 4.3.1, Python 3.12, numpy 2.2.6). The build dir is `build/` in the worktree.
- Build with the `maxset` feature config (`maintainer/configs/maxset.hpp`) which enables `DPD` and `LENNARD_JONES`. Build with **CUDA off** for this work. Select the config with `-D ESPRESSO_MYCONFIG_FILE=<abs path to maxset.hpp>`.
- Use `make -j8` (never `-j$(nproc)`).
- Use `git -C <path> ...`, never `cd <path> && git ...`.
- Run `maintainer/format` (clang-format / autopep8) on changed files before every commit, or CI will not start.
- Running python tests: two equivalent routes. (a) Fast iteration against the source tree — `mpiexec -n <N> build/pypresso testsuite/python/<file>.py` — always uses current sources and finds the helper modules (they sit next to the test in source). (b) CI-equivalent — after editing any registered test file or adding one, run `cmake` (auto re-runs on next `make`) and `make python_test_data`, then `ctest -R "<name>$"`. The build copies test files at configure time, so a plain edit-then-ctest can run a stale copy; prefer route (a) while iterating.
- The LE offset stays in `LeesEdwardsBC::distance()`. Ghost positions carry only periodic-image shifts (multiples of the box length), never the LE offset.
- The dynamic path requires `node_grid[shear_dir] == 1` (splitting along the shear direction is out of scope; keep and re-word the existing guard for the dynamic path). Allowed grids: serial `[1,1,1]`, split-normal `[1,N,1]`, split-neutral `[1,1,N]`, and any mix that keeps `node_grid[shear_dir]==1`.
- Bitwise identity vs the old code is defined per fixed decomposition (dynamic vs `fully_connected`, same grid). Across different decompositions only the pair set (trace) and the average shear stress are expected to match.
- Through Phases 1-4 the dynamic path engages only when `m_box.type() == BoxType::LEES_EDWARDS` **and** `fully_connected_boundary()` is unset; when it is set, the exact current behavior is kept as the bitwise reference. Phase 5 then removes `fully_connected_boundary` entirely (user decision 2026-08-30): the setter throws, the C++ implementation is deleted, and the gate simplifies to "active iff the box is Lees-Edwards". Non-LE behavior is unchanged throughout.
- Notation used throughout: `sn = shear_plane_normal`, `sd = shear_direction`, `cs = cell_size[sd]`, `O = pos_offset`, `s = lround(O / cs)` (integer column shift), `frac = O - s*cs` (residual in `[-cs/2, cs/2]`).

---

## Phase 0: Build and baseline

### Task 0.1: Configure and build with maxset, CUDA off

**Files:**
- Create: `build/` (out-of-source build dir in the worktree)

- [ ] **Step 0: Activate the build environment**

All build/test/pypresso commands in this plan run inside the `es-env`
virtualenv (cmake 4.3.1, Python 3.12, numpy 2.2.6). At the start of every
shell:
```bash
source /tikhome/weeber/es-env/bin/activate
```

- [ ] **Step 1: Configure**

The myconfig header is selected by the `ESPRESSO_MYCONFIG_FILE` cmake
variable (the cmake script also honors the `ESPRESSO_MYCONFIG` *environment*
variable; a `-D ESPRESSO_MYCONFIG=` cache variable is NOT read — do not use
it). Run:
```bash
source /tikhome/weeber/es-env/bin/activate
cmake -S /tikhome/weeber/es/.claude/worktrees/comm_le \
      -B /tikhome/weeber/es/.claude/worktrees/comm_le/build \
      -D ESPRESSO_BUILD_WITH_CUDA=OFF \
      -D ESPRESSO_BUILD_WITH_CALIPER=OFF \
      -D CMAKE_BUILD_TYPE=Release \
      -D ESPRESSO_TEST_NP=4 \
      -D ESPRESSO_MYCONFIG_FILE=/tikhome/weeber/es/.claude/worktrees/comm_le/maintainer/configs/maxset.hpp
```
Expected: configuration succeeds. Verify the chosen config:
`grep -i "ESPRESSO_MYCONFIG_FILE" build/CMakeCache.txt` should point at
`maxset.hpp`.

- [ ] **Step 2: Build**

Run: `make -j8 -C /tikhome/weeber/es/.claude/worktrees/comm_le/build`
Expected: `build/pypresso` exists and the build completes without error.
(A cold build takes 15-30 min; do not interrupt it.)

- [ ] **Step 2b: Populate the python test tree**

`make -j8` builds `pypresso` but does NOT copy the test-helper modules
(`tests_common.py`, `unittest_decorators.py`, `thermostats_common.py`,
`data/`) into `build/testsuite/python/`. Without them `ctest` fails with
`ModuleNotFoundError: No module named 'unittest_decorators'`. Build the
data target:
```bash
make -C /tikhome/weeber/es/.claude/worktrees/comm_le/build python_test_data
```
Expected: `build/testsuite/python/unittest_decorators.py` and
`tests_common.py` now exist.

- [ ] **Step 3: Smoke-check pypresso**

Run: `/tikhome/weeber/es/.claude/worktrees/comm_le/build/pypresso -c "import espressomd; espressomd.assert_features(['DPD','LENNARD_JONES']); print('ok')"`
Expected: prints `ok`.

- [ ] **Step 4: Commit** (build dir is git-ignored; nothing to commit — record success in the task log instead.)

### Task 0.2: Confirm the LE testsuite baseline is green

**Files:**
- Test: `testsuite/python/lees_edwards.py` (existing, unmodified)

- [ ] **Step 1: Run the existing LE test on 1 and 4 ranks**

Run:
```bash
cd /tikhome/weeber/es/.claude/worktrees/comm_le/build
ctest -R "lees_edwards$" --output-on-failure
```
Expected: PASS. This confirms the `fully_connected` reference path works before any change.

- [ ] **Step 2: Record baseline** (no commit; note the pass in the task log.)

---

## Phase 1: DPD benchmark

### Task 1.1: Create `maintainer/benchmarks/dpd.py`

**Files:**
- Create: `maintainer/benchmarks/dpd.py`
- Reference: `maintainer/benchmarks/lj.py`, `maintainer/benchmarks/benchmarks.py`, `samples/dpd.py`, `~/dpd/lin_shear.py`

**Interfaces:**
- Consumes: `benchmarks.add_common_args`, `benchmarks.validate_mode`, `benchmarks.resolve_n_part`, `benchmarks.minimize`, `benchmarks.get_timings`, `benchmarks.get_average_time`, `benchmarks.get_omp_num_threads`, `benchmarks.save_state`, `benchmarks.load_state`, `benchmarks.verify_topology`, `benchmarks.topology_meta`, `benchmarks.write_report`, `benchmarks.tune_skin_unless_fixed` (all defined in `benchmarks.py`, read in this session).
- Produces: an executable benchmark script with `tune`/`run`/default modes matching `lj.py`, plus a `--shear_velocity`, `--density`, and `--fully_connected` CLI option. Later verification scripts import nothing from it; it is standalone.

- [ ] **Step 1: Write the benchmark script**

Create `maintainer/benchmarks/dpd.py` with the GPL header (copy the 18-line header verbatim from `lj.py`) followed by:

```python
import sys
import espressomd
import espressomd.lees_edwards
import benchmarks
import numpy as np
import argparse

DPD_GAMMA = 4.5
DPD_R_CUT = 1.5
KT = 1.0
SEED = 42
DEFAULT_PARTICLES_PER_CORE = 1000
INITIAL_SKIN = 0.4
SHEAR_DIRECTION = "x"
SHEAR_PLANE_NORMAL = "y"

# Overridable by the importlib smoke test (see testsuite/benchmarks).
measurement_steps = None
n_iterations = 30
min_skin = 0.2
max_skin = 1.0

parser = argparse.ArgumentParser(
    description="Benchmark a sheared DPD fluid (Lees-Edwards). "
                "Save the results to a CSV file.")
benchmarks.add_common_args(parser, DEFAULT_PARTICLES_PER_CORE)
parser.add_argument("--density", metavar="RHO", action="store", type=float,
                    default=2.0, required=False,
                    help="Number density of DPD particles (default: 2.0)")
parser.add_argument("--shear_velocity", metavar="V", action="store",
                    type=float, default=1.0, required=False,
                    help="Lees-Edwards shear velocity (default: 1.0)")
parser.add_argument("--fully_connected", action="store_true",
                    help="Use the legacy fully_connected_boundary decomposition "
                         "(reference). Default: dynamic sheared halo.")
parser.add_argument("--retune_skin_after", metavar="N", action="store",
                    type=int, default=None, required=False,
                    help="Retune skin every N timing iterations "
                    "(0 disables, default: 5).")
parser.add_argument("--output", metavar="FILEPATH", action="store",
                    type=str, required=False, default="benchmarks.csv",
                    help="Output file (default: benchmarks.csv)")

args = parser.parse_args()
benchmarks.validate_mode(args)

espressomd.assert_features(["DPD"])
np.random.seed(SEED)

system = espressomd.System(box_l=[1, 1, 1])
system.time_step = 0.01


def configure_dpd(system):
    system.non_bonded_inter[0, 0].dpd.set_params(
        weight_function=1, gamma=DPD_GAMMA, r_cut=DPD_R_CUT,
        trans_weight_function=1, trans_gamma=DPD_GAMMA, trans_r_cut=DPD_R_CUT)


def configure_decomposition(system, args):
    """Select the cell decomposition for LE shear."""
    if args.fully_connected:
        normal_axis = {"x": 0, "y": 1, "z": 2}[SHEAR_PLANE_NORMAL]
        n_nodes = system.cell_system.get_state()["n_nodes"]
        system.cell_system.node_grid = [
            n_nodes if i == normal_axis else 1 for i in range(3)]
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True,
            fully_connected_boundary={"boundary": SHEAR_PLANE_NORMAL,
                                      "direction": SHEAR_DIRECTION})
    else:
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True, fully_connected_boundary=None)


def configure_lees_edwards(system, args):
    protocol = espressomd.lees_edwards.LinearShear(
        initial_pos_offset=0.0, shear_velocity=args.shear_velocity)
    system.lees_edwards.set_boundary_conditions(
        shear_direction=SHEAR_DIRECTION,
        shear_plane_normal=SHEAR_PLANE_NORMAL,
        protocol=protocol)


def resolve_retune(args, meta=None):
    if args.retune_skin_after is not None:
        value = args.retune_skin_after
    elif args.skin is not None:
        value = 0
    elif meta is not None:
        value = int(meta["retune_skin_after"])
    else:
        value = 5
    return None if value <= 0 else value


def build_and_tune(system, args):
    n_part = benchmarks.resolve_n_part(system, args)
    steps = measurement_steps
    if steps is None:
        steps = max(20, int(np.round(
            5e6 / n_part * system.cell_system.get_state()["n_nodes"], -2)))
    box_l = (n_part / args.density)**(1. / 3.)
    system.box_l = 3 * (box_l,)
    system.cell_system.skin = INITIAL_SKIN if args.skin is None else args.skin
    configure_dpd(system)
    system.part.add(pos=np.random.random((n_part, 3)) * system.box_l)

    benchmarks.minimize(system, 1 / system.time_step)
    system.integrator.set_vv()
    system.thermostat.set_dpd(kT=KT, seed=SEED)

    print("Equilibration")
    system.time_step /= 10.
    system.integrator.run(100)
    system.time_step *= 10.
    system.integrator.run(min(5 * steps, 60000))
    print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
        system, args, min_skin, max_skin, tol=0.025, int_steps=200)))

    configure_decomposition(system, args)
    configure_lees_edwards(system, args)
    print("Equilibration")
    system.integrator.run(min(10 * steps, 60000))

    return {"n_part": n_part, "measurement_steps": steps}


def save_dpd_state(system, args, ctx):
    resolved_retune = resolve_retune(args)
    p = system.part.all()
    meta = {
        "skin": float(system.cell_system.skin),
        "box_l": [float(x) for x in system.box_l],
        "n_part": int(ctx["n_part"]),
        "time_step": float(system.time_step),
        "measurement_steps": int(ctx["measurement_steps"]),
        "n_iterations": n_iterations,
        "density": float(args.density),
        "shear_velocity": float(args.shear_velocity),
        "fully_connected": bool(args.fully_connected),
        "retune_skin_after": 0 if resolved_retune is None else resolved_retune,
        "kT": KT, "gamma": DPD_GAMMA, "r_cut": DPD_R_CUT, "seed": SEED,
    }
    meta.update(benchmarks.topology_meta(system))
    benchmarks.save_state(
        args.state_file, meta,
        pos=np.copy(p.pos), vel=np.copy(p.v))
    print(f"Saved state to {args.state_file}")


def run_from_state(system, args):
    meta, handle = benchmarks.load_state(args.state_file)
    benchmarks.verify_topology(system, meta)
    system.box_l = meta["box_l"]
    system.time_step = meta["time_step"]
    system.cell_system.skin = args.skin if args.skin is not None else meta["skin"]
    configure_dpd(system)
    system.part.add(pos=handle["pos"], v=handle["vel"])
    system.integrator.set_vv()
    system.thermostat.set_dpd(kT=meta["kT"], seed=int(meta["seed"]))
    args.fully_connected = bool(meta["fully_connected"])
    args.shear_velocity = float(meta["shear_velocity"])
    configure_decomposition(system, args)
    configure_lees_edwards(system, args)
    retune = resolve_retune(args, meta)
    return meta["measurement_steps"], meta["n_iterations"], retune


def time_and_report(system, args, measurement_steps, n_iterations, retune):
    timings = benchmarks.get_timings(
        system, measurement_steps, n_iterations,
        retune_skin_after_steps=retune)
    avg, ci = benchmarks.get_average_time(timings)
    print(f"average: {avg:.3e} +/- {ci:.3e} (95% C.I.)")
    n_proc = system.cell_system.get_state()["n_nodes"]
    benchmarks.write_report(args.output, n_proc, timings, measurement_steps,
                            n_threads=benchmarks.get_omp_num_threads(system))


if args.mode == "run":
    measurement_steps, n_iterations, retune = run_from_state(system, args)
    time_and_report(system, args, measurement_steps, n_iterations, retune)
    sys.exit(0)

ctx = build_and_tune(system, args)

if args.state_file:
    save_dpd_state(system, args, ctx)

if args.mode == "tune":
    sys.exit(0)

time_and_report(system, args, ctx["measurement_steps"], n_iterations,
                resolve_retune(args))
```

- [ ] **Step 2: Format**

Run: `/tikhome/weeber/es/.claude/worktrees/comm_le/maintainer/format/autopep8.sh maintainer/benchmarks/dpd.py` (or the repo's autopep8 wrapper). Expected: no changes needed or file reformatted.

- [ ] **Step 3: Smoke-run serial, dynamic path**

Run:
```bash
/tikhome/weeber/es/.claude/worktrees/comm_le/build/pypresso \
  maintainer/benchmarks/dpd.py --n_particles 200 --output /tmp/dpd_bench.csv
```
Expected: runs to completion, prints an `average: ... +/- ...` line, writes `/tmp/dpd_bench.csv`.

- [ ] **Step 4: Smoke-run the fully_connected reference on 2 ranks**

Run:
```bash
mpiexec -n 2 /tikhome/weeber/es/.claude/worktrees/comm_le/build/pypresso \
  maintainer/benchmarks/dpd.py --n_particles 400 --fully_connected \
  --output /tmp/dpd_bench_fc.csv
```
Expected: runs to completion (node grid set to split along the shear-plane-normal `y`, which `fully_connected` allows).

- [ ] **Step 5: Commit**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/comm_le add maintainer/benchmarks/dpd.py
git -C /tikhome/weeber/es/.claude/worktrees/comm_le commit -m "Add sheared DPD Lees-Edwards benchmark"
```

---

## Phase 2: Verification harness

This phase builds the tooling for the three proofs and validates it against the **existing** `fully_connected` code (so the tooling is trusted before the refactor). Two artifacts: a standalone reproducibility script (trajectory + shear stress, run manually across MPI/OMP) and a testsuite extension (trace pair-set, CI-grade).

### Task 2.1: Standalone reproducibility script

**Files:**
- Create: `maintainer/benchmarks/dpd_le_verify.py`

**Design note (why this shape):** The physically meaningful, decomposition-invariant object is the set of **within-cutoff** interacting pairs on a **given configuration** — NOT the raw `non_bonded_loop_trace` set (which also contains beyond-cutoff *candidate* pairs whose count varies with the ghost/cell structure: e.g. `fully_connected` reports ~472k raw candidates serial vs ~603k on 2 ranks, while the within-cutoff subset is 28247 on both and equals a Python brute-force enumeration exactly). Two independently-integrated runs on different decompositions also **diverge** (DPD trajectories are not bitwise-equal across decompositions), so their late-time pair sets are not comparable. Therefore: (1) compare within-cutoff pair sets only on the **identical step-0 configuration** (deterministic IC ⇒ same config on every decomposition ⇒ this is the cross-decomposition criterion-2 proof); (2) self-check each run's within-cutoff trace against Python brute force (proves that decomposition finds exactly the interacting pairs); (3) compare shear stress **statistically** (z-score over the two time series), never with an arbitrary absolute tolerance.

**Interfaces:**
- Produces: a script invoked as `pypresso dpd_le_verify.py {dump,compare} ...`. `dump` builds a deterministic DPD+LE system with a NONZERO `--initial_pos_offset` (so the shear boundary is offset at step 0), self-checks the within-cutoff trace against Python brute force on the step-0 config (asserts missing=0, extra=0), integrates `--steps` accumulating the shear stress each step, self-checks the trace again on the final config, and writes `DIR/NAME.npz` with: `pos` (final folded positions sorted by id), `pairs0` (sorted within-cutoff pair set on the step-0 config), `stress_mean`, `stress_std`, `stress_n`, `trace_ok` (both self-checks passed). `compare` loads two `.npz`, asserts both `trace_ok`, asserts equal `pairs0` sets (same-config cross-decomposition proof), optionally (`--bitwise`) asserts bitwise-equal final positions (same-decomposition only), and compares shear stress by z-score (`--stress_z`, default 5).

- [ ] **Step 1: Write the script**

Create `maintainer/benchmarks/dpd_le_verify.py` with the GPL header (copy the 18-line header verbatim from `maintainer/benchmarks/lj.py`), then:

```python
import argparse
import itertools
import numpy as np
import espressomd
import espressomd.lees_edwards

SHEAR_DIRECTION = "x"
SHEAR_PLANE_NORMAL = "y"
AXIS = {"x": 0, "y": 1, "z": 2}
DPD_GAMMA = 4.5
DPD_R_CUT = 1.5
KT = 1.0
SEED = 42


def build_system(node_grid, fully_connected, density, box_l, shear_velocity,
                 initial_pos_offset):
    system = espressomd.System(box_l=3 * [box_l])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.non_bonded_inter[0, 0].dpd.set_params(
        weight_function=1, gamma=DPD_GAMMA, r_cut=DPD_R_CUT,
        trans_weight_function=1, trans_gamma=DPD_GAMMA, trans_r_cut=DPD_R_CUT)
    # Deterministic identical initial condition on every decomposition.
    rng = np.random.default_rng(SEED)
    n_part = int(round(density * box_l**3))
    pos = rng.random((n_part, 3)) * box_l
    vel = (rng.random((n_part, 3)) - 0.5)
    system.part.add(pos=pos, v=vel)
    system.integrator.set_vv()
    system.thermostat.set_dpd(kT=KT, seed=SEED)

    system.cell_system.node_grid = list(node_grid)
    if fully_connected:
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True,
            fully_connected_boundary={"boundary": SHEAR_PLANE_NORMAL,
                                      "direction": SHEAR_DIRECTION})
    else:
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True, fully_connected_boundary=None)
    protocol = espressomd.lees_edwards.LinearShear(
        initial_pos_offset=initial_pos_offset, shear_velocity=shear_velocity)
    system.lees_edwards.set_boundary_conditions(
        shear_direction=SHEAR_DIRECTION,
        shear_plane_normal=SHEAR_PLANE_NORMAL, protocol=protocol)
    return system


def trace_selfcheck(system):
    """Within-cutoff pairs the core reports vs Python brute force on the SAME
    folded positions. Returns (ok, within_cut_pairs, n_missing, n_extra)."""
    system.integrator.run(0)
    folded = {p.id: p.pos_folded for p in system.part.all()}
    core = set()
    for id1, id2, _p1, _p2, vec21, _node in \
            system.cell_system.non_bonded_loop_trace():
        if np.linalg.norm(vec21) <= DPD_R_CUT:
            core.add((id1, id2) if id1 < id2 else (id2, id1))
    py = set()
    for i, j in itertools.combinations(folded.keys(), 2):
        if np.linalg.norm(system.distance_vec(folded[i], folded[j])) <= DPD_R_CUT:
            py.add((i, j) if i < j else (j, i))
    missing = py - core
    extra = core - py
    return (len(missing) == 0 and len(extra) == 0), core, len(missing), len(extra)


def shear_stress(system):
    t = system.analysis.pressure_tensor()["total"]
    return float(t[AXIS[SHEAR_DIRECTION]][AXIS[SHEAR_PLANE_NORMAL]])


def dump(args):
    system = build_system(args.node_grid, args.fully_connected, args.density,
                          args.box_l, args.shear_velocity,
                          args.initial_pos_offset)
    # Step-0 within-cutoff self-check on the identical deterministic config:
    # this pair set is the cross-decomposition criterion-2 proof.
    ok0, pairs0, miss0, extra0 = trace_selfcheck(system)
    assert ok0, f"step-0 trace mismatch: missing={miss0} extra={extra0}"
    stresses = []
    for _ in range(args.steps):
        system.integrator.run(1)
        stresses.append(shear_stress(system))
    # Per-run correctness on the (now diverged) final config.
    okN, _pairsN, missN, extraN = trace_selfcheck(system)
    stresses = np.array(stresses)
    p = system.part.all()
    order = np.argsort(p.id)
    np.savez(f"{args.out}/{args.tag}.npz",
             pos=np.copy(p.pos_folded)[order],
             pairs0=np.array(sorted(pairs0), dtype=int),
             stress_mean=float(np.mean(stresses)),
             stress_std=float(np.std(stresses)),
             stress_n=int(len(stresses)),
             trace_ok=bool(ok0 and okN))
    print(f"dumped {args.tag}: step0_within_cut={len(pairs0)} "
          f"trace_ok={ok0 and okN} finalcheck(missing={missN},extra={extraN})")


def compare(args):
    a = np.load(args.a, allow_pickle=True)
    b = np.load(args.b, allow_pickle=True)
    assert bool(a["trace_ok"]) and bool(b["trace_ok"]), \
        "a per-run within-cutoff trace self-check failed"
    print("trace: both runs find exactly the within-cutoff pairs (vs Python)")
    sa = {tuple(r) for r in a["pairs0"]}
    sb = {tuple(r) for r in b["pairs0"]}
    assert sa == sb, \
        f"step-0 within-cutoff pair sets differ: {len(sa ^ sb)} symmetric-diff"
    print(f"trace: identical step-0 within-cutoff pair sets ({len(sa)})")
    if args.bitwise:
        assert np.array_equal(a["pos"], b["pos"]), "positions not bitwise equal"
        print("pos: bitwise identical")
    ma, mb = float(a["stress_mean"]), float(b["stress_mean"])
    sea = float(a["stress_std"]) / max(1, int(a["stress_n"]))**0.5
    seb = float(b["stress_std"]) / max(1, int(b["stress_n"]))**0.5
    denom = (sea**2 + seb**2)**0.5
    z = abs(ma - mb) / denom if denom > 0 else 0.0
    print(f"shear stress: {ma:.6e} vs {mb:.6e}  z={z:.2f} "
          f"(se_a={sea:.2e}, se_b={seb:.2e})")
    assert z < args.stress_z, f"shear stress differs beyond {args.stress_z} sigma"


p = argparse.ArgumentParser()
sub = p.add_subparsers(dest="cmd", required=True)
d = sub.add_parser("dump")
d.add_argument("--tag", required=True)
d.add_argument("--out", default="/tmp")
d.add_argument("--fully_connected", action="store_true")
d.add_argument("--node_grid", type=int, nargs=3, default=[1, 1, 1])
d.add_argument("--steps", type=int, default=200)
d.add_argument("--density", type=float, default=2.0)
d.add_argument("--box_l", type=float, default=10.0)
d.add_argument("--shear_velocity", type=float, default=1.0)
d.add_argument("--initial_pos_offset", type=float, default=1.3)
d.set_defaults(func=dump)
c = sub.add_parser("compare")
c.add_argument("--a", required=True)
c.add_argument("--b", required=True)
c.add_argument("--bitwise", action="store_true")
c.add_argument("--stress_z", type=float, default=5.0)
c.set_defaults(func=compare)
args = p.parse_args()
args.func(args)
```

Note the step-0 self-check uses a nonzero `--initial_pos_offset` so the shear boundary is genuinely offset at step 0. On the pre-refactor dynamic path (`fully_connected=None`) this self-check is EXPECTED to fail (the dynamic halo does not exist yet); the Phase-2 self-consistency check below therefore uses `--fully_connected`, which works today. After Phase 3, dynamic dumps pass their self-check.

- [ ] **Step 2: Format** — run the autopep8 wrapper in `maintainer/format/` on the file.

- [ ] **Step 3: Self-consistency check against the reference (serial vs y-split fully_connected)**

Run:
```bash
B=/tikhome/weeber/es/.claude/worktrees/comm_le/build/pypresso
V=/tikhome/weeber/es/.claude/worktrees/comm_le/maintainer/benchmarks/dpd_le_verify.py
$B $V dump --tag ref1 --fully_connected --node_grid 1 1 1 --out /tmp
mpiexec -n 2 $B $V dump --tag ref2 --fully_connected --node_grid 1 2 1 --out /tmp
$B $V compare --a /tmp/ref1.npz --b /tmp/ref2.npz
```
Expected: "both runs find exactly the within-cutoff pairs", "identical step-0 within-cutoff pair sets (…)", and a shear-stress z-score below 5. This proves the harness works on known-good code. (`node_grid 1 2 1` splits along the shear-plane normal `y`, which `fully_connected` supports; do NOT use a grid that splits along the shear direction `x` here — `fully_connected` throws on that.)

- [ ] **Step 4: Commit**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/comm_le add maintainer/benchmarks/dpd_le_verify.py
git -C /tikhome/weeber/es/.claude/worktrees/comm_le commit -m "Add DPD Lees-Edwards decomposition verification script"
```

### Task 2.2: Add a DPD pair-visibility test to the LE testsuite

**Files:**
- Modify: `testsuite/python/lees_edwards.py` (add a method mirroring `run_lj_pair_visibility`, near line 937)

**Interfaces:**
- Consumes: `tests_common.check_non_bonded_loop_trace(ut_obj, system, cutoff)` (read this session; asserts every within-cutoff pair is found exactly once and vec21 matches the minimum-image distance).
- Produces: `LeesEdwards.run_dpd_pair_visibility(self, shear_direction, shear_plane_normal)`, used by a later test task.

- [ ] **Step 1: Add the method**

Insert after `run_lj_pair_visibility` (before `test_zz_lj_pair_visibility`):

```python
    def run_dpd_pair_visibility(self, shear_direction, shear_plane_normal):
        """
        DPD particles under linear shear: verify the short-range loop finds
        every within-cutoff pair across the shear boundary for a range of
        offsets, using the pair trace.
        """
        assert espressomd.has_features(["DPD"])
        shear_axis, normal_axis = axis(
            shear_direction), axis(shear_plane_normal)
        system = self.system
        system.part.clear()
        system.time = 0
        system.time_step = 0.05
        cutoff = 1.5
        protocol = espressomd.lees_edwards.LinearShear(
            shear_velocity=2., initial_pos_offset=0.)
        system.lees_edwards.set_boundary_conditions(
            shear_direction=shear_direction,
            shear_plane_normal=shear_plane_normal, protocol=protocol)
        system.cell_system.skin = 0.2
        system.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=1, gamma=4.5, r_cut=cutoff,
            trans_weight_function=1, trans_gamma=4.5, trans_r_cut=cutoff)
        rng = np.random.default_rng(42)
        system.part.add(pos=rng.random((50, 3)) * system.box_l,
                        v=(rng.random((50, 3)) - 0.5))
        system.thermostat.set_dpd(kT=1., seed=42)
        tests_common.check_non_bonded_loop_trace(
            self, system, cutoff=cutoff + system.cell_system.skin)
        for _ in range(50):
            system.integrator.run(3)
            tests_common.check_non_bonded_loop_trace(
                self, system, cutoff=cutoff + system.cell_system.skin)
```

- [ ] **Step 2: Run the method under the existing (fully_connected) reference to confirm it's a valid check**

Add a temporary test that exercises it with `fully_connected` on `y`-split, run it, confirm PASS, then remove the temporary test (the permanent tests are added in Task 4.2). Concretely:

Temporarily append inside the class:
```python
    @utx.skipIfMissingFeatures(["DPD"])
    def test_tmp_dpd_reference(self):
        system = self.system
        system.box_l = [10, 10, 10]
        system.cell_system.node_grid = [1, self.n_nodes, 1]
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True,
            fully_connected_boundary={"boundary": "y", "direction": "x"})
        self.run_dpd_pair_visibility("x", "y")
```
Run:
```bash
cd /tikhome/weeber/es/.claude/worktrees/comm_le/build && make -j8 && \
  ctest -R "lees_edwards$" --output-on-failure
```
Expected: PASS (the DPD trace check is satisfied by the reference). Then delete `test_tmp_dpd_reference`.

- [ ] **Step 3: Format** — run autopep8 on `testsuite/python/lees_edwards.py`.

- [ ] **Step 4: Commit**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/comm_le add testsuite/python/lees_edwards.py
git -C /tikhome/weeber/es/.claude/worktrees/comm_le commit -m "Add DPD pair-visibility helper to Lees-Edwards testsuite"
```

---

## Phase 3: Core refactor

All C++ changes live in `src/core/cell_system/RegularDecomposition.{hpp,cpp}`. The refactor is staged so each task has an independently runnable trace test: Task 3.2 covers `node_grid[sd]==1` (local self-copy halo), Task 3.4 covers splitting along `sd` (cross-rank halo). Task 3.3 keeps the stencil valid as the offset grows.

### Task 3.1: Add an LE-shear geometry helper

**Files:**
- Modify: `src/core/cell_system/RegularDecomposition.hpp` (add a private helper declaration near line 213)
- Modify: `src/core/cell_system/RegularDecomposition.cpp` (add the helper definition; add `#include <cmath>` and `#include "lees_edwards/LeesEdwardsBC.hpp"` if not already present)

**Interfaces:**
- Produces: a private member
  ```cpp
  struct LeShear { bool active; unsigned sn; unsigned sd; int shift; double frac; };
  LeShear le_shear() const;
  ```
  `active` is true iff `m_box.type() == BoxType::LEES_EDWARDS` and `fully_connected_boundary()` is unset. `sn`/`sd` are the shear-plane-normal / shear-direction coordinates. `shift = lround(pos_offset / cell_size[sd])`. `frac = pos_offset - shift * cell_size[sd]`. When inactive, `active=false` and the other fields are 0.

- [ ] **Step 1: Write a C++ unit test**

Add to `src/core/unit_tests/` a small test (or extend the existing `RegularDecomposition`-related unit test if present; search `src/core/unit_tests/CMakeLists.txt` for an existing target). If no suitable target exists, defer the assertion to the Python trace tests and skip this step — but first check:
Run: `grep -rn "RegularDecomposition\|le_shear\|LeShear" src/core/unit_tests/` and reuse the nearest cell-system unit test file. The test asserts, for a box with `cell_size[sd]=1.0` and `pos_offset=2.6`, that `le_shear().shift == 3` and `le_shear().frac` is close to `-0.4`.

- [ ] **Step 2: Run it, expect fail** (undefined `le_shear`).

- [ ] **Step 3: Implement the helper**

In `RegularDecomposition.hpp`, add to the private section:
```cpp
  struct LeShear {
    bool active = false;
    unsigned sn = 0u;
    unsigned sd = 0u;
    int shift = 0;
    double frac = 0.;
  };
  LeShear le_shear() const;
```
In `RegularDecomposition.cpp`:
```cpp
RegularDecomposition::LeShear RegularDecomposition::le_shear() const {
  LeShear r;
  if (m_box.type() != BoxType::LEES_EDWARDS or fully_connected_boundary())
    return r;
  auto const &le = m_box.lees_edwards_bc();
  r.active = true;
  r.sn = le.shear_plane_normal;
  r.sd = le.shear_direction;
  auto const cs = cell_size[r.sd];
  r.shift = static_cast<int>(std::lround(le.pos_offset / cs));
  r.frac = le.pos_offset - static_cast<double>(r.shift) * cs;
  return r;
}
```

- [ ] **Step 4: Run the unit test, expect pass** (or, if deferred, build only: `make -j8 -C build`).

- [ ] **Step 5: Commit**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/comm_le add src/core/cell_system/RegularDecomposition.hpp src/core/cell_system/RegularDecomposition.cpp
git -C /tikhome/weeber/es/.claude/worktrees/comm_le commit -m "Add Lees-Edwards shear geometry helper to RegularDecomposition"
```

### Task 3.2: Shifted sheared stencil with shear-axis fold

**Design (read before coding — this supersedes any earlier "shift-in-source /
make_halo_plan" approach):** With `node_grid[shear_dir] == 1` enforced, the
shear direction (`sd`) is fully local to every rank, so the single
shear-plane-normal (`sn`) ghost layer already contains EVERY `sd` column
(filled geometrically). Therefore all dynamic behavior lives in
`init_cell_interactions()` and `make_halo_plan()` is NOT touched:

1. At the `sn` boundary, SHIFT the `sd` neighbor window to the sheared
   columns: `[cur + shift - reach, cur + shift + reach]` (this is the
   `fully_connected` trick pruned to the sheared window). The ghost fill
   stays geometric; `distance()` applies the offset as today.
2. Fold the `sd` neighbor index into the local `sd` row `[0, cell_grid[sd])`
   for every LE-active neighbor (both the `one_mpi_rank` path and the
   multi-rank path). Because `node_grid[sd]==1`, the global `sd` range equals
   the local range, so the fold maps the (possibly far, shifted) column to a
   valid local/ghost cell — no second ghost layer, no out-of-bounds.

This works uniformly for serial `[1,1,1]`, split-normal `[1,N,1]`, and
split-neutral `[1,1,N]`. Bitwise identity vs `fully_connected` on a fixed
grid is preserved: same minimal cells, same geometric fill, and the sheared
in-range columns get the same `folded_linear_index` order (the columns
`fully_connected` additionally wires hold no in-range partner).

**Sign:** the sign of `shift` in the window is determined empirically by the
trace tests below — two choices only. Start with `+ le.shift`; if a test
reports missed pairs, flip to `- le.shift` and rebuild. The tests are the
arbiter.

**Note:** the current HEAD already contains an earlier geometric-stencil +
serial-source-shift attempt (commit 7f328a8e44) that passes serial but
crashes multi-rank. REPLACE that logic with the shifted-window + fold below;
remove the `one_mpi_rank`-only source shift it added.

**Files:**
- Modify: `src/core/cell_system/RegularDecomposition.cpp` — `init_cell_interactions()` (the LE window branch and the neighbor-fold). Do NOT modify `make_halo_plan()`.

**Interfaces:**
- Consumes: `le_shear()` (Task 3.1); `max_range()` for the reach; existing `at_boundary`, `folded_linear_index`, the neighbor loop, red/black machinery.
- Produces: for the LE-active case, a shifted `±reach` `sd` window at the `sn` boundary and a `sd`-into-local-row fold. The fully-connected branch is left untouched (reference path).

- [ ] **Step 1: Write the failing static trace tests (serial + split-normal)**

Both fix the offset (shear_velocity=0), set it BEFORE building the
decomposition, and do no integration — isolating the stencil. Append to
`testsuite/python/lees_edwards.py`:
```python
    @utx.skipIfMissingFeatures(["DPD"])
    def _dpd_static_visibility(self, node_grid):
        system = self.system
        system.part.clear()
        system.box_l = [10, 10, 10]
        cutoff = 1.5
        system.cell_system.skin = 0.2
        protocol = espressomd.lees_edwards.LinearShear(
            shear_velocity=0., initial_pos_offset=3.7)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=protocol)
        system.cell_system.node_grid = node_grid
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True, fully_connected_boundary=None)
        system.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=1, gamma=4.5, r_cut=cutoff,
            trans_weight_function=1, trans_gamma=4.5, trans_r_cut=cutoff)
        rng = np.random.default_rng(7)
        system.part.add(pos=rng.random((60, 3)) * system.box_l,
                        v=(rng.random((60, 3)) - 0.5))
        system.thermostat.set_dpd(kT=1., seed=7)
        system.integrator.run(0)
        tests_common.check_non_bonded_loop_trace(
            self, system, cutoff=cutoff + system.cell_system.skin)

    @utx.skipIfMissingFeatures(["DPD"])
    def test_dpd_dynamic_serial_static(self):
        """Serial dynamic sheared halo at a fixed nonzero offset."""
        if self.n_nodes > 1:
            self.skipTest("serial (node_grid [1,1,1]) test")
        self._dpd_static_visibility([1, 1, 1])

    @utx.skipIfMissingFeatures(["DPD"])
    def test_dpd_dynamic_split_normal_static(self):
        """Dynamic sheared halo split along the shear-plane normal, fixed
        offset. node_grid[sd]==1 keeps the shear axis local."""
        if self.n_nodes < 2:
            self.skipTest("needs >= 2 MPI ranks")
        self._dpd_static_visibility([1, self.n_nodes, 1])
```

- [ ] **Step 2: Run them, expect fail**

Run:
```bash
source /tikhome/weeber/es-env/bin/activate && cd build && make -j8
mpiexec -n 1 ./pypresso ../testsuite/python/lees_edwards.py LeesEdwards.test_dpd_dynamic_serial_static
mpiexec -n 2 ./pypresso ../testsuite/python/lees_edwards.py LeesEdwards.test_dpd_dynamic_split_normal_static
```
Expected: FAIL ("Pair not found by the core") or the current out-of-bounds crash on the multi-rank case — both are the pre-implementation state.

- [ ] **Step 3: Implement the shifted window + shear-axis fold**

In `init_cell_interactions()`, after `global_size`, add `auto const le = le_shear();` and:
```cpp
  auto const le_reach =
      le.active
          ? static_cast<int>(std::ceil(max_range()[le.sd] / cell_size[le.sd])) + 1
          : 0;
```
Replace any earlier LE window branch with the SHIFTED window at the `sn`
boundary (keep the fully-connected branch untouched):
```cpp
        if (le.active and at_boundary(le.sn, {m, n, o})) {
          auto const cur = Utils::Vector3i{m, n, o}[le.sd];
          lower_index[le.sd] = cur + le.shift - le_reach;  // flip sign if test fails
          upper_index[le.sd] = cur + le.shift + le_reach;
        }
```
Remove the `one_mpi_rank`-only source shift added by the prior attempt. In
the neighbor loop, after the existing `one_mpi_rank` per-coordinate fold
block, add an UNCONDITIONAL (for LE-active) fold of `sd` into the local row —
it must run in BOTH the serial and multi-rank paths:
```cpp
          /* node_grid[sd]==1 is enforced, so the shear direction is fully
           * local: fold the (shifted) sd neighbor into [0, cell_grid[sd]) and
           * connect to that local/ghost cell directly.  No second ghost layer
           * is needed. */
          if (le.active) {
            neighbor[le.sd] =
                ((neighbor[le.sd] % cell_grid[le.sd]) + cell_grid[le.sd]) %
                cell_grid[le.sd];
          }
```
Do NOT apply `fcb_is_inner_connection` in the LE branch. Do NOT modify
`make_halo_plan()`.

- [ ] **Step 4: Run the tests, expect pass (flip the sign if needed)**

Run:
```bash
cd build && make -j8
mpiexec -n 1 ./pypresso ../testsuite/python/lees_edwards.py LeesEdwards.test_dpd_dynamic_serial_static
mpiexec -n 2 ./pypresso ../testsuite/python/lees_edwards.py LeesEdwards.test_dpd_dynamic_split_normal_static
mpiexec -n 3 ./pypresso ../testsuite/python/lees_edwards.py LeesEdwards.test_dpd_dynamic_split_normal_static
```
Expected: PASS on 1/2/3 ranks. If missed pairs, flip the window sign to
`cur - le.shift ± le_reach` and rebuild; one sign passes all.

- [ ] **Step 5: Regression + format**

Run `mpiexec -n 1 ./pypresso ../testsuite/python/lees_edwards.py` and confirm all tests pass (the fully_connected and non-LE tests are untouched). Then run the clang-format wrapper on `src/core/cell_system/RegularDecomposition.cpp`.

- [ ] **Step 6: Commit**

```bash
git -C /tikhome/weeber/es/.claude/worktrees/comm_le add src/core/cell_system/RegularDecomposition.cpp testsuite/python/lees_edwards.py
git -C /tikhome/weeber/es/.claude/worktrees/comm_le commit -m "Shifted sheared stencil with shear-axis fold for Lees-Edwards"
```

### Task 3.3: Rebuild neighborships and halo plan on resort

**Files:**
- Modify: `src/core/cell_system/RegularDecomposition.hpp` (add member `int m_le_shift_at_last_build = 0;` in the private section)
- Modify: `src/core/cell_system/RegularDecomposition.cpp` — `resort()` (lines 187-...); factor the constructor's build sequence (lines 829-935: `init_cell_interactions(); mark_cells(); make_halo_plan(); ...boundary marking...`) into a private `rebuild_topology()` method and call it from both the constructor and `resort()`.

**Interfaces:**
- Consumes: `le_shear()` (3.1).
- Produces: a private `void rebuild_topology();` that runs `init_cell_interactions()`, `mark_cells()`, `m_halo_plan = make_halo_plan()`, and the boundary-marking block. `resort()` calls it when `le_shear().active && le_shear().shift != m_le_shift_at_last_build`. `m_le_shift_at_last_build` is updated inside `rebuild_topology()`.

- [ ] **Step 1: Write the failing serial dynamic (offset-growth) test**

The 3.2 static test fixed the offset. This test integrates under shear on 1
rank so the offset GROWS through integer-cell crossings, which requires a
topology rebuild. Append to `testsuite/python/lees_edwards.py`:
```python
    @utx.skipIfMissingFeatures(["DPD"])
    def test_dpd_dynamic_serial(self):
        """Serial dynamic sheared halo under integration (offset grows past
        cell boundaries): exercises rebuild-on-resort."""
        if self.n_nodes > 1:
            self.skipTest("serial (node_grid [1,1,1]) test")
        system = self.system
        system.box_l = [10, 10, 10]
        system.cell_system.node_grid = [1, 1, 1]
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True, fully_connected_boundary=None)
        self.run_dpd_pair_visibility("x", "y")

    @utx.skipIfMissingFeatures(["DPD"])
    def test_dpd_dynamic_split_normal(self):
        """Dynamic sheared halo split along the shear-plane normal, under
        integration (offset grows): exercises rebuild-on-resort on >1 rank."""
        if self.n_nodes < 2:
            self.skipTest("needs >= 2 MPI ranks")
        system = self.system
        system.box_l = [10, 10, 10]
        system.cell_system.node_grid = [1, self.n_nodes, 1]
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True, fully_connected_boundary=None)
        self.run_dpd_pair_visibility("x", "y")
```

Run:
```bash
source /tikhome/weeber/es-env/bin/activate && cd build && make -j8
mpiexec -n 1 ./pypresso ../testsuite/python/lees_edwards.py LeesEdwards.test_dpd_dynamic_serial
mpiexec -n 2 ./pypresso ../testsuite/python/lees_edwards.py LeesEdwards.test_dpd_dynamic_split_normal
```
Expected: FAIL after several integrator runs — the stencil built at the
current offset goes stale once the offset drifts past a cell boundary and a
resort fires without a topology rebuild.

- [ ] **Step 2: Extract `rebuild_topology()`**

Move the constructor body from `init_cell_interactions();` (line 829) through the end of the boundary-marking block (line ~935, before the closing brace of the constructor) into:
```cpp
void RegularDecomposition::rebuild_topology() {
  init_cell_interactions();
  mark_cells();
  m_halo_plan = make_halo_plan();
  // ... (the existing boundary-marking block verbatim) ...
  m_le_shift_at_last_build = le_shear().shift;
}
```
and replace the constructor block with a single `rebuild_topology();` call (after `create_cell_grid(range);`). Declare `void rebuild_topology();` and `int m_le_shift_at_last_build = 0;` in the header.

- [ ] **Step 3: Hook it into `resort()`**

At the end of `RegularDecomposition::resort()` (after the particle-move sweep, before returning), add:
```cpp
  auto const le = le_shear();
  if (le.active and le.shift != m_le_shift_at_last_build) {
    rebuild_topology();
  }
```
Note: `resort()` runs after particles are folded and moved to their target cells; `rebuild_topology()` depends only on the cell grid and the LE offset, not on particle positions, so the ordering is safe. `CellStructure::resort_particles` already sets `m_rebuild_verlet_list*` so the Verlet list is rebuilt from the new neighborships, and every `ghosts_update` re-reads `decomposition().halo_plan()`, so the fresh plan is picked up automatically.

- [ ] **Step 4: Run the tests, expect pass**

Run:
```bash
cd build && make -j8
mpiexec -n 1 ./pypresso ../testsuite/python/lees_edwards.py LeesEdwards.test_dpd_dynamic_serial
mpiexec -n 2 ./pypresso ../testsuite/python/lees_edwards.py LeesEdwards.test_dpd_dynamic_split_normal
mpiexec -n 3 ./pypresso ../testsuite/python/lees_edwards.py LeesEdwards.test_dpd_dynamic_split_normal
mpiexec -n 1 ./pypresso ../testsuite/python/lees_edwards.py   # full regression
```
Expected: PASS across all 50 integrator runs at 1/2/3 ranks (offset now triggers a topology rebuild whenever the integer column shift changes; between rebuilds the `cutoff+skin`-sized window and the resort trigger's `skin/2` drift bound keep every pair visible).

- [ ] **Step 5: Format and commit**

```bash
maintainer/format/clang-format.sh src/core/cell_system/RegularDecomposition.cpp
git -C ... add src/core/cell_system/RegularDecomposition.hpp src/core/cell_system/RegularDecomposition.cpp testsuite/python/lees_edwards.py
git -C ... commit -m "Rebuild Lees-Edwards topology on integer-offset crossing"
```

### Task 3.4: Enforce the shear-axis node-grid guard + full verification

The dynamic behavior is complete after Tasks 3.2-3.3 (shifted stencil +
shear-axis fold + rebuild-on-resort); `make_halo_plan` is deliberately
unchanged. This task adds the safety guard that the dynamic path requires
`node_grid[shear_dir] == 1`, and verifies the full allowed decomposition
matrix.

**Where the guard goes — read this first.** Task 3.2 made `le_shear()` return
`active=false` when `node_grid[shear_dir] != 1` (to avoid an out-of-bounds
crash during transient decomposition rebuilds, e.g. `test_zz` reconfigures
`node_grid` while a previous iteration's LE still targets the now-split
axis). Consequently a throw inside `init_cell_interactions()` gated on
`le.active` would NEVER fire, and a throw NOT gated on `le.active` would fire
on those harmless transient rebuilds. So the loud guard must live at
**integration time**, not decomposition-build time: it catches the genuine
"user runs with `node_grid[shear]>1` + LE" config while ignoring transient
rebuilds that never integrate.

**Files:**
- Modify: `src/core/integrate.cpp` — `System::integrate()` prologue (near `on_integration_start()`, ~line 626): throw when the box is Lees-Edwards and `node_grid[shear_direction] != 1`.
- Modify: `testsuite/python/lees_edwards.py` — add the guard test and the split-neutral test.

**Interfaces:**
- Consumes: `box_geo->lees_edwards_bc().shear_direction`; `::communicator.node_grid`.
- Produces: a clear `std::runtime_error` at the start of integration when a Lees-Edwards run has `node_grid[shear_dir] != 1`; a green full matrix (serial, split-normal, split-neutral) at 1/2/3/4 ranks.

- [ ] **Step 1: Write the guard + split-neutral tests**

Append to `testsuite/python/lees_edwards.py`:
```python
    @utx.skipIfMissingFeatures(["DPD"])
    def test_dpd_dynamic_split_neutral(self):
        """Dynamic sheared halo split along the neutral axis (z)."""
        if self.n_nodes < 2:
            self.skipTest("needs >= 2 MPI ranks")
        system = self.system
        system.box_l = [10, 10, 10]
        system.cell_system.node_grid = [1, 1, self.n_nodes]  # split along z
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True, fully_connected_boundary=None)
        self.run_dpd_pair_visibility("x", "y")

    @utx.skipIfMissingFeatures(["DPD"])
    def test_dpd_split_shear_rejected(self):
        """Splitting along the shear direction is not supported by the
        dynamic path and must raise a clear error."""
        if self.n_nodes < 2:
            self.skipTest("needs >= 2 MPI ranks")
        system = self.system
        system.box_l = [10, 10, 10]
        system.cell_system.node_grid = [self.n_nodes, 1, 1]  # split along x = sd
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True, fully_connected_boundary=None)
        protocol = espressomd.lees_edwards.LinearShear(
            shear_velocity=1., initial_pos_offset=0.)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=protocol)
        with self.assertRaises(Exception):
            system.integrator.run(0)  # run-time guard rejects the bad config
```

- [ ] **Step 2: Run them, observe current state**

Run:
```bash
source /tikhome/weeber/es-env/bin/activate && cd build && make -j8
mpiexec -n 2 ./pypresso ../testsuite/python/lees_edwards.py LeesEdwards.test_dpd_dynamic_split_neutral LeesEdwards.test_dpd_split_shear_rejected
```
Expected: `split_neutral` PASSES already (z-split does not touch the shear
axis). `split_shear_rejected` FAILS (no guard yet — `run(0)` completes
silently instead of raising, because `le_shear()` fell back to the standard
stencil).

- [ ] **Step 3: Add the run-time guard**

In `src/core/integrate.cpp`, in `System::integrate()` prologue (just after
`on_integration_start();`, ~line 626), add:
```cpp
  if (box_geo->type() == BoxType::LEES_EDWARDS) {
    auto const sd =
        static_cast<unsigned>(box_geo->lees_edwards_bc().shear_direction);
    if (sd < 3u and ::communicator.node_grid[sd] != 1) {
      throw std::runtime_error(
          "Lees-Edwards requires the MPI node grid to be 1 along the shear "
          "direction; splitting the domain along the shear direction is not "
          "supported.");
    }
  }
```
`node_grid` and the box type are identical on every rank, so all ranks
evaluate the condition the same way and throw together — no collective
deadlock. Add includes if needed (`communicator.hpp` for `::communicator`,
`BoxGeometry.hpp` for `BoxType` — both are likely already visible). This
guard is intentionally general (it also covers the retained `fully_connected`
path, which already requires `node_grid[shear]==1`), so Phase 5 can delete
`fully_connected`'s own build-time throw and rely on this one.

- [ ] **Step 4: Run the full allowed matrix, expect pass**

Run:
```bash
cd build && make -j8 && make python_test_data
for n in 1 2 3 4; do
  echo "=== $n ranks ==="
  mpiexec --oversubscribe -n $n ./pypresso ../testsuite/python/lees_edwards.py
done
```
Expected: PASS at 1/2/3/4 ranks. `split_shear_rejected` now raises and passes; serial, split-normal, split-neutral all find every pair.

- [ ] **Step 5: Verify no regression in the non-LE and fully_connected paths**

Run:
```bash
cd build && ctest -R "lees_edwards$|regular_decomposition|random_pairs|hybrid_decomposition" --output-on-failure
```
Expected: all PASS.

- [ ] **Step 6: Format and commit**

```bash
maintainer/format/clang-format.sh src/core/integrate.cpp
git -C ... add src/core/integrate.cpp testsuite/python/lees_edwards.py
git -C ... commit -m "Guard Lees-Edwards to node_grid[shear]==1 at integration; verify full matrix"
```

---

## Phase 4: Prove the three criteria

### Task 4.1: Bitwise trajectory identity, dynamic vs fully_connected, per fixed decomposition

**Files:**
- Use: `maintainer/benchmarks/dpd_le_verify.py` (Task 2.1)

**Interfaces:**
- Consumes: `dump`/`compare` modes; the `--bitwise` flag; `compare` asserts equal step-0 within-cutoff pair sets, both `trace_ok`, bitwise positions, and a shear-stress z-score.

- [ ] **Step 1: Serial, dynamic vs fully_connected**

Run:
```bash
B=.../build/pypresso ; V=.../maintainer/benchmarks/dpd_le_verify.py
$B $V dump --tag dyn_1rank --node_grid 1 1 1 --steps 200 --out /tmp
$B $V dump --tag fc_1rank  --fully_connected --node_grid 1 1 1 --steps 200 --out /tmp
$B $V compare --a /tmp/dyn_1rank.npz --b /tmp/fc_1rank.npz --bitwise
```
Expected: "identical step-0 within-cutoff pair sets", "pos: bitwise identical", stress z-score below 5.

- [ ] **Step 2: y-split (node_grid[sd]==1), dynamic vs fully_connected, per rank count**

Run (2 ranks; repeat for 4):
```bash
mpiexec -n 2 $B $V dump --tag dyn_y2 --node_grid 1 2 1 --steps 200 --out /tmp
mpiexec -n 2 $B $V dump --tag fc_y2  --fully_connected --node_grid 1 2 1 --steps 200 --out /tmp
$B $V compare --a /tmp/dyn_y2.npz --b /tmp/fc_y2.npz --bitwise
```
Expected: bitwise identical (same decomposition, dynamic vs reference).

- [ ] **Step 3: OpenMP thread variation, bitwise per (ranks, threads)**

Run each `dump` twice — once with `OMP_NUM_THREADS=1`, once with `=4` — pairing dynamic vs fully_connected at the *same* thread count, and compare `--bitwise`. Expected: bitwise identical at each fixed (ranks, threads).

- [ ] **Step 4: Record results** in the task log (no code change). If any bitwise comparison fails, STOP and debug with `superpowers:systematic-debugging` before proceeding — bitwise identity per fixed decomposition is a hard acceptance criterion.

### Task 4.2: Trace pair-set equality across all decompositions (permanent tests)

**Files:**
- Modify: `testsuite/python/lees_edwards.py` — replace the obsolete `test_zz_lj_pair_visibility` expectation (lines 979-1013) and add the permanent DPD matrix.

- [ ] **Step 1: Confirm the inverted negative assertion (likely already done in 3.2)**

The old `test_zz_lj_pair_visibility` used `assertRaises(AssertionError)` to document that a regular decomposition *without* `fully_connected` FAILS to find the pair. Task 3.2 already had to invert this (the dynamic serial path now succeeds), changing it to a positive check on `node_grid = [1, self.n_nodes, 1]` (split-normal). Verify that inversion is present and uses **split-normal** (`[1, self.n_nodes, 1]`), NOT split-shear. If for any reason it still references `assertRaises` or `[self.n_nodes, 1, 1]`, fix it to:
```python
        # Dynamic sheared halo (no fully_connected) now finds the pair on the
        # allowed grids (shear axis unsplit).
        system.cell_system.set_regular_decomposition(
            fully_connected_boundary=None)
        system.cell_system.node_grid = [1, self.n_nodes, 1]  # split along normal
        self.run_lj_pair_visibility("x", "y")
```

- [ ] **Step 2: Add the permanent DPD matrix test (allowed grids only)**

```python
    @utx.skipIfMissingFeatures(["DPD"])
    def test_dpd_pair_visibility_matrix(self):
        system = self.system
        system.box_l = [10, 10, 10]
        grids = {"y": [1, self.n_nodes, 1],   # split along normal
                 "z": [1, 1, self.n_nodes]}   # split along neutral
        for verlet in (False, True):
            for split, grid in grids.items():
                system.cell_system.node_grid = grid
                system.cell_system.set_regular_decomposition(
                    use_verlet_lists=verlet, fully_connected_boundary=None)
                self.run_dpd_pair_visibility("x", "y")
```
Splitting along the shear direction is out of scope and covered by the
`test_dpd_split_shear_rejected` guard test (Task 3.4).

- [ ] **Step 3: Run on 1, 2, 4 ranks**

Run:
```bash
cd build && make -j8 && make python_test_data
ctest -R "lees_edwards$" --output-on-failure           # default rank count
```
Also run the file directly under specific rank counts:
```bash
for n in 1 2 4; do mpiexec --oversubscribe -n $n build/pypresso testsuite/python/lees_edwards.py; done
```
Expected: all PASS.

- [ ] **Step 4: Remove the leftover DPD helper temporary test** if any remain; format; commit.

```bash
maintainer/format/autopep8.sh testsuite/python/lees_edwards.py
git -C ... add testsuite/python/lees_edwards.py
git -C ... commit -m "Prove Lees-Edwards pair visibility across all decompositions"
```

### Task 4.3: Average shear-stress agreement

**Files:**
- Use: `maintainer/benchmarks/dpd_le_verify.py`

- [ ] **Step 1: Compare average shear stress across decompositions**

Run (all with the same seed/box, longer run for statistics, e.g. `--steps 2000`):
```bash
$B $V dump --tag st_1rank --node_grid 1 1 1 --steps 2000 --out /tmp
mpiexec -n 2 $B $V dump --tag st_y2 --node_grid 1 2 1 --steps 2000 --out /tmp
mpiexec -n 2 $B $V dump --tag st_z2 --node_grid 1 1 2 --steps 2000 --out /tmp
$B $V compare --a /tmp/st_1rank.npz --b /tmp/st_y2.npz
$B $V compare --a /tmp/st_1rank.npz --b /tmp/st_z2.npz
```
Expected: both `trace_ok`; identical step-0 within-cutoff pair sets; shear-stress z-score below `--stress_z` (default 5). Grids are the allowed ones (shear axis `x` unsplit): serial, split-normal `[1,2,1]`, split-neutral `[1,1,2]`. All dumps use the dynamic path — no `--fully_connected` — so each must pass its own step-0 self-check, which only holds after Phase 3.

- [ ] **Step 2: Record results** in the task log.

---

## Phase 5: Remove `fully_connected_boundary`

**User decision (2026-08-30):** hard-remove `fully_connected_boundary`. The
Python/script-interface setter no longer accepts it — passing a non-`None`
`fully_connected_boundary` raises a clear error. The C++ implementation
(member, ctor parameter, stencil/halo branches, getter) is deleted, and the
`le_shear()` gate simplifies to "active iff the box is Lees-Edwards".

**Ordering:** this phase runs only after Phase 4 is complete and its bitwise
proof is recorded, because Phase 4 uses `fully_connected` as the reference.

### Task 5.1: Make the script-interface setter reject `fully_connected_boundary`

**Files:**
- Modify: `src/script_interface/cell_system/CellSystem.cpp` — the setter block (lines 308-322) and the read-only accessor (lines 142-155).
- Modify: `src/python/espressomd/cell_system.py` — remove the `fully_connected_boundary` docstring entry (line ~110) and any argument plumbing in `set_regular_decomposition`.

**Interfaces:**
- Consumes: nothing new.
- Produces: `set_regular_decomposition(fully_connected_boundary=<not None>)` raises `RuntimeError` with a message pointing to the automatic dynamic behavior; `system.cell_system.fully_connected_boundary` is removed.

- [ ] **Step 1: Write the failing test**

Append to `testsuite/python/lees_edwards.py`:
```python
    def test_fully_connected_removed(self):
        system = self.system
        with self.assertRaisesRegex(Exception, "no longer"):
            system.cell_system.set_regular_decomposition(
                fully_connected_boundary={"boundary": "y", "direction": "x"})
```

- [ ] **Step 2: Run it, expect fail** (the old setter silently accepts the dict).

Run: `cd build && ctest -R "lees_edwards$" --output-on-failure`
Expected: FAIL (no exception raised).

- [ ] **Step 3: Make the setter throw**

Replace the setter block in `CellSystem.cpp` (lines 308-322) with:
```cpp
  } else if (cs_type == CellStructureType::REGULAR) {
    if (params.contains("fully_connected_boundary") and
        not is_none(params.at("fully_connected_boundary"))) {
      throw std::runtime_error(
          "fully_connected_boundary is no longer supported: Lees-Edwards "
          "shear now updates the ghost neighborships dynamically, so it is "
          "not needed. Remove this argument.");
    }
    context()->parallel_try_catch([this]() {
      m_cell_structure->set_regular_decomposition(
          get_system().get_interaction_range());
    });
  } else {
```
Remove the read-only `"fully_connected_boundary"` accessor (lines 142-155) from the parameter list. Remove the docstring entry in `cell_system.py`.

- [ ] **Step 4: Run it, expect pass** (build first).

Run: `cd build && make -j8 && ctest -R "lees_edwards$" --output-on-failure`
Expected: `test_fully_connected_removed` PASSES. (Other LE tests still reference the old API at this point; they are fixed in Task 5.3. Run just the new test in isolation: append `-- --gtest_filter` is not applicable; instead run `build/pypresso testsuite/python/lees_edwards.py LeesEdwards.test_fully_connected_removed`.)

- [ ] **Step 5: Format and commit**

```bash
maintainer/format/clang-format.sh src/script_interface/cell_system/CellSystem.cpp
maintainer/format/autopep8.sh src/python/espressomd/cell_system.py testsuite/python/lees_edwards.py
git -C ... add src/script_interface/cell_system/CellSystem.cpp src/python/espressomd/cell_system.py testsuite/python/lees_edwards.py
git -C ... commit -m "Reject fully_connected_boundary at the script interface"
```

### Task 5.2: Delete the C++ `fully_connected` implementation

**Files:**
- Modify: `src/core/cell_system/RegularDecomposition.hpp` — remove `m_fully_connected_boundary` (line 84), the `fully_connected_boundary()` getter (line 121), and the ctor's `fully_connected` parameter (line ~97).
- Modify: `src/core/cell_system/RegularDecomposition.cpp` — remove the ctor param and initializer (lines 818-823); delete the fully-connected sanity-check block (lines 544-559), the `fcb_is_inner_connection` lambda (lines 508-521), the fully-connected stencil branch (lines 579-587) and its filter use (lines 609-614) in `init_cell_interactions()`; simplify `le_shear()` (Task 3.1) to gate on `m_box.type() == BoxType::LEES_EDWARDS` only. **KEEP** the dynamic-path `node_grid[shear_dir] != 1` guard added in Task 3.4 — it is a separate, still-required check (do not delete it along with the fully-connected sanity block).
- Modify: `src/core/cell_system/CellStructure.hpp` (lines 829-833) and `CellStructure.cpp` (lines 692-697) — remove the `fully_connected_boundary` parameter from `set_regular_decomposition`.
- Modify: `src/core/system/System.cpp` (lines 182-194) — drop the fully-connected preservation; call `set_regular_decomposition(get_interaction_range())` in both branches.

**Interfaces:**
- Consumes: `le_shear()` now activates purely on box type.
- Produces: `set_regular_decomposition(double range)` (no second parameter).

- [ ] **Step 1: Remove the getter's only remaining consumer**

Confirm `le_shear()` is the only place still calling `fully_connected_boundary()` after Task 5.1:
Run: `grep -rn "fully_connected_boundary()" src/core`
Expected: only `RegularDecomposition.cpp` (the `le_shear()` gate). Edit `le_shear()`:
```cpp
  if (m_box.type() != BoxType::LEES_EDWARDS)
    return r;
```

- [ ] **Step 2: Delete the members, params, and branches** listed in Files. After edits:
Run: `grep -rn "fully_connected\|m_fully_connected\|fcb_is_inner_connection" src/core`
Expected: no matches.

- [ ] **Step 3: Build**

Run: `cd build && make -j8`
Expected: compiles (all call sites updated).

- [ ] **Step 4: Run the core-affected tests**

Run: `ctest -R "lees_edwards$|regular_decomposition|random_pairs|hybrid_decomposition" --output-on-failure`
Expected: `lees_edwards` may still fail on old-API references in the test files (fixed in 5.3); `random_pairs`/`hybrid_decomposition` PASS. `regular_decomposition` is fixed in 5.3.

- [ ] **Step 5: Format and commit**

```bash
maintainer/format/clang-format.sh src/core/cell_system/RegularDecomposition.hpp src/core/cell_system/RegularDecomposition.cpp src/core/cell_system/CellStructure.hpp src/core/cell_system/CellStructure.cpp src/core/system/System.cpp
git -C ... add src/core/cell_system/ src/core/system/System.cpp
git -C ... commit -m "Delete the fully_connected_boundary implementation"
```

### Task 5.3: Clean up tests and scripts that used `fully_connected`

**Files:**
- Modify: `testsuite/python/regular_decomposition.py` — remove/convert the fully-connected tests (the `non_bonded_loop_trace` fully-connected cases around lines 130-331).
- Modify: `testsuite/python/lees_edwards.py` — remove the fully-connected loop in `test_zz_lj_pair_visibility` (lines 998-1013; the positive dynamic check from Task 4.2 replaces it).
- Modify: `maintainer/benchmarks/dpd.py` — remove the `--fully_connected` option and `configure_decomposition`'s fully-connected branch (now always dynamic).
- Modify: `maintainer/benchmarks/dpd_le_verify.py` — remove the `--fully_connected` option and its branch.

**Interfaces:**
- Consumes: nothing new.
- Produces: a test/benchmark suite with no `fully_connected` references.

- [ ] **Step 1: Convert `regular_decomposition.py`**

For each test that set `fully_connected_boundary` to verify cross-boundary visibility, either delete it (if it duplicates the LE dynamic coverage now in `lees_edwards.py`) or convert it to assert the setter now raises. Keep the non-LE regular-decomposition coverage intact. After editing:
Run: `grep -n "fully_connected" testsuite/python/regular_decomposition.py`
Expected: no matches (or only an assert-raises test).

- [ ] **Step 2: Simplify the benchmark scripts**

In `dpd.py`, delete the `--fully_connected` argument and collapse `configure_decomposition` to the dynamic branch only. In `dpd_le_verify.py`, delete the `--fully_connected` argument and the `fully_connected` parameter of `build_system` (always dynamic). After editing:
Run: `grep -rn "fully_connected" maintainer/benchmarks/`
Expected: no matches.

- [ ] **Step 3: Full affected test run on 1, 2, 4 ranks**

Run:
```bash
cd build && make -j8
ctest -R "lees_edwards$|regular_decomposition|random_pairs|hybrid_decomposition" --output-on-failure
mpiexec -n 2 build/pypresso testsuite/python/lees_edwards.py
mpiexec -n 4 build/pypresso testsuite/python/lees_edwards.py
```
Expected: all PASS.

- [ ] **Step 4: Confirm no references remain anywhere**

Run: `grep -rln "fully_connected\|fully connected" src/ testsuite/ maintainer/ doc/`
Expected: no matches (docs handled in 5.4).

- [ ] **Step 5: Format and commit**

```bash
maintainer/format/autopep8.sh testsuite/python/regular_decomposition.py testsuite/python/lees_edwards.py maintainer/benchmarks/dpd.py maintainer/benchmarks/dpd_le_verify.py
git -C ... add testsuite/python/ maintainer/benchmarks/
git -C ... commit -m "Remove fully_connected usage from tests and benchmarks"
```

### Task 5.4: Documentation and changelog

**Files:**
- Modify: `doc/sphinx/` (search for the Lees-Edwards and cell-system sections) and `CHANGELOG.md`.

- [ ] **Step 1: Find the doc references**

Run: `grep -rln "fully_connected\|Lees-Edwards\|lees_edwards" doc/`
Expected: a Lees-Edwards user-guide section and a cell-system section.

- [ ] **Step 2: Update the prose** to state: Lees-Edwards shear updates ghost neighborships dynamically as the offset evolves; `fully_connected_boundary` has been removed and is no longer needed; the MPI node grid may now be split along the shear direction. Add a `CHANGELOG.md` entry under the current development version noting the removal and the new split-along-shear capability.

- [ ] **Step 3: Verify docs build** (if a docs target exists):
Run: `grep -rn "sphinx\|doc" build/CMakeCache.txt | head` then `make -C build sphinx` if available; otherwise skip.

- [ ] **Step 4: Commit**

```bash
git -C ... add doc/ CHANGELOG.md
git -C ... commit -m "Document dynamic Lees-Edwards halo and fully_connected removal"
```

---

## Self-review notes

- **Spec coverage (revised scope):** shifted sheared stencil + shear-axis fold (Task 3.2), rebuild-on-resort (3.3), shear-axis node-grid guard + full verification (3.4), skin-sized reach (`le_reach` uses `max_range`, which includes skin), LE-driven activation (3.1 gate), dpd.py (1.1), verification of all three criteria (4.1 bitwise, 4.2 trace, 4.3 stress) across serial/split-normal/split-neutral + OpenMP, removal of `fully_connected` (Phase 5: 5.1 setter throws, 5.2 C++ deleted but the 3.4 shear-axis guard kept, 5.3 tests/scripts cleaned, 5.4 docs). `make_halo_plan` is intentionally NOT modified.
- **Scope:** splitting along the shear direction is out of scope (user decision 2026-08-31); `node_grid[shear_dir]==1` is enforced by the 3.4 guard. This matches `fully_connected`'s existing limit, so removing the attribute is not a regression.
- **Removal ordering:** Phase 5 deletes the `fully_connected` reference path, so it must run strictly after Phase 4 records the bitwise proof. The plan states this in the Phase 5 header and Global Constraints.
- **Bitwise-only-where-valid:** Task 4.1 compares dynamic vs fully_connected only for `node_grid[sd]==1` grids (serial, y-split), matching the spec's constraint. Split-along-shear is proven by trace + stress only (4.2, 4.3).
- **Coexistence:** the dynamic path is gated off when `fully_connected_boundary()` is set (Global Constraints + 3.1), so the reference and the new path coexist through Phase 4.
- **Signature consistency:** `le_shear()`/`LeShear` (3.1) is consumed by 3.2, 3.3, 3.4; `rebuild_topology()`/`m_le_shift_at_last_build` (3.3) are used only within the class; `run_dpd_pair_visibility` (2.2) is consumed by 3.2, 3.4, 4.2; `dump`/`compare` (2.1) by 4.1, 4.3.
- **Open risk flagged in-plan:** the exact `ESPRESSO_MYCONFIG` cmake variable name (0.1 Step 1) and the presence of a C++ unit-test target (3.1 Step 1) are verified at execution time with the given `grep`/`cmake -L` commands rather than assumed.
