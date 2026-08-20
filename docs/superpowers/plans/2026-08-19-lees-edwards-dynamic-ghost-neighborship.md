# Lees-Edwards Dynamic Ghost Neighborship Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the physical-cell -> ghost-cell wiring for Lees-Edwards (LE) shear a function of the current position offset, rebuilt as the offset evolves, so the short-range loop sees every in-range pair across the shear boundary without `fully_connected_boundary` and even when the MPI grid splits along the shear direction.

**Architecture:** Keep the LE offset in `LeesEdwardsBC::distance()` (ghost positions carry only periodic-image shifts). Change three things in `RegularDecomposition`: (1) `make_halo_plan()` sources each shear-crossing ghost from the offset-shifted source column and owning rank; (2) `init_cell_interactions()` wires a narrow offset-centered window instead of the full shear row; (3) `resort()` rebuilds neighborships and the halo plan when the integer offset shift changes. Correctness is proven by bitwise-identical trajectories (per fixed decomposition, vs the retained `fully_connected` reference), by `non_bonded_loop_trace` pair-set equality across decompositions, and by average shear stress.

**Tech Stack:** C++20 (ESPResSo core, Kokkos, boost::mpi), Python 3 (espressomd, numpy), CMake, ctest, mpiexec.

**Spec:** `docs/superpowers/specs/2026-08-19-lees-edwards-dynamic-ghost-neighborship-design.md`

## Global Constraints

- Build with the `maxset` feature config (`maintainer/configs/maxset.hpp`) which enables `DPD` and `LENNARD_JONES`. Build with **CUDA off** for this work.
- Use `make -j8` (never `-j$(nproc)`).
- Use `git -C <path> ...`, never `cd <path> && git ...`.
- Run `maintainer/format` (clang-format / autopep8) on changed files before every commit, or CI will not start.
- The LE offset stays in `LeesEdwardsBC::distance()`. Ghost positions carry only periodic-image shifts (multiples of the box length), never the LE offset.
- Bitwise identity vs the old code is only defined where `fully_connected` runs, i.e. `node_grid[shear_dir] == 1`. Across different decompositions only the pair set (trace) and the average shear stress are expected to match.
- The dynamic path engages only when `m_box.type() == BoxType::LEES_EDWARDS` **and** `fully_connected_boundary()` is unset. When `fully_connected_boundary()` is set, keep the exact current behavior (it is the bitwise reference during Phases 1-4). Non-LE behavior is unchanged.
- Notation used throughout: `sn = shear_plane_normal`, `sd = shear_direction`, `cs = cell_size[sd]`, `O = pos_offset`, `s = lround(O / cs)` (integer column shift), `frac = O - s*cs` (residual in `[-cs/2, cs/2]`).

---

## Phase 0: Build and baseline

### Task 0.1: Configure and build with maxset, CUDA off

**Files:**
- Create: `build/` (out-of-source build dir in the worktree)

- [ ] **Step 1: Configure**

Run:
```bash
cmake -S /tikhome/weeber/es/.claude/worktrees/comm_le \
      -B /tikhome/weeber/es/.claude/worktrees/comm_le/build \
      -D ESPRESSO_BUILD_WITH_CUDA=OFF \
      -D ESPRESSO_BUILD_WITH_CALIPER=OFF \
      -D CMAKE_BUILD_TYPE=Release \
      -D ESPRESSO_TEST_NP=4 \
      -D CMAKE_CXX_FLAGS="-DMYCONFIG_H" \
      -D ESPRESSO_MYCONFIG=/tikhome/weeber/es/.claude/worktrees/comm_le/maintainer/configs/maxset.hpp
```
Expected: configuration succeeds, `DPD` and `LENNARD_JONES` reported enabled. (If the exact `-D ESPRESSO_MYCONFIG` flag name differs in this tree, run `cmake -LAH build | grep -i myconfig` to find the correct variable and re-run.)

- [ ] **Step 2: Build**

Run: `make -j8 -C /tikhome/weeber/es/.claude/worktrees/comm_le/build`
Expected: `build/pypresso` exists and the build completes without error.

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

**Interfaces:**
- Produces: a script invoked as `pypresso dpd_le_verify.py --mode {dump,compare} --tag NAME [--fully_connected] --node_grid a b c --steps N --out DIR`. In `dump` mode it builds a deterministic DPD+LE system, integrates `--steps`, and writes `DIR/NAME.npz` containing `pos` (final folded positions sorted by particle id), `trace` (sorted pair-set array), and `stress` (time-averaged shear-stress scalar). In `compare` mode it loads two `.npz` files and asserts bitwise pos equality, exact trace set equality, and stress closeness.

- [ ] **Step 1: Write the script**

Create `maintainer/benchmarks/dpd_le_verify.py` with the GPL header, then:

```python
import argparse
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


def build_system(node_grid, fully_connected, density, box_l, shear_velocity):
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
        initial_pos_offset=0.0, shear_velocity=shear_velocity)
    system.lees_edwards.set_boundary_conditions(
        shear_direction=SHEAR_DIRECTION,
        shear_plane_normal=SHEAR_PLANE_NORMAL, protocol=protocol)
    return system


def collect_trace(system):
    pairs = system.cell_system.non_bonded_loop_trace()
    rows = []
    for id1, id2, _p1, _p2, vec21, _node in pairs:
        a, b = (id1, id2) if id1 < id2 else (id2, id1)
        rows.append((a, b))
    rows.sort()
    return np.array(rows, dtype=int)


def shear_stress(system):
    t = system.analysis.pressure_tensor()["total"]
    return float(t[AXIS[SHEAR_DIRECTION]][AXIS[SHEAR_PLANE_NORMAL]])


def dump(args):
    system = build_system(args.node_grid, args.fully_connected, args.density,
                          args.box_l, args.shear_velocity)
    stresses = []
    for _ in range(args.steps):
        system.integrator.run(1)
        stresses.append(shear_stress(system))
    p = system.part.all()
    order = np.argsort(p.id)
    np.savez(f"{args.out}/{args.tag}.npz",
             pos=np.copy(p.pos_folded)[order],
             trace=collect_trace(system),
             stress=np.mean(stresses))
    print(f"dumped {args.tag}")


def compare(args):
    a = np.load(args.a, allow_pickle=True)
    b = np.load(args.b, allow_pickle=True)
    ta = {tuple(r) for r in a["trace"]}
    tb = {tuple(r) for r in b["trace"]}
    assert ta == tb, f"trace mismatch: {ta ^ tb}"
    print("trace: identical pair sets")
    if args.bitwise:
        assert np.array_equal(a["pos"], b["pos"]), "positions not bitwise equal"
        print("pos: bitwise identical")
    ds = abs(float(a["stress"]) - float(b["stress"]))
    print(f"shear stress: |{float(a['stress']):.6e} - "
          f"{float(b['stress']):.6e}| = {ds:.3e}")
    assert ds < args.stress_tol, "shear stress mismatch"


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
d.set_defaults(func=dump)
c = sub.add_parser("compare")
c.add_argument("--a", required=True)
c.add_argument("--b", required=True)
c.add_argument("--bitwise", action="store_true")
c.add_argument("--stress_tol", type=float, default=5e-2)
c.set_defaults(func=compare)
args = p.parse_args()
args.func(args)
```

- [ ] **Step 2: Format** — run autopep8 on the file.

- [ ] **Step 3: Self-consistency check against the reference (serial vs y-split fully_connected)**

Run:
```bash
B=/tikhome/weeber/es/.claude/worktrees/comm_le/build/pypresso
V=/tikhome/weeber/es/.claude/worktrees/comm_le/maintainer/benchmarks/dpd_le_verify.py
$B $V dump --tag ref1 --fully_connected --node_grid 1 1 1 --out /tmp
mpiexec -n 2 $B $V dump --tag ref2 --fully_connected --node_grid 1 2 1 --out /tmp
$B $V compare --a /tmp/ref1.npz --b /tmp/ref2.npz
```
Expected: "trace: identical pair sets" and a small shear-stress difference (no `--bitwise`, since the decompositions differ). This proves the harness works on known-good code.

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

### Task 3.2: Dynamic neighbor stencil for the local (node_grid[sd]==1) case

**Files:**
- Modify: `src/core/cell_system/RegularDecomposition.cpp` — `init_cell_interactions()` (lines 486-658)

**Interfaces:**
- Consumes: `le_shear()` (Task 3.1); `max_range()` for the interaction reach; existing `at_boundary`, `folded_linear_index`, `fcb_is_inner_connection`, red/black machinery.
- Produces: neighbor lists that, for the LE-active case, wire a narrow offset-centered window along `sd` at the `sn` boundary, instead of the fully-connected full row. The fully-connected branch is left intact for the reference path.

- [ ] **Step 1: Write the failing Python trace test**

Append to `testsuite/python/lees_edwards.py` inside the class:
```python
    @utx.skipIfMissingFeatures(["DPD"])
    def test_dpd_dynamic_local(self):
        """Dynamic sheared halo, no fully_connected, split along the
        shear-plane-normal only (node_grid[sd]==1)."""
        system = self.system
        system.box_l = [10, 10, 10]
        system.cell_system.node_grid = [1, self.n_nodes, 1]
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True, fully_connected_boundary=None)
        self.run_dpd_pair_visibility("x", "y")
```

- [ ] **Step 2: Run it, expect fail**

Run:
```bash
cd /tikhome/weeber/es/.claude/worktrees/comm_le/build && make -j8 && \
  ctest -R "lees_edwards$" --output-on-failure
```
Expected: FAIL — without the dynamic stencil, cross-shear-boundary pairs are missed ("Pair not found by the core").

- [ ] **Step 3: Implement the dynamic window**

In `init_cell_interactions()`, compute `auto const le = le_shear();` near the top (after `global_size`). Compute the shear-direction reach in cells from the interaction range:
```cpp
  auto const le_reach =
      le.active
          ? static_cast<int>(std::ceil(max_range()[le.sd] / cell_size[le.sd])) + 1
          : 0;
```
Then, in the per-cell loop, add an LE branch parallel to the existing fully-connected branch (keep the latter unchanged):
```cpp
        if (le.active) {
          // Narrow, offset-centered window along the shear direction, only
          // at the shear-plane-normal boundary. Center on the integer-shifted
          // column; width covers the cutoff+skin reach plus the fractional
          // offset (the +1 in le_reach).
          if (at_boundary(le.sn, {m, n, o})) {
            lower_index[le.sd] = (m_or_n_or_o_along_sd(m, n, o, le.sd)) - le.shift - le_reach;
            upper_index[le.sd] = (m_or_n_or_o_along_sd(m, n, o, le.sd)) - le.shift + le_reach;
          }
        }
```
where the current index along `sd` is selected directly (replace the pseudo-call with the literal component): use `{m, n, o}[le.sd]`, i.e.
```cpp
            auto const cur = Utils::Vector3i{m, n, o}[le.sd];
            lower_index[le.sd] = cur - le.shift - le_reach;
            upper_index[le.sd] = cur - le.shift + le_reach;
```
No `fcb_is_inner_connection` filter is applied in the LE branch (the window is already narrow). The existing `folded_linear_index` set and red/black split then order the neighbors exactly as before.

- [ ] **Step 4: Run the test, expect pass**

Run: `cd build && make -j8 && ctest -R "lees_edwards$" --output-on-failure`
Expected: PASS at the initial offset. (Motion over many steps is covered once Task 3.3 lands; if the 50-step loop in `run_dpd_pair_visibility` still fails here due to offset growth, that is expected and fixed next — temporarily reduce the loop to 1 iteration to confirm the stencil itself is correct, then restore.)

- [ ] **Step 5: Format C++**

Run: `/tikhome/weeber/es/.claude/worktrees/comm_le/maintainer/format/clang-format.sh src/core/cell_system/RegularDecomposition.cpp` (or the repo wrapper).

- [ ] **Step 6: Commit**

```bash
git -C ... add src/core/cell_system/RegularDecomposition.cpp testsuite/python/lees_edwards.py
git -C ... commit -m "Dynamic Lees-Edwards neighbor stencil for local shear halo"
```

### Task 3.3: Rebuild neighborships and halo plan on resort

**Files:**
- Modify: `src/core/cell_system/RegularDecomposition.hpp` (add member `int m_le_shift_at_last_build = 0;` in the private section)
- Modify: `src/core/cell_system/RegularDecomposition.cpp` — `resort()` (lines 187-...); factor the constructor's build sequence (lines 829-935: `init_cell_interactions(); mark_cells(); make_halo_plan(); ...boundary marking...`) into a private `rebuild_topology()` method and call it from both the constructor and `resort()`.

**Interfaces:**
- Consumes: `le_shear()` (3.1).
- Produces: a private `void rebuild_topology();` that runs `init_cell_interactions()`, `mark_cells()`, `m_halo_plan = make_halo_plan()`, and the boundary-marking block. `resort()` calls it when `le_shear().active && le_shear().shift != m_le_shift_at_last_build`. `m_le_shift_at_last_build` is updated inside `rebuild_topology()`.

- [ ] **Step 1: Restore the full 50-step loop** in `run_dpd_pair_visibility` if it was shortened, and confirm `test_dpd_dynamic_local` now FAILS at a later offset.

Run: `ctest -R "lees_edwards$" --output-on-failure`
Expected: FAIL after several integrator runs (the stencil built at offset 0 goes stale as the offset grows past `cell_size[sd]`).

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

- [ ] **Step 4: Run the test, expect pass**

Run: `cd build && make -j8 && ctest -R "lees_edwards$" --output-on-failure`
Expected: PASS across all 50 integrator runs (offset now triggers a topology rebuild whenever the integer column shift changes; between rebuilds the `cutoff+skin`-sized window and the resort trigger's `skin/2` drift bound keep every pair visible).

- [ ] **Step 5: Format and commit**

```bash
maintainer/format/clang-format.sh src/core/cell_system/RegularDecomposition.cpp
git -C ... add src/core/cell_system/RegularDecomposition.hpp src/core/cell_system/RegularDecomposition.cpp
git -C ... commit -m "Rebuild Lees-Edwards topology on integer-offset crossing"
```

### Task 3.4: Offset-driven halo sourcing for the split-along-shear case

**Files:**
- Modify: `src/core/cell_system/RegularDecomposition.cpp` — `make_halo_plan()` (lines 660-816), the ghost-enumeration loop (lines 737-787) and `SendRegion` construction (line 803)
- Modify: `src/core/cell_system/RegularDecomposition.cpp` — remove the `node_grid[fc_dir] != 1` throw's applicability to the LE-dynamic path (the throw at lines 550-553 is guarded by `fully_connected_boundary()`, so it already does not fire for the dynamic path — verify and add a comment).

**Interfaces:**
- Consumes: `le_shear()` (3.1); existing `owner_of`, `global_key`, `side`, `mc`, `SendRegion`, and `m_box.length()` for the periodic image shift.
- Produces: for shear-plane-normal-crossing ghosts, a `ghost_global[sd]` shifted by `-side[sn] * shift`, resolving to the sheared source column and its owning rank; any `sd` wrap past the box becomes an image shift stored in `SendRegion.shift`.

- [ ] **Step 1: Write the failing split-along-shear trace test**

Append to `testsuite/python/lees_edwards.py`:
```python
    @utx.skipIfMissingFeatures(["DPD"])
    def test_dpd_dynamic_split_shear(self):
        """Dynamic sheared halo, MPI grid split ALONG the shear direction
        (impossible with fully_connected)."""
        if self.n_nodes < 2:
            self.skipTest("needs >= 2 MPI ranks")
        system = self.system
        system.box_l = [10, 10, 10]
        system.cell_system.node_grid = [self.n_nodes, 1, 1]  # split along x = sd
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True, fully_connected_boundary=None)
        self.run_dpd_pair_visibility("x", "y")
```

- [ ] **Step 2: Run it, expect fail**

Run: `ctest -R "lees_edwards$" --output-on-failure` (with 2 or 4 ranks)
Expected: FAIL — the halo along `sd` currently sources the geometrically-adjacent column, not the sheared one, so cross-boundary pairs on split ranks are missed.

- [ ] **Step 3: Implement the sheared sourcing**

In `make_halo_plan()`, compute `auto const le = le_shear();` after `global_origin`. In the ghost loop, after computing `side` and before `ghost_global`, apply the shear shift to the mirrored global index when this ghost crosses the shear-plane-normal boundary:
```cpp
        auto ghost_global = global_origin + (nc - one);
        Utils::Vector3d region_shift{};
        if (le.active and side[le.sn] != 0) {
          ghost_global[le.sd] += -side[le.sn] * le.shift;
          // A shifted column that leaves [0, global_size[sd]) wraps the
          // x-periodic boundary; record the corresponding image shift so the
          // packed ghost positions land in the right periodic image. This is
          // a multiple of the box length, NOT the LE offset.
          auto const gs = global_size[le.sd];
          auto const wraps = static_cast<int>(
              std::floor(static_cast<double>(ghost_global[le.sd]) / gs));
          region_shift[le.sd] = static_cast<double>(wraps) * m_box.length()[le.sd];
        }
        auto const peer = owner_of(ghost_global);
        auto const recv_key = global_key(ghost_global);
```
Mirror the same shift when building the dual send cell so `recv[k]` still aligns with the peer's `send[k]`: apply `-side[le.sn] * le.shift` to `send_global[le.sd]` (via `mc`) as well. Then construct the send region with the image shift:
```cpp
        peers[peer].send.emplace_back(global_key(send_global),
                                      SendRegionInfo{list_at(mc), region_shift});
```
Adjust the `PeerBucket::send` element type to carry the shift alongside the key (a small local struct `SendRegionInfo{ParticleList*, Utils::Vector3d}`), and in the emission loop (lines 793-804) push `SendRegion{info.cell, info.shift}` instead of `SendRegion{cell, {}}`.

Heed the `HaloExchange.cpp` note: if shifts differ within one `NeighborComm`, `pack_regions` must still write all bonds into one shared archive. For DPD there are no bonds; but keep regions with distinct shifts grouped so the existing per-neighbor common-shift assumption holds — if a single peer receives columns with different image shifts, split them into separate `NeighborComm`s keyed by (peer, shift). Implement the split: key `peers` by `std::pair<int, Utils::Vector3d-as-ivec>` when `le.active`, else by `peer` as today.

- [ ] **Step 4: Run the test, expect pass**

Run: `cd build && make -j8 && ctest -R "lees_edwards$" --output-on-failure` (2 and 4 ranks)
Expected: PASS — split-along-shear now finds every pair.

- [ ] **Step 5: Verify no regression in the non-LE and fully_connected paths**

Run:
```bash
cd build && ctest -R "lees_edwards$|regular_decomposition|random_pairs|hybrid_decomposition" --output-on-failure
```
Expected: all PASS.

- [ ] **Step 6: Format and commit**

```bash
maintainer/format/clang-format.sh src/core/cell_system/RegularDecomposition.cpp
git -C ... add src/core/cell_system/RegularDecomposition.cpp testsuite/python/lees_edwards.py
git -C ... commit -m "Offset-driven Lees-Edwards halo sourcing across the shear boundary"
```

---

## Phase 4: Prove the three criteria

### Task 4.1: Bitwise trajectory identity, dynamic vs fully_connected, per fixed decomposition

**Files:**
- Use: `maintainer/benchmarks/dpd_le_verify.py` (Task 2.1)

**Interfaces:**
- Consumes: `dump`/`compare` modes; the `--bitwise` flag.

- [ ] **Step 1: Serial, dynamic vs fully_connected**

Run:
```bash
B=.../build/pypresso ; V=.../maintainer/benchmarks/dpd_le_verify.py
$B $V dump --tag dyn_1rank --node_grid 1 1 1 --steps 200 --out /tmp
$B $V dump --tag fc_1rank  --fully_connected --node_grid 1 1 1 --steps 200 --out /tmp
$B $V compare --a /tmp/dyn_1rank.npz --b /tmp/fc_1rank.npz --bitwise
```
Expected: "pos: bitwise identical", "trace: identical pair sets".

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

- [ ] **Step 1: Invert the obsolete negative assertion**

The current `test_zz_lj_pair_visibility` asserts (lines 985-990) that a regular decomposition *without* `fully_connected` FAILS to find the pair. After the refactor that must succeed. Replace that `assertRaises` block with a positive check that the dynamic decomposition finds all pairs:
```python
        # Dynamic sheared halo (no fully_connected) must now find the pair,
        # including when the grid splits along the shear direction.
        system.cell_system.set_regular_decomposition(
            fully_connected_boundary=None)
        system.cell_system.node_grid = [self.n_nodes, 1, 1]  # split along shear
        self.run_lj_pair_visibility("x", "y")
```

- [ ] **Step 2: Add the permanent DPD matrix test**

```python
    @utx.skipIfMissingFeatures(["DPD"])
    def test_dpd_pair_visibility_matrix(self):
        system = self.system
        system.box_l = [10, 10, 10]
        grids = {"y": [1, self.n_nodes, 1],   # split along normal
                 "x": [self.n_nodes, 1, 1],   # split along shear
                 "z": [1, 1, self.n_nodes]}   # split along neutral
        for verlet in (False, True):
            for split, grid in grids.items():
                system.cell_system.node_grid = grid
                system.cell_system.set_regular_decomposition(
                    use_verlet_lists=verlet, fully_connected_boundary=None)
                self.run_dpd_pair_visibility("x", "y")
```

- [ ] **Step 3: Run on 1, 2, 4 ranks**

Run:
```bash
cd build && make -j8
ctest -R "lees_edwards$" --output-on-failure           # default rank count
```
Also run the file directly under specific rank counts:
```bash
mpiexec -n 1 build/pypresso testsuite/python/lees_edwards.py
mpiexec -n 2 build/pypresso testsuite/python/lees_edwards.py
mpiexec -n 4 build/pypresso testsuite/python/lees_edwards.py
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
mpiexec -n 2 $B $V dump --tag st_x2 --node_grid 2 1 1 --steps 2000 --out /tmp
$B $V compare --a /tmp/st_1rank.npz --b /tmp/st_y2.npz
$B $V compare --a /tmp/st_1rank.npz --b /tmp/st_x2.npz
```
Expected: trace pair sets identical; shear stress within `--stress_tol`. (Split-along-shear `st_x2` is the new capability; it has no fully_connected counterpart, so stress + trace are its proof.)

- [ ] **Step 2: Record results** in the task log.

---

## Phase 5: Retire `fully_connected_boundary`

### Task 5.1: Make the dynamic path the default and deprecate the attribute

**Files:**
- Modify: `src/script_interface/cell_system/CellSystem.cpp` (lines ~310-322: `fully_connected_boundary` handling) — emit a deprecation warning; keep it functional as the reference until removal is approved.
- Modify: `doc/` and `CHANGELOG.md` — document that LE no longer needs `fully_connected_boundary` and that splitting along the shear direction is now supported.
- Modify: `testsuite/python/regular_decomposition.py` — update any test that assumes `fully_connected` is required for LE.

**Interfaces:**
- Consumes: nothing new.
- Produces: user-visible deprecation of `fully_connected_boundary`; the dynamic path is automatic when LE is active.

- [ ] **Step 1: Decide removal vs deprecation** — this is a user gate. Present to the user: (a) hard-remove the `fully_connected_boundary` API now, or (b) keep it as a deprecated alias for one release. Do not proceed until answered.

- [ ] **Step 2: Implement the chosen option**

For (b): in `CellSystem.cpp`, when `fully_connected_boundary` is passed, log a deprecation notice via the standard warning mechanism (search `grep -rn "deprecat" src/script_interface` for the pattern in use) and keep current behavior.

- [ ] **Step 3: Update docs and changelog** with the concrete sentences: LE shear now uses an offset-driven dynamic halo; `fully_connected_boundary` is no longer required; the MPI grid may split along the shear direction.

- [ ] **Step 4: Run the full affected test set**

Run:
```bash
cd build && make -j8
ctest -R "lees_edwards$|regular_decomposition|random_pairs|hybrid_decomposition" --output-on-failure
```
Expected: all PASS.

- [ ] **Step 5: Format and commit**

```bash
maintainer/format/clang-format.sh <changed .cpp>
maintainer/format/autopep8.sh <changed .py>
git -C ... add -A
git -C ... commit -m "Deprecate fully_connected_boundary; dynamic Lees-Edwards halo is default"
```

---

## Self-review notes

- **Spec coverage:** halo sourcing (Task 3.4), narrow stencil (3.2), rebuild-on-resort (3.3), skin-sized window (`le_reach` uses `max_range`, which includes skin), LE-driven activation (3.1 gate), dpd.py (1.1), verification of all three criteria (4.1 bitwise, 4.2 trace, 4.3 stress), split-along-shear (3.4, 4.2, 4.3), retirement of `fully_connected` (Phase 5). All spec sections map to tasks.
- **Bitwise-only-where-valid:** Task 4.1 compares dynamic vs fully_connected only for `node_grid[sd]==1` grids (serial, y-split), matching the spec's constraint. Split-along-shear is proven by trace + stress only (4.2, 4.3).
- **Coexistence:** the dynamic path is gated off when `fully_connected_boundary()` is set (Global Constraints + 3.1), so the reference and the new path coexist through Phase 4.
- **Signature consistency:** `le_shear()`/`LeShear` (3.1) is consumed by 3.2, 3.3, 3.4; `rebuild_topology()`/`m_le_shift_at_last_build` (3.3) are used only within the class; `run_dpd_pair_visibility` (2.2) is consumed by 3.2, 3.4, 4.2; `dump`/`compare` (2.1) by 4.1, 4.3.
- **Open risk flagged in-plan:** the exact `ESPRESSO_MYCONFIG` cmake variable name (0.1 Step 1) and the presence of a C++ unit-test target (3.1 Step 1) are verified at execution time with the given `grep`/`cmake -L` commands rather than assumed.
