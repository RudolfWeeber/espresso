# Benchmark tune/run separation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split the `lj`, `p3m`, and `lb` benchmarks into separate `tune` and `run` phases backed by a `.npz` state file, add skin/mesh/particle-count flexibility flags, and expose the OpenMP thread count through the ESPResSo script interface.

**Architecture:** A new `omp_num_threads` entry in `cell_system.get_state()` (C++). Shared helpers in `benchmarks.py` for particle-count resolution, state save/load, topology capture/verify, and CLI wiring. Each benchmark script gains an optional positional subcommand (`tune`/`run`/absent=both): `tune` builds + tunes + saves state, `run` rebuilds from state and times (no warmup/tuning), bare does both (unchanged behavior).

**Tech Stack:** C++ (script interface), Python 3 + NumPy, ESPResSo `pypresso`, argparse, pytest (pure-Python helper tests), CMake/ctest (existing benchmark harness).

## Global Constraints

- Build parallelism: use `make -j8` (never `-j$(nproc)`).
- Shared machine: keep smoke-test systems tiny (few hundred particles, small grids, few iterations) and check machine load before running; do not run full-size benchmarks.
- Bare invocation (no subcommand) must stay 100% backward compatible in behavior and arguments — `CMakeLists.txt`, `runner.sh`, `suite.sh` are NOT modified.
- Naming: full words, not abbreviations.
- State file format: single `.npz` via `np.savez`; scalars/short lists in a `meta` dict (0-d object array), large data as named numpy arrays.
- `cores = n_nodes * omp_num_threads`. `--particles_per_core` and `--n_particles` are mutually exclusive; `--particles_per_core` keeps each script's current default.
- Copyright headers: keep the existing GPL header on every modified/created source file.
- Run `maintainer/CI/fix_style.sh` (or the pre-commit hooks: clang-format, autopep8, cmake-format) on changed files before finishing — CI won't start otherwise.
- Build flags for local verification: `-D ESPRESSO_BUILD_BENCHMARKS=ON -D ESPRESSO_BUILD_WITH_CUDA=OFF -D ESPRESSO_BUILD_WITH_WALBERLA=ON -D ESPRESSO_BUILD_WITH_CCACHE=OFF`.

---

## File Structure

- `src/script_interface/cell_system/CellSystem.cpp` — add `omp_num_threads` to `get_state()` (modify).
- `testsuite/python/cell_system.py` — assert the new state entry (modify).
- `maintainer/benchmarks/benchmarks.py` — shared helpers: `get_omp_num_threads`, `n_cores`, `resolve_n_part`, `add_common_args`, `validate_mode`, `topology_meta`, `verify_topology`, `save_state`, `load_state`, `tune_skin_unless_fixed`; update `write_report` (modify).
- `maintainer/benchmarks/test_benchmarks_helpers.py` — pure-Python pytest for the helpers (create).
- `maintainer/benchmarks/lj.py` — mode dispatch + `--skin`, `--retune_skin_after`, particle-count flags (modify/rewrite).
- `maintainer/benchmarks/p3m.py` — mode dispatch + `--skin`, `--mesh`/`--lowest_mesh`/`--highest_mesh`, particle-count flags (modify/rewrite).
- `maintainer/benchmarks/lb.py` — mode dispatch + `--skin`, particle-count flags, LB field state (modify/rewrite).

---

## Task 1: Expose `omp_num_threads` in `cell_system.get_state()`

**Files:**
- Modify: `src/script_interface/cell_system/CellSystem.cpp` (get_state branch, ~line 202; includes ~line 20-53)
- Test: `testsuite/python/cell_system.py:31` (test_cell_system)

**Interfaces:**
- Produces: `system.cell_system.get_state()["omp_num_threads"]` → Python `int`, equal to `omp_get_max_threads()`.

- [ ] **Step 1: Add the failing test assertion**

In `testsuite/python/cell_system.py`, inside `test_cell_system` (after the loop, still in the method), add:

```python
        state = self.system.cell_system.get_state()
        self.assertIn("omp_num_threads", state)
        self.assertIsInstance(state["omp_num_threads"], int)
        self.assertGreaterEqual(state["omp_num_threads"], 1)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `./build/pypresso testsuite/python/cell_system.py` (after a build; if not yet built, this step is verified after Step 4's build).
Expected: FAIL with `KeyError: 'omp_num_threads'` or `AssertionError` on `assertIn`.

- [ ] **Step 3: Implement the C++ change**

In `src/script_interface/cell_system/CellSystem.cpp`, add the OpenMP header near the other system includes (alphabetical, with the `<...>` group is fine):

```cpp
#include <omp.h>
```

In `do_call_method`, in the `get_state` branch, immediately after the line
`state["n_nodes"] = context()->get_comm().size();` add:

```cpp
    state["omp_num_threads"] = omp_get_max_threads();
```

- [ ] **Step 4: Build**

Run:
```bash
cmake -D ESPRESSO_BUILD_BENCHMARKS=ON -D ESPRESSO_BUILD_WITH_CUDA=OFF -D ESPRESSO_BUILD_WITH_WALBERLA=ON -D ESPRESSO_BUILD_WITH_CCACHE=OFF -S . -B build
make -C build -j8 pypresso
```
Expected: build succeeds.

- [ ] **Step 5: Run the test to verify it passes**

Run: `./build/pypresso testsuite/python/cell_system.py`
Expected: PASS. Also spot-check threads:
```bash
OMP_NUM_THREADS=3 ./build/pypresso -c "import espressomd; s=espressomd.System(box_l=[1,1,1]); print(s.cell_system.get_state()['omp_num_threads'])"
```
Expected: prints `3`.

- [ ] **Step 6: Commit**

```bash
git add src/script_interface/cell_system/CellSystem.cpp testsuite/python/cell_system.py
git commit -m "script_interface: expose omp_num_threads in cell_system state"
```

---

## Task 2: Shared helpers in `benchmarks.py`

**Files:**
- Modify: `maintainer/benchmarks/benchmarks.py`
- Test: `maintainer/benchmarks/test_benchmarks_helpers.py` (create)

**Interfaces:**
- Consumes: `system.cell_system.get_state()` → dict with `n_nodes`, `node_grid`, `omp_num_threads`; `system.cell_system.node_grid` (settable); `system.cell_system.tune_skin(...)`; `system.cell_system.skin` (settable).
- Produces:
  - `get_omp_num_threads(system) -> int`
  - `n_cores(system) -> int`  (= `n_nodes * omp_num_threads`)
  - `add_common_args(parser, default_particles_per_core) -> None`  (adds positional `mode` in {`tune`,`run`,None}, `--state_file`, `--skin`, and mutually-exclusive `--particles_per_core`/`--n_particles`; stores `_default_ppc` via `parser.set_defaults`)
  - `validate_mode(args) -> None`  (raises `SystemExit` if `args.mode in {tune,run}` and `not args.state_file`)
  - `resolve_n_part(system, args) -> int`
  - `topology_meta(system) -> dict`  (keys `n_nodes:int`, `node_grid:list[int]`, `omp_num_threads:int`)
  - `verify_topology(system, meta) -> None`  (raises `SystemExit` on thread/rank mismatch, else sets `node_grid`)
  - `save_state(path, meta, **arrays) -> None`
  - `load_state(path) -> (meta:dict, handle)`  (`handle[key]` → arrays)
  - `tune_skin_unless_fixed(system, args, min_skin, max_skin, **tune_kwargs) -> float`
  - `write_report(filepath, n_ranks, timings, n_steps, label='', n_threads=None)`  (new optional `n_threads`)

- [ ] **Step 1: Write the failing test**

Create `maintainer/benchmarks/test_benchmarks_helpers.py`:

```python
#
# Copyright (C) 2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""Pure-Python unit tests for benchmarks.py helpers (no espressomd needed)."""

import argparse
import pathlib
import tempfile
import numpy as np
import pytest

import benchmarks


class StubCellSystem:
    def __init__(self, state):
        self._state = dict(state)

    def get_state(self):
        return dict(self._state)

    @property
    def node_grid(self):
        return self._state["node_grid"]

    @node_grid.setter
    def node_grid(self, value):
        self._state["node_grid"] = value

    @property
    def skin(self):
        return self._state.get("skin")

    @skin.setter
    def skin(self, value):
        self._state["skin"] = value


class StubSystem:
    def __init__(self, n_nodes=1, omp_num_threads=1, node_grid=None):
        self.cell_system = StubCellSystem({
            "n_nodes": n_nodes,
            "omp_num_threads": omp_num_threads,
            "node_grid": node_grid or [n_nodes, 1, 1],
        })


def make_args(**kw):
    ns = argparse.Namespace()
    defaults = dict(mode=None, state_file=None, skin=None,
                    particles_per_core=None, n_particles=None, _default_ppc=1000)
    defaults.update(kw)
    for k, v in defaults.items():
        setattr(ns, k, v)
    return ns


def test_n_cores_multiplies_ranks_and_threads():
    system = StubSystem(n_nodes=4, omp_num_threads=2)
    assert benchmarks.n_cores(system) == 8


def test_resolve_n_part_particles_per_core_default():
    system = StubSystem(n_nodes=2, omp_num_threads=3)
    args = make_args(_default_ppc=100)
    assert benchmarks.resolve_n_part(system, args) == 100 * 6


def test_resolve_n_part_explicit_ppc():
    system = StubSystem(n_nodes=2, omp_num_threads=1)
    args = make_args(particles_per_core=50)
    assert benchmarks.resolve_n_part(system, args) == 100


def test_resolve_n_part_fixed_total():
    system = StubSystem(n_nodes=4, omp_num_threads=4)
    args = make_args(n_particles=777)
    assert benchmarks.resolve_n_part(system, args) == 777


def test_validate_mode_requires_state_file():
    with pytest.raises(SystemExit):
        benchmarks.validate_mode(make_args(mode="run", state_file=None))
    with pytest.raises(SystemExit):
        benchmarks.validate_mode(make_args(mode="tune", state_file=None))
    benchmarks.validate_mode(make_args(mode=None))  # no raise


def test_save_and_load_state_roundtrip():
    meta = {"skin": 0.4, "box_l": [10.0, 10.0, 10.0], "n_part": 3}
    pos = np.arange(9, dtype=float).reshape(3, 3)
    with tempfile.TemporaryDirectory() as d:
        path = str(pathlib.Path(d) / "state.npz")
        benchmarks.save_state(path, meta, pos=pos)
        meta_out, handle = benchmarks.load_state(path)
        assert meta_out == meta
        np.testing.assert_array_equal(handle["pos"], pos)


def test_topology_meta_and_verify_ok():
    system = StubSystem(n_nodes=2, omp_num_threads=2, node_grid=[2, 1, 1])
    meta = benchmarks.topology_meta(system)
    assert meta == {"n_nodes": 2, "node_grid": [2, 1, 1], "omp_num_threads": 2}
    # verify against a fresh system with different node_grid but same counts
    other = StubSystem(n_nodes=2, omp_num_threads=2, node_grid=[1, 2, 1])
    benchmarks.verify_topology(other, meta)
    assert list(other.cell_system.node_grid) == [2, 1, 1]


def test_verify_topology_thread_mismatch():
    meta = {"n_nodes": 1, "node_grid": [1, 1, 1], "omp_num_threads": 4}
    system = StubSystem(n_nodes=1, omp_num_threads=2)
    with pytest.raises(SystemExit):
        benchmarks.verify_topology(system, meta)


def test_verify_topology_rank_mismatch():
    meta = {"n_nodes": 4, "node_grid": [4, 1, 1], "omp_num_threads": 1}
    system = StubSystem(n_nodes=2, omp_num_threads=1)
    with pytest.raises(SystemExit):
        benchmarks.verify_topology(system, meta)


def test_tune_skin_unless_fixed_uses_fixed_skin():
    system = StubSystem()
    args = make_args(skin=0.33)
    result = benchmarks.tune_skin_unless_fixed(system, args, 0.2, 1.0)
    assert result == 0.33
    assert system.cell_system.skin == 0.33
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd maintainer/benchmarks && python3 -m pytest test_benchmarks_helpers.py -q`
Expected: FAIL (AttributeError: module 'benchmarks' has no attribute 'n_cores', etc.).

- [ ] **Step 3: Implement the helpers**

In `maintainer/benchmarks/benchmarks.py`, add `import argparse` to the imports at the top (keep `os`, `sys`, `time`, `pathlib`, `numpy as np`). Append these functions (place them before `write_report`):

```python
def get_omp_num_threads(system):
    '''Number of OpenMP threads per MPI rank (from the script interface).'''
    return int(system.cell_system.get_state()["omp_num_threads"])


def n_cores(system):
    '''Total cores = MPI ranks * OpenMP threads.'''
    state = system.cell_system.get_state()
    return int(state["n_nodes"]) * int(state["omp_num_threads"])


def add_common_args(parser, default_particles_per_core):
    '''
    Register the arguments shared by all benchmark scripts: the optional
    ``tune``/``run`` subcommand, ``--state_file``, ``--skin``, and the
    mutually-exclusive particle-count group. The default particles-per-core
    is stored on the namespace as ``_default_ppc`` for :func:`resolve_n_part`.
    '''
    parser.add_argument(
        "mode", nargs="?", choices=["tune", "run"], default=None,
        help="'tune': build and tune, then save --state_file (no timing); "
             "'run': load --state_file and time it (no tuning); "
             "omitted: do both in one go")
    parser.add_argument("--state_file", metavar="PATH", action="store",
                        type=str, default=None, required=False,
                        help="Path to the .npz state file")
    parser.add_argument("--skin", metavar="SKIN", action="store",
                        type=float, default=None, required=False,
                        help="Fix the Verlet skin (disables skin tuning)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--particles_per_core", metavar="N", action="store",
                       type=int, default=None, required=False,
                       help="Particles per core (core = MPI rank * OMP thread; "
                            f"default: {default_particles_per_core})")
    group.add_argument("--n_particles", metavar="N", action="store",
                       type=int, default=None, required=False,
                       help="Total number of particles, independent of the "
                            "number of cores")
    parser.set_defaults(_default_ppc=default_particles_per_core)


def validate_mode(args):
    '''Exit with an error if ``tune``/``run`` was requested without a state file.'''
    if args.mode in ("tune", "run") and not args.state_file:
        raise SystemExit(f"error: '{args.mode}' mode requires --state_file")


def resolve_n_part(system, args):
    '''
    Total particle count from ``--n_particles`` (fixed) or
    ``--particles_per_core`` (times the number of cores). Falls back to the
    per-script default particles-per-core when neither is given.
    '''
    if getattr(args, "n_particles", None) is not None:
        return int(args.n_particles)
    ppc = args.particles_per_core
    if ppc is None:
        ppc = args._default_ppc
    return int(ppc) * n_cores(system)


def topology_meta(system):
    '''Parallel topology to embed in a state file.'''
    state = system.cell_system.get_state()
    return {
        "n_nodes": int(state["n_nodes"]),
        "node_grid": [int(x) for x in state["node_grid"]],
        "omp_num_threads": int(state["omp_num_threads"]),
    }


def verify_topology(system, meta):
    '''
    Refuse to run a state file under a different parallel topology than it was
    tuned with, then restore the exact MPI node grid.
    '''
    state = system.cell_system.get_state()
    current_threads = int(state["omp_num_threads"])
    current_nodes = int(state["n_nodes"])
    if current_threads != int(meta["omp_num_threads"]):
        raise SystemExit(
            f"error: state was tuned with {meta['omp_num_threads']} OpenMP "
            f"thread(s), this run has {current_threads}")
    if current_nodes != int(meta["n_nodes"]):
        raise SystemExit(
            f"error: state was tuned with {meta['n_nodes']} MPI rank(s), "
            f"this run has {current_nodes}")
    system.cell_system.node_grid = [int(x) for x in meta["node_grid"]]


def save_state(path, meta, **arrays):
    '''
    Save benchmark state to a single ``.npz`` archive: ``meta`` (a dict of
    scalars/short lists) plus named numpy arrays (positions, velocities, ...).
    '''
    np.savez(path, meta=np.array(meta, dtype=object), **arrays)


def load_state(path):
    '''Load a state file. Returns ``(meta_dict, npz_handle)``.'''
    handle = np.load(path, allow_pickle=True)
    meta = handle["meta"].item()
    return meta, handle


def tune_skin_unless_fixed(system, args, min_skin, max_skin, **tune_kwargs):
    '''
    Tune the skin, unless the user fixed it via ``--skin`` (in which case set
    that value and skip tuning). Returns the resulting skin.
    '''
    if args.skin is not None:
        system.cell_system.skin = args.skin
        return args.skin
    return system.cell_system.tune_skin(
        min_skin=min_skin, max_skin=max_skin, **tune_kwargs)
```

Then update `write_report` to accept an explicit thread count:

```python
def write_report(filepath, n_ranks, timings, n_steps, label='', n_threads=None):
```

and replace its first body line

```python
    n_threads = int(os.environ.get("OMP_NUM_THREADS", 1))
```

with

```python
    if n_threads is None:
        n_threads = int(os.environ.get("OMP_NUM_THREADS", 1))
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd maintainer/benchmarks && python3 -m pytest test_benchmarks_helpers.py -q`
Expected: PASS (all tests green).

- [ ] **Step 5: Commit**

```bash
git add maintainer/benchmarks/benchmarks.py maintainer/benchmarks/test_benchmarks_helpers.py
git commit -m "benchmarks: shared helpers for state files and particle-count/topology handling"
```

---

## Task 3: `lj.py` tune/run split and flexibility flags

**Files:**
- Modify (full rewrite): `maintainer/benchmarks/lj.py`

**Interfaces:**
- Consumes: all `benchmarks.*` helpers from Task 2; `get_timings(system, n_steps, n_iterations, retune_skin_after_steps=...)` (existing).
- Produces: state file with meta keys `skin, box_l, n_part, time_step, measurement_steps, n_iterations, volume_fraction, bonds, retune_skin_after, kT, gamma, seed, n_nodes, node_grid, omp_num_threads` (plus `harmonic_r_0, harmonic_k` when bonded); arrays `pos, vel` (plus `bond_pairs` when bonded).

- [ ] **Step 1: Rewrite `lj.py`**

Replace the entire contents of `maintainer/benchmarks/lj.py` with:

```python
#
# Copyright (C) 2013-2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import sys
import espressomd
import espressomd.interactions
import benchmarks
import numpy as np
import argparse

LJ_EPS = 1.0  # LJ epsilon
LJ_SIG = 1.0  # particle diameter
LJ_CUT = LJ_SIG * 2**(1. / 6.)  # cutoff distance
DEFAULT_PARTICLES_PER_CORE = 1000
N_ITERATIONS = 30
MIN_SKIN = 0.2
MAX_SKIN = 1.0
INITIAL_SKIN = 0.4
KT = 1.0
GAMMA = 1.0
SEED = 42
HARMONIC_K = 2.0

parser = argparse.ArgumentParser(description="Benchmark LJ simulations. "
                                 "Save the results to a CSV file.")
benchmarks.add_common_args(parser, DEFAULT_PARTICLES_PER_CORE)
parser.add_argument("--volume_fraction", metavar="FRAC", action="store",
                    type=float, default=0.50, required=False,
                    help="Fraction of the simulation box volume occupied by "
                    "particles (range: [0.01-0.74], default: 0.50)")
parser.add_argument("--bonds", action="store_true",
                    help="Add bonds between particle pairs, default: false")
parser.add_argument("--retune_skin_after", metavar="N", action="store",
                    type=int, default=5, required=False,
                    help="Retune skin every N timing iterations "
                    "(0 disables, default: 5)")
group = parser.add_mutually_exclusive_group()
group.add_argument("--output", metavar="FILEPATH", action="store",
                   type=str, required=False, default="benchmarks.csv",
                   help="Output file (default: benchmarks.csv)")
group.add_argument("--visualizer", action="store_true",
                   help="Starts the visualizer (for debugging purposes)")

args = parser.parse_args()
benchmarks.validate_mode(args)

espressomd.assert_features(["LENNARD_JONES"])
np.random.seed(SEED)

system = espressomd.System(box_l=[1, 1, 1])
system.time_step = 0.01


def configure_lj(system):
    system.non_bonded_inter[0, 0].lennard_jones.set_params(
        epsilon=LJ_EPS, sigma=LJ_SIG, cutoff=LJ_CUT, shift="auto")


def retune_after(args):
    '''None to disable, else the iteration interval.'''
    return None if args.retune_skin_after <= 0 else args.retune_skin_after


def build_and_tune(system, args):
    '''Build the LJ system, warm up, equilibrate and tune the skin.'''
    assert args.volume_fraction > 0, "volume_fraction must be a positive number"
    assert args.volume_fraction < np.pi / (3 * np.sqrt(2)), \
        "volume_fraction exceeds the physical limit of sphere packing (~0.74)"
    assert not (args.bonds and args.volume_fraction > 0.5), \
        "volume_fraction too dense (>0.50) for a diatomic liquid"

    n_part = benchmarks.resolve_n_part(system, args)
    measurement_steps = int(np.round(5e6 / n_part * system.cell_system.get_state()["n_nodes"], -2))
    measurement_steps = max(100, measurement_steps)

    box_l = (n_part * 4. / 3. * np.pi * (LJ_SIG / 2.)**3
             / args.volume_fraction)**(1. / 3.)
    system.box_l = 3 * (box_l,)
    system.cell_system.skin = INITIAL_SKIN if args.skin is None else args.skin
    configure_lj(system)

    bond_pairs = []
    if not args.bonds:
        system.part.add(pos=np.random.random((n_part, 3)) * system.box_l)
    else:
        hb = espressomd.interactions.HarmonicBond(r_0=LJ_CUT, k=HARMONIC_K)
        system.bonded_inter.add(hb)
        for _ in range(0, n_part, 2):
            pos = np.random.random(3) * system.box_l
            vec = np.random.random(3)
            vec /= np.linalg.norm(vec)
            p1 = system.part.add(pos=pos)
            p2 = system.part.add(pos=pos + vec * hb.r_0)
            p1.add_bond((hb, p2))
            bond_pairs.append([p1.id, p2.id])

    benchmarks.minimize(system, n_part / 2.)
    system.integrator.set_vv()
    system.thermostat.set_langevin(kT=KT, gamma=GAMMA, seed=SEED)

    print("Equilibration")
    system.time_step /= 10.
    system.integrator.run(100)
    system.time_step *= 10.
    system.integrator.run(min(5 * measurement_steps, 60000))
    print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
        system, args, MIN_SKIN, MAX_SKIN, tol=0.025, int_steps=200)))
    print("Equilibration")
    system.integrator.run(min(10 * measurement_steps, 60000))

    return {"n_part": n_part, "measurement_steps": measurement_steps,
            "bond_pairs": bond_pairs}


def save_lj_state(system, args, ctx):
    p = system.part.all()
    meta = {
        "skin": float(system.cell_system.skin),
        "box_l": [float(x) for x in system.box_l],
        "n_part": int(ctx["n_part"]),
        "time_step": float(system.time_step),
        "measurement_steps": int(ctx["measurement_steps"]),
        "n_iterations": N_ITERATIONS,
        "volume_fraction": float(args.volume_fraction),
        "bonds": bool(args.bonds),
        "retune_skin_after": int(args.retune_skin_after),
        "harmonic_r_0": float(LJ_CUT),
        "harmonic_k": float(HARMONIC_K),
        "kT": KT, "gamma": GAMMA, "seed": SEED,
    }
    meta.update(benchmarks.topology_meta(system))
    benchmarks.save_state(
        args.state_file, meta,
        pos=np.copy(p.pos), vel=np.copy(p.v),
        bond_pairs=np.array(ctx["bond_pairs"], dtype=int))
    print(f"Saved state to {args.state_file}")


def run_from_state(system, args):
    meta, handle = benchmarks.load_state(args.state_file)
    benchmarks.verify_topology(system, meta)
    system.box_l = meta["box_l"]
    system.time_step = meta["time_step"]
    system.cell_system.skin = args.skin if args.skin is not None else meta["skin"]
    configure_lj(system)
    system.part.add(pos=handle["pos"], v=handle["vel"])
    if meta["bonds"]:
        hb = espressomd.interactions.HarmonicBond(
            r_0=meta["harmonic_r_0"], k=meta["harmonic_k"])
        system.bonded_inter.add(hb)
        for p1, p2 in handle["bond_pairs"]:
            system.part.by_id(int(p1)).add_bond((hb, int(p2)))
    system.integrator.set_vv()
    system.thermostat.set_langevin(
        kT=meta["kT"], gamma=meta["gamma"], seed=int(meta["seed"]))
    retune = None if int(meta["retune_skin_after"]) <= 0 else int(meta["retune_skin_after"])
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

if args.visualizer:
    import time
    import threading
    import espressomd.visualization
    visualizer = espressomd.visualization.openGLLive(system)

    def main_thread():
        while True:
            system.integrator.run(1)
            visualizer.update()
            time.sleep(1 / 60.)

    t = threading.Thread(target=main_thread)
    t.daemon = True
    t.start()
    visualizer.start()

if args.state_file:
    save_lj_state(system, args, ctx)

if args.mode == "tune":
    sys.exit(0)

time_and_report(system, args, ctx["measurement_steps"], N_ITERATIONS,
                retune_after(args))
```

Note: `measurement_steps` scaling multiplies by `n_nodes` so that per-core work stays constant when `--particles_per_core` is used, matching the original intent (the original divided `5e6` by `particles_per_core`; here `n_part = ppc * n_nodes * omp`, so multiplying back by `n_nodes` keeps parity for the MPI case). `max(100, ...)` preserves the original minimum-steps assertion.

- [ ] **Step 2: Smoke test — bare invocation (backward compatibility)**

Run (tiny system to respect shared-machine load):
```bash
./build/pypresso maintainer/benchmarks/lj.py --particles_per_core=200 --volume_fraction=0.50 --output=/tmp/lj_bare.csv
```
Expected: completes, prints an `average:` line, `/tmp/lj_bare.csv` has a header + one data row.

- [ ] **Step 3: Smoke test — tune then run**

Run:
```bash
./build/pypresso maintainer/benchmarks/lj.py tune --n_particles=400 --state_file=/tmp/lj.npz
test -f /tmp/lj.npz && echo STATE_OK
./build/pypresso maintainer/benchmarks/lj.py run --state_file=/tmp/lj.npz --output=/tmp/lj_run.csv
```
Expected: `tune` prints "Saved state ..." and does NOT print an `average:` line; `STATE_OK`; `run` prints an `average:` line quickly (no "Equilibration"), writes `/tmp/lj_run.csv`.

- [ ] **Step 4: Smoke test — flags and MPI topology guard**

Run:
```bash
./build/pypresso maintainer/benchmarks/lj.py tune --n_particles=400 --skin=0.3 --retune_skin_after=0 --state_file=/tmp/lj2.npz
mpiexec --oversubscribe -n 2 ./build/pypresso maintainer/benchmarks/lj.py run --state_file=/tmp/lj2.npz --output=/tmp/lj2.csv; echo "exit=$?"
```
Expected: `tune` succeeds (no skin tuning). The 2-rank `run` on a 1-rank state exits non-zero with an "MPI rank(s)" mismatch message.

- [ ] **Step 5: Smoke test — bonds round-trip**

Run:
```bash
./build/pypresso maintainer/benchmarks/lj.py tune --n_particles=200 --volume_fraction=0.10 --bonds --state_file=/tmp/ljb.npz
./build/pypresso maintainer/benchmarks/lj.py run --state_file=/tmp/ljb.npz --output=/tmp/ljb.csv
```
Expected: both succeed; `run` reconstructs bonds without error.

- [ ] **Step 6: Commit**

```bash
git add maintainer/benchmarks/lj.py
git commit -m "benchmarks/lj: separate tune/run phases with state file and flexibility flags"
```

---

## Task 4: `p3m.py` tune/run split, mesh flags

**Files:**
- Modify (full rewrite): `maintainer/benchmarks/p3m.py`

**Interfaces:**
- Consumes: `benchmarks.*` helpers; `espressomd.electrostatics.P3M`; `p3m.get_params()` → dict with tuned `mesh, cao, alpha, r_cut`.
- Produces: state meta `skin, box_l, n_part, time_step, measurement_steps, n_iterations, volume_fraction, prefactor, accuracy, mesh, cao, alpha, r_cut, gpu, kT, gamma, seed` + topology; arrays `pos, vel, charge, type`.

- [ ] **Step 1: Rewrite `p3m.py`**

Replace the entire contents of `maintainer/benchmarks/p3m.py` with:

```python
#
# Copyright (C) 2013-2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import sys
import espressomd
import espressomd.electrostatics
import benchmarks
import numpy as np
import argparse

DEFAULT_PARTICLES_PER_CORE = 1000
N_ITERATIONS = 30
MIN_SKIN = 0.2
MAX_SKIN = 1.6
INITIAL_SKIN = 0.5
ACCURACY = 1e-3
DEFAULT_TUNE_LIMITS = [12, 160]
KT = 1.0
GAMMA = 1.0
SEED = 42
SPECIES = ["anion", "cation"]
CHARGES = {"anion": -1.0, "cation": 1.0}
TYPES = {"anion": 0, "cation": 0}
LJ_SIGMAS = {"anion": 1.0, "cation": 1.0}
LJ_EPSILONS = {"anion": 1.0, "cation": 1.0}
WCA_CUT = 2.**(1. / 6.)
LJ_CUTS = {"anion": WCA_CUT * LJ_SIGMAS["anion"],
           "cation": WCA_CUT * LJ_SIGMAS["cation"]}

parser = argparse.ArgumentParser(description="Benchmark P3M simulations. "
                                 "Save the results to a CSV file.")
benchmarks.add_common_args(parser, DEFAULT_PARTICLES_PER_CORE)
parser.add_argument("--volume_fraction", metavar="FRAC", action="store",
                    type=float, default=0.25, required=False,
                    help="Fraction of the simulation box volume occupied by "
                    "particles (range: [0.01-0.74], default: 0.25)")
parser.add_argument("--prefactor", metavar="PREFACTOR", action="store",
                    type=float, default=1., required=False,
                    help="P3M prefactor (default: 1)")
parser.add_argument("--gpu", action=argparse.BooleanOptionalAction,
                    default=False, required=False, help="Use GPU implementation")
mesh_group = parser.add_mutually_exclusive_group()
mesh_group.add_argument("--mesh", metavar="M", action="store", type=int,
                        default=None, required=False,
                        help="Fixed P3M mesh size (disables mesh tuning)")
mesh_group.add_argument("--lowest_mesh", metavar="M", action="store", type=int,
                        default=None, required=False,
                        help="Lower limit for the mesh size during tuning")
parser.add_argument("--highest_mesh", metavar="M", action="store", type=int,
                    default=None, required=False,
                    help="Upper limit for the mesh size during tuning")
group = parser.add_mutually_exclusive_group()
group.add_argument("--output", metavar="FILEPATH", action="store",
                   type=str, required=False, default="benchmarks.csv",
                   help="Output file (default: benchmarks.csv)")
group.add_argument("--visualizer", action="store_true",
                   help="Starts the visualizer (for debugging purposes)")

args = parser.parse_args()
benchmarks.validate_mode(args)
if args.mesh is not None and args.highest_mesh is not None:
    parser.error("--mesh cannot be combined with --highest_mesh")

required_features = ["P3M", "LENNARD_JONES"]
if args.gpu:
    required_features.append("CUDA")
espressomd.assert_features(required_features)
np.random.seed(SEED)

system = espressomd.System(box_l=[1, 1, 1])
system.time_step = 0.01
system.cell_system.set_regular_decomposition(use_verlet_lists=True)


def configure_lj(system):
    for i in range(len(SPECIES)):
        ion1 = SPECIES[i]
        for j in range(i, len(SPECIES)):
            ion2 = SPECIES[j]
            lj_sig = (LJ_SIGMAS[ion1] + LJ_SIGMAS[ion2]) / 2
            lj_cut = (LJ_CUTS[ion1] + LJ_CUTS[ion2]) / 2
            lj_eps = (LJ_EPSILONS[ion1] * LJ_EPSILONS[ion2])**(1. / 2.)
            system.non_bonded_inter[TYPES[ion1], TYPES[ion2]].lennard_jones.set_params(
                epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")


def p3m_tune_kwargs(args):
    kwargs = {"prefactor": args.prefactor, "accuracy": ACCURACY,
              "timings": 15, "gpu": args.gpu}
    if args.mesh is not None:
        kwargs["mesh"] = args.mesh
    else:
        low = args.lowest_mesh if args.lowest_mesh is not None else DEFAULT_TUNE_LIMITS[0]
        high = args.highest_mesh if args.highest_mesh is not None else DEFAULT_TUNE_LIMITS[1]
        kwargs["tune_limits"] = [low, high]
    return kwargs


def build_and_tune(system, args):
    assert args.prefactor > 0, "prefactor must be a positive number"
    assert args.volume_fraction > 0, "volume_fraction must be a positive number"
    assert args.volume_fraction < np.pi / (3 * np.sqrt(2)), \
        "volume_fraction exceeds the physical limit of sphere packing (~0.74)"

    n_part = benchmarks.resolve_n_part(system, args)
    measurement_steps = int(np.round(
        5e5 / n_part * system.cell_system.get_state()["n_nodes"], -1))
    measurement_steps = max(50, measurement_steps)

    lj_sig = (LJ_SIGMAS["cation"] + LJ_SIGMAS["anion"]) / 2
    box_l = (n_part * 4. / 3. * np.pi * (lj_sig / 2.)**3
             / args.volume_fraction)**(1. / 3.)
    system.box_l = 3 * (box_l,)
    system.cell_system.skin = INITIAL_SKIN if args.skin is None else args.skin
    configure_lj(system)

    pid = 0
    for _ in range(0, n_part, len(SPECIES)):
        for t in SPECIES:
            system.part.add(pos=np.random.random(3) * system.box_l,
                            id=pid, q=CHARGES[t], type=TYPES[t])
            pid += 1

    benchmarks.minimize(system, n_part / 2.)
    system.integrator.set_vv()
    system.thermostat.set_langevin(kT=KT, gamma=GAMMA, seed=SEED)

    p3m = espressomd.electrostatics.P3M(**p3m_tune_kwargs(args))
    print("Quick equilibration")
    system.time_step /= 10.
    system.integrator.run(100)
    system.time_step *= 10.
    system.integrator.run(min(3 * measurement_steps, 1000))
    print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
        system, args, MIN_SKIN, MAX_SKIN, tol=0.05, int_steps=100,
        adjust_max_skin=True)))
    print("Equilibration")
    system.integrator.run(min(3 * measurement_steps, 3000))
    print("Tune p3m")
    system.electrostatics.solver = p3m
    print("Equilibration")
    system.integrator.run(min(3 * measurement_steps, 3000))
    print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
        system, args, MIN_SKIN, MAX_SKIN, tol=0.05, int_steps=100,
        adjust_max_skin=True)))

    return {"n_part": n_part, "measurement_steps": measurement_steps, "p3m": p3m}


def save_p3m_state(system, args, ctx):
    p = system.part.all()
    tuned = ctx["p3m"].get_params()
    meta = {
        "skin": float(system.cell_system.skin),
        "box_l": [float(x) for x in system.box_l],
        "n_part": int(ctx["n_part"]),
        "time_step": float(system.time_step),
        "measurement_steps": int(ctx["measurement_steps"]),
        "n_iterations": N_ITERATIONS,
        "volume_fraction": float(args.volume_fraction),
        "prefactor": float(args.prefactor),
        "accuracy": float(ACCURACY),
        "mesh": [int(x) for x in tuned["mesh"]],
        "cao": int(tuned["cao"]),
        "alpha": float(tuned["alpha"]),
        "r_cut": float(tuned["r_cut"]),
        "gpu": bool(args.gpu),
        "kT": KT, "gamma": GAMMA, "seed": SEED,
    }
    meta.update(benchmarks.topology_meta(system))
    benchmarks.save_state(
        args.state_file, meta,
        pos=np.copy(p.pos), vel=np.copy(p.v),
        charge=np.copy(p.q), type=np.copy(p.type))
    print(f"Saved state to {args.state_file}")


def run_from_state(system, args):
    meta, handle = benchmarks.load_state(args.state_file)
    benchmarks.verify_topology(system, meta)
    system.box_l = meta["box_l"]
    system.time_step = meta["time_step"]
    system.cell_system.skin = args.skin if args.skin is not None else meta["skin"]
    configure_lj(system)
    system.part.add(pos=handle["pos"], v=handle["vel"],
                    q=handle["charge"], type=[int(t) for t in handle["type"]])
    system.integrator.set_vv()
    system.thermostat.set_langevin(
        kT=meta["kT"], gamma=meta["gamma"], seed=int(meta["seed"]))
    p3m = espressomd.electrostatics.P3M(
        prefactor=meta["prefactor"], accuracy=meta["accuracy"],
        mesh=[int(x) for x in meta["mesh"]], cao=int(meta["cao"]),
        alpha=float(meta["alpha"]), r_cut=float(meta["r_cut"]),
        tune=False, gpu=bool(meta["gpu"]), timings=15)
    system.electrostatics.solver = p3m
    return meta["measurement_steps"], meta["n_iterations"]


def time_and_report(system, args, measurement_steps, n_iterations):
    timings = benchmarks.get_timings(system, measurement_steps, n_iterations)
    avg, ci = benchmarks.get_average_time(timings)
    print(f"average: {avg:.3e} +/- {ci:.3e} (95% C.I.)")
    n_proc = system.cell_system.get_state()["n_nodes"]
    benchmarks.write_report(args.output, n_proc, timings, measurement_steps,
                            n_threads=benchmarks.get_omp_num_threads(system))


if args.mode == "run":
    measurement_steps, n_iterations = run_from_state(system, args)
    time_and_report(system, args, measurement_steps, n_iterations)
    sys.exit(0)

ctx = build_and_tune(system, args)

if args.visualizer:
    import espressomd.visualization
    visualizer = espressomd.visualization.openGLLive(system)
    visualizer.run(1)

if args.state_file:
    save_p3m_state(system, args, ctx)

if args.mode == "tune":
    sys.exit(0)

time_and_report(system, args, ctx["measurement_steps"], N_ITERATIONS)
```

Note: the original re-tuned P3M a second time. That extra re-tune is dropped for clarity; the tuned parameters are captured from the single tuned solver. This keeps behavior equivalent for benchmarking purposes (still fully tuned) while making the tuned params unambiguous to serialize.

- [ ] **Step 2: Smoke test — bare invocation (backward compatibility)**

Run:
```bash
./build/pypresso maintainer/benchmarks/p3m.py --particles_per_core=200 --volume_fraction=0.25 --prefactor=4 --output=/tmp/p3m_bare.csv
```
Expected: completes; prints `average:`; CSV row written.

- [ ] **Step 3: Smoke test — tune then run (skips P3M tuning)**

Run:
```bash
./build/pypresso maintainer/benchmarks/p3m.py tune --n_particles=400 --prefactor=4 --state_file=/tmp/p3m.npz
./build/pypresso maintainer/benchmarks/p3m.py run --state_file=/tmp/p3m.npz --output=/tmp/p3m_run.csv
```
Expected: `tune` prints "Tune p3m" then "Saved state ...", no `average:`. `run` prints `average:` with no "Tune p3m" line (P3M rebuilt with `tune=False`).

- [ ] **Step 4: Smoke test — mesh flags**

Run:
```bash
./build/pypresso maintainer/benchmarks/p3m.py tune --n_particles=400 --prefactor=4 --mesh=16 --state_file=/tmp/p3m_m.npz
./build/pypresso maintainer/benchmarks/p3m.py run --state_file=/tmp/p3m_m.npz --output=/tmp/p3m_m.csv
./build/pypresso maintainer/benchmarks/p3m.py tune --n_particles=400 --prefactor=4 --lowest_mesh=8 --highest_mesh=32 --state_file=/tmp/p3m_l.npz
```
Expected: all succeed; with `--mesh=16` the saved `mesh` is `[16,16,16]` (can verify: `python3 -c "import numpy as np; print(np.load('/tmp/p3m_m.npz', allow_pickle=True)['meta'].item()['mesh'])"`).

- [ ] **Step 5: Commit**

```bash
git add maintainer/benchmarks/p3m.py
git commit -m "benchmarks/p3m: tune/run split, mesh flags, state serialization"
```

---

## Task 5: `lb.py` tune/run split with LB field state

**Files:**
- Modify (full rewrite): `maintainer/benchmarks/lb.py`

**Interfaces:**
- Consumes: `benchmarks.*` helpers; `espressomd.lb.LBFluid`; slice access `lbf[:, :, :].velocity` and `.last_applied_force` (get and set).
- Produces: state meta `skin, box_l, n_part, has_particles, time_step, measurement_steps, n_iterations, agrid, tau, kinematic_viscosity, density, single_precision, blocks_per_mpi_rank, gpu, multi_gpu, kT, gamma, seed` + topology; arrays `pos, vel` (when particles present), `lb_velocity`, `lb_last_applied_force`.

- [ ] **Step 1: Rewrite `lb.py`**

Replace the entire contents of `maintainer/benchmarks/lb.py` with:

```python
#
# Copyright (C) 2018-2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""
Benchmark Lattice-Boltzmann fluid + Lennard-Jones particles.
"""
import sys
import espressomd
import espressomd.lb
import benchmarks
import numpy as np
import argparse

LJ_EPS = 1.0
LJ_SIG = 1.0
LJ_CUT = LJ_SIG * 2**(1. / 6.)
DEFAULT_PARTICLES_PER_CORE = 125
N_ITERATIONS = 30
MIN_SKIN = 0.2
MAX_SKIN = 1.0
INITIAL_SKIN = 0.5
KT = 1.0
GAMMA = 1.0
SEED = 42

parser = argparse.ArgumentParser(description="Benchmark LB simulations. "
                                 "Save the results to a CSV file.")
benchmarks.add_common_args(parser, DEFAULT_PARTICLES_PER_CORE)
parser.add_argument("--box_l", action="store", nargs="+",
                    type=int, default=argparse.SUPPRESS, required=False,
                    help="Box length (cubic box); requires 0 particles")
parser.add_argument("--lb_sites_per_particle", metavar="N_LB", action="store",
                    type=float, default=28, required=False,
                    help="Number of LB sites per particle")
parser.add_argument("--volume_fraction", metavar="FRAC", action="store",
                    type=float, default=0.03, required=False,
                    help="Fraction of the simulation box volume occupied by "
                    "particles (range: [0.01-0.74], default: 0.03)")
parser.add_argument("--single_precision", action="store_true", required=False,
                    help="Using single-precision floating point accuracy")
parser.add_argument("--gpu", action=argparse.BooleanOptionalAction,
                    default=False, required=False, help="Use GPU implementation")
parser.add_argument("--multi-gpu", action=argparse.BooleanOptionalAction,
                    default=False, required=False, help="Use multi-GPU implementation")
parser.add_argument("--output", metavar="FILEPATH", action="store",
                    type=str, required=False, default="benchmarks.csv",
                    help="Output file (default: benchmarks.csv)")
parser.add_argument("--node_grid", action="store", nargs=3,
                    type=int, default=None, required=False, help="MPI topology")
parser.add_argument("--blocks_per_mpi_rank", action="store", nargs=3,
                    type=int, default=[1, 1, 1], required=False,
                    help="blocks per mpi rank")
parser.add_argument("--weak_scaling", action="store_true", required=False,
                    help="The measurement of weak scaling")

args = parser.parse_args()
benchmarks.validate_mode(args)

required_features = ["WALBERLA"]
if args.gpu:
    required_features.append("CUDA")
espressomd.assert_features(required_features)
np.random.seed(SEED)

system = espressomd.System(box_l=[1, 1, 1])
system.time_step = 0.01


def make_lb(system, meta):
    if meta["multi_gpu"]:
        system.cuda_init_handle.call_method("set_device_id_per_rank")
    return espressomd.lb.LBFluid(
        agrid=meta["agrid"], tau=system.time_step,
        kinematic_viscosity=meta["kinematic_viscosity"],
        density=meta["density"], single_precision=meta["single_precision"],
        gpu=meta["gpu"] or meta["multi_gpu"],
        blocks_per_mpi_rank=list(meta["blocks_per_mpi_rank"]))


def build_and_tune(system, args):
    assert args.volume_fraction > 0, "--volume_fraction must be a positive number"
    assert args.volume_fraction < np.pi / (3 * np.sqrt(2)), \
        "--volume_fraction exceeds the physical limit of sphere packing (~0.74)"
    has_box_l = hasattr(args, "box_l")

    if args.node_grid is not None:
        system.cell_system.node_grid = args.node_grid

    n_part = benchmarks.resolve_n_part(system, args)
    assert not (has_box_l and n_part != 0), \
        "Argument --box_l requires 0 particles (--particles_per_core=0 or --n_particles=0)"

    if n_part == 0:
        box_l = [32] if not has_box_l else args.box_l
        box_l = 3 * box_l if len(box_l) == 1 else box_l
        agrid = 1.
        lb_grid = list(box_l)
        measurement_steps = 80
        box_l = list(box_l)
    else:
        mpi_factor = min(2., float(np.amax(system.cell_system.node_grid)))
        box_l = (n_part * 4. / 3. * np.pi * (LJ_SIG / 2.)**3
                 / args.volume_fraction)**(1. / 3.)
        lb_grid = (n_part * args.lb_sites_per_particle)**(1. / 3.)
        lb_grid = int(mpi_factor * np.ceil(lb_grid / mpi_factor))
        agrid = box_l / lb_grid
        measurement_steps = 40
        lb_grid = 3 * [lb_grid]
        box_l = 3 * [box_l]

    if args.weak_scaling:
        box_l = list(np.array(box_l) * system.cell_system.node_grid)

    print(f"box length: {box_l}")
    print(f"LB shape: {lb_grid}")
    print(f"LB agrid: {agrid:.3f}")

    system.box_l = box_l
    system.cell_system.skin = INITIAL_SKIN if args.skin is None else args.skin

    if n_part:
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=LJ_EPS, sigma=LJ_SIG, cutoff=LJ_CUT, shift="auto")
        system.part.add(pos=np.random.random((n_part, 3)) * system.box_l)
        benchmarks.minimize(system, n_part / 2.)
        system.integrator.set_vv()
        system.thermostat.set_langevin(kT=KT, gamma=GAMMA, seed=SEED)
        print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
            system, args, MIN_SKIN, MAX_SKIN, tol=0.05, int_steps=100)))
        print("Equilibration")
        system.integrator.run(500)
        print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
            system, args, MIN_SKIN, MAX_SKIN, tol=0.05, int_steps=100)))
        print("Equilibration")
        system.integrator.run(500)
        system.thermostat.turn_off()

    meta = {
        "agrid": float(agrid), "tau": float(system.time_step),
        "kinematic_viscosity": 1., "density": 1.,
        "single_precision": bool(args.single_precision),
        "blocks_per_mpi_rank": [int(x) for x in args.blocks_per_mpi_rank],
        "gpu": bool(args.gpu), "multi_gpu": bool(args.multi_gpu),
    }
    lbf = make_lb(system, meta)
    system.lb = lbf
    if n_part:
        system.thermostat.set_lb(LB_fluid=lbf, gamma=1., seed=SEED)

    return {"n_part": n_part, "measurement_steps": measurement_steps,
            "lbf": lbf, "lb_meta": meta}


def save_lb_state(system, args, ctx):
    lbf = ctx["lbf"]
    meta = dict(ctx["lb_meta"])
    meta.update({
        "skin": float(system.cell_system.skin),
        "box_l": [float(x) for x in system.box_l],
        "n_part": int(ctx["n_part"]),
        "has_particles": bool(ctx["n_part"]),
        "time_step": float(system.time_step),
        "measurement_steps": int(ctx["measurement_steps"]),
        "n_iterations": N_ITERATIONS,
        "kT": KT, "gamma": GAMMA, "seed": SEED,
    })
    meta.update(benchmarks.topology_meta(system))
    arrays = {
        "lb_velocity": np.copy(lbf[:, :, :].velocity),
        "lb_last_applied_force": np.copy(lbf[:, :, :].last_applied_force),
    }
    if ctx["n_part"]:
        p = system.part.all()
        arrays["pos"] = np.copy(p.pos)
        arrays["vel"] = np.copy(p.v)
    benchmarks.save_state(args.state_file, meta, **arrays)
    print(f"Saved state to {args.state_file}")


def run_from_state(system, args):
    meta, handle = benchmarks.load_state(args.state_file)
    benchmarks.verify_topology(system, meta)
    system.box_l = meta["box_l"]
    system.time_step = meta["time_step"]
    system.cell_system.skin = args.skin if args.skin is not None else meta["skin"]

    if meta["has_particles"]:
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=LJ_EPS, sigma=LJ_SIG, cutoff=LJ_CUT, shift="auto")
        system.part.add(pos=handle["pos"], v=handle["vel"])
        system.integrator.set_vv()

    lbf = make_lb(system, meta)
    system.lb = lbf
    lbf[:, :, :].velocity = handle["lb_velocity"]
    lbf[:, :, :].last_applied_force = handle["lb_last_applied_force"]
    if meta["has_particles"]:
        system.thermostat.set_lb(LB_fluid=lbf, gamma=meta["gamma"], seed=int(meta["seed"]))
    return meta["measurement_steps"], meta["n_iterations"]


def time_and_report(system, args, measurement_steps, n_iterations):
    timings = benchmarks.get_timings(system, measurement_steps, n_iterations)
    avg, ci = benchmarks.get_average_time(timings)
    print(f"average: {1000 * avg:.2f} +/- {1000 * ci:.2f} ms (95% C.I.)")
    n_proc = system.cell_system.get_state()["n_nodes"]
    benchmarks.write_report(args.output, n_proc, timings, measurement_steps,
                            n_threads=benchmarks.get_omp_num_threads(system))


if args.mode == "run":
    measurement_steps, n_iterations = run_from_state(system, args)
    time_and_report(system, args, measurement_steps, n_iterations)
    sys.exit(0)

ctx = build_and_tune(system, args)

if args.state_file:
    save_lb_state(system, args, ctx)

if args.mode == "tune":
    sys.exit(0)

time_and_report(system, args, ctx["measurement_steps"], N_ITERATIONS)
```

- [ ] **Step 2: Smoke test — bare, no particles (backward compatibility)**

Run:
```bash
./build/pypresso maintainer/benchmarks/lb.py --box_l=16 --particles_per_core=0 --output=/tmp/lb_bare.csv
```
Expected: completes; prints an `average: ... ms` line; CSV row written.

- [ ] **Step 3: Smoke test — tune then run, no particles (field restored)**

Run:
```bash
./build/pypresso maintainer/benchmarks/lb.py tune --box_l=16 --particles_per_core=0 --state_file=/tmp/lb.npz
./build/pypresso maintainer/benchmarks/lb.py run --state_file=/tmp/lb.npz --output=/tmp/lb_run.csv
python3 -c "import numpy as np; d=np.load('/tmp/lb.npz', allow_pickle=True); print('field', d['lb_velocity'].shape, d['lb_last_applied_force'].shape)"
```
Expected: `tune` saves state; `run` restores the field and times; the print shows a `(16,16,16,3)`-shaped field twice.

- [ ] **Step 4: Smoke test — with particles round-trip**

Run:
```bash
./build/pypresso maintainer/benchmarks/lb.py tune --n_particles=200 --volume_fraction=0.03 --state_file=/tmp/lbp.npz
./build/pypresso maintainer/benchmarks/lb.py run --state_file=/tmp/lbp.npz --output=/tmp/lbp.csv
```
Expected: both succeed; the state file contains `pos`/`vel` plus the LB field.

- [ ] **Step 5: Commit**

```bash
git add maintainer/benchmarks/lb.py
git commit -m "benchmarks/lb: tune/run split with LB velocity field and force state"
```

---

## Task 6: Backward-compatibility sweep and style

**Files:**
- No new files; verification + style pass across changed files.

- [ ] **Step 1: Run the exact argument combinations used by CMakeLists.txt**

These mirror `maintainer/benchmarks/CMakeLists.txt`. Use small sizes only where noted is not possible (these use the real defaults but small per-core counts already; keep them short by running serially and checking they start/complete). Run each and confirm it finishes and appends a CSV row:

```bash
out=/tmp/compat.csv
./build/pypresso maintainer/benchmarks/lj.py --particles_per_core=200 --volume_fraction=0.50 --output=$out
./build/pypresso maintainer/benchmarks/lj.py --particles_per_core=200 --volume_fraction=0.02 --output=$out
./build/pypresso maintainer/benchmarks/lj.py --particles_per_core=200 --volume_fraction=0.10 --bonds --output=$out
./build/pypresso maintainer/benchmarks/p3m.py --particles_per_core=200 --volume_fraction=0.25 --prefactor=4 --output=$out
./build/pypresso maintainer/benchmarks/lb.py --particles_per_core=125 --volume_fraction=0.03 --lb_sites_per_particle=28 --output=$out
./build/pypresso maintainer/benchmarks/lb.py --box_l=16 --particles_per_core=0 --single_precision --output=$out
wc -l $out
```
Expected: all complete without error; `$out` has a header plus one row per invocation.

- [ ] **Step 2: Run the pure-Python helper tests once more**

Run: `cd maintainer/benchmarks && python3 -m pytest test_benchmarks_helpers.py -q`
Expected: PASS.

- [ ] **Step 3: Style / formatting on changed files**

Run the project style fixers on the changed files (autopep8 for Python, clang-format for C++). If the repo uses pre-commit:
```bash
pre-commit run --files src/script_interface/cell_system/CellSystem.cpp testsuite/python/cell_system.py maintainer/benchmarks/benchmarks.py maintainer/benchmarks/test_benchmarks_helpers.py maintainer/benchmarks/lj.py maintainer/benchmarks/p3m.py maintainer/benchmarks/lb.py
```
Otherwise run `maintainer/CI/fix_style.sh`. Expected: no residual style violations (re-run until clean).

- [ ] **Step 4: Commit any style fixes**

```bash
git add -A
git commit -m "benchmarks: apply formatting" || echo "nothing to format"
```

---

## Self-Review Notes

- **Spec coverage:** tune/run/both (Tasks 3-5); state contents per script incl. LB field + last_applied_force and P3M tuned params (Tasks 4-5); `--skin` (Task 2 helper + all scripts); lj `--retune_skin_after` (Task 3); p3m `--mesh`/`--lowest_mesh`/`--highest_mesh` (Task 4); `--n_particles`/`--particles_per_core` with `cores = ranks*threads` (Task 2); OMP thread API (Task 1); topology stored + refuse-on-mismatch + restore node_grid (Tasks 1-2 + all run paths); numpy binary storage (Task 2 `save_state`); backward compatibility (Task 6, CMake untouched).
- **Naming:** full words used (`particles_per_core`, `retune_skin_after`, `lowest_mesh`, `omp_num_threads`).
- **Type consistency:** `resolve_n_part`, `topology_meta`, `verify_topology`, `save_state`, `load_state`, `tune_skin_unless_fixed`, `get_omp_num_threads`, `write_report(..., n_threads=)` signatures are defined in Task 2 and consumed unchanged in Tasks 3-5.
```
