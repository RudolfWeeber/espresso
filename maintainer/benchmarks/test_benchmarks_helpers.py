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


def test_save_and_load_state_appends_npz_suffix():
    meta = {"n": 1}
    pos = np.zeros((2, 3))
    with tempfile.TemporaryDirectory() as d:
        stem = str(pathlib.Path(d) / "state")  # no .npz
        benchmarks.save_state(stem, meta, pos=pos)
        assert pathlib.Path(stem + ".npz").exists()
        meta_out, handle = benchmarks.load_state(stem)  # also no .npz
        assert meta_out == meta
        np.testing.assert_array_equal(handle["pos"], pos)
