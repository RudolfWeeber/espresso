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

import itertools
import unittest as ut
import unittest_decorators as utx
import numpy as np
import scipy.spatial

import espressomd
import espressomd.electrostatics


@utx.skipIfMissingFeatures(["HAT", "ELECTROSTATICS"])
class ShortRangePairSearch(ut.TestCase):
    """
    Check the short-range forces of the link-cell pair search (with and
    without Verlet lists, regular, N-square and hybrid decomposition)
    against a brute-force sum over all pairs. The system mixes three
    particle types with different cutoffs, inactive type pairs, and charged
    and neutral particles with a Debye-Hueckel cutoff larger than most type
    cutoffs, so every branch of the Verlet criterion decides some
    candidates. The particle density gives cells with more than 64
    particles.

    Forces are checked twice: right after the Verlet list is built, and
    after a few integration steps in which the particles move by less than
    half the skin, so that the Verlet list is reused. Pairs that were in
    the skin shell (between the cutoff and the cutoff plus skin) when the
    list was built and have moved inside the cutoff only contribute in the
    second check, so it fails if the pair search drops skin-shell pairs.
    """

    system = espressomd.System(box_l=[1., 1., 1.])
    system.time_step = 0.05
    system.cell_system.skin = 0.4

    # Hat cutoffs per type pair, None = no interaction
    hat_cutoffs = {(0, 0): 2.0, (0, 1): 2.5, (1, 1): 3.0, (0, 2): None,
                   (1, 2): 1.5, (2, 2): 1.0}
    hat_fmax = 1.3
    dh_params = {"prefactor": 0.7, "kappa": 0.8, "r_cut": 2.8}
    # interaction range 3.0 + skin 0.4
    min_local_box = 3.4
    # integration steps with the Verlet list of the first check
    n_steps = 2
    # initial speed: 0.1 displacement in n_steps without forces, well
    # below half the skin even with the force contribution
    speed = 1.0

    def tearDown(self):
        self.reset_system()

    def reset_system(self):
        self.system.part.clear()
        self.system.electrostatics.clear()
        self.system.non_bonded_inter.reset()
        self.system.periodicity = [True, True, True]
        self.system.cell_system.set_regular_decomposition()

    def setup_system(self, box_l, periodicity):
        system = self.system
        system.box_l = box_l
        system.periodicity = periodicity
        for (t1, t2), cut in self.hat_cutoffs.items():
            if cut is not None:
                system.non_bonded_inter[t1, t2].hat.set_params(
                    F_max=self.hat_fmax, cutoff=cut)
        # jittered lattice (no close contacts), random types, charges and
        # velocity directions
        rng = np.random.default_rng(42)
        spacing = 0.7
        grid = [np.arange(0.5 * spacing, length, spacing) for length in box_l]
        pos = np.array(np.meshgrid(*grid, indexing="ij")).reshape(3, -1).T
        pos += rng.uniform(-0.15, 0.15, pos.shape)
        pos = np.clip(pos, 0., np.array(box_l) * (1. - 1e-12))
        types = rng.integers(0, 3, len(pos))
        q = rng.choice([-1., 0., 1.], len(pos))
        v = rng.normal(size=pos.shape)
        v *= self.speed / np.linalg.norm(v, axis=1)[:, np.newaxis]
        # particles must not leave the box through non-periodic faces
        v[:, np.logical_not(periodicity)] = 0.
        system.part.add(pos=pos, type=types, q=q)
        system.electrostatics.solver = espressomd.electrostatics.DH(
            check_neutrality=False, **self.dh_params)
        return pos, types, q, v

    def reference_forces(self, pos, types, q, box_l, periodicity):
        # boxes of non-periodic directions are inflated, so that no pair
        # interacts through them
        box_l = np.array(box_l)
        box = np.where(periodicity, box_l, 10. * box_l)
        pos = np.where(periodicity, np.mod(pos, box_l), pos + 4.5 * box_l)
        tree = scipy.spatial.cKDTree(pos, boxsize=box)
        pairs = tree.query_pairs(3.0, output_type="ndarray")
        i, j = pairs[:, 0], pairs[:, 1]
        d = pos[i] - pos[j]
        d -= np.rint(d / box) * box
        r = np.linalg.norm(d, axis=1)
        cut_table = np.zeros((3, 3))
        for (t1, t2), cut in self.hat_cutoffs.items():
            if cut is not None:
                cut_table[t1, t2] = cut_table[t2, t1] = cut
        cut = cut_table[types[i], types[j]]
        mask = r < cut
        fac = np.zeros(len(r))
        fac[mask] = self.hat_fmax * (1. - r[mask] / cut[mask]) / r[mask]
        mask = r < self.dh_params["r_cut"]
        kappa = self.dh_params["kappa"]
        fac[mask] += (self.dh_params["prefactor"] * q[i][mask] * q[j][mask] /
                      r[mask]**3 * np.exp(-kappa * r[mask]) *
                      (1. + kappa * r[mask]))
        f_pair = fac[:, np.newaxis] * d
        n = len(pos)
        return np.array([np.bincount(i, f_pair[:, k], minlength=n) -
                         np.bincount(j, f_pair[:, k], minlength=n)
                         for k in range(3)]).T

    def check_forces(self, f_ref0, pos0, types, q, v, box_l, periodicity):
        system = self.system
        partcls = system.part.all()
        partcls.pos = pos0
        partcls.v = v
        system.integrator.run(0, recalc_forces=True)
        np.testing.assert_allclose(np.copy(partcls.f), f_ref0, rtol=0.,
                                   atol=1e-10)
        # move the particles with the Verlet list built above
        system.integrator.run(self.n_steps)
        pos = np.copy(partcls.pos)
        # no particle moved by more than half the skin: the list was reused
        self.assertLess(np.max(np.linalg.norm(pos - pos0, axis=1)),
                        0.5 * system.cell_system.skin)
        f_ref = self.reference_forces(pos, types, q, box_l, periodicity)
        np.testing.assert_allclose(np.copy(partcls.f), f_ref, rtol=0.,
                                   atol=1e-10)

    def check_box(self, box_l, periodicity):
        args = self.setup_system(box_l, periodicity)
        f_ref0 = self.reference_forces(*args[:3], box_l, periodicity)
        self.assertGreater(np.max(np.linalg.norm(f_ref0, axis=1)), 0.)
        cell_system = self.system.cell_system
        cell_systems = [
            ("regular", True, lambda: cell_system.set_regular_decomposition(
                use_verlet_lists=True)),
            ("regular", False, lambda: cell_system.set_regular_decomposition(
                use_verlet_lists=False)),
            ("n_square", True, lambda: cell_system.set_n_square(
                use_verlet_lists=True)),
            ("hybrid", True, lambda: cell_system.set_hybrid_decomposition(
                n_square_types={2}, cutoff_regular=3.0,
                use_verlet_lists=True)),
        ]
        for name, use_verlet_lists, set_cell_system in cell_systems:
            with self.subTest(box_l=box_l, periodicity=periodicity,
                              cell_system=name,
                              use_verlet_lists=use_verlet_lists):
                set_cell_system()
                self.check_forces(f_ref0, *args, box_l, periodicity)
        self.reset_system()

    def test_pair_search(self):
        node_grid = np.array(self.system.cell_system.node_grid)
        boxes = [[11., 14., 17.5],  # 3 x 4 x 5 cells on one MPI rank
                 [7., 7., 10.5],    # 2 x 2 x 3 cells, several images
                 [3.5, 7., 14.]]    # 1 x 2 x 4 cells
        n_checks = 0
        for box_l in boxes:
            if np.any(np.array(box_l) / node_grid < self.min_local_box):
                continue
            for periodicity in itertools.product((True, False), repeat=3):
                if periodicity.count(False) == 1 or all(periodicity):
                    self.check_box(box_l, list(periodicity))
                    n_checks += 1
        self.assertGreater(n_checks, 0)


if __name__ == "__main__":
    ut.main()
