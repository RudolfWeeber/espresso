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
import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.electrostatics


@utx.skipIfMissingFeatures(["P3M"])
class P3MChargeAssignment(ut.TestCase):

    """
    Check that the per-thread charge assignment buffers of P3M are fully
    cleared between force calculations: the forces of a solver that has
    already assigned charges for other particle configurations must be
    bitwise identical to the forces of a freshly created solver.

    """

    system = espressomd.System(box_l=[10., 11., 12.])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    rng = np.random.default_rng(seed=42)
    n_part = 300

    def tearDown(self):
        self.system.part.clear()
        self.system.electrostatics.clear()
        self.system.cell_system.set_regular_decomposition()

    def make_solver(self):
        return espressomd.electrostatics.P3M(
            prefactor=1., accuracy=1e-3, mesh=[20, 22, 24], cao=7, r_cut=2.,
            alpha=1.5, tune=False, check_neutrality=False)

    def check_against_fresh_solver(self, partcls):
        self.system.integrator.run(0, recalc_forces=True)
        forces = np.copy(partcls.f)
        self.system.electrostatics.solver = self.make_solver()
        self.system.integrator.run(0, recalc_forces=True)
        np.testing.assert_array_equal(np.copy(partcls.f), forces)

    def check_configurations(self):
        box_l = np.copy(self.system.box_l)
        partcls = self.system.part.add(
            pos=self.rng.random((self.n_part, 3)) * box_l,
            q=self.rng.choice([-1., 1.], self.n_part),
            type=self.rng.integers(0, 2, self.n_part))
        self.system.electrostatics.solver = self.make_solver()
        positions = [
            self.rng.random((self.n_part, 3)) * box_l,
            # all particles in one corner of the box
            self.rng.random((self.n_part, 3)) * 2.,
            self.rng.random((self.n_part, 3)) * box_l,
        ]
        for pos in positions:
            partcls.pos = pos
            self.check_against_fresh_solver(partcls)
        # fewer particles than threads, then none
        self.system.part.clear()
        partcls = self.system.part.add(pos=[[1., 1., 1.], [5., 5., 5.]],
                                       q=[1., -1.])
        self.check_against_fresh_solver(partcls)
        self.system.part.clear()
        self.system.integrator.run(0, recalc_forces=True)

    def test_regular_decomposition(self):
        self.system.cell_system.set_regular_decomposition()
        self.check_configurations()

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] > 1,
               "P3M only supports the hybrid decomposition on one MPI rank")
    def test_hybrid_decomposition(self):
        # particles of the N-square part are not ordered by position
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=2.)
        self.check_configurations()


if __name__ == "__main__":
    ut.main()
