#
# Copyright (C) 2013-2019 The ESPResSo project
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
import numpy as np
import itertools
import unittest as ut
import unittest_decorators as utx
import tests_common
import espressomd.electrostatics


class ElectrostaticInteractionsTests:
    MMM1D = None

    # Handle to espresso system
    system = espressomd.System(box_l=[10.0] * 3)
    system.periodicity = [0, 0, 1]
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.cell_system.set_n_square()
    system.thermostat.set_langevin(kT=0, gamma=1, seed=8)

    data = np.loadtxt(tests_common.data_path("mmm1d_data.txt"))
    p_pos = data[:, 1:4]
    p_q = data[:, 4]
    forces_target = data[:, 5:8]
    energy_target = -7.156365298205383

    allowed_error = 2e-5

    def setUp(self):
        self.system.periodicity = [0, 0, 1]
        self.system.cell_system.set_n_square()

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    def test_forces_and_energy(self):
        self.system.part.add(pos=self.p_pos, q=self.p_q)
        mmm1d = self.MMM1D(prefactor=1.0, maxPWerror=1e-20)
        self.system.actors.add(mmm1d)
        self.system.integrator.run(steps=0)
        measured_f = np.copy(self.system.part.all().f)
        np.testing.assert_allclose(measured_f, self.forces_target,
                                   atol=self.allowed_error)
        measured_el_energy = self.system.analysis.energy()["coulomb"]
        self.assertAlmostEqual(
            measured_el_energy, self.energy_target, delta=self.allowed_error,
            msg="Measured energy deviates too much from stored result")

    def check_with_analytical_result(self, prefactor, accuracy):
        p = self.system.part.by_id(0)
        f_measured = p.f
        energy_measured = self.system.analysis.energy()["total"]
        target_energy_config = -1.00242505606 * prefactor
        target_force_z_config = 0.99510759 * prefactor

        self.assertAlmostEqual(
            f_measured[0], 0, delta=self.allowed_error,
            msg="Measured force in x deviates too much from analytical result")
        self.assertAlmostEqual(
            f_measured[1], 0, delta=self.allowed_error,
            msg="Measured force in y deviates too much from analytical result")
        self.assertAlmostEqual(
            f_measured[2], target_force_z_config, delta=accuracy,
            msg="Measured force in z deviates too much from analytical result")
        self.assertAlmostEqual(
            energy_measured, target_energy_config, delta=self.allowed_error,
            msg="Measured energy deviates too much from analytical result")

    def test_with_analytical_result(self):
        self.system.part.add(pos=[0, 0, 0], q=1)
        self.system.part.add(pos=[0, 0, 1], q=-1)
        mmm1d = self.MMM1D(prefactor=1.0, maxPWerror=1e-20)
        self.system.actors.add(mmm1d)
        self.system.integrator.run(steps=0)
        self.check_with_analytical_result(prefactor=1.0, accuracy=0.0004)

    def test_bjerrum_length_change(self):
        self.system.part.add(pos=[0, 0, 0], q=1)
        self.system.part.add(pos=[0, 0, 1], q=-1)
        mmm1d = self.MMM1D(prefactor=2.0, maxPWerror=1e-20)
        self.system.actors.add(mmm1d)
        self.system.integrator.run(steps=0)
        self.check_with_analytical_result(prefactor=2.0, accuracy=0.0017)

    def test_exceptions(self):
        # check periodicity exceptions
        for periodicity in itertools.product(range(2), range(2), range(2)):
            if periodicity == (0, 0, 1):
                continue
            self.system.periodicity = periodicity
            with self.assertRaisesRegex(Exception, r"MMM1D requires periodicity \(0, 0, 1\)"):
                mmm1d = self.MMM1D(prefactor=1.0, maxPWerror=1e-2)
                self.system.actors.add(mmm1d)
            self.system.periodicity = (0, 0, 1)
            self.system.actors.clear()
        if self.MMM1D is espressomd.electrostatics.MMM1D:
            with self.assertRaisesRegex(Exception, "MMM1D requires the N-square cellsystem"):
                mmm1d = self.MMM1D(prefactor=1.0, maxPWerror=1e-2)
                self.system.cell_system.set_regular_decomposition()
                self.system.actors.add(mmm1d)
            self.system.cell_system.set_n_square()
            self.system.actors.clear()


@utx.skipIfMissingFeatures(["ELECTROSTATICS"])
class MMM1D_Test(ElectrostaticInteractionsTests, ut.TestCase):

    def setUp(self):
        self.MMM1D = espressomd.electrostatics.MMM1D
        super().setUp()


if __name__ == "__main__":
    ut.main()
