#
# Copyright (C) 2020 The ESPResSo project
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
import espressomd
import unittest as ut
import unittest_decorators as utx


class Integrator_test(ut.TestCase):

    system = espressomd.System(box_l=[1., 1., 1.])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    msg = 'Encountered errors during integrate: ERROR: '

    def setUp(self):
        self.system.part.add(pos=(0, 0, 0))
        self.system.integrator.set_vv()
        self.system.periodicity = 3 * [True]

    def tearDown(self):
        self.system.thermostat.turn_off()
        self.system.part.clear()

    def test_vv_integrator(self):
        self.system.thermostat.set_brownian(kT=1.0, gamma=1.0, seed=42)
        self.system.integrator.set_vv()
        with self.assertRaisesRegex(Exception, self.msg + 'The VV integrator is incompatible with the currently active combination of thermostats'):
            self.system.integrator.run(0)

    def test_brownian_integrator(self):
        self.system.integrator.set_brownian_dynamics()
        with self.assertRaisesRegex(Exception, self.msg + 'The BD integrator requires the BD thermostat'):
            self.system.integrator.run(0)

    @utx.skipIfMissingFeatures("NPT")
    def test_npt_integrator(self):
        self.system.thermostat.set_brownian(kT=1.0, gamma=1.0, seed=42)
        self.system.integrator.set_isotropic_npt(ext_pressure=1.0, piston=1.0)
        with self.assertRaisesRegex(Exception, self.msg + 'The NpT integrator requires the NpT thermostat'):
            self.system.integrator.run(0)

    @utx.skipIfMissingFeatures("STOKESIAN_DYNAMICS")
    def test_stokesian_integrator(self):
        self.system.periodicity = 3 * [False]
        self.system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
        self.system.integrator.set_stokesian_dynamics(
            viscosity=1.0, radii={0: 1.0})
        with self.assertRaisesRegex(Exception, self.msg + 'The SD integrator requires the SD thermostat'):
            self.system.integrator.run(0)

    def test_steepest_descent_integrator(self):
        with self.assertRaisesRegex(ValueError, r"The following keys have to be given as keyword arguments: "
                                                r"\['f_max', 'gamma', 'max_displacement'\], got "
                                                r"\['f_max', 'gamma', 'max_d'\] \(missing \['max_displacement'\]\)"):
            self.system.integrator.set_steepest_descent(
                f_max=0, gamma=0.1, max_d=5)
        self.system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
        self.system.integrator.set_steepest_descent(
            f_max=0, gamma=0.1, max_displacement=0.1)
        with self.assertRaisesRegex(Exception, self.msg + 'The steepest descent integrator is incompatible with thermostats'):
            self.system.integrator.run(0)


if __name__ == "__main__":
    ut.main()
