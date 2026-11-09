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
import espressomd
import espressomd.interactions
import numpy as np


@utx.skipIfMissingFeatures(["LENNARD_JONES"])
class GhostExchangeThreads(ut.TestCase):
    """
    Forces with the regular decomposition against the N-square cell system,
    in a system large enough that, with several MPI ranks and at least four
    OpenMP threads, the ghost communication runs its threaded loops: every
    message and the self-image copies (directions with a single MPI rank)
    carry several hundred particles. With fewer threads, the serial loops
    are tested.
    """
    system = espressomd.System(box_l=3 * [1.])

    def setUp(self):
        n_part = 8000
        box_l = (n_part / 0.5)**(1. / 3.)
        self.system.box_l = 3 * [box_l]
        self.system.time_step = 0.005
        self.system.cell_system.skin = 0.4
        self.system.cell_system.set_regular_decomposition(
            use_verlet_lists=False)
        rng = np.random.default_rng(seed=42)
        m = int(np.ceil(n_part**(1. / 3.)))
        lattice = (np.indices((m, m, m)).reshape(3, -1).T[:n_part] + 0.5)
        pos = lattice * box_l / m + rng.uniform(-0.1, 0.1, (n_part, 3))
        self.partcls = self.system.part.add(pos=pos)
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1., sigma=1., cutoff=2.5, shift="auto")
        harmonic = espressomd.interactions.HarmonicBond(k=10., r_0=1.2)
        self.system.bonded_inter.add(harmonic)
        for i in range(0, n_part - 1, 2):
            self.system.part.by_id(i).add_bond((harmonic, i + 1))

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.non_bonded_inter[0, 0].lennard_jones.deactivate()
        self.system.bonded_inter.clear()

    def forces(self):
        self.system.integrator.run(0, recalc_forces=True)
        return np.copy(self.partcls.f)

    def test_forces_against_nsquare(self):
        system = self.system
        # move the particles: resorts, position pushes and force reductions
        system.thermostat.set_langevin(kT=1., gamma=1., seed=7)
        system.integrator.run(50)
        system.thermostat.turn_off()
        f_regular = self.forces()
        system.cell_system.set_n_square(use_verlet_lists=False)
        f_nsquare = self.forces()
        np.testing.assert_allclose(
            f_regular, f_nsquare, rtol=0.,
            atol=1e-10 * np.max(np.abs(f_nsquare)))


if __name__ == "__main__":
    ut.main()
