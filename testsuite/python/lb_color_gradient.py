#
# Copyright (C) 2024 The ESPResSo project
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
Basic tests for two-component color gradient LB.

Tests:
  - Construction with two viscosities triggers two-component mode
  - Setting/getting per-component densities at nodes
  - Running integration steps without crash
  - Mass conservation (total rho_a and rho_b are conserved)
"""

import unittest as ut
import numpy as np

import espressomd
import espressomd.lb


DOMAIN_SIZE = 12
AGRID = 1.0
RHO_0 = 1.0
EPSILON = 1e-6
VISCOSITY = 1.0 / 6.0
RADIUS = 4.0
SMOOTHING_WIDTH = 2.0


def tanh_interpolation(distances, radius, smoothing_width,
                       rho_outer, rho_inner):
    """Smooth tanh interface profile between two density values."""
    return (rho_outer - rho_inner) / 2.0 * (
        np.tanh((distances - radius) / smoothing_width * 2.64665) - 1.0
    ) + rho_outer


def droplet_densities(grid_size, radius, smoothing_width, rho_0, epsilon):
    """
    Generate initial density fields for a spherical droplet.

    Returns:
        rho_a: solvent density (rho_0 outside, epsilon*rho_0 inside)
        rho_b: droplet density (epsilon*rho_0 outside, rho_0 inside)
    """
    x = np.arange(grid_size) + 0.5
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    center = grid_size / 2.0
    distances = np.sqrt((xx - center)**2 + (yy - center)**2
                        + (zz - center)**2)

    rho_a = tanh_interpolation(distances, radius, smoothing_width,
                               rho_0, epsilon * rho_0)
    rho_b = tanh_interpolation(distances, radius, smoothing_width,
                               epsilon * rho_0, rho_0)
    return rho_a, rho_b


class ColorGradientLBTest(ut.TestCase):
    """Test the two-component color gradient LB method."""

    system = espressomd.System(box_l=[DOMAIN_SIZE] * 3)
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.periodicity = [True, True, True]

    def tearDown(self):
        self.system.lb = None

    def _create_lbf(self):
        """Create a two-component LB fluid and attach it."""
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=self.system.time_step,
            kinematic_viscosity=[VISCOSITY, VISCOSITY])
        self.system.lb = lbf
        return lbf

    def _init_droplet(self, lbf):
        """Set per-node densities to a spherical droplet profile and
        initialize the PDFs."""
        rho_a, rho_b = droplet_densities(
            DOMAIN_SIZE, RADIUS, SMOOTHING_WIDTH, RHO_0, EPSILON)
        for x in range(DOMAIN_SIZE):
            for y in range(DOMAIN_SIZE):
                for z in range(DOMAIN_SIZE):
                    lbf[x, y, z].density = [rho_a[x, y, z], rho_b[x, y, z]]
        lbf.init_two_component()

    def test_construction(self):
        """Two viscosities should create a two-component LB."""
        lbf = self._create_lbf()
        visc = lbf.kinematic_viscosity
        self.assertEqual(len(visc), 2)
        self.assertAlmostEqual(visc[0], VISCOSITY, places=10)
        self.assertAlmostEqual(visc[1], VISCOSITY, places=10)

    def test_node_density_two_component(self):
        """Getting density on a node should return two values in CG mode."""
        lbf = self._create_lbf()
        dens = lbf[0, 0, 0].density
        self.assertEqual(len(dens), 2)

    def test_set_density(self):
        """Setting per-node density with two values should round-trip."""
        lbf = self._create_lbf()
        lbf[0, 0, 0].density = [0.7, 0.3]
        dens = lbf[0, 0, 0].density
        self.assertAlmostEqual(dens[0], 0.7, places=10)
        self.assertAlmostEqual(dens[1], 0.3, places=10)

    def test_run_steps(self):
        """Two-component LB with droplet should run without crash."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)
        self.system.integrator.run(10)

    def test_mass_conservation(self):
        """Total mass of each component should be conserved."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        # Measure initial total density per component
        rho_a_init = 0.0
        rho_b_init = 0.0
        for x in range(DOMAIN_SIZE):
            for y in range(DOMAIN_SIZE):
                for z in range(DOMAIN_SIZE):
                    d = lbf[x, y, z].density
                    rho_a_init += d[0]
                    rho_b_init += d[1]

        self.assertGreater(rho_a_init, 0.0, "rho_a should be > 0 after init")
        self.assertGreater(rho_b_init, 0.0, "rho_b should be > 0 after init")

        # Run some steps
        self.system.integrator.run(20)

        # Measure final total density per component
        rho_a_final = 0.0
        rho_b_final = 0.0
        for x in range(DOMAIN_SIZE):
            for y in range(DOMAIN_SIZE):
                for z in range(DOMAIN_SIZE):
                    d = lbf[x, y, z].density
                    rho_a_final += d[0]
                    rho_b_final += d[1]

        self.assertAlmostEqual(rho_a_init, rho_a_final, places=8,
                               msg="rho_a not conserved")
        self.assertAlmostEqual(rho_b_init, rho_b_final, places=8,
                               msg="rho_b not conserved")


if __name__ == "__main__":
    ut.main()
