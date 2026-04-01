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
  - Density set/get for two-component mode (node and slice)
  - Population get/set for two-component mode (node and slice)
  - Setting velocity raises RuntimeError in two-component mode
  - Pressure tensor raises RuntimeError in two-component mode
  - Adding boundaries raises RuntimeError in two-component mode
  - Velocity getter returns sensible values after integration
  - init_two_component produces PDFs consistent with set densities
  - Viscosity setter works on live CG fluid
  - Bulk slice init matches per-node init
  - kT > 0 is rejected for two-component mode
  - Single viscosity does not create two-component mode
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
        """Set densities to a spherical droplet profile and
        initialize the PDFs."""
        rho_a, rho_b = droplet_densities(
            DOMAIN_SIZE, RADIUS, SMOOTHING_WIDTH, RHO_0, EPSILON)
        lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)
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

    def test_slice_density_set_get(self):
        """Bulk slice density set/get should round-trip for two components."""
        lbf = self._create_lbf()
        N = DOMAIN_SIZE
        rho_a, rho_b = droplet_densities(
            N, RADIUS, SMOOTHING_WIDTH, RHO_0, EPSILON)

        # Set via bulk slice
        lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)

        # Read back via bulk slice
        densities = np.copy(lbf[:, :, :].density)
        self.assertEqual(densities.shape, (N, N, N, 2))
        np.testing.assert_allclose(densities[:, :, :, 0], rho_a, rtol=1e-10)
        np.testing.assert_allclose(densities[:, :, :, 1], rho_b, rtol=1e-10)

        # Verify consistency with per-node access
        for x in range(0, N, N // 3):
            for y in range(0, N, N // 3):
                for z in range(0, N, N // 3):
                    d = lbf[x, y, z].density
                    self.assertAlmostEqual(d[0], rho_a[x, y, z], places=10)
                    self.assertAlmostEqual(d[1], rho_b[x, y, z], places=10)

    def test_node_population_set_get(self):
        """Node population get/set should round-trip for two components."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        # Read populations (2 * stencil_size values)
        pop = np.copy(lbf[3, 3, 3]._population)
        stencil_size = len(pop) // 2
        self.assertEqual(len(pop), 2 * stencil_size)

        # Modify and write back
        pop_modified = pop * 0.99
        lbf[3, 3, 3]._population = pop_modified
        pop_read = np.copy(lbf[3, 3, 3]._population)
        np.testing.assert_allclose(pop_read, pop_modified, rtol=1e-10)

    def test_slice_population_set_get(self):
        """Slice population get/set should round-trip for two components."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        N = DOMAIN_SIZE
        pop = np.copy(lbf[:, :, :]._population)
        stencil_size = pop.shape[-1] // 2
        self.assertEqual(pop.shape, (N, N, N, 2 * stencil_size))

        # Verify consistency with per-node access
        node_pop = np.copy(lbf[2, 3, 4]._population)
        np.testing.assert_allclose(pop[2, 3, 4], node_pop, rtol=1e-10)

    def test_set_velocity_raises(self):
        """Setting velocity should raise in two-component mode."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        with self.assertRaisesRegex(RuntimeError, "not supported for two-component"):
            lbf[0, 0, 0].velocity = [0.1, 0.0, 0.0]

        with self.assertRaisesRegex(RuntimeError, "not supported for two-component"):
            lbf[:, :, :].velocity = np.zeros((DOMAIN_SIZE,) * 3 + (3,))

    def test_pressure_tensor_raises(self):
        """Pressure tensor should raise in two-component mode."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        with self.assertRaisesRegex(RuntimeError, "not implemented for two-component"):
            _ = lbf[0, 0, 0].pressure_tensor

        with self.assertRaisesRegex(RuntimeError, "not implemented for two-component"):
            _ = lbf[:, :, :].pressure_tensor

    def test_boundary_raises(self):
        """Adding boundaries should raise in two-component mode."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        with self.assertRaisesRegex(RuntimeError, "not implemented for two-component"):
            lbf[0, 0, 0].boundary = espressomd.lb.VelocityBounceBack([0, 0, 0])

    def test_velocity_getter(self):
        """Velocity getter should return finite values after integration."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)
        self.system.integrator.run(5)

        # Node velocity
        vel = np.copy(lbf[3, 3, 3].velocity)
        self.assertEqual(len(vel), 3)
        self.assertTrue(np.all(np.isfinite(vel)))

        # Slice velocity
        vel_slice = np.copy(lbf[:, :, :].velocity)
        N = DOMAIN_SIZE
        self.assertEqual(vel_slice.shape, (N, N, N, 3))
        self.assertTrue(np.all(np.isfinite(vel_slice)))

        # Node and slice should be consistent
        np.testing.assert_allclose(vel_slice[3, 3, 3], vel, atol=1e-10)

    def test_init_two_component_consistency(self):
        """After init_two_component, densities should match what was set."""
        lbf = self._create_lbf()
        rho_a, rho_b = droplet_densities(
            DOMAIN_SIZE, RADIUS, SMOOTHING_WIDTH, RHO_0, EPSILON)
        lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)
        lbf.init_two_component()

        # Densities should still match after init
        densities = np.copy(lbf[:, :, :].density)
        np.testing.assert_allclose(densities[:, :, :, 0], rho_a, rtol=1e-5)
        np.testing.assert_allclose(densities[:, :, :, 1], rho_b, rtol=1e-5)

    def test_viscosity_setter(self):
        """Changing viscosity on a live CG fluid should work."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        new_visc = [VISCOSITY * 2, VISCOSITY * 3]
        lbf.kinematic_viscosity = new_visc
        visc = lbf.kinematic_viscosity
        self.assertEqual(len(visc), 2)
        self.assertAlmostEqual(visc[0], new_visc[0], places=10)
        self.assertAlmostEqual(visc[1], new_visc[1], places=10)

    def test_bulk_slice_init_matches_node_init(self):
        """Bulk slice density init should match per-node init."""
        lbf = self._create_lbf()
        N = DOMAIN_SIZE
        rho_a, rho_b = droplet_densities(
            N, RADIUS, SMOOTHING_WIDTH, RHO_0, EPSILON)

        # Init via bulk slice
        lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)
        lbf.init_two_component()
        densities_slice = np.copy(lbf[:, :, :].density)

        # Reset and init via per-node loop (to verify both paths)
        self.system.lb = None
        lbf2 = self._create_lbf()
        for x in range(N):
            for y in range(N):
                for z in range(N):
                    lbf2[x, y, z].density = [rho_a[x, y, z], rho_b[x, y, z]]
        lbf2.init_two_component()
        densities_node = np.copy(lbf2[:, :, :].density)

        np.testing.assert_allclose(
            densities_slice, densities_node, rtol=1e-10)

    def test_thermalized_cg_rejected(self):
        """kT > 0 should be rejected for two-component mode."""
        with self.assertRaisesRegex(ValueError, "not supported.*two-component"):
            espressomd.lb.LBFluid(
                agrid=AGRID, density=RHO_0, tau=self.system.time_step,
                kinematic_viscosity=[VISCOSITY, VISCOSITY],
                kT=1.0, seed=42)

    def test_single_viscosity_not_two_component(self):
        """Single viscosity should not create a two-component LB."""
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=self.system.time_step,
            kinematic_viscosity=VISCOSITY)
        self.system.lb = lbf
        # Single-component: density should be a scalar
        dens = lbf[0, 0, 0].density
        self.assertIsInstance(dens, float)

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
        densities = np.copy(lbf[:, :, :].density)
        rho_a_init = np.sum(densities[:, :, :, 0])
        rho_b_init = np.sum(densities[:, :, :, 1])

        self.assertGreater(rho_a_init, 0.0, "rho_a should be > 0 after init")
        self.assertGreater(rho_b_init, 0.0, "rho_b should be > 0 after init")

        # Run some steps
        self.system.integrator.run(20)

        # Measure final total density per component
        densities = np.copy(lbf[:, :, :].density)
        rho_a_final = np.sum(densities[:, :, :, 0])
        rho_b_final = np.sum(densities[:, :, :, 1])

        self.assertAlmostEqual(rho_a_init, rho_a_final, places=8,
                               msg="rho_a not conserved")
        self.assertAlmostEqual(rho_b_init, rho_b_final, places=8,
                               msg="rho_b not conserved")


if __name__ == "__main__":
    ut.main()
