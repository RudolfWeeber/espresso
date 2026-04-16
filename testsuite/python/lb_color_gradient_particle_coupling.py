#
# Copyright (C) 2024-2026 The ESPResSo project
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
Tests for particle-fluid friction coupling in two-component color gradient LB.

Tests:
  - Drag force on particle matches gamma * (u_fluid - v_particle)
  - Mass conservation is maintained with particle coupling active
  - Simulation runs without crash for multiple steps with coupling
"""

import unittest as ut
import numpy as np

import espressomd
import espressomd.lb


DOMAIN_SIZE = 12
AGRID = 1.0
TAU = 1.0
RHO_0 = 1.0
EPSILON = 1e-6
VISCOSITY = 1.0 / 6.0
RADIUS = 4.0
SMOOTHING_WIDTH = 2.0
GAMMA = 0.5


def tanh_interpolation(distances, radius, smoothing_width,
                       rho_outer, rho_inner):
    """Smooth tanh interface profile between two density values."""
    return (rho_outer - rho_inner) / 2.0 * (
        np.tanh((distances - radius) / smoothing_width * 2.64665) - 1.0
    ) + rho_outer


def droplet_densities(grid_size, radius, smoothing_width, rho_0, epsilon):
    """Generate initial density fields for a spherical droplet."""
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


class ColorGradientParticleCouplingTest(ut.TestCase):
    """Test particle-fluid friction coupling for two-component LB."""

    system = espressomd.System(box_l=[DOMAIN_SIZE] * 3)
    system.time_step = TAU
    system.cell_system.skin = 0.4
    system.periodicity = [True, True, True]

    def tearDown(self):
        self.system.lb = None
        self.system.thermostat.turn_off()
        self.system.part.clear()

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

    def _test_drag_force_at_pos(self, pos, label):
        """Drag force on a particle in quiescent two-component fluid
        should be gamma * (u_fluid - v_particle) = -gamma * v_particle."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        v_particle = np.array([0.01, -0.02, 0.015])
        p = self.system.part.add(pos=pos, v=v_particle)

        # Enable LB thermostat (kT=0, no noise)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)

        # Run one step to trigger the coupling
        self.system.integrator.run(1)

        # Expected drag force: gamma * (u_fluid - v_particle)
        # Fluid is at rest, so F = -gamma * v_particle
        expected_force = -GAMMA * v_particle

        np.testing.assert_allclose(
            np.copy(p.f), expected_force, atol=1e-10,
            err_msg=f"Drag force mismatch for particle at {label}")

    def test_drag_force_inside_droplet(self):
        """Drag force at droplet center (component b dominant)."""
        center = DOMAIN_SIZE / 2.0
        self._test_drag_force_at_pos(
            [center + 0.5] * 3, "droplet center")

    def test_drag_force_outside_droplet(self):
        """Drag force outside droplet (component a dominant)."""
        self._test_drag_force_at_pos(
            [1.5, 1.5, 1.5], "outside droplet")

    def test_drag_force_at_interface(self):
        """Drag force at the droplet interface (mixed densities)."""
        center = DOMAIN_SIZE / 2.0
        self._test_drag_force_at_pos(
            [center + RADIUS, center, center], "interface")

    def test_mass_conservation_with_particle(self):
        """Total mass of each component should be conserved even with
        particle coupling active."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        # Place particle inside the droplet
        center = [DOMAIN_SIZE / 2.0] * 3
        p = self.system.part.add(pos=center, v=[0.01, 0.0, 0.0])
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)

        # Measure initial masses
        densities = np.copy(lbf[:, :, :].density)
        mass_a_init = np.sum(densities[:, :, :, 0])
        mass_b_init = np.sum(densities[:, :, :, 1])

        # Run with coupling
        self.system.integrator.run(20)

        # Measure final masses
        densities = np.copy(lbf[:, :, :].density)
        mass_a_final = np.sum(densities[:, :, :, 0])
        mass_b_final = np.sum(densities[:, :, :, 1])

        self.assertAlmostEqual(mass_a_init, mass_a_final, places=8,
                               msg="rho_a not conserved with particle coupling")
        self.assertAlmostEqual(mass_b_init, mass_b_final, places=8,
                               msg="rho_b not conserved with particle coupling")

    def test_run_with_coupling(self):
        """Two-component LB with particle coupling should run
        without crash for multiple steps."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        # Add particle and enable coupling
        center = [DOMAIN_SIZE / 2.0] * 3
        p = self.system.part.add(pos=center, v=[0.005, -0.005, 0.003])
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)

        # Run 50 steps — should not crash
        self.system.integrator.run(50)

        # Particle should still have finite position and velocity
        self.assertTrue(np.all(np.isfinite(np.copy(p.pos))))
        self.assertTrue(np.all(np.isfinite(np.copy(p.v))))

        # Fluid velocities should be finite
        vel = np.copy(lbf[:, :, :].velocity)
        self.assertTrue(np.all(np.isfinite(vel)))


if __name__ == "__main__":
    ut.main()
