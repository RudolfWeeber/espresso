#
# Copyright (C) 2013-2022 The ESPResSo project
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
import numpy as np
import unittest as ut


class SymplecticEuler(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()
        self.system.non_bonded_inter.reset()

    def test_integrator_activation(self):
        """Test that symplectic Euler integrator can be activated."""
        # Set symplectic Euler integrator
        self.system.integrator.set_symplectic_euler()
        
        # Should not raise any exception
        p = self.system.part.add(pos=[0, 0, 0], v=[1, 0, 0], mass=1.0)
        self.system.integrator.run(1)
        
        # Particle should have moved
        self.assertNotEqual(p.pos[0], 0.0)

    def test_free_particle_propagation(self):
        """Test free particle propagation with symplectic Euler."""
        # Add a free particle
        p = self.system.part.add(pos=[0, 0, 0], v=[1, 0, 0], mass=1.0)
        
        self.system.integrator.set_symplectic_euler()
        
        # Expected position after time dt with symplectic Euler:
        # v_new = v_old + dt * F / m (F = 0 for free particle, so v_new = v_old)
        # x_new = x_old + dt * v_new = x_old + dt * v_old
        dt = self.system.time_step
        expected_pos = dt * 1.0  # initial velocity is [1, 0, 0]
        
        self.system.integrator.run(1)
        
        # Check position
        self.assertAlmostEqual(p.pos[0], expected_pos, places=10)
        self.assertAlmostEqual(p.pos[1], 0.0, places=10)
        self.assertAlmostEqual(p.pos[2], 0.0, places=10)
        
        # Check velocity (should remain unchanged for free particle)
        self.assertAlmostEqual(p.v[0], 1.0, places=10)
        self.assertAlmostEqual(p.v[1], 0.0, places=10)
        self.assertAlmostEqual(p.v[2], 0.0, places=10)

    def test_harmonic_oscillator_energy_conservation(self):
        """Test energy conservation for harmonic oscillator with symplectic Euler."""
        # Set up harmonic oscillator with LJ potential acting as spring
        # Use soft parameters to get harmonic behavior
        p1 = self.system.part.add(pos=[4.9, 5.0, 5.0], v=[0.1, 0, 0], mass=1.0)
        p2 = self.system.part.add(pos=[5.1, 5.0, 5.0], v=[-0.1, 0, 0], mass=1.0)
        
        # Use WCA potential (repulsive part of LJ) as harmonic spring
        self.system.non_bonded_inter[0, 0].wca.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2**(1/6))
        
        self.system.integrator.set_symplectic_euler()
        
        # Calculate initial energy
        def total_energy():
            kinetic = sum(0.5 * p.mass * np.sum(p.v**2) for p in self.system.part)
            potential = self.system.analysis.energy()["total"]
            return kinetic + potential
        
        initial_energy = total_energy()
        
        # Run simulation
        energies = []
        for _ in range(100):
            self.system.integrator.run(1)
            energies.append(total_energy())
        
        # Symplectic methods should conserve energy better than non-symplectic ones
        # Energy should not drift systematically (small oscillations are OK)
        energy_drift = np.abs(energies[-1] - initial_energy)
        max_energy_variation = max(energies) - min(energies)
        
        # Energy conservation should be good (within numerical precision)
        self.assertLess(energy_drift, 0.1 * initial_energy)
        self.assertLess(max_energy_variation, 0.2 * initial_energy)

    def test_velocity_update_order(self):
        """Test that velocity is updated before position in symplectic Euler."""
        # Add particle with force (use external force for simplicity)
        p = self.system.part.add(pos=[0, 0, 0], v=[0, 0, 0], mass=1.0)
        p.ext_force = [1.0, 0, 0]  # Constant external force
        
        self.system.integrator.set_symplectic_euler()
        
        # In symplectic Euler:
        # v(n+1) = v(n) + dt * F(n) / m
        # x(n+1) = x(n) + dt * v(n+1)
        
        dt = self.system.time_step
        expected_v = dt * 1.0 / 1.0  # dt * F / m
        expected_x = dt * expected_v  # dt * v(n+1)
        
        self.system.integrator.run(1)
        
        # Check that velocity was updated first, then position
        self.assertAlmostEqual(p.v[0], expected_v, places=10)
        self.assertAlmostEqual(p.pos[0], expected_x, places=10)

    def test_comparison_with_velocity_verlet(self):
        """Compare symplectic Euler with velocity Verlet for a simple system."""
        # Set up identical initial conditions
        pos_init = [1.0, 1.0, 1.0]
        v_init = [0.1, 0.2, 0.3]
        mass = 2.0
        
        # Test with symplectic Euler
        self.system.part.clear()
        p_se = self.system.part.add(pos=pos_init, v=v_init, mass=mass)
        p_se.ext_force = [0.5, -0.3, 0.1]  # Some external force
        self.system.integrator.set_symplectic_euler()
        self.system.integrator.run(10)
        pos_se = np.copy(p_se.pos)
        v_se = np.copy(p_se.v)
        
        # Test with velocity Verlet
        self.system.part.clear()
        p_vv = self.system.part.add(pos=pos_init, v=v_init, mass=mass)
        p_vv.ext_force = [0.5, -0.3, 0.1]  # Same external force
        self.system.integrator.set_vv()
        self.system.integrator.run(10)
        pos_vv = np.copy(p_vv.pos)
        v_vv = np.copy(p_vv.v)
        
        # Results should be different (different integration schemes)
        self.assertFalse(np.allclose(pos_se, pos_vv))
        self.assertFalse(np.allclose(v_se, v_vv))

    def test_multiple_particles(self):
        """Test symplectic Euler with multiple particles."""
        # Add several particles
        particles = []
        for i in range(5):
            p = self.system.part.add(
                pos=[i, 0, 0], 
                v=[0.1 * i, 0.2 * i, 0], 
                mass=1.0 + 0.1 * i
            )
            particles.append(p)
        
        self.system.integrator.set_symplectic_euler()
        
        # Record initial state
        initial_positions = [np.copy(p.pos) for p in particles]
        initial_velocities = [np.copy(p.v) for p in particles]
        
        # Run integration
        self.system.integrator.run(5)
        
        # Check that all particles moved according to their initial velocities
        for i, p in enumerate(particles):
            # For free particles, position should change according to velocity
            self.assertFalse(np.allclose(p.pos, initial_positions[i]))
            # Velocity should remain the same for free particles
            self.assertTrue(np.allclose(p.v, initial_velocities[i], atol=1e-10))


if __name__ == "__main__":
    ut.main()