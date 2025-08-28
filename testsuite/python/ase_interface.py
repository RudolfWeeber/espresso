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

import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.plugins.ase
import numpy as np
import ase


@utx.skipIfMissingFeatures("EXTERNAL_FORCES")
class ASEInterfaceComprehensiveTest(ut.TestCase):
    """Comprehensive test suite for the ASE interface focusing on update_ase() method."""

    system = espressomd.System(box_l=[10., 10., 10.])

    def setUp(self):
        """Set up system with particles having various properties."""
        self.system.time_step = 0.01
        
        # Add particles with positions, charges, masses, velocities
        self.system.part.add(pos=[1., 2., 3.], q=1.0, mass=2.0, v=[0.1, 0.2, 0.3], type=0)
        self.system.part.add(pos=[4., 5., 6.], q=-0.5, mass=1.5, v=[0.4, 0.5, 0.6], type=1)
        self.system.part.add(pos=[7., 8., 9.], q=2.0, mass=3.0, v=[0.7, 0.8, 0.9], type=0)
        
        self.type_mapping = {0: "H", 1: "O"}
        
    def tearDown(self):
        """Clean up after each test."""
        self.system.part.clear()
        
    def test_ase_interface_instantiation_basic(self):
        """Test basic ASE interface creation without optional exports."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=self.system,
            type_mapping=self.type_mapping,
            particle_slice=self.system.part.all()
        )
        
        # Check that interface is properly initialized
        self.assertIsNotNone(ase_interface.atoms)
        self.assertEqual(ase_interface.export_charges, False)
        self.assertEqual(ase_interface.export_masses, False) 
        self.assertEqual(ase_interface.export_momenta, False)
        
        # Verify atoms object has correct basic properties
        atoms = ase_interface.atoms
        self.assertEqual(len(atoms), 3)
        np.testing.assert_array_equal(atoms.get_chemical_symbols(), ["H", "O", "H"])
        np.testing.assert_allclose(atoms.positions, [[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]])

    def test_ase_interface_instantiation_with_charges(self):
        """Test ASE interface creation with charge export enabled."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=self.system,
            type_mapping=self.type_mapping,
            particle_slice=self.system.part.all(),
            export_charges=True
        )
        
        self.assertEqual(ase_interface.export_charges, True)
        # Verify charges are set in ASE atoms
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.get_initial_charges(), [1.0, -0.5, 2.0])

    def test_ase_interface_instantiation_with_masses(self):
        """Test ASE interface creation with mass export enabled."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=self.system,
            type_mapping=self.type_mapping,
            particle_slice=self.system.part.all(),
            export_masses=True
        )
        
        self.assertEqual(ase_interface.export_masses, True)
        # Verify masses are set in ASE atoms
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.get_masses(), [2.0, 1.5, 3.0])

    def test_ase_interface_instantiation_with_momenta(self):
        """Test ASE interface creation with momentum export enabled.""" 
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=self.system,
            type_mapping=self.type_mapping,
            particle_slice=self.system.part.all(),
            export_momenta=True
        )
        
        self.assertEqual(ase_interface.export_momenta, True)
        # Verify momenta are set in ASE atoms (momentum = mass * velocity)
        atoms = ase_interface.atoms
        expected_momenta = np.array([[0.1*2.0, 0.2*2.0, 0.3*2.0],
                                    [0.4*1.5, 0.5*1.5, 0.6*1.5],
                                    [0.7*3.0, 0.8*3.0, 0.9*3.0]])
        np.testing.assert_allclose(atoms.get_momenta(), expected_momenta)

    def test_ase_interface_instantiation_all_exports(self):
        """Test ASE interface creation with all exports enabled."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=self.system,
            type_mapping=self.type_mapping,
            particle_slice=self.system.part.all(),
            export_charges=True,
            export_masses=True,
            export_momenta=True
        )
        
        # Check all export flags
        self.assertEqual(ase_interface.export_charges, True)
        self.assertEqual(ase_interface.export_masses, True)
        self.assertEqual(ase_interface.export_momenta, True)
        
        # Verify all properties are set
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.get_initial_charges(), [1.0, -0.5, 2.0])
        np.testing.assert_allclose(atoms.get_masses(), [2.0, 1.5, 3.0])
        expected_momenta = np.array([[0.2, 0.4, 0.6], [0.6, 0.75, 0.9], [2.1, 2.4, 2.7]])
        np.testing.assert_allclose(atoms.get_momenta(), expected_momenta)

    def test_update_ase_basic_no_exports(self):
        """Test update_ase() with no exports enabled - only positions should update."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=self.system,
            type_mapping=self.type_mapping,
            particle_slice=self.system.part.all()
        )
        
        # Change particle positions
        self.system.part.by_id(0).pos = [10., 11., 12.]
        self.system.part.by_id(1).pos = [13., 14., 15.]
        self.system.part.by_id(2).pos = [16., 17., 18.]
        
        # Update ASE
        ase_interface.update_ase()
        
        # Check that positions are updated
        atoms = ase_interface.atoms
        expected_positions = self.system.part.all().pos
        np.testing.assert_allclose(atoms.positions, expected_positions)

    def test_update_ase_charges_no_skip(self):
        """Test update_ase() with charges enabled and no skip."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=self.system,
            type_mapping=self.type_mapping,
            particle_slice=self.system.part.all(),
            export_charges=True
        )
        
        # Change particle charges 
        self.system.part.all().q = [5.,-2.,3.5] 
        
        # Update ASE without skipping charge update
        ase_interface.update_ase(skip_charge_update=False)
        
        # Check that charges are updated
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.get_initial_charges(), np.copy(self.system.part.all().q))

    def test_update_ase_charges_with_skip(self):
        """Test update_ase() with charges enabled but skip_charge_update=True."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=self.system,
            type_mapping=self.type_mapping,
            particle_slice=self.system.part.all(),
            export_charges=True
        )
        
        # Store original charges
        original_charges = ase_interface.atoms.get_initial_charges().copy()
        
        # Change particle charges
        self.system.part.all().q = [5.,-2.,3.5] 
        # Update ASE with skip_charge_update=True
        ase_interface.update_ase(skip_charge_update=True)
        
        # Check that charges are NOT updated (should still be original)
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.get_initial_charges(), original_charges)

    def test_update_ase_masses_no_skip(self):
        """Test update_ase() with masses enabled and no skip."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=self.system,
            type_mapping=self.type_mapping,
            particle_slice=self.system.part.all(),
            export_masses=True
        )
        
        # Change particle masses
        self.system.part.all().mass =  [10.0, 5.5, 7.2]
        
        # Update ASE without skipping mass update
        ase_interface.update_ase(skip_mass_update=False)
        
        # Check that masses are updated
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.get_masses(), np.copy(self.system.part.all().mass))

    def test_update_ase_masses_with_skip(self):
        """Test update_ase() with masses enabled but skip_mass_update=True."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=self.system,
            type_mapping=self.type_mapping,
            particle_slice=self.system.part.all(),
            export_masses=True
        )
        
        # Store original masses
        original_masses = ase_interface.atoms.get_masses().copy()
        
        # Change particle masses
        self.system.part.by_id(2).mass = 7.2
        
        # Update ASE with skip_mass_update=True
        ase_interface.update_ase(skip_mass_update=True)
        
        # Check that masses are NOT updated (should still be original)
        atoms = ase_interface.atoms
        np.testing.assert_allclose(atoms.get_masses(), original_masses)

    def test_update_ase_momenta_always_updates(self):
        """Test update_ase() with momenta enabled - should always update (no skip option)."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=self.system,
            type_mapping=self.type_mapping,
            particle_slice=self.system.part.all(),
            export_momenta=True
        )
        
        # Change particle velocities
        self.system.part.by_id(0).v = [1.0, 1.1, 1.2]
        self.system.part.by_id(1).v = [2.0, 2.1, 2.2]
        self.system.part.by_id(2).v = [3.0, 3.1, 3.2]
        
        # Update ASE
        ase_interface.update_ase()
        
        # Check that momenta are updated (momentum = mass * velocity)
        atoms = ase_interface.atoms
        expected_momenta = np.array([[1.0*2.0, 1.1*2.0, 1.2*2.0],
                                    [2.0*1.5, 2.1*1.5, 2.2*1.5],
                                    [3.0*3.0, 3.1*3.0, 3.2*3.0]])
        np.testing.assert_allclose(atoms.get_momenta(), expected_momenta)

    def test_error_handling_no_atoms(self):
        """Test that update_ase() raises error when atoms is None."""
        ase_interface = espressomd.plugins.ase.ASEInterface(
            system=self.system,
            type_mapping=self.type_mapping,
            particle_slice=self.system.part.all()
        )
        
        # Set atoms to None to simulate uninitialized state
        ase_interface.atoms = None
        
        # Should raise RuntimeError
        with self.assertRaisesRegex(RuntimeError, "atoms object not initialized"):
            ase_interface.update_ase()

if __name__ == "__main__":
    ut.main()
