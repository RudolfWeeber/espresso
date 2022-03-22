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

"""
Testmodule for the H5MD interface.
"""
import os
import sys
import unittest as ut
import unittest_decorators as utx
import numpy as np
import espressomd
import espressomd.interactions
import espressomd.io.writer
import tempfile
try:
    import h5py  # h5py has to be imported *after* espressomd (MPI)
    skipIfMissingPythonPackage = utx.no_skip
except ImportError:
    skipIfMissingPythonPackage = ut.skip(
        "Python module h5py not available, skipping test!")


N_PART = 26


@utx.skipIfMissingFeatures(['H5MD'])
@skipIfMissingPythonPackage
class H5mdTests(ut.TestCase):
    """
    Test the core implementation of writing hdf5 files.

    """
    box_l = N_PART // 2
    box_l = [box_l, box_l + 1, box_l + 2]
    system = espressomd.System(box_l=box_l)
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    # set particles outside the main box to verify position folding
    for i in range(N_PART):
        p = system.part.add(id=i, pos=np.array(3 * [i - 4], dtype=float),
                            v=np.array([1.0, 2.0, 3.0]), type=23)
        if espressomd.has_features(['MASS']):
            p.mass = 2.3
        if espressomd.has_features(['EXTERNAL_FORCES']):
            p.ext_force = [0.1, 0.2, 0.3]
        if espressomd.has_features(['ELECTROSTATICS']):
            p.q = i

    vb = espressomd.interactions.Virtual()
    system.bonded_inter.add(vb)

    for i in range(N_PART - 1):
        system.part.by_id(i).add_bond((vb, i + 1))

    system.integrator.run(steps=0)
    system.time = 12.3
    n_nodes = system.cell_system.get_state()["n_nodes"]

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.TemporaryDirectory()
        cls.temp_file = os.path.join(cls.temp_dir.name, 'test.h5')
        h5_units = espressomd.io.writer.h5md.UnitSystem(
            time='ps', mass='u', length='m', charge='e')
        h5 = espressomd.io.writer.h5md.H5md(
            file_path=cls.temp_file, unit_system=h5_units)
        h5.write()
        h5.write()
        h5.flush()
        h5.close()
        cls.h5_params = h5.get_params()
        cls.py_file = h5py.File(cls.temp_file, 'r')
        cls.py_pos = cls.py_file['particles/atoms/position/value'][1]
        cls.py_img = cls.py_file['particles/atoms/image/value'][1]
        cls.py_mass = cls.py_file['particles/atoms/mass/value'][1]
        cls.py_vel = cls.py_file['particles/atoms/velocity/value'][1]
        cls.py_charge = cls.py_file['particles/atoms/charge/value'][1]
        cls.py_f = cls.py_file['particles/atoms/force/value'][1]
        cls.py_id = cls.py_file['particles/atoms/id/value'][1]
        cls.py_id_time = cls.py_file['particles/atoms/id/time'][1]
        cls.py_id_step = cls.py_file['particles/atoms/id/step'][1]
        cls.py_bonds = cls.py_file['connectivity/atoms/value'][1]
        cls.py_box = cls.py_file['particles/atoms/box/edges/value'][1]

    @classmethod
    def tearDownClass(cls):
        cls.temp_dir.cleanup()

    def test_opening(self):
        h5 = espressomd.io.writer.h5md.H5md(file_path=self.temp_file)
        h5.close()

    @ut.skipIf(n_nodes > 1, "only runs for 1 MPI rank")
    def test_exceptions(self):
        h5md = espressomd.io.writer.h5md
        h5_units = h5md.UnitSystem(time='ps', mass='u', length='m', charge='e')
        temp_file = os.path.join(self.temp_dir.name, 'exceptions.h5')
        # write a non-compliant file
        with open(temp_file, 'w') as _:
            pass
        with self.assertRaisesRegex(RuntimeError, 'not a valid HDF5 file'):
            h5md.H5md(file_path=temp_file, unit_system=h5_units)
        # write a leftover backup file
        os.remove(temp_file)
        with open(temp_file + '.bak', 'w') as _:
            pass
        with self.assertRaisesRegex(RuntimeError, 'A backup of the .h5 file exists'):
            h5md.H5md(file_path=temp_file, unit_system=h5_units)

    def test_box(self):
        np.testing.assert_allclose(self.py_box, self.box_l)

    def test_metadata(self):
        """Test if the H5MD metadata has been written properly."""
        self.assertEqual(self.py_file['h5md'].attrs['version'][0], 1)

        self.assertEqual(self.py_file['h5md'].attrs['version'][1], 1)
        self.assertIn('creator', self.py_file['h5md'])
        self.assertIn('name', self.py_file['h5md/creator'].attrs)
        self.assertIn('version', self.py_file['h5md/creator'].attrs)
        self.assertEqual(
            self.py_file['h5md/creator'].attrs['name'][:], b'ESPResSo')
        self.assertIn('author', self.py_file['h5md'])
        self.assertIn('name', self.py_file['h5md/author'].attrs)

    def test_pos(self):
        """Test if positions have been written properly."""
        pos_ref = np.outer(np.arange(N_PART) - 4, [1, 1, 1])
        pos_ref = np.mod(pos_ref, self.box_l)
        pos_read = [x for (_, x) in sorted(zip(self.py_id, self.py_pos))]
        np.testing.assert_allclose(pos_read, pos_ref)

    def test_time(self):
        """Test for time dataset."""
        self.assertEqual(self.py_id_time, 12.3)

    def test_img(self):
        """Test if images have been written properly."""
        pos_ref = np.outer(np.arange(N_PART) - 4, [1, 1, 1])
        images_ref = np.floor_divide(pos_ref, self.box_l)
        images_read = [x for (_, x) in sorted(zip(self.py_id, self.py_img))]
        np.testing.assert_allclose(images_read, images_ref)

    @utx.skipIfMissingFeatures("MASS")
    def test_mass(self):
        """Test if masses have been written correct."""
        np.testing.assert_allclose(self.py_mass, 2.3)

    @utx.skipIfMissingFeatures(['ELECTROSTATICS'])
    def test_charge(self):
        """Test if charges have been written properly."""
        charges = np.arange(N_PART)
        np.testing.assert_allclose(
            [x for (_, x) in sorted(zip(self.py_id, self.py_charge))], charges)

    def test_vel(self):
        """Test if velocities have been written properly."""
        np.testing.assert_allclose(
            np.array([[1.0, 2.0, 3.0] for _ in range(N_PART)]),
            np.array([x for (_, x) in sorted(zip(self.py_id, self.py_vel))]),
            err_msg="Velocities not written correctly by H5md!")

    @utx.skipIfMissingFeatures(['EXTERNAL_FORCES'])
    def test_f(self):
        """Test if forces have been written properly."""
        np.testing.assert_allclose(
            np.array([[0.1, 0.2, 0.3] for _ in range(N_PART)]),
            np.array([x for (_, x) in sorted(zip(self.py_id, self.py_f))]),
            err_msg="Forces not written correctly by H5md!")

    def test_bonds(self):
        """Test if bonds have been written properly."""
        self.assertEqual(len(self.py_bonds), N_PART - 1)

        for i in range(N_PART - 1):
            bond = [x for x in self.py_bonds if x[0] == i][0]
            self.assertEqual(bond[0], i + 0)
            self.assertEqual(bond[1], i + 1)

    def test_script(self):
        with open(sys.argv[0], 'r') as f:
            ref = f.read()
        data = self.py_file['parameters/files'].attrs['script'].decode('utf-8')
        self.assertEqual(data, ref)

    def test_units(self):
        def get_unit(path):
            return self.py_file[path].attrs['unit']
        self.assertEqual(get_unit('particles/atoms/id/time'), b'ps')
        self.assertEqual(get_unit('particles/atoms/position/value'), b'm')
        if espressomd.has_features(['ELECTROSTATICS']):
            self.assertEqual(get_unit('particles/atoms/charge/value'), b'e')
        if espressomd.has_features(['MASS']):
            self.assertEqual(get_unit('particles/atoms/mass/value'), b'u')
        self.assertEqual(get_unit('particles/atoms/force/value'), b'm u ps-2')
        self.assertEqual(get_unit('particles/atoms/velocity/value'), b'm ps-1')

    def test_getters(self):
        self.assertEqual(self.h5_params['file_path'], self.temp_file)
        self.assertEqual(os.path.abspath(self.h5_params['script_path']),
                         os.path.abspath(__file__))
        self.assertEqual(self.h5_params['time_unit'], 'ps')
        if espressomd.has_features(['ELECTROSTATICS']):
            self.assertEqual(self.h5_params['charge_unit'], 'e')
        if espressomd.has_features(['MASS']):
            self.assertEqual(self.h5_params['mass_unit'], 'u')
        self.assertEqual(self.h5_params['force_unit'], 'm u ps-2')
        self.assertEqual(self.h5_params['velocity_unit'], 'm ps-1')

    def test_links(self):
        time_ref = self.py_id_time
        step_ref = self.py_id_step
        for group in "position", "velocity", "force", "charge", "mass", "image":
            time = self.py_file['particles/atoms/' + group + '/time'][1]
            step = self.py_file['particles/atoms/' + group + '/step'][1]
            self.assertEqual(time, time_ref)
            self.assertEqual(step, step_ref)

        bond_time = self.py_file['connectivity/atoms/time'][1]
        self.assertEqual(bond_time, time_ref)
        bond_step = self.py_file['connectivity/atoms/step'][1]
        self.assertEqual(bond_step, step_ref)
        box_time = self.py_file['particles/atoms/box/edges/time'][1]
        self.assertEqual(box_time, time_ref)
        box_step = self.py_file['particles/atoms/box/edges/step'][1]
        self.assertEqual(box_step, step_ref)


if __name__ == "__main__":
    ut.main()
