#
# Copyright (C) 2010-2022 The ESPResSo project
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

import pathlib
import tempfile
import contextlib
import numpy as np

with contextlib.suppress(ImportError):
    import vtk
    import vtk.util.numpy_support

import espressomd
import espressomd.lb
if espressomd.has_features('LB_BOUNDARIES'):
    import espressomd.lbboundaries
    import espressomd.shapes


class TestLBWrite:
    system = espressomd.System(box_l=[10, 11, 12])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def set_lbf(self):
        # setup LB system
        lbf = self.lb_class(
            kT=1, agrid=1.0, dens=1.0, visc=1.0, tau=0.1, seed=42,
            ext_force_density=[0, 0.03, 0])
        self.system.actors.add(lbf)
        if espressomd.has_features('LB_BOUNDARIES'):
            self.system.lbboundaries.add(espressomd.lbboundaries.LBBoundary(
                shape=espressomd.shapes.Wall(normal=[1, 0, 0], dist=1.5)))
            self.system.lbboundaries.add(espressomd.lbboundaries.LBBoundary(
                shape=espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-8.5)))
        return lbf

    def parse_vtk(self, filepath, name, shape):
        reader = vtk.vtkStructuredPointsReader()
        reader.SetFileName(str(filepath))
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()

        data = reader.GetOutput()
        points = data.GetPointData()

        return vtk.util.numpy_support.vtk_to_numpy(
            points.GetArray(name)).reshape(shape, order='F')

    def test_vtk(self):
        '''
        Check VTK files.
        '''

        with tempfile.TemporaryDirectory() as tmp_directory:
            path_vtk_root = pathlib.Path(tmp_directory)
            path_vtk_boundary = path_vtk_root / 'boundary.vtk'
            path_vtk_velocity = path_vtk_root / 'velocity.vtk'
            path_vtk_velocity_bb = path_vtk_root / 'velocity_bb.vtk'
            path_vtk_skip = path_vtk_root / 'skip.vtk'
            path_vtk_invalid = path_vtk_root / 'non_existent_folder' / 'file'

            shape = [10, 11, 12]
            lbf = self.set_lbf()
            self.system.integrator.run(100)

            # write VTK files
            with self.assertRaises(RuntimeError):
                lbf.write_vtk_velocity(str(path_vtk_invalid))
            with self.assertRaises(RuntimeError):
                lbf.write_vtk_boundary(str(path_vtk_invalid))
            lbf.write_vtk_boundary(str(path_vtk_boundary))
            lbf.write_vtk_velocity(str(path_vtk_velocity))
            with self.assertRaises(ValueError):
                lbf.write_vtk_velocity(str(path_vtk_skip), 3 * [0], None)
            with self.assertRaises(ValueError):
                lbf.write_vtk_velocity(str(path_vtk_skip), None, 3 * [0])
            with self.assertRaises(RuntimeError):
                lbf.write_vtk_velocity(str(path_vtk_skip), [-2, 1, 1], 3 * [1])
            with self.assertRaises(RuntimeError):
                lbf.write_vtk_velocity(str(path_vtk_skip), 3 * [0], [1, 2, 16])
            with self.assertRaises(ValueError):
                lbf.write_vtk_velocity(str(path_vtk_skip), [1, 1], 3 * [1])
            with self.assertRaises(ValueError):
                lbf.write_vtk_velocity(
                    str(path_vtk_skip), 3 * [1], np.array([2, 3]))
            bb1, bb2 = ([1, 2, 3], [9, 10, 11])
            lbf.write_vtk_velocity(str(path_vtk_velocity_bb), bb1, bb2)

            # check VTK values match node values
            node_velocity = np.zeros(shape + [3])
            node_boundary = np.zeros(shape, dtype=int)
            for i in range(shape[0]):
                for j in range(shape[1]):
                    for k in range(shape[2]):
                        node = lbf[i, j, k]
                        node_velocity[i, j, k] = node.velocity
                        node_boundary[i, j, k] = node.boundary
            node_velocity_bb = node_velocity[bb1[0]:bb2[0],
                                             bb1[1]:bb2[1],
                                             bb1[2]:bb2[2]]

            vtk_velocity = self.parse_vtk(path_vtk_velocity, 'velocity',
                                          node_velocity.shape)
            np.testing.assert_allclose(vtk_velocity, node_velocity, atol=5e-7)

            vtk_velocity_bb = self.parse_vtk(path_vtk_velocity_bb, 'velocity',
                                             node_velocity_bb.shape)
            np.testing.assert_allclose(
                vtk_velocity_bb, node_velocity_bb, atol=5e-7)

            vtk_boundary = self.parse_vtk(path_vtk_boundary, 'boundary', shape)
            np.testing.assert_equal(vtk_boundary, node_boundary.astype(int))

    def test_print(self):
        '''
        Check data files.
        '''

        with tempfile.TemporaryDirectory() as tmp_directory:
            path_dat_root = pathlib.Path(tmp_directory)
            path_dat_boundary = path_dat_root / 'boundary.vtk'
            path_dat_velocity = path_dat_root / 'velocity.vtk'
            path_dat_invalid = path_dat_root / 'non_existent_folder' / 'file'

            shape = [10, 11, 12]
            lbf = self.set_lbf()
            self.system.integrator.run(100)

            # write data files
            with self.assertRaises(RuntimeError):
                lbf.write_velocity(str(path_dat_invalid))
            with self.assertRaises(RuntimeError):
                lbf.write_boundary(str(path_dat_invalid))
            lbf.write_boundary(str(path_dat_boundary))
            lbf.write_velocity(str(path_dat_velocity))

            # check data values match node values
            node_velocity = np.zeros(shape + [3])
            node_boundary = np.zeros(shape, dtype=int)
            for i in range(shape[0]):
                for j in range(shape[1]):
                    for k in range(shape[2]):
                        node = lbf[i, j, k]
                        node_velocity[i, j, k] = node.velocity
                        node_boundary[i, j, k] = node.boundary

            ref_coord = np.array([
                np.tile(np.arange(shape[0]), shape[1] * shape[2]),
                np.tile(np.repeat(np.arange(shape[1]), shape[0]), shape[2]),
                np.repeat(np.arange(shape[2]), shape[0] * shape[1])]).T

            dat_velocity = np.loadtxt(path_dat_velocity)
            dat_coord = (dat_velocity[:, 0:3] - 0.5).astype(int)
            np.testing.assert_equal(dat_coord, ref_coord)
            dat_vel = dat_velocity[:, 3:]
            ref_vel = np.swapaxes(node_velocity, 0, 2).reshape((-1, 3))
            np.testing.assert_allclose(dat_vel, ref_vel, atol=5e-7)

            dat_boundary = np.loadtxt(path_dat_boundary)
            dat_coord = (dat_boundary[:, 0:3] - 0.5).astype(int)
            np.testing.assert_equal(dat_coord, ref_coord)
            dat_bound = dat_boundary[:, 3].astype(int)
            ref_bound = np.swapaxes(node_boundary, 0, 2).reshape(-1)
            if isinstance(lbf, espressomd.lb.LBFluid):
                ref_bound = (ref_bound != 0).astype(int)
            np.testing.assert_equal(dat_bound, ref_bound)


@utx.skipIfMissingModules("vtk")
class TestLBWriteCPU(TestLBWrite, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluid


@utx.skipIfMissingGPU()
@utx.skipIfMissingModules("vtk")
class TestLBWriteGPU(TestLBWrite, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluidGPU


if __name__ == '__main__':
    ut.main()
