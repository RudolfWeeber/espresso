# Copyright (C) 2010-2019 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.lb
import espressomd.shapes
import espressomd.lbboundaries
import itertools
import numpy as np


class LBBoundariesBase:
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.cell_system.skin = 0.1

    wall_shape1 = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=2.5)
    wall_shape2 = espressomd.shapes.Wall(normal=[-1., 0., 0.], dist=-7.5)

    def setUp(self):
        self.lbf = self.lb_class(
            viscosity=1.0, density=1.0, agrid=0.5, tau=1.0, **self.lb_params)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        self.system.lbboundaries.clear()
        self.system.actors.clear()

    def test_add(self):
        boundary = espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1)

        self.system.lbboundaries.add(boundary)
        self.assertEqual(boundary, self.system.lbboundaries[0])

    def test_remove(self):
        lbb = self.system.lbboundaries

        b1 = lbb.add(
            espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        b2 = lbb.add(
            espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))

        lbb.remove(b1)

        self.assertNotIn(b1, lbb)
        self.assertIn(b2, lbb)

    def test_size(self):
        lbb = self.system.lbboundaries
        self.assertEqual(lbb.size(), 0)

        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        self.assertEqual(lbb.size(), 1)

        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        self.assertEqual(lbb.size(), 2)

    def test_empty(self):
        lbb = self.system.lbboundaries
        self.assertTrue(lbb.empty())

        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        self.assertFalse(lbb.empty())

    def test_clear(self):
        lbb = self.system.lbboundaries

        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))

        lbb.clear()

        self.assertTrue(lbb.empty())

    def check_boundary_flags(self, slip_velocity1, slip_velocity2):
        rng = range(20)

        for i in itertools.product(range(0, 5), rng, rng):
            self.assertTrue(self.lbf[i].is_boundary)
            np.testing.assert_allclose(
                np.copy(self.lbf[i].velocity), slip_velocity1)

        for i in itertools.product(range(5, 15), rng, rng):
            self.assertFalse(self.lbf[i].is_boundary)

        for i in itertools.product(range(15, 20), rng, rng):
            self.assertTrue(self.lbf[i].is_boundary)
            np.testing.assert_allclose(
                np.copy(self.lbf[i].velocity), slip_velocity2)

        self.system.lbboundaries.clear()
        self.lbf.clear_boundaries()
        for i in itertools.product(rng, rng, rng):
            self.assertFalse(self.lbf[i].is_boundary)

    def test_boundary_flags(self):
        lbb = self.system.lbboundaries

        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape1))
        lbb.add(espressomd.lbboundaries.LBBoundary(shape=self.wall_shape2))

        self.check_boundary_flags(np.copy(lbb[0].velocity),
                                  np.copy(lbb[1].velocity))

        slip_velocity1 = 1e-3 * np.array([1., 2., 3.])
        slip_velocity2 = 1e-3 * np.array([4., 5., 6.])
        self.lbf.add_boundary_from_shape(self.wall_shape1, slip_velocity1)
        self.lbf.add_boundary_from_shape(self.wall_shape2, slip_velocity2)
        self.check_boundary_flags(slip_velocity1, slip_velocity2)

        new_slip_velocity1 = 1e-3 * np.array([-1., 0., 1.])
        new_slip_velocity2 = 1e-3 * np.array([-2., 2., 2.])
        bitmask1 = self.lbf.get_shape_bitmask(self.wall_shape1)
        bitmask2 = self.lbf.get_shape_bitmask(self.wall_shape2)
        nodes1 = np.array(np.nonzero(bitmask1)).T
        nodes2 = np.array(np.nonzero(bitmask2)).T
        self.lbf.add_boundary_from_list(nodes1, new_slip_velocity1)
        self.lbf.add_boundary_from_list(
            nodes2, len(nodes2) * [new_slip_velocity2])
        self.check_boundary_flags(new_slip_velocity1, new_slip_velocity2)

    def test_union(self):
        union = espressomd.shapes.Union()
        union.add([self.wall_shape1, self.wall_shape2])

        lbb = self.system.lbboundaries.add(
            espressomd.lbboundaries.LBBoundary(shape=union))
        self.check_boundary_flags(np.copy(lbb.velocity), np.copy(lbb.velocity))

        slip_velocity = 1e-3 * np.array([1., 2., 3.])
        self.lbf.add_boundary_from_shape(union, slip_velocity)
        self.check_boundary_flags(slip_velocity, slip_velocity)


@utx.skipIfMissingFeatures(["LB_BOUNDARIES"])
@utx.skipIfMissingFeatures(["LB_WALBERLA"])
class LBBoundariesWalberla(LBBoundariesBase, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}


@utx.skipIfMissingFeatures(["LB_BOUNDARIES"])
@utx.skipIfMissingFeatures(["LB_WALBERLA"])
class LBBoundariesWalberlaSinglePrecision(LBBoundariesBase, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}


if __name__ == "__main__":
    ut.main()
