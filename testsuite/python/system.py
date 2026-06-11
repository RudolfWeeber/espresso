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

import unittest as ut
import numpy as np
import espressomd


class Test(ut.TestCase):

    def test_isotropic_rescale_noncubic_box(self):
        """Regression test for bug #16: isotropic ('xyz') box rescale silently
        corrupts particle positions on a non-cubic box.

        For coord==3 the old code computed a single scalar scale from the
        x-axis only and applied it to all three particle coordinates.  Each
        coordinate must instead scale by new_box_l[axis] / old_box_l[axis].
        """
        old_box_l = np.array([10., 20., 30.])
        system = espressomd.System(box_l=old_box_l)
        pos = np.array([1., 2., 3.])
        p = system.part.add(pos=pos)
        d_new = 5.
        system.change_volume_and_rescale_particles(d_new, dir="xyz")
        new_box_l = np.copy(system.box_l)
        np.testing.assert_allclose(
            new_box_l, [d_new, d_new, d_new], rtol=0., atol=1e-12)
        expected = pos * (new_box_l / old_box_l)
        np.testing.assert_allclose(
            np.copy(p.pos), expected, rtol=0., atol=1e-12,
            err_msg="isotropic rescale must scale each coordinate by "
                    "new_box_l[axis]/old_box_l[axis] (bug #16: y/z scaled by "
                    "the x-based factor instead)")

    def test_system(self):
        with self.assertRaisesRegex(ValueError, "Required argument 'box_l' not provided"):
            espressomd.System()
        with self.assertRaisesRegex(ValueError, "Property 'unknown' cannot be set via argument to System class"):
            espressomd.System(box_l=[1., 1., 1.], unknown=1)
        system = espressomd.System(box_l=[1., 1., 1.], min_global_cut=0.01)
        self.assertEqual(system.min_global_cut, 0.01)
        self.assertEqual(system.time_step, -1.)
        with self.assertRaisesRegex(RuntimeError, "You can only have one instance of the system class at a time"):
            espressomd.System(box_l=[1., 1., 1.])


if __name__ == "__main__":
    ut.main()
