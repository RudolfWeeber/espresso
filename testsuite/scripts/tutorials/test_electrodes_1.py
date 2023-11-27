#
# Copyright (C) 2023 The ESPResSo project
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
import importlib_wrapper
import numpy as np

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/electrodes/electrodes_part1.py", N_AXIAL_POINTS=6,
    script_suffix="@TEST_SUFFIX@")


@skipIfMissingFeatures
class Tutorial(ut.TestCase):

    def test_force(self):
        ref_force = tutorial.analytic_force_centered(
            tutorial.r, tutorial.box_l_z)
        msg = "The force for particle 1 doesn't agree with analytical expression"
        np.testing.assert_allclose(np.log(tutorial.elc_forces_axial[:, 0]),
                                   np.log(ref_force), rtol=0.05, err_msg=msg)
        msg = "The force for particle 2 doesn't agree with analytical expression"
        np.testing.assert_allclose(np.log(-tutorial.elc_forces_axial[:, 1]),
                                   np.log(ref_force), rtol=0.05, err_msg=msg)


if __name__ == "__main__":
    ut.main()
