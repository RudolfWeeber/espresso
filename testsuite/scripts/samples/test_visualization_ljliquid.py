# Copyright (C) 2019-2022 The ESPResSo project
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
import importlib_wrapper


def disable_GUI(code):
    # integrate without visualizer
    breakpoint = "visualizer.update()"
    assert breakpoint in code
    code = code.replace(breakpoint, "")
    breakpoint = "t = threading.Thread(target=main_thread)"
    assert breakpoint in code
    code = code.split(breakpoint, 1)[0] + "main_thread()"
    return code


sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "/home/thilo/code/espresso2/espresso/testsuite/scripts/samples/local_samples/visualization_ljliquid.py",
    substitutions=disable_GUI, int_n_times=5)


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system


if __name__ == "__main__":
    ut.main()
