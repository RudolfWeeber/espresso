#
# Copyright (C) 2025 The ESPResSo project
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

import numpy as np
import unittest as ut
import scipy.optimize
import importlib_wrapper

sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@SAMPLES_DIR@/lj_liquid_ase.py", int_steps=10, int_n_times=60)


class Sample(ut.TestCase):

    @staticmethod
    def exp_kernel(xdata, a, b, c):
        return a * np.exp(-b * xdata) + c

    def test(self):
        system = sample.system
        ydata = np.array(system.instantaneous_temperatures)
        xdata = np.arange(len(ydata))
        popt, pcov = scipy.optimize.curve_fit(
            self.exp_kernel, xdata, ydata, p0=[-1., 0.1, 1.])
        pvar = np.diag(pcov)
        np.testing.assert_allclose(np.sqrt(pvar) / popt, 0., atol=0.1)
        self.assertAlmostEqual(popt[2], system.thermostat.kT, delta=1e-2)
        self.assertAlmostEqual(popt[1], 7.5e-2, delta=1e-2)
        self.assertAlmostEqual(popt[0], -0.9, delta=0.1)


if __name__ == "__main__":
    ut.main()
