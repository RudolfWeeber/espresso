#
# Copyright (C) 2026 The ESPResSo project
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
import numpy as np

import espressomd
import espressomd.electrostatics
import tests_common


@utx.skipIfMissingFeatures(["ELECTROSTATICS"])
class CoulombSilicaMelt(ut.TestCase):

    """
    Compare P3M, ScaFaCoS P3M and ExScaLiCoS P3M electrostatic forces
    and energy against high-precision reference data of the silica melt
    benchmark system from the ScaFaCoS project.

    """

    system = espressomd.System(box_l=3 * [62.05966799605923])
    data = np.genfromtxt(tests_common.data_path(
        "coulomb_silica_melt_12960.data"))

    # parameters tuned for an RMS force accuracy of 1e-4
    p3m_params = {"r_cut": 15.0301, "mesh": 28, "cao": 6,
                  "alpha": 0.193925}

    # reference energy for a unit prefactor
    reference_energy = -16623.245931437410036

    def setUp(self):
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.
        # Input format: id pos q f
        self.system.part.add(pos=self.data[:, 1:4], q=self.data[:, 4])
        self.reference_forces = self.data[:, 5:8]

    def tearDown(self):
        self.system.part.clear()
        self.system.electrostatics.clear()

    def compare(self, method_name, prefactor):
        err_msg = f"difference too large for method {method_name}"
        self.system.integrator.run(0, recalc_forces=True)
        np.testing.assert_allclose(
            np.copy(self.system.part.all().f) / prefactor,
            self.reference_forces, atol=1e-3,
            err_msg=f"Force {err_msg}")
        self.assertAlmostEqual(
            self.system.analysis.energy()["coulomb"] / prefactor,
            self.reference_energy, delta=5e-2,
            msg=f"Energy {err_msg}")

    @utx.skipIfMissingFeatures(["P3M"])
    def test_p3m_cpu(self):
        self.system.electrostatics.solver = espressomd.electrostatics.P3M(
            prefactor=2.2, accuracy=1e-4, tune=False, check_neutrality=False,
            **self.p3m_params)
        self.compare("p3m", prefactor=2.2)

    @utx.skipIfMissingFeatures(["SCAFACOS"])
    @utx.skipIfMissingScafacosMethod("p3m")
    def test_scafacos_p3m(self):
        self.system.electrostatics.solver = espressomd.electrostatics.Scafacos(
            prefactor=2.8, method_name="p3m", check_neutrality=False,
            method_params={"p3m_r_cut": self.p3m_params["r_cut"],
                           "p3m_grid": self.p3m_params["mesh"],
                           "p3m_cao": self.p3m_params["cao"],
                           "p3m_alpha": self.p3m_params["alpha"]})
        self.compare("scafacos_p3m", prefactor=2.8)

    @utx.skipIfMissingFeatures(["EXSCALICOS"])
    def test_exscalicos_p3m(self):
        for delegate in [True, False]:
            with self.subTest(delegate=delegate):
                actor = espressomd.electrostatics.Exscalicos(
                    prefactor=2.6, method_name="p3m", check_neutrality=False,
                    method_params={"accuracy": 1e-4, **self.p3m_params})
                actor.set_near_field_delegation(delegate=delegate)
                self.system.electrostatics.solver = actor
                self.compare("exscalicos_p3m", prefactor=2.6)
                self.system.electrostatics.clear()


if __name__ == "__main__":
    ut.main()
