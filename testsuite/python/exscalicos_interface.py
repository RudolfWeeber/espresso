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

import numpy as np
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.electrostatics


@utx.skipIfMissingFeatures(["EXSCALICOS"])
class ExscalicosInterface(ut.TestCase):

    system = espressomd.System(box_l=3 * [10.])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def setUp(self):
        self.system.box_l = 3 * [10.]
        self.system.periodicity = 3 * [True]

    def tearDown(self):
        self.system.electrostatics.clear()
        self.system.part.clear()

    def add_charges(self, n_pairs=100, seed=42):
        rng = np.random.default_rng(seed=seed)
        pos = rng.random((2 * n_pairs, 3)) * self.system.box_l
        q = np.repeat([[1., -1.]], n_pairs, axis=0).flatten()
        self.system.part.add(pos=pos, q=q)

    def calc_observables(self):
        self.system.integrator.run(0, recalc_forces=True)
        energy = self.system.analysis.energy()["coulomb"]
        forces = np.copy(self.system.part.all().f)
        return energy, forces

    def p3m_reference(self, prefactor):
        p3m = espressomd.electrostatics.P3M(
            prefactor=prefactor, accuracy=1e-6, mesh=48, cao=7, r_cut=3.,
            alpha=1.13, tune=False)
        self.system.electrostatics.solver = p3m
        observables = self.calc_observables()
        self.system.electrostatics.clear()
        return observables

    def check_observables(self, actor, prefactor, rtol_energy, atol_force):
        ref_energy, ref_forces = self.p3m_reference(prefactor)
        self.system.electrostatics.solver = actor
        energy, forces = self.calc_observables()
        np.testing.assert_allclose(energy, ref_energy, rtol=rtol_energy)
        np.testing.assert_allclose(forces, ref_forces, atol=atol_force)
        self.system.electrostatics.clear()

    def test_actor_exceptions(self):
        Exscalicos = espressomd.electrostatics.Exscalicos
        valid_params = {"accuracy": 1e-3}
        with self.assertRaisesRegex(ValueError, "Parameter 'prefactor' must be > 0"):
            Exscalicos(prefactor=-0.8, method_name="p3m",
                       method_params=valid_params)
        with self.assertRaisesRegex(ValueError, "Method 'impossible' is unknown or not supported in ExScaLiCoS"):
            Exscalicos(prefactor=1., method_name="impossible",
                       method_params=valid_params)
        with self.assertRaisesRegex(ValueError, "Parameter 'p3m_cao' is not a valid ExScaLiCoS parameter"):
            Exscalicos(prefactor=1., method_name="p3m",
                       method_params={"accuracy": 1e-3, "p3m_cao": 7})
        with self.assertRaisesRegex(ValueError, "Parameter 'accuracy' is required"):
            Exscalicos(prefactor=1., method_name="p3m",
                       method_params={"cao": 7})
        with self.assertRaisesRegex(ValueError, "Parameter 'accuracy' must be > 0"):
            Exscalicos(prefactor=1., method_name="p3m",
                       method_params={"accuracy": -1e-3})

        # tuning requires charged particles
        actor = Exscalicos(prefactor=1., method_name="p3m",
                           method_params=valid_params)
        with self.assertRaisesRegex(RuntimeError, "CoulombExscalicos: no charged particles in the system"):
            self.system.electrostatics.solver = actor
        self.assertIsNone(self.system.electrostatics.solver)

        # P3M requires full periodicity
        self.add_charges(n_pairs=10)
        self.system.periodicity = [True, True, False]
        actor = Exscalicos(prefactor=1., method_name="p3m",
                           method_params=valid_params)
        with self.assertRaisesRegex(RuntimeError, r"CoulombExscalicos: requires periodicity \(True, True, True\)"):
            self.system.electrostatics.solver = actor
        self.system.periodicity = [True, True, True]
        self.system.electrostatics.solver = actor
        with self.assertRaisesRegex(Exception, r"CoulombExscalicos: requires periodicity \(True, True, True\)"):
            self.system.periodicity = [True, False, True]
        self.system.periodicity = [True, True, True]
        self.system.electrostatics.clear()

        # errors from the library are forwarded
        actor = Exscalicos(prefactor=1., method_name="p3m",
                           method_params={"accuracy": 1e-3, "r_cut": 6.})
        with self.assertRaisesRegex(RuntimeError, "ExScaLiCoS: .*rMax > 0.5 min"):
            self.system.electrostatics.solver = actor
        self.assertIsNone(self.system.electrostatics.solver)

    def test_actor_parameters(self):
        self.add_charges(n_pairs=50)
        actor = espressomd.electrostatics.Exscalicos(
            prefactor=0.5, method_name="p3m",
            method_params={"accuracy": 1e-3, "r_cut": 2.0, "mesh": 32})
        params = actor.get_params()
        self.assertEqual(params["prefactor"], 0.5)
        self.assertEqual(params["method_name"], "p3m")
        self.assertEqual(set(params["method_params"].keys()),
                         {"accuracy", "r_cut", "mesh"})
        self.assertEqual(params["method_params"]["accuracy"], 1e-3)
        self.assertEqual(params["method_params"]["r_cut"], 2.0)
        np.testing.assert_array_equal(params["method_params"]["mesh"],
                                      [32] * 3)
        self.assertEqual(actor.get_available_methods(), ["p3m"])
        self.assertTrue(actor.get_near_field_delegation())

        self.system.electrostatics.solver = actor
        method_params = actor.method_params
        self.assertEqual(set(method_params.keys()),
                         {"accuracy", "r_cut", "alpha", "cao", "mesh"})
        self.assertAlmostEqual(method_params["accuracy"], 1e-3, delta=1e-15)
        self.assertAlmostEqual(method_params["r_cut"], 2.0, delta=1e-12)
        np.testing.assert_array_equal(method_params["mesh"], [32] * 3)
        self.assertGreater(method_params["alpha"], 0.)
        self.assertGreaterEqual(method_params["cao"], 1)
        self.assertLessEqual(method_params["cao"], 7)
        self.assertAlmostEqual(self.system.cell_system.interaction_range,
                               2.0 + self.system.cell_system.skin, delta=1e-12)

        # the ESPResSo cell system doesn't need the cutoff without delegation
        actor.set_near_field_delegation(delegate=False)
        self.assertFalse(actor.get_near_field_delegation())
        self.assertEqual(self.system.cell_system.interaction_range, -1.)
        actor.set_near_field_delegation(delegate=True)
        self.assertTrue(actor.get_near_field_delegation())
        self.assertAlmostEqual(self.system.cell_system.interaction_range,
                               2.0 + self.system.cell_system.skin, delta=1e-12)

        # tuned values reproduce the same solver
        actor_copy = espressomd.electrostatics.Exscalicos(
            **actor.get_params())
        for key, value in actor.method_params.items():
            np.testing.assert_array_equal(
                np.copy(actor_copy.method_params[key]), np.copy(value))

    def test_tuned_r_cut_fits_cell_system(self):
        self.add_charges(n_pairs=20)
        actor = espressomd.electrostatics.Exscalicos(
            prefactor=1., method_name="p3m", method_params={"accuracy": 1e-3})
        self.system.electrostatics.solver = actor
        r_cut = actor.method_params["r_cut"]
        self.assertGreater(r_cut, 0.)
        local_box_l = np.copy(self.system.box_l) / \
            np.array(self.system.cell_system.node_grid)
        r_cut_max = min(np.min(local_box_l), 0.49 * np.min(self.system.box_l))
        self.assertLessEqual(r_cut + self.system.cell_system.skin,
                             r_cut_max + 1e-12)
        self.system.integrator.run(0, recalc_forces=True)

    def test_p3m_accuracy(self):
        prefactor = 1.3
        self.add_charges(n_pairs=200)
        for delegate in [True, False]:
            for method_params in [
                    {"accuracy": 1e-5, "r_cut": 3.0},
                    {"accuracy": 1e-5, "r_cut": 3.0, "alpha": 1.13,
                     "cao": 7, "mesh": [48, 48, 48]}]:
                with self.subTest(delegate=delegate, **method_params):
                    actor = espressomd.electrostatics.Exscalicos(
                        prefactor=prefactor, method_name="p3m",
                        method_params=method_params)
                    actor.set_near_field_delegation(delegate=delegate)
                    self.check_observables(
                        actor, prefactor, rtol_energy=1e-4, atol_force=1e-3)

    def test_deferred_tuning(self):
        # with a fixed cutoff, the solver can be activated without charges
        prefactor = 1.
        actor = espressomd.electrostatics.Exscalicos(
            prefactor=prefactor, method_name="p3m",
            method_params={"accuracy": 1e-5, "r_cut": 3.0})
        self.system.electrostatics.solver = actor
        self.assertEqual(set(actor.method_params.keys()),
                         {"accuracy", "r_cut"})
        self.system.part.add(pos=[1., 1., 1.])
        energy, forces = self.calc_observables()
        self.assertEqual(energy, 0.)
        np.testing.assert_array_equal(forces, 0.)
        self.system.part.clear()
        self.add_charges(n_pairs=100)
        energy, forces = self.calc_observables()
        self.system.electrostatics.clear()
        ref_energy, ref_forces = self.p3m_reference(prefactor)
        np.testing.assert_allclose(energy, ref_energy, rtol=1e-4)
        np.testing.assert_allclose(forces, ref_forces, atol=1e-3)

    def test_box_change(self):
        prefactor = 1.
        self.add_charges(n_pairs=100)
        actor = espressomd.electrostatics.Exscalicos(
            prefactor=prefactor, method_name="p3m",
            method_params={"accuracy": 1e-5, "r_cut": 3.0})
        self.system.electrostatics.solver = actor
        self.calc_observables()
        self.system.change_volume_and_rescale_particles(11., "xyz")
        self.system.electrostatics.clear()
        self.check_observables(
            actor, prefactor, rtol_energy=1e-4, atol_force=1e-3)

    def test_integration(self):
        self.add_charges(n_pairs=50)
        actor = espressomd.electrostatics.Exscalicos(
            prefactor=1., method_name="p3m",
            method_params={"accuracy": 1e-4, "r_cut": 2.5})
        self.system.electrostatics.solver = actor
        self.system.thermostat.set_langevin(kT=1., gamma=1., seed=42)
        self.system.integrator.run(20)
        self.system.thermostat.turn_off()
        energy, forces = self.calc_observables()
        self.assertTrue(np.isfinite(energy))
        self.assertTrue(np.all(np.isfinite(forces)))

    @utx.skipIfMissingFeatures(["SCAFACOS"])
    @utx.skipIfMissingScafacosMethod("p3m")
    def test_coexistence_with_scafacos(self):
        prefactor = 1.
        self.add_charges(n_pairs=100)
        exscalicos = espressomd.electrostatics.Exscalicos(
            prefactor=prefactor, method_name="p3m",
            method_params={"accuracy": 1e-5, "r_cut": 3.0})
        scafacos = espressomd.electrostatics.Scafacos(
            prefactor=prefactor, method_name="p3m",
            method_params={"p3m_r_cut": 3.0, "tolerance_field": 1e-5})
        self.system.electrostatics.solver = scafacos
        ref_energy, ref_forces = self.calc_observables()
        self.system.electrostatics.solver = exscalicos
        energy, forces = self.calc_observables()
        np.testing.assert_allclose(energy, ref_energy, rtol=1e-4)
        np.testing.assert_allclose(forces, ref_forces, atol=1e-3)
        self.system.electrostatics.solver = scafacos
        energy, forces = self.calc_observables()
        np.testing.assert_allclose(energy, ref_energy, rtol=1e-10)
        np.testing.assert_allclose(forces, ref_forces, atol=1e-10)


if __name__ == "__main__":
    ut.main()
