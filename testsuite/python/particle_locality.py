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


system = espressomd.System(box_l=[10., 10., 10.])
n_nodes = system.cell_system.get_state()["n_nodes"]


@ut.skipIf(n_nodes not in (2, 3),
           "Requires 2 or 3 MPI ranks: the box has to be split for the "
           "invariant to be testable, while each local box must stay "
           "wider than the interaction range")
@utx.skipIfMissingFeatures(["P3M", "VIRTUAL_SITES_RELATIVE"])
class ParticleLocality(ut.TestCase):

    """
    A local particle's stored position must stay the periodic image closest to
    the MPI rank that owns it.  Code that folds a position outside a resort
    (virtual sites, the Lees-Edwards push) can break that: the particle keeps
    its owner but its coordinates jump a whole box.  Pair loops are immune,
    because they use the minimum image, but consumers that index rank-local
    storage by absolute position are not -- the P3M charge-assignment mesh
    reads and writes outside its allocation.

    The skin only controls how often the cell system re-sorts, never the
    physics, so forces must not depend on it.  A large skin suppresses the
    re-sort that used to hide the folding; a small skin does not.
    """

    system = system

    def setUp(self):
        self.system.time_step = 0.005
        self.system.min_global_cut = 0.5
        # split along x so that the far boundary is a domain boundary
        self.system.cell_system.node_grid = [n_nodes, 1, 1]

    def tearDown(self):
        self.system.part.clear()
        self.system.electrostatics.clear()
        self.system.min_global_cut = 0.

    def forces_after_wrapping(self, skin, n_steps=20):
        """Drag a charged virtual site across the far box boundary, so that it
        is folded while the rank at the opposite end of the box still owns it,
        then report the forces."""
        system = self.system
        system.part.clear()
        system.cell_system.skin = skin
        # the reference sits just short of the boundary, the virtual site just
        # beyond it, both owned by the rank covering the upper half of x
        ref = system.part.add(pos=[9.68, 5., 5.], v=[1., 0., 0.], q=0.)
        vs = system.part.add(pos=[9.98, 5., 5.], q=1.)
        vs.vs_auto_relate_to(ref)
        system.part.add(pos=[2., 2., 2.], q=-1.)
        system.electrostatics.solver = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=1e-2, mesh=16, cao=3, r_cut=2.0,
            alpha=1.2, tune=False)
        system.integrator.run(n_steps)
        return np.copy(system.part.all().f), np.copy(system.part.all().pos)

    def test_forces_are_skin_independent_across_a_fold(self):
        f_small, pos_small = self.forces_after_wrapping(skin=0.1)
        f_large, pos_large = self.forces_after_wrapping(skin=0.8)
        # the trajectory must be the same: the skin is not a physical parameter
        np.testing.assert_allclose(pos_large, pos_small, rtol=0., atol=1e-9)
        np.testing.assert_allclose(f_large, f_small, rtol=1e-6, atol=1e-9)

    def test_virtual_site_stays_with_its_owner(self):
        """A virtual site that crosses the box boundary together with its
        reference must keep its stored position next to its owner instead of
        being folded.  The stored position is not visible from Python (``pos``
        is unfolded and ``image_box`` is derived from it), but a fold is: it
        puts the site a box length outside its owner's domain, which the
        resort criterion answers by handing the site to another rank."""
        system = self.system
        system.part.clear()
        system.cell_system.skin = 0.8
        ref = system.part.add(pos=[9.68, 5., 5.], v=[1., 0., 0.], q=0.)
        vs = system.part.add(pos=[9.98, 5., 5.], q=0.)
        vs.vs_auto_relate_to(ref)
        # neither particle travels more than half the skin in total, so the
        # cell system has no reason to re-sort and both stay where they are
        owner = ref.node
        self.assertEqual(vs.node, owner)
        for _ in range(12):
            system.integrator.run(4)
            self.assertEqual(
                vs.node, owner, "the virtual site was folded away from its "
                "owner and had to be re-sorted")
            d = np.copy(vs.pos) - np.copy(ref.pos)
            np.testing.assert_allclose(np.linalg.norm(d), 0.3, atol=1e-9)


@ut.skipIf(n_nodes < 2,
           "Requires at least 2 MPI ranks: on one rank every decomposition "
           "owns the whole box")
@utx.skipIfMissingFeatures(["LENNARD_JONES"])
class IdOwnedParticles(ut.TestCase):

    """
    The n-square decomposition, and the n-square part of the hybrid one,
    distribute particles by id, so a local particle may sit anywhere in the
    box.  The resort criterion must not mistake that for a particle that left
    its owner's domain, or the Verlet list is rebuilt on every step.
    """

    system = system

    def tearDown(self):
        self.system.part.clear()
        self.system.non_bonded_inter.reset()
        self.system.thermostat.turn_off()
        self.system.cell_system.set_regular_decomposition()

    def verlet_reuse(self, set_decomposition, **kwargs):
        system = self.system
        system.time_step = 0.01
        system.cell_system.skin = 0.4
        set_decomposition(**kwargs)
        rng = np.random.default_rng(1)
        n = 200
        system.part.add(pos=rng.random((n, 3)) * 10.,
                        type=[0] * (n // 2) + [1] * (n // 2))
        for types in [(0, 0), (0, 1), (1, 1)]:
            system.non_bonded_inter[types].lennard_jones.set_params(
                epsilon=1., sigma=1., cutoff=1.12246, shift="auto")
        system.integrator.set_steepest_descent(
            f_max=0., gamma=1., max_displacement=0.01)
        system.integrator.run(50)
        system.integrator.set_vv()
        system.thermostat.set_langevin(kT=0.1, gamma=1., seed=2)
        system.integrator.run(100)
        return system.cell_system.get_state()["verlet_reuse"]

    def check_verlet_reuse(self, reuse):
        # verlet_reuse is 0 when the list was never rebuilt during the last
        # run, which is the best possible outcome, not a failure
        self.assertTrue(
            reuse == 0. or reuse > 5.,
            f"a locality check that only holds for positional ownership "
            f"rebuilds the Verlet list on every step, got {reuse}")

    def test_n_square(self):
        self.check_verlet_reuse(
            self.verlet_reuse(self.system.cell_system.set_n_square))

    def test_hybrid(self):
        self.check_verlet_reuse(self.verlet_reuse(
            self.system.cell_system.set_hybrid_decomposition,
            n_square_types={1}, cutoff_regular=1.2))


if __name__ == "__main__":
    ut.main()
