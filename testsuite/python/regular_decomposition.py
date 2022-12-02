#
# Copyright (C) 2013-2022 The ESPResSo project
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
import espressomd
import numpy as np
import itertools

np.random.seed(42)


class RegularDecomposition(ut.TestCase):
    system = espressomd.System(box_l=3 * [50.0])
    original_node_grid = tuple(system.cell_system.node_grid)

    def setUp(self):
        self.system.cell_system.set_regular_decomposition(
            use_verlet_lists=False)
        self.system.cell_system.node_grid = self.original_node_grid
        self.system.time_step = 1e-3

    def tearDown(self):
        self.system.part.clear()

    def check_resort(self):
        n_part = 2351

        # Add the particles on node 0, so that they have to be resorted
        particles = self.system.part.add(
            pos=n_part * [(0, 0, 0)], type=n_part * [1])

        # And now change their positions
        particles.pos = self.system.box_l * \
            np.random.random((n_part, 3))

        # Add an interacting particle in a corner of the box
        self.system.part.add(pos=(0.01, 0.01, 0.01), type=0)
        if espressomd.has_features(['LENNARD_JONES']):
            self.system.non_bonded_inter[0, 1].lennard_jones.set_params(
                epsilon=1.0, sigma=3.0, cutoff=6.0, shift=0.1)
            ref_energy = self.system.analysis.energy()['total']
            assert ref_energy > 10.

        # Distribute the particles on the nodes
        part_dist = self.system.cell_system.resort()

        # Check that we did not lose particles
        self.assertEqual(sum(part_dist), n_part + 1)

        # Check that we can still access all the particles
        # This basically checks if part_node and local_particles
        # is still in a valid state after the particle exchange
        self.assertEqual(sum(self.system.part.all().type), n_part)

        # Check that the system is still valid
        if espressomd.has_features(['LENNARD_JONES']):
            # energy calculation
            new_energy = self.system.analysis.energy()['total']
            self.assertEqual(new_energy, ref_energy)
        # force calculation
        self.system.integrator.run(0, recalc_forces=True)

    def test_resort(self):
        self.check_resort()

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] != 4,
               "Skipping test: only runs for n_nodes >= 4")
    def test_resort_alternating(self):
        # check particle resorting when the left and right cells are different
        self.system.cell_system.node_grid = [4, 1, 1]
        self.check_resort()

    def test_position_rounding(self):
        """This places a particle on the box boundary,
           with parameters that could cause problems with
           rounding."""
        self.system.box_l = [50.0, 50.0, 50.0]
        self.system.cell_system.skin = 0.4
        self.system.min_global_cut = 12.0 / 4.25
        self.system.part.add(pos=[25, 25, 0])
        self.assertEqual(1, len(self.system.part))

    def test_00_fully_connected_boundary(self):
        system = self.system
        system.periodic = [True] * 3
        # Check fthat it's initially disabled
#        self.assertEqual(system.cell_system.get_params()[
#                         "fully_connected_boundary"], None)
#
#        # check setting and getting the parameter
#        system.cell_system.set_regular_decomposition(
#            fully_connected_boundary=(2, 1))
#        self.assertEqual(system.cell_system.get_params()[
#                         "fully_connected_boundary"], [2, 1])
#        # Check that the setting survies cell system re-intis
#        system.cell_system.min_global_cut = system.box_l / 4.1
#        self.assertEqual(system.cell_system.get_params()[
#                         "fully_connected_boundary"], [2, 1])
#
        system.cell_system.set_regular_decomposition()
        # Check that particle visibility
#        fc_normal = np.array((0, 0, 1))
#        fc_normal_idx = 2 
#        fc_dir = np.array((0, 1, 0))
        N = 3
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            sigma=1, epsilon=1, cutoff=system.box_l[0] / N + 0.01, shift="auto")
        indices = [np.array((i, j, k)) for i in range(N)
                   for j in range(N) for k in range(N)]

        def id_for_idx(idx): return (
            idx[0] % N) * 100 + (idx[1] % N) * 10 + idx[2] % N
        ids = [id_for_idx(idx) for idx in indices]
        dx = system.box_l / N
        positions = [idx * dx for idx in indices]
        max_range = system.box_l[0] / N
        epsilon = 1E-10
        print(system.cell_system.get_state())
        system.part.add(id=ids, pos=positions)

        def distance(id1, id2):
            return system.distance(
                system.part.by_id(id1), system.part.by_id(id2))
        must_find = [tuple(sorted([i[0], i[1]]))
                     for i in itertools.combinations(ids, 2)
                     if distance(*i) <= max_range + epsilon]
        can_find = [tuple(sorted([i[0], i[1]]))
                    for i in itertools.combinations(ids, 2)
                    if distance(*i) <= 2 * max_range * np.sqrt(3) + epsilon]
#        fully_connected_neighbors = \
#            [tuple(sorted([id_for_idx(idx), id_for_idx(idx - 1 * fc_normal + i * fc_dir)]))
# for idx in indices if idx[fc_normal_idx] == 0 for i in range(1, N)]
        system.cell_system.get_pairs(1)
        cs_pairs = system.cell_system.non_bonded_loop_trace()
        found = []
        for id1, id2, rest1, rest2, rest3, rest4 in cs_pairs:
            print(sorted([id1, id2]), rest1, rest2, rest3, rest4, system.distance(
                system.part.by_id(id1), system.part.by_id(id2))) 
            p = tuple(sorted((id1, id2)))
            assert p in can_find or p in must_find
            if p in must_find: must_find.remove(p)
            found.append(p)
        # Check for double counting of pairs

        assert len(found) == len(set(found))

        # check that all required pairs have been seen
        self.assertEqual(must_find, [])        


if __name__ == "__main__":
    ut.main()
