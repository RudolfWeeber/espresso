/*
 * Copyright (C) 2019-2020 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#define BOOST_TEST_MODULE Walberla point force test
#define BOOST_TEST_DYN_LINK
#include "config.hpp"

#ifdef LB_WALBERLA

#define BOOST_TEST_NO_MAIN

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "tests_common.hpp"

#include "LBWalberlaImpl.hpp"
#include <LBWalberlaBase.hpp>
#include <lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>

#include <mpi.h>

#include <cmath>
#include <functional>
#include <iostream>
#include <math.h>
#include <memory>
#include <vector>

using Utils::hadamard_product;
using Utils::Vector3d;
using Utils::Vector3i;

namespace bdata = boost::unit_test::data;
constexpr const int offset = 0.;

Vector3i mpi_shape{};
BOOST_AUTO_TEST_CASE(test_interpolation) {
  double density = 1;
  double viscosity = 1. / 7.;
  auto lattice =
      std::make_shared<LatticeWalberla>(Vector3i{64, 64, 64}, mpi_shape, 1);
  auto lb = walberla::LBWalberlaImpl<double>(lattice, viscosity, density);
  auto le_pack = std::make_unique<LeesEdwardsPack>(
      0, 1, [=]() { return offset; }, []() { return 0.0; });
  lb.set_collision_model(std::move(le_pack));

  auto const grid_size_y = lattice->get_grid_dimensions()[1];

  auto const force_pos = Vector3d{32., double(grid_size_y), 32.};
  auto const force_node = Vector3i{32, grid_size_y, 32};
  auto const f1 = Vector3d{0.3, -0.2, 0.3};
  lb.add_force_at_pos(force_pos, f1);

  lb.integrate();

  auto const node = Vector3i{force_node[0] - offset, -1, 32};
//  for (auto const &n : all_nodes_incl_ghosts(lb.lattice())) {
//    if (lb.lattice().node_in_local_halo(node)) {
  auto const laf = *(lb.get_node_last_applied_force(node));
  printf("laf: %f %f %f",laf[0], laf[1], laf[2]);
  BOOST_CHECK_SMALL((laf - f1).norm(), 1E-10);
//    }
//  }
}

int main(int argc, char **argv) {
  int n_nodes;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());
  walberla_mpi_init();

  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#else // ifdef LB_WALBERLA
int main(int argc, char **argv) {}
#endif
