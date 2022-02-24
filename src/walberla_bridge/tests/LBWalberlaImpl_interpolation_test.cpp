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
constexpr const int offset = 10;

Vector3i mpi_shape{};
BOOST_AUTO_TEST_CASE(test_interpolation_force) {
  double density = 1;
  double viscosity = 1. / 7.;
  auto lattice =
      std::make_shared<LatticeWalberla>(Vector3i{64, 64, 64}, mpi_shape, 1);
  auto lb = walberla::LBWalberlaImpl<double>(lattice, viscosity, density);
  auto le_pack = std::make_unique<LeesEdwardsPack>(
      0, 1, [=]() { return offset; }, []() { return 0.0; });
  lb.set_collision_model(std::move(le_pack));

  auto const grid_size_y = lattice->get_grid_dimensions()[1];

  auto const force_pos = Vector3d{32+0.5, double(grid_size_y)-0.5, 32.+0.5};
  auto const force_node = Vector3i{32, 63, 32};
  auto const f1 = Vector3d{0.3, -0.2, 0.3};
  lb.add_force_at_pos(force_pos, f1);

  lb.integrate();

  auto const node = Vector3i{force_node[0] - offset, -1, 32};
  auto const laf = *(lb.get_node_last_applied_force(node, true));
  BOOST_CHECK_SMALL((laf - f1).norm(), 1E-10);
}

BOOST_AUTO_TEST_CASE(test_interpolation_velocity) {
  double density = 1;
  double viscosity = 1. / 7.;
  auto lattice =
      std::make_shared<LatticeWalberla>(Vector3i{64, 64, 64}, mpi_shape, 1);
  auto lb = walberla::LBWalberlaImpl<double>(lattice, viscosity, density);
  auto le_pack = std::make_unique<LeesEdwardsPack>(
      0, 1, [=]() { return offset; }, []() { return 0.0; });
  lb.set_collision_model(std::move(le_pack));

  auto const node_up = Vector3i{32, 63, 32};
  
  auto const node_up_neighbor_1 = Vector3i{32-1, 63, 32};
  auto const node_up_neighbor_2 = Vector3i{32, 63, 32-1};
  auto const node_up_neighbor_3 = Vector3i{32+1, 63, 32};
  auto const node_up_neighbor_4 = Vector3i{32, 63, 32+1};

  auto const vel = Vector3d{0.05, 0.1, 0.15};
  lb.set_node_velocity(node_up, vel);

  lb.integrate();
 
  //This is the original node
  //auto const node_down = Vector3i{32 - offset, -1, 32};
  
  auto const node_down_neighbor_1 = Vector3i{32-1-offset, -1, 32};
  auto const node_down_neighbor_2 = Vector3i{32-offset, -1, 32-1};
  auto const node_down_neighbor_3 = Vector3i{32+1-offset, -1, 32};
  auto const node_down_neighbor_4 = Vector3i{32-offset, -1, 32+1};
  
  auto const vel_up_1 = *(lb.get_node_velocity(node_up_neighbor_1));
  auto const vel_up_2 = *(lb.get_node_velocity(node_up_neighbor_2));
  auto const vel_up_3 = *(lb.get_node_velocity(node_up_neighbor_3));
  auto const vel_up_4 = *(lb.get_node_velocity(node_up_neighbor_4));
  
  auto const vel_down_1 = *(lb.get_node_velocity(node_down_neighbor_1));
  auto const vel_down_2 = *(lb.get_node_velocity(node_down_neighbor_2));
  auto const vel_down_3 = *(lb.get_node_velocity(node_down_neighbor_3));
  auto const vel_down_4 = *(lb.get_node_velocity(node_down_neighbor_4));
 
  //printf("node_up %f %f %f \n", vel_up2[0], vel_up2[1], vel_up2[2]);
  //printf("node_down %f %f %f \n", vel_down[0], vel_down[1], vel_down[2]);

  BOOST_CHECK_SMALL((vel_up_1 - vel_down_1).norm(), 1E-10);
  BOOST_CHECK_SMALL((vel_up_2 - vel_down_2).norm(), 1E-10);
  BOOST_CHECK_SMALL((vel_up_3 - vel_down_3).norm(), 1E-10);
  BOOST_CHECK_SMALL((vel_up_4 - vel_down_4).norm(), 1E-10);
}

BOOST_AUTO_TEST_CASE(test_interpolation_pdfs) {
  double density = 1;
  double viscosity = 1. / 7.;
  auto lattice =
      std::make_shared<LatticeWalberla>(Vector3i{64, 64, 64}, mpi_shape, 1);
  auto lb = walberla::LBWalberlaImpl<double>(lattice, viscosity, density);
  auto le_pack = std::make_unique<LeesEdwardsPack>(
      0, 1, [=]() { return 0.0; }, []() { return 0.0; });
  lb.set_collision_model(std::move(le_pack));

  auto const node_up = Vector3i{32, 63, 32};
  auto const pop = std::vector<double>{0.00, 0.01, 0.02, 
                                       0.03, 0.04, 0.05,
                                       0.06, 0.07, 0.08,
                                       0.09, 0.10, 0.11,
                                       0.12, 0.13, 0.14, 
                                       0.15, 0.16, 0.17,
                                       0.18};

  lb.set_node_pop(node_up, pop);
  
  lb.integrate();
 
  auto const node_down_unshifted = Vector3i{32, 63, 32};
  auto const pop_unshifted = *(lb.get_node_pop(node_down_unshifted));

  le_pack = std::make_unique<LeesEdwardsPack>(
      0, 1, [=]() { return offset; }, []() { return 0.0; });
  lb.set_collision_model(std::move(le_pack));
  
  auto const node_up_2 = Vector3i{25, 63, 32};
  lb.set_node_pop(node_up_2, pop);
  
  lb.integrate();
 
  auto const node_down_shifted = Vector3i{25 - offset, -1, 32};
  auto const pop_shifted = *(lb.get_node_pop(node_down_shifted, true));

//  printf("unshifted: %f", pop_unshifted[0]);
//  printf("shifted: %f", pop_shifted[0]);
  
  for (int i = 0; i < pop_unshifted.size(); i++){
    BOOST_CHECK_SMALL((pop_shifted[i] - pop_unshifted[i]), 1E-10);
  }
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
