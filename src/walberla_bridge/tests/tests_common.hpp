#pragma once

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
#include "config.hpp"

#ifdef LB_WALBERLA

#include <LBWalberlaBase.hpp>
#include <LBWalberlaD3Q19FluctuatingMRT.hpp>
#include <LBWalberlaD3Q19MRT.hpp>
#include <lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <functional>
#include <memory>
#include <vector>

class LBTestParameters {
public:
  int seed;
  double agrid;
  double tau;
  double kT;
  double viscosity;
  double density;
  Utils::Vector3d box_dimensions;
  Utils::Vector3i grid_dimensions;
};

using LbGeneratorVector =
    std::vector<std::function<std::shared_ptr<LBWalberlaBase>(
        Utils::Vector3i, LBTestParameters)>>;

// Add all LBs with kT=0 to be tested, here
LbGeneratorVector unthermalized_lbs() {
  LbGeneratorVector lbs;
  // Unthermalized D3Q19 MRT
  lbs.push_back(
      [](const Utils::Vector3i mpi_shape, const LBTestParameters &params) {
        return std::make_shared<walberla::LBWalberlaD3Q19MRT>(
            params.viscosity, params.density, params.agrid, params.tau,
            params.box_dimensions, mpi_shape, 1);
      });

  // Thermalized D3Q19 MRT with kT set to 0
  lbs.push_back([](Utils::Vector3i mpi_shape, const LBTestParameters &params) {
    return std::make_shared<walberla::LBWalberlaD3Q19FluctuatingMRT>(
        params.viscosity, params.density, params.agrid, params.tau,
        params.box_dimensions, mpi_shape, 1, 0.0, params.seed);
  });
  return lbs;
}

// Add all LBs with thermalization to be tested, here
LbGeneratorVector thermalized_lbs() {
  LbGeneratorVector lbs;

  // Thermalized D3Q19 MRT with kT set to 0
  lbs.push_back(
      [](const Utils::Vector3i mpi_shape, const LBTestParameters &params) {
        return std::make_shared<walberla::LBWalberlaD3Q19FluctuatingMRT>(
            params.viscosity, params.density, params.agrid, params.tau,
            params.box_dimensions, mpi_shape, 1, params.kT, params.seed);
      });
  return lbs;
}

LbGeneratorVector all_lbs() {
  auto lbs = unthermalized_lbs();
  auto thermalized = thermalized_lbs();
  lbs.insert(lbs.end(), thermalized.begin(), thermalized.end());
  return lbs;
}

// Disable printing of type which does not support it
BOOST_TEST_DONT_PRINT_LOG_VALUE(LbGeneratorVector::value_type)

std::vector<Utils::Vector3i>
all_nodes_incl_ghosts(const Utils::Vector3i &grid_dimensions,
                      int n_ghost_layers) {
  std::vector<Utils::Vector3i> res;
  for (int x = -n_ghost_layers; x < grid_dimensions[0] + n_ghost_layers; x++) {
    for (int y = -n_ghost_layers; y < grid_dimensions[1] + n_ghost_layers;
         y++) {
      for (int z = -n_ghost_layers; z < grid_dimensions[2] + n_ghost_layers;
           z++) {
        res.push_back(Utils::Vector3i{x, y, z});
      }
    }
  }
  return res;
}

std::vector<Utils::Vector3i> local_nodes_incl_ghosts(
    std::pair<Utils::Vector3d, Utils::Vector3d> local_domain,
    int n_ghost_layers) {
  std::vector<Utils::Vector3i> res;
  auto const left = local_domain.first;
  auto const right = local_domain.second;
  for (int x = left[0] - n_ghost_layers; x < right[0] + n_ghost_layers; x++) {
    for (int y = left[1] - n_ghost_layers; y < right[1] + n_ghost_layers; y++) {
      for (int z = left[2] - n_ghost_layers; z < right[2] + n_ghost_layers;
           z++) {
        res.push_back(Utils::Vector3i{x, y, z});
      }
    }
  }
  return res;
}

#endif