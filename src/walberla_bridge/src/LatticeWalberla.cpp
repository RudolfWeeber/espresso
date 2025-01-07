/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#include <walberla_bridge/BlockAndCell.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>

#include "utils/types_conversion.hpp"

#include <blockforest/Initialization.h>
#include <blockforest/StructuredBlockForest.h>
#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <domain_decomposition/IBlock.h>

#include <utils/Vector.hpp>

#include <cmath>
#include <initializer_list>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

LatticeWalberla::LatticeWalberla(Utils::Vector3i const &grid_dimensions,
                                 Utils::Vector3i const &node_grid,
                                 Utils::Vector3i const &block_grid,
                                 unsigned int n_ghost_layers)
    : m_grid_dimensions{grid_dimensions}, m_n_ghost_layers{n_ghost_layers} {
  using walberla::real_t;
  using walberla::uint_c;

  for (auto const i : {0u, 1u, 2u}) {
    if (m_grid_dimensions[i] % node_grid[i] != 0) {
      throw std::runtime_error(
          "Lattice grid dimensions and MPI node grid are not compatible.");
    }
    if (m_grid_dimensions[i] % block_grid[i] != 0) {
      throw std::runtime_error(
          "Lattice grid dimensions and block grid are not compatible.");
    }
  }

  auto constexpr lattice_constant = real_t{1};
  auto const cells_block = Utils::hadamard_division(grid_dimensions, block_grid);

  m_blocks = walberla::blockforest::createUniformBlockGrid(
      // number of blocks in each direction
      uint_c(block_grid[0]), uint_c(block_grid[1]), uint_c(block_grid[2]),
      // number of cells per block in each direction
      uint_c(cells_block[0]), uint_c(cells_block[1]), uint_c(cells_block[2]),
      lattice_constant,
      // number of cpus per direction
      uint_c(node_grid[0]), uint_c(node_grid[1]), uint_c(node_grid[2]),
      // periodicity
      true, true, true,
      // keep global block information
      false);
  for (IBlock &block : *m_blocks) {
    m_cached_blocks.push_back(&block);
  }
}

[[nodiscard]] std::pair<Utils::Vector3d, Utils::Vector3d>
LatticeWalberla::get_local_domain() const {
  using walberla::to_vector3d;
  // We allocate some blocks per mpi rank
  int64_t const stride_y = m_grid_dimensions[2];
  int64_t const stride_x = m_grid_dimensions[1]*stride_y;
  auto aa = m_blocks->begin()->getAABB();
  auto bb = m_blocks->begin()->getAABB();
  int64_t aa_index = stride_x*static_cast<int>(aa.min()[0]) + stride_y*static_cast<int>(aa.min()[1]) + static_cast<int>(aa.min()[2]);
  int64_t bb_index = stride_x*static_cast<int>(bb.max()[0]) + stride_y*static_cast<int>(bb.max()[1]) + static_cast<int>(bb.max()[2]);
  for (auto b = m_blocks->begin(); b != m_blocks->end(); ++b) {
    auto cc = b->getAABB();
    for (auto const i : {0u, 1u, 2u}) {
      if ((cc.max()[i] - cc.min()[i]) != 0) {
        assert(m_grid_dimensions[i] % static_cast<int>(cc.max()[i] - cc.min()[i]) == 0);
      }
    }
    int64_t min_index = stride_x*static_cast<int>(cc.min()[0]) + stride_y*static_cast<int>(cc.min()[1]) + static_cast<int>(cc.min()[2]);
    int64_t max_index = stride_x*static_cast<int>(cc.max()[0]) + stride_y*static_cast<int>(cc.max()[1]) + static_cast<int>(cc.max()[2]);
    if (min_index < aa_index) {
      aa = cc;
      aa_index = min_index;
    }
    if (max_index > bb_index) {
      bb = cc;
      bb_index = max_index;
    }
  }
  return {to_vector3d(aa.min()), to_vector3d(bb.max())};
}

[[nodiscard]] bool
LatticeWalberla::node_in_local_domain(Utils::Vector3i const &node) const {
  // Note: Lattice constant =1, cell centers offset by .5
  return ::walberla::get_block_and_cell(*this, node, false) != std::nullopt;
}
[[nodiscard]] bool
LatticeWalberla::node_in_local_halo(Utils::Vector3i const &node) const {
  return ::walberla::get_block_and_cell(*this, node, true) != std::nullopt;
}
[[nodiscard]] bool
LatticeWalberla::pos_in_local_domain(Utils::Vector3d const &pos) const {
  return ::walberla::get_block(*this, pos, false) != nullptr;
}
[[nodiscard]] bool
LatticeWalberla::pos_in_local_halo(Utils::Vector3d const &pos) const {
  return ::walberla::get_block(*this, pos, true) != nullptr;
}

[[nodiscard]] Utils::Vector3i
LatticeWalberla::calc_grid_dimensions(Utils::Vector3d const &box_size,
                                      double agrid) {
  auto const grid_dimensions =
      Utils::Vector3i{{static_cast<int>(std::round(box_size[0] / agrid)),
                       static_cast<int>(std::round(box_size[1] / agrid)),
                       static_cast<int>(std::round(box_size[2] / agrid))}};
  for (auto const i : {0u, 1u, 2u}) {
    if (std::abs(grid_dimensions[i] * agrid - box_size[i]) / box_size[i] >
        std::numeric_limits<double>::epsilon()) {
      throw std::runtime_error(
          "Box length not commensurate with agrid in direction " +
          std::to_string(i) + " length " + std::to_string(box_size[i]) +
          " agrid " + std::to_string(agrid));
    }
  }
  return grid_dimensions;
}
