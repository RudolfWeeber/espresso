/*
 * Copyright (C) 2019-2026 The ESPResSo project
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

#pragma once

/**
 * @file
 * Out-of-class position-based interpolation definitions for
 * @ref walberla::LBWalberlaImplSingleComponent. Read-side helpers and the
 * velocity/density batch accessors have moved into
 * @ref walberla::LBWalberlaCommon; this file retains the leaf-local
 * force-deposition path (layout-dependent).
 */

#include <utils/Vector.hpp>
#include <utils/interpolation/bspline_3d.hpp>
#include <utils/interpolation/bspline_3d_gradient.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace walberla {

/**
 * @brief Distribute forces to the lattice at given positions.
 * Uses B-spline interpolation to spread each force over the surrounding
 * lattice nodes. On GPU, positions are transformed to block-local
 * coordinates and the operation is performed in a single kernel launch.
 */
template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplSingleComponent<FloatType, Architecture>::add_forces_at_pos(
    std::vector<Utils::Vector3d> const &pos,
    std::vector<Utils::Vector3d> const &forces) {
  assert(pos.size() == forces.size());
  if (pos.empty()) {
    return;
  }
  if constexpr (Architecture == lbmpy::Arch::CPU) {
    auto const kernel = make_force_interpolation_kernel();
    for (std::size_t i = 0ul; i < pos.size(); ++i) {
      kernel(pos[i], forces[i]);
    }
  }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  if constexpr (Architecture == lbmpy::Arch::GPU) {
    auto const &lattice = get_lattice();
    auto const &block = *(lattice.get_blocks()->begin());
    auto const origin = block.getAABB().min();
    std::vector<FloatType> host_pos;
    std::vector<FloatType> host_force;
    host_pos.reserve(3ul * pos.size());
    host_force.reserve(3ul * forces.size());
    assert(lattice.get_blocks()->getNumberOfBlocks() == 1u);
    for (auto const &vec : pos) {
#pragma unroll
      for (std::size_t i : {0ul, 1ul, 2ul}) {
        host_pos.emplace_back(static_cast<FloatType>(vec[i] - origin[i]));
      }
    }
    for (auto const &vec : forces) {
#pragma unroll
      for (std::size_t i : {0ul, 1ul, 2ul}) {
        host_force.emplace_back(static_cast<FloatType>(vec[i]));
      }
    }
    zero_centered_to_lb_in_place(host_force);
    auto const gl = lattice.get_ghost_layers();
    auto field = block.template uncheckedFastGetData<VectorField>(
        m_force_to_be_applied_id);
    lbm::accessor::Interpolation::add_force(field, host_pos, host_force, gl);
  }
#endif
}

template <typename FloatType, lbmpy::Arch Architecture>
auto LBWalberlaImplSingleComponent<
    FloatType, Architecture>::make_force_interpolation_kernel() const {
  auto const &lattice = *m_lattice;
  auto const &blocks = *lattice.get_blocks();
  assert(lattice.get_ghost_layers() == 1u);
  return [&](Utils::Vector3d const &pos, Utils::Vector3d const &force) {
    if (not get_block_extended(lattice, pos, 1u)) {
      return;
    }
    interpolate_bspline_at_pos(
        pos, [&, conv = m_zc_to_lb, field_id = m_force_to_be_applied_id](
                 std::array<int, 3> const node, double weight) {
          auto block = get_block_extended(lattice, node, 0u);
          if (!block)
            block = get_block_extended(lattice, node, 1u);
          if (block) {
            auto cell = to_cell(node);
            blocks.transformGlobalToBlockLocalCell(cell, *block);
            weight *= conv;
            auto const weighted_force = to_vector3<FloatType>(weight * force);
            auto field =
                block->template uncheckedFastGetData<VectorField>(field_id);
            lbm::accessor::Vector::add(field, weighted_force, cell);
          }
        });
  };
}

/**
 * @brief add_density_weighted_forces_at_pos is not implemented for
 * single-component LB. Throws std::runtime_error.
 */
template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplSingleComponent<FloatType, Architecture>::
    add_density_weighted_forces_at_pos(
        std::vector<Utils::Vector3d> const &pos,
        std::vector<Utils::Vector3d> const &forces) {
  (void)pos;
  (void)forces;
  throw std::runtime_error(
      "add_density_weighted_forces_at_pos is not implemented for "
      "single-component LB");
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImplSingleComponent<FloatType, Architecture>::add_force_at_pos(
    Utils::Vector3d const &pos, Utils::Vector3d const &force) {
  if (!m_lattice->pos_in_local_halo(pos))
    return false;
  auto const kernel = make_force_interpolation_kernel();
  kernel(pos, force);
  return true;
}

} // namespace walberla
