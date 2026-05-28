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
 * @ref walberla::LBWalberlaImplColorGradient.
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
 * lattice nodes.
 */
template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplColorGradient<FloatType, Architecture>::add_forces_at_pos(
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
auto LBWalberlaImplColorGradient<
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
 * @brief Distribute forces density-weighted to the two lattices
 * for the two components at given positions.
 * Uses B-spline interpolation to spread each force over the surrounding
 * lattice nodes. GPU not implemented.
 */
template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplColorGradient<FloatType, Architecture>::
    add_density_weighted_forces_at_pos(
        std::vector<Utils::Vector3d> const &pos,
        std::vector<Utils::Vector3d> const &forces) {
  assert(pos.size() == forces.size());
  if (pos.empty()) {
    return;
  }
  if constexpr (Architecture == lbmpy::Arch::CPU) {
    auto const kernel = make_density_weighted_force_interpolation_kernel();
    for (std::size_t i = 0ul; i < pos.size(); ++i) {
      kernel(pos[i], forces[i]);
    }
  }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  if constexpr (Architecture == lbmpy::Arch::GPU) {
    throw std::runtime_error(
        "Density-weighted force interpolation not implemented on GPU");
  }
#endif
}

/**
 * @brief Distribute solvation force on the two lattices
 * for the two components at given positions.
 * GPU not implemented.
 */
template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplColorGradient<FloatType, Architecture>::
    add_solvation_forces_at_pos(std::vector<Utils::Vector3d> const &pos,
                                std::vector<double> const &delta_mus) {
  assert(pos.size() == delta_mus.size());
  if (pos.empty())
    return;
  if constexpr (Architecture == lbmpy::Arch::CPU) {
    auto const kernel = make_solvation_force_interpolation_kernel();
    for (std::size_t i = 0ul; i < pos.size(); ++i) {
      kernel(pos[i], delta_mus[i]);
    }
  }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  if constexpr (Architecture == lbmpy::Arch::GPU) {
    throw std::runtime_error(
        "Solvation force interpolation not implemented on GPU");
  }
#endif
}

template <typename FloatType, lbmpy::Arch Architecture>
auto LBWalberlaImplColorGradient<
    FloatType, Architecture>::make_solvation_force_interpolation_kernel()
    /* Only for two-component LB, adds solvation force to fluid*/
    const {
  auto const &lattice = *m_lattice;
  auto const &blocks = *lattice.get_blocks();
  assert(lattice.get_ghost_layers() == 1u);
  return [&](Utils::Vector3d const &pos, double delta_mu) {
    if (not get_block_extended(lattice, pos, 1u)) {
      return;
    }
    auto const grid_spacing = Utils::Vector3d::broadcast(1.);
    auto const offset = Utils::Vector3d::broadcast(.5);
    /* Accumulate grad(rho_a) and grad(rho_b) at particle position using
     * bspline_3d_gradient*/
    Utils::Vector3d grad_rho_a{}, grad_rho_b{};
    Utils::Interpolation::bspline_3d_gradient<2>(
        pos,
        [&](std::array<int, 3> const node, Utils::Vector3d const &grad_weight) {
          auto block = get_block_extended(lattice, node, 1u);
          if (!block)
            return;
          auto cell = to_cell(node);
          blocks.transformGlobalToBlockLocalCell(cell, *block);
          auto const cell_rho_a = static_cast<double>(
              block
                  ->template uncheckedFastGetData<ScalarField>(
                      m_rho_field_id[0])
                  ->get(cell));
          auto const cell_rho_b = static_cast<double>(
              block
                  ->template uncheckedFastGetData<ScalarField>(
                      m_rho_field_id[1])
                  ->get(cell));
          grad_rho_a += cell_rho_a * grad_weight;
          grad_rho_b += cell_rho_b * grad_weight;
        },
        grid_spacing, offset);

    /* Accumulate rho_a and rho_b at particle position using bspline_3d*/
    double rho_a{}, rho_b{};
    Utils::Interpolation::bspline_3d<2>(
        pos,
        [&](std::array<int, 3> const node, double weight) {
          auto block = get_block_extended(lattice, node, 1u);
          if (!block)
            return;
          auto cell = to_cell(node);
          blocks.transformGlobalToBlockLocalCell(cell, *block);
          auto const cell_rho_a = static_cast<double>(
              block
                  ->template uncheckedFastGetData<ScalarField>(
                      m_rho_field_id[0])
                  ->get(cell));
          auto const cell_rho_b = static_cast<double>(
              block
                  ->template uncheckedFastGetData<ScalarField>(
                      m_rho_field_id[1])
                  ->get(cell));
          rho_a += cell_rho_a * weight;
          rho_b += cell_rho_b * weight;
        },
        grid_spacing, offset);

    /*spread force back onto fluid using bspline_3d*/
    interpolate_bspline_at_pos(pos, [&](std::array<int, 3> const node,
                                        double weight) {
      auto block = get_block_extended(lattice, node, 0u);
      if (!block)
        block = get_block_extended(lattice, node, 1u);
      if (!block)
        return;
      auto cell = to_cell(node);
      blocks.transformGlobalToBlockLocalCell(cell, *block);
      auto const total_rho = rho_a + rho_b;
      if (total_rho == 0.)
        /* No fluid to couple to */
        return;
      auto const inv_rho_sq = 1. / (total_rho * total_rho);
      // f^a = +delta_mu * rho_b * inv_rho_sq * grad_rho_a * weight
      // f^b = -delta_mu * rho_a * inv_rho_sq * grad_rho_b * weight
      auto const f_a = to_vector3<FloatType>(+delta_mu * rho_b * inv_rho_sq *
                                             grad_rho_a * weight);
      auto const f_b = to_vector3<FloatType>(-delta_mu * rho_a * inv_rho_sq *
                                             grad_rho_b * weight);
      auto field_a = block->template uncheckedFastGetData<VectorField>(
          m_force_color_gradient_field_id[0]);
      auto field_b = block->template uncheckedFastGetData<VectorField>(
          m_force_color_gradient_field_id[1]);
      lbm::accessor::Vector::add(field_a, f_a, cell);
      lbm::accessor::Vector::add(field_b, f_b, cell);
    });
  };
}

template <typename FloatType, lbmpy::Arch Architecture>
auto LBWalberlaImplColorGradient<
    FloatType, Architecture>::make_density_weighted_force_interpolation_kernel()
    /* Only for two-component LB, e.g. friction coupling*/
    const {
  auto const &lattice = *m_lattice;
  auto const &blocks = *lattice.get_blocks();
  assert(lattice.get_ghost_layers() == 1u);
  return [&](Utils::Vector3d const &pos, Utils::Vector3d const &force) {
    if (not get_block_extended(lattice, pos, 1u)) {
      return;
    }
    auto const grid_spacing = Utils::Vector3d::broadcast(1.);
    auto const offset = Utils::Vector3d::broadcast(.5);
    /* Accumulate rho_a and rho_b at particle position using bspline_3d*/
    double rho_a{}, rho_b{};
    Utils::Interpolation::bspline_3d<2>(
        pos,
        [&](std::array<int, 3> const node, double weight) {
          auto block = get_block_extended(lattice, node, 1u);
          if (!block)
            return;
          auto cell = to_cell(node);
          blocks.transformGlobalToBlockLocalCell(cell, *block);
          auto const cell_rho_a = static_cast<double>(
              block
                  ->template uncheckedFastGetData<ScalarField>(
                      m_rho_field_id[0])
                  ->get(cell));
          auto const cell_rho_b = static_cast<double>(
              block
                  ->template uncheckedFastGetData<ScalarField>(
                      m_rho_field_id[1])
                  ->get(cell));
          rho_a += cell_rho_a * weight;
          rho_b += cell_rho_b * weight;
        },
        grid_spacing, offset);

    /*spread force back onto fluid using bspline_3d*/
    interpolate_bspline_at_pos(
        pos, [&](std::array<int, 3> const node, double weight) {
          /* no conversion factor m_zc_to_lb since densities for two-component
           * lbs are not saved zero_centered*/
          auto block = get_block_extended(lattice, node, 0u);
          if (!block)
            block = get_block_extended(lattice, node, 1u);
          if (block) {
            auto cell = to_cell(node);
            blocks.transformGlobalToBlockLocalCell(cell, *block);
            auto const total_rho = rho_a + rho_b;
            if (total_rho == 0.)
              /* No fluid to couple to */
              return;
            auto const inv_rho = 1. / total_rho;
            auto const weighted_force_a =
                to_vector3<FloatType>(weight * (rho_a * inv_rho) * force);
            auto const weighted_force_b =
                to_vector3<FloatType>(weight * (rho_b * inv_rho) * force);
            auto field_a = block->template uncheckedFastGetData<VectorField>(
                m_force_color_gradient_field_id[0]);
            auto field_b = block->template uncheckedFastGetData<VectorField>(
                m_force_color_gradient_field_id[1]);
            lbm::accessor::Vector::add(field_a, weighted_force_a, cell);
            lbm::accessor::Vector::add(field_b, weighted_force_b, cell);
          }
        });
  };
}

template <typename FloatType, lbmpy::Arch Architecture>
auto LBWalberlaImplColorGradient<
    FloatType, Architecture>::make_color_gradient_interpolation_kernel() const {
  auto const &lattice = *m_lattice;
  auto const &blocks = *lattice.get_blocks();
  assert(lattice.get_ghost_layers() == 1u);
  return [&](Utils::Vector3d const &pos) {
    Utils::Vector3d acc{0., 0., 0.};
    interpolate_bspline_at_pos(
        pos, [&, field_id = m_color_gradient_field_id](
                 std::array<int, 3> const node, double weight) {
          // Nodes with zero weight might not be accessible, because they can be
          // outside ghost layers
          if (weight != 0.) {
            auto block = get_block_extended(lattice, node, 1u);
            if (!block)
              throw interpolation_illegal_access("color gradient", pos, node,
                                                 weight);
            auto cell = to_cell(node);
            blocks.transformGlobalToBlockLocalCell(cell, *block);
            auto field =
                block->template uncheckedFastGetData<VectorField>(field_id);
            auto const cg = lbm::accessor::Vector::get(field, cell);
            acc += to_vector3d(cg) * weight;
          }
        });
    return acc;
  };
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<Utils::Vector3d>
LBWalberlaImplColorGradient<FloatType, Architecture>::
    get_color_gradients_at_pos(std::vector<Utils::Vector3d> const &pos) {
  if (pos.empty()) {
    return {};
  }
  std::vector<Utils::Vector3d> color_gradient{};
  color_gradient.reserve(pos.size());
  if constexpr (Architecture == lbmpy::Arch::CPU) {
    auto const kernel = make_color_gradient_interpolation_kernel();
    std::ranges::transform(pos, std::back_inserter(color_gradient), kernel);
  }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  if constexpr (Architecture == lbmpy::Arch::GPU) {
    throw std::runtime_error(
        "Density-weighted force interpolation not implemented on GPU");
  }
#endif
  return color_gradient;
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImplColorGradient<FloatType, Architecture>::add_force_at_pos(
    Utils::Vector3d const &pos, Utils::Vector3d const &force) {
  throw std::runtime_error(
      "add_force_at_pos is not implemented for two-component LB");
}

} // namespace walberla
