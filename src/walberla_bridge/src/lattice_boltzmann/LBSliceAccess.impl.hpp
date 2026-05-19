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
 * Out-of-class slice access definitions for
 * @ref walberla::LBWalberlaImpl.
 */

#include <utils/Vector.hpp>

#include <span>
#include <vector>

namespace walberla {

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double> LBWalberlaImpl<FloatType, Architecture>::get_slice_velocity(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(3u * ci.numCells());
        auto const field =
            block.template getData<VectorField>(m_velocity_field_id);
        auto values = lbm::accessor::Vector::get(field, bci);

        auto kernel = [&values, &out, this](unsigned const block_index,
                                            unsigned const local_index,
                                            Utils::Vector3i const &node) {
          if (m_boundary->node_is_boundary(node)) {
            auto const &vec = m_boundary->get_node_value_at_boundary(node);
            for (uint_t f = 0u; f < 3u; ++f) {
              out[3u * local_index + f] = vec[f];
            }
          } else {
            for (uint_t f = 0u; f < 3u; ++f) {
              out[3u * local_index + f] =
                  double_c(values[3u * block_index + f]);
            }
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::set_slice_velocity(
    Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
    std::vector<double> const &velocity) {
  if (has_two_components()) {
    throw std::runtime_error(
        "set_slice_velocity is not supported for two-component LB. "
        "Set densities and populations instead to control the barycentric "
        "velocity.");
  }
  m_pending_ghost_comm.set(GhostComm::PDF);
  m_pending_ghost_comm.set(GhostComm::VEL);
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        assert(velocity.size() == 3u * ci.numCells());
        auto pdf_field = block.template getData<PdfField>(m_pdf_field_id[0]);
        auto force_field =
            block.template getData<VectorField>(m_last_applied_force_field_id);
        auto vel_field =
            block.template getData<VectorField>(m_velocity_field_id);
        std::vector<FloatType> values(3u * bci.numCells());

        auto kernel = [&values, &velocity](unsigned const block_index,
                                           unsigned const local_index,
                                           Utils::Vector3i const &) {
          for (uint_t f = 0u; f < 3u; ++f) {
            values[3u * block_index + f] =
                numeric_cast<FloatType>(velocity[3u * local_index + f]);
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        lbm::accessor::Velocity::set(pdf_field, vel_field, force_field, values,
                                     bci);
      });
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double>
LBWalberlaImpl<FloatType, Architecture>::get_slice_last_applied_force(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(3u * ci.numCells());
        auto const field =
            block.template getData<VectorField>(m_last_applied_force_field_id);
        auto const values = lbm::accessor::Vector::get(field, bci);

        auto kernel = [&values, &out](unsigned const block_index,
                                      unsigned const local_index,
                                      Utils::Vector3i const &) {
          for (uint_t f = 0u; f < 3u; ++f) {
            out[3u * local_index + f] = values[3u * block_index + f];
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  zero_centered_to_md_in_place(out);
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::set_slice_last_applied_force(
    Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
    std::vector<double> const &force) {
  m_pending_ghost_comm.set(GhostComm::VEL);
  m_pending_ghost_comm.set(GhostComm::LAF);
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        assert(force.size() == 3u * ci.numCells());
        auto pdf_field = block.template getData<PdfField>(m_pdf_field_id[0]);
        auto force_field =
            block.template getData<VectorField>(m_last_applied_force_field_id);
        auto vel_field =
            block.template getData<VectorField>(m_velocity_field_id);
        std::vector<FloatType> values(3u * bci.numCells());

        auto kernel = [&values, &force](unsigned const block_index,
                                        unsigned const local_index,
                                        Utils::Vector3i const &) {
          for (uint_t f = 0u; f < 3u; ++f) {
            values[3u * block_index + f] =
                numeric_cast<FloatType>(force[3u * local_index + f]);
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        lbm::accessor::Force::set(pdf_field, vel_field, force_field, values,
                                  m_density, bci);
      });
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double>
LBWalberlaImpl<FloatType, Architecture>::get_slice_population(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  auto const pop_per_node =
      has_two_components() ? 2u * this->stencil_size() : this->stencil_size();
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(pop_per_node * ci.numCells());

        auto const pdf_field_a =
            block.template getData<PdfField>(m_pdf_field_id[0]);
        auto const values_a = lbm::accessor::Population::get(pdf_field_a, bci);

        if (has_two_components()) {
          auto const pdf_field_b =
              block.template getData<PdfField>(m_pdf_field_id[1]);
          auto const values_b =
              lbm::accessor::Population::get(pdf_field_b, bci);

          auto kernel = [&values_a, &values_b, &out,
                         this](unsigned const block_index,
                               unsigned const local_index,
                               Utils::Vector3i const &) {
            auto const ss = this->stencil_size();
            for (uint_t f = 0u; f < ss; ++f) {
              out[2u * ss * local_index + f] = values_a[ss * block_index + f];
              out[2u * ss * local_index + ss + f] =
                  values_b[ss * block_index + f];
            }
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        } else {
          auto kernel = [&values_a, &out, this](unsigned const block_index,
                                                unsigned const local_index,
                                                Utils::Vector3i const &) {
            for (uint_t f = 0u; f < this->stencil_size(); ++f) {
              out[this->stencil_size() * local_index + f] =
                  values_a[this->stencil_size() * block_index + f];
            }
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        }
      });
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::set_slice_population(
    Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
    std::vector<double> const &population) {
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        auto const pop_per_node = has_two_components()
                                      ? 2u * this->stencil_size()
                                      : this->stencil_size();
        assert(population.size() == pop_per_node * ci.numCells());
        (void)pop_per_node;
        auto pdf_field_a = block.template getData<PdfField>(m_pdf_field_id[0]);
        auto force_field =
            block.template getData<VectorField>(m_last_applied_force_field_id);
        auto vel_field =
            block.template getData<VectorField>(m_velocity_field_id);
        std::vector<FloatType> values_a(this->stencil_size() * bci.numCells());

        if (has_two_components()) {
          auto pdf_field_b =
              block.template getData<PdfField>(m_pdf_field_id[1]);
          std::vector<FloatType> values_b(this->stencil_size() *
                                          bci.numCells());

          auto kernel = [&values_a, &values_b, &population,
                         this](unsigned const block_index,
                               unsigned const local_index,
                               Utils::Vector3i const &) {
            auto const ss = this->stencil_size();
            for (uint_t f = 0u; f < ss; ++f) {
              values_a[ss * block_index + f] = numeric_cast<FloatType>(
                  population[2u * ss * local_index + f]);
              values_b[ss * block_index + f] = numeric_cast<FloatType>(
                  population[2u * ss * local_index + ss + f]);
            }
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
          lbm::accessor::Population::set(pdf_field_a, vel_field, force_field,
                                         values_a, bci);
          lbm::accessor::Population::set(pdf_field_b, vel_field, force_field,
                                         values_b, bci);
        } else {
          auto kernel = [&values_a, &population,
                         this](unsigned const block_index,
                               unsigned const local_index,
                               Utils::Vector3i const &) {
            for (uint_t f = 0u; f < this->stencil_size(); ++f) {
              values_a[this->stencil_size() * block_index + f] =
                  numeric_cast<FloatType>(
                      population[this->stencil_size() * local_index + f]);
            }
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
          lbm::accessor::Population::set(pdf_field_a, vel_field, force_field,
                                         values_a, bci);
        }
      });
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double> LBWalberlaImpl<FloatType, Architecture>::get_slice_density(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  auto const n_components = has_two_components() ? 2u : 1u;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(n_components * ci.numCells());

        if (has_two_components()) {
          auto const rho_a_field =
              block.template getData<ScalarField>(m_rho_field_id[0]);
          auto const rho_b_field =
              block.template getData<ScalarField>(m_rho_field_id[1]);
          // Pre-collect scalar field values into flat arrays
          std::vector<double> values_a(bci.numCells());
          std::vector<double> values_b(bci.numCells());
          unsigned idx = 0u;
          for (auto x = bci.xMin(); x <= bci.xMax(); ++x) {
            for (auto y = bci.yMin(); y <= bci.yMax(); ++y) {
              for (auto z = bci.zMin(); z <= bci.zMax(); ++z) {
                values_a[idx] = double_c(rho_a_field->get(x, y, z));
                values_b[idx] = double_c(rho_b_field->get(x, y, z));
                ++idx;
              }
            }
          }

          auto kernel = [&values_a, &values_b, &out](unsigned const block_index,
                                                     unsigned const local_index,
                                                     Utils::Vector3i const &) {
            out[2u * local_index + 0u] = values_a[block_index];
            out[2u * local_index + 1u] = values_b[block_index];
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        } else {
          auto const pdf_field =
              block.template getData<PdfField>(m_pdf_field_id[0]);
          auto const values =
              lbm::accessor::Density::get(pdf_field, m_density, bci);

          auto kernel = [&values, &out](unsigned const block_index,
                                        unsigned const local_index,
                                        Utils::Vector3i const &) {
            out[local_index] = values[block_index];
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        }
      });
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::set_slice_density(
    Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
    std::vector<double> const &density) {
  m_pending_ghost_comm.set(GhostComm::PDF);
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        auto const n_components = has_two_components() ? 2u : 1u;
        assert(density.size() == n_components * ci.numCells());
        (void)n_components;

        if (has_two_components()) {
          auto rho_a_field =
              block.template getData<ScalarField>(m_rho_field_id[0]);
          auto rho_b_field =
              block.template getData<ScalarField>(m_rho_field_id[1]);
          // Collect values from input into block-local arrays
          std::vector<FloatType> values_a(bci.numCells());
          std::vector<FloatType> values_b(bci.numCells());

          auto kernel = [&values_a, &values_b,
                         &density](unsigned const block_index,
                                   unsigned const local_index,
                                   Utils::Vector3i const &) {
            values_a[block_index] =
                numeric_cast<FloatType>(density[2u * local_index + 0u]);
            values_b[block_index] =
                numeric_cast<FloatType>(density[2u * local_index + 1u]);
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);

          // Write to scalar fields
          unsigned idx = 0u;
          for (auto x = bci.xMin(); x <= bci.xMax(); ++x) {
            for (auto y = bci.yMin(); y <= bci.yMax(); ++y) {
              for (auto z = bci.zMin(); z <= bci.zMax(); ++z) {
                rho_a_field->get(x, y, z) = values_a[idx];
                rho_b_field->get(x, y, z) = values_b[idx];
                ++idx;
              }
            }
          }
        } else {
          auto pdf_field = block.template getData<PdfField>(m_pdf_field_id[0]);
          std::vector<FloatType> values(bci.numCells());

          auto kernel = [&values, &density](unsigned const block_index,
                                            unsigned const local_index,
                                            Utils::Vector3i const &) {
            values[block_index] = numeric_cast<FloatType>(density[local_index]);
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
          lbm::accessor::Density::set(pdf_field, values, m_density, bci);
        }
      });
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double>
LBWalberlaImpl<FloatType, Architecture>::get_slice_pressure_tensor(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  if (has_two_components()) {
    throw std::runtime_error(
        "pressure tensor not implemented for two-component LB");
  }
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(9u * ci.numCells());
        auto const pdf_field =
            block.template getData<PdfField>(m_pdf_field_id[0]);
        auto values =
            lbm::accessor::PressureTensor::get(pdf_field, m_density, bci);

        auto kernel = [&values, &out, this](unsigned const block_index,
                                            unsigned const local_index,
                                            Utils::Vector3i const &) {
          pressure_tensor_correction(
              std::span<FloatType, 9ul>(&values[9u * block_index], 9ul));
          for (uint_t f = 0u; f < 9u; ++f) {
            out[9u * local_index + f] = values[9u * block_index + f];
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

} // namespace walberla
