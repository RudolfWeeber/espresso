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
 * Out-of-class VTK writer registration definition for
 * @ref walberla::LBWalberlaImplColorGradient.
 */

#include "LBVTKHelpers.hpp"

#include <memory>
#include <optional>
#include <string>

namespace walberla {

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplColorGradient<FloatType, Architecture>::
    register_vtk_field_writers(walberla::vtk::VTKOutput &vtk_obj,
                               LatticeModel::units_map const &units,
                               int flag_observables) {
  // The writers hold a host copy of the field, filled before each write.
  auto const add_scalar_writer = [&](BlockDataID const &field_id,
                                     std::string const &name,
                                     FloatType unit_conversion) {
    auto const blocks = m_lattice->get_blocks();
    WALBERLA_ASSERT_NOT_NULLPTR(blocks);
    auto scalar_writer =
        std::make_shared<ScalarVTKWriter<FloatType, ScalarField, float>>(
            field_id, name, unit_conversion);
    auto before_function = [blocks, field_id, scalar_writer]() {
      for (auto &block : *blocks) {
        auto *scalar_field = block.template getData<ScalarField>(field_id);
        auto const bci = scalar_field->xyzSize();
        scalar_writer->set_content(
            ek::accessor::Scalar::get(scalar_field, bci));
        scalar_writer->set_dims(
            Vector3<uint_t>(bci.xSize(), bci.ySize(), bci.zSize()));
      }
    };
    vtk_obj.addBeforeFunction(std::move(before_function));
    vtk_obj.addCellDataWriter(scalar_writer);
  };

  if (flag_observables & static_cast<int>(OutputVTK::phasefield)) {
    add_scalar_writer(m_phasefield_id, "phasefield", FloatType{1});
  }
  if (flag_observables & static_cast<int>(OutputVTK::density_a)) {
    add_scalar_writer(m_rho_field_id[0], "density_a",
                      FloatType_c(zero_centered_to_md(units.at("density"))));
  }
  if (flag_observables & static_cast<int>(OutputVTK::density_b)) {
    add_scalar_writer(m_rho_field_id[1], "density_b",
                      FloatType_c(zero_centered_to_md(units.at("density"))));
  }
  if (flag_observables & static_cast<int>(OutputVTK::velocity_vector)) {
    auto const unit_conversion = FloatType_c(units.at("velocity"));
    auto const blocks = m_lattice->get_blocks();
    WALBERLA_ASSERT_NOT_NULLPTR(blocks);
    auto velocity_writer =
        std::make_shared<VelocityVTKWriter<FloatType, float>>(
            m_velocity_field_id, "velocity_vector", unit_conversion);
    auto before_function = [this, blocks, velocity_writer]() {
      for (auto &block : *blocks) {
        auto *velocity_field =
            block.template getData<VectorField>(m_velocity_field_id);
        auto const bci = velocity_field->xyzSize();
        velocity_writer->set_content(
            lbm::accessor::Vector::get(velocity_field, bci));
        velocity_writer->set_dims(
            Vector3<uint_t>(bci.xSize(), bci.ySize(), bci.zSize()));
      }
    };
    vtk_obj.addBeforeFunction(std::move(before_function));
    vtk_obj.addCellDataWriter(velocity_writer);
  }
}

} // namespace walberla
