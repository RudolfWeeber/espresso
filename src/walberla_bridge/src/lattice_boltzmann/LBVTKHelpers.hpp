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
 * VTK writer helper class templates shared by LBWalberlaImpl and
 * LBWalberlaImplSingleComponent.
 */

#include <memory>
#include <string>

namespace walberla {

/**
 * @brief Base class for LB field VTK writers.
 * Provides unit conversion and field access for cell-based VTK output.
 * On GPU builds, the GPU field is copied to a CPU mirror before writing.
 * @tparam FloatType   Internal LB precision (float or double).
 * @tparam Field_T     waLBerla field type to read from.
 * @tparam F_SIZE_ARG  Number of components per cell (1, 3, or 9).
 * @tparam OutputType  VTK output precision (default: float).
 */
template <typename FloatType, typename Field_T, uint_t F_SIZE_ARG,
          typename OutputType>
class VTKWriter : public vtk::BlockCellDataWriter<OutputType, F_SIZE_ARG> {
public:
  VTKWriter(ConstBlockDataID const &block_id, std::string const &id,
            FloatType unit_conversion)
      : vtk::BlockCellDataWriter<OutputType, F_SIZE_ARG>(id),
        m_block_id(block_id), m_field(nullptr), m_conversion(unit_conversion) {}

protected:
  void configure() override {
    WALBERLA_ASSERT_NOT_NULLPTR(this->block_);
    m_field = this->block_->template getData<Field_T>(m_block_id);
  }

  ConstBlockDataID const m_block_id;
  Field_T const *m_field;
  FloatType const m_conversion;
};

template <typename FloatType, typename PdfField, typename OutputType = float>
class DensityVTKWriter : public VTKWriter<FloatType, PdfField, 1u, OutputType> {
public:
  using Base = VTKWriter<FloatType, PdfField, 1u, OutputType>;
  using Base::Base;
  using Base::evaluate;

protected:
  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const) override {
    WALBERLA_ASSERT_NOT_NULLPTR(this->m_field);
    auto const density =
        lbm::accessor::Density::get(this->m_field, 1., {x, y, z});
    return numeric_cast<OutputType>(this->m_conversion * density);
  }
};

template <typename FloatType, typename ScalarField, typename OutputType = float>
class PhasefieldVTKWriter
    : public VTKWriter<FloatType, ScalarField, 1u, OutputType> {
public:
  using Base = VTKWriter<FloatType, ScalarField, 1u, OutputType>;
  using Base::Base;
  using Base::evaluate;

protected:
  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const) override {
    WALBERLA_ASSERT_NOT_NULLPTR(this->m_field);
    auto const phasefield = this->m_field->get(x, y, z);
    return numeric_cast<OutputType>(phasefield);
  }
};

template <typename FloatType, typename VectorField, typename OutputType = float>
class VelocityVTKWriter
    : public VTKWriter<FloatType, VectorField, 3u, OutputType> {
public:
  using Base = VTKWriter<FloatType, VectorField, 3u, OutputType>;
  using Base::Base;
  using Base::evaluate;

protected:
  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const f) override {
    WALBERLA_ASSERT_NOT_NULLPTR(this->m_field);
    auto const velocity = lbm::accessor::Vector::get(this->m_field, {x, y, z});
    return numeric_cast<OutputType>(this->m_conversion * velocity[uint_c(f)]);
  }
};

template <typename FloatType, typename PdfField, typename OutputType = float>
class PressureTensorVTKWriter
    : public VTKWriter<FloatType, PdfField, 9u, OutputType> {
public:
  using Base = VTKWriter<FloatType, PdfField, 9u, OutputType>;
  using Base::Base;
  using Base::evaluate;

  PressureTensorVTKWriter(ConstBlockDataID const &block_id,
                          std::string const &id, FloatType unit_conversion,
                          FloatType off_diag_factor)
      : Base(block_id, id, unit_conversion),
        m_off_diag_factor(off_diag_factor) {}

protected:
  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const f) override {
    WALBERLA_ASSERT_NOT_NULLPTR(this->m_field);
    auto const pressure =
        lbm::accessor::PressureTensor::get(this->m_field, 1., {x, y, z});
    auto const revert_factor =
        (f == 0 or f == 4 or f == 8) ? FloatType{1} : m_off_diag_factor;
    return numeric_cast<OutputType>(this->m_conversion * revert_factor *
                                    pressure[uint_c(f)]);
  }
  FloatType const m_off_diag_factor;
};

} // namespace walberla
