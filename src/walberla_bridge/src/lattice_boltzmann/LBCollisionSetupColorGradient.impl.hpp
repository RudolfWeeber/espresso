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
 * Out-of-class collision model setup definitions for
 * @ref walberla::LBWalberlaImplColorGradient.
 * Only the relaxation-rate helpers are here; the collision model
 * setup methods themselves are in the header (they throw inline).
 */

#include <memory>
#include <stdexcept>
#include <utility>

namespace walberla {

template <typename FloatType, lbmpy::Arch Architecture>
FloatType LBWalberlaImplColorGradient<FloatType, Architecture>::
    shear_mode_relaxation_rate(std::size_t component) const {
  return FloatType{2} / (FloatType{6} * m_viscosity[component] + FloatType{1});
}

template <typename FloatType, lbmpy::Arch Architecture>
FloatType
LBWalberlaImplColorGradient<FloatType, Architecture>::odd_mode_relaxation_rate(
    FloatType shear_relaxation, FloatType magic_number) const {
  return (FloatType{4} - FloatType{2} * shear_relaxation) /
         (FloatType{4} * magic_number * shear_relaxation + FloatType{2} -
          shear_relaxation);
}

} // namespace walberla
