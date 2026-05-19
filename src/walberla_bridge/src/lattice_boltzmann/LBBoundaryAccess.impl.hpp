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
 * Out-of-class boundary access definitions for
 * @ref walberla::LBWalberlaImpl. Most boundary methods have moved
 * into @ref walberla::LBWalberlaCommon; this file retains the leaf-
 * specific helper(s) that depend on the PDF field id.
 */

#include <stdexcept>

namespace walberla {

/**
 * @brief Lazily enable boundary mode on first boundary addition.
 * Switches the streaming communicator to the generic pack info,
 * which correctly handles boundary-adjacent cells.
 */
template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::on_boundary_add() {
  if (has_two_components()) {
    throw std::runtime_error(
        "boundaries are not implemented for two-component LB");
  }
  if (not m_has_boundaries) {
    m_has_boundaries = true;
    setup_streaming_communicator();
  }
  m_has_boundaries = true;
}

} // namespace walberla
