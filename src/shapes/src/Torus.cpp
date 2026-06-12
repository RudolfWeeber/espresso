/*
 * Copyright (C) 2010-2026 The ESPResSo project
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
#include <cmath>
#include <shapes/Torus.hpp>
#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

namespace Shapes {

void Torus::calculate_dist(const Utils::Vector3d &pos, double &dist,
                           Utils::Vector3d &vec) const {
  /* Coordinate transform to cylinder coords
     with origin at m_center. */
  auto const c_dist = pos - m_center;
  auto const z = e_z * c_dist;
  auto const r_vec = c_dist - z * e_z;
  auto const r = r_vec.norm();

  dist = (std::sqrt(Utils::sqr(r - m_rad) + z * z) - m_tube_rad) * m_direction;

  /* If exactly on the symmetry axis (r == 0), r_vec / r is undefined (0/0).
     Use the precomputed e_r_axis as a fallback, matching the convention in
     Cylinder::calculate_dist. */
  auto const e_r = (r == 0.) ? e_r_axis : r_vec / r;
  auto const dir_vec = c_dist - e_r * m_rad;
  auto const dir_vec_norm = dir_vec / dir_vec.norm();
  vec = dir_vec_norm * std::abs(dist);
}
} // namespace Shapes
