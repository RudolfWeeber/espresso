/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "config/config.hpp"

#include "TabulatedPotential.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/optional.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cassert>
#include <cmath>
#include <memory>
#include <tuple>
#include <vector>

struct VelDependentTabulated {
  static const int num = 1; // bond partner
  TabulatedPotential approach_tab;
  TabulatedPotential reced_tab;
  VelDependentTabulated() = default;
  VelDependentTabulated(double min_val, double max_val,
                        const std::vector<double> &approach_force,
                        const std::vector<double> &reced_force) {
    // energy is not defined for a history depenent force
    if (approach_force.size() != reced_force.size()) {
      throw std::domain_error(
          "approach_force and reced_force must have the same length");
    }
    std::vector<double> energy(approach_force.size());
    approach_tab = {min_val, max_val, approach_force, energy};
    reced_tab = {min_val, max_val, reced_force, energy};
  }
  double cutoff() const { return approach_tab.cutoff(); };
  boost::optional<Utils::Vector3d> force(const Utils::Vector3d &dx,
                                         const Utils::Vector3d &dv) const;
  boost::optional<double> energy(const Utils::Vector3d &dx,
                                 const Utils::Vector3d &dv) {
    return 0.0;
  };

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &approach_tab;
    ar &reced_tab;
  }
};
inline boost::optional<Utils::Vector3d>
VelDependentTabulated::force(Utils::Vector3d const &dx,
                             const Utils::Vector3d &dv) const {
  auto const dist = dx.norm();
  if (dist >= cutoff())
    return {};

  double rel_motion = dx * dv;

  if (rel_motion >= 0) {
    return reced_tab.force(dist) / dist * dx;
  } else {
    return approach_tab.force(dist) / dist * dx;
  }
}
