/*
 * Copyright (C) 2010-2021 The ESPResSo project
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

#include "ibm_volcons.hpp"
#include "errorhandling.hpp"

/** Set parameters of volume conservation */
IBMVolCons::IBMVolCons(const int softID, const double kappaV) {
  // Specific stuff
  if (softID > IBM_MAX_NUM) {
    runtimeErrorMsg() << "Error: softID (" << softID
                      << ") is larger than IBM_MAX_NUM (" << IBM_MAX_NUM << ")";
  }
  if (softID < 0) {
    runtimeErrorMsg() << "Error: softID (" << softID
                      << ") must be non-negative";
  }

  this->softID = softID;
  this->kappaV = kappaV;
  volRef = 0;
  // NOTE: We cannot compute the reference volume here because not all
  // interactions are setup and thus we do not know which triangles belong to
  // this softID. Calculate it later in the init function of \ref
  // ImmersedBoundaries
}
