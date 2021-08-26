/*
 * Copyright (C) 2015-2019 The ESPResSo project
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

#include "initialize.hpp"

#include "LBBoundaries.hpp"
#include "LBBoundary.hpp"

#include "EKBoundaries.hpp"
#include "EKBoundary.hpp"

namespace ScriptInterface {
namespace LBBoundaries {

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<LBBoundaries>("LBBoundaries::LBBoundaries");
  om->register_new<LBBoundary>("LBBoundaries::LBBoundary");
}
} /* namespace LBBoundaries */
} /* namespace ScriptInterface */

namespace ScriptInterface {
namespace EKBoundaries {

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<EKBoundaries>("EKBoundaries::EKBoundaries");
  om->register_new<EKBoundary>("EKBoundaries::EKBoundary");
}

} // namespace EKBoundaries
} // namespace ScriptInterface
