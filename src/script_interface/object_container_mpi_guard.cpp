/*
 * Copyright (C) 2021 The ESPResSo project
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

#include "script_interface/object_container_mpi_guard.hpp"

#include "core/communication.hpp"

#include <boost/utility/string_ref.hpp>

#include <cstddef>
#include <sstream>
#include <stdexcept>

void object_container_mpi_guard(boost::string_ref const &name,
                                std::size_t n_elements) {
  if (comm_cart.size() > 1 and n_elements) {
    std::stringstream error_msg;
    error_msg << "Non-empty object containers do not support checkpointing in "
              << "MPI environments. Container " << name << " contains "
              << n_elements << " elements.";
    throw std::runtime_error(error_msg.str());
  }
}
