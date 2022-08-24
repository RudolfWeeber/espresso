/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "rotate_system.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "config/config.hpp"
#include "event.hpp"
#include "rotation.hpp"

#include <utils/Vector.hpp>
#include <utils/math/vec_rotate.hpp>

#include <boost/mpi/collectives.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>

void rotate_system(double phi, double theta, double alpha) {
  auto const particles = cell_structure.local_particles();

  // Calculate center of mass
  Utils::Vector3d local_com{};
  double local_mass = 0.0;

  for (auto const &p : particles) {
    if (not p.is_virtual()) {
      local_com += p.mass() * p.pos();
      local_mass += p.mass();
    }
  }

  auto const total_mass =
      boost::mpi::all_reduce(comm_cart, local_mass, std::plus<>());
  auto const com =
      boost::mpi::all_reduce(comm_cart, local_com, std::plus<>()) / total_mass;

  // Rotation axis in Cartesian coordinates
  Utils::Vector3d axis;
  axis[0] = std::sin(theta) * std::cos(phi);
  axis[1] = std::sin(theta) * std::sin(phi);
  axis[2] = std::cos(theta);

  // Rotate particle coordinates
  for (auto &p : particles) {
    // Move the center of mass of the system to the origin
    p.pos() = com + Utils::vec_rotate(axis, alpha, p.pos() - com);
#ifdef ROTATION
    local_rotate_particle(p, axis, alpha);
#endif
  }

  cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
  on_particle_change();
  update_dependent_particles();
}
