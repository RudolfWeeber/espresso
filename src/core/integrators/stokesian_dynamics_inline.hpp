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

#pragma once

#include <config/config.hpp>

#ifdef ESPRESSO_STOKESIAN_DYNAMICS

#include "rotation.hpp"
#include "stokesian_dynamics/sd_interface.hpp"
#include "thermostat.hpp"

inline void stokesian_dynamics_step_1(ParticleRangeStokesian const &particles,
                                      StokesianDynamics const &integrator,
                                      StokesianThermostat const &stokesian,
                                      double time_step, double kT) {
  integrator.propagate_vel_pos(particles, stokesian, time_step, kT);

  for (auto &p : particles) {
    // translate
    p.pos() += p.v() * time_step;
    // rotate around the space-frame angular velocity, without the
    // components along blocked body-frame axes
    auto const dphi = Utils::mask(
        p.rotation(), convert_vector_space_to_body(p, time_step * p.omega()));
    auto const angle = dphi.norm();
    if (angle != 0.) {
      p.quat() = local_rotate_particle_body(p, dphi / angle, angle);
    }
  }
}

#endif // ESPRESSO_STOKESIAN_DYNAMICS
