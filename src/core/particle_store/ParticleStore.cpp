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

#include "particle_store/ParticleStore.hpp"

#include "Particle.hpp"

#include <cstddef>

void ParticleStore::begin_rebuild(std::size_t const number_of_local_particles,
                                  std::size_t const number_of_ghost_particles) {
  m_old_force = m_force;
#ifdef ESPRESSO_ROTATION
  m_old_torque = m_torque;
#endif
  m_old_number_of_particles =
      m_number_of_local_particles + m_number_of_ghost_particles;
  m_number_of_local_particles = number_of_local_particles;
  m_number_of_ghost_particles = number_of_ghost_particles;
  auto const total = number_of_local_particles + number_of_ghost_particles;
  // fresh zero-initialized allocation (Kokkos default-initializes to zero)
  m_force = Column("particle_store::force", total);
#ifdef ESPRESSO_ROTATION
  m_torque = Column("particle_store::torque", total);
#endif
}

void ParticleStore::assign_row(Particle &particle, int const row) {
  auto const old_row = particle.store_row();
  if (particle.store() == this and old_row >= 0 and
      static_cast<std::size_t>(old_row) < m_old_number_of_particles) {
    auto const &old_force = m_old_force.view_host();
    auto &new_force = m_force.view_host();
    for (std::size_t j = 0u; j < 3u; ++j) {
      new_force(row, j) = old_force(old_row, j);
    }
#ifdef ESPRESSO_ROTATION
    auto const &old_torque = m_old_torque.view_host();
    auto &new_torque = m_torque.view_host();
    for (std::size_t j = 0u; j < 3u; ++j) {
      new_torque(row, j) = old_torque(old_row, j);
    }
#endif
  }
  particle.attach_to_store(*this, row);
}

void ParticleStore::finish_rebuild() {
  m_old_force = Column{};
#ifdef ESPRESSO_ROTATION
  m_old_torque = Column{};
#endif
  m_old_number_of_particles = 0u;
  m_dirty = false;
}
