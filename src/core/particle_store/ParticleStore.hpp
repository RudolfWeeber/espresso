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

#include <cstddef>

/**
 * @brief Array-based particle storage (structure of arrays).
 *
 * Owns every per-particle quantity in a single index space: local
 * particles first, contiguously sorted by cell, ghost particles appended
 * after the locals. Fields are grouped into parameters, state, and
 * observables. See
 * docs/superpowers/specs/2026-07-03-array-based-particle-storage-design.md
 *
 * This is scaffolding (migration phase 0): the columns are added in
 * migration phase 2.
 */
class ParticleStore {
public:
  std::size_t number_of_local_particles() const {
    return m_number_of_local_particles;
  }
  std::size_t number_of_ghost_particles() const {
    return m_number_of_ghost_particles;
  }

private:
  std::size_t m_number_of_local_particles = 0u;
  std::size_t m_number_of_ghost_particles = 0u;
};
