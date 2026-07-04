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

#include "particle_store/VectorReference.hpp"

#include <utils/Vector.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include <cassert>
#include <cstddef>

class Particle; // attach_to_store is defined in Particle.hpp

/**
 * @brief Array-based particle storage (structure of arrays).
 *
 * Owns per-particle quantities in a single index space: local particles
 * first (cell-traversal order), ghosts appended. Fields are component-major
 * (LayoutLeft) Kokkos columns; see
 * docs/superpowers/specs/2026-07-03-array-based-particle-storage-design.md
 *
 * Migration phase 2: force and torque columns (observables). Rebuild
 * protocol: mark_dirty() on any topology change; the owner (CellStructure)
 * later runs begin_rebuild / assign_row-per-particle / finish_rebuild.
 * Rebuild preserves values by old row and seeds new rows from the particle's
 * migration carrier (zero for genuinely new particles).
 * Rebuilds are purely rank-local (no MPI).
 */
class ParticleStore {
public:
  using Column = Kokkos::DualView<double *[3], Kokkos::LayoutLeft>;

  std::size_t number_of_local_particles() const {
    return m_number_of_local_particles;
  }
  std::size_t number_of_ghost_particles() const {
    return m_number_of_ghost_particles;
  }
  std::size_t number_of_particles() const {
    return m_number_of_local_particles + m_number_of_ghost_particles;
  }

  void mark_dirty() { m_dirty = true; }
  bool is_dirty() const { return m_dirty; }

  void begin_rebuild(std::size_t number_of_local_particles,
                     std::size_t number_of_ghost_particles);
  void assign_row(Particle &particle, int row);
  void finish_rebuild();

  /** Release all Kokkos-backed columns. Must be called while the Kokkos
   *  runtime is still alive (e.g. before Kokkos::finalize); afterwards the
   *  store is empty and dirty. */
  void release_columns() {
    m_force = Column{};
#ifdef ESPRESSO_ROTATION
    m_torque = Column{};
#endif
    m_old_force = Column{};
#ifdef ESPRESSO_ROTATION
    m_old_torque = Column{};
#endif
    m_number_of_local_particles = 0u;
    m_number_of_ghost_particles = 0u;
    m_old_number_of_particles = 0u;
    m_dirty = true;
  }

  VectorReference force_reference(int const row) {
    return column_reference(m_force, row);
  }
  Utils::Vector3d force_value(int const row) const {
    return column_value(m_force, row);
  }
#ifdef ESPRESSO_ROTATION
  VectorReference torque_reference(int const row) {
    return column_reference(m_torque, row);
  }
  Utils::Vector3d torque_value(int const row) const {
    return column_value(m_torque, row);
  }
#endif

private:
  VectorReference column_reference(Column &column, int const row) {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    auto &view = column.view_host();
    // stride between the three components of one row (LayoutLeft: number of
    // rows). stride(1) is the non-deprecated spelling of stride_1() in the
    // installed Kokkos version.
    return VectorReference(view.data() + row, view.stride(1));
  }
  Utils::Vector3d column_value(Column const &column, int const row) const {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    auto const &view = column.view_host();
    return {view(row, 0), view(row, 1), view(row, 2)};
  }

  std::size_t m_number_of_local_particles = 0u;
  std::size_t m_number_of_ghost_particles = 0u;
  bool m_dirty = false;
  Column m_force;
#ifdef ESPRESSO_ROTATION
  Column m_torque;
#endif
  // previous-generation columns, alive between begin_rebuild/finish_rebuild
  Column m_old_force;
#ifdef ESPRESSO_ROTATION
  Column m_old_torque;
#endif
  std::size_t m_old_number_of_particles = 0u;
};
