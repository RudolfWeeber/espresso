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

#include <utils/Vector.hpp>
#include <utils/quaternion.hpp>

#include <cstddef>

namespace {
// Generalized per-column rebuild copy: fill row @p row of the freshly
// allocated @p new_column either from @p old_column at @p old_row (row
// survived a rank-local rebuild) or from @p seed (the migration carrier /
// default for a detached or brand-new particle). @tparam N is the number of
// component columns (1 for scalars, 3 for vectors, 4 for quaternions); the
// seed is an object indexable with operator[] over [0, N).
template <std::size_t N, class ColumnType, class SeedType>
void preserve_or_seed(ColumnType &new_column, ColumnType const &old_column,
                      int const row, int const old_row, bool const preserve,
                      SeedType const &seed) {
  auto &new_view = new_column.view_host();
  if (preserve) {
    auto const &old_view = old_column.view_host();
    for (std::size_t j = 0u; j < N; ++j) {
      new_view(row, j) = old_view(old_row, j);
    }
  } else {
    for (std::size_t j = 0u; j < N; ++j) {
      new_view(row, j) = seed[j];
    }
  }
}

// Scalar (single-component) counterpart of preserve_or_seed.
template <class ColumnType, class SeedType>
void preserve_or_seed_scalar(ColumnType &new_column,
                             ColumnType const &old_column, int const row,
                             int const old_row, bool const preserve,
                             SeedType const &seed) {
  auto &new_view = new_column.view_host();
  if (preserve) {
    new_view(row) = old_column.view_host()(old_row);
  } else {
    new_view(row) = seed;
  }
}
} // namespace

void ParticleStore::begin_rebuild(std::size_t const number_of_local_particles,
                                  std::size_t const number_of_ghost_particles) {
  m_old_force = m_force;
#ifdef ESPRESSO_ROTATION
  m_old_torque = m_torque;
#endif
  m_old_position = m_position;
  m_old_image_box = m_image_box;
#ifdef ESPRESSO_ROTATION
  m_old_quaternion = m_quaternion;
#endif
  m_old_position_at_last_verlet_update = m_position_at_last_verlet_update;
#ifdef ESPRESSO_BOND_CONSTRAINT
  m_old_position_last_time_step = m_position_last_time_step;
#endif
  m_old_lees_edwards_offset = m_lees_edwards_offset;
  m_old_lees_edwards_flag = m_lees_edwards_flag;
  m_old_number_of_particles =
      m_number_of_local_particles + m_number_of_ghost_particles;
  m_number_of_local_particles = number_of_local_particles;
  m_number_of_ghost_particles = number_of_ghost_particles;
  auto const total = number_of_local_particles + number_of_ghost_particles;
  // fresh zero-initialized allocations (Kokkos default-initializes to zero).
  // NOTE: quaternion rows land as (0,0,0,0), an INVALID quaternion; assign_row
  // seeds each row (identity for genuinely new rows), see below.
  m_force = Column("particle_store::force", total);
#ifdef ESPRESSO_ROTATION
  m_torque = Column("particle_store::torque", total);
#endif
  m_position = Column("particle_store::position", total);
  m_image_box = IntegerColumn("particle_store::image_box", total);
#ifdef ESPRESSO_ROTATION
  m_quaternion = QuaternionColumn("particle_store::quaternion", total);
#endif
  m_position_at_last_verlet_update =
      Column("particle_store::position_at_last_verlet_update", total);
#ifdef ESPRESSO_BOND_CONSTRAINT
  m_position_last_time_step =
      Column("particle_store::position_last_time_step", total);
#endif
  m_lees_edwards_offset =
      ScalarColumn("particle_store::lees_edwards_offset", total);
  m_lees_edwards_flag = ShortColumn("particle_store::lees_edwards_flag", total);
}

void ParticleStore::assign_row(Particle &particle, int const row) {
  assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
  auto const old_row = particle.store_row();
  // A row "survives" a rank-local rebuild when the particle was already
  // attached to THIS store and its old row is within the previous generation.
  // Otherwise the particle is detached (brand-new or just migrated from
  // another rank): its columns are seeded from the migration carriers, which
  // hold the values ferried through Particle serialization (see Particle.hpp).
  // Reading a raw carrier never touches the fresh (zero-initialized) columns.
  bool const preserve =
      particle.store() == this and old_row >= 0 and
      static_cast<std::size_t>(old_row) < m_old_number_of_particles;

  // Observable columns (phase 2).
  preserve_or_seed<3u>(m_force, m_old_force, row, old_row, preserve,
                       particle.migration_force());
#ifdef ESPRESSO_ROTATION
  preserve_or_seed<3u>(m_torque, m_old_torque, row, old_row, preserve,
                       particle.migration_torque());
#endif

  // State columns (phase 3). Genuinely-new rows (detached, carriers at their
  // defaults) must match Particle's member defaults: position/p_old zero,
  // image box zero, quaternion IDENTITY (1,0,0,0), Lees-Edwards offset/flag 0.
  // The migration carriers return exactly those defaults for a freshly
  // constructed Particle, so seeding from them yields the correct new-row
  // values (in particular an identity quaternion, never the zero-init).
  preserve_or_seed<3u>(m_position, m_old_position, row, old_row, preserve,
                       particle.migration_position());
  preserve_or_seed<3u>(m_image_box, m_old_image_box, row, old_row, preserve,
                       particle.migration_image_box());
#ifdef ESPRESSO_ROTATION
  preserve_or_seed<4u>(m_quaternion, m_old_quaternion, row, old_row, preserve,
                       particle.migration_quaternion());
#endif
  preserve_or_seed<3u>(m_position_at_last_verlet_update,
                       m_old_position_at_last_verlet_update, row, old_row,
                       preserve,
                       particle.migration_position_at_last_verlet_update());
#ifdef ESPRESSO_BOND_CONSTRAINT
  preserve_or_seed<3u>(m_position_last_time_step, m_old_position_last_time_step,
                       row, old_row, preserve,
                       particle.migration_position_last_time_step());
#endif
  preserve_or_seed_scalar(m_lees_edwards_offset, m_old_lees_edwards_offset, row,
                          old_row, preserve,
                          particle.migration_lees_edwards_offset());
  preserve_or_seed_scalar(m_lees_edwards_flag, m_old_lees_edwards_flag, row,
                          old_row, preserve,
                          particle.migration_lees_edwards_flag());

  particle.attach_to_store(*this, row);
}

void ParticleStore::finish_rebuild() {
  m_old_force = Column{};
#ifdef ESPRESSO_ROTATION
  m_old_torque = Column{};
#endif
  m_old_position = Column{};
  m_old_image_box = IntegerColumn{};
#ifdef ESPRESSO_ROTATION
  m_old_quaternion = QuaternionColumn{};
#endif
  m_old_position_at_last_verlet_update = Column{};
#ifdef ESPRESSO_BOND_CONSTRAINT
  m_old_position_last_time_step = Column{};
#endif
  m_old_lees_edwards_offset = ScalarColumn{};
  m_old_lees_edwards_flag = ShortColumn{};
  m_old_number_of_particles = 0u;
  m_dirty = false;
}
