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
#include <string>
#include <utility>

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

namespace {
// Grow @p column to hold at least @p total rows, WITHOUT zero-initializing the
// new storage. Only reallocates when the current extent is too small; existing
// capacity is reused verbatim. Safe because assign_row's preserve_or_seed
// writes every logical row in [0, total) before finish_rebuild, so no reader
// ever observes uninitialized storage (see begin_rebuild's contract below).
template <class ColumnType>
void grow_without_init(ColumnType &column, std::size_t const total,
                       char const *label) {
  if (column.extent(0) < total) {
    // Kokkos::view_alloc treats a bare char const* as a pointer-to-memory
    // (unmanaged view); wrap it in a std::string so it is taken as the label.
    column = ColumnType(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string{label}),
        total);
  }
}
} // namespace

// Capacity-cached double buffering (phase 3.5 perf reclamation). Instead of
// freshly allocating all columns on every rebuild (every resort), keep two
// generations and SWAP them: the previous-generation columns become the write
// target for the new generation, and the just-current columns become the
// read source ("old") for preserve_or_seed. A (re)allocation happens ONLY when
// the swapped-in write target is too small for the new particle count; growth
// uses WithoutInitializing since every logical row is overwritten by assign_row
// before it is read. Column extents therefore track CAPACITY (a high-water
// mark), never the logical count -- accessors bound by number_of_particles().
void ParticleStore::begin_rebuild(std::size_t const number_of_local_particles,
                                  std::size_t const number_of_ghost_particles) {
  // Swap current <-> spare generations. After the swap, m_old_* holds the
  // just-current data (the read source for surviving rows) and m_* holds the
  // older spare buffer (the write target we grow / overwrite in place).
  using std::swap;
  swap(m_force, m_old_force);
#ifdef ESPRESSO_ROTATION
  swap(m_torque, m_old_torque);
#endif
  swap(m_position, m_old_position);
  swap(m_image_box, m_old_image_box);
#ifdef ESPRESSO_ROTATION
  swap(m_quaternion, m_old_quaternion);
#endif
  swap(m_position_at_last_verlet_update, m_old_position_at_last_verlet_update);
#ifdef ESPRESSO_BOND_CONSTRAINT
  swap(m_position_last_time_step, m_old_position_last_time_step);
#endif
  swap(m_lees_edwards_offset, m_old_lees_edwards_offset);
  swap(m_lees_edwards_flag, m_old_lees_edwards_flag);
  swap(m_velocity, m_old_velocity);
#ifdef ESPRESSO_ROTATION
  swap(m_angular_velocity, m_old_angular_velocity);
#endif

  m_old_number_of_particles =
      m_number_of_local_particles + m_number_of_ghost_particles;
  m_number_of_local_particles = number_of_local_particles;
  m_number_of_ghost_particles = number_of_ghost_particles;
  auto const total = number_of_local_particles + number_of_ghost_particles;

  // Grow the (swapped-in) write target only when it cannot hold the new count.
  // NOTE: rows are NOT zero-initialized; assign_row seeds/preserves every row.
  // A quaternion row that survives is preserved; a genuinely-new row is seeded
  // to identity (1,0,0,0) from the migration carrier (never the zero-init).
  grow_without_init(m_force, total, "particle_store::force");
#ifdef ESPRESSO_ROTATION
  grow_without_init(m_torque, total, "particle_store::torque");
#endif
  grow_without_init(m_position, total, "particle_store::position");
  grow_without_init(m_image_box, total, "particle_store::image_box");
#ifdef ESPRESSO_ROTATION
  grow_without_init(m_quaternion, total, "particle_store::quaternion");
#endif
  grow_without_init(m_position_at_last_verlet_update, total,
                    "particle_store::position_at_last_verlet_update");
#ifdef ESPRESSO_BOND_CONSTRAINT
  grow_without_init(m_position_last_time_step, total,
                    "particle_store::position_last_time_step");
#endif
  grow_without_init(m_lees_edwards_offset, total,
                    "particle_store::lees_edwards_offset");
  grow_without_init(m_lees_edwards_flag, total,
                    "particle_store::lees_edwards_flag");
  grow_without_init(m_velocity, total, "particle_store::velocity");
#ifdef ESPRESSO_ROTATION
  grow_without_init(m_angular_velocity, total,
                    "particle_store::angular_velocity");
#endif
}

void ParticleStore::assign_row(Particle &particle, int const row) {
  assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
  auto const old_row = particle.store_row();
  // A row "survives" a rank-local rebuild when the particle was already
  // attached to THIS store and its old row is within the previous generation.
  // Otherwise the particle is detached (brand-new or just migrated from
  // another rank): its columns are seeded from the migration carriers, which
  // hold the values ferried through Particle serialization (see Particle.hpp).
  // Reading a raw carrier never touches the fresh columns (allocated
  // WithoutInitializing; every row is written by assign_row before
  // finish_rebuild).
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

  // Momentum columns (phase 4). Genuinely-new rows default to zero velocity,
  // matching Particle's ParticleMomentum member defaults (m.v = {0,0,0},
  // m.omega = {0,0,0}).
  preserve_or_seed<3u>(m_velocity, m_old_velocity, row, old_row, preserve,
                       particle.migration_velocity());
#ifdef ESPRESSO_ROTATION
  preserve_or_seed<3u>(m_angular_velocity, m_old_angular_velocity, row, old_row,
                       preserve, particle.migration_angular_velocity());
#endif

  particle.attach_to_store(*this, row);
}

void ParticleStore::finish_rebuild() {
  ++m_generation;
  // Keep the old-generation columns alive as the spare buffer for the next
  // rebuild's swap (capacity-cached double buffering). Only the bookkeeping is
  // cleared; release_columns drops BOTH generations at teardown.
  m_old_number_of_particles = 0u;
  m_dirty = false;
}
