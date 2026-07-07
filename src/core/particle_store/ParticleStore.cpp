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

#include <cassert>
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

// Host-sidecar (std::vector of POD) counterpart of preserve_or_seed_scalar:
// copy the whole POD from the old vector at @p old_row (row survived) or from
// the migration carrier @p seed (detached / brand-new particle).
template <class SidecarVector, class SeedType>
void preserve_or_seed_sidecar(SidecarVector &new_sidecar,
                              SidecarVector const &old_sidecar, int const row,
                              int const old_row, bool const preserve,
                              SeedType const &seed) {
  if (preserve) {
    new_sidecar[static_cast<std::size_t>(row)] =
        old_sidecar[static_cast<std::size_t>(old_row)];
  } else {
    new_sidecar[static_cast<std::size_t>(row)] = seed;
  }
}

// Ragged host-sidecar (std::vector of heap-owning element, e.g. BondList /
// compact_vector) counterpart of preserve_or_seed_sidecar. Unlike the POD
// sidecars, the element owns heap storage, so on preserve the surviving element
// is MOVED out of the old vector (transfers the buffer -- no reallocation /
// deep copy of a ragged run) instead of copied. The old vector is discarded
// after the rebuild (retired as the spare), so leaving a moved-from element
// behind is harmless. On seed the carrier value is copied (it is still owned by
// the migrating particle). @p old_sidecar is a non-const reference here because
// the move reads (and empties) its element.
template <class SidecarVector, class SeedType>
void preserve_or_move_sidecar(SidecarVector &new_sidecar,
                              SidecarVector &old_sidecar, int const row,
                              int const old_row, bool const preserve,
                              SeedType const &seed) {
  if (preserve) {
    new_sidecar[static_cast<std::size_t>(row)] =
        std::move(old_sidecar[static_cast<std::size_t>(old_row)]);
  } else {
    new_sidecar[static_cast<std::size_t>(row)] = seed;
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
#ifdef ESPRESSO_BOND_CONSTRAINT
  swap(m_rattle_correction, m_old_rattle_correction);
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
  // Parameter columns (phase 5).
  swap(m_id, m_old_id);
  swap(m_mol_id, m_old_mol_id);
  swap(m_type, m_old_type);
  swap(m_propagation, m_old_propagation);
#ifdef ESPRESSO_ROTATION
  swap(m_rotation, m_old_rotation);
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  swap(m_ext_flag, m_old_ext_flag);
#endif
#ifdef ESPRESSO_MASS
  swap(m_mass, m_old_mass);
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  swap(m_q, m_old_q);
#endif
#ifdef ESPRESSO_DIPOLES
  swap(m_dipm, m_old_dipm);
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  swap(m_rinertia, m_old_rinertia);
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  swap(m_mu_E, m_old_mu_E);
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  swap(m_dip_fld, m_old_dip_fld);
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  swap(m_ext_force, m_old_ext_force);
#ifdef ESPRESSO_ROTATION
  swap(m_ext_torque, m_old_ext_torque);
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  swap(m_gamma, m_old_gamma);
#ifdef ESPRESSO_ROTATION
  swap(m_gamma_rot, m_old_gamma_rot);
#endif
#endif
  // Host sidecars: swap current <-> spare (spare now holds the old-row values,
  // the preserve source; the swapped-in current vector is resized below).
#ifdef ESPRESSO_ENGINE
  swap(m_swimming, m_old_swimming);
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  swap(m_magnetodynamics, m_old_magnetodynamics);
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  swap(m_vs_relative, m_old_vs_relative);
#endif
  // Ragged host sidecars: swap current <-> spare (spare now holds the old-row
  // elements, the preserve/move source; the swapped-in current vector is
  // resized below).
  swap(m_bonds_sidecar, m_old_bonds_sidecar);
#ifdef ESPRESSO_EXCLUSIONS
  swap(m_exclusions_sidecar, m_old_exclusions_sidecar);
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
#ifdef ESPRESSO_BOND_CONSTRAINT
  grow_without_init(m_rattle_correction, total,
                    "particle_store::rattle_correction");
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
  // Parameter columns (phase 5).
  grow_without_init(m_id, total, "particle_store::id");
  grow_without_init(m_mol_id, total, "particle_store::mol_id");
  grow_without_init(m_type, total, "particle_store::type");
  grow_without_init(m_propagation, total, "particle_store::propagation");
#ifdef ESPRESSO_ROTATION
  grow_without_init(m_rotation, total, "particle_store::rotation");
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  grow_without_init(m_ext_flag, total, "particle_store::ext_flag");
#endif
#ifdef ESPRESSO_MASS
  grow_without_init(m_mass, total, "particle_store::mass");
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  grow_without_init(m_q, total, "particle_store::q");
#endif
#ifdef ESPRESSO_DIPOLES
  grow_without_init(m_dipm, total, "particle_store::dipm");
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  grow_without_init(m_rinertia, total, "particle_store::rinertia");
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  grow_without_init(m_mu_E, total, "particle_store::mu_E");
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  grow_without_init(m_dip_fld, total, "particle_store::dip_fld");
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  grow_without_init(m_ext_force, total, "particle_store::ext_force");
#ifdef ESPRESSO_ROTATION
  grow_without_init(m_ext_torque, total, "particle_store::ext_torque");
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  grow_without_init(m_gamma, total, "particle_store::gamma");
#ifdef ESPRESSO_ROTATION
  grow_without_init(m_gamma_rot, total, "particle_store::gamma_rot");
#endif
#endif
  // Host sidecars: resize the (swapped-in) current vector to the new count.
  // Every row is overwritten by assign_row's preserve_or_seed_sidecar before
  // finish_rebuild, so the value-initialized filler resize inserts is never
  // observed; it only ensures the storage exists.
#ifdef ESPRESSO_ENGINE
  m_swimming.resize(total);
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  m_magnetodynamics.resize(total);
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  m_vs_relative.resize(total);
#endif
  // Ragged host sidecars: resize the (swapped-in) current vector to the new
  // count. Every row is overwritten by assign_row's preserve_or_move_sidecar
  // before finish_rebuild; the default-constructed (empty) filler that resize
  // inserts is only observed for a row that assign_row never reaches, which
  // never happens (all [0,total) rows are assigned). A ghost/new row is seeded
  // to an empty element.
  m_bonds_sidecar.resize(total);
#ifdef ESPRESSO_EXCLUSIONS
  m_exclusions_sidecar.resize(total);
#endif
}

// The per-field coverage below is the CANONICAL field list of the store. Two
// consumers must stay in sync with it: (1) ParticleStore::copy_row (row->row
// full copy, same file), and (2) the MigrationPack per-field wire pack
// (particle_store/MigrationPack.cpp). Any field added here must be added there
// (and vice versa); the maximal-population round-trip unit test enforces this.
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

  // RATTLE correction observable column (phase 6). Preserve-by-old-row for
  // uniformity with the other observables (SHAKE zeroes it each iteration, so a
  // preserved value is never actually relied upon); a genuinely-new row
  // defaults to zero. There is NO migration carrier -- the correction is a
  // per-iteration scratch that is never persisted nor migrated -- so the seed
  // is a literal zero vector, not a carrier read.
#ifdef ESPRESSO_BOND_CONSTRAINT
  preserve_or_seed<3u>(m_rattle_correction, m_old_rattle_correction, row,
                       old_row, preserve, Utils::Vector3d{0., 0., 0.});
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

  // Parameter columns (phase 5). Genuinely-new rows are seeded from the
  // migration carriers, whose defaults match the old ParticleProperties member
  // defaults (id -1, mol_id 0, type 0, propagation SYSTEM_DEFAULT, bitfields 0,
  // mass 1, rinertia {1,1,1}, q/dipm 0, mu_E/dip_fld/ext_force/ext_torque zero,
  // gamma/gamma_rot -1 or {-1,-1,-1}). Post-flip the migration carriers are the
  // authoritative parameter store for a detached particle and are serialized
  // with it (like the state/momentum carriers), so a migrated particle carries
  // its parameters through the carriers and this seed reproduces them.
  preserve_or_seed_scalar(m_id, m_old_id, row, old_row, preserve,
                          particle.migration_id());
  preserve_or_seed_scalar(m_mol_id, m_old_mol_id, row, old_row, preserve,
                          particle.migration_mol_id());
  preserve_or_seed_scalar(m_type, m_old_type, row, old_row, preserve,
                          particle.migration_type());
  preserve_or_seed_scalar(m_propagation, m_old_propagation, row, old_row,
                          preserve, particle.migration_propagation());
#ifdef ESPRESSO_ROTATION
  preserve_or_seed_scalar(m_rotation, m_old_rotation, row, old_row, preserve,
                          particle.migration_rotation());
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  preserve_or_seed_scalar(m_ext_flag, m_old_ext_flag, row, old_row, preserve,
                          particle.migration_ext_flag());
#endif
#ifdef ESPRESSO_MASS
  preserve_or_seed_scalar(m_mass, m_old_mass, row, old_row, preserve,
                          particle.migration_mass());
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  preserve_or_seed_scalar(m_q, m_old_q, row, old_row, preserve,
                          particle.migration_q());
#endif
#ifdef ESPRESSO_DIPOLES
  preserve_or_seed_scalar(m_dipm, m_old_dipm, row, old_row, preserve,
                          particle.migration_dipm());
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  preserve_or_seed<3u>(m_rinertia, m_old_rinertia, row, old_row, preserve,
                       particle.migration_rinertia());
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  preserve_or_seed<3u>(m_mu_E, m_old_mu_E, row, old_row, preserve,
                       particle.migration_mu_E());
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  preserve_or_seed<3u>(m_dip_fld, m_old_dip_fld, row, old_row, preserve,
                       particle.migration_dip_fld());
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  preserve_or_seed<3u>(m_ext_force, m_old_ext_force, row, old_row, preserve,
                       particle.migration_ext_force());
#ifdef ESPRESSO_ROTATION
  preserve_or_seed<3u>(m_ext_torque, m_old_ext_torque, row, old_row, preserve,
                       particle.migration_ext_torque());
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  preserve_or_seed<3u>(m_gamma, m_old_gamma, row, old_row, preserve,
                       particle.migration_gamma());
#ifdef ESPRESSO_ROTATION
  preserve_or_seed<3u>(m_gamma_rot, m_old_gamma_rot, row, old_row, preserve,
                       particle.migration_gamma_rot());
#endif
#else // ESPRESSO_PARTICLE_ANISOTROPY
  preserve_or_seed_scalar(m_gamma, m_old_gamma, row, old_row, preserve,
                          particle.migration_gamma());
#ifdef ESPRESSO_ROTATION
  preserve_or_seed_scalar(m_gamma_rot, m_old_gamma_rot, row, old_row, preserve,
                          particle.migration_gamma_rot());
#endif
#endif // ESPRESSO_PARTICLE_ANISOTROPY
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE

  // Host sidecars (phase 5): whole-POD preserve-by-old-row / seed-from-carrier.
#ifdef ESPRESSO_ENGINE
  preserve_or_seed_sidecar(m_swimming, m_old_swimming, row, old_row, preserve,
                           particle.migration_swimming());
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  preserve_or_seed_sidecar(m_magnetodynamics, m_old_magnetodynamics, row,
                           old_row, preserve,
                           particle.migration_magnetodynamics());
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  preserve_or_seed_sidecar(m_vs_relative, m_old_vs_relative, row, old_row,
                           preserve, particle.migration_vs_relative());
#endif

  // Ragged host sidecars (phase 6): a surviving row is MOVED out of the old
  // vector element (transfers the heap buffer -- no deep copy of the ragged
  // run); a genuinely-new / migrated row is seeded (copied) from the migration
  // carrier, which for the dormant phase reads the Particle's bl/el members
  // directly (empty for a fresh ghost). Nothing reads these sidecars yet.
  preserve_or_move_sidecar(m_bonds_sidecar, m_old_bonds_sidecar, row, old_row,
                           preserve, particle.migration_bonds());
#ifdef ESPRESSO_EXCLUSIONS
  preserve_or_move_sidecar(m_exclusions_sidecar, m_old_exclusions_sidecar, row,
                           old_row, preserve, particle.migration_exclusions());
#endif

  particle.attach_to_store(*this, row);
}

Particle ParticleStore::make_view(int const row) {
  assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
  Particle view;
  view.attach_to_store(*this, row);
  return view;
}

Particle ParticleStore::snapshot_row(int const row) {
  assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
  Particle snapshot;
  // Attach so the detached_*() getters read THIS store's columns/sidecars at
  // `row`, sync those live values into the carriers, then detach so the result
  // owns its data (all accessors read the carriers). Also carries the ghost
  // flag (stored in `l`, not a column) forward from the view.
  snapshot.attach_to_store(*this, row);
  snapshot.detach_from_store();
  return snapshot;
}

// Row-to-row full copy. Field coverage IDENTICAL to assign_row (see the sync
// note above assign_row). Reads each column/sidecar of `source` at `src` by
// value and writes it into THIS store at `dst` through the element reference.
// The store owns no ghost flag (that lives in Particle::l, not a column), so
// there is nothing to copy for it here -- consistent with assign_row, which
// also never touches a ghost flag.
void ParticleStore::copy_row(ParticleStore const &source, int const src,
                             int const dst) {
  assert(src >= 0 and
         static_cast<std::size_t>(src) < source.number_of_particles());
  assert(dst >= 0 and static_cast<std::size_t>(dst) < number_of_particles());

  // Observable columns (phase 2).
  force_reference(dst) = source.force_value(src);
#ifdef ESPRESSO_ROTATION
  torque_reference(dst) = source.torque_value(src);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  rattle_correction_reference(dst) = source.rattle_correction_value(src);
#endif

  // State columns (phase 3).
  position_reference(dst) = source.position_value(src);
  image_box_reference(dst) = source.image_box_value(src);
#ifdef ESPRESSO_ROTATION
  quaternion_reference(dst) = source.quaternion_value(src);
#endif
  position_at_last_verlet_update_reference(dst) =
      source.position_at_last_verlet_update_value(src);
#ifdef ESPRESSO_BOND_CONSTRAINT
  position_last_time_step_reference(dst) =
      source.position_last_time_step_value(src);
#endif
  lees_edwards_offset(dst) = source.lees_edwards_offset(src);
  lees_edwards_flag(dst) = source.lees_edwards_flag(src);

  // Momentum columns (phase 4).
  velocity_reference(dst) = source.velocity_value(src);
#ifdef ESPRESSO_ROTATION
  angular_velocity_reference(dst) = source.angular_velocity_value(src);
#endif

  // Parameter columns (phase 5).
  id(dst) = source.id(src);
  mol_id(dst) = source.mol_id(src);
  type(dst) = source.type(src);
  propagation(dst) = source.propagation(src);
#ifdef ESPRESSO_ROTATION
  rotation(dst) = source.rotation(src);
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  ext_flag(dst) = source.ext_flag(src);
#endif
#ifdef ESPRESSO_MASS
  mass(dst) = source.mass(src);
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  q(dst) = source.q(src);
#endif
#ifdef ESPRESSO_DIPOLES
  dipm(dst) = source.dipm(src);
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  rinertia_reference(dst) = source.rinertia_value(src);
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  mu_E_reference(dst) = source.mu_E_value(src);
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  dip_fld_reference(dst) = source.dip_fld_value(src);
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  ext_force_reference(dst) = source.ext_force_value(src);
#ifdef ESPRESSO_ROTATION
  ext_torque_reference(dst) = source.ext_torque_value(src);
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  gamma_reference(dst) = source.gamma_value(src);
#ifdef ESPRESSO_ROTATION
  gamma_rot_reference(dst) = source.gamma_rot_value(src);
#endif
#endif

  // Host POD sidecars (phase 5).
#ifdef ESPRESSO_ENGINE
  swimming(dst) = source.swimming(src);
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  magnetodynamics(dst) = source.magnetodynamics(src);
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  vs_relative(dst) = source.vs_relative(src);
#endif

  // Ragged host sidecars (phase 6): copied by value (deep copy of the run).
  bonds_sidecar_reference(dst) = source.bonds_sidecar_reference(src);
#ifdef ESPRESSO_EXCLUSIONS
  exclusions_sidecar_reference(dst) = source.exclusions_sidecar_reference(src);
#endif
}

void ParticleStore::finish_rebuild() {
  ++m_generation;
  // Keep the old-generation columns alive as the spare buffer for the next
  // rebuild's swap (capacity-cached double buffering). Only the bookkeeping is
  // cleared; release_columns drops BOTH generations at teardown.
  m_old_number_of_particles = 0u;
  m_dirty = false;
}
