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

#include "BondList.hpp"
#include "PropagationMode.hpp"
#include "particle_store/ParticleParameters.hpp"
#include "particle_store/ParticleStore.hpp"

#include <utils/Vector.hpp>
#include <utils/compact_vector.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/quaternion.hpp>

#include <boost/container/vector.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/serialization/level.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

namespace detail {
inline bool get_nth_bit(uint8_t const bitfield, unsigned int const bit_idx) {
  return bitfield & (1u << bit_idx);
}
} // namespace detail

// ParticleParametersSwimming, ThermalStonerWohlfarthParameters and
// VirtualSitesRelativeParameters (the cold parameter PODs) are defined in
// particle_store/ParticleParameters.hpp so that both this header and
// ParticleStore can name the complete types (migration phase 5).

/** Properties of a particle which are not supposed to
 *  change during the integration, but have to be known
 *  for all ghosts. Ghosts are particles which are
 *  needed in the interaction calculation, but are just copies of
 *  particles stored on different nodes.
 *
 *  Migration phase 5: these parameter fields have left the @ref Particle
 *  struct; they live only in the @ref ParticleStore columns / host sidecars
 *  (single ownership, spec section 4). @ref ParticleProperties is retained
 *  purely as a type anchor: several call sites name its field types via
 *  @c decltype(ParticleProperties::identity) etc., the
 *  @c ParticleProperties::VirtualSitesRelativeParameters member typedef keeps
 *  the historical script-interface spelling compiling, and the standalone
 *  properties-serialization unit test exercises it. It is no longer a member
 *  of @ref Particle nor part of @ref Particle serialization.
 */
struct ParticleProperties {
  /** unique identifier for the particle. */
  int identity = -1;
  /** Molecule identifier. */
  int mol_id = 0;
  /** particle type, used for non-bonded interactions. */
  int type = 0;
  /** which propagation schemes should be applied to the particle **/
  int propagation = PropagationMode::SYSTEM_DEFAULT;

#ifdef ESPRESSO_ROTATION
  /** Bitfield for the particle axes of rotation.
   *  Values:
   *  - 0: no rotation
   *  - 1: allow rotation around the x axis
   *  - 2: allow rotation around the y axis
   *  - 4: allow rotation around the z axis
   *  By default, the particle cannot rotate.
   */
  uint8_t rotation = static_cast<uint8_t>(0b000u);
#else
  /** Bitfield for the particle axes of rotation. Particle cannot rotate. */
  static constexpr uint8_t rotation = static_cast<uint8_t>(0b000u);
#endif

#ifdef ESPRESSO_EXTERNAL_FORCES
  /** Flag for fixed particle coordinates.
   *  Values:
   *  - 0: no fixed coordinates
   *  - 1: fix translation along the x axis
   *  - 2: fix translation along the y axis
   *  - 4: fix translation along the z axis
   */
  uint8_t ext_flag = static_cast<uint8_t>(0b000u);
#else  // ESPRESSO_EXTERNAL_FORCES
  /** Bitfield for fixed particle coordinates. Coordinates cannot be fixed. */
  static constexpr uint8_t ext_flag = static_cast<uint8_t>(0b000u);
#endif // ESPRESSO_EXTERNAL_FORCES

  /** particle mass */
#ifdef ESPRESSO_MASS
  double mass = 1.0;
#else
  constexpr static double mass{1.0};
#endif

  /** rotational inertia */
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Utils::Vector3d rinertia = {1., 1., 1.};
#else
  static constexpr Utils::Vector3d rinertia = {1., 1., 1.};
#endif

  /** charge. */
#ifdef ESPRESSO_ELECTROSTATICS
  double q = 0.0;
#else
  constexpr static double q{0.0};
#endif

#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  /** electrophoretic mobility times E-field: mu_0 * E */
  Utils::Vector3d mu_E = {0., 0., 0.};
#endif

#ifdef ESPRESSO_DIPOLES
  /** dipole moment (absolute value) */
  double dipm = 0.;
#endif

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  /** total dipole field */
  Utils::Vector3d dip_fld = {0., 0., 0.};
#endif

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  /** @see ::VirtualSitesRelativeParameters (defined in
   *  particle_store/ParticleParameters.hpp). Kept as a member typedef so that
   *  the historical spelling @c
   * ParticleProperties::VirtualSitesRelativeParameters used by the script
   * interface keeps compiling. */
  using VirtualSitesRelativeParameters = ::VirtualSitesRelativeParameters;
  VirtualSitesRelativeParameters vs_relative;
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE

#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
/** Friction coefficient for translation */
#ifndef ESPRESSO_PARTICLE_ANISOTROPY
  double gamma = -1.;
#else
  Utils::Vector3d gamma = {-1., -1., -1.};
#endif // ESPRESSO_PARTICLE_ANISOTROPY
#ifdef ESPRESSO_ROTATION
/** Friction coefficient for rotation */
#ifndef ESPRESSO_PARTICLE_ANISOTROPY
  double gamma_rot = -1.;
#else
  Utils::Vector3d gamma_rot = {-1., -1., -1.};
#endif // ESPRESSO_PARTICLE_ANISOTROPY
#endif // ESPRESSO_ROTATION
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE

#ifdef ESPRESSO_EXTERNAL_FORCES
  /** External force. */
  Utils::Vector3d ext_force = {0., 0., 0.};
#ifdef ESPRESSO_ROTATION
  /** External torque. */
  Utils::Vector3d ext_torque = {0., 0., 0.};
#endif // ESPRESSO_ROTATION
#endif // ESPRESSO_EXTERNAL_FORCES

#ifdef ESPRESSO_ENGINE
  ParticleParametersSwimming swim;
#endif

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  ThermalStonerWohlfarthParameters magnetodynamics;
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & identity;
    ar & mol_id;
    ar & type;
    ar & propagation;
#ifdef ESPRESSO_MASS
    ar & mass;
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    ar & rinertia;
#endif
#ifdef ESPRESSO_ROTATION
    ar & rotation;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    ar & q;
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
    ar & mu_E;
#endif
#ifdef ESPRESSO_DIPOLES
    ar & dipm;
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    ar & dip_fld;
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    ar & vs_relative;
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    ar & gamma;
#ifdef ESPRESSO_ROTATION
    ar & gamma_rot;
#endif
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_EXTERNAL_FORCES
    ar & ext_flag;
    ar & ext_force;
#ifdef ESPRESSO_ROTATION
    ar & ext_torque;
#endif
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_ENGINE
    ar & swim;
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
    ar & magnetodynamics;
#endif
  }
};

/** Positional information on a particle. Information that is
 *  communicated to calculate interactions with ghost particles.
 *
 *  Migration phase 3: the position, image box, quaternion and
 *  position-at-last-time-step fields have left this struct; they live only in
 *  the @ref ParticleStore columns (single ownership, spec section 4). The
 *  struct is retained (empty) for the ghost-transfer documentation reference in
 *  ghosts.hpp; it is no longer a member of @ref Particle nor serialized.
 */
struct ParticlePosition {};

/** Force information on a particle. Forces of ghost particles are
 *  collected and added up to the force of the original particle.
 */
struct ParticleForce {
  ParticleForce() = default;
  ParticleForce(ParticleForce const &) = default;
  ParticleForce &operator=(ParticleForce const &) = default;
  ParticleForce(const Utils::Vector3d &f) : f(f) {}
#ifdef ESPRESSO_ROTATION
  ParticleForce(const Utils::Vector3d &f, const Utils::Vector3d &torque)
      : f(f), torque(torque) {}
#endif

  friend ParticleForce operator+(ParticleForce const &lhs,
                                 ParticleForce const &rhs) {
    ParticleForce result = lhs;
    result += rhs;
    return result;
  }

  ParticleForce &operator+=(ParticleForce const &rhs) {
    f += rhs.f;
#ifdef ESPRESSO_ROTATION
    torque += rhs.torque;
#endif
    return *this;
  }

  /** force. */
  Utils::Vector3d f = {0., 0., 0.};

#ifdef ESPRESSO_ROTATION
  /** torque. */
  Utils::Vector3d torque = {0., 0., 0.};
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & f;
#ifdef ESPRESSO_ROTATION
    ar & torque;
#endif
  }
};

/** Momentum information on a particle. Information communicated to calculate
 *  velocity-dependent interactions with ghost particles.
 *
 *  Migration phase 4: the velocity and angular-velocity fields have left this
 *  struct; they live only in the @ref ParticleStore columns (single ownership,
 *  spec section 4). The struct is retained (empty) for the ghost-transfer
 *  documentation reference in ghosts.hpp; it is no longer a member of @ref
 *  Particle nor serialized.
 */
struct ParticleMomentum {};

/** Information on a particle that is needed only on the
 *  node the particle belongs to.
 */
struct ParticleLocal {
  /** is particle a ghost particle. */
  bool ghost = false;

  // Migration phase 3: the Lees-Edwards flag/offset and the
  // position-from-the-last-Verlet-list-update (p_old) fields have left this
  // struct; they live only in the @ref ParticleStore columns (single
  // ownership, spec section 4). Only the ghost flag remains here.

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & ghost;
  }
};

#ifdef ESPRESSO_BOND_CONSTRAINT
/** Migration phase 6: the position/velocity RATTLE correction has left this
 *  struct; it lives only in a @ref ParticleStore observable column (single
 *  ownership, spec section 4). The struct is retained (empty) purely as a type
 *  anchor for the @ref GHOSTTRANS_RATTLE documentation reference in ghosts.hpp;
 *  it is no longer a member of @ref Particle nor part of @ref Particle
 *  serialization (it was never serialized), and the RATTLE ghost wire now
 *  archives one @ref Utils::Vector3d directly via @ref
 * Particle::rattle_correction().
 */
struct ParticleRattle {};
#endif

/** Struct holding all information for one particle. */
struct Particle { // NOLINT(bugprone-exception-escape)
private:
  ParticleLocal l;

  /** Transitional (migration phase 6): dual-role bonds/exclusions storage.
   *  These were the owned @c bl / @c el members; phase 6 evicts the
   *  authoritative copy into the @ref ParticleStore ragged host sidecars
   *  (single ownership, spec section 4). The struct members survive with the
   *  @c m_migration_ prefix in their SECOND role: they are the DETACHED storage
   *  (returned by @ref bonds() / @ref exclusions() when the particle is not
   *  attached to a store) AND the migration/fetch envelope carried by @ref
   *  Particle::serialize across the boost-serialized inter-rank exchange (which
   *  carries neither Kokkos columns nor host sidecars). The @c detached_*()
   *  getters read the sidecar when attached and the member when detached; @ref
   *  Particle::serialize syncs the members from them on SAVE so the envelope
   *  ferries the LIVE value, and @ref ParticleStore::assign_row seeds a
   *  new/migrated row's sidecar from the member via @ref migration_bonds() /
   *  @ref migration_exclusions(). Removed in phase 7 when the inter-rank
   *  exchange switches to per-field packing. */
  BondList m_migration_bonds;
#ifdef ESPRESSO_EXCLUSIONS
  /** list of particles, with which this particle has no non-bonded
   *  interactions (dual-role: detached storage + migration envelope, see
   *  @c m_migration_bonds above)
   */
  Utils::compact_vector<int> m_migration_exclusions;
#endif

  /** Transitional (migration phase 2): row of this particle in the
   *  ParticleStore, -1 while detached. Rank-local; never serialized. */
  ParticleStore *m_particle_store = nullptr;
  int m_store_row = -1;

  /** Transitional (migration phase 2): force/torque carrier used only while the
   *  particle is detached from a store (freshly constructed or just
   *  deserialized after a cross-rank migration). The @ref ParticleStore columns
   *  are the source of truth once attached; these members ferry the observable
   *  values across the boost-serialized inter-rank particle exchange (which
   * does not carry Kokkos columns) so that @ref ParticleStore::assign_row can
   * seed the new row. This preserves the pre-migration behavior where the force
   *  travelled with the migrating particle. Removed in phase 7, when the
   *  inter-rank exchange switches to per-field column packing. */
  Utils::Vector3d m_detached_force = {0., 0., 0.};
#ifdef ESPRESSO_ROTATION
  Utils::Vector3d m_detached_torque = {0., 0., 0.};
#endif

  /** Transitional (migration phase 3): STATE migration carriers, DORMANT until
   *  phase 8. Mirroring the force/torque carriers above, these members ferry
   *  the state fields (position, image box, quaternion, position-at-last-
   *  verlet-update, position-at-last-time-step, Lees-Edwards offset and flag)
   *  across the boost-serialized inter-rank exchange so @ref
   *  ParticleStore::assign_row can seed a migrated particle's new row. Until
   *  the phase-8 flip the corresponding sub-struct members remain authoritative
   *  and the migration_*() getters read them directly (see below); these
   *  members are NOT yet added to @ref Particle::serialize (phase 8 adds the
   *  @c ar & lines when the sub-structs leave, so the serialization format
   *  changes exactly once). Defaults match the sub-struct member defaults. */
  Utils::Vector3d m_migration_position = {0., 0., 0.};
  Utils::Vector3i m_migration_image_box = {0, 0, 0};
#ifdef ESPRESSO_ROTATION
  Utils::Quaternion<double> m_migration_quaternion =
      Utils::Quaternion<double>::identity();
#endif
  Utils::Vector3d m_migration_position_at_last_verlet_update = {0., 0., 0.};
#ifdef ESPRESSO_BOND_CONSTRAINT
  Utils::Vector3d m_migration_position_last_time_step = {0., 0., 0.};
#endif
  double m_migration_lees_edwards_offset = 0.;
  short int m_migration_lees_edwards_flag = 0;

  /** Transitional (migration phase 4): MOMENTUM migration carriers (now LIVE).
   *  Mirror of the phase-3 state carriers for velocity and angular velocity,
   *  which live in the @ref ParticleStore columns. These members ferry the
   *  velocity/angular-velocity across the boost-serialized inter-rank exchange
   *  (which does not carry Kokkos columns) so @ref ParticleStore::assign_row
   *  can seed a migrated particle's new row. @ref Particle::serialize fills
   * them from the detached_*() getters on SAVE; the migration_*() getters read
   * them raw (safe on a detached particle). Defaults match the pre-migration
   *  m.v/m.omega defaults (zero velocity for brand-new rows). Removed in phase
   * 7 when the inter-rank exchange switches to per-field column packing. */
  Utils::Vector3d m_migration_velocity = {0., 0., 0.};
#ifdef ESPRESSO_ROTATION
  Utils::Vector3d m_migration_angular_velocity = {0., 0., 0.};
#endif

  /** Transitional (migration phase 5): PARAMETER migration carriers (now LIVE).
   *  Mirror of the phase-3/4 state/momentum carriers for every parameter field
   *  that leaves the (now-deleted) @ref ParticleProperties member. These
   *  members ferry the parameters across the boost-serialized inter-rank
   *  exchange (which does not carry Kokkos columns / host sidecars) so @ref
   *  ParticleStore::assign_row can seed a migrated particle's new row. @ref
   *  Particle::serialize fills them from the detached_*() getters on SAVE; the
   *  migration_*() getters read them raw (safe on a detached particle).
   *  Defaults match the pre-flip ParticleProperties member defaults EXACTLY.
   *  The constexpr-when-disabled fields (mass/q/rinertia/rotation/ext_flag
   *  under their #else branches) have NO carrier: their accessors keep the
   *  static fallback (no column, nothing to ferry when the feature is off).
   *  Removed in phase 7 when the inter-rank exchange switches to per-field
   *  column packing. */
  int m_migration_id = -1;
  int m_migration_mol_id = 0;
  int m_migration_type = 0;
  int m_migration_propagation = PropagationMode::SYSTEM_DEFAULT;
#ifdef ESPRESSO_ROTATION
  std::uint8_t m_migration_rotation = static_cast<std::uint8_t>(0b000u);
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  std::uint8_t m_migration_ext_flag = static_cast<std::uint8_t>(0b000u);
#endif
#ifdef ESPRESSO_MASS
  double m_migration_mass = 1.0;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  double m_migration_q = 0.0;
#endif
#ifdef ESPRESSO_DIPOLES
  double m_migration_dipm = 0.0;
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Utils::Vector3d m_migration_rinertia = {1., 1., 1.};
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  Utils::Vector3d m_migration_mu_E = {0., 0., 0.};
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  Utils::Vector3d m_migration_dip_fld = {0., 0., 0.};
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  Utils::Vector3d m_migration_ext_force = {0., 0., 0.};
#ifdef ESPRESSO_ROTATION
  Utils::Vector3d m_migration_ext_torque = {0., 0., 0.};
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifndef ESPRESSO_PARTICLE_ANISOTROPY
  double m_migration_gamma = -1.;
#else
  Utils::Vector3d m_migration_gamma = {-1., -1., -1.};
#endif
#ifdef ESPRESSO_ROTATION
#ifndef ESPRESSO_PARTICLE_ANISOTROPY
  double m_migration_gamma_rot = -1.;
#else
  Utils::Vector3d m_migration_gamma_rot = {-1., -1., -1.};
#endif
#endif // ESPRESSO_ROTATION
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_ENGINE
  ParticleParametersSwimming m_migration_swimming{};
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  ThermalStonerWohlfarthParameters m_migration_magnetodynamics{};
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  VirtualSitesRelativeParameters m_migration_vs_relative{};
#endif

  /** Static fallbacks for the constexpr-when-disabled parameter accessors.
   *  When the feature is off there is no ParticleStore column and no migration
   *  carrier; the accessor returns a reference to this immutable default so the
   *  read-only accessor keeps its pre-migration constexpr semantics. */
#ifndef ESPRESSO_MASS
  static constexpr double mass_fallback{1.0};
#endif
#ifndef ESPRESSO_ROTATIONAL_INERTIA
  static constexpr Utils::Vector3d rinertia_fallback = {1., 1., 1.};
#endif
#ifndef ESPRESSO_ELECTROSTATICS
  static constexpr double q_fallback{0.0};
#endif

public:
  void attach_to_store(ParticleStore &store, int const row) {
    m_particle_store = &store;
    m_store_row = row;
  }
  auto store() const { return m_particle_store; }
  auto store_row() const { return m_store_row; }

  /** @brief Observable force as a plain value, valid whether the particle is
   *  attached to a store (reads the column) or detached (reads the migration
   *  carrier). Used by Particle serialization to capture the value to ferry. */
  Utils::Vector3d detached_force() const {
    return (m_particle_store != nullptr) ? force() : m_detached_force;
  }
#ifdef ESPRESSO_ROTATION
  Utils::Vector3d detached_torque() const {
    return (m_particle_store != nullptr) ? torque() : m_detached_torque;
  }
#endif
  /** @brief Raw migration carrier (never touches columns). Used by
   *  ParticleStore::assign_row to seed a migrated/new particle's row. */
  Utils::Vector3d const &migration_force() const { return m_detached_force; }
#ifdef ESPRESSO_ROTATION
  Utils::Vector3d const &migration_torque() const { return m_detached_torque; }
#endif

  /** @brief Phase-3 STATE migration carriers (now LIVE, mirroring the phase-2
   *  force carriers).
   *
   *  The @c detached_*() getters return the current value whether the particle
   *  is attached to a store (read the column) or detached (read the migration
   *  carrier). @ref Particle serialization fills the carriers from these
   *  getters on SAVE. The @c migration_*() getters return the raw carrier that
   *  @ref ParticleStore::assign_row seeds a new/migrated row from (never
   * touches a column, so it is safe on a detached particle whose column row is
   *  invalid). Removed in phase 7 when the inter-rank exchange switches to
   *  per-field column packing.
   *  @{ */
  Utils::Vector3d detached_position() const {
    return (m_particle_store != nullptr) ? pos() : m_migration_position;
  }
  Utils::Vector3d const &migration_position() const {
    return m_migration_position;
  }
  Utils::Vector3i detached_image_box() const {
    return (m_particle_store != nullptr) ? image_box() : m_migration_image_box;
  }
  Utils::Vector3i const &migration_image_box() const {
    return m_migration_image_box;
  }
#ifdef ESPRESSO_ROTATION
  Utils::Quaternion<double> detached_quaternion() const {
    return (m_particle_store != nullptr) ? quat() : m_migration_quaternion;
  }
  Utils::Quaternion<double> const &migration_quaternion() const {
    return m_migration_quaternion;
  }
#endif
  Utils::Vector3d detached_position_at_last_verlet_update() const {
    return (m_particle_store != nullptr)
               ? pos_at_last_verlet_update()
               : m_migration_position_at_last_verlet_update;
  }
  Utils::Vector3d const &migration_position_at_last_verlet_update() const {
    return m_migration_position_at_last_verlet_update;
  }
#ifdef ESPRESSO_BOND_CONSTRAINT
  Utils::Vector3d detached_position_last_time_step() const {
    return (m_particle_store != nullptr) ? pos_last_time_step()
                                         : m_migration_position_last_time_step;
  }
  Utils::Vector3d const &migration_position_last_time_step() const {
    return m_migration_position_last_time_step;
  }
#endif
  double detached_lees_edwards_offset() const {
    return (m_particle_store != nullptr) ? lees_edwards_offset()
                                         : m_migration_lees_edwards_offset;
  }
  double migration_lees_edwards_offset() const {
    return m_migration_lees_edwards_offset;
  }
  short int detached_lees_edwards_flag() const {
    return (m_particle_store != nullptr) ? lees_edwards_flag()
                                         : m_migration_lees_edwards_flag;
  }
  short int migration_lees_edwards_flag() const {
    return m_migration_lees_edwards_flag;
  }
  /** @} */

  /** @brief Phase-4 MOMENTUM migration carriers (now LIVE, mirroring the
   *  phase-3 state carriers above).
   *
   *  The @c detached_*() getters return the current value whether the particle
   *  is attached to a store (read the column) or detached (read the migration
   *  carrier). @ref Particle serialization fills the carriers from these
   * getters on SAVE. The @c migration_*() getters return the raw carrier that
   * @ref ParticleStore::assign_row seeds a new/migrated row from (never touches
   * a column, so it is safe on a detached particle whose column row is
   * invalid). Removed in phase 7 when the inter-rank exchange switches to
   * per-field column packing.
   *  @{ */
  Utils::Vector3d detached_velocity() const {
    return (m_particle_store != nullptr) ? v() : m_migration_velocity;
  }
  Utils::Vector3d const &migration_velocity() const {
    return m_migration_velocity;
  }
#ifdef ESPRESSO_ROTATION
  Utils::Vector3d detached_angular_velocity() const {
    return (m_particle_store != nullptr) ? omega()
                                         : m_migration_angular_velocity;
  }
  Utils::Vector3d const &migration_angular_velocity() const {
    return m_migration_angular_velocity;
  }
#endif
  /** @} */

  /** @brief Phase-5 PARAMETER migration carriers (now LIVE, mirroring the
   *  phase-3/4 state/momentum carriers above).
   *
   *  The @c detached_*() getters return the current value whether the particle
   *  is attached to a store (read the column / sidecar) or detached (read the
   *  migration carrier). @ref Particle serialization fills the carriers from
   *  these getters on SAVE. The @c migration_*() getters return the raw carrier
   *  that @ref ParticleStore::assign_row seeds a new/migrated row from (never
   *  touches a column, so it is safe on a detached particle whose column row is
   *  invalid). The constexpr-when-disabled fields keep their static fallbacks
   *  under #else and have no carrier. Removed in phase 7 when the inter-rank
   *  exchange switches to per-field column packing.
   *  @{ */
  int detached_id() const {
    return (m_particle_store != nullptr) ? id() : m_migration_id;
  }
  int migration_id() const { return m_migration_id; }
  int detached_mol_id() const {
    return (m_particle_store != nullptr) ? mol_id() : m_migration_mol_id;
  }
  int migration_mol_id() const { return m_migration_mol_id; }
  int detached_type() const {
    return (m_particle_store != nullptr) ? type() : m_migration_type;
  }
  int migration_type() const { return m_migration_type; }
  int detached_propagation() const {
    return (m_particle_store != nullptr) ? propagation()
                                         : m_migration_propagation;
  }
  int migration_propagation() const { return m_migration_propagation; }
#ifdef ESPRESSO_ROTATION
  std::uint8_t detached_rotation() const {
    return (m_particle_store != nullptr) ? rotation() : m_migration_rotation;
  }
  std::uint8_t migration_rotation() const { return m_migration_rotation; }
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  std::uint8_t detached_ext_flag() const {
    return (m_particle_store != nullptr) ? fixed() : m_migration_ext_flag;
  }
  std::uint8_t migration_ext_flag() const { return m_migration_ext_flag; }
#endif
#ifdef ESPRESSO_MASS
  double detached_mass() const {
    return (m_particle_store != nullptr) ? mass() : m_migration_mass;
  }
  double migration_mass() const { return m_migration_mass; }
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  double detached_q() const {
    return (m_particle_store != nullptr) ? q() : m_migration_q;
  }
  double migration_q() const { return m_migration_q; }
#endif
#ifdef ESPRESSO_DIPOLES
  double detached_dipm() const {
    return (m_particle_store != nullptr) ? dipm() : m_migration_dipm;
  }
  double migration_dipm() const { return m_migration_dipm; }
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Utils::Vector3d detached_rinertia() const {
    return (m_particle_store != nullptr) ? rinertia() : m_migration_rinertia;
  }
  Utils::Vector3d const &migration_rinertia() const {
    return m_migration_rinertia;
  }
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  Utils::Vector3d detached_mu_E() const {
    return (m_particle_store != nullptr) ? mu_E() : m_migration_mu_E;
  }
  Utils::Vector3d const &migration_mu_E() const { return m_migration_mu_E; }
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  Utils::Vector3d detached_dip_fld() const {
    return (m_particle_store != nullptr) ? dip_fld() : m_migration_dip_fld;
  }
  Utils::Vector3d const &migration_dip_fld() const {
    return m_migration_dip_fld;
  }
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  Utils::Vector3d detached_ext_force() const {
    return (m_particle_store != nullptr) ? ext_force() : m_migration_ext_force;
  }
  Utils::Vector3d const &migration_ext_force() const {
    return m_migration_ext_force;
  }
#ifdef ESPRESSO_ROTATION
  Utils::Vector3d detached_ext_torque() const {
    return (m_particle_store != nullptr) ? ext_torque()
                                         : m_migration_ext_torque;
  }
  Utils::Vector3d const &migration_ext_torque() const {
    return m_migration_ext_torque;
  }
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  auto detached_gamma() const {
    return (m_particle_store != nullptr) ? gamma() : m_migration_gamma;
  }
  auto const &migration_gamma() const { return m_migration_gamma; }
#ifdef ESPRESSO_ROTATION
  auto detached_gamma_rot() const {
    return (m_particle_store != nullptr) ? gamma_rot() : m_migration_gamma_rot;
  }
  auto const &migration_gamma_rot() const { return m_migration_gamma_rot; }
#endif
#endif
#ifdef ESPRESSO_ENGINE
  ParticleParametersSwimming detached_swimming() const {
    return (m_particle_store != nullptr) ? swimming() : m_migration_swimming;
  }
  ParticleParametersSwimming const &migration_swimming() const {
    return m_migration_swimming;
  }
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  ThermalStonerWohlfarthParameters detached_magnetodynamics() const {
    return (m_particle_store != nullptr) ? magnetodynamics()
                                         : m_migration_magnetodynamics;
  }
  ThermalStonerWohlfarthParameters const &migration_magnetodynamics() const {
    return m_migration_magnetodynamics;
  }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  VirtualSitesRelativeParameters detached_vs_relative() const {
    return (m_particle_store != nullptr) ? vs_relative()
                                         : m_migration_vs_relative;
  }
  VirtualSitesRelativeParameters const &migration_vs_relative() const {
    return m_migration_vs_relative;
  }
#endif

  /** @brief Phase-6 RAGGED migration seed / detached getters (now LIVE).
   *
   *  Bonds/exclusions are the last owned non-POD members; phase 6 evicts the
   *  authoritative copy into ParticleStore ragged sidecars. Unlike the POD
   *  carriers, no separate carrier member is introduced (the dual-role design
   *  of exploration finding 3 keeps the struct members themselves as the
   *  detached storage / the migration envelope). The @c detached_*() getters
   *  return the CURRENT value whether attached (read the sidecar) or detached
   *  (read the member); @ref Particle::serialize fills the members from them on
   *  SAVE. The @c migration_*() getters return the raw member that @ref
   *  ParticleStore::assign_row seeds a new/migrated row's sidecar from (never
   *  touches a sidecar, so it is safe on a detached particle whose store row is
   *  invalid). Removed in phase 7 when the inter-rank exchange switches to
   *  per-field packing.
   *  @{ */
  BondList detached_bonds() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->bonds_sidecar_reference(m_store_row)
               : m_migration_bonds;
  }
  BondList const &migration_bonds() const { return m_migration_bonds; }
#ifdef ESPRESSO_EXCLUSIONS
  Utils::compact_vector<int> detached_exclusions() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->exclusions_sidecar_reference(m_store_row)
               : m_migration_exclusions;
  }
  Utils::compact_vector<int> const &migration_exclusions() const {
    return m_migration_exclusions;
  }
#endif
  /** @} */

  int const &id() const {
    return (m_particle_store != nullptr) ? m_particle_store->id(m_store_row)
                                         : m_migration_id;
  }
  int &id() {
    return (m_particle_store != nullptr) ? m_particle_store->id(m_store_row)
                                         : m_migration_id;
  }
  int const &mol_id() const {
    return (m_particle_store != nullptr) ? m_particle_store->mol_id(m_store_row)
                                         : m_migration_mol_id;
  }
  int &mol_id() {
    return (m_particle_store != nullptr) ? m_particle_store->mol_id(m_store_row)
                                         : m_migration_mol_id;
  }
  int const &type() const {
    return (m_particle_store != nullptr) ? m_particle_store->type(m_store_row)
                                         : m_migration_type;
  }
  int &type() {
    return (m_particle_store != nullptr) ? m_particle_store->type(m_store_row)
                                         : m_migration_type;
  }

  int const &propagation() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->propagation(m_store_row)
               : m_migration_propagation;
  }
  int &propagation() {
    return (m_particle_store != nullptr)
               ? m_particle_store->propagation(m_store_row)
               : m_migration_propagation;
  }

  bool operator==(Particle const &rhs) const { return id() == rhs.id(); }

  bool operator!=(Particle const &rhs) const { return id() != rhs.id(); }

  BondList const &bonds() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->bonds_sidecar_reference(m_store_row)
               : m_migration_bonds;
  }
  BondList &bonds() {
    return (m_particle_store != nullptr)
               ? m_particle_store->bonds_sidecar_reference(m_store_row)
               : m_migration_bonds;
  }

  Utils::Vector3d pos() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->position_value(m_store_row)
               : m_migration_position;
  }
  VectorReference pos() {
    return (m_particle_store != nullptr)
               ? m_particle_store->position_reference(m_store_row)
               : VectorReference(m_migration_position.data(), 1u);
  }
  Utils::Vector3d v() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->velocity_value(m_store_row)
               : m_migration_velocity;
  }
  VectorReference v() {
    return (m_particle_store != nullptr)
               ? m_particle_store->velocity_reference(m_store_row)
               : VectorReference(m_migration_velocity.data(), 1u);
  }
  auto force() {
    assert(m_particle_store != nullptr);
    return m_particle_store->force_reference(m_store_row);
  }
  Utils::Vector3d force() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->force_value(m_store_row);
  }

  bool is_ghost() const { return l.ghost; }
  void set_ghost(bool const ghost_flag) { l.ghost = ghost_flag; }
  VectorReference pos_at_last_verlet_update() {
    return (m_particle_store != nullptr)
               ? m_particle_store->position_at_last_verlet_update_reference(
                     m_store_row)
               : VectorReference(
                     m_migration_position_at_last_verlet_update.data(), 1u);
  }
  Utils::Vector3d pos_at_last_verlet_update() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->position_at_last_verlet_update_value(
                     m_store_row)
               : m_migration_position_at_last_verlet_update;
  }
  Utils::Vector3i image_box() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->image_box_value(m_store_row)
               : m_migration_image_box;
  }
  IntegerVectorReference image_box() {
    return (m_particle_store != nullptr)
               ? m_particle_store->image_box_reference(m_store_row)
               : IntegerVectorReference(m_migration_image_box.data(), 1u);
  }
  double lees_edwards_offset() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->lees_edwards_offset(m_store_row)
               : m_migration_lees_edwards_offset;
  }
  double &lees_edwards_offset() {
    return (m_particle_store != nullptr)
               ? m_particle_store->lees_edwards_offset(m_store_row)
               : m_migration_lees_edwards_offset;
  }
  short lees_edwards_flag() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->lees_edwards_flag(m_store_row)
               : m_migration_lees_edwards_flag;
  }
  short &lees_edwards_flag() {
    return (m_particle_store != nullptr)
               ? m_particle_store->lees_edwards_flag(m_store_row)
               : m_migration_lees_edwards_flag;
  }

#ifdef ESPRESSO_MASS
  double const &mass() const {
    return (m_particle_store != nullptr) ? m_particle_store->mass(m_store_row)
                                         : m_migration_mass;
  }
  double &mass() {
    return (m_particle_store != nullptr) ? m_particle_store->mass(m_store_row)
                                         : m_migration_mass;
  }
#else
  constexpr auto &mass() const { return mass_fallback; }
#endif
#ifdef ESPRESSO_ROTATION
  std::uint8_t const &rotation() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->rotation(m_store_row)
               : m_migration_rotation;
  }
  std::uint8_t &rotation() {
    return (m_particle_store != nullptr)
               ? m_particle_store->rotation(m_store_row)
               : m_migration_rotation;
  }
  bool can_rotate() const { return static_cast<bool>(rotation()); }
  bool can_rotate_around(unsigned int const axis) const {
    assert(axis <= 2u);
    return detail::get_nth_bit(rotation(), axis);
  }
  void set_can_rotate_around(unsigned int const axis, bool const rot_flag) {
    assert(axis <= 2u);
    if (rot_flag) {
      rotation() |= static_cast<uint8_t>(1u << axis);
    } else {
      rotation() &= static_cast<uint8_t>(~(1u << axis));
    }
  }
  void set_can_rotate_all_axes() { rotation() = static_cast<uint8_t>(0b111u); }
  void set_cannot_rotate_all_axes() {
    rotation() = static_cast<uint8_t>(0b000u);
  }
  Utils::Quaternion<double> quat() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->quaternion_value(m_store_row)
               : m_migration_quaternion;
  }
  QuaternionReference quat() {
    return (m_particle_store != nullptr)
               ? m_particle_store->quaternion_reference(m_store_row)
               : QuaternionReference(m_migration_quaternion.data(), 1u);
  }
  auto torque() {
    assert(m_particle_store != nullptr);
    return m_particle_store->torque_reference(m_store_row);
  }
  Utils::Vector3d torque() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->torque_value(m_store_row);
  }
  Utils::Vector3d omega() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->angular_velocity_value(m_store_row)
               : m_migration_angular_velocity;
  }
  VectorReference omega() {
    return (m_particle_store != nullptr)
               ? m_particle_store->angular_velocity_reference(m_store_row)
               : VectorReference(m_migration_angular_velocity.data(), 1u);
  }
#ifdef ESPRESSO_EXTERNAL_FORCES
  Utils::Vector3d ext_torque() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->ext_torque_value(m_store_row)
               : m_migration_ext_torque;
  }
  VectorReference ext_torque() {
    return (m_particle_store != nullptr)
               ? m_particle_store->ext_torque_reference(m_store_row)
               : VectorReference(m_migration_ext_torque.data(), 1u);
  }
#endif // ESPRESSO_EXTERNAL_FORCES
  auto calc_director() const {
    return Utils::convert_quaternion_to_director(quat());
  }
#else  // ESPRESSO_ROTATION
  auto can_rotate() const { return false; }
  auto can_rotate_around(unsigned int const) const { return false; }
#endif // ESPRESSO_ROTATION
#ifdef ESPRESSO_DIPOLES
  double const &dipm() const {
    return (m_particle_store != nullptr) ? m_particle_store->dipm(m_store_row)
                                         : m_migration_dipm;
  }
  double &dipm() {
    return (m_particle_store != nullptr) ? m_particle_store->dipm(m_store_row)
                                         : m_migration_dipm;
  }
  auto calc_dip() const { return calc_director() * dipm(); }
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  // The Stoner-Wohlfarth parameters live in the ParticleStore host sidecar (a
  // whole POD indexed by store row) when attached, or in the migration carrier
  // when detached. The accessors return a reference into whichever holds the
  // POD; the individual-field accessors then reference into that.
  ThermalStonerWohlfarthParameters &magnetodynamics() {
    return (m_particle_store != nullptr)
               ? m_particle_store->magnetodynamics(m_store_row)
               : m_migration_magnetodynamics;
  }
  ThermalStonerWohlfarthParameters const &magnetodynamics() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->magnetodynamics(m_store_row)
               : m_migration_magnetodynamics;
  }
  auto const &stoner_wohlfarth_is_enabled() const {
    return magnetodynamics().is_enabled;
  }
  auto &stoner_wohlfarth_is_enabled() { return magnetodynamics().is_enabled; }
  auto const &stoner_wohlfarth_phi_0() const { return magnetodynamics().phi0; }
  auto &stoner_wohlfarth_phi_0() { return magnetodynamics().phi0; }
  auto const &saturation_magnetization() const {
    return magnetodynamics().sat_mag;
  }
  auto &saturation_magnetization() { return magnetodynamics().sat_mag; }
  auto const &magnetic_anisotropy_field_inv() const {
    return magnetodynamics().ani_fld_inv;
  }
  auto &magnetic_anisotropy_field_inv() {
    return magnetodynamics().ani_fld_inv;
  }
  auto const &magnetic_anisotropy_energy() const {
    return magnetodynamics().ani_energy;
  }
  auto &magnetic_anisotropy_energy() { return magnetodynamics().ani_energy; }
  auto const &stoner_wohlfarth_tau0_inv() const {
    return magnetodynamics().tau0_inv;
  }
  auto &stoner_wohlfarth_tau0_inv() { return magnetodynamics().tau0_inv; }
  auto const &stoner_wohlfarth_dt_incr() const {
    return magnetodynamics().dt_incr;
  }
  auto &stoner_wohlfarth_dt_incr() { return magnetodynamics().dt_incr; }
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  Utils::Vector3d dip_fld() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->dip_fld_value(m_store_row)
               : m_migration_dip_fld;
  }
  VectorReference dip_fld() {
    return (m_particle_store != nullptr)
               ? m_particle_store->dip_fld_reference(m_store_row)
               : VectorReference(m_migration_dip_fld.data(), 1u);
  }
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Utils::Vector3d rinertia() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->rinertia_value(m_store_row)
               : m_migration_rinertia;
  }
  VectorReference rinertia() {
    return (m_particle_store != nullptr)
               ? m_particle_store->rinertia_reference(m_store_row)
               : VectorReference(m_migration_rinertia.data(), 1u);
  }
#else
  constexpr auto &rinertia() const { return rinertia_fallback; }
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  double const &q() const {
    return (m_particle_store != nullptr) ? m_particle_store->q(m_store_row)
                                         : m_migration_q;
  }
  double &q() {
    return (m_particle_store != nullptr) ? m_particle_store->q(m_store_row)
                                         : m_migration_q;
  }
#else
  constexpr auto &q() const { return q_fallback; }
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  Utils::Vector3d mu_E() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->mu_E_value(m_store_row)
               : m_migration_mu_E;
  }
  VectorReference mu_E() {
    return (m_particle_store != nullptr)
               ? m_particle_store->mu_E_reference(m_store_row)
               : VectorReference(m_migration_mu_E.data(), 1u);
  }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES
  auto is_virtual() const {
    return (propagation() & (PropagationMode::TRANS_VS_RELATIVE |
                             PropagationMode::TRANS_VS_CENTER_OF_MASS |
                             PropagationMode::ROT_VS_RELATIVE |
                             PropagationMode::ROT_VS_INDEPENDENT |
                             PropagationMode::TRANS_LB_TRACER)) != 0;
  }
#else
  constexpr auto is_virtual() const { return false; }
#endif // ESPRESSO_VIRTUAL_SITES
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  VirtualSitesRelativeParameters &vs_relative() {
    return (m_particle_store != nullptr)
               ? m_particle_store->vs_relative(m_store_row)
               : m_migration_vs_relative;
  }
  VirtualSitesRelativeParameters const &vs_relative() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->vs_relative(m_store_row)
               : m_migration_vs_relative;
  }
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  // gamma/gamma_rot: element reference (double&) when isotropic, or a
  // VectorReference when ESPRESSO_PARTICLE_ANISOTROPY selects per-axis
  // friction; const returns a value. Detached falls back to the carrier.
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  Utils::Vector3d gamma() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->gamma_value(m_store_row)
               : m_migration_gamma;
  }
  VectorReference gamma() {
    return (m_particle_store != nullptr)
               ? m_particle_store->gamma_reference(m_store_row)
               : VectorReference(m_migration_gamma.data(), 1u);
  }
#else
  double gamma() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->gamma_value(m_store_row)
               : m_migration_gamma;
  }
  double &gamma() {
    return (m_particle_store != nullptr)
               ? m_particle_store->gamma_reference(m_store_row)
               : m_migration_gamma;
  }
#endif // ESPRESSO_PARTICLE_ANISOTROPY
#ifdef ESPRESSO_ROTATION
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  Utils::Vector3d gamma_rot() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->gamma_rot_value(m_store_row)
               : m_migration_gamma_rot;
  }
  VectorReference gamma_rot() {
    return (m_particle_store != nullptr)
               ? m_particle_store->gamma_rot_reference(m_store_row)
               : VectorReference(m_migration_gamma_rot.data(), 1u);
  }
#else
  double gamma_rot() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->gamma_rot_value(m_store_row)
               : m_migration_gamma_rot;
  }
  double &gamma_rot() {
    return (m_particle_store != nullptr)
               ? m_particle_store->gamma_rot_reference(m_store_row)
               : m_migration_gamma_rot;
  }
#endif // ESPRESSO_PARTICLE_ANISOTROPY
#endif // ESPRESSO_ROTATION
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_EXTERNAL_FORCES
  std::uint8_t const &fixed() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->ext_flag(m_store_row)
               : m_migration_ext_flag;
  }
  std::uint8_t &fixed() {
    return (m_particle_store != nullptr)
               ? m_particle_store->ext_flag(m_store_row)
               : m_migration_ext_flag;
  }
  bool has_fixed_coordinates() const { return static_cast<bool>(fixed()); }
  /** @brief Raw fixed-coordinate bitfield, by value.
   *  Hot integrator loops read this ONCE per particle and bit-test locally
   *  instead of calling @ref is_fixed_along (which re-reads the ParticleStore
   *  ext_flag column on every axis). Bit @c axis set means the coordinate is
   *  fixed. */
  std::uint8_t fixed_flags_byte() const { return fixed(); }
  bool is_fixed_along(unsigned int const axis) const {
    assert(axis <= 2u);
    return detail::get_nth_bit(fixed(), axis);
  }
  void set_fixed_along(int const axis, bool const fixed_flag) {
    // set new flag
    if (fixed_flag) {
      fixed() |= static_cast<uint8_t>(1u << axis);
    } else {
      fixed() &= static_cast<uint8_t>(~(1u << axis));
    }
  }
  Utils::Vector3d ext_force() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->ext_force_value(m_store_row)
               : m_migration_ext_force;
  }
  VectorReference ext_force() {
    return (m_particle_store != nullptr)
               ? m_particle_store->ext_force_reference(m_store_row)
               : VectorReference(m_migration_ext_force.data(), 1u);
  }
#else  // ESPRESSO_EXTERNAL_FORCES
  constexpr bool has_fixed_coordinates() const { return false; }
  constexpr std::uint8_t fixed_flags_byte() const {
    return static_cast<std::uint8_t>(0u);
  }
  constexpr bool is_fixed_along(unsigned int const) const { return false; }
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_ENGINE
  ParticleParametersSwimming &swimming() {
    return (m_particle_store != nullptr)
               ? m_particle_store->swimming(m_store_row)
               : m_migration_swimming;
  }
  ParticleParametersSwimming const &swimming() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->swimming(m_store_row)
               : m_migration_swimming;
  }
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  Utils::Vector3d pos_last_time_step() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->position_last_time_step_value(m_store_row)
               : m_migration_position_last_time_step;
  }
  VectorReference pos_last_time_step() {
    return (m_particle_store != nullptr)
               ? m_particle_store->position_last_time_step_reference(
                     m_store_row)
               : VectorReference(m_migration_position_last_time_step.data(),
                                 1u);
  }
  // RATTLE/SHAKE correction (phase 6): evicted to a ParticleStore observable
  // column (structurally like force). It is per-iteration SHAKE scratch --
  // never persisted, never migrated, never checkpointed -- so there is no
  // migration carrier and the accessor is attached-only (asserts a store),
  // exactly like force().
  VectorReference rattle_correction() {
    assert(m_particle_store != nullptr);
    return m_particle_store->rattle_correction_reference(m_store_row);
  }
  Utils::Vector3d rattle_correction() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->rattle_correction_value(m_store_row);
  }
#endif

#ifdef ESPRESSO_EXCLUSIONS
  Utils::compact_vector<int> &exclusions() {
    return (m_particle_store != nullptr)
               ? m_particle_store->exclusions_sidecar_reference(m_store_row)
               : m_migration_exclusions;
  }
  Utils::compact_vector<int> const &exclusions() const {
    return (m_particle_store != nullptr)
               ? m_particle_store->exclusions_sidecar_reference(m_store_row)
               : m_migration_exclusions;
  }
  bool has_exclusion(int pid) const {
    return std::ranges::find(exclusions(), pid) != exclusions().end();
  }
#endif

private:
  friend boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & l;
    // Migration phase 6: bonds/exclusions live authoritatively in the
    // ParticleStore ragged host sidecars when the particle is attached. The
    // struct members are dual-role (detached storage + migration/fetch
    // envelope). On SAVE, sync the members from the sidecar (detached_*() reads
    // the sidecar when attached) BEFORE serializing them, so the envelope
    // carries the LIVE value -- this must run before these legs, unlike the
    // other carriers whose legs sit after the shared is_saving block below. The
    // leg form/order is byte-identical to the pre-flip `ar & bl` / `ar & el`.
    if (Archive::is_saving::value) {
      m_migration_bonds = detached_bonds();
#ifdef ESPRESSO_EXCLUSIONS
      m_migration_exclusions = detached_exclusions();
#endif
    }
    ar & m_migration_bonds;
#ifdef ESPRESSO_EXCLUSIONS
    ar & m_migration_exclusions;
#endif
    // Migration phases 2, 3, 4 & 5: force/torque (phase 2), the STATE fields
    // (position, image box, quaternion, position-at-last-verlet-update,
    // position-at-last-time-step, Lees-Edwards offset and flag; phase 3), the
    // MOMENTUM fields (velocity and angular velocity; phase 4), and ALL the
    // PARAMETER fields (id, mol_id, type, propagation, the rotation/ext_flag
    // bitfields, mass, rinertia, q, mu_E, dipm, dip_fld, gamma/gamma_rot,
    // ext_force/ext_torque and the three cold PODs swim/magnetodynamics/
    // vs_relative; phase 5) live in ParticleStore columns/sidecars, which are
    // not carried by this (boost) serializer used for cross-rank particle
    // exchange. Ferry the values so they survive a global resort that moves the
    // particle to another rank (matching pre-migration behavior). SAVE reads
    // the current value (column/sidecar if attached, carrier otherwise); LOAD
    // lands in the detached carrier, from which ParticleStore::assign_row seeds
    // the rebuilt row. The constexpr-when-disabled fields have no carrier and
    // are not serialized (their static fallback is the only value when the
    // feature is off).
    if (Archive::is_saving::value) {
      m_detached_force = detached_force();
#ifdef ESPRESSO_ROTATION
      m_detached_torque = detached_torque();
#endif
      m_migration_position = detached_position();
      m_migration_image_box = detached_image_box();
#ifdef ESPRESSO_ROTATION
      m_migration_quaternion = detached_quaternion();
#endif
      m_migration_position_at_last_verlet_update =
          detached_position_at_last_verlet_update();
#ifdef ESPRESSO_BOND_CONSTRAINT
      m_migration_position_last_time_step = detached_position_last_time_step();
#endif
      m_migration_lees_edwards_offset = detached_lees_edwards_offset();
      m_migration_lees_edwards_flag = detached_lees_edwards_flag();
      m_migration_velocity = detached_velocity();
#ifdef ESPRESSO_ROTATION
      m_migration_angular_velocity = detached_angular_velocity();
#endif
      m_migration_id = detached_id();
      m_migration_mol_id = detached_mol_id();
      m_migration_type = detached_type();
      m_migration_propagation = detached_propagation();
#ifdef ESPRESSO_ROTATION
      m_migration_rotation = detached_rotation();
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
      m_migration_ext_flag = detached_ext_flag();
#endif
#ifdef ESPRESSO_MASS
      m_migration_mass = detached_mass();
#endif
#ifdef ESPRESSO_ELECTROSTATICS
      m_migration_q = detached_q();
#endif
#ifdef ESPRESSO_DIPOLES
      m_migration_dipm = detached_dipm();
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
      m_migration_rinertia = detached_rinertia();
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
      m_migration_mu_E = detached_mu_E();
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
      m_migration_dip_fld = detached_dip_fld();
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
      m_migration_ext_force = detached_ext_force();
#ifdef ESPRESSO_ROTATION
      m_migration_ext_torque = detached_ext_torque();
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
      m_migration_gamma = detached_gamma();
#ifdef ESPRESSO_ROTATION
      m_migration_gamma_rot = detached_gamma_rot();
#endif
#endif
#ifdef ESPRESSO_ENGINE
      m_migration_swimming = detached_swimming();
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
      m_migration_magnetodynamics = detached_magnetodynamics();
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
      m_migration_vs_relative = detached_vs_relative();
#endif
    }
    ar & m_detached_force;
#ifdef ESPRESSO_ROTATION
    ar & m_detached_torque;
#endif
    ar & m_migration_position;
    ar & m_migration_image_box;
#ifdef ESPRESSO_ROTATION
    ar & m_migration_quaternion;
#endif
    ar & m_migration_position_at_last_verlet_update;
#ifdef ESPRESSO_BOND_CONSTRAINT
    ar & m_migration_position_last_time_step;
#endif
    ar & m_migration_lees_edwards_offset;
    ar & m_migration_lees_edwards_flag;
    ar & m_migration_velocity;
#ifdef ESPRESSO_ROTATION
    ar & m_migration_angular_velocity;
#endif
    ar & m_migration_id;
    ar & m_migration_mol_id;
    ar & m_migration_type;
    ar & m_migration_propagation;
#ifdef ESPRESSO_ROTATION
    ar & m_migration_rotation;
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    ar & m_migration_ext_flag;
#endif
#ifdef ESPRESSO_MASS
    ar & m_migration_mass;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    ar & m_migration_q;
#endif
#ifdef ESPRESSO_DIPOLES
    ar & m_migration_dipm;
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    ar & m_migration_rinertia;
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
    ar & m_migration_mu_E;
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    ar & m_migration_dip_fld;
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    ar & m_migration_ext_force;
#ifdef ESPRESSO_ROTATION
    ar & m_migration_ext_torque;
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    ar & m_migration_gamma;
#ifdef ESPRESSO_ROTATION
    ar & m_migration_gamma_rot;
#endif
#endif
#ifdef ESPRESSO_ENGINE
    ar & m_migration_swimming;
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
    ar & m_migration_magnetodynamics;
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    ar & m_migration_vs_relative;
#endif
  }
};

BOOST_CLASS_IMPLEMENTATION(Particle, object_serializable)
#ifdef ESPRESSO_ENGINE
BOOST_CLASS_IMPLEMENTATION(ParticleParametersSwimming, object_serializable)
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
BOOST_CLASS_IMPLEMENTATION(ThermalStonerWohlfarthParameters,
                           object_serializable)
#endif
BOOST_CLASS_IMPLEMENTATION(ParticleProperties, object_serializable)
BOOST_CLASS_IMPLEMENTATION(ParticlePosition, object_serializable)
BOOST_CLASS_IMPLEMENTATION(ParticleMomentum, object_serializable)
BOOST_CLASS_IMPLEMENTATION(ParticleForce, object_serializable)
BOOST_CLASS_IMPLEMENTATION(ParticleLocal, object_serializable)
#ifdef ESPRESSO_BOND_CONSTRAINT
BOOST_CLASS_IMPLEMENTATION(ParticleRattle, object_serializable)
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
BOOST_CLASS_IMPLEMENTATION(decltype(ParticleProperties::vs_relative),
                           object_serializable)
#endif

#ifdef ESPRESSO_ENGINE
BOOST_IS_BITWISE_SERIALIZABLE(ParticleParametersSwimming)
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
BOOST_IS_BITWISE_SERIALIZABLE(ThermalStonerWohlfarthParameters)
#endif
BOOST_IS_BITWISE_SERIALIZABLE(ParticleProperties)
BOOST_IS_BITWISE_SERIALIZABLE(ParticlePosition)
BOOST_IS_BITWISE_SERIALIZABLE(ParticleMomentum)
BOOST_IS_BITWISE_SERIALIZABLE(ParticleForce)
BOOST_IS_BITWISE_SERIALIZABLE(ParticleLocal)
#ifdef ESPRESSO_BOND_CONSTRAINT
BOOST_IS_BITWISE_SERIALIZABLE(ParticleRattle)
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
BOOST_IS_BITWISE_SERIALIZABLE(decltype(ParticleProperties::vs_relative))
#endif
