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
struct ParticleRattle {
  /** position/velocity correction */
  Utils::Vector3d correction = {0., 0., 0.};

  friend ParticleRattle operator+(ParticleRattle const &lhs,
                                  ParticleRattle const &rhs) {
    return {lhs.correction + rhs.correction};
  }

  ParticleRattle &operator+=(ParticleRattle const &rhs) {
    return *this = *this + rhs;
  }

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & correction;
  }
};
#endif

/** Struct holding all information for one particle. */
struct Particle { // NOLINT(bugprone-exception-escape)
private:
  ParticleProperties p;
  ParticleLocal l;
#ifdef ESPRESSO_BOND_CONSTRAINT
  ParticleRattle rattle;
#endif
  BondList bl;
#ifdef ESPRESSO_EXCLUSIONS
  /** list of particles, with which this particle has no non-bonded
   *  interactions
   */
  Utils::compact_vector<int> el;
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

  /** @brief Phase-5 PARAMETER migration carriers, DORMANT until the phase-4
   *  parameter flip (this migration's THE FLIP task).
   *
   *  Pre-flip the @ref ParticleProperties member @c p (and the cold parameter
   *  PODs it holds) is still authoritative and is carried across the
   *  boost-serialized inter-rank exchange by @c ar & p, so no dedicated carrier
   *  members are needed yet. These getters read the struct fields directly and
   *  are safe whether the particle is attached or detached (parameters do not
   *  live in the store columns yet). @ref ParticleStore::assign_row uses them
   *  to seed a new/migrated row so the new parameter columns/sidecars already
   *  carry the correct values before the flip makes them authoritative. At the
   *  flip these turn into real detached carriers wired into serialize().
   *  @{ */
  int migration_id() const { return p.identity; }
  int migration_mol_id() const { return p.mol_id; }
  int migration_type() const { return p.type; }
  int migration_propagation() const { return p.propagation; }
#ifdef ESPRESSO_ROTATION
  std::uint8_t migration_rotation() const { return p.rotation; }
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  std::uint8_t migration_ext_flag() const { return p.ext_flag; }
#endif
#ifdef ESPRESSO_MASS
  double migration_mass() const { return p.mass; }
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  double migration_q() const { return p.q; }
#endif
#ifdef ESPRESSO_DIPOLES
  double migration_dipm() const { return p.dipm; }
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Utils::Vector3d const &migration_rinertia() const { return p.rinertia; }
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  Utils::Vector3d const &migration_mu_E() const { return p.mu_E; }
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  Utils::Vector3d const &migration_dip_fld() const { return p.dip_fld; }
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  Utils::Vector3d const &migration_ext_force() const { return p.ext_force; }
#ifdef ESPRESSO_ROTATION
  Utils::Vector3d const &migration_ext_torque() const { return p.ext_torque; }
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  auto const &migration_gamma() const { return p.gamma; }
#ifdef ESPRESSO_ROTATION
  auto const &migration_gamma_rot() const { return p.gamma_rot; }
#endif
#endif
#ifdef ESPRESSO_ENGINE
  ParticleParametersSwimming const &migration_swimming() const {
    return p.swim;
  }
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  ThermalStonerWohlfarthParameters const &migration_magnetodynamics() const {
    return p.magnetodynamics;
  }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  VirtualSitesRelativeParameters const &migration_vs_relative() const {
    return p.vs_relative;
  }
#endif
  /** @} */

  auto const &id() const { return p.identity; }
  auto &id() { return p.identity; }
  auto const &mol_id() const { return p.mol_id; }
  auto &mol_id() { return p.mol_id; }
  auto const &type() const { return p.type; }
  auto &type() { return p.type; }

  auto const &propagation() const { return p.propagation; }
  auto &propagation() { return p.propagation; }

  bool operator==(Particle const &rhs) const { return id() == rhs.id(); }

  bool operator!=(Particle const &rhs) const { return id() != rhs.id(); }

  auto const &bonds() const { return bl; }
  auto &bonds() { return bl; }

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
  auto const &mass() const { return p.mass; }
  auto &mass() { return p.mass; }
#else
  constexpr auto &mass() const { return p.mass; }
#endif
#ifdef ESPRESSO_ROTATION
  auto const &rotation() const { return p.rotation; }
  auto &rotation() { return p.rotation; }
  bool can_rotate() const { return static_cast<bool>(p.rotation); }
  bool can_rotate_around(unsigned int const axis) const {
    assert(axis <= 2u);
    return detail::get_nth_bit(p.rotation, axis);
  }
  void set_can_rotate_around(unsigned int const axis, bool const rot_flag) {
    assert(axis <= 2u);
    if (rot_flag) {
      p.rotation |= static_cast<uint8_t>(1u << axis);
    } else {
      p.rotation &= static_cast<uint8_t>(~(1u << axis));
    }
  }
  void set_can_rotate_all_axes() { p.rotation = static_cast<uint8_t>(0b111u); }
  void set_cannot_rotate_all_axes() {
    p.rotation = static_cast<uint8_t>(0b000u);
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
  auto const &ext_torque() const { return p.ext_torque; }
  auto &ext_torque() { return p.ext_torque; }
#endif // ESPRESSO_EXTERNAL_FORCES
  auto calc_director() const {
    return Utils::convert_quaternion_to_director(quat());
  }
#else  // ESPRESSO_ROTATION
  auto can_rotate() const { return false; }
  auto can_rotate_around(unsigned int const) const { return false; }
#endif // ESPRESSO_ROTATION
#ifdef ESPRESSO_DIPOLES
  auto const &dipm() const { return p.dipm; }
  auto &dipm() { return p.dipm; }
  auto calc_dip() const { return calc_director() * dipm(); }
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  auto const &stoner_wohlfarth_is_enabled() const {
    return p.magnetodynamics.is_enabled;
  }
  auto &stoner_wohlfarth_is_enabled() { return p.magnetodynamics.is_enabled; }
  auto const &stoner_wohlfarth_phi_0() const { return p.magnetodynamics.phi0; }
  auto &stoner_wohlfarth_phi_0() { return p.magnetodynamics.phi0; }
  auto const &saturation_magnetization() const {
    return p.magnetodynamics.sat_mag;
  }
  auto &saturation_magnetization() { return p.magnetodynamics.sat_mag; }
  auto const &magnetic_anisotropy_field_inv() const {
    return p.magnetodynamics.ani_fld_inv;
  }
  auto &magnetic_anisotropy_field_inv() {
    return p.magnetodynamics.ani_fld_inv;
  }
  auto const &magnetic_anisotropy_energy() const {
    return p.magnetodynamics.ani_energy;
  }
  auto &magnetic_anisotropy_energy() { return p.magnetodynamics.ani_energy; }
  auto const &stoner_wohlfarth_tau0_inv() const {
    return p.magnetodynamics.tau0_inv;
  }
  auto &stoner_wohlfarth_tau0_inv() { return p.magnetodynamics.tau0_inv; }
  auto const &stoner_wohlfarth_dt_incr() const {
    return p.magnetodynamics.dt_incr;
  }
  auto &stoner_wohlfarth_dt_incr() { return p.magnetodynamics.dt_incr; }
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  auto const &dip_fld() const { return p.dip_fld; }
  auto &dip_fld() { return p.dip_fld; }
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  auto const &rinertia() const { return p.rinertia; }
  auto &rinertia() { return p.rinertia; }
#else
  constexpr auto &rinertia() const { return p.rinertia; }
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  auto const &q() const { return p.q; }
  auto &q() { return p.q; }
#else
  constexpr auto &q() const { return p.q; }
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  auto const &mu_E() const { return p.mu_E; }
  auto &mu_E() { return p.mu_E; }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES
  auto is_virtual() const {
    return (p.propagation & (PropagationMode::TRANS_VS_RELATIVE |
                             PropagationMode::TRANS_VS_CENTER_OF_MASS |
                             PropagationMode::ROT_VS_RELATIVE |
                             PropagationMode::ROT_VS_INDEPENDENT |
                             PropagationMode::TRANS_LB_TRACER)) != 0;
  }
#else
  constexpr auto is_virtual() const { return false; }
#endif // ESPRESSO_VIRTUAL_SITES
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  auto const &vs_relative() const { return p.vs_relative; }
  auto &vs_relative() { return p.vs_relative; }
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  auto const &gamma() const { return p.gamma; }
  auto &gamma() { return p.gamma; }
#ifdef ESPRESSO_ROTATION
  auto const &gamma_rot() const { return p.gamma_rot; }
  auto &gamma_rot() { return p.gamma_rot; }
#endif // ESPRESSO_ROTATION
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_EXTERNAL_FORCES
  auto const &fixed() const { return p.ext_flag; }
  auto &fixed() { return p.ext_flag; }
  bool has_fixed_coordinates() const { return static_cast<bool>(p.ext_flag); }
  bool is_fixed_along(unsigned int const axis) const {
    assert(axis <= 2u);
    return detail::get_nth_bit(p.ext_flag, axis);
  }
  void set_fixed_along(int const axis, bool const fixed_flag) {
    // set new flag
    if (fixed_flag) {
      p.ext_flag |= static_cast<uint8_t>(1u << axis);
    } else {
      p.ext_flag &= static_cast<uint8_t>(~(1u << axis));
    }
  }
  auto const &ext_force() const { return p.ext_force; }
  auto &ext_force() { return p.ext_force; }
#else  // ESPRESSO_EXTERNAL_FORCES
  constexpr bool has_fixed_coordinates() const { return false; }
  constexpr bool is_fixed_along(unsigned int const) const { return false; }
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_ENGINE
  auto const &swimming() const { return p.swim; }
  auto &swimming() { return p.swim; }
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
  auto const &rattle_params() const { return rattle; }
  auto &rattle_params() { return rattle; }
  auto const &rattle_correction() const { return rattle.correction; }
  auto &rattle_correction() { return rattle.correction; }
#endif

#ifdef ESPRESSO_EXCLUSIONS
  Utils::compact_vector<int> &exclusions() { return el; }
  Utils::compact_vector<int> const &exclusions() const { return el; }
  bool has_exclusion(int pid) const {
    return std::ranges::find(el, pid) != el.end();
  }
#endif

private:
  friend boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & p;
    ar & l;
    ar & bl;
#ifdef ESPRESSO_EXCLUSIONS
    ar & el;
#endif
    // Migration phases 2, 3 & 4: force/torque (phase 2), the STATE fields
    // (position, image box, quaternion, position-at-last-verlet-update,
    // position-at-last-time-step, Lees-Edwards offset and flag; phase 3), and
    // the MOMENTUM fields (velocity and angular velocity; phase 4) live in
    // ParticleStore columns, which are not carried by this (boost) serializer
    // used for cross-rank particle exchange. Ferry the values so they survive a
    // global resort that moves the particle to another rank (matching
    // pre-migration behavior). SAVE reads the current value (column if
    // attached, carrier otherwise); LOAD lands in the detached carrier, from
    // which ParticleStore::assign_row seeds the rebuilt row.
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
