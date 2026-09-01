/*
 * Copyright (C) 2026 The ESPResSo project
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

#include "system/Leaf.hpp"

struct Particle;

namespace System {

/**
 * @brief Centralized runtime feature state.
 *
 * Answers "is feature X in use right now?".
 * Particle-derived predicates read a cached bitmask that is filled by a
 * single sweep over local particles and reduced across all MPI ranks with
 * a bitwise OR. System-derived predicates read live system state through
 * @ref Leaf::get_system and store no copies.
 *
 * Queries return the state as of the last update. They never trigger the
 * reduction lazily: the reduction is collective, and query sites are not
 * guaranteed to be collective.
 *
 * The cached state is invalidated via Propagation::recalc_active_features
 * (set on particle changes, integrator changes, and thermostat changes)
 * and refreshed at collective update points: the start of integrate(),
 * on_observable_calc(), update_dependent_particles(), and after
 * collision-detection topology changes.
 */
class ActiveFeatures : public Leaf<ActiveFeatures> {
public:
  /** @brief Recompute the particle-derived state if invalidated.
   *  Collective call: all ranks must enter together. */
  void update_if_needed();
  /** @brief Recompute the particle-derived state unconditionally.
   *  Collective call: all ranks must enter together. */
  void update();

  bool particles_are_virtual() const {
    return (m_particle_features & IS_VIRTUAL) != 0u;
  }
  bool particles_can_rotate() const {
    return (m_particle_features & CAN_ROTATE) != 0u;
  }
  bool particles_are_fixed() const {
    return (m_particle_features & IS_FIXED) != 0u;
  }
  bool particles_have_charge() const {
    return (m_particle_features & HAS_CHARGE) != 0u;
  }
#ifdef ESPRESSO_EXCLUSIONS
  bool particles_have_exclusions() const {
    return (m_particle_features & HAS_EXCLUSIONS) != 0u;
  }
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  bool particles_have_stoner_wohlfarth() const {
    return (m_particle_features & HAS_STONER_WOHLFARTH) != 0u;
  }
#endif
#ifdef ESPRESSO_DIPOLES
  bool particles_have_dipole_moment() const {
    return (m_particle_features & HAS_DIPOLE_MOMENT) != 0u;
  }
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  bool particles_have_ext_force() const {
    return (m_particle_features & HAS_EXT_FORCE) != 0u;
  }
#ifdef ESPRESSO_ROTATION
  bool particles_have_ext_torque() const {
    return (m_particle_features & HAS_EXT_TORQUE) != 0u;
  }
#endif
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_ENGINE
  bool particles_are_swimmers() const {
    return (m_particle_features & IS_SWIMMER) != 0u;
  }
#endif

private:
  enum ParticleFeature : unsigned {
    IS_VIRTUAL = 1u << 0u,
    CAN_ROTATE = 1u << 1u,
    IS_FIXED = 1u << 2u,
    HAS_CHARGE = 1u << 3u,
#ifdef ESPRESSO_EXCLUSIONS
    HAS_EXCLUSIONS = 1u << 4u,
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
    HAS_STONER_WOHLFARTH = 1u << 5u,
#endif
#ifdef ESPRESSO_DIPOLES
    HAS_DIPOLE_MOMENT = 1u << 6u,
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    HAS_EXT_FORCE = 1u << 7u,
#ifdef ESPRESSO_ROTATION
    HAS_EXT_TORQUE = 1u << 8u,
#endif
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_ENGINE
    IS_SWIMMER = 1u << 9u,
#endif
  };

  static unsigned particle_bits(Particle const &p);

  /** Bitwise OR of @ref ParticleFeature over all particles on all ranks. */
  unsigned m_particle_features = 0u;
};

} // namespace System
