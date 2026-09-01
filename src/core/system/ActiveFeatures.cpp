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

#include <config/config.hpp>

#include "ActiveFeatures.hpp"

#include "Particle.hpp"
#include "PropagationMode.hpp"
#include "cell_system/CellStructure.hpp"
#include "collision_detection/CollisionDetection.hpp"
#include "communication.hpp"
#include "integrators/Propagation.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_reduction.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <functional>

namespace System {

unsigned ActiveFeatures::particle_bits(Particle const &p) {
  auto features = 0u;
  if (p.is_virtual()) {
    features |= IS_VIRTUAL;
  }
  if (p.can_rotate()) {
    features |= CAN_ROTATE;
  }
  if (p.has_fixed_coordinates()) {
    features |= IS_FIXED;
  }
  if (p.q() != 0.) {
    features |= HAS_CHARGE;
  }
#ifdef ESPRESSO_EXCLUSIONS
  if (not p.exclusions().empty()) {
    features |= HAS_EXCLUSIONS;
  }
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  if (p.stoner_wohlfarth_is_enabled()) {
    features |= HAS_STONER_WOHLFARTH;
  }
#endif
#ifdef ESPRESSO_DIPOLES
  if (p.dipm() != 0.) {
    features |= HAS_DIPOLE_MOMENT;
  }
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  if (p.ext_force() != Utils::Vector3d{0., 0., 0.}) {
    features |= HAS_EXT_FORCE;
  }
#ifdef ESPRESSO_ROTATION
  if (p.ext_torque() != Utils::Vector3d{0., 0., 0.}) {
    features |= HAS_EXT_TORQUE;
  }
#endif
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_ENGINE
  if (p.swimming().swimming) {
    features |= IS_SWIMMER;
  }
#endif
  return features;
}

void ActiveFeatures::update() {
  auto &system = get_system();
  auto &propagation = *system.propagation;
  propagation.update_default_propagation(system.thermostat->thermo_switch);
  struct SweepResult {
    int propagations = PropagationMode::NONE;
    unsigned features = 0u;
  };
  auto const local = reduce_over_local_particles<SweepResult>(
      *system.cell_structure,
      [](SweepResult &acc, Particle const &p) {
        acc.propagations |= p.propagation();
        acc.features |= particle_bits(p);
      },
      [](SweepResult &acc, SweepResult const &other) {
        acc.propagations |= other.propagations;
        acc.features |= other.features;
      });
  int const local_masks[2] = {local.propagations,
                              static_cast<int>(local.features)};
  int global_masks[2];
  boost::mpi::all_reduce(::comm_cart, local_masks, 2, global_masks,
                         std::bit_or<int>());
  auto used_propagations = global_masks[0];
  if (used_propagations & PropagationMode::SYSTEM_DEFAULT) {
    used_propagations |= propagation.default_propagation;
  }
  propagation.used_propagations = used_propagations;
  m_particle_features = static_cast<unsigned>(global_masks[1]);
  propagation.recalc_active_features = false;
}

void ActiveFeatures::update_if_needed() {
  if (get_system().propagation->recalc_active_features) {
    update();
  }
}

bool ActiveFeatures::has_gay_berne() const {
  return get_system().nonbonded_ias->pair_potential_active(
      PairPotential::GayBerne);
}

#ifdef ESPRESSO_ROTATION
/**
 * Return true when any active physics requires orientation of ghost
 * particles.  Used by both @c System::get_global_ghost_flags (QUAT push) and
 * @c System::get_force_reduce_ghost_flags (TORQUE reduce).
 *
 * Conservative: include a reader when in doubt.  Non-rotating plain LJ
 * must yield false so that both bits stay OFF for that common case; the
 * rotational-propagation arm therefore additionally requires a particle
 * that can actually rotate, since default_propagation carries ROT_* bits
 * for every SYSTEM_DEFAULT particle regardless of its rotation flags.
 *
 * Readers of ghost particle quat / director / calc_dip:
 * - short_range_cabana commit_particle: quat -> director for Gay-Berne /
 *   Dipoles pair kernels (ghost particles in AoSoA)
 * - vs_relative_update_particles: p_ref.quat() where p_ref may be a ghost
 * - LB particle coupling (swimmer): p.calc_director() on ghost particles
 *   when ESPRESSO_ENGINE and LB are both active
 * - HomogeneousMagneticField constraint: p.calc_dip() (local particles only,
 *   does NOT read ghosts; included conservatively via dipole check)
 * - dipolar solvers (dp3m, dds, dlc, scafacos): p.calc_dip() on
 *   unique_particles which includes ghost particles that have no local
 *   counterpart (see @ref CellStructure::set_index_map); covered by
 *   the dipolar-solver-set condition
 * - ShapeBasedConstraint / calc_non_central_force: reads p.quat() on LOCAL
 *   particles only (Constraints iterate local_particles); safe without the
 *   bit
 * - Stoner-Wohlfarth integrate_magnetodynamics: calls p_ref->calc_director()
 *   where p_ref comes from get_reference_particle / get_local_particle, which
 *   CAN return a ghost; covered by the HAS_STONER_WOHLFARTH particle bit,
 *   not by a thermostat proxy
 * - rotational integrators: operate on local particles only (OK)
 * - ICC blocking reduce: resets force_and_torque on local particles and may
 *   carry a zero-valued TORQUE payload on the wire when orientation physics is
 *   concurrently active — harmless (bytes only, value is zero)
 * - calculate_vs_relate_to_params (virtual_sites.cpp, reached from
 *   collision_detection/utils.hpp place_vs_and_relate_to_particle): reads
 *   p_relate_to.quat() where p_relate_to can be a ghost when collision
 *   detection creates a virtual site mid-step; covered by the
 *   collision-detection-active condition below, not by a particle bit
 */
bool ActiveFeatures::orientation_ghosts_needed() const {
  auto const &system = get_system();

  // Rotational propagation active on a particle that can actually rotate:
  // ghost quat needed for any kernel that uses the director of a ghost
  // particle. The can-rotate check matters because default_propagation
  // carries ROT_* bits for every SYSTEM_DEFAULT particle, while all
  // rotational integrators skip particles with rotation flags 0b000.
  if ((system.propagation->used_propagations &
       (PropagationMode::ROT_EULER | PropagationMode::ROT_LANGEVIN |
        PropagationMode::ROT_BROWNIAN | PropagationMode::ROT_STOKESIAN |
        PropagationMode::ROT_VS_RELATIVE |
        PropagationMode::ROT_VS_INDEPENDENT)) and
      particles_can_rotate()) {
    return true;
  }

  // Virtual sites relative uses p_ref.quat() where p_ref may be a ghost.
  if (system.propagation->used_propagations &
      (PropagationMode::TRANS_VS_RELATIVE | PropagationMode::ROT_VS_RELATIVE |
       PropagationMode::ROT_VS_INDEPENDENT)) {
    return true;
  }

#ifdef ESPRESSO_DIPOLES
  // Dipolar solver set and at least one particle carries a dipole moment:
  // short_range_cabana reads quat -> director for dipole pair kernels on
  // ghosts, and the dipolar solvers read p.calc_dip() on unique_particles.
  if (system.dipoles.impl and system.dipoles.impl->solver and
      particles_have_dipole_moment()) {
    return true;
  }
#endif

#ifdef ESPRESSO_GAY_BERNE
  // Gay-Berne anisotropic nonbonded interaction configured: short_range_cabana
  // commits ghost quat -> director for the Cabana pair kernel.
  if (has_gay_berne()) {
    return true;
  }
#endif

#if defined(ESPRESSO_COLLISION_DETECTION) and                                  \
    defined(ESPRESSO_VIRTUAL_SITES_RELATIVE)
  // Collision detection can create virtual sites mid-step;
  // calculate_vs_relate_to_params reads the quaternion of the reference
  // particle, which can be a ghost (collision_detection/utils.hpp).
  if (not system.collision_detection->is_off()) {
    return true;
  }
#endif

#ifdef ESPRESSO_ENGINE
  // LB active with a swimmer present: the coupling loop includes ghost
  // particles and calls p.calc_director() for swimmer_force_on_fluid
  // particles.
  if (system.lb.is_solver_set() and particles_are_swimmers()) {
    return true;
  }
#endif

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  // Stoner-Wohlfarth integrate_magnetodynamics calls
  // get_reference_particle / get_local_particle which can return a ghost,
  // then reads p_ref->calc_director() on it.
  if (particles_have_stoner_wohlfarth()) {
    return true;
  }
#endif

  return false;
}
#endif // ESPRESSO_ROTATION

} // namespace System
