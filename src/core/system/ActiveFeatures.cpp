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
#include "communication.hpp"
#include "integrators/Propagation.hpp"
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

} // namespace System
