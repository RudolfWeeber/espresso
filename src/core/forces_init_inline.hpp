/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

// Per-particle force initialization (external forces, Langevin friction and
// noise). Used by the separate pass init_forces_and_thermostat()
// (forces_init.cpp) and by the AoSoA commit sweep (short_range_verlet.cpp),
// which does the same work while it visits the particles anyway.

#include <config/config.hpp>

#include "Particle.hpp"
#include "PropagationMode.hpp"
#include "integrators/Propagation.hpp"
#include "rotation.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"
#include "thermostats/langevin_inline.hpp"

/** External particle forces */
inline ParticleForce external_force(Particle const &p) {
  ParticleForce f = {};

#ifdef ESPRESSO_EXTERNAL_FORCES
  f.f += p.ext_force();
#ifdef ESPRESSO_ROTATION
  f.torque += p.ext_torque();
#endif
#endif

#ifdef ESPRESSO_ENGINE
  // apply a swimming force in the direction of
  // the particle's orientation axis
  if (p.swimming().swimming and !p.swimming().is_engine_force_on_fluid) {
    f.f += p.swimming().f_swim * p.calc_director();
  }
#endif

  return f;
}

/**
 * @brief Initial force and torque of a local particle: external force and
 * torque, plus Langevin friction and noise if the particle is propagated
 * with the Langevin thermostat.
 */
class LocalForceInit {
  Propagation const &m_propagation;
  LangevinThermostat const *m_langevin; ///< nullptr if Langevin is inactive
  double m_kT;
  double m_time_step;

public:
  explicit LocalForceInit(System::System const &system)
      : m_propagation{*system.propagation}, m_langevin{nullptr},
        m_kT{system.thermostat->kT}, m_time_step{system.get_time_step()} {
    auto const &thermostat = *system.thermostat;
    if (thermostat.langevin and
        (m_propagation.used_propagations &
         (PropagationMode::TRANS_LANGEVIN | PropagationMode::ROT_LANGEVIN))) {
      m_langevin = thermostat.langevin.get();
    }
  }

  void operator()(Particle &p) const {
    // Initialize force with external forces
    p.force_and_torque() = external_force(p);

    // Apply Langevin noise if thermostat is active
    if (m_langevin) {
      auto const &langevin = *m_langevin;
      if (m_propagation.should_propagate_with(p,
                                              PropagationMode::TRANS_LANGEVIN))
        p.force() += friction_thermo_langevin(langevin, p, m_time_step, m_kT);
#ifdef ESPRESSO_ROTATION
      if (m_propagation.should_propagate_with(p, PropagationMode::ROT_LANGEVIN))
        p.torque() += convert_vector_body_to_space(
            p,
            friction_thermo_langevin_rotation(langevin, p, m_time_step, m_kT));
#endif
    }
  }
};
