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

// Force initialization and Langevin thermostat live in their own translation
// unit. Keeping the init tree out of forces.cpp leaves the whole-TU
// inline-growth budget of forces.cpp for the hot pair kernel
// (forces_cabana.hpp) and the Verlet-list build.

#include "forces_init.hpp"
#include "forces_init_inline.hpp"

#include <config/config.hpp>

#include "cell_system/CellStructure.hpp"
#include "cell_system/for_each_particle.hpp"
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
#include "magnetostatics/dipoles.hpp"
#endif
#include "system/System.hpp"

#ifdef ESPRESSO_CALIPER
#include "caliper_utils.hpp"

#include <caliper/cali.h>
#endif

/** Combined force initialization and Langevin noise application */
void init_forces_and_thermostat(System::System const &system) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif

  // Single pass over all local particles
  system.cell_structure->for_each_local_particle(LocalForceInit{system});
  // The local force/torque scatter buffers were already zeroed this force
  // call by prepare_verlet_list_cabana() (via update_verlet_state).

  init_ghost_forces(system);
}

void init_ghost_forces(System::System const &system) {
  auto &cell_structure = *system.cell_structure;
  cell_structure.ghosts_reset_forces();
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (system.dipoles.impl->solver.has_value()) {
    cell_structure.ghosts_reset_dipole_fields();
  }
#endif
}
