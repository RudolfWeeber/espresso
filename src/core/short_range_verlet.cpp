/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

// Sole instantiation site of the update_cabana_state AoSoA-commit /
// Verlet-list-build template for VerletCriterion<>. forces.cpp, energy.cpp and
// pressure.cpp used to each instantiate this giant (identical each time); they
// now call the non-template wrapper below so the build-lambda instantiations
// exist ONLY here.

#include "short_range_verlet.hpp"

#include <config/config.hpp>

#include "cell_system/CellStructure.hpp"
#include "nonbonded_interactions/VerletCriterion.hpp"

#include "short_range_cabana.hpp"

void update_verlet_state(CellStructure &cell_structure,
                         System::System const &system, double coulomb_cut,
                         double dipolar_cut, double collision_cut,
                         double pair_cutoff, int integ_switch) {
  auto const skin = cell_structure.get_verlet_skin();
  // When no electrostatics/dipolar/collision cutoff is active, build the
  // criterion variant that compiles those dead per-candidate branches out of
  // the build loop; otherwise the full one. pair_cutoff is the interaction
  // range, which is also the criterion's maximum cutoff. Both
  // update_cabana_state instantiations live in this single TU.
  bool const short_range_only = coulomb_cut == inactive_cutoff and
                                dipolar_cut == inactive_cutoff and
                                collision_cut == inactive_cutoff;
  if (short_range_only) {
    VerletCriterion<GetNonbondedCutoff, true> const criterion{
        system, skin, pair_cutoff, coulomb_cut, dipolar_cut, collision_cut};
    update_cabana_state(cell_structure, criterion, pair_cutoff, integ_switch);
  } else {
    VerletCriterion<GetNonbondedCutoff, false> const criterion{
        system, skin, pair_cutoff, coulomb_cut, dipolar_cut, collision_cut};
    update_cabana_state(cell_structure, criterion, pair_cutoff, integ_switch);
  }
}
