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
                         VerletCriterion<> const &criterion, double pair_cutoff,
                         int integ_switch) {
  update_cabana_state(cell_structure, criterion, pair_cutoff, integ_switch);
}
