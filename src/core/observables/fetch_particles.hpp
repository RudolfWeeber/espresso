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

#ifndef FETCH_PARTICLES_HPP
#define FETCH_PARTICLES_HPP

#include "PidObservable.hpp"
#include "cell_system/CellStructure.hpp"
#include "system/System.hpp"

#include <algorithm>
#include <functional>
#include <iterator>
#include <set>
#include <vector>

/** Fetch a group of particles.
 *
 *  @param ids particle identifiers
 *  @return array of particle copies, with positions in the current box.
 */
inline auto fetch_particles(std::vector<int> const &ids) {
  auto const &system = System::get_system();
  auto &cell_structure = *system.cell_structure;
  Observables::ParticleReferenceRange local_particle_refs;
  // Phase 7a: local_particles() hands out transient cached VIEWS, so iterating
  // + copying references (the pre-flip std::copy_if approach) would store
  // references to a view that is overwritten on the next increment (the
  // cached-view multipass hazard). Instead resolve each requested id through
  // get_local_particle, which returns a pointer into the STABLE view pool --
  // valid for the observable's lifetime (no rebuild during evaluation). Skip
  // ids that are not a local (owned) particle here (nullptr or ghost),
  // preserving the pre-flip filter (only local, non-ghost particles).
  for (auto const id : ids) {
    auto const *p = cell_structure.get_local_particle(id);
    if (p != nullptr and not p->is_ghost()) {
      local_particle_refs.emplace_back(std::cref(*p));
    }
  }
  return local_particle_refs;
}
#endif
