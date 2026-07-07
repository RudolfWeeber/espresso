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

/** Owning result of @ref fetch_particles (phase 7e).
 *
 *  Phase 7e: get_local_particle returns by-value @ref Particle views (16-byte
 *  handles aliasing the store), so the reference range that observables consume
 *  can no longer point into a persistent view pool -- it must point into an
 *  owned buffer that lives as long as the range. This struct owns that buffer
 *  (@ref owned) and exposes the @ref ParticleReferenceRange (@ref refs) that
 *  references into it. Callers keep the returned object alive for the duration
 *  of the observable evaluation and pass @c .refs to @c evaluate. */
struct FetchedParticles {
  std::vector<Particle> owned;
  Observables::ParticleReferenceRange refs;
};

/** Fetch a group of particles.
 *
 *  @param ids particle identifiers
 *  @return owned particle views (with positions in the current box) plus a
 *          reference range into them (phase 7e).
 */
inline FetchedParticles fetch_particles(std::vector<int> const &ids) {
  auto const &system = System::get_system();
  auto &cell_structure = *system.cell_structure;
  FetchedParticles result;
  // Reserve so `owned` never reallocates while we take references into it
  // (ids.size() is the upper bound; ghosts/absent ids are filtered out).
  result.owned.reserve(ids.size());
  // Resolve each requested id through get_local_particle, which returns a
  // by-value view over the store row -- valid for the observable's lifetime (no
  // rebuild during evaluation). Skip ids that are not a local (owned) particle
  // here (absent or ghost), preserving the pre-flip filter (only local,
  // non-ghost particles).
  for (auto const id : ids) {
    auto const p = cell_structure.get_local_particle(id);
    if (p and not p->is_ghost()) {
      result.owned.emplace_back(*p);
    }
  }
  for (auto const &p : result.owned) {
    result.refs.emplace_back(std::cref(p));
  }
  return result;
}
#endif
