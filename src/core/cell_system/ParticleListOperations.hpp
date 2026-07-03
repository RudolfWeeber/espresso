/*
 * Copyright (C) 2017-2026 The ESPResSo project
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

#include "Particle.hpp"
#include "ParticleList.hpp"

#include <cstddef>
#include <utility>

/**
 * @brief Primitives for mutating the particle storage of cells.
 *
 * Every insertion, removal, or move of a @ref Particle in *cell storage*
 * (a @ref ParticleList owned by a @ref Cell, including ghost cells) MUST
 * go through these functions. They are the single hook point for
 * mirroring rows into the ParticleStore in later migration phases; see
 * docs/superpowers/specs/2026-07-03-array-based-particle-storage-design.md
 * section 3, phase 1.
 *
 * Plain communication buffers (send/receive @ref ParticleList instances
 * that do not belong to a cell) are exempt and use the @ref Utils::Bag
 * API directly.
 *
 * maintainer/CI/check_cell_storage_mutations.sh enforces this rule.
 */
namespace CellParticleStorage {

/**
 * @brief Insert a particle into a cell's particle storage.
 * May reallocate the storage, invalidating pointers and iterators into it.
 * @return Reference to the stored particle.
 */
inline Particle &insert_particle(ParticleList &storage, Particle &&particle) {
  return storage.insert(std::move(particle));
}

/**
 * @brief Move a particle out of a cell's particle storage.
 * Swap-with-back removal: element order is not preserved.
 * @return The extracted particle and the iterator past the removed element.
 */
inline std::pair<Particle, ParticleList::iterator>
extract_particle(ParticleList &storage, ParticleList::iterator position) {
  auto particle = std::move(*position);
  auto next = storage.erase(position);
  return {std::move(particle), next};
}

/**
 * @brief Erase a particle from a cell's particle storage, destroying it.
 * Swap-with-back removal: element order is not preserved.
 * @return Iterator past the removed element.
 */
inline ParticleList::iterator erase_particle(ParticleList &storage,
                                             ParticleList::iterator position) {
  return storage.erase(position);
}

/** @brief Remove all particles from a cell's particle storage. */
inline void clear_particles(ParticleList &storage) { storage.clear(); }

/**
 * @brief Resize ghost-cell storage to exactly @p count particles.
 * Newly created particles are default-constructed; all particles in the
 * storage are (re)marked as ghosts.
 */
inline void resize_ghost_storage(ParticleList &storage, std::size_t count) {
  storage.resize(count);
  for (auto &particle : storage) {
    particle.set_ghost(true);
  }
}

} // namespace CellParticleStorage
