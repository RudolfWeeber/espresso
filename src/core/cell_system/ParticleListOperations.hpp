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
#include "cell_system/Cell.hpp"

#include <cassert>
#include <cstddef>
#include <utility>

/**
 * @brief Primitives for mutating the particle storage of cells (phase 7a; 7b
 * row-ref staging).
 *
 * Since the cell flip, a @ref Cell no longer owns @c Particle objects: it owns
 * @ref ParticleStore ROW INDICES (@ref Cell::rows) for its committed particles
 * plus a small staging buffer (@ref Cell::staged) of @ref StagedParticle row
 * references not yet committed. Every insertion, removal, or move of a particle
 * in *cell storage* (including ghost cells) MUST go through these functions.
 * They are the single choke point that keeps the store's row bookkeeping
 * consistent; the actual column data is (re)materialized by the next store
 * rebuild (@ref CellStructure::ensure_particle_store_synchronized), so all of
 * these mutations only stage/unstage rows and mark the store dirty -- the
 * caller is responsible for the dirty mark (via
 * @ref CellStructure::mark_particle_store_dirty) and for triggering the rebuild
 * before the changes become visible through @ref Cell::particles().
 *
 * Phase 7b (Task 4): the migration envelope is dead, so a cell can no longer
 * stage a detached, data-carrying @c Particle. It stages a @ref StagedParticle
 * (a reference to a SOURCE-store row, or a fresh-default marker); the rebuild
 * COPIES that source row into a fresh committed row (@ref
 * ParticleStore::copy_row) or seeds it to defaults (@ref
 * ParticleStore::seed_default_row). Callers that need to move a live row (local
 * mis-cell move, migration receive) first copy it into a source store
 * (@ref CellStructure::stage_row) and then stage a reference to that row here.
 *
 * maintainer/CI/check_cell_storage_mutations.sh enforces this rule.
 */
namespace CellParticleStorage {

/**
 * @brief Stage a source-store row for insertion into a cell (phase 7b).
 *
 * Appends a @ref StagedParticle referencing @p source_row of @p source_store to
 * the cell's staging buffer; the next store rebuild copies that row into a
 * fresh committed row (@ref ParticleStore::copy_row). The source row (and its
 * store) must remain valid until the rebuild runs. The caller must mark the
 * store dirty.
 */
inline void insert_staged_row(Cell &cell, ParticleStore &source_store,
                              int const source_row) {
  cell.staged().push_back(StagedParticle{&source_store, source_row});
}

/**
 * @brief Drop the committed row at bag position @p index from a cell (phase
 * 7b).
 *
 * Removes the cell's row-index bag entry via a constant-time swap-with-back
 * erase (so the row-bag churn is bitwise-equivalent to the pre-flip
 * @c Bag<Particle> erase: element order not preserved -- surviving rows are
 * renumbered on the next rebuild anyway). The store row itself is dropped by
 * the next rebuild; the caller must mark the store dirty. @p index refers to a
 * position in @ref Cell::rows (0-based).
 */
inline void drop_row(Cell &cell, std::size_t index) {
  auto &rows = cell.rows();
  assert(index < rows.size());
  rows.erase(rows.begin() + static_cast<std::ptrdiff_t>(index));
}

/** @brief Remove all committed rows and staged particles from a cell. */
inline void clear_particles(Cell &cell) {
  cell.rows().clear();
  cell.staged().clear();
}

/**
 * @brief Resize a ghost cell to exactly @p count particles (phase 7a; 7b).
 *
 * Drops the cell's committed rows and staging buffer and stages @p count
 * fresh-default ghost particles (a @ref StagedParticle with a null source
 * store); the next store rebuild seeds each as a fresh default ghost row
 * (@ref ParticleStore::seed_default_row). Ghost-ness is STRUCTURAL: these rows
 * land in the store's ghost suffix @c [n_local, n_total) by construction, so no
 * flag is carried. The caller must mark the store dirty. (Pre-flip this resized
 * the cell's @ref ParticleList in place; the staging + rebuild reproduces the
 * same "count default ghosts" outcome for the row-based cell.)
 */
inline void resize_ghost_storage(Cell &cell, std::size_t count) {
  cell.rows().clear();
  auto &staged = cell.staged();
  staged.clear();
  staged.resize(count); // default StagedParticle: source_store == nullptr
}

} // namespace CellParticleStorage
