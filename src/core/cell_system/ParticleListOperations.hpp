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
 * @brief Primitives for mutating the particle storage of cells (phase 7a).
 *
 * Since the cell flip, a @ref Cell no longer owns @c Particle objects: it owns
 * @ref ParticleStore ROW INDICES (@ref Cell::rows) for its committed particles
 * plus a small staging buffer (@ref Cell::staged) of not-yet-committed detached
 * particles. Every insertion, removal, or move of a particle in *cell storage*
 * (including ghost cells) MUST go through these functions. They are the single
 * choke point that keeps the store's row bookkeeping consistent; the actual
 * column data is (re)materialized by the next store rebuild
 * (@ref CellStructure::ensure_particle_store_synchronized), so all of these
 * mutations only stage/unstage rows and mark the store dirty -- the caller is
 * responsible for the dirty mark (via @ref CellStructure::mark_particle_store_
 * dirty) and for triggering the rebuild before the changes become visible
 * through @ref Cell::particles().
 *
 * Plain communication buffers (send/receive @ref ParticleList instances that
 * do not belong to a cell) are exempt and use the @ref Utils::Bag API directly.
 *
 * Exception (decomposition swap/teardown): @ref
 * CellStructure::set_particle_decomposition re-inserts particles into the new
 * decomposition through the routed @ref CellStructure::add_particle, then
 * destroys the old @ref ParticleDecomposition wholesale. The old cells are torn
 * down by their destructors, which never touch these primitives; the whole
 * store is rebuilt (mark-dirty), not per-row.
 *
 * maintainer/CI/check_cell_storage_mutations.sh enforces this rule.
 */
namespace CellParticleStorage {

/**
 * @brief Stage a particle for insertion into a cell.
 *
 * The incoming detached, carrier-laden @p particle is appended to the cell's
 * staging buffer; it becomes a committed row (with a store row and live
 * columns seeded from its carriers) only at the next store rebuild. The caller
 * must mark the store dirty. Returns a reference to the staged particle (valid
 * until the buffer reallocates or the rebuild consumes it).
 */
inline Particle &insert_particle(Cell &cell, Particle &&particle) {
  cell.staged().push_back(std::move(particle));
  return cell.staged().back();
}

/**
 * @brief Snapshot the committed row at bag position @p index out of a cell.
 *
 * Returns a detached, carrier-laden @c Particle (via @ref
 * ParticleStore::snapshot_row) holding the row's current column/sidecar values,
 * and removes that row index from the cell's row bag (constant-time
 * swap-with-back, so row-bag order is not preserved -- the surviving rows are
 * renumbered on the next rebuild anyway). The store row itself is dropped by
 * the next rebuild; the caller must mark the store dirty. @p index refers to a
 * position in @ref Cell::rows (0-based).
 */
inline Particle extract_row(Cell &cell, std::size_t index) {
  auto &rows = cell.rows();
  assert(index < rows.size());
  auto const row = rows.begin()[index];
  auto snapshot = cell.store().snapshot_row(row);
  // Constant-time swap-with-back removal of the row index (mirrors the pre-flip
  // Bag<Particle> erase: element order not preserved).
  rows.erase(rows.begin() + static_cast<std::ptrdiff_t>(index));
  return snapshot;
}

/** @brief Remove all committed rows and staged particles from a cell. */
inline void clear_particles(Cell &cell) {
  cell.rows().clear();
  cell.staged().clear();
}

/**
 * @brief Resize a ghost cell to exactly @p count particles (phase 7a).
 *
 * Drops the cell's committed rows and staging buffer and stages @p count
 * default-constructed ghost particles; the next store rebuild commits them as
 * fresh ghost rows. All staged particles are marked as ghosts. The caller must
 * mark the store dirty. (Pre-flip this resized the cell's @ref ParticleList in
 * place; the staging + rebuild reproduces the same "count default ghosts"
 * outcome for the row-based cell.)
 */
inline void resize_ghost_storage(Cell &cell, std::size_t count) {
  cell.rows().clear();
  auto &staged = cell.staged();
  staged.clear();
  staged.resize(count);
  for (auto &particle : staged) {
    particle.set_ghost(true);
  }
}

} // namespace CellParticleStorage
