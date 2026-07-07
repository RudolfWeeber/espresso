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

#pragma once

#include <utils/Bag.hpp>

/**
 * @brief Per-cell container of @ref ParticleStore row indices.
 *
 * Each @ref Cell carries a @ref CellRows bag as its authoritative content
 * (phase 7a flip, Task 4): every particle in the cell is represented by
 * exactly one row index in this bag (one @c int per particle, in cell-traversal
 * order). @ref Cell::particles() constructs a @ref RowParticleRange over this
 * bag and the cell's associated @ref ParticleStore, yielding live @ref Particle
 * views; it is the primary way to iterate a cell's particles.
 *
 * A @ref Utils::Bag is used (rather than a plain @c std::vector) so that
 * @ref CellParticleStorage::extract_row can remove a row via constant-time
 * swap-with-back erase, renumbering the cell's row sequence the same way the
 * old @c Bag<Particle> removed particles.  @c int is trivially
 * swappable/movable so the Bag has no per-element overhead.
 *
 * Rebuilt from scratch during every
 * @ref CellStructure::ensure_particle_store_synchronized; the row-bag churn
 * during extract is bitwise-equivalent to the pre-flip Bag<Particle> churn.
 * The collapse to a contiguous @c (offset,count) range is deferred to phase 7c.
 */
using CellRows = Utils::Bag<int>;
