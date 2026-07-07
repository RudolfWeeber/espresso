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
 * Migration phase 7a: each @ref Cell carries a @ref CellRows bag in parallel to
 * its @ref ParticleList, recording the store row assigned to every particle in
 * the cell (one @c int per particle, in @ref ParticleList iteration order).
 * A @ref Utils::Bag is used (rather than a plain @c std::vector) so that the
 * eventual flip (Task 4) can reuse the Bag's constant-time swap-with-back
 * erase to renumber a cell's rows the same way particles are removed today;
 * @c int is trivially swappable/movable so the Bag has no per-element cost.
 *
 * The bag is DORMANT in phase 7a: it is refilled from scratch during every
 * @ref CellStructure::ensure_particle_store_synchronized rebuild and read by
 * no production code yet (only cross-checked against the @ref ParticleList
 * under @c ESPRESSO_ADDITIONAL_CHECKS). The collapse to a contiguous
 * @c (offset,count) range happens in phase 7c.
 */
using CellRows = Utils::Bag<int>;
