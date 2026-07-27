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

#include "cell_system/Cell.hpp"

#include <span>
#include <unordered_set>

namespace GhostComm {

/**
 * @brief Classify each local cell as interior or boundary.
 *
 * A local cell is *boundary* iff at least one neighbor returned by
 * `cell->neighbors().all()` is a ghost cell (i.e. its `ParticleList *` is in
 * the set of ghost cells).  All other local cells are *interior*.
 *
 * The function first resets every local cell to interior (idempotent on
 * repeated calls / plan rebuild), then marks boundary cells.
 *
 * Call this right after `m_halo_plan = make_halo_plan()` in each
 * decomposition, where `local_cells()` and `ghost_cells()` are already
 * populated.
 */
inline void mark_boundary_cells(std::span<Cell *const> local_cells,
                                std::span<Cell *const> ghost_cells) {
  // Build a set of ghost ParticleList pointers for O(1) lookup.
  std::unordered_set<ParticleList const *> ghost_set;
  ghost_set.reserve(ghost_cells.size());
  for (Cell *c : ghost_cells) {
    ghost_set.insert(&c->particles());
  }

  // Reset all local cells to interior first (idempotent on rebuild).
  for (Cell *c : local_cells) {
    c->m_is_boundary = false;
  }

  // Mark a cell as boundary iff any of its neighbors is a ghost.
  for (Cell *c : local_cells) {
    for (Cell *n : c->neighbors().all()) {
      if (ghost_set.count(&n->particles())) {
        c->m_is_boundary = true;
        break;
      }
    }
  }
}

} // namespace GhostComm
