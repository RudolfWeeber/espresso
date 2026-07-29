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

#include <functional>
#include <span>
#include <unordered_set>

namespace GhostComm {

/**
 * @brief Classify each local cell as interior or boundary.
 *
 * A local cell is *boundary* iff:
 *   (a) at least one neighbor returned by `cell->neighbors().all()` is a
 *       ghost cell (i.e. its `ParticleList *` is in the set of ghost cells),
 *       OR
 *   (b) the optional @p wrap_predicate returns `true` for the (cell, neighbor)
 *       pair — even when both cells are local.  This captures periodic
 *       wrap-around neighborships that exist on a single MPI rank where the
 *       local domain spans the whole box and no ghost cells are generated for
 *       that axis.  The predicate defaults to always-false (no wrapping).
 *
 * All other local cells are *interior*.
 *
 * The function first resets every local cell to interior (idempotent on
 * repeated calls / plan rebuild), then marks boundary cells.
 *
 * Call this right after `m_halo_plan = make_halo_plan()` in each
 * decomposition, where `local_cells()` and `ghost_cells()` are already
 * populated.
 *
 * @param local_cells  Local cell pointer span.
 * @param ghost_cells  Ghost cell pointer span.
 * @param wrap_predicate  Optional predicate `(Cell const *a, Cell const *b)
 *   -> bool` that returns `true` when the neighbor relation a→b crosses a
 *   periodic box boundary (both cells are local).  When omitted, no
 *   additional wrapping boundary is detected.
 */
inline void mark_boundary_cells(
    std::span<Cell *const> local_cells, std::span<Cell *const> ghost_cells,
    std::function<bool(Cell const *, Cell const *)> wrap_predicate = nullptr) {
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

  // Mark a cell as boundary iff any of its neighbors is a ghost (rule a) or
  // the wrap predicate fires for that neighbor pair (rule b).
  for (Cell *c : local_cells) {
    for (Cell *n : c->neighbors().all()) {
      if (ghost_set.count(&n->particles()) ||
          (wrap_predicate && wrap_predicate(c, n))) {
        c->m_is_boundary = true;
        break;
      }
    }
  }
}

} // namespace GhostComm
