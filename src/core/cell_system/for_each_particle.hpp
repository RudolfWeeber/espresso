
/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "cell_system/CellStructure.hpp"

#include <Kokkos_Core.hpp>

#include <span>

template <typename Callable>
inline void
CellStructure::parallel_for_each_particle_impl(std::span<Cell *const> cells,
                                               Callable &f) const {
  auto &store = const_cast<ParticleStore &>(m_particle_store);
  if (cells.size() > 1) {
    Kokkos::parallel_for( // loop over cells
        "for_each_local_particle", cells.size(), [&](auto cell_idx) {
          // One row-range iterator per cell (per thread) reuses its own cached
          // view across the cell's particles (phase 7a lifetime contract).
          for (auto &p : cells[cell_idx]->particles())
            f(p);
        });
  } else if (cells.size() == 1) {
    // Single-cell parallel-over-particles: each index gets its OWN view (a
    // shared cached-view iterator would not be thread-safe here).
    auto const &rows = cells.front()->rows();
    Kokkos::parallel_for( // loop over particles
        "for_each_local_particle", rows.size(), [&](auto part_idx) {
          auto view = store.make_view(rows.begin()[part_idx]);
          f(view);
        });
  }
}
