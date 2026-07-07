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

#include "Particle.hpp"
#include "ParticleList.hpp"
#include "cell_system/CellRows.hpp"

#include <boost/range/iterator_range.hpp>

#include <algorithm>
#include <span>
#include <utility>
#include <vector>

template <class CellRef> class Neighbors {
  using storage_type = std::vector<CellRef>;

public:
  using value_type = storage_type::value_type;
  using iterator = storage_type::iterator;
  using const_iterator = storage_type::const_iterator;
  using cell_range = boost::iterator_range<iterator>;

private:
  void copy(const Neighbors &rhs) {
    m_neighbors = rhs.m_neighbors;
    m_red_black_divider =
        m_neighbors.begin() +
        std::distance(rhs.m_neighbors.begin(),
                      const_iterator(rhs.m_red_black_divider));
  }

public:
  Neighbors() = default;
  Neighbors(const Neighbors &rhs) { copy(rhs); }
  Neighbors &operator=(const Neighbors &rhs) {
    copy(rhs);
    return *this;
  }

  Neighbors(std::span<const CellRef> red_neighbors,
            std::span<const CellRef> black_neighbors) {
    m_neighbors.resize(red_neighbors.size() + black_neighbors.size());
    auto const res = std::ranges::copy(red_neighbors, m_neighbors.begin());
    m_red_black_divider = res.out;
    std::ranges::copy(black_neighbors, m_red_black_divider);
  }

  /**
   * @brief All neighbors.
   */
  cell_range all() { return {m_neighbors.begin(), m_neighbors.end()}; }
  /**
   * @brief Red partition of neighbors.
   *
   * An partition of the neighbors so that iterating over all
   * neighbors_red of all cells visits every pair exactly once.
   * Complement of neighbors_black.
   */
  cell_range red() { return {m_neighbors.begin(), m_red_black_divider}; }
  /**
   * @brief Black partition of neighbors.
   *
   * An partition of the neighbors so that iterating over all
   * neighbors_black of all cells visits every pair exactly once.
   * Complement of neighbors_red.
   */
  cell_range black() { return {m_red_black_divider, m_neighbors.end()}; }

private:
  /** Container with all the neighbors.
      Red neighbors are first, black second. */
  storage_type m_neighbors;
  /** Iterator pointing to the first black neighbor
      in the container. */
  iterator m_red_black_divider;
};

class Cell {
  using neighbors_type = Neighbors<Cell *>;

  ParticleList m_particles;

  /** Store row indices of the particles in @ref m_particles (phase 7a).
   *  Parallel to @ref m_particles: one entry per particle, in the same
   *  iteration order, refilled during every store rebuild. DORMANT -- read by
   *  no production code yet (Task 4 will flip @ref particles() to hand out
   *  views over these rows). */
  CellRows m_rows;

public:
  /** Particles */
  auto &particles() { return m_particles; }
  auto const &particles() const { return m_particles; }

  /** @brief Store row indices of this cell's particles (phase 7a, dormant). */
  auto &rows() { return m_rows; }
  auto const &rows() const { return m_rows; }

  neighbors_type m_neighbors;

  /** Interaction pairs as @ref ParticleStore ROW indices (phase 7a). Cells no
   *  longer own stable @c Particle addresses, so pairs are recorded as store
   *  rows rather than pointers, mirroring @ref CellStructure::m_verlet_list.
   *  NOTE: this member has no reader or writer anywhere in the code base (the
   *  live Verlet list is @ref CellStructure::m_verlet_list); the type is kept
   *  in sync with the pointer-lifetime hardening for consistency and is a
   *  candidate for removal in a later phase. */
  std::vector<std::pair<int, int>> m_verlet_list;

  /**
   * @brief All neighbors of the cell.
   */
  neighbors_type &neighbors() { return m_neighbors; }
};
