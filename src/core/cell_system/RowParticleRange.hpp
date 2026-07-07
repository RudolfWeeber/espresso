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
#include "cell_system/CellRows.hpp"
#include "particle_store/ParticleStore.hpp"

#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range.hpp>

#include <cstddef>

/**
 * @brief Iterator over a @ref CellRows bag yielding @ref Particle views.
 *
 * Migration phase 7a groundwork (DORMANT). Walks the row indices stored in a
 * @ref CellRows bag and, for each, hands out a @ref Particle view over the
 * paired @ref ParticleStore row (built via @ref ParticleStore::make_view).
 * Its contract mirrors @ref ParticleIterator / @ref ParticleRange so that the
 * flip task (Task 4) can substitute it for the cell-Bag iterator with minimal
 * churn:
 *   - @c value_type is @ref Particle and @c reference is @ref Particle& (as in
 *     @ref ParticleIterator, whose @c dereference() returns @c Particle&);
 *   - it is a forward iterator (@c boost::forward_traversal_tag).
 *
 * LIFETIME CONTRACT (the load-bearing design point): @c operator* returns a
 * reference to a @ref Particle view CACHED INSIDE THE ITERATOR (the
 * @ref m_view member), refreshed on construction and on every @c increment().
 * That reference stays valid for as long as the iterator object lives and is
 * not incremented; incrementing or destroying the iterator invalidates it.
 * This lets a caller bind @c Particle &p = *it and keep using @c p across a
 * loop body (the pattern the bond handlers rely on -- see the phase 7a task-1
 * audit), which a by-value @c operator* returning a temporary could not
 * support. The referenced view itself aliases the store; it is invalidated by
 * a store rebuild just like any other view.
 */
class RowParticleIterator
    : public boost::iterator_facade<RowParticleIterator, Particle,
                                    boost::forward_traversal_tag, Particle &> {
  CellRows::const_iterator m_row;
  ParticleStore *m_store;
  /** Cached view over the current row; the target of @c operator*. Its address
   *  is stable for the lifetime of the iterator (per the class lifetime
   *  contract), refreshed to the current row on each dereference. */
  Particle m_view;

public:
  /** @brief Iterator at @p row over @p store. */
  RowParticleIterator(CellRows::const_iterator row, ParticleStore &store)
      : m_row(row), m_store(&store) {}

  /** @brief Past-the-end iterator (no store dereference happens). */
  explicit RowParticleIterator(CellRows::const_iterator end)
      : m_row(end), m_store(nullptr) {}

private:
  friend class boost::iterator_core_access;

  void increment() { ++m_row; }

  bool equal(RowParticleIterator const &rhs) const {
    return m_row == rhs.m_row;
  }

  // Lazily build the view on dereference so a valid past-the-end iterator (with
  // no store) is never dereferenced. The view is stored in m_view so the
  // returned reference outlives the expression, honouring the lifetime
  // contract documented on the class.
  Particle &dereference() const {
    // Refresh the mutable cache to the current row. `const` on iterator_facade
    // dereference is a boost requirement; the cache is logically transient.
    auto &self = const_cast<RowParticleIterator &>(*this);
    self.m_view = m_store->make_view(*m_row);
    return self.m_view;
  }
};

/**
 * @brief Range of @ref Particle views over a @ref CellRows bag + store.
 *
 * DORMANT phase 7a groundwork. Modelled on @ref ParticleRange (a
 * @c boost::iterator_range with a cached size). Iterating it yields the same
 * particle sequence, in the same order, as iterating the paired cell's
 * @ref ParticleList Bag on a synchronized store.
 */
class RowParticleRange : public boost::iterator_range<RowParticleIterator> {
  using base_type = boost::iterator_range<RowParticleIterator>;

public:
  RowParticleRange(CellRows const &rows, ParticleStore &store)
      : base_type(rows.empty() ? RowParticleIterator(rows.end())
                               : RowParticleIterator(rows.begin(), store),
                  RowParticleIterator(rows.end())),
        m_size(rows.size()) {}

  base_type::size_type size() const { return m_size; }

private:
  std::size_t m_size;
};
