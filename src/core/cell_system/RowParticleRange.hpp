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
#include <optional>

/**
 * @brief Iterator over a @ref CellRows bag yielding @ref Particle views.
 *
 * The live iterator backing @ref Cell::particles() since the phase-7a flip
 * (Task 4). Walks the row indices stored in a @ref CellRows bag and, for each,
 * hands out a @ref Particle view over the paired @ref ParticleStore row (a
 * cached view rebound via @ref Particle::attach_to_store; see the PERFORMANCE
 * note below).  Its contract mirrors @ref ParticleIterator / @ref
 * ParticleRange:
 *   - @c value_type is @ref Particle and @c reference is @ref Particle& (as in
 *     @ref ParticleIterator, whose @c dereference() returns @c Particle&);
 *   - it is a forward iterator (@c boost::forward_traversal_tag).
 *
 * LIFETIME CONTRACT (the load-bearing design point): @c operator* returns a
 * reference to a @ref Particle view CACHED INSIDE THE ITERATOR (the
 * @ref m_view member), REBOUND to the current row on each dereference. That
 * reference stays valid for as long as the iterator object lives and is not
 * incremented; incrementing or destroying the iterator invalidates it. This
 * lets a caller bind @c Particle &p = *it and keep using @c p across a loop
 * body (the pattern the bond handlers rely on -- see the phase 7a task-1
 * audit), which a by-value @c operator* returning a temporary could not
 * support. The referenced view itself aliases the store; it is invalidated by
 * a store rebuild just like any other view.
 *
 * PERFORMANCE (phase 7a perf fix): dereference REBINDS the cached view (sets
 * the two store-handle fields via @ref Particle::attach_to_store) instead of
 * CONSTRUCTING + copy-assigning a fresh @ref Particle. The old @c m_view =
 * make_view(...) built a full Particle (all migration carriers, incl. the
 * heap-owning BondList/exclusion members) and copy-assigned it -- two full
 * Particle materialisations per dereference, per element, in every core
 * particle loop (the measured 2-3x regression). Rebinding touches only the two
 * handle fields; the carriers stay default-constructed and are never read
 * while the view is attached (every accessor takes the store-attached branch).
 *
 * The cache is a @c std::optional<Particle> so that DEFAULT / COPY / end
 * construction and @c begin() (all frequent: @ref ParticleIterator reassigns
 * its inner iterator once per cell, @ref Algorithm::link_cell copies iterators
 * via @c std::next and builds a fresh range per neighbour) cost NOTHING beyond
 * the two index fields -- no Particle is materialised until the iterator is
 * actually dereferenced. Exactly one Particle shell is constructed per iterator
 * that is dereferenced at least once, and it is reused (rebound) across every
 * subsequent dereference of that iterator.
 */
class RowParticleIterator
    : public boost::iterator_facade<RowParticleIterator, Particle,
                                    boost::forward_traversal_tag, Particle &> {
  CellRows::const_iterator m_row;
  ParticleStore *m_store;
  /** Cached view over the current row; the target of @c operator*. Lazily
   *  materialised on first dereference (see the class PERFORMANCE note) and
   *  then rebound to the current row on each subsequent dereference. Its
   * address is stable for the lifetime of the iterator once materialised (per
   * the class lifetime contract). Declared @c mutable because @c dereference()
   * is
   *  @c const per the boost iterator_facade contract while the cache is
   *  logically transient. */
  mutable std::optional<Particle> m_view;

public:
  /** @brief Default (singular) iterator, needed by @ref ParticleIterator's
   *  end/singular states. Never dereferenced. */
  RowParticleIterator() : m_row(nullptr), m_store(nullptr) {}

  /** @brief Iterator at @p row over @p store. */
  RowParticleIterator(CellRows::const_iterator row, ParticleStore &store)
      : m_row(row), m_store(&store) {}

  /** @brief Past-the-end iterator (no store dereference happens). */
  explicit RowParticleIterator(CellRows::const_iterator end)
      : m_row(end), m_store(nullptr) {}

  /** @brief The store + row this iterator currently addresses. Lets a composing
   *  iterator (@ref ParticleIterator) rebind its OWN reused view without
   *  materialising this iterator's cache -- avoiding a Particle construction
   * per cell boundary. Valid only on a non-end iterator. */
  ParticleStore *current_store() const { return m_store; }
  int current_row() const { return *m_row; }

private:
  friend class boost::iterator_core_access;

  void increment() { ++m_row; }

  bool equal(RowParticleIterator const &rhs) const {
    return m_row == rhs.m_row;
  }

  // Materialise the cached view on first dereference (so a valid past-the-end
  // iterator, which is never dereferenced, never builds one), then REBIND it to
  // the current row. Rebinding sets ONLY the two store-handle fields
  // (attach_to_store) -- it never copy-assigns a Particle, so the heap-owning
  // migration carriers are neither reallocated nor touched. The view lives in
  // m_view so the returned reference outlives the expression, honouring the
  // lifetime contract documented on the class.
  Particle &dereference() const {
    if (not m_view) {
      m_view.emplace();
      m_view->attach_to_store(*m_store, *m_row);
    } else if (m_view->store() != m_store or m_view->store_row() != *m_row) {
      // Rebind only if not already bound to (m_store, *m_row); a repeat
      // dereference of the same row is then a pure comparison.
      m_view->attach_to_store(*m_store, *m_row);
    }
    return *m_view;
  }
};

/**
 * @brief Range of @ref Particle views over a @ref CellRows bag + store.
 *
 * The live return type of @ref Cell::particles() since the phase-7a flip.
 * Modelled on @ref ParticleRange (a @c boost::iterator_range with a cached
 * size). Iterating it yields the cell's particles in cell-traversal order,
 * each as a live @ref Particle view aliasing the @ref ParticleStore columns.
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
