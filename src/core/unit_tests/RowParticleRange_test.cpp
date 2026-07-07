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

// Phase-7a groundwork tests for the DORMANT view machinery: the per-cell
// row-index bag (CellRows) fill order, and the RowParticleRange iterator
// adaptor that yields Particle views from a CellRows + ParticleStore. All are
// standalone (no MPI, no decomposition): a ParticleStore is built by hand and a
// Cell's Bag + row bag are populated the way ensure_particle_store_synchronized
// does, so the tests pin the machinery's contract without the integration path
// (which the physics ctest battery exercises end-to-end under
// ADDITIONAL_CHECKS).

#define BOOST_TEST_MODULE RowParticleRange test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "cell_system/Cell.hpp"
#include "cell_system/CellRows.hpp"
#include "cell_system/RowParticleRange.hpp"

#include "Particle.hpp"
#include "particle_store/ParticleStore.hpp"

#include <Kokkos_Core.hpp>

#include <cstddef>
#include <span>
#include <vector>

// ParticleStore allocates Kokkos Views, which requires an initialized runtime.
struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

namespace {
// Build a store with `n` rows and set row r's id to `first_id + r`. Returns the
// detached seed particles kept alive so their carriers back the columns; the
// caller only needs the store afterwards.
void build_store(ParticleStore &store, std::vector<Particle> &seeds,
                 std::size_t n_local, std::size_t n_ghost, int first_id) {
  auto const n = n_local + n_ghost;
  seeds.clear();
  seeds.resize(n);
  store.mark_dirty();
  store.begin_rebuild(n_local, n_ghost);
  int row = 0;
  for (auto &p : seeds) {
    store.assign_row(p, row++);
  }
  store.finish_rebuild();
  for (int r = 0; r < static_cast<int>(n); ++r) {
    store.id(r) = first_id + r;
  }
}

// Fill a cell's row bag and wire its store exactly as
// ensure_particle_store_synchronized does (phase 7a: the cell no longer holds a
// Bag<Particle>; it holds row indices + a store pointer, and hands out views).
// `rows` are the store rows assigned to this cell.
void fill_cell(Cell &cell, ParticleStore &store, std::vector<int> const &rows) {
  cell.set_store(store);
  auto &row_bag = cell.rows();
  row_bag.clear();
  for (auto const r : rows) {
    row_bag.insert(r);
  }
}
} // namespace

// The per-cell row bag must record the exact store rows of the cell's Bag
// contents, in Bag iteration order. Mirrors the ensure_particle_store_
// synchronized loop for both a "local" and a "ghost" cell block.
BOOST_AUTO_TEST_CASE(cell_rows_match_bag_traversal_order) {
  ParticleStore store{};
  std::vector<Particle> seeds;
  // 3 local rows [0,3) then 2 ghost rows [3,5), ids 100..104.
  build_store(store, seeds, 3u, 2u, 100);

  Cell local_cell, ghost_cell;
  fill_cell(local_cell, store, {0, 1, 2});
  fill_cell(ghost_cell, store, {3, 4});

  auto check = [&store](Cell const &cell) {
    auto const &bag = cell.particles();
    auto const &rows = cell.rows();
    BOOST_REQUIRE_EQUAL(rows.size(), bag.size());
    std::size_t k = 0u;
    for (auto const &p : bag) {
      // store.id(rows[k]) == bag[k].id(), the ADDITIONAL_CHECKS invariant.
      BOOST_CHECK_EQUAL(store.id(rows.begin()[k]), p.id());
      ++k;
    }
  };
  check(local_cell);
  check(ghost_cell);

  // Spell out the recorded rows to guard the helper itself.
  BOOST_REQUIRE_EQUAL(local_cell.rows().size(), 3u);
  BOOST_CHECK_EQUAL(local_cell.rows().begin()[0], 0);
  BOOST_CHECK_EQUAL(local_cell.rows().begin()[2], 2);
  BOOST_REQUIRE_EQUAL(ghost_cell.rows().size(), 2u);
  BOOST_CHECK_EQUAL(ghost_cell.rows().begin()[0], 3);
  BOOST_CHECK_EQUAL(ghost_cell.rows().begin()[1], 4);
}

// The iterator adaptor must yield the same particle sequence (by id) as
// iterating the cell's Bag, on a built store.
BOOST_AUTO_TEST_CASE(row_range_sequence_equals_bag_iteration) {
  ParticleStore store{};
  std::vector<Particle> seeds;
  build_store(store, seeds, 4u, 0u, 200);

  Cell cell;
  fill_cell(cell, store, {0, 1, 2, 3});

  std::vector<int> bag_ids;
  for (auto const &p : cell.particles()) {
    bag_ids.push_back(p.id());
  }

  std::vector<int> range_ids;
  RowParticleRange range{cell.rows(), store};
  BOOST_CHECK_EQUAL(range.size(), cell.rows().size());
  for (auto const &p : range) {
    range_ids.push_back(p.id());
  }

  BOOST_CHECK(range_ids == bag_ids);
  BOOST_REQUIRE_EQUAL(range_ids.size(), 4u);
  BOOST_CHECK_EQUAL(range_ids[0], 200);
  BOOST_CHECK_EQUAL(range_ids[3], 203);
}

// An EMPTY row bag yields an empty range (begin() == end()).
BOOST_AUTO_TEST_CASE(row_range_empty) {
  ParticleStore store{};
  std::vector<Particle> seeds;
  build_store(store, seeds, 1u, 0u, 0);

  Cell cell; // no rows filled
  RowParticleRange range{cell.rows(), store};
  BOOST_CHECK_EQUAL(range.size(), 0u);
  BOOST_CHECK(range.begin() == range.end());
  std::size_t count = 0u;
  for (auto const &p : range) {
    static_cast<void>(p);
    ++count;
  }
  BOOST_CHECK_EQUAL(count, 0u);
}

// LIFETIME CONTRACT: `Particle &p = *it` stays valid while the iterator lives
// and is not incremented. Binding a reference to the view, then reading it
// after other statements, must still refer to the same particle. This is the
// pattern the bond handlers rely on (Particle &p2 = *ptr held across a call).
BOOST_AUTO_TEST_CASE(
    row_range_dereference_reference_stable_while_iterator_lives) {
  ParticleStore store{};
  std::vector<Particle> seeds;
  build_store(store, seeds, 3u, 0u, 300);

  Cell cell;
  fill_cell(cell, store, {0, 1, 2});

  RowParticleRange range{cell.rows(), store};
  auto it = range.begin();

  // Bind a reference to the current view and keep using it across other
  // statements without re-dereferencing the iterator.
  Particle &p = *it;
  int const id_at_bind = p.id();
  BOOST_CHECK_EQUAL(id_at_bind, 300);

  // Do unrelated work; the reference must remain valid and refer to the same
  // particle (row 0), i.e. the referent was not a destroyed temporary.
  volatile int scratch = 0;
  for (int i = 0; i < 1000; ++i) {
    scratch += i;
  }
  static_cast<void>(scratch);

  BOOST_CHECK_EQUAL(p.id(), 300);
  BOOST_CHECK_EQUAL(p.store_row(), 0);
  BOOST_CHECK(p.store() == &store);
  // A write through the still-live reference reaches the store row it names.
  p.pos() = Utils::Vector3d{1., 2., 3.};
  BOOST_CHECK_EQUAL(store.position_value(0)[1], 2.);

  // After incrementing, the SAME iterator now names the next row; a fresh
  // dereference reflects the advance (the cache is refreshed, per the
  // documented contract that incrementing invalidates the previous referent).
  ++it;
  Particle &p_next = *it;
  BOOST_CHECK_EQUAL(p_next.id(), 301);
  BOOST_CHECK_EQUAL(p_next.store_row(), 1);
}
