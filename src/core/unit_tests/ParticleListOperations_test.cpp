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

// Phase 7a: cells hold ParticleStore ROW indices + a staging buffer, not a
// Bag<Particle>. The CellParticleStorage choke points therefore operate on a
// Cell (staging inserts, snapshot+remove extracts, clear, ghost resize) rather
// than on a ParticleList. These tests pin that row/staging behaviour with a
// hand-built store and cell (no MPI / decomposition).

#define BOOST_TEST_MODULE ParticleListOperations test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Particle.hpp"
#include "cell_system/Cell.hpp"
#include "cell_system/ParticleListOperations.hpp"
#include "particle_store/ParticleStore.hpp"

#include <Kokkos_Core.hpp>

#include <cstddef>
#include <vector>

// ParticleStore allocates Kokkos Views, which requires an initialized runtime.
struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

static Particle make_particle(int id) {
  Particle p{};
  p.id() = id;
  return p;
}

namespace {
// Commit a cell's staged particles into store rows, exactly as
// ensure_particle_store_synchronized does for a single-cell system: assign one
// row per staged particle (seeded from carriers) and record it in the row bag.
void commit(Cell &cell, ParticleStore &store) {
  auto const n = cell.staged().size();
  store.mark_dirty();
  store.begin_rebuild(n, 0u);
  int row = 0;
  cell.rows().clear();
  for (auto &p : cell.staged()) {
    store.assign_row(p, row);
    cell.rows().insert(row);
    ++row;
  }
  store.finish_rebuild();
  cell.staged().clear();
  cell.set_store(store);
}
} // namespace

BOOST_AUTO_TEST_CASE(insert_particle_stages_and_returns_reference) {
  Cell cell;
  auto &staged = CellParticleStorage::insert_particle(cell, make_particle(7));
  // Not committed yet: no rows, one staged particle.
  BOOST_CHECK_EQUAL(cell.rows().size(), 0ul);
  BOOST_CHECK_EQUAL(cell.staged().size(), 1ul);
  BOOST_CHECK_EQUAL(staged.id(), 7);
}

BOOST_AUTO_TEST_CASE(extract_row_snapshots_and_removes_with_swap_with_back) {
  ParticleStore store{};
  Cell cell;
  CellParticleStorage::insert_particle(cell, make_particle(1));
  CellParticleStorage::insert_particle(cell, make_particle(2));
  CellParticleStorage::insert_particle(cell, make_particle(3));
  commit(cell, store);
  BOOST_REQUIRE_EQUAL(cell.rows().size(), 3ul);

  auto extracted = CellParticleStorage::extract_row(cell, 0u);
  BOOST_CHECK_EQUAL(extracted.id(), 1);
  BOOST_CHECK_EQUAL(cell.rows().size(), 2ul);
  // The snapshot is detached and carries the row's data via its carriers.
  BOOST_CHECK(extracted.store() == nullptr);
  // swap-with-back: the row bag's freed position now holds the last row.
  BOOST_CHECK_EQUAL(store.id(cell.rows().begin()[0]),
                    store.id(2)); // row 2 held id 3
}

BOOST_AUTO_TEST_CASE(extract_row_on_single_element_empties_row_bag) {
  ParticleStore store{};
  Cell cell;
  CellParticleStorage::insert_particle(cell, make_particle(42));
  commit(cell, store);

  auto extracted = CellParticleStorage::extract_row(cell, 0u);
  BOOST_CHECK_EQUAL(extracted.id(), 42);
  BOOST_CHECK_EQUAL(cell.rows().size(), 0ul);
}

BOOST_AUTO_TEST_CASE(clear_particles_empties_rows_and_staging) {
  ParticleStore store{};
  Cell cell;
  CellParticleStorage::insert_particle(cell, make_particle(1));
  commit(cell, store);
  CellParticleStorage::insert_particle(cell, make_particle(2)); // staged
  CellParticleStorage::clear_particles(cell);
  BOOST_CHECK_EQUAL(cell.rows().size(), 0ul);
  BOOST_CHECK_EQUAL(cell.staged().size(), 0ul);
}

BOOST_AUTO_TEST_CASE(resize_ghost_storage_stages_ghosts) {
  Cell cell;
  CellParticleStorage::insert_particle(cell, make_particle(1));
  CellParticleStorage::resize_ghost_storage(cell, 3ul);
  // Stages exactly `count` default ghost particles; drops any prior content.
  BOOST_CHECK_EQUAL(cell.rows().size(), 0ul);
  BOOST_CHECK_EQUAL(cell.staged().size(), 3ul);
  for (auto const &p : cell.staged()) {
    BOOST_CHECK(p.is_ghost());
  }
}

BOOST_AUTO_TEST_CASE(resize_ghost_storage_shrinks) {
  Cell cell;
  CellParticleStorage::insert_particle(cell, make_particle(1));
  CellParticleStorage::insert_particle(cell, make_particle(2));
  CellParticleStorage::insert_particle(cell, make_particle(3));

  CellParticleStorage::resize_ghost_storage(cell, 1ul);
  BOOST_CHECK_EQUAL(cell.staged().size(), 1ul);
  for (auto const &p : cell.staged()) {
    BOOST_CHECK(p.is_ghost());
  }
}
