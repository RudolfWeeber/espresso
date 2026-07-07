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

// Phase 7a/7b: cells hold ParticleStore ROW indices + a staging buffer, not a
// Bag<Particle>. Phase 7b (Task 4): the migration envelope died, so a cell no
// longer stages a detached, data-carrying Particle -- it stages a
// StagedParticle (a reference to a SOURCE-store row, or a fresh-default
// marker). The CellParticleStorage choke points therefore operate on a Cell
// (row-ref staging inserts, row drops, clear, ghost resize) rather than on a
// ParticleList. These tests pin that row/staging behaviour with a hand-built
// store and cell (no MPI / decomposition).

#define BOOST_TEST_MODULE ParticleListOperations test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Particle.hpp"
#include "cell_system/Cell.hpp"
#include "cell_system/ParticleListOperations.hpp"
#include "particle_store/ParticleStore.hpp"

#include <Kokkos_Core.hpp>

#include <cstddef>

// ParticleStore allocates Kokkos Views, which requires an initialized runtime.
struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

namespace {
// A source store that outlives the commit, holding one seeded row per particle
// a test wants to stage. Rows are handed out monotonically; each is seeded to
// the new-particle defaults and its id set, so a committed copy_row preserves
// it.
struct SourceStore {
  ParticleStore store{};
  int next_row = 0;

  explicit SourceStore(std::size_t capacity) {
    store.begin_rebuild(capacity, 0u);
    for (std::size_t i = 0u; i < capacity; ++i) {
      store.seed_default_row(static_cast<int>(i));
    }
    store.finish_rebuild();
  }
  // Seed the next source row with the given id and return it.
  int make_source_row(int id) {
    auto const row = next_row++;
    store.id(row) = id;
    return row;
  }
};

// Commit a cell's staged rows into store rows, mirroring the staged loop of
// CellStructure::ensure_particle_store_synchronized for a single-cell system:
// one committed row per staged entry, copied from its source row (copy_row) or
// seeded to defaults (fresh-default ghost, source_store == nullptr).
void commit(Cell &cell, ParticleStore &store) {
  auto const n = cell.rows().size() + cell.staged().size();
  store.mark_dirty();
  store.begin_rebuild(n, 0u);
  int row = 0;
  cell.rows().clear();
  for (auto const &staged : cell.staged()) {
    if (staged.source_store != nullptr) {
      store.copy_row(*staged.source_store, staged.source_row, row);
    } else {
      store.seed_default_row(row);
    }
    cell.rows().insert(row);
    ++row;
  }
  store.finish_rebuild();
  cell.staged().clear();
  cell.set_store(store);
}
} // namespace

BOOST_AUTO_TEST_CASE(insert_staged_row_stages_a_row_reference) {
  SourceStore src{1u};
  auto const row = src.make_source_row(7);

  Cell cell;
  CellParticleStorage::insert_staged_row(cell, src.store, row);
  // Not committed yet: no rows, one staged entry that references the source
  // row.
  BOOST_CHECK_EQUAL(cell.rows().size(), 0ul);
  BOOST_REQUIRE_EQUAL(cell.staged().size(), 1ul);
  BOOST_CHECK(cell.staged()[0].source_store == &src.store);
  BOOST_CHECK_EQUAL(cell.staged()[0].source_row, row);
}

BOOST_AUTO_TEST_CASE(drop_row_removes_with_swap_with_back) {
  SourceStore src{3u};
  ParticleStore store{};
  Cell cell;
  CellParticleStorage::insert_staged_row(cell, src.store,
                                         src.make_source_row(1));
  CellParticleStorage::insert_staged_row(cell, src.store,
                                         src.make_source_row(2));
  CellParticleStorage::insert_staged_row(cell, src.store,
                                         src.make_source_row(3));
  commit(cell, store);
  BOOST_REQUIRE_EQUAL(cell.rows().size(), 3ul);
  // Committed rows 0,1,2 hold ids 1,2,3 respectively.
  BOOST_REQUIRE_EQUAL(store.id(0), 1);
  BOOST_REQUIRE_EQUAL(store.id(2), 3);

  CellParticleStorage::drop_row(cell, 0u);
  BOOST_CHECK_EQUAL(cell.rows().size(), 2ul);
  // swap-with-back: the row bag's freed position now holds the former last row
  // (row 2, which held id 3).
  BOOST_CHECK_EQUAL(store.id(cell.rows().begin()[0]), store.id(2));
}

BOOST_AUTO_TEST_CASE(drop_row_on_single_element_empties_row_bag) {
  SourceStore src{1u};
  ParticleStore store{};
  Cell cell;
  CellParticleStorage::insert_staged_row(cell, src.store,
                                         src.make_source_row(42));
  commit(cell, store);
  BOOST_REQUIRE_EQUAL(cell.rows().size(), 1ul);

  CellParticleStorage::drop_row(cell, 0u);
  BOOST_CHECK_EQUAL(cell.rows().size(), 0ul);
}

BOOST_AUTO_TEST_CASE(clear_particles_empties_rows_and_staging) {
  SourceStore src{2u};
  ParticleStore store{};
  Cell cell;
  CellParticleStorage::insert_staged_row(cell, src.store,
                                         src.make_source_row(1));
  commit(cell, store);
  CellParticleStorage::insert_staged_row(cell, src.store,
                                         src.make_source_row(2)); // staged
  CellParticleStorage::clear_particles(cell);
  BOOST_CHECK_EQUAL(cell.rows().size(), 0ul);
  BOOST_CHECK_EQUAL(cell.staged().size(), 0ul);
}

BOOST_AUTO_TEST_CASE(resize_ghost_storage_stages_ghosts) {
  SourceStore src{1u};
  Cell cell;
  CellParticleStorage::insert_staged_row(cell, src.store,
                                         src.make_source_row(1));
  CellParticleStorage::resize_ghost_storage(cell, 3ul);
  // Stages exactly `count` fresh-default ghost entries; drops any prior
  // content.
  BOOST_CHECK_EQUAL(cell.rows().size(), 0ul);
  BOOST_REQUIRE_EQUAL(cell.staged().size(), 3ul);
  // Fresh-default ghost: null source store (the rebuild seeds defaults).
  for (auto const &staged : cell.staged()) {
    BOOST_CHECK(staged.source_store == nullptr);
  }
}

BOOST_AUTO_TEST_CASE(resize_ghost_storage_shrinks) {
  SourceStore src{3u};
  Cell cell;
  CellParticleStorage::insert_staged_row(cell, src.store,
                                         src.make_source_row(1));
  CellParticleStorage::insert_staged_row(cell, src.store,
                                         src.make_source_row(2));
  CellParticleStorage::insert_staged_row(cell, src.store,
                                         src.make_source_row(3));

  CellParticleStorage::resize_ghost_storage(cell, 1ul);
  BOOST_REQUIRE_EQUAL(cell.staged().size(), 1ul);
  for (auto const &staged : cell.staged()) {
    BOOST_CHECK(staged.source_store == nullptr);
  }
}
