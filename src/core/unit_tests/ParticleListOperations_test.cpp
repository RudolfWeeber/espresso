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

#define BOOST_TEST_MODULE ParticleListOperations test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Particle.hpp"
#include "ParticleList.hpp"
#include "cell_system/ParticleListOperations.hpp"

#include <utility>

static Particle make_particle(int id) {
  Particle p{};
  p.id() = id;
  return p;
}

BOOST_AUTO_TEST_CASE(insert_particle_appends_and_returns_reference) {
  ParticleList storage;
  auto &stored =
      CellParticleStorage::insert_particle(storage, make_particle(7));
  BOOST_CHECK_EQUAL(storage.size(), 1ul);
  BOOST_CHECK_EQUAL(stored.id(), 7);
}

BOOST_AUTO_TEST_CASE(extract_particle_moves_out_with_swap_with_back) {
  ParticleList storage;
  CellParticleStorage::insert_particle(storage, make_particle(1));
  CellParticleStorage::insert_particle(storage, make_particle(2));
  CellParticleStorage::insert_particle(storage, make_particle(3));

  auto [extracted, next] =
      CellParticleStorage::extract_particle(storage, storage.begin());
  BOOST_CHECK_EQUAL(extracted.id(), 1);
  BOOST_CHECK_EQUAL(storage.size(), 2ul);
  // swap-with-back: the last element (id 3) now occupies the freed slot
  BOOST_CHECK_EQUAL(next->id(), 3);
}

BOOST_AUTO_TEST_CASE(extract_particle_on_single_element_returns_end) {
  ParticleList storage;
  CellParticleStorage::insert_particle(storage, make_particle(42));

  auto [extracted, next] =
      CellParticleStorage::extract_particle(storage, storage.begin());
  // erase-last path of the underlying Bag: nothing to swap in
  BOOST_CHECK_EQUAL(extracted.id(), 42);
  BOOST_CHECK_EQUAL(storage.size(), 0ul);
  BOOST_CHECK(next == storage.end());
}

BOOST_AUTO_TEST_CASE(erase_particle_destroys_element) {
  ParticleList storage;
  CellParticleStorage::insert_particle(storage, make_particle(1));
  CellParticleStorage::insert_particle(storage, make_particle(2));
  auto it = CellParticleStorage::erase_particle(storage, storage.begin());
  BOOST_CHECK_EQUAL(storage.size(), 1ul);
  BOOST_CHECK_EQUAL(it->id(), 2);
}

BOOST_AUTO_TEST_CASE(clear_particles_empties_storage) {
  ParticleList storage;
  CellParticleStorage::insert_particle(storage, make_particle(1));
  CellParticleStorage::clear_particles(storage);
  BOOST_CHECK_EQUAL(storage.size(), 0ul);
}

BOOST_AUTO_TEST_CASE(resize_ghost_storage_marks_all_particles_as_ghosts) {
  ParticleList storage;
  CellParticleStorage::insert_particle(storage, make_particle(1));
  CellParticleStorage::resize_ghost_storage(storage, 3ul);
  BOOST_CHECK_EQUAL(storage.size(), 3ul);
  for (auto const &p : storage) {
    BOOST_CHECK(p.is_ghost());
  }
}

BOOST_AUTO_TEST_CASE(resize_ghost_storage_shrinks_and_marks_remaining_ghost) {
  ParticleList storage;
  CellParticleStorage::insert_particle(storage, make_particle(1));
  CellParticleStorage::insert_particle(storage, make_particle(2));
  CellParticleStorage::insert_particle(storage, make_particle(3));

  CellParticleStorage::resize_ghost_storage(storage, 1ul);
  BOOST_CHECK_EQUAL(storage.size(), 1ul);
  for (auto const &p : storage) {
    BOOST_CHECK(p.is_ghost());
  }
}
