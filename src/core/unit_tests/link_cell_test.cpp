/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE link_cell test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "algorithm/link_cell.hpp"

#include "Particle.hpp"
#include "cell_system/Cell.hpp"

#include <unordered_set>
#include <utility>
#include <vector>

// Custom hash function for pair<int, int>
struct pair_hash {
  std::size_t operator()(const std::pair<int, int> &p) const {
    auto h1 = std::hash<int>{}(p.first);
    auto h2 = std::hash<int>{}(p.second);
    return h1 ^ (h2 << 1);
  }
};

BOOST_AUTO_TEST_CASE(link_cell) {
  auto const n_cells = 10u;
  auto const n_part_per_cell = 10u;
  auto const n_part = n_cells * n_part_per_cell;

  std::vector<Cell> cells(n_cells);

  auto id = 0;
  for (auto &c : cells) {
    std::vector<Cell *> neighbors;

    for (auto &n : cells) {
      if (&c != &n)
        neighbors.push_back(&n);
    }

    c.m_neighbors = Neighbors<Cell *>(neighbors, {});

    c.particles().resize(n_part_per_cell);

    for (auto &p : c.particles()) {
      p.id() = id++;
    }
  }

  // Collect all pairs found by the link_cell algorithm
  std::unordered_set<std::pair<int, int>, pair_hash> found_pairs;

  Algorithm::link_cell(cells.begin(), cells.end(),
                       [&found_pairs](Particle const &p1, Particle const &p2) {
                         // Store pairs in normalized order (smaller id first)
                         if (p1.id() < p2.id())
                           found_pairs.emplace(p1.id(), p2.id());
                         else if (p2.id() < p1.id())
                           found_pairs.emplace(p2.id(), p1.id());
                         // Skip if p1.id() == p2.id()
                       });

  // Generate expected pairs
  std::unordered_set<std::pair<int, int>, pair_hash> expected_pairs;
  for (auto i = 0; i < static_cast<int>(n_part); i++)
    for (auto j = i + 1; j < static_cast<int>(n_part); j++) {
      expected_pairs.emplace(i, j);
    }

  // Check that we found exactly the expected number of pairs
  BOOST_CHECK_EQUAL(found_pairs.size(), expected_pairs.size());

  // Check that the sets are equal (all expected pairs were found)
  BOOST_CHECK(found_pairs == expected_pairs);
}
