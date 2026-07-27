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

/**
 * @file
 * TDD for Task 3.1: interior/boundary cell classification.
 *
 * Tests:
 *  1. mark_boundary_cells: a local cell with only local neighbors → interior.
 *  2. mark_boundary_cells: a local cell with a ghost neighbor → boundary.
 *  3. mark_boundary_cells is idempotent (second call gives same result).
 *  4. validate_halo_plan consistency check: fires when a cell is mis-marked
 *     interior but has a ghost neighbor.
 *  5. validate_halo_plan consistency check: silent for a correctly-marked
 *     boundary cell.
 */

#define BOOST_TEST_MODULE Cell_boundary test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "cell_system/Cell.hpp"
#include "ghosts/HaloPlan.hpp"
#include "ghosts/HaloPlanValidator.hpp"
#include "ghosts/mark_boundary_cells.hpp"

#include <span>
#include <vector>

// ── mark_boundary_cells ──────────────────────────────────────────────────────

// A local cell whose only neighbors are other local cells must be interior.
BOOST_AUTO_TEST_CASE(local_only_neighbors_is_interior) {
  Cell local_a, local_b, local_c;
  // local_a has two neighbors, both local (not ghost)
  local_a.m_neighbors =
      Neighbors<Cell *>(std::vector<Cell *>{&local_b, &local_c}, {});

  std::vector<Cell *> locals{&local_a, &local_b, &local_c};
  std::vector<Cell *> ghosts{};

  GhostComm::mark_boundary_cells(locals, ghosts);

  BOOST_CHECK_MESSAGE(!local_a.is_boundary(),
                      "cell with only local neighbors must be interior");
  BOOST_CHECK_MESSAGE(!local_b.is_boundary(),
                      "local_b (no neighbors listed) must be interior");
  BOOST_CHECK_MESSAGE(!local_c.is_boundary(),
                      "local_c (no neighbors listed) must be interior");
}

// A local cell that has at least one ghost neighbor must be boundary.
BOOST_AUTO_TEST_CASE(ghost_neighbor_makes_boundary) {
  Cell local_a, local_b, ghost;
  // local_a has one local neighbor and one ghost neighbor
  local_a.m_neighbors =
      Neighbors<Cell *>(std::vector<Cell *>{&local_b, &ghost}, {});
  // local_b has only a local neighbor
  local_b.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&local_a}, {});

  std::vector<Cell *> locals{&local_a, &local_b};
  std::vector<Cell *> ghosts{&ghost};

  GhostComm::mark_boundary_cells(locals, ghosts);

  BOOST_CHECK_MESSAGE(local_a.is_boundary(),
                      "cell with a ghost neighbor must be boundary");
  BOOST_CHECK_MESSAGE(!local_b.is_boundary(),
                      "cell without ghost neighbors must be interior");
}

// Calling mark_boundary_cells a second time (rebuild) must yield the same
// result (reset + re-mark = idempotent).
BOOST_AUTO_TEST_CASE(mark_boundary_cells_is_idempotent) {
  Cell local, ghost;
  local.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&ghost}, {});
  std::vector<Cell *> locals{&local};
  std::vector<Cell *> ghosts{&ghost};

  GhostComm::mark_boundary_cells(locals, ghosts);
  bool const first_result = local.is_boundary();

  GhostComm::mark_boundary_cells(locals, ghosts);
  bool const second_result = local.is_boundary();

  BOOST_CHECK_EQUAL(first_result, second_result);
  BOOST_CHECK(first_result); // should be boundary
}

// ── validator consistency check ──────────────────────────────────────────────

// A local cell whose is_boundary()==false but that has a ghost neighbor must
// trigger the "interior cell has a ghost neighbor" consistency violation.
BOOST_AUTO_TEST_CASE(validator_fires_for_mismarked_interior_cell) {
  using namespace GhostComm;
  Cell local, ghost;
  // Wire ghost as a neighbor of local, but do NOT call mark_boundary_cells,
  // so local.is_boundary() stays false (the default).
  local.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&ghost}, {});

  std::vector<Cell *> locals{&local};
  std::vector<Cell *> ghosts{&ghost};

  // Build a minimal valid plan that covers the ghost so that the coverage
  // checks pass — we only want to test the consistency check in isolation.
  HaloPlan plan;
  plan.neighbors.push_back(NeighborComm{
      1, {SendRegion{&local.particles(), {}}}, {&ghost.particles()}});

  auto violations = validate_halo_plan(plan, locals, ghosts);
  bool found_consistency =
      std::ranges::any_of(violations, [](std::string const &v) {
        return v.find("interior cell has a ghost neighbor") !=
               std::string::npos;
      });
  BOOST_CHECK_MESSAGE(
      found_consistency,
      "expected 'interior cell has a ghost neighbor' violation; "
      "got: " +
          (violations.empty() ? "(none)" : violations.front()));
}

// A correctly-marked boundary cell must NOT trigger the consistency violation.
BOOST_AUTO_TEST_CASE(validator_silent_for_correctly_marked_boundary_cell) {
  using namespace GhostComm;
  Cell local, ghost;
  local.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&ghost}, {});

  std::vector<Cell *> locals{&local};
  std::vector<Cell *> ghosts{&ghost};

  // Properly mark the cell as boundary first.
  GhostComm::mark_boundary_cells(locals, ghosts);
  BOOST_REQUIRE(local.is_boundary()); // pre-condition for this test

  HaloPlan plan;
  plan.neighbors.push_back(NeighborComm{
      1, {SendRegion{&local.particles(), {}}}, {&ghost.particles()}});

  auto violations = validate_halo_plan(plan, locals, ghosts);
  bool found_consistency =
      std::ranges::any_of(violations, [](std::string const &v) {
        return v.find("interior cell has a ghost neighbor") !=
               std::string::npos;
      });
  BOOST_CHECK_MESSAGE(
      !found_consistency,
      "expected no consistency violation for a correctly-marked "
      "boundary cell; got: " +
          (violations.empty() ? "(none)" : violations.front()));
  BOOST_CHECK_MESSAGE(violations.empty(),
                      "expected no violations at all for a correct plan");
}
