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

#define BOOST_TEST_MODULE HaloPlanValidator test
#define BOOST_TEST_DYN_LINK
#include "EspressoCoreGlobalConfig.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "cell_system/AtomDecomposition.hpp"
#include "cell_system/Cell.hpp"
#include "cell_system/HybridDecomposition.hpp"
#include "cell_system/RegularDecomposition.hpp"
#include "communication.hpp"
#include "ghosts/HaloPlan.hpp"
#include "ghosts/HaloPlanValidator.hpp"

#include <boost/mpi.hpp>
#include <boost/test/unit_test.hpp>

#include <functional>
#include <optional>
#include <set>
#include <span>
#include <vector>

BOOST_TEST_GLOBAL_FIXTURE(EspressoCoreGlobalConfig);

namespace {
BoxGeometry make_periodic_box(double box_l) {
  BoxGeometry box;
  box.set_length(Utils::Vector3d::broadcast(box_l));
  box.set_periodic(0u, true);
  box.set_periodic(1u, true);
  box.set_periodic(2u, true);
  return box;
}

RegularDecomposition make_dd(Utils::Vector3i const &node_grid, double box_l,
                             double range) {
  ::communicator.set_node_grid(node_grid);
  auto const box = make_periodic_box(box_l);
  auto const local_box = LocalBox::make_regular_decomposition(
      box.length(), ::communicator.calc_node_index(), ::communicator.node_grid);
  return RegularDecomposition(::communicator.comm, range, box, local_box,
                              std::nullopt);
}

// The AtomDecomposition (n-square) covers every ghost via its collective
// broadcast/reduce section, so the validator must recognise collective-covered
// ghosts as covered. The box lifetime must outlive the returned decomposition;
// AtomDecomposition holds a BoxGeometry const& (see AtomDecomposition.hpp).
AtomDecomposition make_atom_dd(BoxGeometry const &box) {
  return AtomDecomposition(::communicator.comm, box);
}

// The HybridDecomposition combines the regular child's point-to-point section
// with the n-square child's collective section, so its plan exercises BOTH
// coverage paths at once.
HybridDecomposition make_hybrid_dd(Utils::Vector3i const &node_grid,
                                   BoxGeometry const &box,
                                   double cutoff_regular, double skin,
                                   std::set<int> const &n_square_types) {
  ::communicator.set_node_grid(node_grid);
  auto const local_box = LocalBox::make_regular_decomposition(
      box.length(), ::communicator.calc_node_index(), ::communicator.node_grid);
  return HybridDecomposition(
      ::communicator.comm, cutoff_regular, skin, []() { return false; }, box,
      local_box, n_square_types);
}
} // namespace

namespace utf = boost::unit_test;
static bool has_4_mpi_ranks() { return boost::mpi::communicator{}.size() == 4; }
static bool has_1_mpi_rank() { return boost::mpi::communicator{}.size() == 1; }

// ── Task 2.1: single-rank local checks ──────────────────────────────────────
// These use hand-made cells and are rank-independent, so they pass at NUM_PROC
// 4 without modification.

BOOST_AUTO_TEST_CASE(detects_defects) {
  using namespace GhostComm;
  // 1 local cell whose only ghost neighbor is `ghost`.
  Cell local, ghost;
  local.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&ghost}, {});
  std::vector<Cell *> locals{&local}, ghosts{&ghost};

  HaloPlan good;
  good.neighbors.push_back(NeighborComm{
      1, {SendRegion{&local.particles(), {}}}, {&ghost.particles()}});
  BOOST_CHECK(validate_halo_plan(good, locals, ghosts).empty());

  // uncovered ghost (neighborship-match + coverage failure)
  HaloPlan uncovered; // empty plan, ghost never filled
  BOOST_CHECK(!validate_halo_plan(uncovered, locals, ghosts).empty());

  // duplicate peer
  HaloPlan duppeer = good;
  duppeer.neighbors.push_back(NeighborComm{
      1, {SendRegion{&local.particles(), {}}}, {&ghost.particles()}});
  BOOST_CHECK(!validate_halo_plan(duppeer, locals, ghosts).empty());

  // double-filled ghost
  HaloPlan dbl = good;
  dbl.neighbors.push_back(NeighborComm{
      2, {SendRegion{&local.particles(), {}}}, {&ghost.particles()}});
  BOOST_CHECK(!validate_halo_plan(dbl, locals, ghosts).empty());

  // send/recv size mismatch
  HaloPlan badshape = good;
  badshape.neighbors[0].recv.clear();
  BOOST_CHECK(!validate_halo_plan(badshape, locals, ghosts).empty());

  // recv target outside ghost set
  Cell alien;
  HaloPlan alien_recv = good;
  alien_recv.neighbors[0].recv[0] = &alien.particles();
  auto violations = validate_halo_plan(alien_recv, locals, ghosts);
  BOOST_CHECK_MESSAGE(!violations.empty(),
                      "Expected recv-outside-ghost violation");

  // isolated shape mismatch (send.size() != recv.size())
  HaloPlan shape_only = good;
  shape_only.neighbors[0].send.push_back(SendRegion{&local.particles(), {}});
  // send.size()==2, recv.size()==1 -> shape violation; ghost still covered once
  auto shape_violations = validate_halo_plan(shape_only, locals, ghosts);
  BOOST_CHECK_MESSAGE(!shape_violations.empty(),
                      "Expected shape-mismatch violation");
}

// ── Task 2.3: real-plan local checks ────────────────────────────────────────
// Build a real RegularDecomposition on 4 ranks and assert that the local
// validator (coverage + neighborship-match + peer-uniqueness + shape) passes.
// This is the guard that the wiring in RegularDecomposition.cpp won't
// false-trip in CI (which builds with ESPRESSO_ADDITIONAL_CHECKS on).
BOOST_AUTO_TEST_CASE(local_checks_pass_real_decomposition,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_4_mpi_ranks();
                     })) {
  using namespace GhostComm;
  auto const dd = make_dd({2, 2, 1}, 8., 1.2);
  auto const *plan = dd.halo_plan();
  BOOST_REQUIRE(plan != nullptr);
  auto violations =
      validate_halo_plan(*plan, dd.local_cells(), dd.ghost_cells());
  BOOST_CHECK_MESSAGE(
      violations.empty(),
      "Expected no local-check violations for real decomposition");
}

// Single-rank real RegularDecomposition. make_halo_plan() returns an empty
// plan on one rank (periodic neighbourships are wired directly between local
// cells in init_cell_interactions), yet mark_cells() still classifies the
// halo-layer cells as ghosts. Those ghosts are never referenced by a local
// cell's neighbor stencil, so the validator must NOT flag them as uncovered.
// This is the regression that aborted particle_reduction_test,
// p3m_test_1_mpi_ranks, EspressoSystemStandAlone_test_1_mpi_ranks etc. under a
// checks-on (RelWithAssert) build. Runs only when launched on a single rank.
BOOST_AUTO_TEST_CASE(local_checks_pass_single_rank_decomposition,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_1_mpi_rank();
                     })) {
  using namespace GhostComm;
  auto const dd = make_dd({1, 1, 1}, 8., 1.2);
  auto const *plan = dd.halo_plan();
  BOOST_REQUIRE(plan != nullptr);
  // Sanity: the single-rank plan is intentionally empty (no p2p, no local),
  // but there are ghost cells present.
  BOOST_CHECK(plan->neighbors.empty());
  BOOST_CHECK(plan->local.empty());
  BOOST_CHECK(not dd.ghost_cells().empty());
  auto violations =
      validate_halo_plan(*plan, dd.local_cells(), dd.ghost_cells());
  BOOST_CHECK_MESSAGE(violations.empty(),
                      "Expected no local-check violations for a single-rank "
                      "RegularDecomposition (unreferenced halo-layer ghosts)");
}

// ── Fix 2 guard: collective-covered ghosts ──────────────────────────────────
// A real AtomDecomposition (n-square) covers ALL of its ghosts via the
// collective broadcast/reduce section, not via point-to-point recv/dst
// targets. Before the fix the validator false-positived on every such ghost
// ("never filled"). These cases run the validator directly in the checks-OFF
// build, so they catch the false-positives without needing a checks-on build.
BOOST_AUTO_TEST_CASE(atom_decomposition_collective_coverage,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_4_mpi_ranks();
                     })) {
  using namespace GhostComm;
  ::communicator.set_node_grid({4, 1, 1});
  auto const box = make_periodic_box(8.);
  auto const dd = make_atom_dd(box);
  auto const *plan = dd.halo_plan();
  BOOST_REQUIRE(plan != nullptr);
  auto violations =
      validate_halo_plan(*plan, dd.local_cells(), dd.ghost_cells());
  BOOST_CHECK_MESSAGE(violations.empty(),
                      "Expected no local-check violations for real "
                      "AtomDecomposition (collective-covered ghosts)");
  BOOST_CHECK(validate_halo_plan_symmetry(*plan).empty());
}

// A real HybridDecomposition covers the regular ghosts via point-to-point and
// the n-square ghosts via the collective section — both paths at once.
BOOST_AUTO_TEST_CASE(hybrid_decomposition_mixed_coverage,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_4_mpi_ranks();
                     })) {
  using namespace GhostComm;
  auto const box = make_periodic_box(8.);
  auto const dd = make_hybrid_dd({2, 2, 1}, box, 1.2, 0.1, std::set<int>{0});
  auto const *plan = dd.halo_plan();
  BOOST_REQUIRE(plan != nullptr);
  auto violations =
      validate_halo_plan(*plan, dd.local_cells(), dd.ghost_cells());
  BOOST_CHECK_MESSAGE(violations.empty(),
                      "Expected no local-check violations for real "
                      "HybridDecomposition (mixed p2p + collective coverage)");
  BOOST_CHECK(validate_halo_plan_symmetry(*plan).empty());
}

// ── Task 2.2: cross-rank symmetry checks ────────────────────────────────────

// Good case: real RegularDecomposition on 4 ranks must be symmetric.
BOOST_AUTO_TEST_CASE(symmetry_good_real_decomposition,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_4_mpi_ranks();
                     })) {
  using namespace GhostComm;
  auto const dd = make_dd({2, 2, 1}, 8., 1.2);
  auto const *plan = dd.halo_plan();
  BOOST_REQUIRE(plan != nullptr);
  auto violations = validate_halo_plan_symmetry(*plan);
  BOOST_CHECK_MESSAGE(violations.empty(),
                      "Expected no symmetry violations for real decomposition");
}

// Asymmetric case: every rank posts a NeighborComm to peer=(me+1)%n with
// send=[1 cell] and recv=[] (expects 0 from peer), but peer sends 1 to me.
// Each rank should detect a mismatch for peer (me-1+n)%n.
BOOST_AUTO_TEST_CASE(symmetry_detects_mismatch,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_4_mpi_ranks();
                     })) {
  using namespace GhostComm;
  boost::mpi::communicator comm;
  int const n = comm.size();
  int const me = comm.rank();

  // Dummy cell for the send region.
  Cell dummy;

  HaloPlan asym;
  asym.comm = comm;
  // Every rank sends 1 item to (me+1)%n but declares recv=0 from that peer.
  // So rank (me+1)%n will receive my send count=1 via all_to_all but see that
  // my declared recv count from (me+1)%n is 0, which doesn't match what
  // (me+1)%n actually sends.  Each rank therefore detects a mismatch.
  int const send_peer = (me + 1) % n;
  asym.neighbors.push_back(
      NeighborComm{send_peer,
                   {SendRegion{&dummy.particles(), {}}}, // send.size() == 1
                   {}});                                 // recv.size() == 0

  auto violations = validate_halo_plan_symmetry(asym);
  BOOST_CHECK_MESSAGE(!violations.empty(),
                      "Expected symmetry violation: recv=0 but peer sends 1");
}
