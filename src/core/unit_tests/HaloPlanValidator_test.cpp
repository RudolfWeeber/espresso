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
#include "cell_system/Cell.hpp"
#include "ghosts/HaloPlan.hpp"
#include "ghosts/HaloPlanValidator.hpp"
#include <boost/test/unit_test.hpp>
#include <span>
#include <vector>

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
