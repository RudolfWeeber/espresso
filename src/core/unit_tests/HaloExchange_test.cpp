/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE HaloExchange test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/mpi.hpp>

#include "BondList.hpp"
#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "ParticleList.hpp"
#include "ghosts.hpp"
#include "ghosts/HaloExchange.hpp"
#include "ghosts/HaloPlan.hpp"

#include <array>

namespace utf = boost::unit_test;

/*
 * On two ranks, each rank owns a single particle tagged by its rank. A Push
 * of GHOSTTRANS_POSITION must overwrite the local ghost cell with the peer's
 * position. This exercises the deadlock-free irecv/isend protocol and the
 * a-priori recv sizing.
 */
BOOST_AUTO_TEST_CASE(push_positions_between_two_ranks,
                     *utf::precondition([](utf::test_unit_id) {
                       return boost::mpi::communicator{}.size() == 2;
                     })) {
  using namespace GhostComm;
  boost::mpi::communicator world;
  int const me = world.rank();
  int const other = 1 - me;

  BoxGeometry box;
  box.set_length({10., 10., 10.});

  ParticleList local, ghost;
  local.resize(1);
  local.begin()->id() = me; // tag by owner rank
  local.begin()->pos() = {double(me), 2. * double(me), 3. * double(me)};
  ghost.resize(1); // sized as if ghosts_count() had already run

  HaloPlan plan;
  plan.comm = world;
  plan.neighbors.push_back(
      NeighborComm{other, {SendRegion{&local, {}}}, {&ghost}});

  halo_exchange(plan, box, GHOSTTRANS_POSITION,
                {Direction::Push, Combine::Overwrite});

  // ghost cell now holds the peer's (folded) position.
  BOOST_CHECK_CLOSE(ghost.begin()->pos()[0], double(other), 1e-12);
  BOOST_CHECK_CLOSE(ghost.begin()->pos()[1], 2. * double(other), 1e-12);
  BOOST_CHECK_CLOSE(ghost.begin()->pos()[2], 3. * double(other), 1e-12);
  // the local (owned) particle must be untouched by a Push.
  BOOST_CHECK_CLOSE(local.begin()->pos()[0], double(me), 1e-12);
}

/*
 * A Reduce with Combine::Add of GHOSTTRANS_FORCE folds ghost forces back onto
 * the owner. The send/recv roles swap: we send from the *ghost* cell and add
 * into the *owned* cell. Each rank seeds a distinct force on its ghost, so the
 * owner's force must end up as (its own pre-existing force + the peer's ghost
 * force).
 */
BOOST_AUTO_TEST_CASE(reduce_forces_between_two_ranks,
                     *utf::precondition([](utf::test_unit_id) {
                       return boost::mpi::communicator{}.size() == 2;
                     })) {
  using namespace GhostComm;
  boost::mpi::communicator world;
  int const me = world.rank();
  int const other = 1 - me;

  BoxGeometry box;
  box.set_length({10., 10., 10.});

  ParticleList local, ghost;
  local.resize(1);
  ghost.resize(1);
  // Pre-existing force on the owned particle; Reduce must add to it.
  local.begin()->force() = {1., 0., 0.};
  // Force accumulated on this rank's ghost (from short-range interactions).
  // Encode the ghost's owner (== other) so the sum is checkable.
  ghost.begin()->force() = {0., 10. * double(me + 1), 100. * double(me + 1)};

  HaloPlan plan;
  plan.comm = world;
  plan.neighbors.push_back(
      NeighborComm{other, {SendRegion{&local, {}}}, {&ghost}});

  halo_exchange(plan, box, GHOSTTRANS_FORCE, {Direction::Reduce, Combine::Add});

  // Owner force == own pre-existing force + peer's ghost force.
  BOOST_CHECK_CLOSE(local.begin()->force()[0], 1., 1e-12);
  BOOST_CHECK_CLOSE(local.begin()->force()[1], 10. * double(other + 1), 1e-12);
  BOOST_CHECK_CLOSE(local.begin()->force()[2], 100. * double(other + 1), 1e-12);
}

/*
 * The cold, resort-only BONDS path is exchanged as a *second* per-neighbor
 * message (the CommBuf bond buffer) alongside the flat particle buffer. This
 * checks that a bond placed on the owner round-trips byte-correct into the
 * peer's ghost cell over the async engine.
 */
BOOST_AUTO_TEST_CASE(push_bonds_between_two_ranks,
                     *utf::precondition([](utf::test_unit_id) {
                       return boost::mpi::communicator{}.size() == 2;
                     })) {
  using namespace GhostComm;
  boost::mpi::communicator world;
  int const me = world.rank();
  int const other = 1 - me;

  BoxGeometry box;
  box.set_length({10., 10., 10.});

  ParticleList local, ghost;
  local.resize(1);
  local.begin()->id() = me;
  std::array<int, 2> partners{me, me + 5};
  local.begin()->bonds().insert(BondView{me + 1, partners});
  ghost.resize(1);

  HaloPlan plan;
  plan.comm = world;
  plan.neighbors.push_back(
      NeighborComm{other, {SendRegion{&local, {}}}, {&ghost}});

  halo_exchange(plan, box, GHOSTTRANS_PROPRTS | GHOSTTRANS_BONDS,
                {Direction::Push, Combine::Overwrite});

  BOOST_REQUIRE_EQUAL(ghost.begin()->bonds().size(), 1);
  auto it = ghost.begin()->bonds().begin();
  BOOST_CHECK_EQUAL((*it).bond_id(), other + 1);
  BOOST_REQUIRE_EQUAL((*it).partner_ids().size(), 2u);
  BOOST_CHECK_EQUAL((*it).partner_ids()[0], other);
  BOOST_CHECK_EQUAL((*it).partner_ids()[1], other + 5);
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
