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
/** \file
 *  Asynchronous, split-phase ghost-communication engine. See HaloExchange.hpp.
 */

#include "ghosts/HaloExchange.hpp"

#include "BoxGeometry.hpp"
#include "ghosts.hpp"
#include "ghosts/HaloPlan.hpp"
#include "ghosts/particle_packing.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/nonblocking.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <span>
#include <unordered_set>
#include <vector>

namespace GhostComm {

namespace {

/**
 * @brief MPI tags for ghost communications.
 *
 * One tag per exchange-kind keeps concurrent exchanges (a future goal) from
 * aliasing: within a single exchange each message is already uniquely matched
 * by @c (peer, tag). The BONDS transfer needs its own tag because it is sent
 * as a *second* message to the same peer, mirroring the legacy two-transfer.
 */
enum GhostTag : int {
  TAG_POSITION = 100,
  TAG_FORCE = 101,
  TAG_PARTNUM = 102,
  TAG_DEFAULT = 103,
  TAG_BONDS = 104,
};

int main_tag(unsigned data_parts) {
  if (data_parts & GHOSTTRANS_PARTNUM)
    return TAG_PARTNUM;
  if (data_parts & GHOSTTRANS_FORCE)
    return TAG_FORCE;
  if (data_parts & GHOSTTRANS_POSITION)
    return TAG_POSITION;
  return TAG_DEFAULT;
}
// NOTE: the (peer, tag) pair is what makes each MPI message unambiguous.
// This scheme is safe only when each peer appears at most once in the plan
// (the "one NeighborComm per peer" invariant). See halo_exchange_start @pre.

/** @brief View a list of cell pointers as a span for the packing routines. */
std::span<ParticleList *const> as_span(std::vector<ParticleList *> const &v) {
  return {v.data(), v.size()};
}

/** @brief Extract the plain cell pointers from a list of send regions. */
std::vector<ParticleList *> region_cells(std::vector<SendRegion> const &send) {
  std::vector<ParticleList *> cells;
  cells.reserve(send.size());
  for (auto const &r : send)
    cells.push_back(r.cell);
  return cells;
}

/**
 * @brief Pack a list of send regions into one buffer, honoring per-region
 *        shifts, producing bytes identical to a flat pack over the region
 *        cells in order.
 *
 * @c pack_cells resizes/overwrites the buffer it is given, so each region is
 * packed into a scratch buffer and its bytes appended here. All shifts are
 * currently zero (R5), but the per-region shift field is honored.
 */
void pack_regions(CommBuf &buf, std::vector<SendRegion> const &regions,
                  BoxGeometry const &box, unsigned data_parts) {
  buf.resize(0);
  buf.bonds().clear();

  CommBuf tmp;
  for (auto const &region : regions) {
    std::array<ParticleList *, 1> one{region.cell};
    pack_cells(tmp, std::span<ParticleList *const>{one}, region.shift, box,
               data_parts);
    auto const *first = tmp.data();
    auto const old_size = buf.size();
    buf.resize(old_size + tmp.size());
    std::copy(first, first + tmp.size(), buf.data() + old_size);
    if (data_parts & GHOSTTRANS_BONDS) {
      auto &dst = buf.bonds();
      dst.insert(dst.end(), tmp.bonds().begin(), tmp.bonds().end());
    }
  }
}

/**
 * @brief Collective (broadcast/reduce-sum) section.
 *
 * Task 1.6 will implement real MPI collectives here. For now no plan produces
 * a non-trivial collective section, so we only assert that.
 */
void run_collective(HaloPlan const &plan, BoxGeometry const &,
                    unsigned /*data_parts*/, ExchangeOp /*op*/) {
  // Task 1.6: implement CollectivePattern::Broadcast / ReduceSum here.
  assert(!plan.collective ||
         plan.collective->pattern == CollectivePattern::None);
  static_cast<void>(plan);
}

} // namespace

GhostExchange halo_exchange_start(HaloPlan const &plan, BoxGeometry const &box,
                                  unsigned data_parts, ExchangeOp op) {
  GhostExchange st;
  st.op = op;
  st.data_parts = data_parts;
  st.box = &box;
  st.plan = &plan;
  if (data_parts == GHOSTTRANS_NONE)
    return st;

#ifdef ESPRESSO_ADDITIONAL_CHECKS
  // Invariant: each peer appears at most once in plan.neighbors.
  // Multiple send/recv regions to the same peer must be folded into a single
  // NeighborComm; see halo_exchange_start @pre for why this is required.
  {
    std::unordered_set<int> seen_peers;
    seen_peers.reserve(plan.neighbors.size());
    for (auto const &nc : plan.neighbors) {
      assert(seen_peers.insert(nc.peer).second &&
             "halo_exchange_start: duplicate peer in plan.neighbors — "
             "fold multiple regions to the same peer into one NeighborComm");
    }
  }
#endif

  auto const &comm = plan.comm;
  auto const n = plan.neighbors.size();
  st.send.resize(n);
  st.recv.resize(n);
  st.send_cells.resize(n);
  st.recv_cells.resize(n);

  bool const push = op.direction == Direction::Push;
  bool const bonds = (data_parts & GHOSTTRANS_BONDS) != 0u;
  int const tag = main_tag(data_parts);

  // Resolve, per neighbor, which cells we send from and which we receive into.
  // Push:   send from the region cells, receive into the ghost (recv) cells.
  // Reduce: roles swap -- send from ghost cells, receive into region cells,
  //         adding on arrival.
  for (std::size_t i = 0; i < n; ++i) {
    auto const &nc = plan.neighbors[i];
    if (push) {
      st.send_cells[i] = region_cells(nc.send);
      st.recv_cells[i] = nc.recv;
    } else {
      st.send_cells[i] = nc.recv;
      st.recv_cells[i] = region_cells(nc.send);
    }
  }

  // 1) Post ALL receives first (deadlock-free). Sizes are known a-priori from
  //    the packed size of the cells we will receive into, except PARTNUM whose
  //    ghost cells are not sized yet -- there we post a fixed-size recv and let
  //    unpack_cells resize the cells (legacy bootstrap path).
  for (std::size_t i = 0; i < n; ++i) {
    auto const &nc = plan.neighbors[i];
    st.recv[i].resize(
        calc_transmit_size(as_span(st.recv_cells[i]), box, data_parts));
    st.requests.push_back(comm.irecv(nc.peer, tag, st.recv[i].data(),
                                     static_cast<int>(st.recv[i].size())));
  }

  // 2) Pack, then post ALL sends.
  for (std::size_t i = 0; i < n; ++i) {
    auto const &nc = plan.neighbors[i];
    if (push) {
      pack_regions(st.send[i], nc.send, box, data_parts);
    } else {
      // Reduce: pack plain ghost cells, no shift.
      pack_cells(st.send[i], as_span(st.send_cells[i]), {}, box, data_parts);
    }
    st.requests.push_back(comm.isend(nc.peer, tag, st.send[i].data(),
                                     static_cast<int>(st.send[i].size())));
  }

  // BONDS (cold, resort-only path): a second per-neighbor message for the bond
  // buffers, mirroring the legacy two-transfer. Kept entirely off the hot
  // POSITION/FORCE path via the flag guard. boost::mpi serializes the
  // std::vector<char> length + payload, so recv sizing is automatic here.
  if (bonds) {
    for (std::size_t i = 0; i < n; ++i) {
      auto const &nc = plan.neighbors[i];
      st.requests.push_back(comm.irecv(nc.peer, TAG_BONDS, st.recv[i].bonds()));
    }
    for (std::size_t i = 0; i < n; ++i) {
      auto const &nc = plan.neighbors[i];
      st.requests.push_back(comm.isend(nc.peer, TAG_BONDS, st.send[i].bonds()));
    }
  }

  return st;
}

void halo_exchange_finish(GhostExchange &st) {
  if (st.data_parts == GHOSTTRANS_NONE)
    return;

  // Same-rank copies can run while messages are in flight.
  for (auto const &lc : st.plan->local)
    local_cell_copy(*lc.src, *lc.dst, lc.shift, *st.box, st.data_parts);

  if (st.plan->collective)
    run_collective(*st.plan, *st.box, st.data_parts, st.op);

  boost::mpi::wait_all(st.requests.begin(), st.requests.end());

  // Unpack / reduce into the destination cells resolved at start.
  bool const add = st.op.combine == Combine::Add;
  for (std::size_t i = 0; i < st.plan->neighbors.size(); ++i) {
    auto dst = as_span(st.recv_cells[i]);
    if (add && st.data_parts == GHOSTTRANS_FORCE) {
      add_forces(st.recv[i], dst);
    }
#ifdef ESPRESSO_BOND_CONSTRAINT
    else if (add && st.data_parts == GHOSTTRANS_RATTLE) {
      add_rattle(st.recv[i], dst);
    }
#endif
    else {
      unpack_cells(st.recv[i], dst, *st.box, st.data_parts);
    }
  }
}

void halo_exchange(HaloPlan const &plan, BoxGeometry const &box,
                   unsigned data_parts, ExchangeOp op) {
  auto st = halo_exchange_start(plan, box, data_parts, op);
  halo_exchange_finish(st);
}

} // namespace GhostComm
