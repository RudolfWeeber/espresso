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
#include <functional>
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
 * @brief Pack all send regions of one NeighborComm into a single buffer.
 *
 * All regions are packed with a **single** @c pack_cells call so that bonds
 * produce exactly one boost binary archive per neighbor message -- matching
 * the single @c unpack_cells call on the receive side.  Concatenating
 * per-region archives would corrupt the bond stream (the second archive
 * header is misread as bond data by the receiver).
 *
 * Per-region shifts are honored via the common shift. All @c SendRegion.shift
 * values within a @c NeighborComm are equal — @c RegularDecomposition::
 * make_halo_plan() sets them all to @c {}. This lets us pack every region's
 * cells in one @c pack_cells call (a single bond archive, matching the
 * receiver's single @c unpack_cells).
 *
 * @note If a future decomposition (e.g. a Lees-Edwards refactor) ever stores
 * DISTINCT per-region shifts, @c pack_regions must be generalized to pack
 * each region's non-bond (flat) data with its own shift while still writing
 * all bonds into ONE shared archive — do not simply concatenate per-region
 * archives (that was the bond-corruption bug fixed in commit 314e7c366b).
 */
void pack_regions(CommBuf &buf, std::vector<SendRegion> const &regions,
                  BoxGeometry const &box, unsigned data_parts) {
#ifdef ESPRESSO_ADDITIONAL_CHECKS
  // All regions in a NeighborComm must share the same shift.  If per-region
  // shifts ever diverge the packer must be generalized (separate archives).
  if (!regions.empty()) {
    auto const &ref = regions.front().shift;
    for (auto const &r : regions) {
      assert(r.shift == ref &&
             "pack_regions: per-region shifts differ within one NeighborComm "
             "-- packer must be generalized for per-region shift support");
    }
  }
#endif

  Utils::Vector3d const common_shift =
      regions.empty() ? Utils::Vector3d{} : regions.front().shift;

  std::vector<ParticleList *> cells;
  cells.reserve(regions.size());
  for (auto const &r : regions)
    cells.push_back(r.cell);

  pack_cells(buf, as_span(cells), common_shift, box, data_parts);
}

/**
 * @brief Collective (broadcast/reduce-sum) section.
 *
 * Mirrors the legacy GHOST_BCST / GHOST_RDCE loop from ghost_communicator()
 * in ghosts.cpp.  For each root rank in [0, comm.size()):
 *
 *   Broadcast (Push): root packs its owned cell (cells[root]) and broadcasts
 *   the buffer to all ranks; every non-root rank unpacks into cells[root]
 *   (their ghost copy of root's data).
 *
 *   ReduceSum (Reduce): every rank packs cells[root] (ghost-force contribution)
 *   and reduces to root with std::plus on the raw double buffer; root unpacks
 *   into cells[root] (the owned cell, so forces accumulate there).
 *
 * The legacy code uses one GHOST_BCST (or GHOST_RDCE) step per rank, each
 * step addressing a single cell pointer (ghost_comm.part_lists[0] ==
 * &cells.at(n)).  We replicate that loop exactly.
 */
void run_collective(HaloPlan const &plan, BoxGeometry const &box,
                    unsigned data_parts, ExchangeOp op) {
  // Precondition: for non-PARTNUM data parts, ghost cell sizes must already be
  // synced by a prior GHOSTTRANS_PARTNUM exchange (same invariant as the legacy
  // GHOST_BCST/GHOST_RDCE path). The per-root broadcast/reduce byte count is
  // derived from the (already-sized) cells; a size mismatch across ranks would
  // be undefined behavior.

  if (!plan.collective || plan.collective->pattern == CollectivePattern::None) {
    return;
  }

  auto const &cs = *plan.collective;
  auto const &comm = plan.comm;
  int const comm_size = comm.size();
  int const my_rank = comm.rank();

  // Use op.direction to select the actual MPI collective: Push -> broadcast
  // (particles flow from owner to all ghosts); Reduce -> reduce-sum (ghost
  // forces flow back to owner).  cs.pattern == Broadcast marks the section as
  // active; the engine caller sets op correctly for each exchange.
  bool const is_broadcast = (op.direction == Direction::Push);

  // cs.cells has one entry per rank: cs.cells[root] is the ParticleList for
  // that root rank (owned on root, ghost on all others).
  assert(static_cast<int>(cs.cells.size()) == comm_size);

  CommBuf buf;

  for (int root = 0; root < comm_size; ++root) {
    ParticleList *cell = cs.cells[static_cast<std::size_t>(root)];
    auto const cell_span = std::span<ParticleList *const>{&cell, 1};

    if (is_broadcast) {
      // Push: root broadcasts its owned particles to all other ranks.
      if (my_rank == root) {
        pack_cells(buf, cell_span, {}, box, data_parts);
        boost::mpi::broadcast(comm, buf.data(), static_cast<int>(buf.size()),
                              root);
        boost::mpi::broadcast(comm, buf.bonds(), root);
      } else {
        buf.resize(calc_transmit_size(cell_span, box, data_parts));
        buf.bonds().clear();
        boost::mpi::broadcast(comm, buf.data(), static_cast<int>(buf.size()),
                              root);
        boost::mpi::broadcast(comm, buf.bonds(), root);
        unpack_cells(buf, cell_span, box, data_parts);
      }
    } else {
      // ReduceSum: every rank sends its ghost-force contribution; root
      // accumulates into its owned cell.  Mirrors the legacy GHOST_RDCE
      // which reduces the raw double buffer with std::plus<double>.
      pack_cells(buf, cell_span, {}, box, data_parts);
      auto *raw = reinterpret_cast<double *>(buf.data());
      int const count = static_cast<int>(buf.size() / sizeof(double));
      if (my_rank == root) {
        CommBuf recv_buf;
        recv_buf.resize(buf.size());
        auto *recv_raw = reinterpret_cast<double *>(recv_buf.data());
        boost::mpi::reduce(comm, raw, count, recv_raw, std::plus<double>{},
                           root);
        unpack_cells(recv_buf, cell_span, box, data_parts);
      } else {
        boost::mpi::reduce(comm, raw, count, std::plus<double>{}, root);
      }
    }
  }
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
