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
#pragma once

/** \file
 *  Asynchronous, split-phase ghost-communication engine.
 *
 *  Drives a @ref GhostComm::HaloPlan with non-blocking point-to-point MPI
 *  (plus same-rank copies), reusing the byte-identical serialization from
 *  particle_packing.hpp. The engine is split into a @ref halo_exchange_start
 *  phase (post all receives, pack and post all sends) and a
 *  @ref halo_exchange_finish phase (same-rank copies, wait, unpack/reduce) so
 *  that callers can overlap other work with in-flight messages.
 *
 *  Deadlock-freedom is guaranteed by posting every @c irecv before any
 *  @c isend and matching each message by @c (peer, tag); message sizes are
 *  known a-priori from the receiving cells' packed size (see
 *  @ref GhostComm::calc_transmit_size), except for the PARTNUM bootstrap.
 */

#include "BoxGeometry.hpp"
#include "ghosts/HaloPlan.hpp"
#include "ghosts/particle_packing.hpp"

#include <boost/mpi/request.hpp>

#include <vector>

namespace GhostComm {

/**
 * @brief Opaque handle for one in-flight halo exchange.
 *
 * Owns the per-neighbor send/recv buffers and the outstanding MPI requests, so
 * that N neighbor messages can be in flight simultaneously (no function-static
 * buffers). Created by @ref halo_exchange_start and consumed exactly once by
 * @ref halo_exchange_finish.
 */
struct GhostExchange {
  ExchangeOp op{};
  unsigned data_parts = 0u;
  BoxGeometry const *box = nullptr;
  HaloPlan const *plan = nullptr;

  /** Per-neighbor packed send buffers (index-aligned with plan->neighbors). */
  std::vector<CommBuf> send;
  /** Per-neighbor recv buffers (index-aligned with plan->neighbors). */
  std::vector<CommBuf> recv;
  /** Outstanding non-blocking send/recv requests. */
  std::vector<boost::mpi::request> requests;

  /**
   * Scratch storage backing the @c std::span views handed to the packing
   * routines. For the Reduce direction the send/recv roles swap, so we need
   * the plain cell pointers extracted from the @ref SendRegion list; these
   * vectors must outlive the pack/unpack calls, hence they live on the handle.
   */
  std::vector<std::vector<ParticleList *>> send_cells;
  std::vector<std::vector<ParticleList *>> recv_cells;
};

/**
 * @brief Begin a halo exchange: post all receives, then pack and post all
 *        sends. Returns immediately with an in-flight @ref GhostExchange.
 *
 * @param plan        Communication plan (peers, regions, local copies).
 * @param box         Box geometry for position folding.
 * @param data_parts  Bitmask of GHOSTTRANS_* flags to transfer.
 * @param op          Direction (Push/Reduce) and combine mode (Overwrite/Add).
 *
 * @pre Each peer rank appears at most once in @p plan.neighbors. Multiple
 *      send/receive regions to the same peer **must** be folded into that
 *      peer's single @ref NeighborComm (as the plan builders do). This is
 *      what makes the @c (peer, data-part tag) message matching unambiguous:
 *      without this invariant two @c NeighborComm entries sharing the same
 *      peer and tag would cross-match, silently corrupting data or deadlocking.
 *      The plan builders (e.g. RegularDecomposition) always produce unique
 *      peers; the invariant is asserted at runtime when
 *      @c ESPRESSO_ADDITIONAL_CHECKS is defined.
 */
GhostExchange halo_exchange_start(HaloPlan const &plan, BoxGeometry const &box,
                                  unsigned data_parts, ExchangeOp op);

/**
 * @brief Complete a halo exchange: run same-rank copies (overlapping the
 *        in-flight messages), wait for all requests, then unpack or reduce.
 */
void halo_exchange_finish(GhostExchange &state);

/**
 * @brief Blocking convenience wrapper: start + finish.
 */
void halo_exchange(HaloPlan const &plan, BoxGeometry const &box,
                   unsigned data_parts, ExchangeOp op);

} // namespace GhostComm
