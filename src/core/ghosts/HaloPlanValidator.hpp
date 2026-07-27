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
#pragma once

#include "cell_system/Cell.hpp"
#include "ghosts/HaloPlan.hpp"

#include <span>
#include <string>
#include <vector>

namespace GhostComm {

/**
 * @brief Validate a HaloPlan for correctness.
 *
 * Checks:
 * 1. Coverage: every ghost cell is filled. Point-to-point ghosts appear as a
 *    recv/local.dst target exactly once; no p2p target lies outside the ghost
 *    set; none appears twice. Ghosts filled by the collective broadcast/reduce
 *    section (AtomDecomposition, HybridDecomposition) are covered by that
 *    section instead of a recv/dst target.
 * 2. Neighborship-match: every local cell's ghost neighbor is covered, either
 *    as a recv/dst target or by the collective section.
 * 3. Peer-uniqueness: each peer value appears in at most one NeighborComm.
 * 4. Shape: each NeighborComm has send.size() == recv.size().
 *
 * @returns a human-readable violation string per problem; empty = valid.
 */
std::vector<std::string> validate_halo_plan(HaloPlan const &plan,
                                            std::span<Cell *const> local_cells,
                                            std::span<Cell *const> ghost_cells);

/**
 * @brief Cross-rank symmetry check for a HaloPlan.
 *
 * Uses a collective all-to-all exchange of per-rank send-counts so that every
 * rank participates regardless of whether the neighbor sets are symmetric.
 * This avoids the deadlock that a naive per-neighbor isend/irecv pattern
 * would cause when a rank sends to a peer that has no matching recv posted.
 *
 * Invariant: for every peer P, my recv-count from P must equal P's send-count
 * to me (i.e. my_recv_from[P] == peers_send_to_me[P]).
 *
 * @returns a human-readable violation string per mismatch; empty = valid.
 */
std::vector<std::string> validate_halo_plan_symmetry(HaloPlan const &plan);

} // namespace GhostComm
