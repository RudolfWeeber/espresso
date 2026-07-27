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
 * 1. Coverage: every ghost cell appears as a recv/local.dst target exactly
 *    once; no target lies outside the ghost set; none appears twice.
 * 2. Neighborship-match: every local cell's ghost neighbor is a covered
 *    recv/dst target.
 * 3. Peer-uniqueness: each peer value appears in at most one NeighborComm.
 * 4. Shape: each NeighborComm has send.size() == recv.size().
 *
 * @returns a human-readable violation string per problem; empty = valid.
 */
std::vector<std::string> validate_halo_plan(HaloPlan const &plan,
                                            std::span<Cell *const> local_cells,
                                            std::span<Cell *const> ghost_cells);

} // namespace GhostComm
