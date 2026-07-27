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
/** \file
 *  Local (single-rank) validator for GhostComm::HaloPlan.
 *  See HaloPlanValidator.hpp.
 */

#include "ghosts/HaloPlanValidator.hpp"

#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace GhostComm {

std::vector<std::string>
validate_halo_plan(HaloPlan const &plan, std::span<Cell *const> local_cells,
                   std::span<Cell *const> ghost_cells) {
  std::vector<std::string> violations;

  // Build the set of expected ghost ParticleList pointers.
  std::unordered_set<ParticleList const *> ghost_set;
  for (Cell *c : ghost_cells) {
    ghost_set.insert(&c->particles());
  }

  // Check peer-uniqueness and shape; accumulate recv/dst fill counts.
  std::unordered_map<ParticleList const *, int> fill_count;
  std::unordered_set<int> seen_peers;

  for (auto const &nc : plan.neighbors) {
    // Peer-uniqueness
    if (!seen_peers.insert(nc.peer).second) {
      std::ostringstream oss;
      oss << "peer " << nc.peer << " appears in more than one NeighborComm";
      violations.push_back(oss.str());
    }

    // Shape: send.size() == recv.size()
    if (nc.send.size() != nc.recv.size()) {
      std::ostringstream oss;
      oss << "NeighborComm peer=" << nc.peer
          << " has send.size()=" << nc.send.size()
          << " != recv.size()=" << nc.recv.size();
      violations.push_back(oss.str());
    }

    // Accumulate recv fill counts; check targets are in ghost set.
    for (ParticleList const *pl : nc.recv) {
      if (ghost_set.find(pl) == ghost_set.end()) {
        std::ostringstream oss;
        oss << "NeighborComm peer=" << nc.peer
            << " recv target is not a ghost cell";
        violations.push_back(oss.str());
      }
      ++fill_count[pl];
    }
  }

  // Accumulate local.dst fill counts; check targets are in ghost set.
  for (auto const &lc : plan.local) {
    ParticleList const *pl = lc.dst;
    if (ghost_set.find(pl) == ghost_set.end()) {
      violations.push_back("LocalComm dst target is not a ghost cell");
    }
    ++fill_count[pl];
  }

  // Coverage: every ghost must appear exactly once.
  for (Cell *c : ghost_cells) {
    ParticleList const *pl = &c->particles();
    auto it = fill_count.find(pl);
    int count = (it != fill_count.end()) ? it->second : 0;
    if (count == 0) {
      violations.push_back("ghost cell is never filled (missing recv/dst)");
    } else if (count > 1) {
      std::ostringstream oss;
      oss << "ghost cell is filled " << count << " times (expected 1)";
      violations.push_back(oss.str());
    }
  }

  // Neighborship-match: every local cell's ghost neighbor must be covered.
  for (Cell *c : local_cells) {
    for (Cell *n : c->neighbors().all()) {
      ParticleList const *pl = &n->particles();
      if (ghost_set.find(pl) == ghost_set.end()) {
        continue; // not a ghost neighbor
      }
      auto it = fill_count.find(pl);
      int count = (it != fill_count.end()) ? it->second : 0;
      if (count == 0) {
        violations.push_back(
            "local cell has ghost neighbor that is not a covered recv/dst "
            "target (referenced-but-uncommunicated ghost)");
      }
    }
  }

  return violations;
}

} // namespace GhostComm
