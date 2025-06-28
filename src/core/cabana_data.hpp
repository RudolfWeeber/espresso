/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#ifdef SHARED_MEMORY_PARALLELISM

#include "custom_verlet_list.hpp"
#include <Cabana_Core.hpp>
#include <unordered_map>
#include <vector>

using memory_space = Kokkos::SharedSpace;
using execution_space = Kokkos::DefaultExecutionSpace;

using ListAlgorithm = Cabana::HalfNeighborTag;
using ListType = Cabana::CustomVerletList<memory_space, ListAlgorithm,
                                          Cabana::VerletLayout2D>;

class CabanaData {
private:
  ListType verlet_list;
  std::vector<Particle *> unique_particles;
  int max_id;

public:
  CabanaData() = default;
  CabanaData(ListType verlet_list, std::vector<Particle *> unique_particles,
             int max_id)
      : verlet_list(verlet_list), unique_particles(unique_particles),
        max_id(max_id) {}

  ListType get_verlet_list() const { return verlet_list; }
  int get_index() const { return unique_particles.size(); }
  int get_max_id() const { return max_id; }
  std::vector<Particle *> get_unique_particles() const {
    return unique_particles;
  }

  ~CabanaData() {};
};
#endif
