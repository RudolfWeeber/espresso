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

#include "config/config.hpp"

#include "cell_system/CellStructure.hpp"
#include "lees_edwards/lees_edwards.hpp"

#ifdef CALIPER
#include <caliper/cali.h>
#endif

#ifdef SHARED_MEMORY_PARALLELISM

#include "aosoa_pack.hpp"
#include "cabana_data.hpp"
#include "custom_verlet_list.hpp"
#include <Cabana_Core.hpp>
#include <Cabana_NeighborList.hpp>
#include <cassert>
#include <iostream>
#include <stdio.h>
#include <unordered_set>
#include <utility>

using memory_space = Kokkos::HostSpace; // Kokkos::SharedSpace;
using execution_space = Kokkos::DefaultExecutionSpace;

inline void set_offset_and_size_indexed_by_cid(
    int &total_bins, int *cell_num,
    Cabana::LinkedCellList<memory_space, double> &cell_list,
    Kokkos::View<int *, Kokkos::LayoutRight> &bin_offset,
    Kokkos::View<int *, Kokkos::LayoutRight> &bin_size,
    Kokkos::View<int *, Kokkos::LayoutRight> &original_idx) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  // for (int cid = 0; cid < total_bins; ++cid) {
  Kokkos::parallel_for(
      "set_offset", total_bins,
      [&cell_num, &cell_list, &bin_offset, &bin_size](const int cid) {
        int dx[3] = {};
        dx[0] = static_cast<int>(cid / (cell_num[1] * cell_num[2]));
        dx[1] = static_cast<int>((cid - dx[0] * (cell_num[1] * cell_num[2])) /
                                 cell_num[2]);
        dx[2] = cid % cell_num[2];
        bin_offset(cid) = cell_list.binOffset(dx[0], dx[1], dx[2]);
        bin_size(cid) = cell_list.binSize(dx[0], dx[1], dx[2]);

        // Calculate particle_bins
        cell_list(cid);
      });
  Kokkos::parallel_for("set_permutation", original_idx.extent(0),
                       [&cell_list, &original_idx](const int i) {
                         original_idx(i) = cell_list.permutation(i);
                       });
}

using ActiveProtocol = std::variant<LeesEdwards::Off, LeesEdwards::LinearShear,
                                    LeesEdwards::OscillatoryShear>;
inline int set_interacting_pair_cell(
    int &total_bins, int total_pair_cell, int *cell_num, int *delta_lebc,
    int le_direction, int le_normal,
    std::shared_ptr<ActiveProtocol> le_protocol,
    // ActiveProtocol le_protocol,
    Kokkos::View<int *, Kokkos::LayoutRight> &bin_size,
    Cabana::LinkedCellList<memory_space, double> &cell_list,
    //std::vector<std::pair<int, int>> &interacting_pair_cell) {
    Kokkos::View<int **, Kokkos::LayoutRight> &interacting_pair_cell) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  constexpr int ijkIndexes[27][3] = {
      {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0},
      {-1, 0, 1},   {-1, 1, -1}, {-1, 1, 0},  {-1, 1, 1},  {0, -1, -1},
      {0, -1, 0},   {0, -1, 1},  {0, 0, -1},  {0, 0, 0},   {0, 0, 1},
      {0, 1, -1},   {0, 1, 0},   {0, 1, 1},   {1, -1, -1}, {1, -1, 0},
      {1, -1, 1},   {1, 0, -1},  {1, 0, 0},   {1, 0, 1},   {1, 1, -1},
      {1, 1, 0},    {1, 1, 1}};

  int empty_pair_number = 0;
  int pair_cell_id = 0;

  for (int cid_i = 0; cid_i < total_bins; ++cid_i) {
    // Kokkos::parallel_for("set_interacting_pair_cell", total_bins,
    //		  KOKKOS_LAMBDA(const int cid_i) {
    // auto thread_id = omp_get_thread_num();
    //  Obtaining 3 dimentional cell index from cid_i
    int index[3] = {};
    cell_list.ijkBinIndex(cid_i, index[0], index[1], index[2]);
    int dx[3];
    // From 27 neighbor cell, the list of interacting pair cell is created
    for (int n = 0; n < 27; ++n) {

      if (le_protocol == nullptr) {
        if (index[0] != 0 and ijkIndexes[n][0] == -1)
          continue;

        if (index[1] != 0 and ijkIndexes[n][1] == -1 and ijkIndexes[n][0] == 0)
          continue;
      }

      bool duplicate_cell = false;
      // Obtaining 3 dimentional cell index from neighbor cell
      for (int d = 0; d < 3; ++d) {
        dx[d] = (ijkIndexes[n][d] + index[d] + cell_num[d]) % cell_num[d];
        if (cell_num[d] <= 2 and ijkIndexes[n][d] + index[d] != dx[d])
          duplicate_cell = true;
      }
      if (duplicate_cell)
        continue;

      // Lees-Edwards BC
      int le_crossing = 0;
      if (le_protocol != nullptr) {
        le_crossing =
            ijkIndexes[n][le_normal] + index[le_normal] - dx[le_normal];
        if (le_crossing < 0) {
          dx[le_direction] = (dx[le_direction] + delta_lebc[le_direction] +
                              cell_num[le_direction]) %
                             cell_num[le_direction];
        } else if (le_crossing > 0) {
          dx[le_direction] = (dx[le_direction] - delta_lebc[le_direction] +
                              cell_num[le_direction]) %
                             cell_num[le_direction];
        }
      }
      // Additional Cell
      /*
      if (le_crossing != 0 && index[le_direction] == 1) {
        if (le_crossing < 0) {
          dx[le_direction] = (dx[le_direction] + 1 +
                  cell_num[le_direction]) % cell_num[le_direction];
        } else if (le_crossing > 0) {
          dx[le_direction] = (dx[le_direction] - 1 +
            cell_num[le_direction]) % cell_num[le_direction];
        }
        cell_offset = bin_offset(dx[0], dx[1], dx[2]);
        cell_size = bin_size(dx[0], dx[1], dx[2]);
      }
      */

      // Interacting pair cell is registered in the list
      int cid_j = cell_list.cardinalBinIndex(dx[0], dx[1], dx[2]);
      // if (cid_i <= cid_j) {
      if (cid_i < cid_j) {
        if (bin_size(cid_i) != 0 and bin_size(cid_j) != 0) {
          // std::size_t pcid = Kokkos::atomic_fetch_inc(&pair_cell_id());
          // interacting_pair_cell(pcid, 0) = cid_i;
          // interacting_pair_cell(pcid, 1) = cid_j;
          interacting_pair_cell(pair_cell_id, 0) = cid_i;
          interacting_pair_cell(pair_cell_id, 1) = cid_j;
          // interacting_pair_cell.emplace_back(std::pair(cid_i, cid_j));
          ++pair_cell_id;
          // interacting_pair_cell_thread(thread_id, pair_id_thread(thread_id),
          // 0) = cid_i; interacting_pair_cell_thread(thread_id,
          // pair_id_thread(thread_id), 1) = cid_j; pair_id_thread(thread_id) +=
          // 1;
        } else {
          // Kokkos::atomic_inc(&empty_pair_number());
          ++empty_pair_number;
          // empty_thread(thread_id) += 1;
        }
      }
    }
    //});
  }
  return empty_pair_number;
}

using ListAlgorithm = Cabana::HalfNeighborTag;
using ListType = Cabana::CustomVerletList<memory_space, ListAlgorithm,
                                          Cabana::VerletLayout2D>;
template <class VerletCriterion, class Kernel>
ListType create_verlet_list(double const max_cutoff, int const max_counts,
                            AoSoA_pack &aosoa,
                            std::vector<Particle *> &unique_particles,
                            VerletCriterion const &verlet_criterion,
                            Kernel &first_neighbor_kernel,
                            CellStructure &cell_structure) {
  // Creating LinkedCellList and VerletList:
  // Box Properties
  auto const &system = ::System::get_system();
  auto box_geo = *(system.box_geo);
  auto box_l = box_geo.length();
  Cabana::LinkedCellList<memory_space, double> cell_list;
  double grid_min[3] = {0.0, 0.0, 0.0};
  double grid_max[3] = {box_l[0], box_l[1], box_l[2]};
  double grid_delta[3] = {};
  int cell_num[3] = {};
  double eff_cutoff;
  for (int d = 0; d < 3; ++d) {
    eff_cutoff = max_cutoff;
    if (eff_cutoff > box_l[d])
      eff_cutoff = box_l[d];
    cell_num[d] = static_cast<int>(box_l[d] / eff_cutoff);
    grid_delta[d] = std::nextafter(box_l[d] / cell_num[d], 0);
  }
  // For Lees-Edwards boundary condition
  double le_offset;
  int le_direction;
  int le_normal;
  int delta_lebc[3] = {0, 0, 0};
  auto le_protocol = system.lees_edwards->get_protocol();
  if (le_protocol == nullptr) {
    le_offset = 0.;
    le_direction = -1;
    le_normal = -1;
  } else {
    le_offset = box_geo.lees_edwards_bc().pos_offset;
    le_direction = box_geo.lees_edwards_bc().shear_direction;
    le_normal = box_geo.lees_edwards_bc().shear_plane_normal;
    delta_lebc[le_direction] =
        static_cast<int>(std::ceil(le_offset / grid_delta[le_direction])) %
        cell_num[le_direction];
  }
#ifdef CALIPER
  CALI_MARK_BEGIN("Cabana - CellList");
#endif
  cell_list = Cabana::createLinkedCellList<memory_space>(
      aosoa.position, grid_delta, grid_min, grid_max);
#ifdef CALIPER
  CALI_MARK_END("Cabana - CellList");
#endif
  int total_bins = cell_list.totalBins();
  // Now permute the AoSoA (i.e. reorder the data) using the linked cell
  // list.
  // Cabana::permute( cell_list, particle_storage );

  // Number of threads
  //int num_threads = execution_space().concurrency();
  ListType verlet_list = ListType(aosoa.position, 0, aosoa.position.size(),
                                  max_counts);//, num_threads);

  // Offset particle id and the number of particle in specific cell
  Kokkos::View<int *, Kokkos::LayoutRight> bin_offset("bin_offset", total_bins);
  Kokkos::View<int *, Kokkos::LayoutRight> bin_size("bin_size", total_bins);
  Kokkos::View<int *, Kokkos::LayoutRight> original_idx("original_idx",
                                                        aosoa.position.size());
  set_offset_and_size_indexed_by_cid(total_bins, cell_num, cell_list,
                                     bin_offset, bin_size, original_idx);
  auto const particle_bins = cell_list.getParticleBins();

  // Creating Interacting cell
  int total_pair_cell;
  if (total_bins < 27) {
    total_pair_cell = (total_bins - 1) * total_bins / 2 + total_bins;
  } else {
    total_pair_cell = 14 * total_bins;
  }
  Kokkos::View<int **, Kokkos::LayoutRight> interacting_pair_cell(
      "interacting_pair_cell", total_pair_cell - total_bins, 2);
  //std::vector<std::pair<int, int>> interacting_pair_cell;
  //interacting_pair_cell.reserve(total_pair_cell - total_bins);
  int empty_pair_number = set_interacting_pair_cell(
      total_bins, total_pair_cell, cell_num, delta_lebc, le_direction,
      le_normal, le_protocol, bin_size, cell_list, interacting_pair_cell);

  auto const distance_function = detail::MinimalImageDistance{
      std::as_const(cell_structure).decomposition().box()};

  auto aosoa_id = aosoa.id;
  auto aosoa_ghost = aosoa.ghost;

  // This kernel calculate within each cell
  auto kernel_each = [&bin_offset, &bin_size, &original_idx, &aosoa_id,
                      &aosoa_ghost, &unique_particles, &verlet_criterion,
                      &distance_function, &verlet_list](const int cid_i) {
                      //&first_neighbor_kernel](const int cid_i) {
    //auto thread_id = omp_get_thread_num();

    int offset_i = bin_offset(cid_i);
    int size_i = bin_size(cid_i);

    for (int i = offset_i; i < offset_i + size_i; ++i) {
      // int ii = i;
      int ii = original_idx(i); // get previous id
      int id_i = aosoa_id(ii);
      auto p1 = unique_particles.at(ii);
      //  auto p1 = cell_structure.get_local_particle(id_i);
      for (int j = i + 1; j < offset_i + size_i; ++j) {
        // int jj = j;
        int jj = original_idx(j);
        int id_j = aosoa_id(jj);
        if (aosoa_ghost(ii) or aosoa_ghost(jj)) {
          if (((id_i < id_j) and aosoa_ghost(ii)) or
              ((id_i > id_j) and aosoa_ghost(jj))) {
            continue;
          }
        }
        // auto p2 = unique_particles.at(jj);
        //  auto p2 = cell_structure.get_local_particle(id_j);
        if (verlet_criterion(*p1, *unique_particles.at(jj),
                             distance_function(*p1,
                                               *unique_particles.at(jj)))) {
          //verlet_list.addNeighborNonAtomic(thread_id, ii, jj);
          verlet_list.addNeighborNonAtomic(ii, jj);
          //first_neighbor_kernel(ii, jj);
        }
      }
    } // i-loop
  };

  // This kernel used the loop for the pair of interacting cell
  auto kernel_neighbor = [&interacting_pair_cell, &bin_offset, &bin_size,
                          &original_idx, &aosoa_id, &aosoa_ghost,
                          &unique_particles, &verlet_criterion,
                          &distance_function, &verlet_list](const int pair_cell_i) {
                          //&first_neighbor_kernel](const int pair_cell_i) {
    int cid_i = interacting_pair_cell(pair_cell_i, 0);
    int cid_j = interacting_pair_cell(pair_cell_i, 1);
    // int cid_i = interacting_pair_cell.at(pair_cell_i).first;
    // int cid_j = interacting_pair_cell.at(pair_cell_i).second;

    //auto thread_id = omp_get_thread_num();

    int offset_i = bin_offset(cid_i);
    int size_i = bin_size(cid_i);
    int offset_j = bin_offset(cid_j);
    int size_j = bin_size(cid_j);

    for (int i = offset_i; i < offset_i + size_i; ++i) {
      // int ii = i;
      int ii = original_idx(i); // get previous id
      int id_i = aosoa_id(ii);
      auto p1 = unique_particles.at(ii);
      // auto p1 = cell_structure.get_local_particle(id_i);

      for (int j = offset_j; j < offset_j + size_j; ++j) {
        // int jj = j;
        int jj = original_idx(j);
        int id_j = aosoa_id(jj);
        if (aosoa_ghost(ii) or aosoa_ghost(jj)) {
          if (((id_i < id_j) and aosoa_ghost(ii)) or
              ((id_i > id_j) and aosoa_ghost(jj))) {
            continue;
          }
        }
        // auto p2 = unique_particles.at(jj);
        // auto p2 = cell_structure.get_local_particle(id_j);
        if (verlet_criterion(*p1, *unique_particles.at(jj),
                             distance_function(*p1,
                                               *unique_particles.at(jj)))) {
          //verlet_list.addNeighbor(thread_id, ii, jj);
          verlet_list.addNeighbor(ii, jj);
          //first_neighbor_kernel(ii, jj);
        }
      } // i-loop
    }
  };

  Kokkos::RangePolicy<execution_space> policy_each(0, total_bins);
  Kokkos::parallel_for("calc_by_cell_list_each", policy_each, kernel_each);
  Kokkos::fence();

  Kokkos::RangePolicy<execution_space> policy_neighbor(
      0, total_pair_cell - total_bins - empty_pair_number);
  Kokkos::parallel_for("calc_by_cell_list_beighbor", policy_neighbor,
                       kernel_neighbor);
  Kokkos::fence();

  return verlet_list;
}
#endif // SHARED_MEMORY_PARALLELISM
