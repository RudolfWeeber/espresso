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
    Kokkos::View<int *, Kokkos::LayoutRight> &bin_size) {
  for (int cid = 0; cid < total_bins; ++cid) {
    int dx[3] = {};
    dx[0] = static_cast<int>(cid / (cell_num[1] * cell_num[2]));
    dx[1] = static_cast<int>((cid - dx[0] * (cell_num[1] * cell_num[2])) /
			     cell_num[2]);
    dx[2] = cid % cell_num[2];
    bin_offset(cid) = cell_list.binOffset(dx[0], dx[1], dx[2]);
    bin_size(cid) = cell_list.binSize(dx[0], dx[1], dx[2]);

    // Calculate particle_bins
    cell_list(cid);
  }
}

using ActiveProtocol = std::variant<LeesEdwards::Off, LeesEdwards::LinearShear,
				    LeesEdwards::OscillatoryShear>;
inline int set_interacting_pair_cell(
    int &total_bins, int total_pair_cell, int *cell_num, int *delta_lebc,
    int le_direction, int le_normal,
    std::shared_ptr<ActiveProtocol> le_protocol,
    //ActiveProtocol le_protocol,
    Kokkos::View<int *, Kokkos::LayoutRight> &bin_size,
    Cabana::LinkedCellList<memory_space, double> &cell_list,
    Kokkos::View<int **, Kokkos::LayoutRight> &interacting_pair_cell) {

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
    // Obtaining 3 dimentional cell index from cid_i
    int index[3] = {};
    cell_list.ijkBinIndex(cid_i, index[0], index[1], index[2]);
    int dx[3];
    // From 27 neighbor cell, the list of interacting pair cell is created
    for (int n = 0; n < 27; ++n) {
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
      if (cid_i <= cid_j) {
	if (bin_size(cid_i) != 0 and bin_size(cid_j) != 0) {
	  interacting_pair_cell(pair_cell_id, 0) = cid_i;
	  interacting_pair_cell(pair_cell_id, 1) = cid_j;
	  ++pair_cell_id;
	} else {
	  ++empty_pair_number;
	}
      }
    }
  }
  return empty_pair_number;
}

using ListAlgorithm = Cabana::HalfNeighborTag;
using ListType = Cabana::CustomVerletList<memory_space, ListAlgorithm,
					    Cabana::VerletLayout2D>;
//template <class SliceDouble3, class SliceInt, class SliceBool, class VerletCriterion, class Kernel>
template <class VerletCriterion, class Kernel>
ListType create_verlet_list(
		double const max_cutoff, int const max_counts,
		AoSoA_pack aosoa,
		std::vector<Particle*> unique_particles,
		VerletCriterion const &verlet_criterion,
		Kernel first_neighbor_kernel,
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
  //std::shared_ptr<ActiveProtocol> le_protocol = system.lees_edwards->get_protocol();
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
  cell_list = Cabana::createLinkedCellList<memory_space>(
      aosoa.position, grid_delta, grid_min, grid_max);
  int total_bins = cell_list.totalBins();
  // Now permute the AoSoA (i.e. reorder the data) using the linked cell
  // list.
  // Cabana::permute( cell_list, particle_storage );

  ListType verlet_list =
      ListType(aosoa.position, 0, aosoa.position.size(), max_counts);

  // Offset particle id and the number of particle in specific cell
  Kokkos::View<int *, Kokkos::LayoutRight> bin_offset("bin_offset",
						      total_bins);
  Kokkos::View<int *, Kokkos::LayoutRight> bin_size("bin_size", total_bins);
  set_offset_and_size_indexed_by_cid(total_bins, cell_num, cell_list,
				     bin_offset, bin_size);
  auto const particle_bins = cell_list.getParticleBins();

  // Creating Interacting cell
  int total_pair_cell;
  if (total_bins < 27) {
    total_pair_cell = (total_bins - 1) * total_bins / 2 + total_bins;
  } else {
    total_pair_cell = 14 * total_bins;
  }
  Kokkos::View<int **, Kokkos::LayoutRight> interacting_pair_cell(
      "interacting_pair_cell", total_pair_cell, 2);
  int empty_pair_number = set_interacting_pair_cell(
      total_bins, total_pair_cell, cell_num, delta_lebc, le_direction,
      le_normal, le_protocol, bin_size, cell_list, interacting_pair_cell);

  auto const distance_function = detail::MinimalImageDistance{
          std::as_const(cell_structure).decomposition().box()};

  // This kernel used the loop for the pair of interacting cell
  auto kernel = [&](const int pair_cell_i) {
    int cid_i = interacting_pair_cell(pair_cell_i, 0);
    int cid_j = interacting_pair_cell(pair_cell_i, 1);

    auto verlet_kernel = [&](Particle *p1, int ii, int id_i,
			     int cell_offset, int cell_size) {
      for (int j = cell_offset; j < cell_offset + cell_size; ++j) {
	//int ii = cell_list.permutation(i); // debug
	// int jj = j;
	int jj = cell_list.permutation(j);
	int id_j = aosoa.id(jj);
	if (aosoa.ghost(ii) or aosoa.ghost(jj)) {
	  if (((id_i < id_j) and aosoa.ghost(ii)) or
	      ((id_i > id_j) and aosoa.ghost(jj))) {
	    continue;
	  }
	} else if (aosoa.ghost(ii) and aosoa.ghost(jj)) {
	  continue; // reject both ghost
	}
	auto p2 = unique_particles.at(jj);//cell_structure.get_local_particle(id_j);
	//if (p2 == nullptr)
	//  continue;
	if (verlet_criterion(*p1, *p2, distance_function(*p1, *p2))) {
	  verlet_list.addNeighbor(ii, jj);
	  /*std::cout << "*Cabana* "
		    << i << " "
		    << j << " "
		    << id_i << " "
		    << id_j << " "
		    << aosoa.ghost(ii) << " "
		    << aosoa.ghost(jj) << " "
		    << cid_i << " "
		    << cid_j << " "
		    << aosoa.position(ii, 0) << ", "
		    << aosoa.position(ii, 1) << ", "
		    << aosoa.position(ii, 2) << " "
		    << aosoa.position(jj, 0) << ", "
		    << aosoa.position(jj, 1) << ", "
		    << aosoa.position(jj, 2) << "\n";//
	  //std::cout << "CHECK "
		    << n << " "
		    << i << " "
		    << j << " "
		    << dx[0] << " "
		    << dx[1] << " "
		    << dx[2] << "\n";*/
	  first_neighbor_kernel(ii, jj);
	}
      } // j-loop
    };

    int offset_i = bin_offset(cid_i);
    int size_i = bin_size(cid_i);

    for (int i = offset_i; i < offset_i + size_i; ++i) {
      // int ii = i;
      int ii = cell_list.permutation(i); //get previous id
      int id_i = aosoa.id(ii);
      auto p1 = unique_particles.at(ii); //cell_structure.get_local_particle(id_i);
      //if (p1 == nullptr)
      //  continue;

      if (cid_i == cid_j) {
	verlet_kernel(p1, ii, id_i, i + 1,
		      size_i + offset_i - i - 1); // j-loop
	// verlet_kernel(p1, i, id_i, i + 1,
	//               size_i + offset_i - i - 1); // j-loop
      } else {
	int offset_j = bin_offset(cid_j);
	int size_j = bin_size(cid_j);
	verlet_kernel(p1, ii, id_i, offset_j, size_j); // j-loop
	//verlet_kernel(p1, i, id_i, offset_j, size_j); // j-loop
      }
    } // i-loop
  };

  Kokkos::RangePolicy<execution_space> policy(0, total_pair_cell -
						     empty_pair_number);
  Kokkos::parallel_for("calc_by_cell_list", policy, kernel);
  Kokkos::fence();

  return verlet_list;
}
#endif // SHARED_MEMORY_PARALLELISM
