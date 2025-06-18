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

#include "cabana_data.hpp"
#include "custom_verlet_list.hpp"
#include <Cabana_Core.hpp>
#include <Cabana_NeighborList.hpp>
#include <cassert>
#include <iostream>
#include <stdio.h>
#include <unordered_set>
#include <utility>

inline double wrap(double x, double L) {
  auto result = x - std::floor(x / L) * L;
  if (result >= L)
    result -= std::nextafter(L, 0.);
  return result;
}

template <class SliceDouble3, class SliceDouble, class SliceInt,
          class SliceBool>
inline void
write_particle(Particle &p, int const &id, SliceDouble3 &s_position,
               SliceDouble3 &s_force, SliceDouble3 &s_torque,
               SliceDouble &s_charge, SliceInt &s_id, SliceInt &s_type,
               SliceBool &s_ghost, SliceBool &s_has_excl,
               std::vector<Particle *> &s_pointers, Utils::Vector3d &box_l) {
  auto const pos = p.pos();
  s_position(id, 0) = wrap(pos[0], box_l[0]);
  s_position(id, 1) = wrap(pos[1], box_l[1]);
  s_position(id, 2) = wrap(pos[2], box_l[2]);
  s_id(id) = p.id();
  s_charge(id) = p.q();
  s_type(id) = p.type();
  s_ghost(id) = p.is_ghost();
  s_force(id, 0) = 0.0;
  s_force(id, 1) = 0.0;
  s_force(id, 2) = 0.0;
  s_torque(id, 0) = 0.0;
  s_torque(id, 1) = 0.0;
  s_torque(id, 2) = 0.0;
  s_has_excl(id) = (p.exclusions().size() != 0);
  s_pointers[id] = &p;
  assert(s_position(id, 0) >= 0. and s_position(id, 0) < box_l[0]);
  assert(s_position(id, 1) >= 0. and s_position(id, 1) < box_l[1]);
  assert(s_position(id, 2) >= 0. and s_position(id, 2) < box_l[2]);
}

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

template <class BondKernel, class VerletCriterion = detail::True>
void cabana_short_range(
    BondKernel bond_kernel,
    [[maybe_unused]] BondedInteractionsMap const &bonded_ias,
    Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel,
    Dipoles::ShortRangeForceKernel::kernel_type const *dipoles_kernel,
    Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const *elc_kernel,
    Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel,
#ifdef COLLISION_DETECTION
    std::shared_ptr<CollisionDetection::CollisionDetection> collision_detection,
#endif
    CellStructure &cell_structure, double pair_cutoff, double bond_cutoff,
    Thermostat::Thermostat const &thermostat, BoxGeometry const &box_geo,
    InteractionsNonBonded &nonbonded_ias, ParticleRange particles,
    ParticleRange ghost_particles,
    VerletCriterion const &verlet_criterion = {}) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

#ifdef CALIPER
  CALI_MARK_BEGIN("Espresso - Bond Kernel");
#endif

  assert(cell_structure.get_resort_particles() == Cells::RESORT_NONE);

  if (bond_cutoff >= 0.) {
    cell_structure.bond_loop(bond_kernel);
  }

#ifdef CALIPER
  CALI_MARK_END("Espresso - Bond Kernel");
#endif

  // Cabana short range loop
  if (pair_cutoff > 0.) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // ===================================================
    // Setup Cabana Variables
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Setup");
#endif
    // Dont know where to do this better
    using data_types = Cabana::MemberTypes<double[3], // 0
                                           double[3], // 1
                                           double[3], // 2
                                           double,    // 3
                                           int,       // 4
                                           int,       // 5
                                           bool,      // 6
                                           bool>;     // 7
    using memory_space = Kokkos::HostSpace;           // Kokkos::SharedSpace;
    using execution_space = Kokkos::DefaultExecutionSpace;

    using ListAlgorithm = Cabana::HalfNeighborTag;
    using ListType = Cabana::CustomVerletList<memory_space, ListAlgorithm,
                                              Cabana::VerletLayout2D>;

    // Number of threads
    const int num_threads = execution_space().concurrency();

    const int vector_length = 8;
#ifdef CALIPER
    CALI_MARK_END("Cabana - Setup");
#endif

    // ===================================================
    // Count unique particles and create Index map
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Index map");
#endif
    std::unordered_map<int, int> id_to_index{}; // For DEBUG
    std::unordered_set<int> registered_index{};
    std::vector<int> index_to_id{};
    int index = 0;

    bool const rebuild = cell_structure.get_rebuild_cabana_verlet_list();
    // std::cout << "For CABANA rebuild " << rebuild << std::endl;

    CabanaData saved_data;

    // Load saved data if we do not have to rebuild
    if (!rebuild) {
      saved_data = cell_structure.get_cabana_data();
    }

    // If we have to rebuild, we need to count the particles and create a new
    // map
    if (rebuild) {

      for (auto const &p : particles) {
        // if (cell_structure.get_local_particle(p.id())) {
        id_to_index[p.id()] = index;
        index++;
        //}
      }

      for (auto const &p : ghost_particles) {
        if (not id_to_index.contains(p.id())) {
          // if (cell_structure.get_local_particle(p.id())) {
          id_to_index[p.id()] = index;
          index++;
          //}
        }
      }
    } else {
      // If we do not rebuild we can use the saved map
      // id_to_index = saved_data.get_id_to_index();
      index = saved_data.get_index();
    }

    const int number_of_unique_particles = index;
#ifdef CALIPER
    CALI_MARK_END("Cabana - Index map");
#endif

    // ===================================================
    // Create and fill particle storage
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Allocation");
#endif
    Cabana::AoSoA<data_types, memory_space, vector_length> particle_storage(
        "particles", number_of_unique_particles);
    auto slice_position = Cabana::slice<0>(particle_storage);
    auto slice_force = Cabana::slice<1>(particle_storage);
    auto slice_torque = Cabana::slice<2>(particle_storage);
    auto slice_charge = Cabana::slice<3>(particle_storage);
    auto slice_id = Cabana::slice<4>(particle_storage);
    auto slice_type = Cabana::slice<5>(particle_storage);
    auto slice_ghost = Cabana::slice<6>(particle_storage);
    auto slice_has_excl = Cabana::slice<7>(particle_storage);
    std::vector<Particle *> slice_pointers(number_of_unique_particles);
    auto box_l = box_geo.length();
#ifdef CALIPER
    CALI_MARK_END("Cabana - Allocation");
#endif
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Fill particle storage");
#endif
    int p_id = 0;
    registered_index.clear();
    for (auto &p : particles) {
      // if (!cell_structure.get_local_particle(p.id())) continue;
      write_particle(p, p_id, slice_position, slice_force, slice_torque,
                     slice_charge, slice_id, slice_type, slice_ghost,
                     slice_has_excl, slice_pointers, box_l);
      registered_index.insert(p.id());
      ++p_id;
    }
    for (auto &p : ghost_particles) {
      // if the ghost is not in the previous map, but mpi moved it to this rank?
      // it will not have neighbors because we did not rebuild the verlet list.
      if (registered_index.contains(p.id())) {
        continue;
      }
      // if (!cell_structure.get_local_particle(p.id())) continue;
      write_particle(p, p_id, slice_position, slice_force, slice_torque,
                     slice_charge, slice_id, slice_type, slice_ghost,
                     slice_has_excl, slice_pointers, box_l);
      registered_index.insert(p.id());
      ++p_id;
    }

    using TP = decltype(slice_position);
    // using TF = decltype(slice_force);
    // using TR = decltype(slice_torque);
    using TQ = decltype(slice_charge);
    using TE = decltype(slice_has_excl);
    using TI = decltype(slice_id);
    using TT = decltype(slice_type);

    Kokkos::View<double ***, Kokkos::LayoutRight> local_force(
        "local_force", num_threads, number_of_unique_particles, 3);

    Kokkos::View<double ***, Kokkos::LayoutRight> local_torque(
        "local_torque", num_threads, number_of_unique_particles, 3);

    Kokkos::View<double **, Kokkos::LayoutRight> local_virial("local_virial",
                                                              num_threads, 3);

#ifdef CALIPER
    CALI_MARK_END("Cabana - Fill particle storage");
#endif

    // The kernel of calculate force
    struct FirstNeighborKernel {
      const CellStructure *cell;
      [[maybe_unused]] const BondedInteractionsMap &bonded_ias;
      const InteractionsNonBonded &nonbonded_ias;
      const Thermostat::Thermostat &thermostat;
      const BoxGeometry &box_geo;
      // std::vector<int> &index_to_id;
      Kokkos::View<double ***> local_force;
      Kokkos::View<double ***> local_torque;
      Kokkos::View<double **> local_virial;
      TP &slice_position;
      TQ &slice_charge;
      TI &slice_id;
      TT &slice_type;
      TE &slice_has_excl;
#ifdef COLLISION_DETECTION
      // std::shared_ptr<CollisionDetection::CollisionDetection>
      // collision_detection;
      mutable CollisionDetection::CollisionDetection collision_detection;
#endif
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel;
      Dipoles::ShortRangeForceKernel::kernel_type const *dipoles_kernel;
      Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const *elc_kernel;
      Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel;

      int num_threads;
      int mpi_rank;
      int particle_number;

      FirstNeighborKernel(
          const CellStructure *cell_,
          [[maybe_unused]] const BondedInteractionsMap &bonded_ias_,
          const InteractionsNonBonded &nonbonded_ias_,
          const Thermostat::Thermostat &thermostat_,
          const BoxGeometry &box_geo_,
          // std::vector<int> &index_to_id_,
          Kokkos::View<double ***> local_force_,
          Kokkos::View<double ***> local_torque_,
          Kokkos::View<double **> local_virial_, TP &slice_position_,
          TQ &slice_charge_, TI &slice_id_, TT &slice_type_, TE &slice_hsa_excl,
#ifdef COLLISION_DETECTION
          // std::shared_ptr<CollisionDetection::CollisionDetection>
          // collision_detection_,
          CollisionDetection::CollisionDetection collision_detection_,
#endif
          Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel_,
          Dipoles::ShortRangeForceKernel::kernel_type const *dipoles_kernel_,
          Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const
              *elc_kernel_,
          Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel_,
          int num_threads_, int mpi_rank_, int particle_number_)
          : cell(cell_), bonded_ias(bonded_ias_), nonbonded_ias(nonbonded_ias_),
            thermostat(thermostat_), box_geo(box_geo_),
            // index_to_id(index_to_id_),
            local_force(local_force_), local_torque(local_torque_),
            local_virial(local_virial_), slice_position(slice_position_),
            slice_charge(slice_charge_), slice_id(slice_id_),
            slice_type(slice_type_), slice_has_excl(slice_has_excl),
#ifdef COLLISION_DETECTION
            collision_detection(collision_detection_),
#endif
            coulomb_kernel(coulomb_kernel_), dipoles_kernel(dipoles_kernel_),
            elc_kernel(elc_kernel_), coulomb_u_kernel(coulomb_u_kernel_),
            num_threads(num_threads_), mpi_rank(mpi_rank_),
            particle_number(particle_number_) {
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(int i, int j) const {

        Utils::Vector3d const pi = {slice_position(i, 0), slice_position(i, 1),
                                    slice_position(i, 2)};
        Utils::Vector3d const pj = {slice_position(j, 0), slice_position(j, 1),
                                    slice_position(j, 2)};

        Utils::Vector3d const d = box_geo.get_mi_vector(pi, pj);
        auto const dist = d.norm();

        auto const q1q2 = slice_charge(i) * slice_charge(j);

        auto thread_id = omp_get_thread_num();
        // auto thread_id = Kokkos::OpenMP::impl_hardware_thread_id();
        // std::cout << "\nin " << thread_id << "\n"; //" " << i << " " << j <<
        // " " <<

        IA_parameters const &ia_params =
            nonbonded_ias.get_ia_param(slice_type(i), slice_type(j));
        //
        ParticleForce pf{};
        Utils::Vector3d virial{};
        bool do_nonbonded_flag = true;
#ifdef EXCLUSIONS
        if (slice_has_excl(i) and slice_has_excl(j)) {
          //        auto p1 = cell->get_local_particle(slice_id(i));
          //        auto p2 = cell->get_local_particle(slice_id(j));

          //        if (p1 == nullptr or p2 == nullptr)
          //          return;

          //        do_nonbonded_flag = do_nonbonded(*p1, *p2);
        }
#endif

        add_non_bonded_pair_withot_p(pf, d, dist, q1q2, ia_params,
                                     do_nonbonded_flag, coulomb_kernel);

#if defined(THOLE) or defined(ELECTROSTATICS) or defined(P3M) or               \
    defined(DPD) or defined(DIPOLES)
        auto const dist2 = dist * dist;

        auto p1 = cell->get_local_particle(slice_id(i));
        auto p2 = cell->get_local_particle(slice_id(j));

        if (p1 == nullptr or p2 == nullptr)
          return;
        add_non_bonded_pair_force_with_p(
            const_cast<Particle &>(*p1), const_cast<Particle &>(*p2), pf,
            virial, d, dist, dist2, q1q2, ia_params, do_nonbonded_flag,
            thermostat, box_geo, bonded_ias, coulomb_kernel, dipoles_kernel,
            elc_kernel, coulomb_u_kernel);
#endif
        local_force(thread_id, i, 0) += pf.f[0];
        local_force(thread_id, i, 1) += pf.f[1];
        local_force(thread_id, i, 2) += pf.f[2];
#ifdef ROTATION
        local_torque(thread_id, i, 0) += pf.torque[0];
        local_torque(thread_id, i, 1) += pf.torque[1];
        local_torque(thread_id, i, 2) += pf.torque[2];
#endif
        auto opf = calc_opposing_force(pf, d);
        local_force(thread_id, j, 0) += opf.f[0];
        local_force(thread_id, j, 1) += opf.f[1];
        local_force(thread_id, j, 2) += opf.f[2];
#ifdef ROTATION
        local_torque(thread_id, j, 0) += opf.torque[0];
        local_torque(thread_id, j, 1) += opf.torque[1];
        local_torque(thread_id, j, 2) += opf.torque[2];
#endif
#ifdef NPT
        local_virial(thread_id, 0) += virial[0];
        local_virial(thread_id, 1) += virial[1];
        local_virial(thread_id, 2) += virial[2];
#endif

#ifdef COLLISION_DETECTION
        // if (not collision_detection.is_off()) {
        //   collision_detection.detect_collision(*p1, *p2, dist2);
        // }
#endif
      };
    };

    //  START VERLET_LIST
    // ===================================================
    // Get Verlet Pairs and Fill list
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Verlet List1");
#endif
    ListType verlet_list;

    // Rebuild verlet list if needed
    auto const &system = ::System::get_system();
    int max_counts;
    double max_cutoff = pair_cutoff; // system.get_interaction_range();
    if (std::isinf(max_cutoff)) {
      max_counts = number_of_unique_particles;
    } else {
      max_counts =
          static_cast<int>(27 * max_cutoff * max_cutoff * max_cutoff / 3);
    }
    if (max_counts < 256)
      max_counts = 256;
    if (rebuild) { // Legacy Velert List
      /*verlet_list =
          ListType(slice_position, 0, slice_position.size(), max_counts);
      auto kernel = [&](Particle const &p1, Particle const &p2) {
        verlet_list.addNeighbor(id_to_index.at(p1.id()),
                                id_to_index.at(p2.id()));
          //std::cout << "WITHSMP "
                    //<< id_to_index.at(p1.id()) << " "
                    //<< id_to_index.at(p2.id()) << " "
                    //<< p1.is_ghost() << " "
                    //<< p2.is_ghost() << " "
                    //<< p1.id() << " "
                    //<< p2.id() << std::endl;
                    //<< p1.pos() << " "
                    //<< p2.pos() << "\n";
      };

      cell_structure.cabana_verlet_list_loop(kernel, verlet_criterion);*/
    } else {
      // Else use the saved verlet list
      verlet_list = saved_data.get_verlet_list();
    }
#ifdef CALIPER
    CALI_MARK_END("Cabana - Verlet List1");
#endif

    FirstNeighborKernel first_neighbor_kernel(
        &cell_structure, bonded_ias, nonbonded_ias, thermostat, box_geo,
        // index_to_id,
        local_force, local_torque, local_virial, slice_position, slice_charge,
        slice_id, slice_type, slice_has_excl,
#ifdef COLLISION_DETECTION
        *collision_detection,
#endif
        coulomb_kernel, dipoles_kernel, elc_kernel, coulomb_u_kernel,
        num_threads, rank, number_of_unique_particles);

    if (rebuild and max_cutoff != INACTIVE_CUTOFF) { // Shared memory
#ifdef CALIPER
      CALI_MARK_BEGIN("Cabana - Verlet Build");
#endif
      // Creating LinkedCellList and VerletList:
      // Box Properties
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
      cell_list = Cabana::createLinkedCellList<memory_space>(
          slice_position, grid_delta, grid_min, grid_max);
      int total_bins = cell_list.totalBins();
      // Now permute the AoSoA (i.e. reorder the data) using the linked cell
      // list.
      // Cabana::permute( cell_list, particle_storage );

      verlet_list =
          ListType(slice_position, 0, slice_position.size(), max_counts);

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
            // int ii = cell_list.permutation(i); // debug
            // int jj = j;
            int jj = cell_list.permutation(j);
            int id_j = slice_id(jj);
            if (slice_ghost(ii) or slice_ghost(jj)) {
              if (((id_i < id_j) and slice_ghost(ii)) or
                  ((id_i > id_j) and slice_ghost(jj))) {
                continue;
              }
            } else if (slice_ghost(ii) and slice_ghost(jj)) {
              continue; // reject both ghost
            }
            auto p2 = cell_structure.get_local_particle(id_j);
            if (p2 == nullptr)
              continue;
            if (verlet_criterion(*p1, *p2, distance_function(*p1, *p2))) {
              verlet_list.addNeighbor(ii, jj);
              /*std::cout << "*Cabana* "
                        << i << " "
                        << j << " "
                        << id_i << " "
                        << id_j << " "
                        << slice_ghost(ii) << " "
                        << slice_ghost(jj) << " "
                        << cid_i << " "
                        << cid_j << " "
                        << slice_position(ii, 0) << ", "
                        << slice_position(ii, 1) << ", "
                        << slice_position(ii, 2) << " "
                        << slice_position(jj, 0) << ", "
                        << slice_position(jj, 1) << ", "
                        << slice_position(jj, 2) << "\n";//
              //std::cout << "CHECK "
                        << n << " "
                        << i << " "
                        << j << " "
                        << dx[0] << " "
                        << dx[1] << " "
                        << dx[2] << "\n";*/
              // first_neighbor_kernel(ii, jj);
            }
          } // j-loop
        };

        int offset_i = bin_offset(cid_i);
        int size_i = bin_size(cid_i);

        for (int i = offset_i; i < offset_i + size_i; ++i) {
          // int ii = i;
          int ii = cell_list.permutation(i);
          int id_i = slice_id(ii);
          auto p1 = cell_structure.get_local_particle(id_i);
          if (p1 == nullptr)
            continue;

          if (cid_i == cid_j) {
            verlet_kernel(p1, ii, id_i, i + 1,
                          size_i + offset_i - i - 1); // j-loop
            // verlet_kernel(p1, i, id_i, i + 1,
            //               size_i + offset_i - i - 1); // j-loop
          } else {
            int offset_j = bin_offset(cid_j);
            int size_j = bin_size(cid_j);
            verlet_kernel(p1, ii, id_i, offset_j, size_j); // j-loop
            // verlet_kernel(p1, i, id_i, offset_j, size_j); // j-loop
          }
        } // i-loop
      };

      Kokkos::RangePolicy<execution_space> policy(0, total_pair_cell -
                                                         empty_pair_number);
      Kokkos::parallel_for("calc_by_cell_list", policy, kernel);
      Kokkos::fence();
#ifdef CALIPER
      CALI_MARK_END("Cabana - Verlet Build");
#endif
    } // else {
    {
#ifdef CALIPER
      CALI_MARK_BEGIN("Cabana - Verlet Use");
#endif
      Kokkos::RangePolicy<execution_space> policy(0, particle_storage.size());
      Cabana::neighbor_parallel_for(policy, first_neighbor_kernel, verlet_list,
                                    Cabana::FirstNeighborsTag(),
                                    Cabana::TeamOpTag());
      Kokkos::fence();
#ifdef CALIPER
      CALI_MARK_END("Cabana - Verlet Use");
#endif
    }

#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Calc Forces");
#endif

    // Save data for next iteration if we just rebuilt
    if (rebuild) {
      CabanaData new_data(verlet_list, particle_storage.size());
      cell_structure.set_cabana_data(std::make_unique<CabanaData>(new_data));
    }

    // Force and Torque reduction
    Kokkos::RangePolicy<execution_space> policy(0, particle_storage.size());
    Kokkos::parallel_for(
        "reduction", policy, KOKKOS_LAMBDA(const int i) {
          double fx = 0.;
          double fy = 0.;
          double fz = 0.;
          double tx = 0.;
          double ty = 0.;
          double tz = 0.;
          for (int tid = 0; tid < num_threads; ++tid) {
            fx += local_force(tid, i, 0);
            fy += local_force(tid, i, 1);
            fz += local_force(tid, i, 2);
            tx += local_torque(tid, i, 0);
            ty += local_torque(tid, i, 1);
            tz += local_torque(tid, i, 2);
          }
          slice_force(i, 0) = fx;
          slice_force(i, 1) = fy;
          slice_force(i, 2) = fz;
          slice_torque(i, 0) = tx;
          slice_torque(i, 1) = ty;
          slice_torque(i, 2) = tz;
        });
    Kokkos::fence();

#ifdef NPT
    double vx = 0.;
    double vy = 0.;
    double vz = 0.;
    for (int tid = 0; tid < num_threads; ++tid) {
      vx += local_virial(tid, 0);
      vy += local_virial(tid, 1);
      vz += local_virial(tid, 2);
    }
    Utils::Vector3d virial_vec{vx, vy, vz};
    npt_add_virial_force_contribution(virial_vec);
#endif
#ifdef CALIPER
    CALI_MARK_END("Cabana - Calc Forces");
#endif

#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Collision Detection");
#endif
#ifdef COLLISION_DETECTION
    auto collision_kernel = [&](Particle const &p1, Particle const &p2,
                                Distance const &d) {
      collision_detection->detect_collision(p1, p2, d.dist2);
    };
    if (not collision_detection->is_off()) {
      cell_structure.non_bonded_loop(collision_kernel, verlet_criterion);
    }
#endif
#ifdef CALIPER
    CALI_MARK_END("Cabana - Collision Detection");
#endif

    // ===================================================
    // Add forces to particles
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Particle Forces");
#endif
    for (auto id = 0; id < particle_storage.size(); ++id) {
      auto p = cell_structure.get_local_particle(slice_id(id));
      if (p == nullptr) {
        return;
      }
      Utils::Vector3d f_vec{slice_force(id, 0), slice_force(id, 1),
                            slice_force(id, 2)};
      Utils::Vector3d torque_vec{slice_torque(id, 0), slice_torque(id, 1),
                                 slice_torque(id, 2)};
#ifdef ROTATION
      ParticleForce f(f_vec, torque_vec);
#else
      ParticleForce f(f_vec);
#endif

      p->force_and_torque() += f;
    }
#ifdef CALIPER
    CALI_MARK_END("Cabana - Particle Forces");
#endif
  }
}

#endif
