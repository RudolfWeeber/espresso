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

inline double wrap1(double x, double L) {
  auto result = x - std::floor(x / L) * L;
  return result;
}

inline double wrap2(double x, double L) {
  auto result = x - std::floor(x / L) * L;
  if (result >= L)
    result -= std::nextafter(L, 0.);
  return result;
}

inline bool contains(std::vector<int> const &storage, int const value) {
  return (std::find(storage.begin(), storage.end(), value) != storage.end());
}

template <class SliceDouble3, class SliceDouble, class SliceInt,
          class SliceBool>
inline void write_particle(Particle const &p, int const &id,
                           SliceDouble3 &s_position, SliceDouble3 &s_force,
                           SliceDouble3 &s_torque, SliceDouble &s_charge,
                           SliceInt &s_id, SliceInt &s_type, SliceBool &s_ghost,
                           Utils::Vector3d &box_l) {
  // SliceInt &s_id, SliceInt &s_type, SliceBool &s_ghost, BoxGeometry const
  // &box_geo) {
  auto const pos = p.pos();
  s_position(id, 0) = wrap2(pos[0], box_l[0]);
  s_position(id, 1) = wrap2(pos[1], box_l[1]);
  s_position(id, 2) = wrap2(pos[2], box_l[2]);
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
  assert(s_position(id, 0) >= 0. && s_position(id, 0) < box_l[0]);
  assert(s_position(id, 1) >= 0. && s_position(id, 1) < box_l[1]);
  assert(s_position(id, 2) >= 0. && s_position(id, 2) < box_l[2]);
  /*if (p.id() == 325) {
    std::cout << "0 CHECK 325 "
              << pos[0] << " "
              << pos[1] << " "
              << pos[2] << "\n";
  }
  if (p.id() == 512) {
    std::cout << "0 CHECK 512 "
              << pos[0] << " "
              << pos[1] << " "
              << pos[2] << "\n";
  }*/
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
    using data_types = Cabana::MemberTypes<double[3], double[3], double[3],
                                           double, int, int, int, bool>;
    using memory_space = Kokkos::SharedSpace;
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
    std::unordered_map<int, int> id_to_index{};
    std::vector<int> index_to_id{};
    int index = 0;

    bool const rebuild = cell_structure.get_rebuild_cabana_verlet_list();

    CabanaData saved_data;

    // Load saved data if we do not have to rebuild
    if (!rebuild) {
      saved_data = cell_structure.get_cabana_data();
    }

    // If we have to rebuild, we need to count the particles and create a new
    // map
    if (rebuild) {

      for (auto const &p : particles) {
        // id_to_index[p.id()] = index;
        index_to_id.emplace_back(p.id());
        index++;
      }

      for (auto const &p : ghost_particles) {
        // if (not id_to_index.contains(p.id())) {
        //   id_to_index[p.id()] = index;
        if (not contains(index_to_id, p.id())) {
          index_to_id.emplace_back(p.id());
          index++;
        }
      }
    } else {
      // If we do not rebuild we can use the saved map
      // id_to_index = saved_data.get_id_to_index();
      index_to_id = saved_data.get_index_to_id();
      index = id_to_index.size();
    }

    const int number_of_unique_particles = index;
#ifdef CALIPER
    CALI_MARK_END("Cabana - Index map");
#endif

    // ===================================================
    // Create and fill particle storage
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Fill particle storage");
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
    auto box_l = box_geo.length();
    int p_id = 0;
    std::vector<int> registered_pid{};
    // std::vector<int> ghost_pid{};
    for (auto const &p : particles) {
      write_particle(p, p_id, slice_position, slice_force, slice_torque,
                     slice_charge, slice_id, slice_type, slice_ghost, box_l);
      registered_pid.emplace_back(p.id());
      ++p_id;
    }
    for (auto const &p : ghost_particles) {
      // if the ghost is not in the previous map, but mpi moved it to this rank?
      // it will not have neighbors because we did not rebuild the verlet list.
      if (contains(registered_pid, p.id())) {
        continue;
      }
      write_particle(p, p_id, slice_position, slice_force, slice_torque,
                     slice_charge, slice_id, slice_type, slice_ghost, box_l);
      registered_pid.emplace_back(p.id());
      // ghost_pid.emplace_back(p.id());
      ++p_id;
    }

    using TP = decltype(slice_position);
    // using TF = decltype(slice_force);
    // using TR = decltype(slice_torque);
    using TQ = decltype(slice_charge);
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

    // ===================================================
    // Get Verlet Pairs and Fill list
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Verlet List");
#endif
    ListType verlet_list;
    std::vector<std::pair<int, int>> pair_check;

    // Rebuild verlet list if needed
    if (rebuild) {

      verlet_list = ListType(slice_position, 0, slice_position.size(), 64);
      /*
      auto kernel = [&](Particle const &p1, Particle const &p2) {
        verlet_list.addNeighbor(id_to_index.at(p1.id()),
                                id_to_index.at(p2.id()));
          std::cout << "Cell_structure "
                    << id_to_index.at(p1.id()) << " "
                    << id_to_index.at(p2.id()) << " "
                    << p1.id() << " "
                    << p2.id() << " "
                    << p1.pos() << " "
                    << p2.pos() << "\n";
        if (p1.id() < p2.id()) {
          pair_check.emplace_back(std::pair{p1.id(), p2.id()});
        } else {
          pair_check.emplace_back(std::pair{p2.id(), p1.id()});
        }
      };*/

      // cell_structure.cabana_verlet_list_loop(kernel, verlet_criterion);
    } else {
      // Else use the saved verlet list
      verlet_list = saved_data.get_verlet_list();
    }

    // Creating LinkedCellList and VerletList:
    Cabana::LinkedCellList<memory_space, double> cell_list;
    double grid_min[3] = {0.0, 0.0, 0.0};
    double grid_max[3] = {box_l[0], box_l[1], box_l[2]};
    double grid_delta[3] = {};
    int cell_num[3] = {};
    double max_cutoff = System::get_system().get_interaction_range();
    for (int d = 0; d < 3; ++d) {
      cell_num[d] = static_cast<int>(box_l[d] / max_cutoff);
      grid_delta[d] = std::nextafter(box_l[d] / cell_num[d], 0);
    }
    cell_list = Cabana::createLinkedCellList<memory_space>(
        slice_position, grid_delta, grid_min, grid_max);
    // Now permute the AoSoA (i.e. reorder the data) using the linked cell list.
    // Cabana::permute( cell_list, particle_storage );
    // ListType s_verlet_list;
    if (rebuild && max_cutoff != INACTIVE_CUTOFF) {
      // if (max_cutoff != INACTIVE_CUTOFF) {
      verlet_list = ListType(slice_position, 0, slice_position.size(), 64);
      for (int cid = 0; cid < cell_list.totalBins(); ++cid) {
        cell_list(cid);
      }
      auto const particle_bins = cell_list.getParticleBins();
      // std::vector< std::vector< std::vector<int> > > ijkIndexesInCell{};
      // for (int cid = 0; cid < cell_list.totalBins(); ++cid) {
      std::vector<std::vector<int>> ijkIndexes{};
      // int index[3] = {};
      // index[0] = static_cast<int>(cid / (cell_num[1] * cell_num[2]));
      // index[1] = static_cast<int>((cid - index[0] * (cell_num[1] *
      // cell_num[2])) / cell_num[2] ); index[2] = cid % cell_num[2];
      for (int n = 0; n < 27; ++n) {
        std::vector<int> dx = {0, 0, 0};
        dx[0] = static_cast<int>(n / 9);
        dx[1] = static_cast<int>((n - 9 * dx[0]) / 3);
        dx[2] = n % 3;
        // for (int d = 0; d < 3; ++d) {
        //   dx[d] = (-1 + dx[d] + index[d] + cell_num[d]) % cell_num[d];
        // }
        ijkIndexes.emplace_back(dx);
      }
      // ijkIndexesInCell.emplace_back(ijkIndexes);
      //}
      //
      auto const distance_function = detail::MinimalImageDistance{
          std::as_const(cell_structure).decomposition().box()};

      // auto kernel = [&cell_list, &particle_bins, &cell_num, &slice_id,
      // &cell_structure, &verlet_criterion, &verlet_list,
      // &distance_function](const int i) {
      auto kernel = [&](const int i) {
        int index[3] = {};
        cell_list.ijkBinIndex(particle_bins(i), index[0], index[1], index[2]);
        // auto ijkIndexes = ijkIndexesInCell[particle_bins(i)];
        int dx[3];
        for (int n = 0; n < 27; ++n) {
          // auto dx = ijkIndexes[n];
          // auto relative_index = ijkIndexes[n];
          dx[0] = static_cast<int>(n / 9);
          dx[1] = static_cast<int>((n - 9 * dx[0]) / 3);
          dx[2] = n % 3;
          for (int d = 0; d < 3; ++d) {
            dx[d] = (-1 + dx[d] + index[d] + cell_num[d]) % cell_num[d];
          }

          int offset = cell_list.binOffset(dx[0], dx[1], dx[2]);
          int size = cell_list.binSize(dx[0], dx[1], dx[2]);

          for (int j = offset; j < offset + size; j++) {
            // int jj = j;
            int jj = cell_list.permutation(j);
            if (slice_id(i) < slice_id(jj)) {
              auto p1 = cell_structure.get_local_particle(slice_id(i));
              auto p2 = cell_structure.get_local_particle(slice_id(jj));
              if (p1 == nullptr or p2 == nullptr)
                continue;
              if (p1->is_ghost()) {
                // if (slice_ghost(slice_id(i))) {
                // std::cout << slice_id(i) << " is ghost in rank " << rank <<
                // "\n";
              } else {
                if (verlet_criterion(*p1, *p2, distance_function(*p1, *p2))) {
                  verlet_list.addNeighbor(i, jj);
                  /*std::cout << "*Cabana* "
                            << i << " "
                            << j << " "
                            << slice_id(i) << " "
                            << slice_id(j) << " "
                            << slice_position(i, 0) << " "
                            << slice_position(i, 1) << " "
                            << slice_position(i, 2) << " "
                            << slice_position(j, 0) << " "
                            << slice_position(j, 1) << " "
                            << slice_position(j, 2) << "\n";*/
                }
              }
            }
          }
        }
      };

      Kokkos::RangePolicy<execution_space> policy(0, particle_storage.size());
      Kokkos::parallel_for("calc_by_cell_list", policy, kernel);
      Kokkos::fence();
    }

    // Save data for next iteration if we just rebuilt
    if (rebuild) {
      CabanaData new_data(verlet_list, id_to_index);
      cell_structure.set_cabana_data(std::make_unique<CabanaData>(new_data));
    }

    // calculate force with customverletlist
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
          TQ &slice_charge_, TI &slice_id_, TT &slice_type_,
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
          int num_threads_, int mpi_rank_)
          : cell(cell_), bonded_ias(bonded_ias_), nonbonded_ias(nonbonded_ias_),
            thermostat(thermostat_), box_geo(box_geo_),
            // index_to_id(index_to_id_),
            local_force(local_force_), local_torque(local_torque_),
            local_virial(local_virial_), slice_position(slice_position_),
            slice_charge(slice_charge_), slice_id(slice_id_),
            slice_type(slice_type_),
#ifdef COLLISION_DETECTION
            collision_detection(collision_detection_),
#endif
            coulomb_kernel(coulomb_kernel_), dipoles_kernel(dipoles_kernel_),
            elc_kernel(elc_kernel_), coulomb_u_kernel(coulomb_u_kernel_),
            num_threads(num_threads_), mpi_rank(mpi_rank_) {
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
        // std::cout << "in " << thread_id << " " << i << " " << j << " " <<
        // q1q2 << "\n";

        IA_parameters const &ia_params =
            nonbonded_ias.get_ia_param(slice_type(i), slice_type(j));
        /*
        auto p1 = cell->get_local_particle(slice_id(i));
        auto p2 = cell->get_local_particle(slice_id(j));

        if (p1 == nullptr or p2 == nullptr)
          return;

        auto[pf, virial] = add_non_bonded_pair_force(
          const_cast<Particle &>(*p1), const_cast<Particle &>(*p2),
          d, dist, dist2, q1q2, ia_params, thermostat, box_geo, bonded_ias,
          coulomb_kernel, dipoles_kernel, elc_kernel, coulomb_u_kernel);
        */
        //
        ParticleForce pf{};
        Utils::Vector3d virial{};

#ifdef EXCLUSIONS
        auto p1 = cell->get_local_particle(slice_id(i));
        auto p2 = cell->get_local_particle(slice_id(j));

        if (p1 == nullptr or p2 == nullptr)
          return;

        bool do_nonbonded_flag = do_nonbonded(*p1, *p2);
#else
        bool do_nonbonded_flag = true;
#endif

        add_non_bonded_pair_withot_p(pf, d, dist, q1q2, ia_params,
                                     do_nonbonded_flag, coulomb_kernel);

#if defined(THOLE) or defined(ELECTROSTATICS) or defined(P3M) or               \
    defined(DPD) or defined(DIPOLES)
        auto const dist2 = dist * dist;

#ifndef EXCLUSIONS
        auto p1 = cell->get_local_particle(slice_id(i));
        auto p2 = cell->get_local_particle(slice_id(j));

        if (p1 == nullptr or p2 == nullptr)
          return;
#endif
        add_non_bonded_pair_force_with_p(
            const_cast<Particle &>(*p1), const_cast<Particle &>(*p2), pf,
            virial, d, dist, dist2, q1q2, ia_params, do_nonbonded_flag,
            thermostat, box_geo, bonded_ias, coulomb_kernel, dipoles_kernel,
            elc_kernel, coulomb_u_kernel);
#endif
        //
        local_force(thread_id, i, 0) += pf.f[0];
        local_force(thread_id, i, 1) += pf.f[1];
        local_force(thread_id, i, 2) += pf.f[2];
        local_torque(thread_id, i, 0) += pf.torque[0];
        local_torque(thread_id, i, 1) += pf.torque[1];
        local_torque(thread_id, i, 2) += pf.torque[2];

        auto opf = calc_opposing_force(pf, d);
        local_force(thread_id, j, 0) += opf.f[0];
        local_force(thread_id, j, 1) += opf.f[1];
        local_force(thread_id, j, 2) += opf.f[2];
        local_torque(thread_id, j, 0) += opf.torque[0];
        local_torque(thread_id, j, 1) += opf.torque[1];
        local_torque(thread_id, j, 2) += opf.torque[2];

        /*if (p1->id() == 512 || p2->id() == 512) {
          std::cout << "2 CHECK FORCE i "
                    << i << " " << j << " "
                    << p1->id() << " "
                    << p2->id() << " "
                    << dist << " "
                    << p1->is_ghost() << " "
                    << p2->is_ghost() << " "
                    << pf.f[0] << " "
                    << pf.f[1] << " "
                    << pf.f[2] << "\n";
        }*/
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
#ifdef CALIPER
    CALI_MARK_END("Cabana - Verlet List");
#endif

    // ===================================================
    // Execute Kernel
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Execute Kernel");
#endif
    Kokkos::RangePolicy<execution_space> policy(0, particle_storage.size());

    FirstNeighborKernel first_neighbor_kernel(
        &cell_structure, bonded_ias, nonbonded_ias, thermostat, box_geo,
        // index_to_id,
        local_force, local_torque, local_virial, slice_position, slice_charge,
        slice_id, slice_type,
#ifdef COLLISION_DETECTION
        *collision_detection,
#endif
        coulomb_kernel, dipoles_kernel, elc_kernel, coulomb_u_kernel,
        num_threads, rank);

    // std::cout << rank << " " << num_threads << " Execute
    // FirstNeighborKernel\n";
    //   TODO: Add option to switch "SerialOpTag" Between "TeamOpTag"
    //   Feels like TeamOpTag is faster, atleast for large particle numbers
    Cabana::neighbor_parallel_for(policy, first_neighbor_kernel, verlet_list,
                                  Cabana::FirstNeighborsTag(),
                                  Cabana::SerialOpTag()); //, "verlet_list");
    Kokkos::fence();

    // Force and Torque reduction
    Kokkos::parallel_for(
        "reduce", policy, KOKKOS_LAMBDA(const int i) {
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

          /*if (slice_id(i) == 325 || slice_id(i) == 512) {
            std::cout << "3 CHECK FORCE "
                      << slice_id(i) << " "
                      << slice_force(i, 0) << " "
                      << slice_force(i, 1) << " "
                      << slice_force(i, 2) << "\n";
          }*/
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
    CALI_MARK_END("Cabana - Execute Kernel");
#endif

#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Collision Detection");
#endif
#ifdef COLLISION_DETECTION
    auto collision_kernel = [&](Particle const &p1, Particle const &p2,
                                Distance const &d) {
      if (not collision_detection->is_off()) {
        collision_detection->detect_collision(p1, p2, d.dist2);
      }
    };
    cell_structure.non_bonded_loop(collision_kernel, verlet_criterion);
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

      ParticleForce f(f_vec, torque_vec);
      p->force_and_torque() += f;
    }
    /*
    std::unordered_set<int> processed_ids;

    for (auto &p : ghost_particles) {
      int const pid = p.id();
      // Check if the particle has already been processed
      if (processed_ids.find(pid) != processed_ids.end()) {
        continue;
      }

      // Check if the ghost particle is in the map, i.e. was used during force
      // calculation
      if (id_to_index.find(pid) == id_to_index.end()) {
        continue;
      }

      auto const id = id_to_index.at(pid);

      // Only add forces to ghost particles that are not as normal particles in
      // the map, as they have already been added to the force calculation
      if (id < particles.size()) {
        continue;
      }

      processed_ids.insert(pid);

      Utils::Vector3d f_vec{slice_force(id, 0), slice_force(id, 1),
                            slice_force(id, 2)};
      Utils::Vector3d torque_vec{slice_torque(id, 0), slice_torque(id, 1),
                                 slice_torque(id, 2)};

      ParticleForce f(f_vec, torque_vec);
      p.force_and_torque() += f;
    }
    */
#ifdef CALIPER
    CALI_MARK_END("Cabana - Particle Forces");
#endif
  }
}

#endif
