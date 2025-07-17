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

inline double wrap(double x, double L) {
  auto result = x - std::floor(x / L) * L;
  // if (result >= L)
  //   result -= std::nextafter(L, 0.);
  return result;
}

inline void write_particle(Particle const &p, int const &id,
                           AoSoA_pack &aosoa) {
  // Utils::Vector3d const &box_l) {
  aosoa.id(id) = p.id();
  aosoa.charge(id) = p.q();
  aosoa.type(id) = p.type();
  // aosoa.ghost(id) = p.is_ghost();
  auto const pos = p.pos();
  for (int d = 0; d < 3; ++d) {
    aosoa.position(id, d) = pos[d];
  }
  /*double wpos[3] = {};
  for (int d = 0; d < 3; ++d) {
    // aosoa.position(id, d) =
    //     pos[d] - std::floor(pos[d] / box_l[d]) * box_l[d];
    wpos[d] = pos[d] - std::floor(pos[d] / box_l[d]) * box_l[d];
  }
  for (int d = 0; d < 3; ++d) {
    aosoa.position(id, d) = wpos[d];
  }*/
  // assert(aosoa.position(id, 0) >= 0. and aosoa.position(id, 0) < box_l[0]);
  // assert(aosoa.position(id, 1) >= 0. and aosoa.position(id, 1) < box_l[1]);
  // assert(aosoa.position(id, 2) >= 0. and aosoa.position(id, 2) < box_l[2]);
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
    // ===================================================
    // Count unique particles and create Index map
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Index map");
#endif
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Number of threads
    int num_threads = execution_space().concurrency();

    std::vector<Particle *> unique_particles;
    int index = 0;
    int max_id = 0;

    bool const rebuild = cell_structure.get_rebuild_cabana_verlet_list() or
                         (not cell_structure.use_verlet_list);
    // if (rank == 0) {
    // std::cout << "\nFor CABANA rebuild " << rebuild
    //	      << " " << rank << std::endl;
    // }

    CabanaData saved_data;

    // If we have to rebuild, we need to count the particles
    if (rebuild) {

      std::unordered_set<int> registered_index{};
      // std::bitset<1000000> registered_index;
      for (auto &p : particles) {
        if (p.id() > max_id)
          max_id = p.id();
        // registered_index.set(p.id());
        unique_particles.emplace_back(&p);
        index++;
      }

      for (auto &p : ghost_particles) {
        if (not cell_structure.get_local_particle(p.id())) {
          continue;
        }
        if (not cell_structure.get_local_particle(p.id())->is_ghost()) {
          continue;
        }
        if (registered_index.contains(p.id())) {
          continue;
        }
        // if (registered_index.test(p.id())) {
        //   continue;
        // }
        if (p.id() > max_id)
          max_id = p.id();
        registered_index.insert(p.id());
        // registered_index.set(p.id());
        unique_particles.emplace_back(&p);
        index++;
      }
      registered_index.clear();
    } else {
      // If we do not rebuild we can use the saved map
      saved_data = cell_structure.get_cabana_data();
      index = saved_data.get_index();
      unique_particles = saved_data.get_unique_particles();
      max_id = saved_data.get_max_id();
    }

    int number_of_unique_particles = index;

    // ===================================================
    // Create essential variable for MD
    // ===================================================
    Kokkos::View<double **[3], Kokkos::LayoutRight> local_force(
        "local_force", number_of_unique_particles, num_threads);

#ifdef ROTATION
    Kokkos::View<double **[3], Kokkos::LayoutRight> local_torque(
        "local_torque", number_of_unique_particles, num_threads);
#endif
#ifdef NPT
    Kokkos::View<double *[3], Kokkos::LayoutRight> local_virial("local_virial",
                                                                num_threads);
#endif
    Cabana::AoSoA<data_types, memory_space, vector_length> particle_storage(
        "particles", number_of_unique_particles);
    particle_storage.resize(number_of_unique_particles);
    // particle properties are defined in aosoa_pack.hpp
    auto aosoa = AoSoA_pack(particle_storage);

    using ListAlgorithm = Cabana::HalfNeighborTag;
    using ListType = Cabana::CustomVerletList<memory_space, ListAlgorithm,
                                              Cabana::VerletLayout2D>;
    ListType verlet_list;

#ifdef CALIPER
    CALI_MARK_END("Cabana - Index map");
#endif

    // The kernel of calculate force
    struct FirstNeighborKernel {
#if defined(EXCLUSIONS) or defined(THOLE) or defined(ELECTROSTATICS) or        \
    defined(P3M) or defined(DPD) or defined(DIPOLES) or defined(NPT)
      std::vector<Particle *> unique_particles;
#endif
      [[maybe_unused]] const BondedInteractionsMap bonded_ias;
      const InteractionsNonBonded nonbonded_ias;
      const BoxGeometry box_geo;
      Kokkos::View<double **[3]> local_force;
#ifdef ROTATION
      Kokkos::View<double **[3]> local_torque;
#endif
#ifdef NPT
      Kokkos::View<double *[3]> local_virial;
#endif
#ifdef COLLISION_DETECTION
      // std::shared_ptr<CollisionDetection::CollisionDetection>
      // collision_detection;
      mutable CollisionDetection::CollisionDetection collision_detection;
#endif
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel;
#if defined(THOLE) or defined(ELECTROSTATICS) or defined(P3M) or               \
    defined(DPD) or defined(DIPOLES) or defined(NPT)
      Dipoles::ShortRangeForceKernel::kernel_type const *dipoles_kernel;
      Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const *elc_kernel;
      Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel;
      const Thermostat::Thermostat &thermostat;
#endif
      // int num_threads;
      // int mpi_rank;
      // int particle_number;
      const AoSoA_pack aosoa;

      FirstNeighborKernel(
      // const CellStructure *cell_,
#if defined(EXCLUSIONS) or defined(THOLE) or defined(ELECTROSTATICS) or        \
    defined(P3M) or defined(DPD) or defined(DIPOLES) or defined(NPT)
          std::vector<Particle *> &unique_particles_,
#endif
          [[maybe_unused]] const BondedInteractionsMap &bonded_ias_,
          const InteractionsNonBonded &nonbonded_ias_,
          const BoxGeometry &box_geo_, Kokkos::View<double **[3]> local_force_,
#ifdef ROTATION
          Kokkos::View<double **[3]> local_torque_,
#endif
#ifdef NPT
          Kokkos::View<double *[3]> local_virial_,
#endif
#ifdef COLLISION_DETECTION
          // std::shared_ptr<CollisionDetection::CollisionDetection>
          // collision_detection_,
          CollisionDetection::CollisionDetection collision_detection_,
#endif
          Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel_,
#if defined(THOLE) or defined(ELECTROSTATICS) or defined(P3M) or               \
    defined(DPD) or defined(DIPOLES) or defined(NPT)
          Dipoles::ShortRangeForceKernel::kernel_type const *dipoles_kernel_,
          Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const
              *elc_kernel_,
          Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel_,
          const Thermostat::Thermostat &thermostat_,
#endif
          const AoSoA_pack &aosoa_)
          // int num_threads_), int mpi_rank_, int particle_number_)
          : // cell(cell_),
#if defined(EXCLUSIONS) or defined(THOLE) or defined(ELECTROSTATICS) or        \
    defined(P3M) or defined(DPD) or defined(DIPOLES) or defined(NPT)
            unique_particles(unique_particles_),
#endif
            bonded_ias(bonded_ias_), nonbonded_ias(nonbonded_ias_),
            box_geo(box_geo_), local_force(local_force_),
#ifdef ROTATION
            local_torque(local_torque_),
#endif
#ifdef NPT
            local_virial(local_virial_),
#endif
#ifdef COLLISION_DETECTION
            collision_detection(collision_detection_),
#endif
            coulomb_kernel(coulomb_kernel_),
#if defined(THOLE) or defined(ELECTROSTATICS) or defined(P3M) or               \
    defined(DPD) or defined(DIPOLES) or defined(NPT)
            dipoles_kernel(dipoles_kernel_), elc_kernel(elc_kernel_),
            coulomb_u_kernel(coulomb_u_kernel_), thermostat(thermostat_),
#endif
            aosoa(aosoa_) {
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(int i, int j) const {

        auto thread_id = omp_get_thread_num();
        // auto thread_id = Kokkos::OpenMP::impl_hardware_thread_id();
        // std::cout << "\nin " << thread_id << "\n"; //" " << i << " " << j <<
        // " " <<

        IA_parameters const &ia_params =
            nonbonded_ias.get_ia_param(aosoa.type(i), aosoa.type(j));

        ParticleForce pf{};
#ifdef NPT
        Utils::Vector3d virial{};
#endif
        Utils::Vector3d const d = box_geo.get_mi_vector(
            aosoa.position(i, 0), aosoa.position(i, 1), aosoa.position(i, 2),
            aosoa.position(j, 0), aosoa.position(j, 1), aosoa.position(j, 2));
        auto const dist = d.norm();

        auto const q1q2 = aosoa.charge(i) * aosoa.charge(j);

#ifdef EXCLUSIONS
        auto p1 = unique_particles.at(i);
        auto p2 = unique_particles.at(j);
        // auto p1 = cell->get_local_particle(aosoa.id(i));
        // auto p2 = cell->get_local_particle(aosoa.id(j));

        // if (p1 == nullptr or p2 == nullptr)
        //   return;

        bool do_nonbonded_flag = do_nonbonded(*p1, *p2);
#else
        bool do_nonbonded_flag = true;
#endif

        add_non_bonded_pair_withot_p(pf, d, dist, q1q2, ia_params,
                                     do_nonbonded_flag, coulomb_kernel);

#if defined(THOLE) or defined(ELECTROSTATICS) or defined(P3M) or               \
    defined(DPD) or defined(DIPOLES) or defined(NPT)
        auto const dist2 = dist * dist;

#ifndef EXCLUSIONS
        auto p1 = unique_particles.at(i);
        auto p2 = unique_particles.at(j);
        // auto p1 = cell->get_local_particle(aosoa.id(i));
        // auto p2 = cell->get_local_particle(aosoa.id(j));

        // if (p1 == nullptr or p2 == nullptr)
        //   return;
#endif // NOT EXCLUSIONS
        add_non_bonded_pair_force_with_p(
            const_cast<Particle &>(*p1), const_cast<Particle &>(*p2), pf,
#ifdef NPT
            virial,
#endif // NPT
            d, dist, dist2, q1q2, ia_params, do_nonbonded_flag, thermostat,
            box_geo, bonded_ias, coulomb_kernel, dipoles_kernel, elc_kernel,
            coulomb_u_kernel);
#endif // ETC
       //
        local_force(i, thread_id, 0) += pf.f[0];
        local_force(i, thread_id, 1) += pf.f[1];
        local_force(i, thread_id, 2) += pf.f[2];
#ifdef ROTATION
        local_torque(i, thread_id, 0) += pf.torque[0];
        local_torque(i, thread_id, 1) += pf.torque[1];
        local_torque(i, thread_id, 2) += pf.torque[2];
#endif

        auto opf = calc_opposing_force(pf, d);
        local_force(j, thread_id, 0) += opf.f[0];
        local_force(j, thread_id, 1) += opf.f[1];
        local_force(j, thread_id, 2) += opf.f[2];
#ifdef ROTATION
        local_torque(j, thread_id, 0) += opf.torque[0];
        local_torque(j, thread_id, 1) += opf.torque[1];
        local_torque(j, thread_id, 2) += opf.torque[2];
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

    // Fill the essential variable for MD
    {
#ifdef CALIPER
      CALI_MARK_BEGIN("Cabana - Allocation");
#endif
      // ===================================================
      // Fill particle storage
      // ===================================================

      Kokkos::View<int *> id_to_index(
          Kokkos::ViewAllocateWithoutInitializing("id_to_index"), max_id + 1);
      // Kokkos::view_alloc("id_to_index", Kokkos::WithoutInitializing), max_id
      // + 1);
      Kokkos::deep_copy(id_to_index, -1);

      auto box_l = box_geo.length();

      using policy_type = Kokkos::RangePolicy<execution_space>;
      Kokkos::parallel_for(
          "AoSoA write", policy_type(0, particle_storage.size()),
          //[&unique_particles, &aosoa, &box_l, &id_to_index](const int p_id) {
          [&unique_particles, &aosoa, &id_to_index](const int p_id) {
            // write_particle(*unique_particles.at(p_id), p_id, aosoa, box_l);
            write_particle(*unique_particles.at(p_id), p_id, aosoa);
            id_to_index(unique_particles.at(p_id)->id()) = p_id;
          });
      Kokkos::fence();
      // After ONLY JUST creating LinkedCellList, force calculation became
      // slower, even if it is not used and It is explicitly deleted.
      if (0) {
        // Cabana::LinkedCellList<memory_space, double> cell_list;
        double grid_min[3] = {0.0, 0.0, 0.0};
        double grid_max[3] = {box_l[0], box_l[1], box_l[2]};
        double grid_delta[3] = {};
        int cell_num[3] = {};
        double eff_cutoff;
        for (int d = 0; d < 3; ++d) {
          eff_cutoff = pair_cutoff;
          if (eff_cutoff > box_l[d])
            eff_cutoff = box_l[d];
          cell_num[d] = static_cast<int>(box_l[d] / eff_cutoff);
          grid_delta[d] = std::nextafter(box_l[d] / cell_num[d], 0);
        }
        auto *cell_list = new Cabana::LinkedCellList<memory_space>(
            aosoa.position, grid_delta, grid_min, grid_max);
        // Now permute the AoSoA (i.e. reorder the data)
        Cabana::permute(*cell_list, particle_storage);
        unique_particles.clear();
        for (int i = 0; i < aosoa.id.size(); ++i) {
          id_to_index(aosoa.id(i)) = i;
          unique_particles.emplace_back(
              cell_structure.get_local_particle(aosoa.id(i)));
        }
        delete cell_list;
        Kokkos::fence();
        /*Kokkos::parallel_for("AoSoA write", policy_type(0,
           particle_storage.size()),
                       [&unique_particles, &aosoa, &box_l](const int p_id) {
                         write_particle(*unique_particles.at(p_id), p_id,
                                        aosoa, box_l);
                         });*/
      }

#ifdef CALIPER
      CALI_MARK_END("Cabana - Allocation");
#endif

      // ===================================================
      // Get Verlet Pairs and Fill list
      // ===================================================

      // Rebuild verlet list if needed
      double max_cutoff = pair_cutoff; // system.get_interaction_range();
      if (rebuild and max_cutoff != INACTIVE_CUTOFF) { // Shared memory
#ifdef CALIPER
        CALI_MARK_BEGIN("Cabana - Verlet List");
#endif
        /*
        verlet_list = create_verlet_list(
            max_cutoff, max_counts, aosoa, unique_particles, verlet_criterion,
            first_neighbor_kernel, cell_structure);
        */
        bool at_steepest_descent = cell_structure.get_steepest_descent_flag();
        int max_counts;
        if (std::isinf(max_cutoff)) {
          max_counts = number_of_unique_particles;
        } else {
          int max_prefactor;
          if (at_steepest_descent) {
            max_prefactor = 8;
          } else {
            max_prefactor = 5;
          }
          max_counts = static_cast<int>(
              std::ceil(max_prefactor * max_cutoff * max_cutoff * max_cutoff));
        }
        int threshold_num = 8;
#ifdef COLLISION_DETECTION
        threshold_num = 64;
#endif
        if (max_counts < threshold_num) {
          max_counts = std::min(threshold_num, number_of_unique_particles);
        }
        // std::cout << "max_counts:" << max_counts << " " << max_cutoff <<
        // std::endl;
        verlet_list = ListType(0, number_of_unique_particles, max_counts);
        auto const &cells =
            std::as_const(cell_structure).decomposition().local_cells();
        auto const distance_function = detail::MinimalImageDistance{
            std::as_const(cell_structure).decomposition().box()};

        auto kernel_each = [&cells, &distance_function, &verlet_criterion,
                            &id_to_index, &verlet_list, max_id](int i) {
          auto &local_particles = cells[i]->particles();
          for (auto it = local_particles.begin(); it != local_particles.end();
               ++it) {
            auto &p1 = *it;
            if (p1.id() > max_id)
              continue;
            int ii = id_to_index(p1.id());
            if (ii < 0)
              continue;
            /* Pairs in this cell */
            for (auto jt = std::next(it); jt != local_particles.end(); ++jt) {
              if ((*jt).id() > max_id)
                continue;
              if (verlet_criterion(p1, *jt, distance_function(p1, *jt))) {
                int jj = id_to_index((*jt).id());
                if (jj >= 0) {
                  verlet_list.addNeighborNonAtomic(ii, jj);
                }
                // first_neighbor_kernel(ii, jj);
              }
            }
          }
        };

        auto kernel_neighbor = [&cells, &distance_function, &verlet_criterion,
                                &id_to_index, &verlet_list, max_id](int i) {
          auto &local_particles = cells[i]->particles();
          for (auto it = local_particles.begin(); it != local_particles.end();
               ++it) {
            auto const &p1 = *it;
            if (p1.id() > max_id)
              continue;
            int ii = id_to_index(p1.id());
            if (ii < 0)
              continue;
            /* Pairs with neighbors */
            for (auto &neighbor : cells[i]->neighbors().red()) {
              for (auto const &p2 : neighbor->particles()) {
                if (p2.id() > max_id)
                  continue;
                if (verlet_criterion(p1, p2, distance_function(p1, p2))) {
                  int jj = id_to_index(p2.id());
                  if (jj >= 0) {
                    verlet_list.addNeighbor(ii, jj);
                  }
                  // first_neighbor_kernel(ii, jj);
                }
              }
            }
          }
        };

        Kokkos::parallel_for("each", cells.size(), kernel_each);
        Kokkos::fence();

        Kokkos::parallel_for("neighbor", cells.size(), kernel_neighbor);
        Kokkos::fence();

        // verlet_list.get_max_counts();

        // Save data for next iteration if we just rebuilt
        CabanaData new_data(verlet_list, unique_particles, max_id);
        cell_structure.set_cabana_data(std::make_unique<CabanaData>(new_data));
#ifdef CALIPER
        CALI_MARK_END("Cabana - Verlet List");
#endif
      } else if (not rebuild) {
        //  Else use the saved verlet list
        verlet_list = saved_data.get_verlet_list();
      }
    } // else {
    {
#ifdef CALIPER
      CALI_MARK_BEGIN("Cabana - calc Force");
#endif
      FirstNeighborKernel first_neighbor_kernel_o(
#if defined(EXCLUSIONS) or defined(THOLE) or defined(ELECTROSTATICS) or        \
    defined(P3M) or defined(DPD) or defined(DIPOLES) or defined(NPT)
          unique_particles,
#endif
          bonded_ias, nonbonded_ias, box_geo, local_force,
#ifdef ROTATION
          local_torque,
#endif
#ifdef NPT
          local_virial,
#endif
#ifdef COLLISION_DETECTION
          *collision_detection,
#endif
          coulomb_kernel,
#if defined(THOLE) or defined(ELECTROSTATICS) or defined(P3M) or               \
    defined(DPD) or defined(DIPOLES) or defined(NPT)
          dipoles_kernel, elc_kernel, coulomb_u_kernel, thermostat,
#endif
          aosoa); // num_threads, rank, number_of_unique_particles);

      const auto &first_neighbor_kernel = first_neighbor_kernel_o;
      // std::vector<std::pair<int, int>> interaction_pairs;
      // std::vector<std::pair<Particle *, Particle *>> interaction_pairs;
      /*using neighbor_list = Cabana::NeighborList<ListType>;
      if (rank == 3) {
      for (int i = 0; i < number_of_unique_particles; ++i) {
        std::cout << "i:" << i;
        for (int n = 0; n < neighbor_list::numNeighbor(verlet_list, i); ++n) {
          int j = neighbor_list::getNeighbor(verlet_list, i, n);
          std::cout << " " << j;
          //first_neighbor_kernel(i, j);
          //interaction_pairs.emplace_back(i, j);
          //interaction_pairs.emplace_back(unique_particles.at(i),
          //		  		 unique_particles.at(j));
        }
        std::cout << std::endl;
      }
      }*/
      //
      // verlet_list.get_max_counts();
      //
      Kokkos::RangePolicy<execution_space> policy(0, particle_storage.size());
      Cabana::neighbor_parallel_for(policy, first_neighbor_kernel, verlet_list,
                                    Cabana::FirstNeighborsTag(),
                                    Cabana::TeamOpTag());
      //
      /*
      using neighbor_list = Cabana::NeighborList<ListType>;
      using SimdPolicy = Cabana::SimdPolicy<vector_length, execution_space>;
      SimdPolicy simd_policy(0, number_of_unique_particles);
      Cabana::simd_parallel_for(simd_policy,
                      [&number_of_unique_particles, &verlet_list,
                       &first_neighbor_kernel] (const int s, const int a) {
                        int i = s * vector_length + a;
                        if (i > number_of_unique_particles) return;
                        for (int n = 0; n <
      neighbor_list::numNeighbor(verlet_list, i); ++n) { int j =
      neighbor_list::getNeighbor(verlet_list, i, n); first_neighbor_kernel(i,
      j);
                        }
                      });
                      */

      Kokkos::fence();
#ifdef CALIPER
      CALI_MARK_END("Cabana - calc Force");
#endif
    }

#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - reduction Forces");
#endif
    // Force and Torque reduction
    Kokkos::RangePolicy<execution_space> policy(0, particle_storage.size());
    Kokkos::parallel_for("reduction", policy,
                         [&local_force,
#ifdef ROTATION
                          &local_torque,
#endif
                          &unique_particles, num_threads](const int i) {
                           double fx = 0.;
                           double fy = 0.;
                           double fz = 0.;
#ifdef ROTATION
                           double tx = 0.;
                           double ty = 0.;
                           double tz = 0.;
#endif
                           for (int tid = 0; tid < num_threads; ++tid) {
                             fx += local_force(i, tid, 0);
                             fy += local_force(i, tid, 1);
                             fz += local_force(i, tid, 2);
#ifdef ROTATION
                             tx += local_torque(i, tid, 0);
                             ty += local_torque(i, tid, 1);
                             tz += local_torque(i, tid, 2);
#endif
                           }
                           auto &p = unique_particles.at(i);
                           // auto p =
                           // cell_structure.get_local_particle(aosoa.id(i));
                           p->force() += Utils::Vector3d{fx, fy, fz};
#ifdef ROTATION
                           p->torque() += Utils::Vector3d{tx, ty, tz};
#endif
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
    CALI_MARK_END("Cabana - reduction Forces");
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
  }
}

#endif
