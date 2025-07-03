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
#include "verlet_list_loop.hpp"
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

inline void write_particle(Particle const &p, int const &id, AoSoA_pack &aosoa,
                           Utils::Vector3d &box_l) {
  aosoa.id(id) = p.id();
  aosoa.charge(id) = p.q();
  aosoa.type(id) = p.type();
  aosoa.ghost(id) = p.is_ghost();
  auto const pos = p.pos();
  double wpos[3] = {};
  for (int d = 0; d < 3; ++d) {
    // aosoa.position(id, d) =
    //     pos[d] - std::floor(pos[d] / box_l[d]) * box_l[d];
    wpos[d] = pos[d] - std::floor(pos[d] / box_l[d]) * box_l[d];
  }
  for (int d = 0; d < 3; ++d) {
    aosoa.position(id, d) = wpos[d];
  }
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
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // ===================================================
    // Setup Cabana Variables
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Setup");
#endif

    using ListAlgorithm = Cabana::HalfNeighborTag;
    using ListType = Cabana::CustomVerletList<memory_space, ListAlgorithm,
                                              Cabana::VerletLayout2D>;

    // Number of threads
    int num_threads = execution_space().concurrency();

#ifdef CALIPER
    CALI_MARK_END("Cabana - Setup");
#endif

    // ===================================================
    // Count unique particles and create Index map
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Index map");
#endif
    std::unordered_set<int> registered_index{};
    std::vector<Particle *> unique_particles;
    // std::vector<Particle> sequential_particles;
    int index = 0;
    // int max_id = 0;

    bool const rebuild = cell_structure.get_rebuild_cabana_verlet_list();
    // if (rank == 0) {
    //   std::cout << "For CABANA rebuild " << rebuild
    //	        << " " << Kokkos::OpenMP::concurrency() << std::endl;
    // }

    CabanaData saved_data;

    // Load saved data if we do not have to rebuild
    if (!rebuild) {
      saved_data = cell_structure.get_cabana_data();
    }

    // If we have to rebuild, we need to count the particles and create a new
    // map
    if (rebuild) {

      for (auto &p : particles) {
        if (cell_structure.get_local_particle(p.id())) {
          // if (p.id() > max_id) max_id = p.id();
          registered_index.insert(p.id());
          unique_particles.emplace_back(&p);
          // sequential_particles.emplace_back(p);
          index++;
        }
      }

      for (auto &p : ghost_particles) {
        if (not registered_index.contains(p.id())) {
          if (cell_structure.get_local_particle(p.id())) {
            // if (p.id() > max_id) max_id = p.id();
            registered_index.insert(p.id());
            unique_particles.emplace_back(&p);
            // sequential_particles.emplace_back(p);
            index++;
          }
        }
      }
    } else {
      // If we do not rebuild we can use the saved map
      index = saved_data.get_index();
      unique_particles = saved_data.get_unique_particles();
      // max_id = saved_data.get_max_id();
    }

    int number_of_unique_particles = index;
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
    particle_storage.resize(number_of_unique_particles);
    auto slice_position = Cabana::slice<0>(particle_storage);
    auto slice_charge = Cabana::slice<1>(particle_storage);
    auto slice_id = Cabana::slice<2>(particle_storage);
    auto slice_type = Cabana::slice<3>(particle_storage);
    auto slice_ghost = Cabana::slice<4>(particle_storage);
    // particle properties are defined in aosoa_pack.hpp
    auto aosoa = AoSoA_pack(particle_storage);
    auto box_l = box_geo.length();
    /*
    using policy_type = Kokkos::RangePolicy<execution_space>;
    Kokkos::parallel_for("AoSoA write", policy_type(0, particle_storage.size()),
                         [&unique_particles, &aosoa, &box_l](
                             const int p_id) {

                           Particle p = *unique_particles.at(p_id);

                           write_particle(p, p_id, aosoa, box_l);
                         });
    */

    using policy_type = Kokkos::TeamPolicy<execution_space>;
    int league_size = unique_particles.size();
    Kokkos::parallel_for(
        "AoSoA Write", policy_type(league_size, Kokkos::AUTO),
        [&unique_particles, &aosoa, &box_l, number_of_unique_particles](
            const policy_type::member_type &team_member) {
          int p_id = team_member.league_rank();
          if (p_id >= number_of_unique_particles)
            return;

          Particle p = *unique_particles.at(p_id);

          write_particle(p, p_id, aosoa, box_l);
        });

    /*using policy_type = Cabana::SimdPolicy<vector_length, execution_space>;
    //int league_size = (particle_storage.size() + v_length - 1) / v_length;
    //Kokkos::parallel_for("SIMD AoSoA Write",
    Cabana::simd_parallel_for(
                         policy_type(0, particle_storage.size()),
                         [&unique_particles, &box_l, vector_length,
                          number_of_unique_particles,
                          &slice_position, &slice_charge, &slice_id,
                          &slice_type, &slice_ghost] (const int s, const int a)
    {

                           int p_id = s * vector_length + a;
                           if (p_id >= number_of_unique_particles) return;

                           Particle p = *unique_particles.at(p_id);

                           auto pos = p.pos();
                           for (int d = 0; d < 3; ++d) {
                             double wrapped = pos[d] - std::floor(pos[d] /
    box_l[d]) * box_l[d]; slice_position.access(s, a, d) = wrapped;
                           }

                           slice_charge.access(s, a) = p.q();         // charge
                           slice_id.access(s, a)     = p.id();        // id
                           slice_type.access(s, a)   = p.type();      // type
                           slice_ghost.access(s, a)  = p.is_ghost();  // ghost
                         }, "SIMD AoSoA Write");*/
    Kokkos::fence();

    Kokkos::View<double ***, Kokkos::LayoutRight> local_force(
        "local_force", num_threads, number_of_unique_particles, 3);

    Kokkos::View<double ***, Kokkos::LayoutRight> local_torque(
        "local_torque", num_threads, number_of_unique_particles, 3);

    Kokkos::View<double **, Kokkos::LayoutRight> local_virial("local_virial",
                                                              num_threads, 3);

#ifdef CALIPER
    CALI_MARK_END("Cabana - Allocation");
#endif

    // The kernel of calculate force
    struct FirstNeighborKernel {
#if defined(EXCLUSIONS) or defined(THOLE) or defined(ELECTROSTATICS) or        \
    defined(P3M) or defined(DPD) or defined(DIPOLES) or defined(NPT)
      std::vector<Particle *> unique_particles;
#endif
      [[maybe_unused]] const BondedInteractionsMap &bonded_ias;
      const InteractionsNonBonded &nonbonded_ias;
      const BoxGeometry &box_geo;
      AoSoA_pack aosoa;
      Kokkos::View<double ***> local_force;
#ifdef ROTATION
      Kokkos::View<double ***> local_torque;
#endif
#ifdef NPT
      Kokkos::View<double **> local_virial;
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

      int num_threads;
      int mpi_rank;
      int particle_number;

      FirstNeighborKernel(
      // const CellStructure *cell_,
#if defined(EXCLUSIONS) or defined(THOLE) or defined(ELECTROSTATICS) or        \
    defined(P3M) or defined(DPD) or defined(DIPOLES) or defined(NPT)
          std::vector<Particle *> &unique_particles_,
#endif
          [[maybe_unused]] const BondedInteractionsMap &bonded_ias_,
          const InteractionsNonBonded &nonbonded_ias_,
          const BoxGeometry &box_geo_, AoSoA_pack &aosoa_,
          Kokkos::View<double ***> local_force_,
#ifdef ROTATION
          Kokkos::View<double ***> local_torque_,
#endif
#ifdef NPT
          Kokkos::View<double **> local_virial_,
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
          int &num_threads_, int &mpi_rank_, int &particle_number_)
          : // cell(cell_),
#if defined(EXCLUSIONS) or defined(THOLE) or defined(ELECTROSTATICS) or        \
    defined(P3M) or defined(DPD) or defined(DIPOLES) or defined(NPT)
            unique_particles(unique_particles_),
#endif
            bonded_ias(bonded_ias_), nonbonded_ias(nonbonded_ias_),
            box_geo(box_geo_), aosoa(aosoa_), local_force(local_force_),
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
            num_threads(num_threads_), mpi_rank(mpi_rank_),
            particle_number(particle_number_) {
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
        Utils::Vector3d const pi = {aosoa.position(i, 0), aosoa.position(i, 1),
                                    aosoa.position(i, 2)};
        Utils::Vector3d const pj = {aosoa.position(j, 0), aosoa.position(j, 1),
                                    aosoa.position(j, 2)};

        Utils::Vector3d const d = box_geo.get_mi_vector(pi, pj);
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
    CALI_MARK_BEGIN("Cabana - Verlet List by ESPRESSO");
#endif
    ListType verlet_list;

    // Rebuild verlet list if needed
    bool at_steepest_descent = cell_structure.get_steepest_descent_flag();
    int max_counts;
    double max_cutoff = pair_cutoff; // system.get_interaction_range();
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
    if (rebuild) { // Legacy Velert List
      /*verlet_list =
          ListType(aosoa.position, 0, aosoa.position.size(), max_counts,
      num_threads); auto kernel = [&](Particle const &p1, Particle const &p2) {
        auto thread_id = omp_get_thread_num();
        verlet_list.addNeighbor(thread_id,
                                id_to_index(p1.id()),
                                id_to_index(p2.id()));
        //std::cout << "WITHSMP "
                  //<< id_to_index(p1.id()) << " "
                  //<< id_to_index(p2.id()) << " "
                  //<< p1.is_ghost() << " "
                  //<< p2.is_ghost() << " "
                  //<< p1.id() << " "
                  //<< p2.id() << std::endl;
                    //<< p1.pos() << " "
                    //<< p2.pos() << "\n";
        };

      cell_structure.cabana_verlet_list_loop(kernel, verlet_criterion);
      */
    } else {
      // Else use the saved verlet list
      verlet_list = saved_data.get_verlet_list();
    }
#ifdef CALIPER
    CALI_MARK_END("Cabana - Verlet List by ESPRESSO");
#endif

    FirstNeighborKernel first_neighbor_kernel(
#if defined(EXCLUSIONS) or defined(THOLE) or defined(ELECTROSTATICS) or        \
    defined(P3M) or defined(DPD) or defined(DIPOLES) or defined(NPT)
        unique_particles,
#endif
        bonded_ias, nonbonded_ias, box_geo, aosoa, local_force,
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
        num_threads, rank, number_of_unique_particles);

    if (rebuild and max_cutoff != INACTIVE_CUTOFF) { // Shared memory
#ifdef CALIPER
      CALI_MARK_BEGIN("Cabana - Verlet List by Cabana");
#endif
      if (1) {
        verlet_list = create_verlet_list(
            max_cutoff, max_counts, aosoa, unique_particles, verlet_criterion,
            // sequential_particles, verlet_criterion,
            first_neighbor_kernel, cell_structure);
        /*using neighbor_list = Cabana::NeighborList<ListType>;
        std::vector<std::pair<int, int>> interaction_pairs;

        for (int i = 0; i < number_of_unique_particles; ++i) {
          for (int n = 0; n < neighbor_list::numNeighbor(verlet_list, i); ++n) {
            int j = neighbor_list::getNeighbor(verlet_list, i, n);
            //interaction_pairs.emplace_back(i, j);
            std::cout << "*Cabana* "
                          << i << " "
                          << j << " "
                          << aosoa.ghost(i) << " "
                          << aosoa.ghost(j) << " "
                          << aosoa.id(i) << " "
                          << aosoa.id(j) << "\n";
          }
        }*/
        // verlet_list.get_max_counts();
      }
#ifdef CALIPER
      CALI_MARK_END("Cabana - Verlet List by Cabana");
#endif
    } else {
      //{
#ifdef CALIPER
      CALI_MARK_BEGIN("Cabana - calc Force");
#endif
      /*
      Kokkos::parallel_for(
        "ForceLoop", Kokkos::RangePolicy<>(0,interaction_pairs.size()),
        KOKKOS_LAMBDA(int idx) {
          auto i = interaction_pairs[idx].first;
          auto j = interaction_pairs[idx].second;
          first_neighbor_kernel(i, j);
      });
      */
      // verlet_list.get_max_counts();

      Kokkos::RangePolicy<execution_space> policy(0, particle_storage.size());
      Cabana::neighbor_parallel_for(policy, first_neighbor_kernel, verlet_list,
                                    Cabana::FirstNeighborsTag(),
                                    Cabana::TeamOpTag());

      Kokkos::fence();
#ifdef CALIPER
      CALI_MARK_END("Cabana - calc Force");
#endif
    }

    // Save data for next iteration if we just rebuilt
    if (rebuild) {
      // CabanaData new_data(verlet_list, unique_particles,
      //                     unique_particles.size());
      CabanaData new_data(verlet_list, unique_particles);
      cell_structure.set_cabana_data(std::make_unique<CabanaData>(new_data));
    }

#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - reduction Forces");
#endif
    // Force and Torque reduction
    Kokkos::RangePolicy<execution_space> policy(0, particle_storage.size());
    Kokkos::parallel_for("reduction", policy,
                         [&local_force, &local_torque, &unique_particles,
                          num_threads](const int i) {
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
