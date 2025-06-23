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
  if (result >= L)
    result -= std::nextafter(L, 0.);
  return result;
}

inline void write_particle(Particle const &p, int const &id, AoSoA_pack &aosoa,
                           Utils::Vector3d &box_l) {
  auto const &pos = p.pos();
  aosoa.position(id, 0) = wrap(pos[0], box_l[0]);
  aosoa.position(id, 1) = wrap(pos[1], box_l[1]);
  aosoa.position(id, 2) = wrap(pos[2], box_l[2]);
  aosoa.id(id) = p.id();
  aosoa.charge(id) = p.q();
  aosoa.type(id) = p.type();
  aosoa.ghost(id) = p.is_ghost();
  aosoa.force(id, 0) = 0.0;
  aosoa.force(id, 1) = 0.0;
  aosoa.force(id, 2) = 0.0;
  aosoa.torque(id, 0) = 0.0;
  aosoa.torque(id, 1) = 0.0;
  aosoa.torque(id, 2) = 0.0;
  assert(aosoa.position(id, 0) >= 0. and aosoa.position(id, 0) < box_l[0]);
  assert(aosoa.position(id, 1) >= 0. and aosoa.position(id, 1) < box_l[1]);
  assert(aosoa.position(id, 2) >= 0. and aosoa.position(id, 2) < box_l[2]);
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

    // const int vector_length = 1;
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
    // std::vector<int> index_to_id{};
    std::vector<Particle *> unique_particles;
    int index = 0;

    bool const rebuild = cell_structure.get_rebuild_cabana_verlet_list();
    // std::cout << "For CABANA rebuild " << rebuild << " " <<
    // Kokkos::OpenMP::concurrency() << std::endl;

    CabanaData saved_data;

    // Load saved data if we do not have to rebuild
    if (!rebuild) {
      saved_data = cell_structure.get_cabana_data();
    }

    // If we have to rebuild, we need to count the particles and create a new
    // map
    if (rebuild) {

      for (auto &p : particles) {
        //        if (cell_structure.get_local_particle(p.id())) {
        // id_to_index[p.id()] = index;
        unique_particles.emplace_back(&p);
        index++;
        //        }
      }

      for (auto &p : ghost_particles) {
        if (not registered_index.contains(p.id())) {
          if (cell_structure.get_local_particle(p.id())->is_ghost()) {
            // id_to_index[p.id()] = index;
            registered_index.insert(p.id());
            unique_particles.emplace_back(&p);
            index++;
          }
        }
      }
    } else {
      // If we do not rebuild we can use the saved map
      // id_to_index = saved_data.get_id_to_index();
      index = saved_data.get_index();
      unique_particles = saved_data.get_unique_particles();
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
    auto aosoa = AoSoA_pack(
        particle_storage); // particle properties are defined in aosoa_pack.hpp
    auto box_l = box_geo.length();
    // int p_id = 0;
    // registered_index.clear();
    // Kokkos::View<Particle*> particle_view("particle_pointer",
    // unique_particles.size());
    Kokkos::RangePolicy<execution_space> allocation_policy(
        0, unique_particles.size());
    Kokkos::parallel_for("allocation", allocation_policy, [&](int p_id) {
      const auto &p = *unique_particles[p_id];
      // if (!cell_structure.get_local_particle(p.id())) continue;
      write_particle(p, p_id, aosoa, box_l);
    });
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
      // std::vector<int> &index_to_id;
      Kokkos::View<double ***> local_force;
      Kokkos::View<double ***> local_torque;
      Kokkos::View<double **> local_virial;
      AoSoA_pack aosoa;
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
          const BoxGeometry &box_geo_,
          // std::vector<int> &index_to_id_,
          Kokkos::View<double ***> local_force_,
          Kokkos::View<double ***> local_torque_,
          Kokkos::View<double **> local_virial_, AoSoA_pack &aosoa_,
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
            box_geo(box_geo_),
            // index_to_id(index_to_id_),
            local_force(local_force_), local_torque(local_torque_),
            local_virial(local_virial_), aosoa(aosoa_),
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
        /*
        auto p1 = unique_particles.at(i); //
cell->get_local_particle(aosoa.id(i)); auto p2 = unique_particles.at(j);//
cell->get_local_particle(aosoa.id(j));

        Utils::Vector3d const d = box_geo.get_mi_vector(p1->pos(), p2->pos());
        auto const dist = d.norm();

        auto const q1q2 =aosoa.charge(i) *aosoa.charge(j);

        auto const dist2 = dist * dist;
#ifdef NPT
        auto [pf, virial]
#else
        auto pf
#endif
        = add_non_bonded_pair_force(
            const_cast<Particle &>(*p1), const_cast<Particle &>(*p2), d, dist,
            dist2, q1q2, ia_params, thermostat, box_geo, bonded_ias,
            coulomb_kernel, dipoles_kernel, elc_kernel, coulomb_u_kernel);
        */

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
        auto p1 =
            unique_particles.at(i); // cell->get_local_particle(aosoa.id(i));
        auto p2 =
            unique_particles.at(j); // cell->get_local_particle(aosoa.id(j));

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
        auto p1 =
            unique_particles.at(i); // cell->get_local_particle(aosoa.id(i));
        auto p2 =
            unique_particles.at(j); // cell->get_local_particle(aosoa.id(j));

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
    // auto const &system = ::System::get_system();
    int max_counts;
    double max_cutoff = pair_cutoff; // system.get_interaction_range();
    if (std::isinf(max_cutoff)) {
      max_counts = number_of_unique_particles;
    } else {
      max_counts = 64;
      static_cast<int>(27 * max_cutoff * max_cutoff * max_cutoff / 3);
    }
    //    if (max_counts < 256)
    //      max_counts = 56;
    if (rebuild) { // Legacy Velert List
      /*verlet_list =
          ListType(aosoa.position, 0,aosoa.position.size(), max_counts);
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
    CALI_MARK_END("Cabana - Verlet List by ESPRESSO");
#endif

    FirstNeighborKernel first_neighbor_kernel(
#if defined(EXCLUSIONS) or defined(THOLE) or defined(ELECTROSTATICS) or        \
    defined(P3M) or defined(DPD) or defined(DIPOLES) or defined(NPT)
        unique_particles,
#endif
        bonded_ias, nonbonded_ias, box_geo, local_force, local_torque,
        local_virial, aosoa,
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
      verlet_list = create_verlet_list(max_cutoff, max_counts, aosoa,
                                       unique_particles, verlet_criterion,
                                       first_neighbor_kernel, cell_structure);
#ifdef CALIPER
      CALI_MARK_END("Cabana - Verlet List by Cabana");
#endif
    } else {
      //{
#ifdef CALIPER
      CALI_MARK_BEGIN("Cabana - calc Force");
#endif
      /*
      using neighbor_list = Cabana::NeighborList<ListType>;
      std::vector<std::pair<int, int>> interaction_pairs;

      for (int i = 0; i < number_of_unique_particles; ++i) {
        for (int n = 0; n < neighbor_list::numNeighbor(verlet_list, i); ++n) {
          int j = neighbor_list::getNeighbor(verlet_list, i, n);
          interaction_pairs.emplace_back(i, j);
        }
      }
      */
      /*
      Kokkos::parallel_for("ForceLoop", Kokkos::RangePolicy<>(0,
      interaction_pairs.size()), KOKKOS_LAMBDA(int idx) { auto i =
      interaction_pairs[idx].first; auto j = interaction_pairs[idx].second;
                        first_neighbor_kernel(i, j);
                      });
      */

      Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>
          policy(0, particle_storage.size());
      Cabana::neighbor_parallel_for(policy, first_neighbor_kernel, verlet_list,
                                    Cabana::FirstNeighborsTag(),
                                    Cabana::SerialOpTag());

      Kokkos::fence();
#ifdef CALIPER
      CALI_MARK_END("Cabana - calc Force");
#endif
    }

    // Save data for next iteration if we just rebuilt
    if (rebuild) {
      CabanaData new_data(verlet_list, unique_particles,
                          particle_storage.size());
      cell_structure.set_cabana_data(std::make_unique<CabanaData>(new_data));
    }

#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - reduction Forces");
#endif
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
          aosoa.force(i, 0) = fx;
          aosoa.force(i, 1) = fy;
          aosoa.force(i, 2) = fz;
          aosoa.torque(i, 0) = tx;
          aosoa.torque(i, 1) = ty;
          aosoa.torque(i, 2) = tz;
        });
    Kokkos::parallel_for(
        "add_force", policy, KOKKOS_LAMBDA(const int i) {
          auto &p = unique_particles[i];
          p->force() += Utils::Vector3d{aosoa.force(i, 0), aosoa.force(i, 1),
                                        aosoa.force(i, 2)};
          p->torque() += Utils::Vector3d{aosoa.torque(i, 0), aosoa.torque(i, 1),
                                         aosoa.torque(i, 2)};
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

    // ===================================================
    // Add forces to particles
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Particle Forces");
#endif
/*
    for (auto id = 0; id < particle_storage.size(); ++id) {
      auto p = unique_particles[id];
      if (p == nullptr) {
        return;
      }
      Utils::Vector3d f_vec{aosoa.force(id, 0), aosoa.force(id, 1),
                            aosoa.force(id, 2)};
      Utils::Vector3d torque_vec{aosoa.torque(id, 0), aosoa.torque(id, 1),
                                 aosoa.torque(id, 2)};

#ifdef ROTATION
      ParticleForce f(f_vec, torque_vec);
#else
      ParticleForce f(f_vec);
#endif
      p->force_and_torque() += f;
    }
*/
#ifdef CALIPER
    CALI_MARK_END("Cabana - Particle Forces");
#endif
  }
}

#endif
