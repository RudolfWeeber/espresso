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
#include <unordered_set>
#include <utility>
#include <stdio.h>

template <class SliceDouble3, class SliceInt>
inline void write_particle(Particle const &p,
                           std::unordered_map<int, int> const &id_to_index,
                           SliceDouble3 &s_position, SliceDouble3 &s_force,
                           SliceDouble3 &s_torque, SliceInt &s_id,
                           SliceInt &s_type) {
  auto const pos = p.pos();
  auto const id = id_to_index.at(p.id());
  s_position(id, 0) = pos[0];
  s_position(id, 1) = pos[1];
  s_position(id, 2) = pos[2];
  s_id(id) = p.id();
  s_type(id) = p.type();
  s_force(id, 0) = 0.0;
  s_force(id, 1) = 0.0;
  s_force(id, 2) = 0.0;
  s_torque(id, 0) = 0.0;
  s_torque(id, 1) = 0.0;
  s_torque(id, 2) = 0.0;
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
    using data_types =
        Cabana::MemberTypes<double[3], double[3], double[3], int, int, int>;
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
        id_to_index[p.id()] = index;
        index_to_id.emplace_back(p.id());
        index++;
      }

      for (auto const &p : ghost_particles) {
        if (not id_to_index.contains(p.id())) {
          id_to_index[p.id()] = index;
          index_to_id.emplace_back(p.id());
          index++;
        }
      }
    } else {
      // If we do not rebuild we can use the saved map
      id_to_index = saved_data.get_id_to_index();
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
    auto slice_id = Cabana::slice<3>(particle_storage);
    auto slice_type = Cabana::slice<4>(particle_storage);
    for (auto const &p : particles) {
      write_particle(p, id_to_index, slice_position, slice_force, slice_torque,
                     slice_id, slice_type);
    }
    //using TP = decltype(slice_position);
    //using TF = decltype(slice_force);
    //using TR = decltype(slice_torque);
    using TT = decltype(slice_type);

    Kokkos::View<double[3], Kokkos::DefaultExecutionSpace> virial_all(
        "virial_all");

    Kokkos::View<double ***, Kokkos::LayoutRight> local_force(
        "local_force", num_threads, number_of_unique_particles, 3);

    Kokkos::View<double ***, Kokkos::LayoutRight> local_torque(
        "local_torque", num_threads, number_of_unique_particles, 3);

    Kokkos::View<double **, Kokkos::LayoutRight> local_virial("local_virial", num_threads, 3);

    for (auto const &p : ghost_particles) {
      // if the ghost is not in the previous map, but mpi moved it to this rank?
      // it will not have neighbors because we did not rebuild the verlet list.
      if (not id_to_index.contains(p.id())) {
        continue;
      }
      write_particle(p, id_to_index, slice_position, slice_force, slice_torque,
                     slice_id, slice_type);
    }
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

    // Rebuild verlet list if needed
    if (rebuild) {

      verlet_list = ListType(slice_position, 0, slice_position.size(), 32);

      auto kernel = [&](Particle const &p1, Particle const &p2) {
        verlet_list.addNeighbor(id_to_index.at(p1.id()),
                                id_to_index.at(p2.id()));
      };

      cell_structure.cabana_verlet_list_loop(kernel, verlet_criterion);
    } else {
      // Else use the saved verlet list
      verlet_list = saved_data.get_verlet_list();
    }

    // Save data for next iteration if we just rebuilt
    if (rebuild) {
      CabanaData new_data(verlet_list, id_to_index, index_to_id);
      cell_structure.set_cabana_data(std::make_unique<CabanaData>(new_data));
    }

    // fill customverletlist with pairs
    struct FirstNeighborKernel {
      const CellStructure *cell;
      [[maybe_unused]] const BondedInteractionsMap &bonded_ias;
      const InteractionsNonBonded &nonbonded_ias;
      const Thermostat::Thermostat &thermostat;
      const BoxGeometry &box_geo;
      std::vector<int> &index_to_id;
      Kokkos::View<double ***> local_force;
      Kokkos::View<double ***> local_torque;
      Kokkos::View<double **> local_virial;
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
          const BoxGeometry &box_geo_, std::vector<int> &index_to_id_,
          Kokkos::View<double ***> local_force_,
      	  Kokkos::View<double ***> local_torque_,
      	  Kokkos::View<double **> local_virial_,
	  TT &slice_type_,
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
            index_to_id(index_to_id_),
	    local_force(local_force_),
	    local_torque(local_torque_), local_virial(local_virial_),
	    slice_type(slice_type_),
            collision_detection(collision_detection_),
            coulomb_kernel(coulomb_kernel_), dipoles_kernel(dipoles_kernel_),
            elc_kernel(elc_kernel_), coulomb_u_kernel(coulomb_u_kernel_),
	    num_threads(num_threads_), mpi_rank(mpi_rank_) {
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(int i, int j) const {
	/*
        Utils::Vector3d const pi = {slice_position(i, 0), slice_position(i, 1),
                                    slice_position(i, 2)};
        Utils::Vector3d const pj = {slice_position(j, 0), slice_position(j, 1),
                                    slice_position(j, 2)};
	*/

        auto p1 = cell->get_local_particle(index_to_id.at(i));
        auto p2 = cell->get_local_particle(index_to_id.at(j));
        if (p1 == nullptr or p2 == nullptr)
          return;

        Utils::Vector3d const d = box_geo.get_mi_vector(p1->pos(), p2->pos());
        auto const dist = d.norm();
        auto const dist2 = dist * dist;

	auto thread_id = omp_get_thread_num();
        // auto thread_id = Kokkos::OpenMP::impl_hardware_thread_id();
        //std::cout << "in " << thread_id << " " << i << " " << j << "\n";

        IA_parameters const &ia_params =
            nonbonded_ias.get_ia_param(slice_type(i), slice_type(j));

	auto[pf, virial] = add_non_bonded_pair_force(
	  const_cast<Particle &>(*p1), const_cast<Particle &>(*p2),
	  d, dist, dist2, ia_params, thermostat, box_geo, bonded_ias,
    	  coulomb_kernel, dipoles_kernel, elc_kernel, coulomb_u_kernel);

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
        index_to_id,
	local_force, local_torque, local_virial,
	slice_type,
#ifdef COLLISION_DETECTION
        *collision_detection,
#endif
        coulomb_kernel, dipoles_kernel, elc_kernel, coulomb_u_kernel,
        num_threads, rank);
    /*
    // For using not custom_verletlist but Cabana::VeletList:
    using s_ListType = Cabana::VerletList<memory_space, ListAlgorithm,
                                              Cabana::VerletLayout2D>;
    auto box_l = box_geo.length();
    double grid_min[3] = { 0.0, 0.0, 0.0 };
    double grid_max[3] = { box_l[0], box_l[1], box_l[2] };
    s_ListType s_verlet_list;
    s_verlet_list = s_ListType(slice_position, 0, slice_position.size(), nonbonded_ias.maximal_cutoff(), 1.0, grid_min, grid_max);
    */

    //std::cout << rank << " " << num_threads << " Execute FirstNeighborKernel\n";
    //  TODO: Add option to switch "SerialOpTag" Between "TeamOpTag"
    //  Feels like TeamOpTag is faster, atleast for large particle numbers
    Cabana::neighbor_parallel_for(policy, first_neighbor_kernel, verlet_list,
                                  Cabana::FirstNeighborsTag(),
                                  Cabana::SerialOpTag());//, "verlet_list");
    /*
    // For checking how custom_verlet_list works:
    using TeamPolicy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>;
    using MemberType = TeamPolicy::member_type;

    TeamPolicy team_policy(number_of_unique_particles, 1);

    using neighbor_list_traits = Cabana::NeighborList<ListType>;
    using s_neighbor_list_traits = Cabana::NeighborList<s_ListType>;

    //std::cout << "Number of threads " << num_threads << "\n";
    Kokkos::parallel_for("force_calc_by_team", team_policy,
      KOKKOS_LAMBDA(const MemberType& team_member)
    {
	const int i = team_member.league_rank(); // particle index in verlet_list
	const int num_neighbors = neighbor_list_traits::numNeighbor(verlet_list, i);
	//const int s_num_neighbors = s_neighbor_list_traits::numNeighbor(s_verlet_list, i);
	//std::cout << "team_size " << team_member.team_size() << "\n";
	//std::cout << "neighbor " << i << " " << num_neighbors << " " << s_num_neighbors << "\n";

	Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, num_neighbors),
			[&](const int n) {
	  const int j = neighbor_list_traits::getNeighbor(verlet_list, i, n); // particle index in verlet_list
	  const int thread_id = omp_get_thread_num();
	  std::cout << "tid " << thread_id << "\n";
	  char region_name[64];
	  sprintf(region_name, "work_region_thread_%d", thread_id);
	  cali_begin_region(region_name);
	  first_neighbor_kernel(i, j);
	  cali_end_region(region_name);
	});
    });
    */
    Kokkos::fence();

    Kokkos::parallel_for("reduce", policy,
      KOKKOS_LAMBDA(const int i) {
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
      }
    );

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
    for (auto &p : particles) {
      auto const id = id_to_index.at(p.id());
      Utils::Vector3d f_vec{slice_force(id, 0), slice_force(id, 1),
                            slice_force(id, 2)};
      Utils::Vector3d torque_vec{slice_torque(id, 0), slice_torque(id, 1),
                                 slice_torque(id, 2)};

      ParticleForce f(f_vec, torque_vec);
      p.force_and_torque() += f;
    }

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
#ifdef CALIPER
    CALI_MARK_END("Cabana - Particle Forces");
#endif
  }
}

#endif
