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

#include "aosoa_pack.hpp"
// #include "cabana_data.hpp"
#include "custom_verlet_list.hpp"
#include "forces_cabana.hpp"
#include <Cabana_Core.hpp>
#include <Cabana_NeighborList.hpp>
#include <iostream>

inline void write_particle(Particle const &p, int const &id,
                           AoSoA_pack &aosoa) {
  aosoa.id(id) = p.id();
  aosoa.charge(id) = p.q();
  aosoa.type(id) = p.type();
  auto const pos = p.pos();
  for (int d = 0; d < 3; ++d) {
    aosoa.position(id, d) = pos[d];
  }
}

using ListAlgorithm = Cabana::HalfNeighborTag;
using ListType = Cabana::CustomVerletList<memory_space, ListAlgorithm,
                                          Cabana::VerletLayout2D>;
template <class VerletCriterion>
__attribute__((always_inline)) inline void construct_verlet_list(
    CellStructure &cell_structure, VerletCriterion const &verlet_criterion,
    Kokkos::View<int *> const &id_to_index, const int max_id) {
  auto const &cells =
      std::as_const(cell_structure).decomposition().local_cells();
  auto const distance_function = detail::MinimalImageDistance{
      std::as_const(cell_structure).decomposition().box()};
  auto verlet_list = cell_structure.get_cabana_verlet_list();

  auto intra_kernel = [&cells, &distance_function, &verlet_criterion,
                       &id_to_index, &verlet_list, max_id](const int i) {
    auto &local_particles = cells[i]->particles();
    for (auto it = local_particles.begin(); it != local_particles.end(); ++it) {
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
            verlet_list.addNeighborLoadBalancing(ii, jj);
          }
        }
      }
    }
  };

  auto inter_kernel = [&cells, &distance_function, &verlet_criterion,
                       &id_to_index, &verlet_list, max_id](const int i) {
    auto &local_particles = cells[i]->particles();
    for (auto it = local_particles.begin(); it != local_particles.end(); ++it) {
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
              // verlet_list.addNeighborAtomic(ii, jj);
              verlet_list.addNeighborNonAtomic(ii, jj);
            }
          }
        }
      }
    }
  };

  Kokkos::parallel_for("inter", cells.size(), intra_kernel);
  Kokkos::fence();

  Kokkos::parallel_for("intra", cells.size(), inter_kernel);
  Kokkos::fence();

  // verlet_list.sortNeighbors();
}

template <class VerletCriterion>
__attribute__((always_inline)) inline void update_cabana_state(
    CellStructure &cell_structure, ParticleRange const &particles,
    ParticleRange const &ghost_particles,
    VerletCriterion const &verlet_criterion, double const pair_cutoff) {
#ifdef CALIPER
  CALI_MARK_BEGIN("Cabana - Index map");
#endif
  // Number of threads
  int num_threads = execution_space().concurrency();

  int number_of_unique_particles = 0;

  bool const rebuild = cell_structure.get_rebuild_cabana_verlet_list() or
                       (not cell_structure.use_verlet_list);
  // std::cout << "rebuild:" << rebuild << std::endl;

  if (rebuild) {
    // If we have to rebuild, we need to count the particles
    cell_structure.set_index_map(particles, ghost_particles,
                                 number_of_unique_particles);
    // Create essential variable for MD
    cell_structure.rebuild_local_properties(number_of_unique_particles,
                                            num_threads, pair_cutoff);
  } else {
    // If we do not rebuild we can use the saved map
    number_of_unique_particles = cell_structure.get_unique_particles().size();
    cell_structure.reset_local_properties();
  }
  auto const unique_particles = cell_structure.get_unique_particles();
  auto aosoa = cell_structure.get_aosoa_data();
  int max_id = cell_structure.get_max_id();

#ifdef CALIPER
  CALI_MARK_END("Cabana - Index map");
#endif
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
    Kokkos::deep_copy(id_to_index, -1);

    using policy_type = Kokkos::RangePolicy<execution_space>;
    Kokkos::parallel_for(
        "AoSoA write", policy_type(0, number_of_unique_particles),
        [&unique_particles, &aosoa, &id_to_index](const int p_id) {
          write_particle(*unique_particles.at(p_id), p_id, aosoa);
          id_to_index(unique_particles.at(p_id)->id()) = p_id;
        });
    Kokkos::fence();
#ifdef CALIPER
    CALI_MARK_END("Cabana - Allocation");
#endif

    // ===================================================
    // Get Verlet Pairs and Fill Verlet list
    // ===================================================

    // Rebuild verlet list if needed
    if (rebuild) {
#ifdef CALIPER
      CALI_MARK_BEGIN("Cabana - Verlet List");
#endif
      construct_verlet_list(cell_structure, verlet_criterion, id_to_index,
                            max_id);
      cell_structure.mark_rebuild_cabana_verlet_list_as_UpToDate();
#ifdef CALIPER
      CALI_MARK_END("Cabana - Verlet List");
#endif
    }
  }
}

template <class BondKernel, class PairKernel,
          class VerletCriterion = detail::True>
void cabana_short_range(
    BondKernel const &bond_kernel, PairKernel const &forces_kernel,
#ifdef COLLISION_DETECTION
    std::shared_ptr<CollisionDetection::CollisionDetection> collision_detection,
#endif
    CellStructure &cell_structure, double pair_cutoff, double bond_cutoff,
    ParticleRange const &particles, ParticleRange const &ghost_particles,
    VerletCriterion const &verlet_criterion = {}) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

  int num_threads = execution_space().concurrency();

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
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - calc Force");
#endif
    auto unique_particles = cell_structure.get_unique_particles();
    auto local_force = cell_structure.get_local_force();
#ifdef ROTATION
    auto local_torque = cell_structure.get_local_torque();
#endif
#ifdef NPT
    auto local_virial = cell_structure.get_local_virial();
#endif
    auto cabana_verlet_list = cell_structure.get_cabana_verlet_list();

    // cabana_verlet_list.get_variance_max_counts();
    Kokkos::RangePolicy<execution_space> policy(0, unique_particles.size());
    Cabana::neighbor_parallel_for(policy, forces_kernel, cabana_verlet_list,
                                  Cabana::FirstNeighborsTag(),
                                  // Cabana::TeamOpTag());
                                  Cabana::SerialOpTag());
    Kokkos::fence();
#ifdef CALIPER
    CALI_MARK_END("Cabana - calc Force");
#endif

#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - reduction Forces");
#endif
    // Force and Torque reduction
    // Kokkos::RangePolicy<execution_space> policy(0, unique_particles.size());
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
                           // auto &p = unique_particles.at(i);
                           // p->force() += Utils::Vector3d{fx, fy, fz};
                           unique_particles.at(i)->force() +=
                               Utils::Vector3d{fx, fy, fz};
#ifdef ROTATION
                           // p->torque() += Utils::Vector3d{tx, ty, tz};
                           unique_particles.at(i)->torque() +=
                               Utils::Vector3d{tx, ty, tz};
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
