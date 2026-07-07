/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include <config/config.hpp>

#include "cell_system/CellStructure.hpp"

#include "aosoa_pack.hpp"
#include "bond_forces_kokkos.hpp"
#include "custom_verlet_list.hpp"
#include "forces_cabana.hpp"

#include <Cabana_Core.hpp>
#include <Cabana_NeighborList.hpp>

#ifdef ESPRESSO_CALIPER
#include <caliper/cali.h>
#endif

#include <iterator>
#include <span>
#include <utility>

template <class KokkosRangePolicy = Kokkos::RangePolicy<>>
ESPRESSO_ATTR_ALWAYS_INLINE inline void
kokkos_parallel_range_for(auto const &name, auto start, auto end,
                          auto const &kernel) {
  if (Kokkos::num_threads() > 1) {
    KokkosRangePolicy policy(start, end);
    Kokkos::parallel_for(name, policy, kernel);
  } else {
    for (auto p_index = start; p_index < end; ++p_index) {
      kernel(p_index);
    }
  }
}

ESPRESSO_ATTR_ALWAYS_INLINE inline void
commit_particle(Particle const &p, auto const index,
                CellStructure::AoSoA_pack &aosoa, bool const rebuild) {
  // phase 3.5: position/image/director are no longer copied here. Kernels read
  // them directly from the ParticleStore columns (via the pack-index->store-row
  // translation view) and the store-side derived director view.
  // phase 4: velocity is no longer copied here either. It aliases the
  // ParticleStore velocity column; velocity-dependent kernels read it by *store
  // row* via row(i), just like position.
  // phase 5: id/mass are no longer copied here — they alias the ParticleStore
  // id/mass columns, read by *store row* via row(i) on the cold bond path.
  // charge/dipm are refreshed per step (when a solver is active) in
  // refresh_pack_charges / refresh_pack_dipm, not here.
  // phase-5 perf recovery: `type` is once again PACK-OWNED and written here on
  // rebuild (it is read pack-indexed by the hot pair kernels). The
  // ParticleStore type column stays authoritative; a mid-run type change forces
  // a rebuild (on_particle_type_change -> set_resort_particles), so this cache
  // is refreshed before it is next read.
  // phase 6: the has-exclusion flag (like `type`) is written ONLY on rebuild.
  // Exclusions cannot change between rebuilds: every add/delete goes through
  // on_particle_change -> set_resort_particles(RESORT_LOCAL), which forces the
  // next update_cabana_state onto the full-commit branch. The pack-owned
  // `flags` column persists across partial-update steps (per-step force reset
  // does not touch it; clear_local_properties zeroes it but also forces a
  // resort/rebuild). Reading the exclusions sidecar per particle per step was
  // pure waste -- and, post-eviction, that read is a ParticleStore column
  // indirection rather than the old inline member, so keeping it out of the
  // partial-update hot loop matters.
  if (rebuild) {
    aosoa.type(index) = p.type();
#ifdef ESPRESSO_EXCLUSIONS
    aosoa.set_has_exclusion(index, !p.exclusions().empty());
#else
    aosoa.flags(index) = 0;
#endif
  }
}

ESPRESSO_ATTR_ALWAYS_INLINE inline void
link_cell_kokkos(std::span<Cell *const> cells, BoxGeometry const &box_geo,
                 auto const &verlet_criterion,
                 Kokkos::View<int *> const &id_to_index, int const max_id,
                 auto const &intra_operator, auto const &inter_operator) {

  // implementation detail: max_id refers to the max local particle id,
  // but ghost particles from other ranks may have larger particle ids;
  // -1 is used as a sentinel value for particle ids from other threads

  auto intra_kernel = [&cells, &box_geo, &verlet_criterion, &id_to_index,
                       &intra_operator, max_id](const int i) {
    auto &local_particles = cells[i]->particles();
    for (auto it = local_particles.begin(); it != local_particles.end(); ++it) {
      auto const &p1 = *it;
      if (p1.id() <= max_id) {
        auto const ii = id_to_index(p1.id());
        if (ii >= 0) {
          // Hoist p1's position out of the inner loop: it lives in the
          // ParticleStore column (phase 3), so materializing the proxy once per
          // outer particle instead of once per pair candidate avoids O(pairs)
          // strided column reads on this hot path.
          auto const p1_pos = Utils::Vector3d(p1.pos());
          // pairs in this cell
          for (auto jt = std::next(it); jt != local_particles.end(); ++jt) {
            if ((*jt).id() <= max_id) {
              if (verlet_criterion(p1, *jt,
                                   box_geo.get_mi_dist2(
                                       p1_pos, Utils::Vector3d(jt->pos())))) {
                auto const jj = id_to_index((*jt).id());
                if (jj >= 0) {
                  intra_operator(ii, jj);
                }
              }
            }
          }
        }
      }
    }
  };

  auto inter_kernel = [&cells, &box_geo, &verlet_criterion, &id_to_index,
                       &inter_operator, max_id](const int i) {
    auto &local_particles = cells[i]->particles();
    for (auto const &p1 : local_particles) {
      if (p1.id() <= max_id) {
        auto const ii = id_to_index(p1.id());
        if (ii >= 0) {
          // Hoist p1's position out of the inner loops (see intra_kernel).
          auto const p1_pos = Utils::Vector3d(p1.pos());
          // pairs with neighboring cells
          for (auto &neighbor : cells[i]->neighbors().red()) {
            for (auto const &p2 : neighbor->particles()) {
              if (p2.id() <= max_id) {
                if (verlet_criterion(p1, p2,
                                     box_geo.get_mi_dist2(
                                         p1_pos, Utils::Vector3d(p2.pos())))) {
                  auto const jj = id_to_index(p2.id());
                  if (jj >= 0) {
                    inter_operator(ii, jj);
                  }
                }
              }
            }
          }
        }
      }
    }
  };

  Kokkos::parallel_for("intra", cells.size(), intra_kernel);
  Kokkos::fence();

  Kokkos::parallel_for("inter", cells.size(), inter_kernel);
  Kokkos::fence();
}

ESPRESSO_ATTR_ALWAYS_INLINE inline void
update_cabana_state(CellStructure &cell_structure, auto const &verlet_criterion,
                    double const pair_cutoff, auto const integ_switch) {
#ifdef ESPRESSO_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  using execution_space = Kokkos::DefaultExecutionSpace;
  using policy_type = Kokkos::RangePolicy<execution_space>;
  auto const rebuild = cell_structure.prepare_verlet_list_cabana(pair_cutoff);
  auto const &unique_particles = cell_structure.get_unique_particles();
  auto const n_part = unique_particles.size();
  auto const max_id = cell_structure.get_cached_max_local_particle_id();
  // phase 3.5: recompute the store-side derived director (same cadence as the
  // old commit-loop director write), then point the pack's store-aliased views
  // at the current store columns + translation/director views.
#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
  cell_structure.update_director_view();
#endif
  cell_structure.bind_pack_store_views();
  auto &aosoa = cell_structure.get_aosoa();

  if (rebuild) {
    auto &id_to_index = cell_structure.get_id_to_index();

    // ===================================================
    // Fill particle storage (full commit)
    // ===================================================
#ifdef ESPRESSO_CALIPER
    CALI_MARK_BEGIN("AoSoA commit full");
#endif
    int pair_count = 0;
    int angle_count = 0;
    int dihedral_count = 0;
    kokkos_parallel_range_for<policy_type>(
        "AoSoA write", std::size_t{0}, n_part,
        [&unique_particles, &aosoa, &id_to_index, &cell_structure, &pair_count,
         &angle_count, &dihedral_count](int const index) {
          auto const &p = *unique_particles.at(index);
          commit_particle(p, index, aosoa, true);
          id_to_index(p.id()) = index;
          if (not p.is_ghost()) {
            cell_structure.update_bond_storage(pair_count, angle_count,
                                               dihedral_count, p);
          }
        });
    Kokkos::fence();
    auto &bs = cell_structure.bond_state();
    auto &pair_bond_list = bs.pair_list;
    Kokkos::parallel_for("resolve_pair_bond_indices", pair_count,
                         [&pair_bond_list, &id_to_index](int idx) {
                           for (int col = 0; col < 2; ++col) {
                             pair_bond_list(idx, col) =
                                 id_to_index(pair_bond_list(idx, col));
                           }
                         });
    auto &angle_bond_list = bs.angle_list;
    Kokkos::parallel_for("resolve_angle_bond_indices", angle_count,
                         [&angle_bond_list, &id_to_index](int idx) {
                           for (int col = 0; col < 3; ++col) {
                             angle_bond_list(idx, col) =
                                 id_to_index(angle_bond_list(idx, col));
                           }
                         });
    auto &dihedral_bond_list = bs.dihedral_list;
    Kokkos::parallel_for("resolve_dihedral_bond_indices", dihedral_count,
                         [&dihedral_bond_list, &id_to_index](int idx) {
                           for (int col = 0; col < 4; ++col) {
                             dihedral_bond_list(idx, col) =
                                 id_to_index(dihedral_bond_list(idx, col));
                           }
                         });
    Kokkos::fence();
#ifdef ESPRESSO_CALIPER
    CALI_MARK_END("AoSoA commit full");
#endif

    // ===================================================
    // Get Verlet pairs and fill Verlet list
    // ===================================================
    bool rebuild_vl = (integ_switch != INTEG_METHOD_STEEPEST_DESCENT and
                       cell_structure.use_verlet_list);
#ifdef ESPRESSO_CALIPER
    CALI_MARK_BEGIN("Verlet list creation");
#endif
    cell_structure.rebuild_verlet_list_cabana(
        [&](std::span<Cell *const> cells, BoxGeometry const &box,
            CellStructure::ListType &verlet_list) {
          link_cell_kokkos(
              std::move(cells), box, verlet_criterion, id_to_index, max_id,
              [&](const int i, const int j) {
                // intra cell loop
                verlet_list.addNeighborLB(i, j);
              },
              [&](const int i, const int j) {
                // inter cell loop
                verlet_list.addNeighbor(i, j);
              });

          if (verlet_list.hasOverflow()) {
            cell_structure.use_verlet_list = false;
            runtimeWarningMsg()
                << "Verlet list overflow detected: neighbor count exceeded "
                   "max_counts. Falling back to the link cell algorithm. "
                   "Configured max is "
                << Cabana::NeighborList<CellStructure::ListType>::maxNeighbor(
                       verlet_list);
          }
        },
        rebuild_vl);
#ifdef ESPRESSO_CALIPER
    CALI_MARK_END("Verlet list creation");
#endif
  }
  // Partial-update (no-rebuild) branch: NOTHING to commit. Since phase 3.5 the
  // pack no longer owns per-step particle data (position/velocity/etc. alias
  // the ParticleStore columns, read by store row), and the only remaining
  // pack-owned commit writes -- `type` and the has-exclusion `flags` -- are
  // rebuild-cadence data written on the full-commit branch above (both are
  // refreshed by a forced rebuild whenever they can change). The old per-step
  // commit_particle(rebuild=false) loop therefore had no work left to do; it
  // iterated every unique particle each non-rebuild step doing nothing. Dropped
  // (was an O(N)-per-step no-op). The store-view rebind (bind_pack_store_views
  // above) already repointed the aliased views for this step.
}

#ifdef ESPRESSO_ELECTROSTATICS
// phase-5 perf recovery: refresh the PACK-OWNED contiguous charge column from
// the authoritative ParticleStore q column. This runs once per step (O(N))
// ONLY when a coulomb actor is active; the hot pair loop and the P3M
// gather/spread loops then read `aosoa.charge(pack_index)` contiguously, which
// is O(pairs) >> O(N). Pure-LJ runs never call this and pay zero cost.
ESPRESSO_ATTR_ALWAYS_INLINE inline void
refresh_pack_charges(CellStructure &cell_structure) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using policy_type = Kokkos::RangePolicy<execution_space>;
  auto const &unique_particles = cell_structure.get_unique_particles();
  auto const n_part = unique_particles.size();
  auto &aosoa = cell_structure.get_aosoa();
  kokkos_parallel_range_for<policy_type>(
      "refresh pack charges", std::size_t{0}, n_part,
      [&unique_particles, &aosoa](std::size_t const index) {
        aosoa.pair_charge(index) = unique_particles.at(index)->q();
      });
}
#endif

#ifdef ESPRESSO_DIPOLES
// phase-5 perf recovery: pack-owned dipm column refresh, guarded by an active
// dipolar actor (see refresh_pack_charges for the rationale).
ESPRESSO_ATTR_ALWAYS_INLINE inline void
refresh_pack_dipm(CellStructure &cell_structure) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using policy_type = Kokkos::RangePolicy<execution_space>;
  auto const &unique_particles = cell_structure.get_unique_particles();
  auto const n_part = unique_particles.size();
  auto &aosoa = cell_structure.get_aosoa();
  kokkos_parallel_range_for<policy_type>(
      "refresh pack dipm", std::size_t{0}, n_part,
      [&unique_particles, &aosoa](std::size_t const index) {
        aosoa.pair_dipm(index) = unique_particles.at(index)->dipm();
      });
}
#endif

void cabana_short_range(auto const &pair_bonds_kernel,
                        auto const &angle_bonds_kernel,
                        auto const &dihedral_bonds_kernel,
                        auto const &nonbonded_kernel,
                        CellStructure &cell_structure, double pair_cutoff,
                        double bond_cutoff, auto const &verlet_criterion,
                        auto const integ_switch) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  assert(cell_structure.get_resort_particles() == Cells::RESORT_NONE);

  if (bond_cutoff >= 0.) {
#ifdef ESPRESSO_CALIPER
    CALI_MARK_BEGIN("cabana_bond_loop");
#endif
    auto const n_pair_bonds = cell_structure.get_local_pair_bond_numbers();
    auto const n_angle_bonds = cell_structure.get_local_angle_bond_numbers();
    auto const n_dihedral_bonds =
        cell_structure.get_local_dihedral_bond_numbers();
    if (n_pair_bonds > 0) {
      Kokkos::parallel_for( // loop over bonds
          "for_each_local_pair_bonds", n_pair_bonds, pair_bonds_kernel);
      Kokkos::fence();
    }
    if (n_angle_bonds > 0) {
      Kokkos::parallel_for( // loop over bonds
          "for_each_local_angle_bonds", n_angle_bonds, angle_bonds_kernel);
      Kokkos::fence();
    }
    if (n_dihedral_bonds > 0) {
      Kokkos::parallel_for( // loop over bonds
          "for_each_local_dihedral_bonds", n_dihedral_bonds,
          dihedral_bonds_kernel);
      Kokkos::fence();
    }
#ifdef ESPRESSO_CALIPER
    CALI_MARK_END("cabana_bond_loop");
#endif
  }

  // Cabana short range loop
  if (pair_cutoff > 0.) {
#ifdef ESPRESSO_CALIPER
    CALI_MARK_BEGIN("cabana_pair_loop");
#endif
    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT and
        cell_structure.use_verlet_list) {
      auto const &verlet_list = cell_structure.get_verlet_list_cabana();
      Kokkos::RangePolicy<execution_space> policy(
          std::size_t{0}, cell_structure.get_unique_particles().size());
      Cabana::neighbor_parallel_for(policy, nonbonded_kernel, verlet_list,
                                    Cabana::FirstNeighborsTag(),
                                    Cabana::SerialOpTag());
    } else {
      cell_structure.cell_list_loop(
          [&](std::span<Cell *const> cells, BoxGeometry const &box) {
            link_cell_kokkos(
                std::move(cells), box, verlet_criterion,
                cell_structure.get_id_to_index(),
                cell_structure.get_cached_max_local_particle_id(),
                [&](const int i, const int j) {
                  // intra cell loop
                  nonbonded_kernel(i, j);
                },
                [&](const int i, const int j) {
                  // inter cell loop
                  nonbonded_kernel(i, j);
                });
          });
    }
    Kokkos::fence();
#ifdef ESPRESSO_CALIPER
    CALI_MARK_END("cabana_pair_loop");
#endif
  }
}
