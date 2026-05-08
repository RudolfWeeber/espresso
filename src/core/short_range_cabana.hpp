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
  // Always commit: positions, velocities, charges, directors, dipm
  auto const pos = p.pos();
  aosoa.set_vector_at(aosoa.position, index, pos);
  aosoa.position_x(index) = pos[0];
  aosoa.position_y(index) = pos[1];
  aosoa.position_z(index) = pos[2];
#ifdef ESPRESSO_ELECTROSTATICS
  aosoa.charge(index) = p.q();
#endif
  aosoa.set_vector_at(aosoa.velocity, index, p.v());
#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
  aosoa.set_vector_at(aosoa.director, index,
                      Utils::convert_quaternion_to_director(p.quat()));
#endif
#ifdef ESPRESSO_DIPOLES
  aosoa.dipm(index) = p.dipm();
#endif

  // Only commit on rebuild: id, type
  if (rebuild) {
    aosoa.id(index) = p.id();
    aosoa.type(index) = p.type();
    aosoa.set_vector_at(aosoa.image, index, p.image_box());
#ifdef ESPRESSO_MASS
    aosoa.mass(index) = p.mass();
#endif
  }

  // Always update exclusion flags (they can change during simulation)
#ifdef ESPRESSO_EXCLUSIONS
  aosoa.set_has_exclusion(index, !p.exclusions().empty());
#else
  aosoa.flags(index) = 0;
#endif
}

namespace detail {

/** AoSoA index of a ghost-cell particle via id lookup.
 *  Reads p.id() from the Particle struct only.
 *  Returns -1 if the particle has no AoSoA entry. */
inline int ghost_particle_index(Particle const &p,
                                Kokkos::View<int *> const &id_to_index,
                                int max_id) {
  if (p.id() > max_id)
    return -1;
  return id_to_index(p.id());
}

} // namespace detail

ESPRESSO_ATTR_ALWAYS_INLINE inline void
link_cell_kokkos(std::span<Cell *const> cells, BoxGeometry const &box_geo,
                 auto const &verlet_criterion,
                 Kokkos::View<int *> const &id_to_index, int const max_id,
                 CellStructure::AoSoA_pack const &aosoa,
                 auto const &intra_operator, auto const &inter_operator) {

  // Hoist box constants for SIMD path (CUBOID + fully periodic only).
  bool const simd_eligible = box_geo.type() == BoxType::CUBOID and
                             box_geo.periodic(0u) and box_geo.periodic(1u) and
                             box_geo.periodic(2u);
  auto const lx = box_geo.length()[0];
  auto const ly = box_geo.length()[1];
  auto const lz = box_geo.length()[2];
  auto const lxi = box_geo.length_inv()[0];
  auto const lyi = box_geo.length_inv()[1];
  auto const lzi = box_geo.length_inv()[2];

  auto kernel = [&cells, box_geo, verlet_criterion, &aosoa, &id_to_index,
                 &intra_operator, &inter_operator, max_id, simd_eligible, lx,
                 ly, lz, lxi, lyi, lzi](const int i) {
    auto const base = cells[i]->aosoa_offset();
    auto const n_i = cells[i]->particles().size();

    // Intra-cell pairs.
    for (std::size_t k = 0; k < n_i; ++k) {
      auto const ii = base + k;
      auto const pos1 = aosoa.get_position(ii);
      for (std::size_t l = k + 1; l < n_i; ++l) {
        auto const jj = base + l;
        auto const dist2 = box_geo.get_mi_dist2(pos1, aosoa.get_position(jj));
        if (verlet_criterion(aosoa, ii, jj, dist2))
          intra_operator(static_cast<int>(ii), static_cast<int>(jj));
      }
    }

    // Inter-cell pairs (red neighbors only — visits each pair once).
    constexpr std::size_t MAX_NJ = 64;
    alignas(64) double dist2_buf[MAX_NJ];
    for (auto *neighbor : cells[i]->neighbors().red()) {
      auto const base_j = neighbor->aosoa_offset();
      if (base_j != Cell::no_aosoa_slot) {
        // Local neighbor — no Particle struct access.
        auto const n_j = neighbor->particles().size();
        bool const use_simd = simd_eligible and n_j <= MAX_NJ;
        for (std::size_t k = 0; k < n_i; ++k) {
          auto const ii = base + k;
          auto const pos1 = aosoa.get_position(ii);
          if (use_simd) {
            auto const px = pos1[0];
            auto const py = pos1[1];
            auto const pz = pos1[2];
#pragma omp simd
            for (std::size_t l = 0; l < n_j; ++l) {
              auto const dx0 = px - aosoa.position_x(base_j + l);
              auto const dy0 = py - aosoa.position_y(base_j + l);
              auto const dz0 = pz - aosoa.position_z(base_j + l);
              auto const dx = dx0 - std::rint(dx0 * lxi) * lx;
              auto const dy = dy0 - std::rint(dy0 * lyi) * ly;
              auto const dz = dz0 - std::rint(dz0 * lzi) * lz;
              dist2_buf[l] = dx * dx + dy * dy + dz * dz;
            }
            for (std::size_t l = 0; l < n_j; ++l) {
              auto const jj = base_j + l;
              if (verlet_criterion(aosoa, ii, jj, dist2_buf[l]))
                inter_operator(static_cast<int>(ii), static_cast<int>(jj));
            }
          } else {
            for (std::size_t l = 0; l < n_j; ++l) {
              auto const jj = base_j + l;
              auto const dist2 =
                  box_geo.get_mi_dist2(pos1, aosoa.get_position(jj));
              if (verlet_criterion(aosoa, ii, jj, dist2))
                inter_operator(static_cast<int>(ii), static_cast<int>(jj));
            }
          }
        }
      } else {
        // Ghost neighbor — read p2.id() only, rest from AoSoA.
        for (std::size_t k = 0; k < n_i; ++k) {
          auto const ii = base + k;
          auto const pos1 = aosoa.get_position(ii);
          for (auto const &p2 : neighbor->particles()) {
            auto const jj =
                detail::ghost_particle_index(p2, id_to_index, max_id);
            if (jj < 0)
              continue;
            auto const dist2 =
                box_geo.get_mi_dist2(pos1, aosoa.get_position(jj));
            if (verlet_criterion(aosoa, ii, jj, dist2))
              inter_operator(static_cast<int>(ii), static_cast<int>(jj));
          }
        }
      }
    }
  };

  Kokkos::RangePolicy<Kokkos::Schedule<Kokkos::Dynamic>> policy(0,
                                                                cells.size());
  policy.set_chunk_size(4);
  Kokkos::parallel_for("link_cell", policy, kernel);
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
              cell_structure.get_aosoa(),
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
  } else {
    // ===================================================
    // Fill particle storage (partial update)
    // ===================================================
#ifdef ESPRESSO_CALIPER
    CALI_MARK_BEGIN("AoSoA commit partial");
#endif
    kokkos_parallel_range_for<policy_type>(
        "AoSoA write", std::size_t{0}, n_part,
        [&unique_particles, &aosoa](int const index) {
          auto const &p = *unique_particles.at(index);
          commit_particle(p, index, aosoa, false);
        });
    Kokkos::fence();
#ifdef ESPRESSO_CALIPER
    CALI_MARK_END("AoSoA commit partial");
#endif
  }
}

#ifdef ESPRESSO_ELECTROSTATICS
ESPRESSO_ATTR_ALWAYS_INLINE inline void
update_aosoa_charges(CellStructure &cell_structure) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using policy_type = Kokkos::RangePolicy<execution_space>;
  auto const &unique_particles = cell_structure.get_unique_particles();
  auto const n_part = unique_particles.size();
  auto &aosoa = cell_structure.get_aosoa();

  kokkos_parallel_range_for<policy_type>(
      "Views update charges", std::size_t{0}, n_part,
      [&unique_particles, &aosoa](std::size_t const index) {
        aosoa.charge(index) = unique_particles.at(index)->q();
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
      auto const &counts_view = verlet_list.counts;
      auto const &neighbors_view = verlet_list.neighbors;
      auto const &aosoa = cell_structure.get_aosoa();
      auto const &box_geo = nonbonded_kernel.box_geo;
      double const cutoff_sq = nonbonded_kernel.system_max_cutoff_sq;
      bool const simd_eligible = box_geo.type() == BoxType::CUBOID and
                                 box_geo.periodic(0u) and
                                 box_geo.periodic(1u) and box_geo.periodic(2u);
      auto const lx = box_geo.length()[0];
      auto const ly = box_geo.length()[1];
      auto const lz = box_geo.length()[2];
      auto const lxi = box_geo.length_inv()[0];
      auto const lyi = box_geo.length_inv()[1];
      auto const lzi = box_geo.length_inv()[2];
      Kokkos::RangePolicy<execution_space> policy(
          std::size_t{0}, cell_structure.get_unique_particles().size());
      Kokkos::parallel_for(
          "cabana_pair_loop_simd", policy,
          [&nonbonded_kernel, &aosoa, &counts_view, &neighbors_view,
           simd_eligible, lx, ly, lz, lxi, lyi, lzi,
           cutoff_sq](std::size_t const i) {
            int const n = counts_view(i);
            if (n == 0)
              return;

            constexpr int MAX_NJ = 256;
            if (simd_eligible and n <= MAX_NJ) {
              alignas(64) double dx_buf[MAX_NJ];
              alignas(64) double dy_buf[MAX_NJ];
              alignas(64) double dz_buf[MAX_NJ];
              alignas(64) double dist2_buf[MAX_NJ];
              auto const px = aosoa.position_x(i);
              auto const py = aosoa.position_y(i);
              auto const pz = aosoa.position_z(i);
#pragma omp simd
              for (int k = 0; k < n; ++k) {
                int const j = neighbors_view(i, k);
                auto const dx0 = px - aosoa.position_x(j);
                auto const dy0 = py - aosoa.position_y(j);
                auto const dz0 = pz - aosoa.position_z(j);
                auto const dx = dx0 - std::rint(dx0 * lxi) * lx;
                auto const dy = dy0 - std::rint(dy0 * lyi) * ly;
                auto const dz = dz0 - std::rint(dz0 * lzi) * lz;
                dx_buf[k] = dx;
                dy_buf[k] = dy;
                dz_buf[k] = dz;
                dist2_buf[k] = dx * dx + dy * dy + dz * dz;
              }
              for (int k = 0; k < n; ++k) {
                if (dist2_buf[k] <= cutoff_sq) {
                  int const j = neighbors_view(i, k);
                  Utils::Vector3d const d{dx_buf[k], dy_buf[k], dz_buf[k]};
                  nonbonded_kernel(i, static_cast<std::size_t>(j), dist2_buf[k],
                                   d);
                }
              }
            } else {
              for (int k = 0; k < n; ++k) {
                int const j = neighbors_view(i, k);
                nonbonded_kernel(i, static_cast<std::size_t>(j));
              }
            }
          });
    } else {
      cell_structure.cell_list_loop(
          [&](std::span<Cell *const> cells, BoxGeometry const &box) {
            link_cell_kokkos(
                std::move(cells), box, verlet_criterion,
                cell_structure.get_id_to_index(),
                cell_structure.get_cached_max_local_particle_id(),
                cell_structure.get_aosoa(),
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
