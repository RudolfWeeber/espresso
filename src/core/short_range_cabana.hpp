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
#include "kokkos_helpers.hpp"

#include <Cabana_Core.hpp>
#include <Cabana_NeighborList.hpp>

#ifdef ESPRESSO_CALIPER
#include "caliper_utils.hpp"
#endif

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

ESPRESSO_ATTR_ALWAYS_INLINE inline void
commit_particle(Particle const &p, auto const index,
                CellStructure::AoSoA_pack &aosoa, bool const rebuild) {
  // Always commit: positions, velocities, charges, directors, dipm
  aosoa.set_vector_at(aosoa.position, index, p.pos());
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
  bool const has_exclusion = not p.exclusions().empty();
  aosoa.set_has_exclusion(index, has_exclusion);
  // Record the any-exclusion aggregate on the host. This is intentionally NOT
  // done inside the device-qualified set_has_exclusion: this commit sweep runs
  // on the host execution space, and the aggregate is a host std::atomic.
  if (has_exclusion) {
    aosoa.mark_any_exclusion();
  }
#else
  aosoa.flags(index) = 0;
#endif
}

/**
 * @brief Link-cell pair search for cuboid boxes.
 *
 * Finds the same pairs, in the same order, as the scalar loops of
 * @ref link_cell_kokkos, but splits the candidate test in two stages. Once
 * per call, the positions of all particles in the local cells and their
 * neighbor cells are copied into SoA arrays, cell by cell, together with
 * the particle's local index (-1 for particles that are skipped). For each
 * particle, the squared minimum-image distances to a block of candidates
 * are then computed by a branch-free loop over these arrays that the
 * compiler can vectorize. Only candidates within the global
 * (cutoff + skin) bound of the criterion, a small fraction of the
 * candidates in the 27-cell neighborhood, go on to the index filter and
 * the full criterion, which reads the particle data.
 */
inline void link_cell_cuboid(
    std::span<Cell *const> cells, CuboidMinimumImage const &minimum_image,
    auto const &verlet_criterion,
    Kokkos::View<int *, Kokkos::DefaultHostExecutionSpace> const &id_to_index,
    int const max_id, auto const &intra_operator, auto const &inter_operator) {
  using execution_space = Kokkos::DefaultHostExecutionSpace;

  // slots: the local cells first (slot i = cells[i]), then the neighbor
  // cells that are not local (ghost cells), in order of first appearance
  std::vector<Cell const *> slot_cells(cells.begin(), cells.end());
  std::unordered_map<Cell const *, std::size_t> slot_of;
  slot_of.reserve(2u * cells.size());
  for (std::size_t i = 0; i < cells.size(); ++i) {
    slot_of.emplace(cells[i], i);
  }
  // red neighbor slots of each local cell
  std::vector<std::size_t> neighbor_slots;
  std::vector<std::size_t> neighbor_begin(cells.size() + 1u, 0u);
  for (std::size_t i = 0; i < cells.size(); ++i) {
    for (auto const *neighbor : cells[i]->neighbors().red()) {
      auto const [it, inserted] =
          slot_of.try_emplace(neighbor, slot_cells.size());
      if (inserted) {
        slot_cells.emplace_back(neighbor);
      }
      neighbor_slots.emplace_back(it->second);
    }
    neighbor_begin[i + 1u] = neighbor_slots.size();
  }
  // particle range [offset[k], offset[k + 1]) of slot k
  auto const n_slots = slot_cells.size();
  std::vector<std::size_t> offset(n_slots + 1u, 0u);
  for (std::size_t k = 0; k < n_slots; ++k) {
    offset[k + 1u] = offset[k] + slot_cells[k]->particles().size();
  }
  auto const n_candidates = offset[n_slots];
  std::vector<double> pos_x(n_candidates);
  std::vector<double> pos_y(n_candidates);
  std::vector<double> pos_z(n_candidates);
  std::vector<int> index(n_candidates);
  std::vector<Particle const *> particle(n_candidates);
  kokkos_parallel_range_for<execution_space>(
      "verlet_pack", std::size_t{0}, n_slots, [&](std::size_t const k) {
        auto j = offset[k];
        for (auto const &p : slot_cells[k]->particles()) {
          auto const &pos = p.pos();
          pos_x[j] = pos[0];
          pos_y[j] = pos[1];
          pos_z[j] = pos[2];
          index[j] = (p.id() <= max_id) ? id_to_index(p.id()) : -1;
          particle[j] = &p;
          ++j;
        }
      });
  Kokkos::fence();

  // test particle @p a against the candidates [begin, end)
  auto const max_cut2 = verlet_criterion.eff_max_cut2();
  auto const scan = [&, max_cut2](std::size_t const a, std::size_t const begin,
                                  std::size_t const end, auto const &op) {
    constexpr std::size_t block_size = 64u;
    double dist2[block_size];
    for (auto j0 = begin; j0 < end; j0 += block_size) {
      auto const m = std::min(block_size, end - j0);
      minimum_image.batch_dist2(pos_x[a], pos_y[a], pos_z[a], m,
                                pos_x.data() + j0, pos_y.data() + j0,
                                pos_z.data() + j0, dist2);
      for (std::size_t t = 0; t < m; ++t) {
        auto const j = j0 + t;
        if (dist2[t] <= max_cut2 and index[j] >= 0 and
            verlet_criterion(*particle[a], *particle[j], dist2[t])) {
          op(index[a], index[j]);
        }
      }
    }
  };

  auto intra_kernel = [&](const int i) {
    auto const end = offset[i + 1];
    for (auto a = offset[i]; a < end; ++a) {
      if (index[a] >= 0) {
        // pairs in this cell
        scan(a, a + 1u, end, intra_operator);
      }
    }
  };

  auto inter_kernel = [&](const int i) {
    for (auto a = offset[i]; a < offset[i + 1]; ++a) {
      if (index[a] >= 0) {
        // pairs with neighboring cells
        for (auto n = neighbor_begin[i]; n < neighbor_begin[i + 1]; ++n) {
          auto const k = neighbor_slots[n];
          scan(a, offset[k], offset[k + 1u], inter_operator);
        }
      }
    }
  };

  kokkos_parallel_range_for<execution_space>("intra", std::size_t{0},
                                             cells.size(), intra_kernel);
  Kokkos::fence();

  kokkos_parallel_range_for<execution_space>("inter", std::size_t{0},
                                             cells.size(), inter_kernel);
  Kokkos::fence();
}

/**
 * @brief Whether @ref link_cell_cuboid is expected to be faster than the
 * scalar loops of @ref link_cell_kokkos.
 *
 * The prefilter has a serial setup cost per call that is proportional to
 * the number of cells (slot map, neighbor slots and particle ranges),
 * while its gain is proportional to the number of candidate pairs, which
 * are processed in parallel. The candidate pairs of a cell scale with the
 * square of its occupation, so the prefilter is used when the mean squared
 * occupation of the local cells exceeds a threshold that grows linearly
 * with the number of threads. The coefficients are the measured break-even
 * points (cell-list force loop, Lennard-Jones and WCA liquids and gases)
 * on AMD EPYC 9374F with 1, 4, 8 and 16 threads, rounded up: below them,
 * in dilute systems, the setup dominates and the scalar loops are faster.
 */
inline bool link_cell_prefilter_pays_off(std::span<Cell *const> cells) {
  if (cells.empty()) {
    return false;
  }
  auto const n_threads = Kokkos::num_threads();
  auto const threshold = 8. + 4.5 * static_cast<double>(n_threads);
  auto const occupation2 = [cells](std::size_t const i) {
    auto const n = cells[i]->particles().size();
    return n * n;
  };
  std::size_t sum_n2 = 0u;
  if (n_threads > 1) {
    // one cache miss per cell: parallel, so the test stays negligible
    // in dilute systems with many cells
    Kokkos::parallel_reduce(
        "verlet_prefilter_occupation",
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(std::size_t{0},
                                                               cells.size()),
        [&occupation2](std::size_t const i, std::size_t &acc) {
          acc += occupation2(i);
        },
        sum_n2);
  } else {
    for (std::size_t i = 0u; i < cells.size(); ++i) {
      sum_n2 += occupation2(i);
    }
  }
  return static_cast<double>(sum_n2) >=
         threshold * static_cast<double>(cells.size());
}

ESPRESSO_ATTR_ALWAYS_INLINE inline void link_cell_kokkos(
    std::span<Cell *const> cells, BoxGeometry const &box_geo,
    auto const &verlet_criterion,
    Kokkos::View<int *, Kokkos::DefaultHostExecutionSpace> const &id_to_index,
    int const max_id, auto const &intra_operator, auto const &inter_operator) {

  // implementation detail: max_id refers to the max local particle id,
  // but ghost particles from other ranks may have larger particle ids;
  // -1 is used as a sentinel value for particle ids from other threads

  if (box_geo.type() != BoxType::LEES_EDWARDS and
      link_cell_prefilter_pays_off(cells)) {
    link_cell_cuboid(cells, box_geo.cuboid_minimum_image(), verlet_criterion,
                     id_to_index, max_id, intra_operator, inter_operator);
    return;
  }

  // Hoist the cuboid minimum-image parameters by value: the fold runs once
  // per candidate pair and must not chase the BoxGeometry reference for its
  // box lengths every time. Lees-Edwards boxes take the full BoxGeometry
  // path (shear offset handling).
  bool const has_lees_edwards = box_geo.type() == BoxType::LEES_EDWARDS;
  auto const cuboid_minimum_image = box_geo.cuboid_minimum_image();
  auto const minimum_image_dist2 =
      [&box_geo, has_lees_edwards, cuboid_minimum_image](
          Utils::Vector3d const &a, Utils::Vector3d const &b) {
        return has_lees_edwards ? box_geo.get_mi_dist2(a, b)
                                : cuboid_minimum_image.dist2(a, b);
      };

  auto intra_kernel = [&cells, minimum_image_dist2, &verlet_criterion,
                       &id_to_index, &intra_operator, max_id](const int i) {
    auto &local_particles = cells[i]->particles();
    for (auto it = local_particles.begin(); it != local_particles.end(); ++it) {
      auto const &p1 = *it;
      if (p1.id() <= max_id) {
        auto const ii = id_to_index(p1.id());
        if (ii >= 0) {
          // pairs in this cell
          for (auto jt = std::next(it); jt != local_particles.end(); ++jt) {
            if ((*jt).id() <= max_id) {
              if (verlet_criterion(p1, *jt,
                                   minimum_image_dist2(p1.pos(), jt->pos()))) {
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

  auto inter_kernel = [&cells, minimum_image_dist2, &verlet_criterion,
                       &id_to_index, &inter_operator, max_id](const int i) {
    auto &local_particles = cells[i]->particles();
    for (auto const &p1 : local_particles) {
      if (p1.id() <= max_id) {
        auto const ii = id_to_index(p1.id());
        if (ii >= 0) {
          // pairs with neighboring cells
          for (auto &neighbor : cells[i]->neighbors().red()) {
            for (auto const &p2 : neighbor->particles()) {
              if (p2.id() <= max_id) {
                if (verlet_criterion(p1, p2,
                                     minimum_image_dist2(p1.pos(), p2.pos()))) {
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

  kokkos_parallel_range_for<Kokkos::DefaultHostExecutionSpace>(
      "intra", std::size_t{0}, cells.size(), intra_kernel);
  Kokkos::fence();

  kokkos_parallel_range_for<Kokkos::DefaultHostExecutionSpace>(
      "inter", std::size_t{0}, cells.size(), inter_kernel);
  Kokkos::fence();
}

// @p make_verlet_criterion is a nullary factory: constructing the criterion
// fills an O(n_types^2) cutoff table, so it is only invoked when the Verlet
// list is actually rebuilt, not on every force call.
template <class execution_space = Kokkos::DefaultHostExecutionSpace>
ESPRESSO_ATTR_ALWAYS_INLINE inline void
update_cabana_state(CellStructure &cell_structure,
                    auto const &make_verlet_criterion, double const pair_cutoff,
                    auto const integ_switch) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
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
    ESPRESSO_CALI_MARK_BEGIN("AoSoA commit full");
#endif
    int pair_count = 0;
    int angle_count = 0;
    int dihedral_count = 0;
#ifdef ESPRESSO_EXCLUSIONS
    // commit_particle accumulates the any-exclusion aggregate below; clear it
    // before the sweep repopulates it (read O(1) at the dispatch gate).
    aosoa.reset_any_exclusion();
#endif
    kokkos_parallel_range_for<execution_space>(
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
    using host_space = Kokkos::DefaultHostExecutionSpace;
    auto &bs = cell_structure.bond_state();
    if (pair_count) {
      auto &pair_bond_list = bs.pair_list;
      kokkos_parallel_range_for<host_space>(
          "resolve_pair_bond_indices", std::size_t{0}, pair_count,
          [&pair_bond_list, &id_to_index](int idx) {
            for (int col = 0; col < 2; ++col) {
              pair_bond_list(idx, col) = id_to_index(pair_bond_list(idx, col));
            }
          });
    }
    if (angle_count) {
      auto &angle_bond_list = bs.angle_list;
      kokkos_parallel_range_for<host_space>(
          "resolve_angle_bond_indices", std::size_t{0}, angle_count,
          [&angle_bond_list, &id_to_index](int idx) {
            for (int col = 0; col < 3; ++col) {
              angle_bond_list(idx, col) =
                  id_to_index(angle_bond_list(idx, col));
            }
          });
    }
    if (dihedral_count) {
      auto &dihedral_bond_list = bs.dihedral_list;
      kokkos_parallel_range_for<host_space>(
          "resolve_dihedral_bond_indices", std::size_t{0}, dihedral_count,
          [&dihedral_bond_list, &id_to_index](int idx) {
            for (int col = 0; col < 4; ++col) {
              dihedral_bond_list(idx, col) =
                  id_to_index(dihedral_bond_list(idx, col));
            }
          });
    }
    if (pair_count != 0 or angle_count != 0 or dihedral_count != 0) {
      Kokkos::fence();
    }
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_END("AoSoA commit full");
#endif

    // ===================================================
    // Get Verlet pairs and fill Verlet list
    // ===================================================
    bool rebuild_vl = (integ_switch != INTEG_METHOD_STEEPEST_DESCENT and
                       cell_structure.use_verlet_list);
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_BEGIN("Verlet list creation");
#endif
    cell_structure.rebuild_verlet_list_cabana(
        [&](std::span<Cell *const> cells, BoxGeometry const &box,
            CellStructure::ListType &verlet_list) {
          auto const verlet_criterion = make_verlet_criterion();
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
    ESPRESSO_CALI_MARK_END("Verlet list creation");
#endif
  } else {
    // ===================================================
    // Fill particle storage (partial update)
    // ===================================================
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_BEGIN("AoSoA commit partial");
#endif
#ifdef ESPRESSO_EXCLUSIONS
    // commit_particle accumulates the any-exclusion aggregate below; clear it
    // before the sweep repopulates it (read O(1) at the dispatch gate).
    aosoa.reset_any_exclusion();
#endif
    kokkos_parallel_range_for<execution_space>(
        "AoSoA write", std::size_t{0}, n_part,
        [&unique_particles, &aosoa](int const index) {
          auto const &p = *unique_particles.at(index);
          commit_particle(p, index, aosoa, false);
        });
    Kokkos::fence();
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_END("AoSoA commit partial");
#endif
  }
}

#ifdef ESPRESSO_ELECTROSTATICS
template <class execution_space = Kokkos::DefaultHostExecutionSpace>
ESPRESSO_ATTR_ALWAYS_INLINE inline void
update_aosoa_charges(CellStructure &cell_structure) {
  auto const &unique_particles = cell_structure.get_unique_particles();
  auto const n_part = unique_particles.size();
  auto &aosoa = cell_structure.get_aosoa();

  kokkos_parallel_range_for<execution_space>(
      "Views update charges", std::size_t{0}, n_part,
      [&unique_particles, &aosoa](std::size_t const index) {
        aosoa.charge(index) = unique_particles.at(index)->q();
      });
}
#endif

// Optional replacement for the Verlet-list pair loop. When set, it is invoked
// with the current Verlet list and the pack particle count instead of the
// generic Cabana::neighbor_parallel_for over @p nonbonded_kernel. forces.cpp
// installs the compile-time-specialized pair kernel this way when the active
// feature set allows it; energy/pressure leave it empty and keep the generic
// loop.
using ShortRangeVerletPairLoop =
    std::function<void(CellStructure::ListType const &, std::size_t)>;

// @p make_verlet_criterion is a nullary factory: constructing the criterion
// fills an O(n_types^2) cutoff table, so it is only invoked on the link-cell
// fallback path, which is the only consumer here.
template <class execution_space = Kokkos::DefaultHostExecutionSpace>
void cabana_short_range(auto const &pair_bonds_kernel,
                        auto const &angle_bonds_kernel,
                        auto const &dihedral_bonds_kernel,
                        auto const &nonbonded_kernel,
                        CellStructure &cell_structure, double pair_cutoff,
                        double bond_cutoff, auto const &make_verlet_criterion,
                        auto const integ_switch,
                        ShortRangeVerletPairLoop const &verlet_pair_loop = {}) {
  assert(cell_structure.get_resort_particles() == Cells::RESORT_NONE);

  if (bond_cutoff >= 0.) {
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_BEGIN("cabana_bond_loop");
#endif
    using host_space = Kokkos::DefaultHostExecutionSpace;
    auto const n_pair_bonds = cell_structure.get_local_pair_bond_numbers();
    auto const n_angle_bonds = cell_structure.get_local_angle_bond_numbers();
    auto const n_dihedral_bonds =
        cell_structure.get_local_dihedral_bond_numbers();
    if (n_pair_bonds > 0) {
      kokkos_parallel_range_for<host_space>("for_each_local_pair_bonds",
                                            std::size_t{0}, n_pair_bonds,
                                            pair_bonds_kernel);
    }
    if (n_angle_bonds > 0) {
      kokkos_parallel_range_for<host_space>("for_each_local_angle_bonds",
                                            std::size_t{0}, n_angle_bonds,
                                            angle_bonds_kernel);
    }
    if (n_dihedral_bonds > 0) {
      kokkos_parallel_range_for<host_space>("for_each_local_dihedral_bonds",
                                            std::size_t{0}, n_dihedral_bonds,
                                            dihedral_bonds_kernel);
    }
    if (n_pair_bonds != 0 or n_angle_bonds != 0 or n_dihedral_bonds != 0) {
      Kokkos::fence();
    }
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_END("cabana_bond_loop");
#endif
  }

  // Cabana short range loop
  if (pair_cutoff > 0.) {
#ifdef ESPRESSO_CALIPER
    ESPRESSO_CALI_MARK_BEGIN("cabana_pair_loop");
#endif
    if (integ_switch != INTEG_METHOD_STEEPEST_DESCENT and
        cell_structure.use_verlet_list) {
      auto const &verlet_list = cell_structure.get_verlet_list_cabana();
      auto const n_particles = cell_structure.get_unique_particles().size();
      if (verlet_pair_loop) {
        verlet_pair_loop(verlet_list, n_particles);
      } else if (Kokkos::num_threads() == 1) {
        using NL = Cabana::NeighborList<CellStructure::ListType>;
        for (std::size_t i = 0; i < n_particles; ++i) {
          auto const nn = NL::numNeighbor(verlet_list, i);
          for (std::size_t n = 0; n < nn; ++n) {
            nonbonded_kernel(i, NL::getNeighbor(verlet_list, i, n));
          }
        }
      } else {
        Kokkos::RangePolicy<execution_space> policy(std::size_t{0},
                                                    n_particles);
        Cabana::neighbor_parallel_for(policy, nonbonded_kernel, verlet_list,
                                      Cabana::FirstNeighborsTag(),
                                      Cabana::SerialOpTag());
      }
    } else {
      cell_structure.cell_list_loop(
          [&](std::span<Cell *const> cells, BoxGeometry const &box) {
            auto const verlet_criterion = make_verlet_criterion();
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
    ESPRESSO_CALI_MARK_END("cabana_pair_loop");
#endif
  }
}
