/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "cell_system/CellStructure.hpp"

#include "cell_system/AtomDecomposition.hpp"
#include "cell_system/HybridDecomposition.hpp"
#include "cell_system/ParticleDecomposition.hpp"
#include "cell_system/ParticleListOperations.hpp"
#include "cell_system/RegularDecomposition.hpp"

#include "BoxGeometry.hpp"
#include "LocalBondState.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "aosoa_pack.hpp"
#include "cell_system/CellStructureType.hpp"
#include "communication.hpp"
#include "ghosts.hpp"
#include "integrators/Propagation.hpp"
#include "lees_edwards/lees_edwards.hpp"
#include "particle_enumeration.hpp"
#include "particle_reduction.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/contains.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/math/sqr.hpp>
#include <utils/quaternion.hpp>

#ifdef ESPRESSO_CALIPER
#include <caliper/cali.h>
#endif

#include <boost/mpi/collectives/all_reduce.hpp>

#include <omp.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <numbers>
#include <optional>
#include <ranges>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

CellStructure::~CellStructure() {
  clear_local_properties();
  // Release the ParticleStore columns while Kokkos is still alive; the handle
  // reset below may finalize Kokkos, and leaving the store's Views to
  // member-destruction would free them after Kokkos::finalize.
  m_particle_store.release_columns();
  // Same Kokkos-lifetime constraint for the (phase 7b) migration staging store.
  m_staging_store.release_columns();
  // Kokkos handle can only be freed after all Cabana containers have been freed
  m_kokkos_handle.reset();
}

void CellStructure::clear_local_properties() {
  m_scatter_force.reset();
  m_local_force.reset();
#ifdef ESPRESSO_ROTATION
  m_scatter_torque.reset();
  m_local_torque.reset();
#endif
#ifdef ESPRESSO_NPT
  m_scatter_virial.reset();
  m_local_virial.reset();
#endif
  m_id_to_index.reset();
  m_aosoa.reset();
  m_verlet_list_cabana.reset();
  m_bond_state->clear();
  // Release the phase-3.5 Kokkos-backed views owned here (translation view and
  // derived director). The pack's aliasing views were just dropped by
  // m_aosoa.reset() above, so releasing the owners now is safe and must happen
  // while Kokkos is still alive (see the destructor's release ordering).
  m_pack_index_to_store_row = Kokkos::View<int *, Kokkos::HostSpace>{};
#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
  m_director_view = Kokkos::View<double *[3], ParticleStore::StateVectorLayout,
                                 Kokkos::HostSpace>{};
#endif
  m_rebuild_verlet_list_cabana = true;
}
void CellStructure::clear_bond_properties() { m_bond_state->reset(); }

void CellStructure::set_kokkos_handle(std::shared_ptr<KokkosHandle> handle) {
  m_kokkos_handle = std::move(handle);
  m_bond_state = std::make_unique<LocalBondState>();
}

static auto estimate_max_counts(double pair_cutoff,
                                std::size_t number_of_unique_particles,
                                double local_box_volume,
                                std::size_t num_local_particles) {
  if (std::isinf(pair_cutoff)) {
    return number_of_unique_particles;
  }
  if (pair_cutoff < 0.) {
    pair_cutoff = 0.;
  }
  // Estimate number of neighbors based on local density and cutoff sphere:
  // volume n_neighbors = rho * (4/3) * pi * r^3, where rho = n_particles /
  // volume
  auto const local_density =
      (local_box_volume > 0. && num_local_particles > 0)
          ? static_cast<double>(num_local_particles) / local_box_volume
          : 0.;
  auto const cutoff_sphere_volume =
      (4. / 3.) * std::numbers::pi * Utils::int_pow<3>(pair_cutoff);
  // account for local fluctuations. Empirical.
  auto const fluctuation_factor = 2.;
  auto max_counts = static_cast<std::size_t>(
      std::ceil(fluctuation_factor * local_density * cutoff_sphere_volume));
  std::size_t constexpr threshold_num = 16;
  if (max_counts < threshold_num) {
    max_counts = std::min(threshold_num, number_of_unique_particles);
  }
  return max_counts;
}

void CellStructure::rebuild_local_properties(double const pair_cutoff) {
#ifdef ESPRESSO_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  assert(m_kokkos_handle);
  auto const num_part = get_unique_particles().size();
  auto const &system = get_system();
  auto const local_box_volume = system.local_geo->volume();
  auto max_counts = estimate_max_counts(pair_cutoff, num_part, local_box_volume,
                                        get_num_local_particles_cached());
#ifdef ESPRESSO_COLLISION_DETECTION
  if (system.has_collision_detection_enabled()) {
    // TODO: use other types of Verlet list data structures
    max_counts = num_part * 2ul;
  }
#endif
  if (m_local_force) { // local properties are reallocated
    Kokkos::realloc(get_local_force(), num_part);
    // underlying View extent changed -> scratch buffers must be rebuilt
    m_scatter_force.emplace(
        Kokkos::Experimental::create_scatter_view(get_local_force()));
#ifdef ESPRESSO_ROTATION
    Kokkos::realloc(get_local_torque(), num_part);
    // underlying View extent changed -> scratch buffers must be rebuilt
    m_scatter_torque.emplace(
        Kokkos::Experimental::create_scatter_view(get_local_torque()));
#endif
    Kokkos::realloc(get_id_to_index(), get_cached_max_local_particle_id() + 1);
    Kokkos::deep_copy(get_id_to_index(), -1);
    // Resize particle views using AoSoA_pack's resize method
    m_aosoa->resize(num_part);
    Kokkos::deep_copy(m_aosoa->flags, uint8_t{0});
    m_verlet_list_cabana->reallocData(num_part, max_counts);
  } else { // local properties are initialized
    m_local_force = std::make_unique<ForceType>("local_force", num_part);
    m_scatter_force.emplace(
        Kokkos::Experimental::create_scatter_view(*m_local_force));
#ifdef ESPRESSO_ROTATION
    m_local_torque = std::make_unique<ForceType>("local_torque", num_part);
    m_scatter_torque.emplace(
        Kokkos::Experimental::create_scatter_view(*m_local_torque));
#endif
    m_id_to_index = std::make_unique<Kokkos::View<int *>>(
        Kokkos::ViewAllocateWithoutInitializing("id_to_index"),
        get_cached_max_local_particle_id() + 1);
    Kokkos::deep_copy(get_id_to_index(), -1);
    // Create AoSoA_pack and initialize with resize
    m_aosoa = std::make_unique<AoSoA_pack>();
    m_aosoa->resize(num_part);
    Kokkos::deep_copy(m_aosoa->flags, uint8_t{0});

    m_verlet_list_cabana =
        std::make_unique<ListType>(0ul, num_part, max_counts);
  }
#ifdef ESPRESSO_NPT
  m_local_virial = std::make_unique<VirialType>("local_virial");
  m_scatter_virial.emplace(
      Kokkos::Experimental::create_scatter_view(*m_local_virial));
#endif
}

void CellStructure::reset_local_force_and_torque() {
#ifdef ESPRESSO_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  Kokkos::deep_copy(get_local_force(), 0.);
  m_scatter_force->reset();
#ifdef ESPRESSO_ROTATION
  Kokkos::deep_copy(get_local_torque(), 0.);
  m_scatter_torque->reset();
#endif
}

void CellStructure::reset_local_properties() {
  Kokkos::deep_copy(get_local_force(), 0.);
  m_scatter_force->reset();
#ifdef ESPRESSO_ROTATION
  Kokkos::deep_copy(get_local_torque(), 0.);
  m_scatter_torque->reset();
#endif
#ifdef ESPRESSO_NPT
  Kokkos::deep_copy(get_local_virial(), 0.);
  m_scatter_virial->reset();
#endif
  // phase 6: the has-exclusion `flags` column is NOT zeroed here. This runs on
  // every partial-update (no-rebuild) step; the flag is rebuild-cadence data
  // (commit_particle writes it only on rebuild, and exclusion changes force a
  // rebuild), so it is valid across partial steps and must be preserved. It is
  // (re)zeroed and rewritten on the rebuild path (rebuild_local_properties +
  // the full-commit loop). Zeroing it here would drop live exclusion flags for
  // the rest of the step.
}

void CellStructure::update_bond_storage(int &pair_count, int &angle_count,
                                        int &dihedral_count,
                                        Particle const &p) {
  auto &pair_list = m_bond_state->pair_list;
  auto &pair_ids = m_bond_state->pair_ids;
  auto &angle_list = m_bond_state->angle_list;
  auto &angle_ids = m_bond_state->angle_ids;
  auto &dihedral_list = m_bond_state->dihedral_list;
  auto &dihedral_ids = m_bond_state->dihedral_ids;
  for (auto const bond : p.bonds()) {
    auto const partner_ids = bond.partner_ids();
    try {
      auto const partners = resolve_bond_partners(partner_ids);
      if (partners.size() == 1u) { // pair bonds
        auto p_index = Kokkos::atomic_fetch_add(&pair_count, 1);
        pair_list(p_index, 0) = p.id();
        pair_list(p_index, 1) = partners[0]->id();
        pair_ids(p_index) = bond.bond_id();
      } else if (partners.size() == 2u) { // angle bond
        auto a_index = Kokkos::atomic_fetch_add(&angle_count, 1);
        angle_list(a_index, 0) = p.id();
        angle_list(a_index, 1) = partners[0]->id();
        angle_list(a_index, 2) = partners[1]->id();
        angle_ids(a_index) = bond.bond_id();
      } else if (partners.size() == 3u) { // dihedral bond
        auto d_index = Kokkos::atomic_fetch_add(&dihedral_count, 1);
        dihedral_list(d_index, 0) = p.id();
        dihedral_list(d_index, 1) = partners[0]->id();
        dihedral_list(d_index, 2) = partners[1]->id();
        dihedral_list(d_index, 3) = partners[2]->id();
        dihedral_ids(d_index) = bond.bond_id();
      }
    } catch (BondResolutionError const &) {
      bond_resolution_error(partner_ids);
    }
  }
}

void CellStructure::set_index_map() {
#ifdef ESPRESSO_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  // Phase 7a: cells hand out row VIEWS (transient cached objects), so
  // m_unique_particles can no longer hold pointers into cell iteration. Instead
  // it points into the stable view pool (m_view_pool, rebuilt by
  // rebuild_particle_index): its local prefix is store rows [0, n_local) in row
  // order (== the pack order, so pack index i == store row i on the local
  // prefix) and its tail is the deduped ghosts (first occurrence of each id not
  // owned by a local) -- exactly the point set and order this function built
  // before. Rebuilding those pointers here is a wholesale refresh over the pool
  // the store rebuild just produced.
  assert(not m_particle_store.is_dirty());
  auto &unique_particles = m_unique_particles;
  unique_particles.clear();
  auto const n_local = count_local_particles();
  unique_particles.reserve(m_view_pool_fill);
  int max_id = -1;
  int pair_count = 0;
  int angle_count = 0;
  int dihedral_count = 0;
  m_bond_state->reset_counts();
  // Iterate the LIVE view-pool prefix [0, m_view_pool_fill): local prefix (pool
  // positions [0, n_local)) then the deduped ghost tail. Slots past
  // m_view_pool_fill are idle reusable capacity (phase 7a perf fix) and must be
  // skipped. Count bonds only on the local prefix (matching the pre-flip
  // enumerate_local_particles pass, which never walked ghosts).
  for (std::size_t pool_index = 0u; pool_index < m_view_pool_fill;
       ++pool_index) {
    auto &p = m_view_pool[pool_index];
    unique_particles.emplace_back(std::addressof(p));
    max_id = std::max(p.id(), max_id);
    if (pool_index < n_local) {
      for (auto const bond : p.bonds()) {
        auto const partner_ids = bond.partner_ids();
        if (partner_ids.empty()) {
          continue;
        }
        if (partner_ids.size() == 1u) {
          ++pair_count;
        } else if (partner_ids.size() == 2u) {
          ++angle_count;
        } else if (partner_ids.size() == 3u) {
          ++dihedral_count;
        }
      }
    }
  }
  set_local_bond_numbers(pair_count, angle_count, dihedral_count);
  m_bond_state->allocate();
  m_cached_max_local_particle_id = max_id;
  m_num_local_particles_cached = unique_particles.size();

  // Build the pack-index -> store-row translation view (phase 3.5). The store
  // must already be synchronized (rows assigned) by the caller. The local
  // prefix is the identity: enumerate_local_particles' exclusive-scan order
  // matches ensure_particle_store_synchronized's local loop, so pack index i
  // equals store row i for every local particle. Only the deduped ghost tail
  // needs a real translation.
  assert(not m_particle_store.is_dirty());
  auto const n_unique = unique_particles.size();
  if (m_pack_index_to_store_row.extent(0) != n_unique) {
    m_pack_index_to_store_row = Kokkos::View<int *, Kokkos::HostSpace>(
        Kokkos::ViewAllocateWithoutInitializing("pack_index_to_store_row"),
        n_unique);
  }
  for (std::size_t i = 0u; i < n_unique; ++i) {
    m_pack_index_to_store_row(i) = unique_particles[i]->store_row();
    // local prefix must be the identity (see comment above)
    assert(i >= n_local or m_pack_index_to_store_row(i) == static_cast<int>(i));
  }
  static_cast<void>(n_local);
}

void CellStructure::bind_pack_store_views() {
  auto &aosoa = *m_aosoa;
  aosoa.position = m_particle_store.position_view();
  aosoa.image = m_particle_store.image_box_view();
  aosoa.velocity = m_particle_store.velocity_view();
  // phase 5: id/mass/charge alias the authoritative ParticleStore scalar
  // columns (indexed by store row), like position/velocity, and are read only
  // on cold paths (bond kernels; charge is needed by BondedCoulomb even with no
  // coulomb solver, so it must always be valid). phase-5 perf recovery: type is
  // a pack-owned contiguous column written on rebuild in commit_particle, and
  // pair_charge/pair_dipm are pack-owned contiguous hot-path columns refreshed
  // per step by refresh_pack_charges/refresh_pack_dipm when the corresponding
  // solver is active; those pack-owned columns are not bound here.
  aosoa.id = m_particle_store.id_view();
#ifdef ESPRESSO_ELECTROSTATICS
  aosoa.charge = m_particle_store.q_view();
#endif
#ifdef ESPRESSO_MASS
  aosoa.mass = m_particle_store.mass_view();
#endif
  aosoa.row_map = m_pack_index_to_store_row;
#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
  aosoa.director = m_director_view;
#endif
}

#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
void CellStructure::update_director_view() {
#ifdef ESPRESSO_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  auto const n_total = m_particle_store.number_of_particles();
  if (m_director_view.extent(0) != n_total) {
    m_director_view =
        Kokkos::View<double *[3], ParticleStore::StateVectorLayout,
                     Kokkos::HostSpace>(
            Kokkos::ViewAllocateWithoutInitializing("director"), n_total);
  }
  auto quaternion = m_particle_store.quaternion_view();
  auto director = m_director_view;
  Kokkos::parallel_for(
      "update_director_view", n_total, [quaternion, director](std::size_t row) {
        Utils::Quaternion<double> q;
        q[0] = quaternion(row, 0);
        q[1] = quaternion(row, 1);
        q[2] = quaternion(row, 2);
        q[3] = quaternion(row, 3);
        auto const d = Utils::convert_quaternion_to_director(q);
        director(row, 0) = d[0];
        director(row, 1) = d[1];
        director(row, 2) = d[2];
      });
  Kokkos::fence();
}
#endif

CellStructure::CellStructure(BoxGeometry const &box)
    : m_decomposition{std::make_unique<AtomDecomposition>(box)} {
  mark_particle_store_dirty();
}

void CellStructure::ensure_particle_store_synchronized() {
  if (not m_particle_store.is_dirty()) {
    return;
  }
  // Wire every cell to the store (cheap; idempotent) so that Cell::particles()
  // and Cell::store() work throughout the rebuild and the following step.
  auto wire_stores = [this](std::span<Cell *const> cells) {
    for (auto *cell : cells) {
      cell->set_store(m_particle_store);
    }
  };
  wire_stores(decomposition().local_cells());
  wire_stores(decomposition().ghost_cells());

  // A cell's content after the flip = its surviving committed rows (preserved
  // by old row) + its staged detached particles (seeded from carriers). Count
  // both so begin_rebuild sizes the columns for the new generation.
  auto cell_count = [](Cell const *cell) {
    return cell->rows().size() + cell->staged().size();
  };
  std::size_t n_local = 0u;
  for (auto const *cell : decomposition().local_cells()) {
    n_local += cell_count(cell);
  }
  std::size_t n_ghost = 0u;
  for (auto const *cell : decomposition().ghost_cells()) {
    n_ghost += cell_count(cell);
  }

  // Snapshot the OLD row bags before begin_rebuild swaps the columns: a
  // surviving row is assigned by handing assign_row a Particle attached to
  // THIS store at the OLD row, from which it preserves the (now-spare) column
  // data; the staged particles are detached, so assign_row seeds them from
  // their carriers. begin_rebuild resizes to the new count, so we must capture
  // the old row indices first (a shrink would put them out of the new bounds).
  m_particle_store.begin_rebuild(n_local, n_ghost);
  int row = 0;
  // One reused survivor view rebound per old row (phase 7a perf): constructing
  // a fresh Particle per survivor default-constructs all migration carriers
  // (heap-owning BondList/exclusion members) every iteration; rebinding the two
  // store-handle fields via attach_to_store instead costs two writes and never
  // touches the carriers (assign_row's preserve branch reads only columns by
  // old row). Reused across all cells (assign_cell is called sequentially).
  Particle survivor;
  auto assign_cell = [this, &row, &survivor](Cell *cell,
                                             std::vector<int> const &old_rows) {
    auto &rows = cell->rows();
    rows.clear();
    // Surviving committed rows: preserve column data by old row.
    for (int const old_row : old_rows) {
      survivor.attach_to_store(m_particle_store, old_row);
      m_particle_store.assign_row(survivor, row);
      rows.insert(row);
      ++row;
    }
    // Staged (newly inserted / migrated / fresh ghost) particles: seed from
    // carriers, then clear the staging area.
    for (auto &staged : cell->staged()) {
      m_particle_store.assign_row(staged, row);
      rows.insert(row);
      ++row;
    }
    cell->staged().clear();
  };
  // Reuse ONE scratch buffer across all cells (phase 7a perf): a surviving
  // cell's row indices must be copied out before assign_cell clears/refills the
  // bag, but allocating a fresh std::vector per cell (thousands per resort)
  // churned the heap in the resort hot path. Reusing a single buffer (assign to
  // clear-and-refill without freeing capacity) removes that per-cell alloc.
  std::vector<int> old_rows;
  auto assign_cells = [&assign_cell, &old_rows](std::span<Cell *const> cells) {
    for (auto *cell : cells) {
      auto const &row_bag = cell->rows();
      old_rows.assign(row_bag.begin(), row_bag.end());
      assign_cell(cell, old_rows);
    }
  };
  assign_cells(decomposition().local_cells());
  assign_cells(decomposition().ghost_cells());
  m_particle_store.finish_rebuild();

  // Refresh the id->view index / stable view pool from the freshly assigned
  // rows (locals win over ghost copies of the same id).
  rebuild_particle_index();

#ifdef ESPRESSO_ADDITIONAL_CHECKS
  // Verify each cell's row bag matches the store: the store id at each recorded
  // row equals the id of the view at the same position. Covers local + ghost.
  auto check_cell_rows = [this](std::span<Cell *const> cells) {
    for (auto *cell : cells) {
      auto const &rows = cell->rows();
      std::size_t k = 0u;
      for (auto const &p : cell->particles()) {
        auto const stored_id = m_particle_store.id(rows.begin()[k]);
        if (stored_id != p.id()) {
          throw std::runtime_error(
              "CellRows id mismatch at position " + std::to_string(k) +
              ": store.id(row) = " + std::to_string(stored_id) +
              " but view id = " + std::to_string(p.id()));
        }
        ++k;
      }
      if (k != rows.size()) {
        throw std::runtime_error("CellRows size mismatch");
      }
    }
  };
  check_cell_rows(decomposition().local_cells());
  check_cell_rows(decomposition().ghost_cells());
#endif
}

int CellStructure::stage_row(int const live_row) {
  // Grow the staging store when the next row would overflow its capacity. The
  // staging store's own rows must survive the growth, so a fresh larger store
  // is built and every already-staged row is copied into it via copy_row (the
  // same machinery the flip uses); the small store is then swapped in. Capacity
  // doubles (min 8) to amortize the copies over a batch of stages.
  auto const needed = static_cast<std::size_t>(m_staging_store_next_row) + 1u;
  if (needed > m_staging_store_capacity) {
    auto new_capacity =
        std::max<std::size_t>(m_staging_store_capacity * 2u, 8u);
    while (new_capacity < needed) {
      new_capacity *= 2u;
    }
    ParticleStore grown{};
    grown.begin_rebuild(new_capacity, 0u);
    grown.finish_rebuild();
    for (int r = 0; r < m_staging_store_next_row; ++r) {
      grown.copy_row(m_staging_store, r, r);
    }
    using std::swap;
    swap(m_staging_store, grown);
    m_staging_store_capacity = new_capacity;
  }
  auto const staging_row = m_staging_store_next_row++;
  m_staging_store.copy_row(m_particle_store, live_row, staging_row);
  return staging_row;
}

void CellStructure::rebuild_particle_index() {
  assert(not m_particle_store.is_dirty());
  m_particle_index.clear();
  // Phase 7a perf fix: rewind the fill cursor and REUSE the pool's constructed
  // slots (acquire_view_slot rebinds them in place) instead of clear + N-fold
  // Particle construction. Idle slots past m_view_pool_fill stay as reusable
  // capacity; they are never read.
  m_view_pool_fill = 0u;
  auto const n_local = m_particle_store.number_of_local_particles();
  // Index locals (rows [0, n_local)); their id columns are valid.
  for (std::size_t r = 0u; r < n_local; ++r) {
    auto *view = acquire_view_slot(static_cast<int>(r));
    update_particle_index(view->id(), view);
  }
  // Also index the ALREADY-VALID ghost rows: a bare rebuild (e.g. an
  // add_particle commit, the clean-store forces.cpp sync, or resort_particles'
  // end-of-function sync) must not drop the ghost index entries a prior
  // ghosts_update established, or a ghost-id lookup (LB coupling's
  // is_ghost_for_local_particle, collision detection) would see nullptr. FRESH
  // ghost rows just created by resize_ghost_storage carry a default id (-1) and
  // are skipped by index_ghost_particles; they are picked up by a later
  // index_ghost_particles once ghosts_update fills their id columns.
  index_ghost_particles();
}

void CellStructure::index_ghost_particles() {
  assert(not m_particle_store.is_dirty());
  auto const n_total = m_particle_store.number_of_particles();
  auto const n_local = m_particle_store.number_of_local_particles();
  // Ghost rows [n_local, n_total): add a pool view + index entry for each ghost
  // whose id is VALID (>= 0; a freshly resized ghost carries -1 until
  // ghosts_update fills it) and not already owned by a local (locals win). The
  // local prefix of the pool is left intact; this appends/extends the ghost
  // tail. Idempotent: a ghost id already present (from a local or an earlier
  // ghost of the same id) is skipped, so re-running after ghosts_update only
  // adds the newly-valid ghosts.
  //
  // Phase 7a perf fix: peek at the next idle slot bound to this row, check its
  // id, and only COMMIT (advance m_view_pool_fill via acquire_view_slot) if the
  // ghost is accepted -- a rejected candidate leaves the slot idle for the next
  // row to reuse, so no Particle is constructed per skipped ghost.
  for (std::size_t r = n_local; r < n_total; ++r) {
    // Bind the next idle slot to row r without advancing the cursor. Grow the
    // pool by one only if there is no idle slot to rebind.
    if (m_view_pool_fill >= m_view_pool.size()) {
      m_view_pool.push_back(m_particle_store.make_view(static_cast<int>(r)));
    } else {
      m_view_pool[m_view_pool_fill].attach_to_store(m_particle_store,
                                                    static_cast<int>(r));
    }
    auto const id = m_view_pool[m_view_pool_fill].id();
    if (id < 0 or get_local_particle(id) != nullptr) {
      continue; // leave the slot idle; the next row rebinds it
    }
    update_particle_index(id, std::addressof(m_view_pool[m_view_pool_fill]));
    ++m_view_pool_fill; // commit
  }
}

void CellStructure::check_particle_index() const {
  auto const max_id = get_max_local_particle_id();

  for (auto const &p : local_particles()) {
    auto const id = p.id();

    if (id < 0 || id > max_id) {
      throw std::runtime_error("Particle id out of bounds.");
    }

    // Phase 7a: particles are views over store rows, so the index entry (a
    // pointer into the view pool) and the iterated view have DIFFERENT
    // addresses even for the same particle. Check identity by store row instead
    // of by address: both must resolve to the same row of the same store.
    auto const *indexed = get_local_particle(id);
    if (indexed == nullptr or indexed->store() != p.store() or
        indexed->store_row() != p.store_row()) {
      throw std::runtime_error("Invalid local particle index entry.");
    }
  }

  /* checks: local particle id. Phase 7a: the index also holds ghost entries
   * (for ghost-id lookups), so count only the LOCAL (non-ghost) entries when
   * comparing against local_particles(). */
  std::size_t local_part_cnt = 0u;
  for (int n = 0; n < get_max_local_particle_id() + 1; n++) {
    auto const *indexed = get_local_particle(n);
    if (indexed != nullptr) {
      if (indexed->id() != n) {
        throw std::runtime_error("local_particles part has corrupted id.");
      }
      if (not indexed->is_ghost()) {
        local_part_cnt++;
      }
    }
  }

  if (local_part_cnt != local_particles().size()) {
    throw std::runtime_error(
        std::to_string(local_particles().size()) + " parts in cells but " +
        std::to_string(local_part_cnt) + " parts in local_particles");
  }
}

void CellStructure::check_particle_sorting() const {
  for (auto cell : decomposition().local_cells()) {
    for (auto const &p : cell->particles()) {
      if (particle_to_cell(p) != cell) {
        throw std::runtime_error("misplaced particle with id " +
                                 std::to_string(p.id()));
      }
    }
  }
}

void CellStructure::remove_particle(int id) {
  auto remove_all_bonds_to = [id](BondList &bl) {
    for (auto it = bl.begin(); it != bl.end();) {
      if (Utils::contains(it->partner_ids(), id)) {
        it = bl.erase(it);
      } else {
        std::advance(it, 1);
      }
    }
  };

  for (auto cell : decomposition().local_cells()) {
    cell->set_store(m_particle_store);
    // Iterate the committed rows by position. extract_row removes a row via
    // swap-with-back (so we re-examine the swapped-in position); a kept
    // particle has any bond to the removed id stripped (written through the
    // view into the store's bond sidecar, preserved by the next rebuild).
    for (std::size_t index = 0u; index < cell->rows().size();) {
      auto view = m_particle_store.make_view(cell->rows().begin()[index]);
      if (view.id() == id) {
        CellParticleStorage::extract_row(*cell, index);
      } else {
        remove_all_bonds_to(view.bonds());
        ++index;
      }
    }
  }
  // NOTE: deferred id→view reindex. Between this call and the next
  // ensure_particle_store_synchronized, get_local_particle(id) returns a
  // pool pointer to the now-orphaned row rather than nullptr.  The window is
  // benign in every current call path: System::on_particle_change() forces a
  // resort (which rebuilds the index) before the next force calculation, and
  // particle_node.cpp erases `id` from the particle_node map immediately after
  // remove_particle returns, so the orphaned pointer is never reached through
  // the public API in the interim.
  mark_particle_store_dirty();
}

Particle *CellStructure::add_local_particle(Particle &&p) {
  auto const sort_cell = particle_to_cell(p);
  if (sort_cell) {
    sort_cell->set_store(m_particle_store);
    auto const id = p.id();
    append_staged_particle(*sort_cell, std::move(p));
    mark_particle_store_dirty();
    // Commit immediately so the new particle is live (indexed, iterable,
    // readable) right after this call -- the pre-flip contract callers rely on.
    // Return a stable pointer into the view pool (phase 7a).
    ensure_particle_store_synchronized();
    return get_local_particle(id);
  }

  return {};
}

Particle *CellStructure::add_particle(Particle &&p) {
  auto const sort_cell = particle_to_cell(p);
  /* There is always at least one cell, so if the particle
   * does not belong to a cell on this node we can put it there. */
  auto cell = sort_cell ? sort_cell : decomposition().local_cells()[0];
  cell->set_store(m_particle_store);

  /* If the particle isn't local a global resort may be
   * needed, otherwise a local resort if sufficient. */
  set_resort_particles(sort_cell ? Cells::RESORT_LOCAL : Cells::RESORT_GLOBAL);

  auto const id = p.id();
  append_staged_particle(*cell, std::move(p));
  mark_particle_store_dirty();
  // Commit immediately so the new particle is live right after this call (the
  // pre-flip contract). Return a stable view-pool pointer (phase 7a).
  ensure_particle_store_synchronized();
  return get_local_particle(id);
}

int CellStructure::get_max_local_particle_id() const {
  auto it = std::ranges::find_if(std::ranges::views::reverse(m_particle_index),
                                 [](auto const *p) { return p != nullptr; });

  return (it != m_particle_index.rend()) ? (*it)->id() : -1;
}

int CellStructure::get_local_pair_bond_numbers() const {
  return m_bond_state->pair_count;
}
int CellStructure::get_local_angle_bond_numbers() const {
  return m_bond_state->angle_count;
}
int CellStructure::get_local_dihedral_bond_numbers() const {
  return m_bond_state->dihedral_count;
}
void CellStructure::set_local_bond_numbers(int pair_value, int angle_value,
                                           int dihedral_value) {
  m_bond_state->set_counts(pair_value, angle_value, dihedral_value);
}
#ifdef ESPRESSO_COLLISION_DETECTION
void CellStructure::clear_new_bonds() { m_bond_state->clear_new_bonds(); }
void CellStructure::add_new_bond(int bond_id,
                                 std::vector<int> const &particle_ids) {
  m_bond_state->add_new_bond(bond_id, particle_ids, get_id_to_index());
}
void CellStructure::rebuild_bond_list() { m_bond_state->rebuild(); }
#endif // ESPRESSO_COLLISION_DETECTION

void CellStructure::remove_all_particles() {
  for (auto cell : decomposition().local_cells()) {
    CellParticleStorage::clear_particles(*cell);
  }
  // Also clear the GHOST cells (phase 7a): they hold ghost copies of the
  // now-deleted particles. If left in place, the next store rebuild would
  // re-index those stale ghost rows (rebuild_particle_index indexes valid
  // ghosts), so get_local_particle would still report the deleted ids as
  // present -- e.g. a fresh add of the same id would wrongly see "already
  // exists". Emptying them keeps the index consistent with the emptied system.
  for (auto cell : decomposition().ghost_cells()) {
    CellParticleStorage::clear_particles(*cell);
  }

  clear_particle_index();
  clear_bond_properties();
  mark_particle_store_dirty();
}

/* Map the data parts flags from cells to those used internally
 * by the ghost communication */
unsigned map_data_parts(unsigned data_parts) {
  using namespace Cells;

  /* clang-format off */
  return GHOSTTRANS_NONE
         | ((data_parts & DATA_PART_PROPERTIES) ? GHOSTTRANS_PROPRTS : 0u)
         | ((data_parts & DATA_PART_POSITION) ? GHOSTTRANS_POSITION : 0u)
         | ((data_parts & DATA_PART_MOMENTUM) ? GHOSTTRANS_MOMENTUM : 0u)
         | ((data_parts & DATA_PART_FORCE) ? GHOSTTRANS_FORCE : 0u)
#ifdef ESPRESSO_BOND_CONSTRAINT
         | ((data_parts & DATA_PART_RATTLE) ? GHOSTTRANS_RATTLE : 0u)
#endif
         | ((data_parts & DATA_PART_BONDS) ? GHOSTTRANS_BONDS : 0u);
  /* clang-format on */
}

void CellStructure::ghosts_count() {
  ghost_communicator(decomposition().exchange_ghosts_comm(),
                     *get_system().box_geo, GHOSTTRANS_PARTNUM);
  mark_particle_store_dirty();
}
void CellStructure::ghosts_update(unsigned data_parts) {
  ghost_communicator(decomposition().exchange_ghosts_comm(),
                     *get_system().box_geo, map_data_parts(data_parts));
}
void CellStructure::ghosts_reduce_forces() {
  ghost_communicator(decomposition().collect_ghost_force_comm(),
                     *get_system().box_geo, GHOSTTRANS_FORCE);
}
#ifdef ESPRESSO_BOND_CONSTRAINT
void CellStructure::ghosts_reduce_rattle_correction() {
  ghost_communicator(decomposition().collect_ghost_force_comm(),
                     *get_system().box_geo, GHOSTTRANS_RATTLE);
}
#endif

void CellStructure::resort_particles(bool global_flag, bool commit) {
  // Phase 7a: the id->view index / view pool is rebuilt wholesale from the
  // store by ensure_particle_store_synchronized after the resort (see
  // rebuild_particle_index), so the pre-flip incremental index bookkeeping
  // (invalidate_ghosts + the ModifiedList/RemovedParticle diff visitor) is no
  // longer needed. The `diff` is still collected (the decomposition contract)
  // but only for its "cells touched" side effect; index correctness comes from
  // the rebuild. Clear the stale index now so nothing reads a dangling view
  // pool pointer in the resort window (get_local_particle is not called there).
  clear_particle_index();

  std::vector<ParticleChange> diff;

  // Ensure every cell knows the store before the decomposition mutates cell
  // storage (extract_row / insert / ghost resize all need Cell::store()).
  for (auto *cell : m_decomposition->local_cells()) {
    cell->set_store(m_particle_store);
  }
  for (auto *cell : m_decomposition->ghost_cells()) {
    cell->set_store(m_particle_store);
  }

  // Mark the store dirty BEFORE the decomposition's resort: the dirty flag must
  // be truthful DURING the resort window. HybridDecomposition runs an internal
  // PARTNUM+ghost update inside that window, and the columnar ghost paths must
  // see dirty and fall back to the per-particle path there.
  mark_particle_store_dirty();
  m_decomposition->resort(global_flag, diff);
  static_cast<void>(diff);

  // Commit the staged (migrated/new) particles into store rows and rebuild the
  // id->view index (phase 7a): resort cleared the index and staged the moved
  // particles, so a caller reading the index right after resort (e.g. the
  // ADDITIONAL_CHECKS below, or a direct/unit-test resort_particles call) must
  // see a consistent, populated index -- the pre-flip contract. Rank-local (no
  // MPI), so it does not affect the collective ordering the caller relies on.
  //
  // Phase 7a perf: the hot-path caller (update_ghosts_and_resort_particle)
  // passes commit=false to DEFER this to a single rebuild AFTER ghosts_count,
  // so locals are not committed here only to be re-copied by a second whole-
  // store rebuild moments later. ghosts_count sizes ghost cells from the
  // SOURCE cell's Cell::size (committed rows + staged), so it does not need the
  // locals committed first; and nothing reads the index in the deferred window.
  if (commit) {
    ensure_particle_store_synchronized();
  }

  auto const &lebc = get_system().box_geo->lees_edwards_bc();
  m_rebuild_verlet_list = true;
  m_rebuild_verlet_list_cabana = true;
  m_le_pos_offset_at_last_resort = lebc.pos_offset;

#ifdef ESPRESSO_ADDITIONAL_CHECKS
  // These checks read the id->view index / sorted store, which only exist after
  // the commit above; when the commit is deferred, the caller runs its own sync
  // and the equivalent checks (ensure_particle_store_synchronized's cell-row
  // validation) after ghosts_count.
  if (commit) {
    check_particle_index();
    check_particle_sorting();
  }
#endif
}

void CellStructure::set_atom_decomposition() {
  auto &system = get_system();
  auto &local_geo = *system.local_geo;
  auto const &box_geo = *system.box_geo;
  set_particle_decomposition(
      std::make_unique<AtomDecomposition>(::comm_cart, box_geo));
  m_type = CellStructureType::NSQUARE;
  local_geo.set_cell_structure_type(m_type);
  system.on_cell_structure_change();
}

void CellStructure::set_regular_decomposition(
    double range, std::optional<std::pair<int, int>> fully_connected_boundary) {
  auto &system = get_system();
  auto &local_geo = *system.local_geo;
  auto const &box_geo = *system.box_geo;
  set_particle_decomposition(std::make_unique<RegularDecomposition>(
      ::comm_cart, range, box_geo, local_geo, fully_connected_boundary));
  m_type = CellStructureType::REGULAR;
  local_geo.set_cell_structure_type(m_type);
  system.on_cell_structure_change();
}

void CellStructure::set_hybrid_decomposition(double cutoff_regular,
                                             std::set<int> n_square_types) {
  auto &system = get_system();
  auto &local_geo = *system.local_geo;
  auto const &box_geo = *system.box_geo;
  auto hybrid = std::make_unique<HybridDecomposition>(
      ::comm_cart, cutoff_regular, m_verlet_skin,
      [&system]() { return system.get_global_ghost_flags(); }, box_geo,
      local_geo, n_square_types);
  // Phase 7a: let the hybrid decomposition commit staged particles to store
  // rows between its internal resort and ghost communications (see
  // HybridDecomposition::resort). Mark the store dirty first: the ghost resize
  // (resize_ghost_storage) stages ghosts WITHOUT marking dirty, so a bare
  // ensure_particle_store_synchronized would early-return (clean store) and
  // never commit them; forcing the rebuild guarantees the staged particles
  // (migrated locals AND freshly resized ghosts) become committed rows.
  hybrid->set_commit_store([this]() {
    mark_particle_store_dirty();
    ensure_particle_store_synchronized();
  });
  set_particle_decomposition(std::move(hybrid));
  m_type = CellStructureType::HYBRID;
  local_geo.set_cell_structure_type(m_type);
  system.on_cell_structure_change();
}

void CellStructure::set_verlet_skin(double value) {
  assert(value >= 0.);
  m_verlet_skin = value;
  m_verlet_skin_set = true;
  m_rebuild_verlet_list_cabana = true;
  get_system().on_verlet_skin_change();
}

void CellStructure::set_verlet_skin_heuristic() {
  assert(not is_verlet_skin_set());
  auto const max_cut = get_system().maximal_cutoff();
  if (max_cut <= 0.) {
    throw std::runtime_error(
        "cannot automatically determine skin, please set it manually");
  }
  /* maximal skin that can be used without resorting is the maximal
   * range of the cell system minus what is needed for interactions. */
  auto const max_range = std::ranges::min(max_cutoff());
  auto const new_skin = std::min(0.4 * max_cut, max_range - max_cut);
  set_verlet_skin(new_skin);
}

void CellStructure::update_ghosts_and_resort_particle(unsigned data_parts) {
  /* data parts that are only updated on resort */
  auto constexpr resort_only_parts =
      Cells::DATA_PART_PROPERTIES | Cells::DATA_PART_BONDS;

  auto const global_resort = boost::mpi::all_reduce(
      ::comm_cart, m_resort_particles, std::bit_or<unsigned>());

  if (global_resort != Cells::RESORT_NONE) {
    auto const do_global_resort = (global_resort & Cells::RESORT_GLOBAL) != 0;

    /* Resort cell system. Phase 7a perf: defer the local commit (commit=false)
     * so it is folded into the SINGLE post-ghosts_count store rebuild below,
     * instead of committing locals here and re-copying every local column a
     * second time in that rebuild (the phase-7a double-rebuild regression). */
    resort_particles(do_global_resort, /*commit=*/false);
    /* Phase 7a: after resort the migrated/new locals are STAGED (not yet
     * committed to store rows). ghosts_count sizes every ghost cell from its
     * source (local) cell's Cell::size, which counts committed rows + staged,
     * so staged locals are counted correctly and do not need committing first.
     * Count ghosts (stages ghost slots + marks the store dirty), then rebuild
     * ONCE to commit locals AND the freshly staged ghosts together. */
    ghosts_count();
    /* Rebuild the store rows NOW: resort staged locals and ghosts_count marked
     * the store dirty (staged ghost slots); a dirty store forces the ghost
     * update below onto the slow per-field accessor path (measured +49% in this
     * slot on lj-4rank). After this single rebuild, locals and freshly created
     * ghosts are all attached and the update writes their state straight into
     * the columns by row. */
    ensure_particle_store_synchronized();
    ghosts_update(data_parts);

    /* Index the ghosts now that ghosts_update has filled their id columns
     * (phase 7a: the rebuild inside ensure_particle_store_synchronized indexed
     * only locals, whose ids were valid; freshly created ghost rows carried a
     * default id until this update). Adds one stable view-pool entry per ghost
     * whose id is not owned by a local. */
    index_ghost_particles();

    /* Particles are now sorted */
    clear_resort_particles();
  } else {
    /* Communication step: ghost information */
    ghosts_update(data_parts & ~resort_only_parts);
  }
}

bool CellStructure::check_resort_required(
    Utils::Vector3d const &additional_offset) const {
  auto const lim = Utils::sqr(m_verlet_skin / 2.) - additional_offset.norm2();

  Reduction::AddPartialResultKernel<bool> add_partial =
      [lim](bool &result, Particle const &p) {
        if ((Utils::Vector3d(p.pos()) -
             Utils::Vector3d(p.pos_at_last_verlet_update()))
                .norm2() > lim) {
          result = true;
        }
      };

  Reduction::ReductionOp<bool> reduce_op = [](bool &acc, bool const &val) {
    acc |= val;
  };

  return reduce_over_local_particles(*this, add_partial, reduce_op);
}
