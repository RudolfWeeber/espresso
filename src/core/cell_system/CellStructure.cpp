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
  auto &unique_particles = m_unique_particles;
  unique_particles.clear();
  unique_particles.resize(count_local_particles());
  std::unordered_set<int> registered_index{};
  using execution_space = Kokkos::DefaultExecutionSpace;
  int n_threads = execution_space().concurrency();
  std::vector<int> max_ids(n_threads);

  m_bond_state->reset_counts();
  std::vector<int> pair_counts(n_threads, 0);
  std::vector<int> angle_counts(n_threads, 0);
  std::vector<int> dihedral_counts(n_threads, 0);

  enumerate_local_particles(
      *this, [&unique_particles, &max_ids, &pair_counts, &angle_counts,
              &dihedral_counts](std::size_t index, Particle &p) {
        unique_particles[index] = &p;
        auto const thread_num = omp_get_thread_num();
        max_ids[thread_num] = std::max(p.id(), max_ids[thread_num]);
        for (auto const bond : p.bonds()) {
          if (not bond.partner_ids().empty()) {
            auto const partner_ids = bond.partner_ids();
            if (partner_ids.size() == 1u) {
              pair_counts[thread_num] += 1;
            } else if (partner_ids.size() == 2u) {
              angle_counts[thread_num] += 1;
            } else if (partner_ids.size() == 3u) {
              dihedral_counts[thread_num] += 1;
            }
          }
        }
      });
  Kokkos::fence();
  int pair_count = std::reduce(std::begin(pair_counts), std::end(pair_counts));
  int angle_count =
      std::reduce(std::begin(angle_counts), std::end(angle_counts));
  int dihedral_count =
      std::reduce(std::begin(dihedral_counts), std::end(dihedral_counts));
  set_local_bond_numbers(pair_count, angle_count, dihedral_count);
  m_bond_state->allocate();

  int max_id = *(std::max_element(max_ids.begin(), max_ids.end()));
  for (auto &p : ghost_particles()) {
    auto const *local_particle = get_local_particle(p.id());
    if (not local_particle) {
      continue;
    }
    if (not local_particle->is_ghost()) {
      continue;
    }
    if (registered_index.contains(p.id())) {
      continue;
    }
    registered_index.insert(p.id());
    unique_particles.emplace_back(&p);
    max_id = std::max(p.id(), max_id);
  }
  registered_index.clear();
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
  auto const n_local = count_local_particles();
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
  auto const n_local = count_local_particles();
  std::size_t n_ghost = 0u;
  for (auto const &p : ghost_particles()) {
    static_cast<void>(p);
    ++n_ghost;
  }
  m_particle_store.begin_rebuild(n_local, n_ghost);
  int row = 0;
  // Iterate cells directly (not the flattened particle ranges) so each cell's
  // parallel row-index bag (phase 7a, dormant) can be refilled as rows are
  // assigned: cleared at the start of the cell's block, then one int appended
  // per particle in ParticleList order. This preserves the row-assignment
  // order (local cells then ghost cells, particles in Bag order) exactly.
  auto assign_cell_rows = [this, &row](std::span<Cell *const> cells) {
    for (auto *cell : cells) {
      auto &rows = cell->rows();
      rows.clear();
      for (auto &p : cell->particles()) {
        m_particle_store.assign_row(p, row);
        rows.insert(row);
        ++row;
      }
    }
  };
  assign_cell_rows(decomposition().local_cells());
  assign_cell_rows(decomposition().ghost_cells());
  m_particle_store.finish_rebuild();

#ifdef ESPRESSO_ADDITIONAL_CHECKS
  // Verify the dormant per-cell row bags match the cell Bag contents: same
  // length, and the store id at each recorded row equals the particle id at the
  // same Bag position. Covers both local and ghost cells.
  auto check_cell_rows = [this](std::span<Cell *const> cells) {
    for (auto *cell : cells) {
      auto const &parts = cell->particles();
      auto const &rows = cell->rows();
      if (rows.size() != parts.size()) {
        throw std::runtime_error("CellRows size mismatch: bag has " +
                                 std::to_string(parts.size()) + " parts but " +
                                 std::to_string(rows.size()) + " rows");
      }
      std::size_t k = 0u;
      for (auto const &p : parts) {
        auto const stored_id = m_particle_store.id(rows.begin()[k]);
        if (stored_id != p.id()) {
          throw std::runtime_error(
              "CellRows id mismatch at position " + std::to_string(k) +
              ": store.id(row) = " + std::to_string(stored_id) +
              " but bag id = " + std::to_string(p.id()));
        }
        ++k;
      }
    }
  };
  check_cell_rows(decomposition().local_cells());
  check_cell_rows(decomposition().ghost_cells());
#endif
}

void CellStructure::check_particle_index() const {
  auto const max_id = get_max_local_particle_id();

  for (auto const &p : local_particles()) {
    auto const id = p.id();

    if (id < 0 || id > max_id) {
      throw std::runtime_error("Particle id out of bounds.");
    }

    if (get_local_particle(id) != &p) {
      throw std::runtime_error("Invalid local particle index entry.");
    }
  }

  /* checks: local particle id */
  std::size_t local_part_cnt = 0u;
  for (int n = 0; n < get_max_local_particle_id() + 1; n++) {
    if (get_local_particle(n) != nullptr) {
      local_part_cnt++;
      if (get_local_particle(n)->id() != n) {
        throw std::runtime_error("local_particles part has corrupted id.");
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
    auto &parts = cell->particles();
    for (auto it = parts.begin(); it != parts.end();) {
      if (it->id() == id) {
        it = CellParticleStorage::erase_particle(parts, it);
        update_particle_index(id, nullptr);
        update_particle_index(parts);
      } else {
        remove_all_bonds_to(it->bonds());
        it++;
      }
    }
  }
  mark_particle_store_dirty();
}

Particle *CellStructure::add_local_particle(Particle &&p) {
  auto const sort_cell = particle_to_cell(p);
  if (sort_cell) {
    mark_particle_store_dirty();
    return std::addressof(
        append_indexed_particle(sort_cell->particles(), std::move(p)));
  }

  return {};
}

Particle *CellStructure::add_particle(Particle &&p) {
  auto const sort_cell = particle_to_cell(p);
  /* There is always at least one cell, so if the particle
   * does not belong to a cell on this node we can put it there. */
  auto cell = sort_cell ? sort_cell : decomposition().local_cells()[0];

  /* If the particle isn't local a global resort may be
   * needed, otherwise a local resort if sufficient. */
  set_resort_particles(sort_cell ? Cells::RESORT_LOCAL : Cells::RESORT_GLOBAL);

  mark_particle_store_dirty();
  return std::addressof(
      append_indexed_particle(cell->particles(), std::move(p)));
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
    CellParticleStorage::clear_particles(cell->particles());
  }

  m_particle_index.clear();
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

namespace {
/**
 * @brief Apply a @ref ParticleChange to a particle index.
 */
struct UpdateParticleIndexVisitor {
  CellStructure *cs;

  void operator()(RemovedParticle rp) const {
    cs->update_particle_index(rp.id, nullptr);
  }
  void operator()(ModifiedList mp) const { cs->update_particle_index(mp.pl); }
};
} // namespace

void CellStructure::resort_particles(bool global_flag) {
  invalidate_ghosts();

  std::vector<ParticleChange> diff;

  // Mark the store dirty BEFORE the decomposition's resort: the dirty flag must
  // be truthful DURING the resort window. HybridDecomposition runs an internal
  // PARTNUM+ghost update inside that window, and the columnar ghost paths must
  // see dirty and fall back to the per-particle path there.
  mark_particle_store_dirty();
  m_decomposition->resort(global_flag, diff);

  for (auto d : diff) {
    std::visit(UpdateParticleIndexVisitor{this}, d);
  }

  auto const &lebc = get_system().box_geo->lees_edwards_bc();
  m_rebuild_verlet_list = true;
  m_rebuild_verlet_list_cabana = true;
  m_le_pos_offset_at_last_resort = lebc.pos_offset;

#ifdef ESPRESSO_ADDITIONAL_CHECKS
  check_particle_index();
  check_particle_sorting();
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
  set_particle_decomposition(std::make_unique<HybridDecomposition>(
      ::comm_cart, cutoff_regular, m_verlet_skin,
      [&system]() { return system.get_global_ghost_flags(); }, box_geo,
      local_geo, n_square_types));
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

    /* Resort cell system */
    resort_particles(do_global_resort);
    ghosts_count();
    /* Rebuild the store rows NOW: resort_particles/ghosts_count marked the
     * store dirty, and a dirty store forces the ghost update below onto the
     * slow per-field accessor path (measured +49% in this slot on lj-4rank).
     * After the rebuild, freshly created ghosts are attached and the update
     * writes their state straight into the columns by row. */
    ensure_particle_store_synchronized();
    ghosts_update(data_parts);

    /* Add the ghost particles to the index if we don't already
     * have them. */
    for (auto &p : ghost_particles()) {
      if (get_local_particle(p.id()) == nullptr) {
        update_particle_index(p.id(), &p);
      }
    }

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
