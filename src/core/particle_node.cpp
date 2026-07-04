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

#include "particle_node.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_store/ParticleStore.hpp"
#include "rotation.hpp"
#include "system/System.hpp"

#include <utils/Cache.hpp>
#include <utils/Vector.hpp>
#include <utils/mpi/gatherv.hpp>

#include <boost/mpi/collectives/all_gather.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/gather.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/collectives/scatter.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <memory>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

constexpr auto some_tag = 42;

/** @brief Enable particle type tracking in @ref particle_type_map. */
static bool type_list_enable;

/** @brief Mapping particle types to lists of particle ids. */
static std::unordered_map<int, std::unordered_set<int>> particle_type_map;

/** @brief Mapping particle ids to MPI ranks. */
static std::unordered_map<int, int> particle_node;

static auto &get_cell_structure() {
  return *System::get_system().cell_structure;
}

/**
 * @brief Keep track of the largest particle id.
 * This book-keeping variable is necessary to make particle insertion run
 * in constant time. Traversing the @ref particle_node to find the largest
 * particle id scales with O(N) and traversing the local cells in parallel
 * followed by a reduction scales with O(N^2).
 */
static int max_seen_pid = -1;

static auto rebuild_needed() {
  auto is_rebuild_needed = ::particle_node.empty();
  boost::mpi::broadcast(::comm_cart, is_rebuild_needed, 0);
  return is_rebuild_needed;
}

static void mpi_synchronize_max_seen_pid_local() {
  boost::mpi::broadcast(::comm_cart, ::max_seen_pid, 0);
}

void init_type_map(int type) {
  if (type < 0) {
    throw std::runtime_error("Types may not be negative");
  }
  ::type_list_enable = true;
  auto &nonbonded_ias = *System::get_system().nonbonded_ias;
  nonbonded_ias.make_particle_type_exist(type);

  std::vector<int> local_pids;
  for (auto const &p : get_cell_structure().local_particles()) {
    if (p.type() == type) {
      local_pids.emplace_back(p.id());
    }
  }

  std::vector<std::vector<int>> global_pids;
  boost::mpi::all_gather(::comm_cart, local_pids, global_pids);
  ::particle_type_map[type].clear();
  for (auto const &vec : global_pids) {
    for (auto const &p_id : vec) {
      ::particle_type_map[type].insert(p_id);
    }
  }
}

static void remove_id_from_map(int p_id, int type) {
  auto it = particle_type_map.find(type);
  if (it != particle_type_map.end())
    it->second.erase(p_id);
}

static void add_id_to_type_map(int p_id, int type) {
  auto it = particle_type_map.find(type);
  if (it != particle_type_map.end())
    it->second.insert(p_id);
}

void on_particle_type_change(int p_id, int old_type, int new_type) {
  if (::type_list_enable) {
    if (old_type == type_tracking::any_type) {
      for (auto &kv : ::particle_type_map) {
        auto it = kv.second.find(p_id);
        if (it != kv.second.end()) {
          kv.second.erase(it);
#ifndef NDEBUG
          if (auto p = get_cell_structure().get_local_particle(p_id)) {
            assert(p->type() == kv.first);
          }
#endif
          break;
        }
      }
    } else if (old_type != type_tracking::new_part) {
      if (old_type != new_type) {
        remove_id_from_map(p_id, old_type);
      }
    }
    add_id_to_type_map(p_id, new_type);
  }
}

namespace {
/* Limit cache to 100 MiB */
std::size_t const max_cache_size = (100ul * 1048576ul) / sizeof(Particle);
Utils::Cache<int, Particle> particle_fetch_cache(max_cache_size);

/**
 * @brief Snapshot store backing the fetch cache (migration phase 3).
 *
 * Every detached @ref Particle placed into @ref particle_fetch_cache is a
 * migration-serialized copy from a worker rank: it arrives with no
 * @ref ParticleStore, so its force/torque/state accessors would read the
 * (soon-to-be-removed) sub-structs. To keep ALL accessor reads on cached
 * particles working after the phase-8 flip (cluster analysis, ParticleHandle
 * getters), each cached particle is attached to this head-node-local store,
 * seeded from its migration carriers via @ref ParticleStore::assign_row.
 * Mirrors @ref Constraints::ShapeBasedConstraint::m_part_rep_store.
 *
 * Pointer-stability design (IMPORTANT):
 * @ref Utils::Cache is backed by a @c std::unordered_map<int,const Particle>,
 * whose nodes are stable across rehash but whose @c put() EVICTS a RANDOM
 * other element when the cache is full (@c drop_random_element). We therefore
 * do NOT keep raw @c Particle* addresses to re-assign on growth (an evicted
 * entry's address would dangle). Instead the store has a FIXED capacity per
 * invalidation epoch: it is lazily built once (on the first attach after an
 * invalidation) sized to the cache's capacity, and rows are handed out
 * monotonically. If the row counter would exceed the capacity, we simply drop
 * the whole cache and store (@ref reset_fetch_cache_store) and rebuild lazily
 * on the next attach — a cache is always safe to drop, and this avoids any
 * pointer bookkeeping or store-growth relocation. The capacity is capped so a
 * single epoch never over-allocates. */
ParticleStore fetch_cache_store;
/** Next free row in @ref fetch_cache_store; -1 while the store is unbuilt. */
int fetch_cache_store_next_row = -1;
/** Capacity @ref fetch_cache_store was built with in the current epoch. */
std::size_t fetch_cache_store_capacity = 0u;
/** Co-ownership of the Kokkos runtime (mirrors @ref CellStructure's
 *  @c m_kokkos_handle and @ref ShapeBasedConstraint's). @ref fetch_cache_store
 *  holds Kokkos Views that must be released before @c Kokkos::finalize().
 *  Captured lazily when the store first allocates (particle_node code runs
 *  well after Kokkos init, but we capture defensively the same way). */
std::shared_ptr<KokkosHandle> fetch_cache_store_kokkos_handle;

/** Upper bound on rows allocated per invalidation epoch, so a single epoch
 *  never over-allocates the store. The store never needs more rows than the
 *  cache can hold entries, so this is min(cache capacity, a sane cap). */
std::size_t fetch_cache_store_capacity_target() {
  constexpr std::size_t cap = 65536u;
  return std::min(particle_fetch_cache.max_size(), cap);
}

/** Release the snapshot store's columns and reset the epoch. Kept in lock-step
 *  with @ref particle_fetch_cache invalidation: dropping cached particles whose
 *  store rows we release simultaneously is exactly correct. */
void reset_fetch_cache_store() {
  fetch_cache_store.release_columns();
  fetch_cache_store_next_row = -1;
  fetch_cache_store_capacity = 0u;
  // Drop the Kokkos runtime co-ownership: the store now holds no Views, so it
  // no longer needs to keep Kokkos alive. Holding the handle here would defer
  // Kokkos::finalize() to this translation unit's static teardown, where the
  // store static is destroyed AFTER the handle static (reverse construction
  // order) and its (already released) Views would otherwise be touched post-
  // finalize. Releasing both together keeps them in lock-step.
  fetch_cache_store_kokkos_handle.reset();
}

/** Allocate the fixed-capacity columns for a fresh epoch and reset the row
 *  counter to 0. Co-owns the Kokkos runtime so the columns can be released even
 *  if this translation unit outlives the last CellStructure (Kokkos is
 *  initialized by the time any particle fetch happens; capture defensively). */
void build_fetch_cache_store() {
  fetch_cache_store_kokkos_handle = ::kokkos_handle;
  fetch_cache_store_capacity = fetch_cache_store_capacity_target();
  fetch_cache_store.begin_rebuild(fetch_cache_store_capacity, 0u);
  fetch_cache_store.finish_rebuild();
  fetch_cache_store_next_row = 0;
}

/**
 * @brief Attach a just-cached particle to @ref fetch_cache_store.
 *
 * Grows nothing: the store is a fixed-capacity snapshot for the current
 * invalidation epoch (see the store's doc comment). Builds the store lazily on
 * first use, then hands out one monotonic row per call, seeded from @p p's
 * migration carriers. If the store is exhausted mid-epoch, drops the cache and
 * store and starts a fresh epoch — the caller has already put @p p into the
 * cache, so we ferry it through a copy whose carriers hold the current values
 * (store binding intact), re-insert it into the cleared cache, and attach the
 * re-inserted copy to row 0 of the fresh store.
 *
 * @param p   the (mutable) particle living inside the fetch cache.
 * @returns   a pointer to the now-attached particle in the cache. This is @p p
 *            itself in the common case, but on the exhaustion path the cache is
 *            cleared and @p p is re-inserted, so the returned pointer differs
 *            and callers MUST use it instead of their prior cache pointer.
 */
Particle const *attach_cached_particle(Particle &p) {
  auto const exhausted = fetch_cache_store_next_row >= 0 and
                         static_cast<std::size_t>(fetch_cache_store_next_row) >=
                             fetch_cache_store_capacity;
  if (exhausted) {
    // Dropping the cache is always safe (it is a cache); it also drops @p p, so
    // ferry it through a detached copy (its carriers hold the values) and
    // re-insert into the freshly cleared cache before rebuilding the store.
    Particle kept = p;
    auto const kept_id = kept.id();
    particle_fetch_cache.invalidate();
    reset_fetch_cache_store();
    auto const reinserted = particle_fetch_cache.put(kept_id, std::move(kept));
    build_fetch_cache_store();
    fetch_cache_store.assign_row(const_cast<Particle &>(*reinserted),
                                 fetch_cache_store_next_row++);
    return reinserted;
  }
  if (fetch_cache_store_next_row < 0) {
    build_fetch_cache_store();
  }
  fetch_cache_store.assign_row(p, fetch_cache_store_next_row++);
  return &p;
}
} // namespace

void invalidate_fetch_cache() {
  particle_fetch_cache.invalidate();
  // Cache and snapshot store lifecycles stay in lock-step: rows we release
  // here back cached particles that we are dropping at the same time.
  reset_fetch_cache_store();
}
std::size_t fetch_cache_max_size() { return particle_fetch_cache.max_size(); }

static void mpi_send_particle_data_local(int p_id) {
  // Owner-side serialize reads the particle's accessors to fill the migration
  // carriers, so the owning rank's store must be clean first. O(1) when clean.
  get_cell_structure().ensure_particle_store_synchronized();
  auto const p = get_cell_structure().get_local_particle(p_id);
  auto const found = p and not p->is_ghost();
  assert(1 == boost::mpi::all_reduce(::comm_cart, static_cast<int>(found),
                                     std::plus<>()) &&
         "Particle not found");
  if (found) {
    ::comm_cart.send(0, 42, *p);
  }
}

REGISTER_CALLBACK(mpi_send_particle_data_local)

const Particle &get_particle_data(int p_id) {
  auto const pnode = get_particle_node(p_id);

  if (pnode == this_node) {
    // The local-return path hands out a live particle whose accessor reads need
    // valid ParticleStore rows. O(1) when the store is clean; rank-local.
    get_cell_structure().ensure_particle_store_synchronized();
    auto const p = get_cell_structure().get_local_particle(p_id);
    assert(p != nullptr);
    return *p;
  }

  /* Query the cache */
  auto const p_ptr = particle_fetch_cache.get(p_id);
  if (p_ptr) {
    return *p_ptr;
  }

  /* Cache miss, fetch the particle,
   * put it into the cache and return a pointer into the cache. */
  Communication::mpiCallbacks().call_all(mpi_send_particle_data_local, p_id);
  Particle result{};
  ::comm_cart.recv(boost::mpi::any_source, boost::mpi::any_tag, result);
  auto const p_cached = particle_fetch_cache.put(p_id, std::move(result));
  // Attach the just-cached (detached) copy to the head-node snapshot store so
  // its accessors keep reading valid rows after the phase-8 flip. The cache
  // value is const-qualified; attaching only sets the particle's store row (a
  // rank-local, transitional book-keeping field), so mutating through the cast
  // is intentional and confined here. On the exhaustion path the cache is
  // cleared and the particle re-inserted, so we return the pointer that attach
  // reports (which may differ from p_cached).
  return *attach_cached_particle(const_cast<Particle &>(*p_cached));
}

static auto
get_local_particle_property(int p_id,
                            Utils::Vector3d (*getter)(Particle const &)) {
  // Ensure the owning rank's particle has a valid ParticleStore row before the
  // getter reads force/torque. O(1) when the store is clean; rank-local.
  get_cell_structure().ensure_particle_store_synchronized();
  auto const p = get_cell_structure().get_local_particle(p_id);
  auto const found = (p != nullptr) and not p->is_ghost();
  assert(1 == boost::mpi::all_reduce(::comm_cart, static_cast<int>(found),
                                     std::plus<>()) &&
         "particle not found exactly once");
  auto const local_value = found ? getter(*p) : Utils::Vector3d{};
  return boost::mpi::all_reduce(::comm_cart, local_value, std::plus<>());
}

static void mpi_get_particle_force_local(int p_id) {
  get_local_particle_property(
      p_id, [](Particle const &p) { return Utils::Vector3d(p.force()); });
}

REGISTER_CALLBACK(mpi_get_particle_force_local)

Utils::Vector3d get_particle_force(int p_id) {
  Communication::mpiCallbacks().call(mpi_get_particle_force_local, p_id);
  return get_local_particle_property(
      p_id, [](Particle const &p) { return Utils::Vector3d(p.force()); });
}

#ifdef ESPRESSO_ROTATION
static void mpi_get_particle_torque_lab_local(int p_id) {
  get_local_particle_property(p_id, [](Particle const &p) {
    return convert_vector_body_to_space(p, Utils::Vector3d(p.torque()));
  });
}

REGISTER_CALLBACK(mpi_get_particle_torque_lab_local)

Utils::Vector3d get_particle_torque_lab(int p_id) {
  Communication::mpiCallbacks().call(mpi_get_particle_torque_lab_local, p_id);
  return get_local_particle_property(p_id, [](Particle const &p) {
    return convert_vector_body_to_space(p, Utils::Vector3d(p.torque()));
  });
}
#endif // ESPRESSO_ROTATION

static void mpi_get_particles_local() {
  std::vector<int> local_ids;
  boost::mpi::scatter(comm_cart, local_ids, 0);

  std::vector<Particle> parts(local_ids.size());
  std::ranges::transform(local_ids, parts.begin(), [](int p_id) {
    auto const p = get_cell_structure().get_local_particle(p_id);
    assert(p != nullptr);
    return *p;
  });

  Utils::Mpi::gatherv(comm_cart, parts.data(), static_cast<int>(parts.size()),
                      0);
}

REGISTER_CALLBACK(mpi_get_particles_local)

/**
 * @brief Get multiple particles at once.
 *
 * *WARNING* Particles are returned in an arbitrary order.
 *
 * @param ids The ids of the particles that should be returned.
 *
 * @returns The particle list.
 */
static std::vector<Particle> mpi_get_particles(std::span<const int> ids) {
  Communication::mpiCallbacks().call(mpi_get_particles_local);
  /* Return value */
  std::vector<Particle> parts(ids.size());

  /* Group ids per node */
  static std::vector<std::vector<int>> node_ids(comm_cart.size());
  for (auto &per_node : node_ids) {
    per_node.clear();
  }

  for (auto const &p_id : ids) {
    auto const p_node = get_particle_node(p_id);
    node_ids[p_node].emplace_back(p_id);
  }
  /* We shouldn't be prefetching particles that are already on the head node */
  assert(node_ids[this_node].empty());

  /* Distribute ids to the worker nodes */
  {
    static std::vector<int> ignore;
    boost::mpi::scatter(comm_cart, node_ids, ignore, 0);
    assert(ignore.empty());
  }

  static std::vector<int> node_sizes(comm_cart.size());
  // cannot use range-based transform with GCC 13 + ASAN
  std::transform(node_ids.begin(), node_ids.end(), node_sizes.begin(),
                 std::size<std::vector<int>>);

  Utils::Mpi::gatherv(comm_cart, parts.data(), static_cast<int>(parts.size()),
                      parts.data(), node_sizes.data(), 0);

  return parts;
}

void prefetch_particle_data(std::span<const int> in_ids) {
  /* Nothing to do on a single node. */
  // NOLINTNEXTLINE(clang-analyzer-core.NonNullParamChecker)
  if (comm_cart.size() == 1)
    return;

  static std::vector<int> ids;
  ids.clear();
  auto out_ids = std::back_inserter(ids);

  /* Don't prefetch particles already on the head node or already cached. */
  std::ranges::copy_if(in_ids, out_ids, [](int id) {
    return (get_particle_node(id) != this_node) && particle_fetch_cache.has(id);
  });

  /* Don't prefetch more particles than fit the cache. */
  if (ids.size() > particle_fetch_cache.max_size())
    ids.resize(particle_fetch_cache.max_size());

  /* Fetch the particles... */
  for (auto &p : mpi_get_particles(ids)) {
    auto id = p.id();
    auto const p_cached = particle_fetch_cache.put(id, std::move(p));
    // Attach the just-cached (detached) copy to the snapshot store; see the
    // matching call and const_cast rationale in get_particle_data.
    attach_cached_particle(const_cast<Particle &>(*p_cached));
  }
}

static void mpi_who_has_local() {
  static std::vector<int> sendbuf;

  auto local_particles = get_cell_structure().local_particles();
  auto const n_part = static_cast<int>(local_particles.size());
  boost::mpi::gather(comm_cart, n_part, 0);

  if (n_part == 0) {
    mpi_synchronize_max_seen_pid_local();
    return;
  }

  sendbuf.resize(n_part);

  std::transform(local_particles.begin(), local_particles.end(),
                 sendbuf.begin(), [](Particle const &p) { return p.id(); });

  comm_cart.send(0, some_tag, sendbuf);
  mpi_synchronize_max_seen_pid_local();
}

REGISTER_CALLBACK(mpi_who_has_local)

static void mpi_who_has_head() {
  auto local_particles = get_cell_structure().local_particles();

  static std::vector<int> n_parts;
  boost::mpi::gather(comm_cart, static_cast<int>(local_particles.size()),
                     n_parts, 0);

  static std::vector<int> pdata;
  auto const n_nodes = ::comm_cart.size();
  max_seen_pid = -1;

  /* then fetch particle locations */
  for (int pnode = 0; pnode < n_nodes; pnode++) {
    if (pnode == this_node) {
      for (auto const &p : local_particles) {
        particle_node[p.id()] = this_node;
        max_seen_pid = std::max(max_seen_pid, p.id());
      }
    } else if (n_parts[pnode] > 0) {
      pdata.resize(n_parts[pnode]);
      comm_cart.recv(pnode, some_tag, pdata);
      for (int i = 0; i < n_parts[pnode]; i++) {
        particle_node[pdata[i]] = pnode;
        max_seen_pid = std::max(max_seen_pid, pdata[i]);
      }
    }
  }
  mpi_synchronize_max_seen_pid_local();
}

/**
 * @brief Rebuild the particle index.
 */
static void build_particle_node() {
  Communication::mpiCallbacks().call(mpi_who_has_local);
  mpi_who_has_head();
}

/**
 * @brief Rebuild the particle index.
 */
static void build_particle_node_parallel() {
  if (this_node == 0) {
    mpi_who_has_head();
  } else {
    mpi_who_has_local();
  }
}

int get_particle_node(int p_id) {
  if (p_id < 0) {
    throw std::domain_error("Invalid particle id: " + std::to_string(p_id));
  }

  if (particle_node.empty())
    build_particle_node();

  auto const needle = particle_node.find(p_id);

  // Check if particle has a node, if not, we assume it does not exist.
  if (needle == particle_node.end()) {
    throw std::runtime_error("Particle node for id " + std::to_string(p_id) +
                             " not found!");
  }
  return needle->second;
}

int get_particle_node_parallel(int p_id) {
  if (p_id < 0) {
    throw std::domain_error("Invalid particle id: " + std::to_string(p_id));
  }

  if (rebuild_needed()) {
    build_particle_node_parallel();
  }

  if (this_node != 0) {
    return -1;
  }

  auto const needle = particle_node.find(p_id);

  // Check if particle has a node, if not, we assume it does not exist.
  if (needle == particle_node.end()) {
    throw std::runtime_error("Particle node for id " + std::to_string(p_id) +
                             " not found!");
  }
  return needle->second;
}

void clear_particle_node() {
  ::max_seen_pid = -1;
  particle_node.clear();
}

static void clear_particle_type_map() {
  for (auto &kv : ::particle_type_map) {
    kv.second.clear();
  }
}

/**
 * @brief Calculate the largest particle id.
 * Traversing the @ref particle_node to find the largest particle id
 * scales with O(N). Consider using the cached value in @ref max_seen_pid
 * if possible. This function is only necessary when the cached value is
 * invalidated, for example when removing the particle which has the
 * largest id.
 */
static int calculate_max_seen_id() {
  return std::accumulate(
      particle_node.begin(), particle_node.end(), -1,
      [](int max, auto const &kv) { return std::max(max, kv.first); });
}

/**
 * @brief Create a new particle and attach it to a cell.
 * @param p_id  The identity of the particle to create.
 * @param pos   The particle position.
 * @return Whether the particle was created on that node.
 */
static bool maybe_insert_particle(int p_id, Utils::Vector3d const &pos) {
  auto const &box_geo = *System::get_system().box_geo;
  auto folded_pos = pos;
  auto image_box = Utils::Vector3i{};
  box_geo.fold_position(folded_pos, image_box);

  Particle new_part;
  new_part.id() = p_id;
  new_part.pos() = folded_pos;
  new_part.image_box() = image_box;

  return get_cell_structure().add_local_particle(std::move(new_part)) !=
         nullptr;
}

/**
 * @brief Move particle to a new position.
 * @param p_id  The identity of the particle to move.
 * @param pos   The new particle position.
 * @return Whether the particle was moved from that node.
 */
static bool maybe_move_particle(int p_id, Utils::Vector3d const &pos) {
  auto const &system = System::get_system();
  auto const &box_geo = *system.box_geo;
  auto p = system.cell_structure->get_local_particle(p_id);
  if (p == nullptr) {
    return false;
  }
  auto folded_pos = pos;
  auto image_box = Utils::Vector3i{};
  box_geo.fold_position(folded_pos, image_box);
  p->pos() = folded_pos;
  p->image_box() = image_box;
  return true;
}

void remove_all_particles() {
  get_cell_structure().remove_all_particles();
  System::get_system().on_particle_change();
  clear_particle_node();
  clear_particle_type_map();
}

void remove_particle(int p_id) {
  if (::type_list_enable) {
    auto p = get_cell_structure().get_local_particle(p_id);
    auto p_type = -1;
    if (p != nullptr and not p->is_ghost()) {
      if (this_node == 0) {
        p_type = p->type();
      } else {
        ::comm_cart.send(0, 42, p->type());
      }
    } else if (this_node == 0) {
      ::comm_cart.recv(boost::mpi::any_source, 42, p_type);
    }
    assert(this_node != 0 or p_type != -1);
    boost::mpi::broadcast(::comm_cart, p_type, 0);
    remove_id_from_map(p_id, p_type);
  }

  if (this_node == 0) {
    particle_node[p_id] = -1;
  }
  get_cell_structure().remove_particle(p_id);
  System::get_system().on_particle_change();
  mpi_synchronize_max_seen_pid_local();
  if (this_node == 0) {
    particle_node.erase(p_id);
    if (p_id == ::max_seen_pid) {
      --::max_seen_pid;
      // if there is a gap (i.e. there is no particle with id max_seen_pid - 1,
      // then the cached value is invalidated and has to be recomputed (slow)
      if (not particle_node.contains(::max_seen_pid) or
          particle_node[::max_seen_pid] == -1) {
        ::max_seen_pid = calculate_max_seen_id();
      }
    }
  }
  mpi_synchronize_max_seen_pid_local();
}

void make_new_particle(int p_id, Utils::Vector3d const &pos) {
  if (rebuild_needed()) {
    build_particle_node_parallel();
  }
  auto const has_created = maybe_insert_particle(p_id, pos);
  System::get_system().on_particle_change();

  auto node = -1;
  auto const node_local = (has_created) ? ::comm_cart.rank() : 0;
  boost::mpi::reduce(::comm_cart, node_local, node, std::plus<int>{}, 0);
  if (::this_node == 0) {
    particle_node[p_id] = node;
    max_seen_pid = std::max(max_seen_pid, p_id);
    assert(not has_created or node == 0);
  }
  mpi_synchronize_max_seen_pid_local();
}

void set_particle_pos(int p_id, Utils::Vector3d const &pos) {
  auto const has_moved = maybe_move_particle(p_id, pos);
  get_cell_structure().set_resort_particles(Cells::RESORT_GLOBAL);
  System::get_system().on_particle_change();

  auto success = false;
  boost::mpi::reduce(::comm_cart, has_moved, success, std::plus<bool>{}, 0);
  if (::this_node == 0 and !success) {
    throw std::runtime_error("Particle node for id " + std::to_string(p_id) +
                             " not found!");
  }
}

int get_random_p_id(int type, int random_index_in_type_map) {
  auto it = particle_type_map.find(type);
  if (it == particle_type_map.end()) {
    throw std::runtime_error("The provided particle type " +
                             std::to_string(type) +
                             " is currently not tracked by the system.");
  }

  if (random_index_in_type_map + 1 > it->second.size())
    throw std::runtime_error("The provided index exceeds the number of "
                             "particle types listed in the particle_type_map");
  // there is no guarantee of order across MPI ranks
  auto p_id = *std::next(it->second.begin(), random_index_in_type_map);
  boost::mpi::broadcast(::comm_cart, p_id, 0);
  return p_id;
}

int number_of_particles_with_type(int type) {
  auto it = particle_type_map.find(type);
  if (it == particle_type_map.end()) {
    throw std::runtime_error("The provided particle type " +
                             std::to_string(type) +
                             " is currently not tracked by the system.");
  }

  return static_cast<int>(it->second.size());
}

bool particle_exists(int p_id) {
  if (particle_node.empty())
    build_particle_node();
  return particle_node.contains(p_id);
}

std::vector<int> get_particle_ids() {
  if (particle_node.empty())
    build_particle_node();

  std::vector<int> pids{};
  std::ranges::copy(std::views::keys(particle_node), std::back_inserter(pids));
  std::ranges::sort(pids);

  return pids;
}

std::vector<int> get_particle_ids_parallel() {
  if (rebuild_needed()) {
    build_particle_node_parallel();
  }
  std::vector<int> pids{};
  std::ranges::copy(std::views::keys(particle_node), std::back_inserter(pids));
  boost::mpi::broadcast(::comm_cart, pids, 0);
  return pids;
}

int get_maximal_particle_id() {
  if (rebuild_needed()) {
    build_particle_node_parallel();
  }

  return max_seen_pid;
}

int get_n_part() {
  if (particle_node.empty())
    build_particle_node();

  return static_cast<int>(particle_node.size());
}
