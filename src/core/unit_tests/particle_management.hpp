/*
 * Copyright (C) 2021-2026 The ESPResSo project
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

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "particle_store/ParticleStore.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <optional>
#include <utility>

/**
 * @brief An optional particle that owns the standalone ParticleStore backing
 * its force/torque columns.
 *
 * Migration phase 2: force/torque live in the ParticleStore and are not part
 * of the Particle struct. A particle copied to the head node needs a store so
 * its accessors work; bundling the store here ties its lifetime to the copy,
 * which matters because the store holds Kokkos Views that must be destroyed
 * before Kokkos::finalize(). Mimics the subset of the std::optional<Particle>
 * interface used by the tests (bool, has_value, operator*).
 */
class ParticleWithStore {
  std::optional<Particle> m_particle{};
  ParticleStore m_store{};

  void reattach() {
    if (m_particle) {
      // Re-point the particle at the moved store WITHOUT re-seeding: the values
      // already live in row 0 of the moved-from store's Views. Calling
      // assign_row here would re-seed from the stale migration carrier, so use
      // attach_to_store to only fix up the pointer/row.
      m_particle->attach_to_store(m_store, 0);
    }
  }

public:
  ParticleWithStore() = default;
  ParticleWithStore(ParticleWithStore &&other) noexcept
      : m_particle(std::move(other.m_particle)),
        m_store(std::move(other.m_store)) {
    // the moved store keeps its Kokkos Views; re-point the particle at it
    reattach();
  }
  ParticleWithStore &operator=(ParticleWithStore &&other) noexcept {
    m_particle = std::move(other.m_particle);
    m_store = std::move(other.m_store);
    reattach();
    return *this;
  }
  ParticleWithStore(ParticleWithStore const &) = delete;
  ParticleWithStore &operator=(ParticleWithStore const &) = delete;

  /** Take ownership of @p p, giving it a single-row store, and restore its
   *  force (and torque) plus the position group.
   *
   *  Migration phase 3: assign_row seeds the fresh store's row from @p p's
   *  MIGRATION CARRIERS, which are only refreshed on a serialize SAVE. A
   *  particle handed to us via the head-node fast path (@c received = *p, a
   *  plain copy, no serialization) therefore carries STALE position-group
   *  carriers relative to its (correct) live column. Capture the current
   *  position-group values through the accessors first (they read the column
   *  when attached, or the carrier when detached -- both correct here), then
   *  write them back through the fresh store's proxies after attach so the copy
   *  reflects the source's actual state, not the stale carrier. */
  void assign(Particle const &p, Utils::Vector3d const &force
#ifdef ESPRESSO_ROTATION
              ,
              Utils::Vector3d const &torque
#endif
  ) {
    // Capture the source's current position group through the const accessors
    // (column when attached, carrier when detached -- both correct) BEFORE the
    // copy is re-attached and re-seeded from its (possibly stale) carrier.
    auto const pos = p.pos();
    auto const image_box = p.image_box();
    auto const pos_at_last_verlet_update = p.pos_at_last_verlet_update();
#ifdef ESPRESSO_ROTATION
    auto const quat = p.quat();
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    auto const pos_last_time_step = p.pos_last_time_step();
#endif
    auto const lees_edwards_offset = p.lees_edwards_offset();
    auto const lees_edwards_flag = p.lees_edwards_flag();

    m_particle = p;
    m_store.begin_rebuild(1u, 0u);
    m_store.finish_rebuild();
    m_store.assign_row(*m_particle, 0);

    m_particle->pos() = pos;
    m_particle->image_box() = image_box;
    m_particle->pos_at_last_verlet_update() = pos_at_last_verlet_update;
#ifdef ESPRESSO_ROTATION
    m_particle->quat() = quat;
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    m_particle->pos_last_time_step() = pos_last_time_step;
#endif
    m_particle->lees_edwards_offset() = lees_edwards_offset;
    m_particle->lees_edwards_flag() = lees_edwards_flag;
    m_particle->force() = force;
#ifdef ESPRESSO_ROTATION
    m_particle->torque() = torque;
#endif
  }

  explicit operator bool() const { return m_particle.has_value(); }
  bool has_value() const { return m_particle.has_value(); }
  Particle &operator*() { return *m_particle; }
  Particle const &operator*() const { return *m_particle; }
  Particle *operator->() { return &*m_particle; }
  Particle const *operator->() const { return &*m_particle; }
};

/**
 * @brief Copy a particle (including its force/torque) to the head node.
 *
 * Migration phase 2: force/torque live in the ParticleStore columns and are
 * neither part of the Particle struct nor serialized. To keep them observable
 * on the head node, they are transmitted explicitly and the returned copy owns
 * a standalone store so its accessors work.
 */
inline ParticleWithStore
copy_particle_to_head_node(boost::mpi::communicator const &comm,
                           System::System &system, int p_id) {
  std::optional<Particle> received{};
  Utils::Vector3d force{};
#ifdef ESPRESSO_ROTATION
  Utils::Vector3d torque{};
#endif
  system.cell_structure->ensure_particle_store_synchronized();
  auto p = system.cell_structure->get_local_particle(p_id);
  if (p and not p->is_ghost()) {
    force = Utils::Vector3d(p->force());
#ifdef ESPRESSO_ROTATION
    torque = Utils::Vector3d(p->torque());
#endif
    if (comm.rank() == 0) {
      received = *p;
    } else {
      comm.send(0, p_id, *p);
      comm.send(0, p_id, force);
#ifdef ESPRESSO_ROTATION
      comm.send(0, p_id, torque);
#endif
    }
  }
  if (comm.rank() == 0 and not received) {
    Particle q{};
    comm.recv(boost::mpi::any_source, p_id, q);
    comm.recv(boost::mpi::any_source, p_id, force);
#ifdef ESPRESSO_ROTATION
    comm.recv(boost::mpi::any_source, p_id, torque);
#endif
    received = q;
  }
  ParticleWithStore result{};
  if (comm.rank() == 0 and received) {
    result.assign(*received, force
#ifdef ESPRESSO_ROTATION
                  ,
                  torque
#endif
    );
  }
  return result;
}
