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
   *  force (and torque) plus the position and momentum groups.
   *
   *  Migration phases 3 & 4: assign_row seeds the fresh store's row from @p p's
   *  MIGRATION CARRIERS, which are only refreshed on a serialize SAVE. A
   *  particle handed to us via the head-node fast path (@c received = *p, a
   *  plain copy, no serialization) therefore carries STALE position/momentum
   *  carriers relative to its (correct) live column. Capture the current
   *  position-group and velocity/angular-velocity values through the accessors
   *  first (they read the column when attached, or the carrier when detached --
   *  both correct here), then write them back through the fresh store's proxies
   *  after attach so the copy reflects the source's actual state, not the stale
   *  carrier. */
  void assign(Particle const &p, Utils::Vector3d const &force
#ifdef ESPRESSO_ROTATION
              ,
              Utils::Vector3d const &torque
#endif
  ) {
    // Capture the source's current position/momentum groups through the const
    // accessors (column when attached, carrier when detached -- both correct)
    // BEFORE the copy is re-attached and re-seeded from its (possibly stale)
    // carrier.
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
    auto const velocity = p.v();
#ifdef ESPRESSO_ROTATION
    auto const angular_velocity = p.omega();
#endif
    // Migration phase 5: ALL parameter fields also live in the store now, so
    // the same stale-carrier problem applies to the head-node fast-path copy.
    // Capture them through the accessors before the re-attach + re-seed.
    auto const id = p.id();
    auto const mol_id = p.mol_id();
    auto const type = p.type();
    auto const propagation = p.propagation();
#ifdef ESPRESSO_ROTATION
    auto const rotation = p.rotation();
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    auto const ext_flag = p.fixed();
    auto const ext_force = Utils::Vector3d(p.ext_force());
#ifdef ESPRESSO_ROTATION
    auto const ext_torque = Utils::Vector3d(p.ext_torque());
#endif
#endif
#ifdef ESPRESSO_MASS
    auto const mass = p.mass();
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    auto const q = p.q();
#endif
#ifdef ESPRESSO_DIPOLES
    auto const dipm = p.dipm();
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    auto const rinertia = Utils::Vector3d(p.rinertia());
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
    auto const mu_E = Utils::Vector3d(p.mu_E());
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    auto const dip_fld = Utils::Vector3d(p.dip_fld());
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    auto const gamma = p.gamma();
#ifdef ESPRESSO_ROTATION
    auto const gamma_rot = p.gamma_rot();
#endif
#endif
#ifdef ESPRESSO_ENGINE
    auto const swimming = p.swimming();
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
    auto const magnetodynamics = p.magnetodynamics();
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    auto const vs_relative = p.vs_relative();
#endif

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
    m_particle->v() = velocity;
#ifdef ESPRESSO_ROTATION
    m_particle->omega() = angular_velocity;
#endif
    m_particle->force() = force;
#ifdef ESPRESSO_ROTATION
    m_particle->torque() = torque;
#endif
    m_particle->id() = id;
    m_particle->mol_id() = mol_id;
    m_particle->type() = type;
    m_particle->propagation() = propagation;
#ifdef ESPRESSO_ROTATION
    m_particle->rotation() = rotation;
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    m_particle->fixed() = ext_flag;
    m_particle->ext_force() = ext_force;
#ifdef ESPRESSO_ROTATION
    m_particle->ext_torque() = ext_torque;
#endif
#endif
#ifdef ESPRESSO_MASS
    m_particle->mass() = mass;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    m_particle->q() = q;
#endif
#ifdef ESPRESSO_DIPOLES
    m_particle->dipm() = dipm;
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    m_particle->rinertia() = rinertia;
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
    m_particle->mu_E() = mu_E;
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    m_particle->dip_fld() = dip_fld;
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    m_particle->gamma() = gamma;
#ifdef ESPRESSO_ROTATION
    m_particle->gamma_rot() = gamma_rot;
#endif
#endif
#ifdef ESPRESSO_ENGINE
    m_particle->swimming() = swimming;
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
    m_particle->magnetodynamics() = magnetodynamics;
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    m_particle->vs_relative() = vs_relative;
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
