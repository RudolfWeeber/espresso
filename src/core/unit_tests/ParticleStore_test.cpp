/*
 * Copyright (C) 2017-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE ParticleStore test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "particle_store/ParticleStore.hpp"

#include "Particle.hpp"

#include <utils/Vector.hpp>
#include <utils/quaternion.hpp>

#include <Kokkos_Core.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <cstddef>
#include <sstream>
#include <vector>

// ParticleStore allocates Kokkos Views, which requires an initialized runtime.
struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

BOOST_AUTO_TEST_CASE(default_constructed_store_is_empty) {
  ParticleStore const store{};
  BOOST_CHECK_EQUAL(store.number_of_local_particles(), 0ul);
  BOOST_CHECK_EQUAL(store.number_of_ghost_particles(), 0ul);
}

BOOST_AUTO_TEST_CASE(rebuild_assigns_rows_and_zero_initializes) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.mark_dirty();
  BOOST_CHECK(store.is_dirty());
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();
  BOOST_CHECK(not store.is_dirty());
  BOOST_CHECK_EQUAL(store.number_of_local_particles(), 2ul);
  BOOST_CHECK_EQUAL(p0.store_row(), 0);
  BOOST_CHECK_EQUAL(p1.store_row(), 1);
  Utils::Vector3d const f0 = store.force_value(0);
  BOOST_CHECK_EQUAL(f0.norm2(), 0.);
}

BOOST_AUTO_TEST_CASE(rebuild_preserves_values_by_old_row) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();
  store.force_reference(p0.store_row()) = Utils::Vector3d{1., 2., 3.};
  store.force_reference(p1.store_row()) = Utils::Vector3d{4., 5., 6.};

  // simulate a resort that swaps the two particles' order and adds one
  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

  Utils::Vector3d const f_p1 = store.force_value(p1.store_row());
  Utils::Vector3d const f_p0 = store.force_value(p0.store_row());
  Utils::Vector3d const f_p2 = store.force_value(p2.store_row());
  BOOST_CHECK_EQUAL(f_p1[0], 4.);      // p1 kept its values at its new row
  BOOST_CHECK_EQUAL(f_p0[2], 3.);      // p0 kept its values at its new row
  BOOST_CHECK_EQUAL(f_p2.norm2(), 0.); // new particle zero-initialized
}

BOOST_AUTO_TEST_CASE(ghost_rows_follow_locals) {
  ParticleStore store{};
  Particle p_local{}, p_ghost{};
  store.begin_rebuild(1u, 1u);
  store.assign_row(p_local, 0);
  store.assign_row(p_ghost, 1);
  store.finish_rebuild();
  BOOST_CHECK_EQUAL(store.number_of_local_particles(), 1ul);
  BOOST_CHECK_EQUAL(store.number_of_ghost_particles(), 1ul);
}

// Regression guard (phase-2 hybrid/n_square 4-rank force bug): force/torque
// live only in the store columns, which are not carried by the boost
// serialization used for cross-rank particle migration. A particle that
// migrates to another rank arrives detached from this store; its force must be
// ferried through the Particle serialization carrier and seeded into the new
// row by assign_row. Otherwise a global resort followed by a force-reusing
// integrator step would report a zeroed force (the original bug).
BOOST_AUTO_TEST_CASE(rebuild_seeds_migrated_particle_force_from_carrier) {
  // 1) A particle attached to a source store with a known force (models the
  //    sending rank right before a global resort).
  ParticleStore source{};
  Particle p{};
  p.id() = 7;
  source.begin_rebuild(1u, 0u);
  source.assign_row(p, 0);
  source.finish_rebuild();
  auto const f_ref = Utils::Vector3d{-1.5, 2.25, -3.75};
  source.force_reference(p.store_row()) = f_ref;
#ifdef ESPRESSO_ROTATION
  auto const t_ref = Utils::Vector3d{4.5, -5.5, 6.5};
  source.torque_reference(p.store_row()) = t_ref;
#endif

  // 2) Serialize and deserialize the particle, exactly as the cross-rank
  //    migration does. The receiving particle is detached (no store).
  std::stringstream stream;
  {
    boost::archive::text_oarchive oa{stream};
    oa << p;
  }
  Particle received{};
  {
    boost::archive::text_iarchive ia{stream};
    ia >> received;
  }
  BOOST_REQUIRE(received.store() == nullptr);
  BOOST_CHECK_EQUAL(received.migration_force()[0], f_ref[0]);
  BOOST_CHECK_EQUAL(received.migration_force()[2], f_ref[2]);

  // 3) The receiving rank's store rebuild must seed the new row from the
  //    carrier so the force survives the migration.
  ParticleStore target{};
  target.begin_rebuild(1u, 0u);
  target.assign_row(received, 0);
  target.finish_rebuild();
  auto const f_new = target.force_value(received.store_row());
  BOOST_CHECK_EQUAL(f_new[0], f_ref[0]);
  BOOST_CHECK_EQUAL(f_new[1], f_ref[1]);
  BOOST_CHECK_EQUAL(f_new[2], f_ref[2]);
#ifdef ESPRESSO_ROTATION
  auto const t_new = target.torque_value(received.store_row());
  BOOST_CHECK_EQUAL(t_new[0], t_ref[0]);
  BOOST_CHECK_EQUAL(t_new[2], t_ref[2]);
#endif
}

// Phase-3 state columns: a rank-local rebuild that shuffles the row order must
// preserve each particle's state (position, image box, quaternion, and the
// Lees-Edwards offset) by old row, exactly as the force column does.
// Phase-4 velocity columns are verified alongside.
BOOST_AUTO_TEST_CASE(rebuild_preserves_state_columns_by_old_row) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();

  store.position_reference(p0.store_row()) = Utils::Vector3d{1., 2., 3.};
  store.position_reference(p1.store_row()) = Utils::Vector3d{4., 5., 6.};
  store.image_box_reference(p0.store_row()) = Utils::Vector3i{-1, 0, 2};
  store.image_box_reference(p1.store_row()) = Utils::Vector3i{7, -8, 9};
  store.position_at_last_verlet_update_reference(p0.store_row()) =
      Utils::Vector3d{0.1, 0.2, 0.3};
  store.position_at_last_verlet_update_reference(p1.store_row()) =
      Utils::Vector3d{0.4, 0.5, 0.6};
  store.lees_edwards_offset(p0.store_row()) = 1.25;
  store.lees_edwards_offset(p1.store_row()) = -2.5;
  store.lees_edwards_flag(p0.store_row()) = static_cast<short>(3);
  store.lees_edwards_flag(p1.store_row()) = static_cast<short>(-4);
  store.velocity_reference(p0.store_row()) = Utils::Vector3d{1.1, 2.2, 3.3};
  store.velocity_reference(p1.store_row()) = Utils::Vector3d{-4.4, 5.5, -6.6};
#ifdef ESPRESSO_ROTATION
  auto const q0 = Utils::Quaternion<double>{{0.5, 0.5, 0.5, 0.5}};
  auto const q1 = Utils::Quaternion<double>{{0., 1., 0., 0.}};
  store.quaternion_reference(p0.store_row()) = q0;
  store.quaternion_reference(p1.store_row()) = q1;
  store.angular_velocity_reference(p0.store_row()) =
      Utils::Vector3d{0.1, 0.2, 0.3};
  store.angular_velocity_reference(p1.store_row()) =
      Utils::Vector3d{-0.4, 0.5, -0.6};
#endif

  // resort that swaps the two particles' order and appends a new one
  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

  Utils::Vector3d const pos_p1 = store.position_value(p1.store_row());
  Utils::Vector3d const pos_p0 = store.position_value(p0.store_row());
  BOOST_CHECK_EQUAL(pos_p1[0], 4.);
  BOOST_CHECK_EQUAL(pos_p0[2], 3.);
  Utils::Vector3i const img_p1 = store.image_box_value(p1.store_row());
  Utils::Vector3i const img_p0 = store.image_box_value(p0.store_row());
  BOOST_CHECK_EQUAL(img_p1[0], 7);
  BOOST_CHECK_EQUAL(img_p0[2], 2);
  Utils::Vector3d const pold_p1 =
      store.position_at_last_verlet_update_value(p1.store_row());
  BOOST_CHECK_EQUAL(pold_p1[1], 0.5);
  BOOST_CHECK_EQUAL(store.lees_edwards_offset(p1.store_row()), -2.5);
  BOOST_CHECK_EQUAL(store.lees_edwards_offset(p0.store_row()), 1.25);
  BOOST_CHECK_EQUAL(store.lees_edwards_flag(p1.store_row()),
                    static_cast<short>(-4));
  BOOST_CHECK_EQUAL(store.lees_edwards_flag(p0.store_row()),
                    static_cast<short>(3));
  // Phase-4: velocity preserved by old row; new row zeroed.
  Utils::Vector3d const vel_p1 = store.velocity_value(p1.store_row());
  Utils::Vector3d const vel_p0 = store.velocity_value(p0.store_row());
  Utils::Vector3d const vel_p2 = store.velocity_value(p2.store_row());
  BOOST_CHECK_EQUAL(vel_p1[0], -4.4);
  BOOST_CHECK_EQUAL(vel_p0[2], 3.3);
  BOOST_CHECK_EQUAL(vel_p2.norm2(), 0.);
#ifdef ESPRESSO_ROTATION
  Utils::Quaternion<double> const quat_p1 =
      store.quaternion_value(p1.store_row());
  Utils::Quaternion<double> const quat_p0 =
      store.quaternion_value(p0.store_row());
  BOOST_CHECK_EQUAL(quat_p1[1], 1.);
  BOOST_CHECK_EQUAL(quat_p0[3], 0.5);
  // Phase-4: angular velocity preserved by old row; new row zeroed.
  Utils::Vector3d const av_p1 = store.angular_velocity_value(p1.store_row());
  Utils::Vector3d const av_p0 = store.angular_velocity_value(p0.store_row());
  Utils::Vector3d const av_p2 = store.angular_velocity_value(p2.store_row());
  BOOST_CHECK_EQUAL(av_p1[1], 0.5);
  BOOST_CHECK_EQUAL(av_p0[2], 0.3);
  BOOST_CHECK_EQUAL(av_p2.norm2(), 0.);
#endif
}

// A genuinely new row (e.g. a fresh ghost, before its field update overwrites
// the column) must default to sensible values that match Particle's member
// defaults. The quaternion in particular must be the IDENTITY (1,0,0,0), NOT
// the Kokkos zero-init (0,0,0,0) which is an invalid quaternion.
// Phase-4: velocity and angular velocity must default to zero.
BOOST_AUTO_TEST_CASE(new_row_state_defaults) {
  ParticleStore store{};
  Particle p{};
  store.begin_rebuild(1u, 0u);
  store.assign_row(p, 0);
  store.finish_rebuild();

  Utils::Vector3d const pos = store.position_value(p.store_row());
  BOOST_CHECK_EQUAL(pos.norm2(), 0.);
  Utils::Vector3i const img = store.image_box_value(p.store_row());
  BOOST_CHECK_EQUAL(img[0], 0);
  BOOST_CHECK_EQUAL(img[1], 0);
  BOOST_CHECK_EQUAL(img[2], 0);
  BOOST_CHECK_EQUAL(store.lees_edwards_offset(p.store_row()), 0.);
  BOOST_CHECK_EQUAL(store.lees_edwards_flag(p.store_row()),
                    static_cast<short>(0));
  // Phase-4: velocity defaults to zero (from migration carrier default
  // {0,0,0}).
  Utils::Vector3d const vel = store.velocity_value(p.store_row());
  BOOST_CHECK_EQUAL(vel.norm2(), 0.);
#ifdef ESPRESSO_ROTATION
  Utils::Quaternion<double> const quat = store.quaternion_value(p.store_row());
  BOOST_CHECK_EQUAL(quat[0], 1.); // identity, not zero
  BOOST_CHECK_EQUAL(quat[1], 0.);
  BOOST_CHECK_EQUAL(quat[2], 0.);
  BOOST_CHECK_EQUAL(quat[3], 0.);
  // Phase-4: angular velocity defaults to zero.
  Utils::Vector3d const av = store.angular_velocity_value(p.store_row());
  BOOST_CHECK_EQUAL(av.norm2(), 0.);
#endif
}

// Phase-4: velocity columns are seeded from the migration carrier on
// assign_row. Pre-flip the carrier is always zero (serialization not wired
// yet), so a detached particle produces zero velocity columns for new rows.
// This guards that assign_row seeds the velocity column (not skip it) and that
// the default is exactly zero (not garbage from WithoutInitializing).
BOOST_AUTO_TEST_CASE(rebuild_seeds_velocity_from_carrier) {
  // detached particle — migration carriers at their default {0,0,0}
  Particle p{};
  BOOST_REQUIRE(p.store() == nullptr);
  BOOST_CHECK_EQUAL(p.migration_velocity()[0], 0.);
  BOOST_CHECK_EQUAL(p.migration_velocity()[2], 0.);

  ParticleStore store{};
  store.begin_rebuild(1u, 0u);
  store.assign_row(p, 0);
  store.finish_rebuild();

  Utils::Vector3d const vel = store.velocity_value(p.store_row());
  BOOST_CHECK_EQUAL(vel[0], 0.);
  BOOST_CHECK_EQUAL(vel[1], 0.);
  BOOST_CHECK_EQUAL(vel[2], 0.);
#ifdef ESPRESSO_ROTATION
  BOOST_CHECK_EQUAL(p.migration_angular_velocity()[0], 0.);
  Utils::Vector3d const av = store.angular_velocity_value(p.store_row());
  BOOST_CHECK_EQUAL(av[0], 0.);
  BOOST_CHECK_EQUAL(av[1], 0.);
  BOOST_CHECK_EQUAL(av[2], 0.);
#endif
}

// A detached particle with its migration carriers populated must have all its
// state columns seeded from those carriers on assign_row (models a particle
// that just migrated in from another rank). Post-flip (phase 3) the state
// fields live only in the store columns, so a detached particle's carriers are
// populated the way a real migration does: set the values on an attached
// source particle, serialize it out (SAVE fills the carriers from the columns),
// and deserialize into a detached receiver (LOAD lands them in its carriers).
BOOST_AUTO_TEST_CASE(rebuild_seeds_migrated_particle_state_from_carrier) {
  ParticleStore source{};
  Particle src{};
  src.id() = 11;
  source.begin_rebuild(1u, 0u);
  source.assign_row(src, 0);
  source.finish_rebuild();
  source.position_reference(src.store_row()) = Utils::Vector3d{1.5, -2.5, 3.5};
  source.image_box_reference(src.store_row()) = Utils::Vector3i{4, -5, 6};
  source.position_at_last_verlet_update_reference(src.store_row()) =
      Utils::Vector3d{7.5, 8.5, 9.5};
  source.lees_edwards_offset(src.store_row()) = 12.75;
  source.lees_edwards_flag(src.store_row()) = static_cast<short>(2);
#ifdef ESPRESSO_ROTATION
  source.quaternion_reference(src.store_row()) =
      Utils::Quaternion<double>{{0., 0., 1., 0.}};
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  source.position_last_time_step_reference(src.store_row()) =
      Utils::Vector3d{-1., -2., -3.};
#endif

  std::stringstream stream;
  {
    boost::archive::text_oarchive oa{stream};
    oa << src;
  }
  Particle p{};
  {
    boost::archive::text_iarchive ia{stream};
    ia >> p;
  }
  BOOST_REQUIRE(p.store() == nullptr);
  BOOST_CHECK_EQUAL(p.migration_position()[0], 1.5);
  BOOST_CHECK_EQUAL(p.migration_image_box()[2], 6);

  ParticleStore target{};
  target.begin_rebuild(1u, 0u);
  target.assign_row(p, 0);
  target.finish_rebuild();

  Utils::Vector3d const pos = target.position_value(p.store_row());
  BOOST_CHECK_EQUAL(pos[0], 1.5);
  BOOST_CHECK_EQUAL(pos[1], -2.5);
  BOOST_CHECK_EQUAL(pos[2], 3.5);
  Utils::Vector3i const img = target.image_box_value(p.store_row());
  BOOST_CHECK_EQUAL(img[0], 4);
  BOOST_CHECK_EQUAL(img[2], 6);
  Utils::Vector3d const pold =
      target.position_at_last_verlet_update_value(p.store_row());
  BOOST_CHECK_EQUAL(pold[1], 8.5);
  BOOST_CHECK_EQUAL(target.lees_edwards_offset(p.store_row()), 12.75);
  BOOST_CHECK_EQUAL(target.lees_edwards_flag(p.store_row()),
                    static_cast<short>(2));
#ifdef ESPRESSO_ROTATION
  Utils::Quaternion<double> const quat = target.quaternion_value(p.store_row());
  BOOST_CHECK_EQUAL(quat[2], 1.);
  BOOST_CHECK_EQUAL(quat[0], 0.);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  Utils::Vector3d const plast =
      target.position_last_time_step_value(p.store_row());
  BOOST_CHECK_EQUAL(plast[0], -1.);
  BOOST_CHECK_EQUAL(plast[2], -3.);
#endif
  // Phase-4: velocity carrier is not serialized yet (pre-flip); column seeds
  // to zero from the default carrier.
  Utils::Vector3d const vel = target.velocity_value(p.store_row());
  BOOST_CHECK_EQUAL(vel.norm2(), 0.);
#ifdef ESPRESSO_ROTATION
  Utils::Vector3d const av = target.angular_velocity_value(p.store_row());
  BOOST_CHECK_EQUAL(av.norm2(), 0.);
#endif
}

// Head-node fetch-cache SNAPSHOT-STORE pattern (migration phase 3): detached
// particles (as they arrive from a worker rank, carriers holding their values)
// are attached to a FIXED-capacity store built once per invalidation epoch,
// with rows handed out monotonically. This mirrors attach_cached_particle in
// particle_node.cpp without needing MPI: build a small store once, attach a
// batch of detached particles to consecutive rows, and check each keeps its own
// carrier-seeded state (never another particle's row). The head-node cache
// itself is exercised end-to-end by the multi-rank python gates.
BOOST_AUTO_TEST_CASE(snapshot_store_attaches_batch_of_detached_particles) {
  constexpr std::size_t capacity = 4u;
  ParticleStore store{};
  store.begin_rebuild(capacity, 0u);
  store.finish_rebuild();

  // Each detached particle's carriers are populated the way a real fetch does:
  // an attached source particle is serialized out (SAVE fills the carriers from
  // its columns) and deserialized into a detached receiver (LOAD lands them in
  // its carriers).
  auto const make_detached = [](int id, Utils::Vector3d const &pos,
                                Utils::Vector3i const &image_box) {
    ParticleStore src_store{};
    Particle src{};
    src.id() = id;
    src_store.begin_rebuild(1u, 0u);
    src_store.assign_row(src, 0);
    src_store.finish_rebuild();
    src_store.position_reference(src.store_row()) = pos;
    src_store.image_box_reference(src.store_row()) = image_box;
    std::stringstream stream;
    {
      boost::archive::text_oarchive oa{stream};
      oa << src;
    }
    Particle received{};
    {
      boost::archive::text_iarchive ia{stream};
      ia >> received;
    }
    return received;
  };

  std::vector<Particle> parts(capacity);
  int next_row = 0;
  for (std::size_t i = 0u; i < capacity; ++i) {
    parts[i] = make_detached(
        static_cast<int>(i),
        Utils::Vector3d{double(i), double(i) + 0.5, double(i) + 0.25},
        Utils::Vector3i{int(i), -int(i), 2 * int(i)});
    auto &p = parts[i];
    BOOST_REQUIRE(p.store() ==
                  nullptr);          // detached, like a freshly-fetched copy
    store.assign_row(p, next_row++); // seeds the row from the carriers
  }

  for (std::size_t i = 0u; i < capacity; ++i) {
    auto const &p = parts[i];
    BOOST_CHECK_EQUAL(p.store_row(), static_cast<int>(i));
    auto const pos = store.position_value(p.store_row());
    BOOST_CHECK_EQUAL(pos[0], double(i));
    BOOST_CHECK_EQUAL(pos[1], double(i) + 0.5);
    auto const img = store.image_box_value(p.store_row());
    BOOST_CHECK_EQUAL(img[0], int(i));
    BOOST_CHECK_EQUAL(img[2], 2 * int(i));
    // Reading through the (now-attached) particle accessor matches the column.
    BOOST_CHECK_EQUAL(p.pos()[0], double(i));
  }
}

// Phase-5 PARAMETER columns/sidecars: a rank-local rebuild that shuffles the
// row order must preserve each particle's parameters (representative subset:
// id, type, mass, gamma, ext_force, dip_fld) by old row, and seed a genuinely
// new row from the migration carrier defaults, exactly as the state columns do.
BOOST_AUTO_TEST_CASE(rebuild_preserves_parameter_columns_by_old_row) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();

  store.id(p0.store_row()) = 10;
  store.id(p1.store_row()) = 20;
  store.type(p0.store_row()) = 3;
  store.type(p1.store_row()) = 7;
#ifdef ESPRESSO_MASS
  store.mass(p0.store_row()) = 2.5;
  store.mass(p1.store_row()) = 4.25;
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  store.gamma_reference(p0.store_row()) = Utils::Vector3d{1., 2., 3.};
  store.gamma_reference(p1.store_row()) = Utils::Vector3d{4., 5., 6.};
#else
  store.gamma_reference(p0.store_row()) = 1.5;
  store.gamma_reference(p1.store_row()) = 2.5;
#endif
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  store.ext_force_reference(p0.store_row()) = Utils::Vector3d{-1., -2., -3.};
  store.ext_force_reference(p1.store_row()) = Utils::Vector3d{7., 8., 9.};
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  store.dip_fld_reference(p0.store_row()) = Utils::Vector3d{0.1, 0.2, 0.3};
  store.dip_fld_reference(p1.store_row()) = Utils::Vector3d{0.4, 0.5, 0.6};
#endif

  // resort that swaps the two particles' order and appends a new one
  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

  BOOST_CHECK_EQUAL(store.id(p1.store_row()), 20);
  BOOST_CHECK_EQUAL(store.id(p0.store_row()), 10);
  BOOST_CHECK_EQUAL(store.id(p2.store_row()), -1); // default
  BOOST_CHECK_EQUAL(store.type(p1.store_row()), 7);
  BOOST_CHECK_EQUAL(store.type(p0.store_row()), 3);
  BOOST_CHECK_EQUAL(store.type(p2.store_row()), 0); // default
#ifdef ESPRESSO_MASS
  BOOST_CHECK_EQUAL(store.mass(p1.store_row()), 4.25);
  BOOST_CHECK_EQUAL(store.mass(p0.store_row()), 2.5);
  BOOST_CHECK_EQUAL(store.mass(p2.store_row()), 1.0); // default
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  Utils::Vector3d const g_p1 = store.gamma_value(p1.store_row());
  Utils::Vector3d const g_p0 = store.gamma_value(p0.store_row());
  Utils::Vector3d const g_p2 = store.gamma_value(p2.store_row());
  BOOST_CHECK_EQUAL(g_p1[0], 4.);
  BOOST_CHECK_EQUAL(g_p0[2], 3.);
  BOOST_CHECK_EQUAL(g_p2[0], -1.); // default {-1,-1,-1}
  BOOST_CHECK_EQUAL(g_p2[2], -1.);
#else
  BOOST_CHECK_EQUAL(store.gamma_value(p1.store_row()), 2.5);
  BOOST_CHECK_EQUAL(store.gamma_value(p0.store_row()), 1.5);
  BOOST_CHECK_EQUAL(store.gamma_value(p2.store_row()), -1.); // default
#endif
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  Utils::Vector3d const ef_p1 = store.ext_force_value(p1.store_row());
  Utils::Vector3d const ef_p0 = store.ext_force_value(p0.store_row());
  Utils::Vector3d const ef_p2 = store.ext_force_value(p2.store_row());
  BOOST_CHECK_EQUAL(ef_p1[0], 7.);
  BOOST_CHECK_EQUAL(ef_p0[2], -3.);
  BOOST_CHECK_EQUAL(ef_p2.norm2(), 0.); // default {0,0,0}
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  Utils::Vector3d const df_p1 = store.dip_fld_value(p1.store_row());
  Utils::Vector3d const df_p0 = store.dip_fld_value(p0.store_row());
  Utils::Vector3d const df_p2 = store.dip_fld_value(p2.store_row());
  BOOST_CHECK_EQUAL(df_p1[0], 0.4);
  BOOST_CHECK_EQUAL(df_p0[2], 0.3);
  BOOST_CHECK_EQUAL(df_p2.norm2(), 0.); // default {0,0,0}
#endif
}

// Phase-5 uint8 column (rotation / ext_flag): exercises the NEW element type.
// Bitfield values must round-trip through the DualView<uint8_t*> column, write
// through the element reference, and preserve by old row across a rebuild.
#if defined(ESPRESSO_ROTATION) || defined(ESPRESSO_EXTERNAL_FORCES)
BOOST_AUTO_TEST_CASE(uint8_parameter_column_write_through_and_preserve) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();

#ifdef ESPRESSO_ROTATION
  store.rotation(p0.store_row()) = static_cast<std::uint8_t>(0b101u);
  store.rotation(p1.store_row()) = static_cast<std::uint8_t>(0b010u);
  BOOST_CHECK_EQUAL(store.rotation(p0.store_row()),
                    static_cast<std::uint8_t>(0b101u));
  std::uint8_t &rot_ref = store.rotation(p1.store_row());
  rot_ref = static_cast<std::uint8_t>(0b111u);
  BOOST_CHECK_EQUAL(store.rotation(p1.store_row()),
                    static_cast<std::uint8_t>(0b111u));
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  store.ext_flag(p0.store_row()) = static_cast<std::uint8_t>(0b011u);
  store.ext_flag(p1.store_row()) = static_cast<std::uint8_t>(0b100u);
#endif

  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

#ifdef ESPRESSO_ROTATION
  BOOST_CHECK_EQUAL(store.rotation(p1.store_row()),
                    static_cast<std::uint8_t>(0b111u));
  BOOST_CHECK_EQUAL(store.rotation(p0.store_row()),
                    static_cast<std::uint8_t>(0b101u));
  BOOST_CHECK_EQUAL(store.rotation(p2.store_row()),
                    static_cast<std::uint8_t>(0u)); // default
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  BOOST_CHECK_EQUAL(store.ext_flag(p1.store_row()),
                    static_cast<std::uint8_t>(0b100u));
  BOOST_CHECK_EQUAL(store.ext_flag(p0.store_row()),
                    static_cast<std::uint8_t>(0b011u));
  BOOST_CHECK_EQUAL(store.ext_flag(p2.store_row()),
                    static_cast<std::uint8_t>(0u)); // default
#endif
}
#endif

// Phase-5 host sidecar (POD std::vector): exercises the NEW sidecar machinery.
// A POD written into a sidecar row must preserve by old row across a rebuild
// and a genuinely-new row must default to a value-constructed POD.
#ifdef ESPRESSO_ENGINE
BOOST_AUTO_TEST_CASE(swimming_sidecar_preserve_and_default) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();

  store.swimming(p0.store_row()).f_swim = 1.5;
  store.swimming(p0.store_row()).swimming = true;
  store.swimming(p1.store_row()).f_swim = -2.5;

  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

  BOOST_CHECK_EQUAL(store.swimming(p1.store_row()).f_swim, -2.5);
  BOOST_CHECK_EQUAL(store.swimming(p0.store_row()).f_swim, 1.5);
  BOOST_CHECK(store.swimming(p0.store_row()).swimming);
  // genuinely-new row: value-constructed POD default (matches the carrier /
  // ParticleParametersSwimming member defaults)
  BOOST_CHECK_EQUAL(store.swimming(p2.store_row()).f_swim, 0.);
  BOOST_CHECK(not store.swimming(p2.store_row()).swimming);
}
#endif

// Phase-5: a genuinely new row's parameter columns/sidecars must default to
// ParticleProperties' member defaults, seeded from the migration carriers of a
// freshly-constructed (detached) particle.
BOOST_AUTO_TEST_CASE(new_row_parameter_defaults) {
  ParticleStore store{};
  Particle p{};
  BOOST_REQUIRE(p.store() == nullptr);
  store.begin_rebuild(1u, 0u);
  store.assign_row(p, 0);
  store.finish_rebuild();

  BOOST_CHECK_EQUAL(store.id(p.store_row()), -1);
  BOOST_CHECK_EQUAL(store.mol_id(p.store_row()), 0);
  BOOST_CHECK_EQUAL(store.type(p.store_row()), 0);
  BOOST_CHECK_EQUAL(store.propagation(p.store_row()),
                    static_cast<int>(PropagationMode::SYSTEM_DEFAULT));
#ifdef ESPRESSO_MASS
  BOOST_CHECK_EQUAL(store.mass(p.store_row()), 1.0);
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  BOOST_CHECK_EQUAL(store.q(p.store_row()), 0.0);
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Utils::Vector3d const ri = store.rinertia_value(p.store_row());
  BOOST_CHECK_EQUAL(ri[0], 1.);
  BOOST_CHECK_EQUAL(ri[1], 1.);
  BOOST_CHECK_EQUAL(ri[2], 1.);
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  BOOST_CHECK_EQUAL(store.ext_force_value(p.store_row()).norm2(), 0.);
  BOOST_CHECK_EQUAL(store.ext_flag(p.store_row()),
                    static_cast<std::uint8_t>(0u));
#endif
#ifdef ESPRESSO_ROTATION
  BOOST_CHECK_EQUAL(store.rotation(p.store_row()),
                    static_cast<std::uint8_t>(0u));
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  Utils::Vector3d const g = store.gamma_value(p.store_row());
  BOOST_CHECK_EQUAL(g[0], -1.);
  BOOST_CHECK_EQUAL(g[2], -1.);
#else
  BOOST_CHECK_EQUAL(store.gamma_value(p.store_row()), -1.);
#endif
#endif
}

// Phase-5: a detached particle (freshly constructed) carries its parameters in
// the ParticleProperties member; assign_row seeds the columns/sidecars from the
// dormant migration carriers (which read those struct fields pre-flip). Set the
// values on the detached particle, then attach it and verify the columns.
BOOST_AUTO_TEST_CASE(rebuild_seeds_parameters_from_carrier) {
  Particle p{};
  BOOST_REQUIRE(p.store() == nullptr);
  p.id() = 42;
  p.type() = 5;
#ifdef ESPRESSO_MASS
  p.mass() = 3.5;
#endif
#ifdef ESPRESSO_ENGINE
  p.swimming().f_swim = 9.5;
#endif
  // The dormant carriers read the struct fields directly.
  BOOST_CHECK_EQUAL(p.migration_id(), 42);
  BOOST_CHECK_EQUAL(p.migration_type(), 5);

  ParticleStore store{};
  store.begin_rebuild(1u, 0u);
  store.assign_row(p, 0);
  store.finish_rebuild();

  BOOST_CHECK_EQUAL(store.id(p.store_row()), 42);
  BOOST_CHECK_EQUAL(store.type(p.store_row()), 5);
#ifdef ESPRESSO_MASS
  BOOST_CHECK_EQUAL(store.mass(p.store_row()), 3.5);
#endif
#ifdef ESPRESSO_ENGINE
  BOOST_CHECK_EQUAL(store.swimming(p.store_row()).f_swim, 9.5);
#endif
}

// Scalar column references write through to the stored value.
BOOST_AUTO_TEST_CASE(scalar_column_references_write_through) {
  ParticleStore store{};
  Particle p{};
  store.begin_rebuild(1u, 0u);
  store.assign_row(p, 0);
  store.finish_rebuild();

  store.lees_edwards_offset(p.store_row()) = 3.14;
  BOOST_CHECK_EQUAL(store.lees_edwards_offset(p.store_row()), 3.14);
  double &offset_ref = store.lees_edwards_offset(p.store_row());
  offset_ref = -1.0;
  BOOST_CHECK_EQUAL(store.lees_edwards_offset(p.store_row()), -1.0);

  store.lees_edwards_flag(p.store_row()) = static_cast<short>(5);
  BOOST_CHECK_EQUAL(store.lees_edwards_flag(p.store_row()),
                    static_cast<short>(5));
  short &flag_ref = store.lees_edwards_flag(p.store_row());
  flag_ref = static_cast<short>(-7);
  BOOST_CHECK_EQUAL(store.lees_edwards_flag(p.store_row()),
                    static_cast<short>(-7));
}
