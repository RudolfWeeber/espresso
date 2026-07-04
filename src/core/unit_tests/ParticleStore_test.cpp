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

#include <sstream>

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
#ifdef ESPRESSO_ROTATION
  auto const q0 = Utils::Quaternion<double>{{0.5, 0.5, 0.5, 0.5}};
  auto const q1 = Utils::Quaternion<double>{{0., 1., 0., 0.}};
  store.quaternion_reference(p0.store_row()) = q0;
  store.quaternion_reference(p1.store_row()) = q1;
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
#ifdef ESPRESSO_ROTATION
  Utils::Quaternion<double> const quat_p1 =
      store.quaternion_value(p1.store_row());
  Utils::Quaternion<double> const quat_p0 =
      store.quaternion_value(p0.store_row());
  BOOST_CHECK_EQUAL(quat_p1[1], 1.);
  BOOST_CHECK_EQUAL(quat_p0[3], 0.5);
#endif
}

// A genuinely new row (e.g. a fresh ghost, before its field update overwrites
// the column) must default to sensible values that match Particle's member
// defaults. The quaternion in particular must be the IDENTITY (1,0,0,0), NOT
// the Kokkos zero-init (0,0,0,0) which is an invalid quaternion.
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
#ifdef ESPRESSO_ROTATION
  Utils::Quaternion<double> const quat = store.quaternion_value(p.store_row());
  BOOST_CHECK_EQUAL(quat[0], 1.); // identity, not zero
  BOOST_CHECK_EQUAL(quat[1], 0.);
  BOOST_CHECK_EQUAL(quat[2], 0.);
  BOOST_CHECK_EQUAL(quat[3], 0.);
#endif
}

// A detached particle with its migration carriers populated must have all its
// state columns seeded from those carriers on assign_row (models a particle
// that just migrated in from another rank). Pre-flip the carrier getters read
// the struct fields, so we set the struct via the public setters.
BOOST_AUTO_TEST_CASE(rebuild_seeds_migrated_particle_state_from_carrier) {
  Particle p{};
  p.id() = 11;
  p.pos() = Utils::Vector3d{1.5, -2.5, 3.5};
  p.image_box() = Utils::Vector3i{4, -5, 6};
  p.pos_at_last_verlet_update() = Utils::Vector3d{7.5, 8.5, 9.5};
  p.lees_edwards_offset() = 12.75;
  p.lees_edwards_flag() = static_cast<short>(2);
#ifdef ESPRESSO_ROTATION
  auto const q = Utils::Quaternion<double>{{0., 0., 1., 0.}};
  p.quat() = q;
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  p.pos_last_time_step() = Utils::Vector3d{-1., -2., -3.};
#endif
  BOOST_REQUIRE(p.store() == nullptr);

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
