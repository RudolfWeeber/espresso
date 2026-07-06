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

#define BOOST_TEST_MODULE "Particle struct test"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <config/config.hpp>

#include "Particle.hpp"
#include "ParticleStoreTestFixture.hpp"
#include "PropagationMode.hpp"

#include <utils/compact_vector.hpp>
#include <utils/serialization/memcpy_archive.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <algorithm>
#include <array>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

void check_particle_force(ParticleForce const &out, ParticleForce const &ref) {
  BOOST_TEST(out.f == ref.f, boost::test_tools::per_element());
#ifdef ESPRESSO_ROTATION
  BOOST_TEST(out.torque == ref.torque, boost::test_tools::per_element());
#endif
}

BOOST_AUTO_TEST_CASE(comparison) {
  {
    Particle p, q;

    p.id() = 1;
    q.id() = 2;

    BOOST_CHECK(p != q);
    BOOST_CHECK(not(p == q));
  }

  {
    Particle p, q;

    p.id() = 2;
    q.id() = 2;

    BOOST_CHECK(not(p != q));
    BOOST_CHECK(p == q);
  }
}

BOOST_AUTO_TEST_CASE(serialization) {
  // Migration phase 2: force/torque live in the ParticleStore columns and are
  // no longer part of Particle serialization. Attach the hand-made particles
  // to a standalone store so the accessors work.
  ParticleStoreTestFixture fixture{};
  auto p = Particle();
  fixture.attach(p);

  auto const bond_id = 5;
  auto const bond_partners = std::array<const int, 3>{{12, 13, 14}};

  p.id() = 15;
  p.bonds().insert({bond_id, bond_partners});
  p.force() = {1., -2., 3.};
#ifdef ESPRESSO_ROTATION
  p.torque() = {-4., 5., -6.};
#endif
#ifdef ESPRESSO_EXCLUSIONS
  std::vector<int> el = {5, 6, 7, 8};
  p.exclusions() = Utils::compact_vector<int>{el.begin(), el.end()};
#endif

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << p;

  boost::archive::text_iarchive in_ar(stream);
  // Migration phase 5: id (with all other parameters) lives in the store column
  // and is ferried across the archive via the migration carrier. Deserialize
  // into a DETACHED q so the carrier value is read back (the assign_row seeding
  // path), then attach so the column read returns it -- mirrors the state and
  // momentum carrier round-trips below.
  auto q = Particle();
  in_ar >> q;

  BOOST_CHECK(q.detached_id() == p.detached_id());
  BOOST_CHECK((*q.bonds().begin() == BondView{bond_id, bond_partners}));
  fixture.attach(q);
  BOOST_CHECK(q.id() == p.id());
  // Force/torque ARE serialized via the migration carriers (m_detached_force /
  // m_detached_torque); after attaching q to the fresh store, assign_row seeds
  // row 0 from those carriers, so the column reads return the ferried values.
  BOOST_TEST(Utils::Vector3d(q.force()) == (Utils::Vector3d{1., -2., 3.}),
             boost::test_tools::per_element());
#ifdef ESPRESSO_ROTATION
  BOOST_TEST(Utils::Vector3d(q.torque()) == (Utils::Vector3d{-4., 5., -6.}),
             boost::test_tools::per_element());
#endif
}

// Migration phase 3: the STATE fields (position, image box, quaternion,
// position-at-last-verlet-update, position-at-last-time-step, Lees-Edwards
// offset and flag) live in the ParticleStore columns and are ferried across a
// boost archive via the migration carriers. This round-trip test attaches p to
// a store, sets known state, serializes, and deserializes into a DETACHED q so
// that the carrier values (not a column) are read back through the detached_*()
// getters -- exactly the path ParticleStore::assign_row uses to seed a migrated
// particle's new row on the destination rank.
BOOST_AUTO_TEST_CASE(state_carrier_serialization_round_trip) {
  ParticleStoreTestFixture fixture{};
  auto p = Particle();
  fixture.attach(p);

  auto const ref_pos = Utils::Vector3d{1.25, -2.5, 3.75};
  auto const ref_image = Utils::Vector3i{4, -5, 6};
  p.pos() = ref_pos;
  p.image_box() = ref_image;
  p.pos_at_last_verlet_update() = Utils::Vector3d{0.5, 0.25, -0.125};
#ifdef ESPRESSO_ROTATION
  auto const ref_quat = Utils::Quaternion<double>{{0.5, 0.5, 0.5, 0.5}};
  p.quat() = ref_quat;
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  auto const ref_pos_last = Utils::Vector3d{7.5, 8.5, 9.5};
  p.pos_last_time_step() = ref_pos_last;
#endif
  p.lees_edwards_offset() = 1.5;
  p.lees_edwards_flag() = short{3};

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << p;

  boost::archive::text_iarchive in_ar(stream);
  auto q = Particle(); // intentionally NOT attached to a store
  in_ar >> q;

  // A detached q reads the ferried carriers through detached_*().
  BOOST_TEST(q.detached_position() == ref_pos,
             boost::test_tools::per_element());
  BOOST_TEST(q.detached_image_box() == ref_image,
             boost::test_tools::per_element());
  BOOST_TEST(q.detached_position_at_last_verlet_update() ==
                 (Utils::Vector3d{0.5, 0.25, -0.125}),
             boost::test_tools::per_element());
#ifdef ESPRESSO_ROTATION
  {
    auto const q_quat = q.detached_quaternion();
    for (std::size_t j = 0u; j < 4u; ++j) {
      BOOST_CHECK_EQUAL(q_quat[j], ref_quat[j]);
    }
  }
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  BOOST_TEST(q.detached_position_last_time_step() == ref_pos_last,
             boost::test_tools::per_element());
#endif
  BOOST_CHECK_EQUAL(q.detached_lees_edwards_offset(), 1.5);
  BOOST_CHECK_EQUAL(q.detached_lees_edwards_flag(), short{3});

  // After attaching q to a fresh store, assign_row seeds the row from those
  // carriers, so the column reads now return the ferried values.
  fixture.attach(q);
  BOOST_TEST(Utils::Vector3d(q.pos()) == ref_pos,
             boost::test_tools::per_element());
  BOOST_TEST(Utils::Vector3i(q.image_box()) == ref_image,
             boost::test_tools::per_element());
#ifdef ESPRESSO_ROTATION
  {
    auto const q_quat = Utils::Quaternion<double>(q.quat());
    for (std::size_t j = 0u; j < 4u; ++j) {
      BOOST_CHECK_EQUAL(q_quat[j], ref_quat[j]);
    }
  }
#endif
  BOOST_CHECK_EQUAL(q.lees_edwards_offset(), 1.5);
  BOOST_CHECK_EQUAL(q.lees_edwards_flag(), short{3});
}

// Migration phase 4: the MOMENTUM fields (velocity and angular velocity) live
// in the ParticleStore columns and are ferried across a boost archive via the
// migration carriers. Mirror of state_carrier_serialization_round_trip for the
// momentum fields: attach p to a store, set known velocity/omega, serialize,
// deserialize into a DETACHED q (reads the carriers through detached_*()), then
// attach q to a fresh store and confirm assign_row seeds the columns from the
// ferried carriers -- exactly the cross-rank migration seeding path.
BOOST_AUTO_TEST_CASE(momentum_carrier_serialization_round_trip) {
  ParticleStoreTestFixture fixture{};
  auto p = Particle();
  fixture.attach(p);

  auto const ref_vel = Utils::Vector3d{-1.5, 2.25, -3.5};
  p.v() = ref_vel;
#ifdef ESPRESSO_ROTATION
  auto const ref_omega = Utils::Vector3d{0.75, -1.25, 4.5};
  p.omega() = ref_omega;
#endif

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << p;

  boost::archive::text_iarchive in_ar(stream);
  auto q = Particle(); // intentionally NOT attached to a store
  in_ar >> q;

  // A detached q reads the ferried carriers through detached_*().
  BOOST_TEST(q.detached_velocity() == ref_vel,
             boost::test_tools::per_element());
#ifdef ESPRESSO_ROTATION
  BOOST_TEST(q.detached_angular_velocity() == ref_omega,
             boost::test_tools::per_element());
#endif

  // After attaching q to a fresh store, assign_row seeds the row from those
  // carriers, so the column reads now return the ferried values.
  fixture.attach(q);
  BOOST_TEST(Utils::Vector3d(q.v()) == ref_vel,
             boost::test_tools::per_element());
#ifdef ESPRESSO_ROTATION
  BOOST_TEST(Utils::Vector3d(q.omega()) == ref_omega,
             boost::test_tools::per_element());
#endif
}

namespace Utils {
template <>
struct is_statically_serializable<ParticleProperties> : std::true_type {};
} // namespace Utils

BOOST_AUTO_TEST_CASE(properties_serialization) {
  auto const expected_size =
      Utils::MemcpyOArchive::packing_size<ParticleProperties>();

  BOOST_CHECK_LE(expected_size, sizeof(ParticleProperties));

  std::vector<char> buf(expected_size);

  auto prop = ParticleProperties{};
  prop.identity = 1234;

  {
    auto oa = Utils::MemcpyOArchive{buf};

    oa << prop;

    BOOST_CHECK_EQUAL(oa.bytes_written(), expected_size);
  }

  {
    auto ia = Utils::MemcpyIArchive{buf};
    ParticleProperties out;

    ia >> out;
    BOOST_CHECK_EQUAL(ia.bytes_read(), expected_size);
    BOOST_CHECK_EQUAL(out.identity, prop.identity);
  }
}

namespace Utils {
template <>
struct is_statically_serializable<ParticleForce> : std::true_type {};
} // namespace Utils

BOOST_AUTO_TEST_CASE(force_serialization) {
  auto const expected_size =
      Utils::MemcpyOArchive::packing_size<ParticleForce>();

  BOOST_CHECK_LE(expected_size, sizeof(ParticleForce));

  std::vector<char> buf(expected_size);

  auto pf = ParticleForce{{1, 2, 3}};
#ifdef ESPRESSO_ROTATION
  pf.torque = {4, 5, 6};
#endif

  {
    auto oa = Utils::MemcpyOArchive{buf};

    oa << pf;

    BOOST_CHECK_EQUAL(oa.bytes_written(), expected_size);
  }

  {
    auto ia = Utils::MemcpyIArchive{buf};
    ParticleForce out;

    ia >> out;

    BOOST_CHECK_EQUAL(ia.bytes_read(), expected_size);
    check_particle_force(out, pf);
  }
}

BOOST_AUTO_TEST_CASE(force_constructors) {

  auto pf = ParticleForce{{1, 2, 3}};
#ifdef ESPRESSO_ROTATION
  pf.torque = {4, 5, 6};
#endif

  // check copy constructor
  {
    ParticleForce out(pf);
    check_particle_force(out, pf);
  }

  // check copy assignment operator
  {
    ParticleForce out; // avoid copy elision
    out = pf;
    check_particle_force(out, pf);
  }
}

// Migration phase 6: ParticleRattle no longer carries the correction Vector3d
// (evicted to a ParticleStore observable column); the struct is now an empty
// type anchor. Its standalone serialization/constructor unit tests are removed
// -- the correction round-trip is exercised by the store column tests
// (ParticleStore_test.cpp) and the RATTLE ghost path.

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH

void check_particle_tsw(ThermalStonerWohlfarthParameters const &out,
                        ThermalStonerWohlfarthParameters const &ref) {
  BOOST_TEST(out.is_enabled == ref.is_enabled);
  BOOST_TEST(out.phi0 == ref.phi0);
  BOOST_TEST(out.sat_mag == ref.sat_mag);
  BOOST_TEST(out.ani_fld_inv == ref.ani_fld_inv);
  BOOST_TEST(out.ani_energy == ref.ani_energy);
  BOOST_TEST(out.tau0_inv == ref.tau0_inv);
  BOOST_TEST(out.dt_incr == ref.dt_incr);
}

BOOST_AUTO_TEST_CASE(thermal_stoner_wohlfarth_serialization) {
  auto const expected_size =
      Utils::MemcpyOArchive::packing_size<ThermalStonerWohlfarthParameters>();

  BOOST_CHECK_LE(expected_size, sizeof(ThermalStonerWohlfarthParameters));

  std::vector<char> buf(expected_size);

  auto pr = ThermalStonerWohlfarthParameters{.is_enabled = false,
                                             .phi0 = 1.,
                                             .sat_mag = 2.,
                                             .ani_fld_inv = 3.,
                                             .ani_energy = 4.,
                                             .tau0_inv = 5.,
                                             .dt_incr = 6.};

  {
    auto oa = Utils::MemcpyOArchive{buf};

    oa << pr;

    BOOST_CHECK_EQUAL(oa.bytes_written(), expected_size);
  }

  {
    auto ia = Utils::MemcpyIArchive{buf};
    ThermalStonerWohlfarthParameters out;

    ia >> out;

    BOOST_CHECK_EQUAL(ia.bytes_read(), expected_size);
    check_particle_tsw(out, pr);
  }
}

BOOST_AUTO_TEST_CASE(thermal_stoner_wohlfarth_constructors) {
  auto pr = ThermalStonerWohlfarthParameters{.is_enabled = false,
                                             .phi0 = 1.,
                                             .sat_mag = 2.,
                                             .ani_fld_inv = 3.,
                                             .ani_energy = 4.,
                                             .tau0_inv = 5.,
                                             .dt_incr = 6.};

  // check copy constructor
  {
    ThermalStonerWohlfarthParameters out(pr);
    check_particle_tsw(out, pr);
  }

  // check copy assignment operator
  {
    ThermalStonerWohlfarthParameters out; // avoid copy elision
    out = pr;
    check_particle_tsw(out, pr);
  }
}
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH

BOOST_AUTO_TEST_CASE(particle_bitfields) {
  auto p = Particle();

  // check default values
  BOOST_CHECK(not p.has_fixed_coordinates());
  BOOST_CHECK(not p.can_rotate());
  BOOST_CHECK(not p.is_fixed_along(1));
  BOOST_CHECK(not p.can_rotate_around(1));

  // check setting of one axis
#ifdef ESPRESSO_EXTERNAL_FORCES
  p.set_fixed_along(1, true);
  BOOST_CHECK(p.is_fixed_along(1));
  BOOST_CHECK(p.has_fixed_coordinates());
#endif
#ifdef ESPRESSO_ROTATION
  p.set_can_rotate_around(1, true);
  BOOST_CHECK(p.can_rotate_around(1));
  BOOST_CHECK(p.can_rotate());
#endif

  // check that unsetting is properly registered
#ifdef ESPRESSO_EXTERNAL_FORCES
  p.set_fixed_along(1, false);
  BOOST_CHECK(not p.has_fixed_coordinates());
#endif
#ifdef ESPRESSO_ROTATION
  p.set_can_rotate_around(1, false);
  BOOST_CHECK(not p.can_rotate());
#endif

  // check setting of all flags at once
#ifdef ESPRESSO_ROTATION
  p.set_can_rotate_all_axes();
  BOOST_CHECK(p.can_rotate_around(0));
  BOOST_CHECK(p.can_rotate_around(1));
  BOOST_CHECK(p.can_rotate_around(2));
  p.set_cannot_rotate_all_axes();
  BOOST_CHECK(not p.can_rotate());
#endif

  static_assert(
      std::is_same_v<std::underlying_type_t<PropagationMode::PropagationMode>,
                     decltype(ParticleProperties::propagation)>);
}
