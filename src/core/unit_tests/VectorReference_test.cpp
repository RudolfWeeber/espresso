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
#define BOOST_TEST_MODULE VectorReference test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "particle_store/VectorReference.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cstddef>

// simulate a component-major (LayoutLeft) column with 4 rows:
// memory is [x0 x1 x2 x3 | y0 y1 y2 y3 | z0 z1 z2 z3], stride 4
struct ColumnFixture {
  std::array<double, 12> storage{};
  static constexpr std::size_t stride = 4;
  VectorReference row(std::size_t i) {
    return VectorReference(storage.data() + i, stride);
  }
};

BOOST_FIXTURE_TEST_CASE(assignment_writes_through, ColumnFixture) {
  row(1) = Utils::Vector3d{1., 2., 3.};
  BOOST_CHECK_EQUAL(storage[1], 1.); // x1
  BOOST_CHECK_EQUAL(storage[5], 2.); // y1
  BOOST_CHECK_EQUAL(storage[9], 3.); // z1
}

BOOST_FIXTURE_TEST_CASE(conversion_reads_components, ColumnFixture) {
  storage[2] = 4.;  // x2
  storage[6] = 5.;  // y2
  storage[10] = 6.; // z2
  Utils::Vector3d const value = row(2);
  BOOST_CHECK_EQUAL(value[0], 4.);
  BOOST_CHECK_EQUAL(value[1], 5.);
  BOOST_CHECK_EQUAL(value[2], 6.);
}

BOOST_FIXTURE_TEST_CASE(compound_operators, ColumnFixture) {
  row(0) = Utils::Vector3d{1., 1., 1.};
  row(0) += Utils::Vector3d{1., 2., 3.};
  row(0) -= Utils::Vector3d{0., 1., 0.};
  row(0) *= 2.;
  Utils::Vector3d const value = row(0);
  BOOST_CHECK_EQUAL(value[0], 4.);
  BOOST_CHECK_EQUAL(value[1], 4.);
  BOOST_CHECK_EQUAL(value[2], 8.);
}

BOOST_FIXTURE_TEST_CASE(subscript_and_norms, ColumnFixture) {
  row(3)[0] = 3.;
  row(3)[1] = 0.;
  row(3)[2] = 4.;
  BOOST_CHECK_EQUAL(row(3)[2], 4.);
  BOOST_CHECK_EQUAL(row(3).norm2(), 25.);
  BOOST_CHECK_EQUAL(row(3).norm(), 5.);
}

BOOST_FIXTURE_TEST_CASE(proxy_copies_alias_the_same_row, ColumnFixture) {
  auto reference = row(1);
  auto alias = reference; // copying the proxy copies the reference, not data
  alias = Utils::Vector3d{7., 8., 9.};
  BOOST_CHECK_EQUAL(storage[1], 7.);
}
