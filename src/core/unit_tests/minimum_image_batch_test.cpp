/*
 * Copyright (C) 2026 The ESPResSo project
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

#define BOOST_TEST_MODULE minimum_image_batch tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "BoxGeometry.hpp"
#include "algorithm/minimum_image_batch.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

namespace {

// Bit-for-bit equality of two doubles (NaN not expected here).
bool bit_equal(double a, double b) {
  std::uint64_t ua, ub;
  static_assert(sizeof(double) == sizeof(std::uint64_t));
  __builtin_memcpy(&ua, &a, sizeof(double));
  __builtin_memcpy(&ub, &b, sizeof(double));
  return ua == ub;
}

BoxGeometry make_box(Utils::Vector3d const &l, bool px, bool py, bool pz) {
  BoxGeometry box;
  box.set_type(BoxType::CUBOID);
  box.set_length(l);
  box.set_periodic(0u, px);
  box.set_periodic(1u, py);
  box.set_periodic(2u, pz);
  return box;
}

// Verify the batch primitives reproduce the scalar get_mi_vector /
// get_mi_dist2 BITWISE for a batch of neighbor positions, including the
// exact-boundary case abs(dx) == box_half.
void check_identity(BoxGeometry const &box, Utils::Vector3d const &ref,
                    std::vector<Utils::Vector3d> const &neigh) {
  auto const n = neigh.size();
  std::vector<double> nx(n), ny(n), nz(n);
  for (std::size_t k = 0; k < n; ++k) {
    nx[k] = neigh[k][0];
    ny[k] = neigh[k][1];
    nz[k] = neigh[k][2];
  }
  std::vector<double> ox(n), oy(n), oz(n), od2(n);
  auto const params = Algorithm::OrthoBoxParams::from(box);
  std::array<double, 3> ref_arr{ref[0], ref[1], ref[2]};

  Algorithm::get_mi_vector_batch(params, ref_arr.data(), nx.data(), ny.data(),
                                 nz.data(), ox.data(), oy.data(), oz.data(), n);
  Algorithm::get_mi_dist2_batch(params, ref_arr.data(), nx.data(), ny.data(),
                                nz.data(), od2.data(), n);

  for (std::size_t k = 0; k < n; ++k) {
    auto const d_scalar = box.get_mi_vector(ref, neigh[k]);
    auto const d2_scalar = box.get_mi_dist2(ref, neigh[k]);
    BOOST_TEST(bit_equal(ox[k], d_scalar[0]));
    BOOST_TEST(bit_equal(oy[k], d_scalar[1]));
    BOOST_TEST(bit_equal(oz[k], d_scalar[2]));
    BOOST_TEST(bit_equal(od2[k], d2_scalar));
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(bitwise_identity_periodic) {
  auto const box = make_box({10.0, 8.0, 6.5}, true, true, true);
  std::mt19937 gen(12345);
  std::uniform_real_distribution<double> ux(-15.0, 15.0);
  std::uniform_real_distribution<double> uy(-12.0, 12.0);
  std::uniform_real_distribution<double> uz(-10.0, 10.0);
  Utils::Vector3d ref{ux(gen), uy(gen), uz(gen)};
  std::vector<Utils::Vector3d> neigh;
  for (int i = 0; i < 500; ++i) {
    neigh.push_back({ux(gen), uy(gen), uz(gen)});
  }
  check_identity(box, ref, neigh);
}

BOOST_AUTO_TEST_CASE(bitwise_identity_mixed_periodicity) {
  auto const box = make_box({10.0, 8.0, 6.5}, true, false, true);
  std::mt19937 gen(999);
  std::uniform_real_distribution<double> u(-20.0, 20.0);
  Utils::Vector3d ref{u(gen), u(gen), u(gen)};
  std::vector<Utils::Vector3d> neigh;
  for (int i = 0; i < 500; ++i) {
    neigh.push_back({u(gen), u(gen), u(gen)});
  }
  check_identity(box, ref, neigh);
}

BOOST_AUTO_TEST_CASE(bitwise_identity_boundary) {
  // Exercise dx exactly at +/- box_half, where a naive branchless fold would
  // diverge from the scalar `abs(dx) > box_half` conditional.
  auto const box = make_box({10.0, 8.0, 6.5}, true, true, true);
  Utils::Vector3d ref{0.0, 0.0, 0.0};
  std::vector<Utils::Vector3d> neigh;
  // Neighbor at -box_half in each axis -> dx = +box_half exactly.
  neigh.push_back({-5.0, -4.0, -3.25});
  // Neighbor at +box_half -> dx = -box_half exactly.
  neigh.push_back({5.0, 4.0, 3.25});
  // Just past the boundary in both directions.
  neigh.push_back({-5.000001, -4.000001, -3.250001});
  neigh.push_back({5.000001, 4.000001, 3.250001});
  check_identity(box, ref, neigh);
}

BOOST_AUTO_TEST_CASE(bitwise_identity_round_ties) {
  // Exercise dx at n*box +/- box_half, i.e. dx*box_inv == n +/- 0.5 exactly.
  // Here std::round (half-away) and rint (half-even) disagree, so the batch's
  // tie-corrected round MUST match the scalar std::round bitwise. Covers ties
  // that rint rounds toward zero (n even) and away (n odd), both signs.
  auto const box = make_box({10.0, 8.0, 6.5}, true, true, true);
  Utils::Vector3d ref{0.0, 0.0, 0.0};
  std::vector<Utils::Vector3d> neigh;
  for (int n = -4; n <= 4; ++n) {
    for (double s : {0.5, -0.5}) {
      // neighbor = -(n+s)*box so that dx = ref - neighbor = (n+s)*box.
      neigh.push_back({-(n + s) * 10.0, -(n + s) * 8.0, -(n + s) * 6.5});
      // one ulp on either side of the tie.
      double const y = (n + s);
      double const yx = std::nextafter(y, 1e9);
      double const yn = std::nextafter(y, -1e9);
      neigh.push_back({-yx * 10.0, -yn * 8.0, -yx * 6.5});
      neigh.push_back({-yn * 10.0, -yx * 8.0, -yn * 6.5});
    }
  }
  check_identity(box, ref, neigh);
}
