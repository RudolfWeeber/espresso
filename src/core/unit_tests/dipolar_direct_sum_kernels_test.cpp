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

#define BOOST_TEST_MODULE "dipolar direct sum kernels"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "magnetostatics/dipolar_direct_sum_kernels.hpp"

#include <utils/Vector.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>

#include <cstddef>
#include <vector>

// Kokkos must be initialized before any simd/view use.
struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};
BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

namespace {
// Build a simd_double whose lane l holds vals[l].
simd_double make_simd(double const *vals) {
  return simd_double(vals, Kokkos::Experimental::simd_flag_default);
}
} // namespace

BOOST_AUTO_TEST_SUITE(suite)

// Anchor: two identical dipoles m=(1,0,0) separated by d=(1,0,0).
// Analytically f=(-6,0,0), torque=0 (parallel moments => no torque).
BOOST_AUTO_TEST_CASE(pair_force_analytic) {
  Utils::Vector3d const d{1., 0., 0.};
  Utils::Vector3d const m{1., 0., 0.};
  auto const pf = pair_force<double>(d, m, m);
  BOOST_CHECK_CLOSE(pf.f[0], -6., 1e-10);
  BOOST_CHECK_SMALL(pf.f[1], 1e-12);
  BOOST_CHECK_SMALL(pf.f[2], 1e-12);
  BOOST_CHECK_SMALL(pf.torque[0], 1e-12);
  BOOST_CHECK_SMALL(pf.torque[1], 1e-12);
  BOOST_CHECK_SMALL(pf.torque[2], 1e-12);
}

// Core SIMD-correctness property: each lane of the simd kernel equals the
// scalar kernel evaluated on that lane's inputs, for distinct per-lane data.
BOOST_AUTO_TEST_CASE(pair_force_simd_matches_scalar) {
  constexpr std::size_t w = simd_double::size();
  // Distinct, non-degenerate inputs per lane.
  std::vector<Utils::Vector3d> ds(w), m1s(w), m2s(w);
  for (std::size_t l = 0; l < w; ++l) {
    auto const x = 1. + 0.37 * double(l);
    ds[l] = {x, 0.5 - 0.1 * double(l), 0.2 + 0.05 * double(l)};
    m1s[l] = {0.3 + 0.2 * double(l), 1.1, -0.4};
    m2s[l] = {-0.7, 0.9 - 0.05 * double(l), 0.6};
  }
  // Pack per-component lane buffers into simd inputs.
  auto pack = [&](std::vector<Utils::Vector3d> const &v, int c) {
    std::vector<double> buf(w);
    for (std::size_t l = 0; l < w; ++l)
      buf[l] = v[l][c];
    return make_simd(buf.data());
  };
  Utils::Vector<simd_double, 3> const d_s{pack(ds, 0), pack(ds, 1),
                                          pack(ds, 2)};
  Utils::Vector<simd_double, 3> const m1_s{pack(m1s, 0), pack(m1s, 1),
                                           pack(m1s, 2)};
  Utils::Vector<simd_double, 3> const m2_s{pack(m2s, 0), pack(m2s, 1),
                                           pack(m2s, 2)};

  auto const pf_s = pair_force<simd_double>(d_s, m1_s, m2_s);

  for (std::size_t l = 0; l < w; ++l) {
    auto const pf = pair_force<double>(ds[l], m1s[l], m2s[l]);
    for (int c = 0; c < 3; ++c) {
      BOOST_CHECK_CLOSE(pf_s.f[c][l], pf.f[c], 1e-10);
      BOOST_CHECK_CLOSE(pf_s.torque[c][l], pf.torque[c], 1e-10);
    }
  }
}

BOOST_AUTO_TEST_CASE(pair_potential_simd_matches_scalar) {
  constexpr std::size_t w = simd_double::size();
  std::vector<Utils::Vector3d> ds(w), m1s(w), m2s(w);
  for (std::size_t l = 0; l < w; ++l) {
    ds[l] = {1. + 0.3 * double(l), 0.4, -0.2};
    m1s[l] = {0.5, -0.3 + 0.1 * double(l), 0.8};
    m2s[l] = {0.2, 0.7, -0.6 + 0.05 * double(l)};
  }
  auto pack = [&](std::vector<Utils::Vector3d> const &v, int c) {
    std::vector<double> buf(w);
    for (std::size_t l = 0; l < w; ++l)
      buf[l] = v[l][c];
    return make_simd(buf.data());
  };
  Utils::Vector<simd_double, 3> const d_s{pack(ds, 0), pack(ds, 1),
                                          pack(ds, 2)};
  Utils::Vector<simd_double, 3> const m1_s{pack(m1s, 0), pack(m1s, 1),
                                           pack(m1s, 2)};
  Utils::Vector<simd_double, 3> const m2_s{pack(m2s, 0), pack(m2s, 1),
                                           pack(m2s, 2)};
  auto const u_s = pair_potential<simd_double>(d_s, m1_s, m2_s);
  for (std::size_t l = 0; l < w; ++l) {
    auto const u = pair_potential<double>(ds[l], m1s[l], m2s[l]);
    BOOST_CHECK_CLOSE(u_s[l], u, 1e-10);
  }
}

BOOST_AUTO_TEST_SUITE_END()
