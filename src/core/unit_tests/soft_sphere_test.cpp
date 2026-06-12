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

#define BOOST_TEST_MODULE soft sphere potential
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "config/config.hpp"

#ifdef ESPRESSO_SOFT_SPHERE

#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "nonbonded_interactions/soft_sphere.hpp"

#include <cmath>

/**
 * The soft-sphere force kernel treats a pair with
 * @c r_off = dist - offset <= 0 as non-interacting (it returns a force factor
 * of exactly 0 there, see the @c r_off > 0 guard in
 * @ref soft_pair_force_factor). The energy must agree: wherever the force is
 * 0, the energy must be finite and equal to 0. Without the matching guard the
 * energy evaluates @c a / pow(r_off, n), which is +Inf for @c r_off == 0 and
 * NaN for @c r_off < 0 with a non-integer exponent. This is regression test
 * for bug-sweep #15.
 */
BOOST_AUTO_TEST_CASE(soft_sphere_energy_force_consistency) {
  // offset > 0, so max_cutoff() == cut + offset and a pair closer than the
  // offset still passes the (dist < max_cutoff()) gate inside the energy.
  {
    IA_parameters ia;
    ia.soft_sphere = SoftSphere_Parameters(/*a*/ 1.0, /*n*/ 2.0, /*cut*/ 2.0,
                                           /*offset*/ 1.0);

    // Case (1): dist == offset => r_off == 0. For any exponent n this makes
    // the unguarded energy a / pow(0, n) == +Inf, while the force is 0.
    auto const dist_at_offset = 1.0;
    BOOST_CHECK_EQUAL(soft_pair_force_factor(ia, dist_at_offset), 0.0);
    BOOST_CHECK(std::isfinite(soft_pair_energy(ia, dist_at_offset)));
    BOOST_CHECK_EQUAL(soft_pair_energy(ia, dist_at_offset), 0.0);

    // Case (1b): dist < offset, r_off < 0, even integer exponent: the
    // unguarded energy is finite but nonzero, still disagreeing with the
    // zero force. The energy must be 0 to match the force.
    auto const dist_below_offset = 0.5;
    BOOST_CHECK_EQUAL(soft_pair_force_factor(ia, dist_below_offset), 0.0);
    BOOST_CHECK(std::isfinite(soft_pair_energy(ia, dist_below_offset)));
    BOOST_CHECK_EQUAL(soft_pair_energy(ia, dist_below_offset), 0.0);
  }

  // Case (2): non-integer exponent n with dist < offset => r_off < 0, so the
  // unguarded energy std::pow(negative, non-integer) == NaN.
  {
    IA_parameters ia;
    ia.soft_sphere = SoftSphere_Parameters(/*a*/ 1.0, /*n*/ 2.5, /*cut*/ 2.0,
                                           /*offset*/ 1.0);
    auto const dist_below_offset = 0.5;
    BOOST_CHECK_EQUAL(soft_pair_force_factor(ia, dist_below_offset), 0.0);
    BOOST_CHECK(std::isfinite(soft_pair_energy(ia, dist_below_offset)));
    BOOST_CHECK_EQUAL(soft_pair_energy(ia, dist_below_offset), 0.0);
  }

  // Sanity: in the genuinely interacting region (r_off > 0, dist < max_cutoff)
  // the energy is the expected positive value and the fix must not change it.
  {
    IA_parameters ia;
    ia.soft_sphere = SoftSphere_Parameters(/*a*/ 1.0, /*n*/ 2.0, /*cut*/ 2.0,
                                           /*offset*/ 1.0);
    auto const dist = 1.5; // r_off == 0.5
    auto const r_off = dist - ia.soft_sphere.offset;
    auto const expected = ia.soft_sphere.a / std::pow(r_off, ia.soft_sphere.n);
    BOOST_CHECK_GT(soft_pair_force_factor(ia, dist), 0.0);
    BOOST_CHECK_CLOSE(soft_pair_energy(ia, dist), expected, 1e-12);

    // beyond the cutoff both force and energy are 0
    BOOST_CHECK_EQUAL(soft_pair_force_factor(ia, 5.0), 0.0);
    BOOST_CHECK_EQUAL(soft_pair_energy(ia, 5.0), 0.0);
  }
}

#else // ESPRESSO_SOFT_SPHERE

int main(int, char **) {}

#endif // ESPRESSO_SOFT_SPHERE
