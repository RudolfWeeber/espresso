/*
 * Copyright (C) 2020 The ESPResSo project
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

/* Unit tests for the ReactionMethods utility functions. */

#define BOOST_TEST_MODULE ReactionMethods utility functions test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "reaction_methods/utils.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

BOOST_AUTO_TEST_CASE(find_minimum_non_negative_value_test) {
  using namespace ReactionMethods;

  BOOST_CHECK_EQUAL(find_minimum_non_negative_value({1, 2, 3}), 1);
  BOOST_CHECK_EQUAL(find_minimum_non_negative_value({3, 2, 1}), 1);
  BOOST_CHECK_EQUAL(find_minimum_non_negative_value({-1, 2, 3}), 2);
  BOOST_CHECK_EQUAL(find_minimum_non_negative_value({-1, -2, -3}), -3);
  BOOST_CHECK_THROW(find_minimum_non_negative_value({}), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(average_list_of_allowed_entries_test) {
  using namespace ReactionMethods;
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  BOOST_CHECK_CLOSE(average_list_of_allowed_entries(std::vector<long>{1, 2}),
                    1.5, tol);
  BOOST_CHECK_CLOSE(
      average_list_of_allowed_entries(std::vector<double>{1, 2, -2}), 1.5, tol);
  BOOST_CHECK_CLOSE(
      average_list_of_allowed_entries(std::vector<double>{1.5, -3.}), 1.5, tol);
  BOOST_CHECK_CLOSE(
      average_list_of_allowed_entries(std::vector<int>{-1, -2, -3}), 0.0, tol);
  BOOST_CHECK_CLOSE(average_list_of_allowed_entries(std::vector<double>{}), 0.0,
                    tol);
}

BOOST_AUTO_TEST_CASE(factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i_test) {
  using namespace ReactionMethods;
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  auto const reaction_ensemble_combinations = [](int N, int nu) {
    return (N + nu < 0) ? 0. : std::tgamma(N + 1) / std::tgamma(N + nu + 1);
  };

  for (int N0 = 0; N0 < 6; ++N0) {
    for (int nu = -4; nu <= 4; ++nu) {
      auto const val = factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(N0, nu);
      auto const ref = reaction_ensemble_combinations(N0, nu);
      BOOST_CHECK_CLOSE(val, ref, 5 * tol);
    }
  }
}
