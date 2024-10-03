/*
 * Copyright (C) 2017-2024 The ESPResSo project
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

#define BOOST_TEST_MODULE "Special mathematical functions"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "p3m/field_layout_helpers.hpp"
#include "utils/Vector.hpp"
#include "utils/index.hpp"

BOOST_AUTO_TEST_CASE(add_remove_halo) {

  int n_x = 3;
  int n_y = 7;
  int n_z = 11;
  Utils::Vector3i shape = {n_x, n_y, n_z};
  std::vector<double> without_halo_orig(n_x * n_y * n_z);
  int n_halo = 3;

  // Fill
  for (int i = 0; i < n_x; i++) {
    for (int j = 0; j < n_y; j++) {
      for (int k = 0; k < n_z; k++) {
        int ind = Utils::get_linear_index(i, j, k, shape,
                                          Utils::MemoryOrder::ROW_MAJOR);
        without_halo_orig[ind] = ind;
      }
    }
  }

  // add halo

  std::vector<double> with_halo =
      pad_with_zeros(without_halo_orig, shape, n_halo);
  Utils::Vector3i shape_with_halo =
      shape + 2 * Utils::Vector3i::broadcast(n_halo);

  BOOST_CHECK_EQUAL(with_halo.size(), shape_with_halo[0] * shape_with_halo[1] *
                                          shape_with_halo[2]);

  // verify inner region
  for (int i = 0; i < n_x; i++) {
    for (int j = 0; j < n_y; j++) {
      for (int k = 0; k < n_z; k++) {
        int ind_without_halo = Utils::get_linear_index(
            i, j, k, shape, Utils::MemoryOrder::ROW_MAJOR);
        int ind_with_halo = Utils::get_linear_index(
            i + n_halo, j + n_halo, k + n_halo, shape_with_halo,
            Utils::MemoryOrder::ROW_MAJOR);
        BOOST_CHECK_EQUAL(with_halo[ind_with_halo],
                          without_halo_orig[ind_without_halo]);
      }
    }
  }

  // verify, added halso are 0
  for (int i = 0; i < shape_with_halo[0]; i++) {
    for (int j = 0; j < shape_with_halo[1]; j++) {
      for (int k = 0; k < shape_with_halo[2]; k++) {
        if ((i < n_halo or i >= shape_with_halo[1] - n_halo) or
            (j < n_halo or j >= shape_with_halo[1] - n_halo) or
            (k < n_halo or k >= shape_with_halo[2] - n_halo)) {
          // this is a halo cell, must be zero
          int ind = Utils::get_linear_index(i, j, k, shape_with_halo,
                                            Utils::MemoryOrder::ROW_MAJOR);
          BOOST_CHECK_EQUAL(with_halo[ind], 0);
        }
      }
    }
  }

  // Remove halo again

  std::vector<double> without_halo_new =
      extract_block(with_halo, shape_with_halo, n_halo);
  for (int i = 0; i < without_halo_orig.size(); i++) {
    BOOST_CHECK_EQUAL(without_halo_new[i], without_halo_orig[i]);
  }
  BOOST_CHECK_EQUAL(without_halo_new.size(), without_halo_orig.size());
}
