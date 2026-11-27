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

#define BOOST_TEST_MODULE "ScatterView without replica on one thread"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "cell_system/CellStructure.hpp"
#include "cell_system/scatter_view_single_thread.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include <cstddef>
#include <type_traits>

struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

using ForceType = CellStructure::ForceType;
using VirialType = CellStructure::VirialType;
using ScatterForce = CellStructure::ScatterForce;
using ScatterVirial = CellStructure::ScatterVirial;
using execution_space = CellStructure::execution_space;

/**
 * True if the scatter views of @ref CellStructure are duplicated, i.e. the
 * single-thread alias applies. If a Kokkos update changes the ScatterView
 * template signature, @ref ScatterViewSingleThread::make_scatter_view
 * silently falls back to a separate replica and this test fails instead.
 */
static constexpr bool duplicated_expected() {
#ifdef KOKKOS_ENABLE_OPENMP
  return std::is_same_v<execution_space, Kokkos::OpenMP>;
#else
  return false;
#endif
}

static bool single_thread() { return execution_space().concurrency() == 1; }

BOOST_AUTO_TEST_CASE(scatter_view_types) {
  if constexpr (duplicated_expected()) {
    static_assert(
        ScatterViewSingleThread::is_duplicated_layout_right<ScatterForce>{});
    static_assert(
        ScatterViewSingleThread::is_duplicated_layout_right<ScatterVirial>{});
  }
}

/** Scatter contributions to a few values, several per value. */
static void scatter(ScatterForce const &scatter_view, std::size_t n) {
  Kokkos::parallel_for("scatter",
                       Kokkos::RangePolicy<execution_space>(0, 4 * n),
                       [=](std::size_t const k) {
                         auto access = scatter_view.access();
                         auto const i = (k * 7u) % n;
                         for (int d = 0; d < 3; ++d) {
                           access(i, d) += 0.1 * static_cast<double>(k + d);
                         }
                       });
  execution_space{}.fence();
}

BOOST_AUTO_TEST_CASE(force_scatter_view) {
  for (std::size_t const n : {1ul, 3ul, 1000ul}) {
    ForceType dest("dest", n);
    ForceType reference("reference", n);
    auto scatter_dest =
        ScatterViewSingleThread::make_scatter_view<ScatterForce>(dest);
    ScatterForce scatter_reference(reference);

    auto const aliased =
        ScatterViewSingleThread::writes_directly(scatter_dest, dest);
    if constexpr (duplicated_expected()) {
      BOOST_CHECK_EQUAL(aliased, single_thread());
    }
    BOOST_CHECK(not ScatterViewSingleThread::writes_directly(scatter_reference,
                                                             reference) or
                not duplicated_expected());

    // two rounds, with the zeroing of the force calculation in between
    for (int round = 0; round < 2; ++round) {
      Kokkos::deep_copy(dest, 0.);
      Kokkos::deep_copy(reference, 0.);
      scatter_dest.reset_except(dest);
      scatter_reference.reset_except(reference);
      scatter(scatter_dest, n);
      scatter(scatter_reference, n);
      if (aliased) {
        // the contributions are already in the destination
        BOOST_CHECK_NE(dest(0, 2), 0.);
      }
      Kokkos::Experimental::contribute(dest, scatter_dest);
      Kokkos::Experimental::contribute(reference, scatter_reference);
      for (std::size_t i = 0u; i < n; ++i) {
        for (int d = 0; d < 3; ++d) {
          // bitwise identical with one thread
          if (single_thread()) {
            BOOST_REQUIRE_EQUAL(dest(i, d), reference(i, d));
          } else {
            BOOST_REQUIRE_CLOSE(dest(i, d), reference(i, d), 1e-12);
          }
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(virial_scatter_view) {
  VirialType dest("dest");
  auto scatter_dest =
      ScatterViewSingleThread::make_scatter_view<ScatterVirial>(dest);
  if constexpr (duplicated_expected()) {
    BOOST_CHECK_EQUAL(ScatterViewSingleThread::writes_directly(scatter_dest,
                                                               dest),
                      single_thread());
  }
  for (int round = 0; round < 2; ++round) {
    Kokkos::deep_copy(dest, 0.);
    scatter_dest.reset_except(dest);
    Kokkos::parallel_for("scatter_virial",
                         Kokkos::RangePolicy<execution_space>(0, 100),
                         [=](std::size_t const k) {
                           auto access = scatter_dest.access();
                           access(k % 3) += 1.;
                         });
    execution_space{}.fence();
    Kokkos::Experimental::contribute(dest, scatter_dest);
    BOOST_CHECK_EQUAL(dest(0), 34.);
    BOOST_CHECK_EQUAL(dest(1), 33.);
    BOOST_CHECK_EQUAL(dest(2), 33.);
  }
}
