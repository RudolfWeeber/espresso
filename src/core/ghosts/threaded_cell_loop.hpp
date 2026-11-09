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

#pragma once

/**
 * @file
 * Threaded loops over the cells of a ghost communication.
 *
 * The loops write each particle from exactly one thread, in the same order
 * as the serial loops, so the results do not depend on the number of
 * threads.
 */

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <functional>
#include <numeric>
#include <vector>

namespace GhostComm {

/**
 * @brief Minimal number of threads for a threaded loop.
 * With two threads per MPI rank, the fork/join of the parallel regions and
 * the bookkeeping of the threaded loops cost more than they save.
 */
inline constexpr int threaded_loop_min_threads = 4;

/**
 * @brief Minimal number of particles for a threaded loop.
 * Smaller loops run on the calling thread: the fork/join of the parallel
 * region costs more than the loop itself.
 */
inline constexpr std::size_t threaded_loop_min_particles = 256u;

/**
 * @brief Whether threaded loops are used at all.
 * Checked before any bookkeeping for a threaded loop, so that with fewer
 * threads exactly the serial loops run. The size of the thread pool is
 * fixed once Kokkos is initialized; it is queried once, since every query
 * constructs an execution space instance.
 * @pre Kokkos is initialized.
 */
inline bool threads_available() {
  static bool const value = Kokkos::num_threads() >= threaded_loop_min_threads;
  return value;
}

/** @brief Whether a loop over @p n_particles particles runs threaded. */
inline bool use_threaded_loop(std::size_t n_particles) {
  return threads_available() and n_particles >= threaded_loop_min_particles;
}

/**
 * @brief Run @p kernel(i) for i in [0, n) on the host threads.
 * The first exception thrown by a kernel invocation is rethrown on the
 * calling thread once all threads have finished.
 */
template <class Kernel>
void threaded_loop(char const *name, std::size_t n, Kernel const &kernel) {
  std::exception_ptr error;
  std::atomic_flag has_error = ATOMIC_FLAG_INIT;
  Kokkos::parallel_for(
      name,
      Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(std::size_t{0}, n),
      [&](std::size_t i) {
        try {
          kernel(i);
        } catch (...) {
          if (not has_error.test_and_set()) {
            error = std::current_exception();
          }
        }
      });
  Kokkos::DefaultHostExecutionSpace{}.fence();
  if (error) {
    std::rethrow_exception(error);
  }
}

/**
 * @brief Entries of a list grouped by the cell they write to.
 * The entries <tt>order[k]</tt> with <tt>k</tt> in
 * <tt>[begin[g], begin[g + 1])</tt> write to the same cell, in their order
 * in the list. Distinct groups write to distinct cells.
 */
struct CellGroups {
  std::vector<std::size_t> order;
  std::vector<std::size_t> begin;
  std::size_t size() const { return begin.size() - 1u; }
};

/**
 * @brief Group the entries [0, n) of a list by the cell they write to.
 * @param n       Number of entries.
 * @param target  Cell written by entry i.
 */
template <class Target>
CellGroups group_by_target_cell(std::size_t n, Target const &target) {
  CellGroups groups;
  groups.order.resize(n);
  std::iota(groups.order.begin(), groups.order.end(), std::size_t{0});
  std::ranges::stable_sort(groups.order, [&target](auto a, auto b) {
    return std::less<>{}(target(a), target(b));
  });
  for (std::size_t k = 0u; k < n; ++k) {
    if (k == 0u or target(groups.order[k]) != target(groups.order[k - 1u])) {
      groups.begin.push_back(k);
    }
  }
  groups.begin.push_back(n);
  return groups;
}

} // namespace GhostComm
