/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include <config/config.hpp>

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"

#include <Kokkos_Core.hpp>

#include <concepts>
#include <functional>
#include <utility>

namespace Reduction {

/** @brief Kernel that adds the result from a single particle to a reduction */
template <typename ResultType>
using AddPartialResultKernel =
    std::function<void(ResultType &, Particle const &)>;

/** @brief Join two partial reduction results */
template <typename ResultType>
using ReductionOp = std::function<void(ResultType &, ResultType const &)>;

/** @brief Custom reduction in the form required by Kokkos */
template <typename ResultType, typename Kernel> class KokkosReducer {
public:
  // Kokkos reduction functors need the value_type typedef.
  // This is the type of the result of the reduction.
  using value_type = ResultType;

  // Just like with parallel_for functors, you may specify
  // an execution_space typedef. If not provided, Kokkos
  // will use the default execution space by default.

  // kernels to wrap
  ReductionOp<ResultType> reduction_op;
  Kernel kernel;
  KokkosReducer(Kernel kernel, ReductionOp<ResultType> reduction_op)
      : reduction_op(std::move(reduction_op)), kernel(std::move(kernel)) {}
  KokkosReducer(KokkosReducer const &other)
      : reduction_op(other.reduction_op), kernel(other.kernel) {}

  KOKKOS_INLINE_FUNCTION void operator()(std::integral auto const i,
                                         value_type &update) const {
    kernel(i, update);
  }

  // "Join" intermediate results from different threads.
  // This should normally implement the same reduction
  // operation as operator() above.
  KOKKOS_INLINE_FUNCTION void join(value_type &dst,
                                   value_type const &src) const {
    reduction_op(dst, src);
  }
};

template <typename ResultType, typename Kernel>
KokkosReducer<ResultType, Kernel>
make_kokkos_reducer(Kernel k, ReductionOp<ResultType> reduce_op) {
  return KokkosReducer<ResultType, Kernel>(k, reduce_op);
}

} // namespace Reduction

/** @brief performs a reduction over all particles
 *
 * @param cs  cell structure to iterate over
 * @param add_partial is a function that adds a reduction result from a single
 * particle
 * @param reduce_op is a function that joins two reduction results
 *
 * both functions have to implement the same reduction.
 */
template <typename ResultType>
ResultType reduce_over_local_particles(
    CellStructure const &cs,
    Reduction::AddPartialResultKernel<ResultType> add_partial,
    [[maybe_unused]] Reduction::ReductionOp<ResultType> reduce_op) {

  ResultType result{};

  auto const &cells = cs.decomposition().local_cells();
  if (cells.size() > 1) { // parallel loop over cells
    auto reducer = Reduction::make_kokkos_reducer<ResultType>(
        [&cells, add_partial](std::size_t const c_index, ResultType &res) {
          // One reused view per cell (per thread), REBOUND per row via
          // attach_to_store instead of materialising a Particle per cell
          // through the row-range iterator. The cell's committed rows are the
          // contiguous range [offset, offset+count); this runs on a clean
          // store, so index it directly.
          auto *cell = cells[c_index];
          auto const offset = cell->offset();
          auto const n_part = cell->count();
          auto &store = cell->store();
          Particle p;
          for (std::size_t idx = 0u; idx < n_part; ++idx) {
            p.attach_to_store(store, static_cast<int>(offset + idx));
            add_partial(res, std::as_const(p));
          }
        },
        reduce_op);
    Kokkos::parallel_reduce( // loop over cells
        "reduce_on_local_particle", cells.size(), reducer, result);
    return result;
  }
  // single cell case: parallel over particles, each index building its OWN
  // view (a shared cached-view iterator would not be thread-safe). The
  // committed rows are the contiguous range [offset, offset+count).
  auto const offset = cells.front()->offset();
  auto const n_part = cells.front()->count();
  auto &store = cells.front()->store();
  auto reducer = Reduction::make_kokkos_reducer<ResultType>(
      [offset, &store, add_partial](std::size_t const p_index,
                                    ResultType &res) {
        auto const view = store.make_view(static_cast<int>(offset + p_index));
        add_partial(res, std::as_const(view));
      },
      reduce_op);
  Kokkos::parallel_reduce( // loop over particles
      "reduce_on_local_particle", n_part, reducer, result);
  return result;
}
