#pragma once
#include <functional>
#include <type_traits>

#include "cell_system/CellStructure.hpp"
#include "config/config.hpp"
#ifdef SHARED_MEMORY_PARALLELISM
#include "Kokkos_Core.hpp"
#endif

#include <functional>

namespace Reduction {

/** @brief Kernel that adds the result from a single particle to a reduction */
template <typename ResultType>
using AddPartialResultKernel =
    std::function<void(const Particle &, ResultType &)>;

/** @brief Join two partial reduciton results */
template <typename ResultType>
using ReductionOp = std::function<void(ResultType &, const ResultType &)>;

#ifdef SHARED_MEMORY_PARALLELISM

/** @brief Implements a custom reduction in the form required by  Kokkos */
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
      : reduction_op(reduction_op), kernel(kernel) {}
  KokkosReducer(const KokkosReducer &other)
      : reduction_op(other.reduction_op), kernel(other.kernel) {};

  KOKKOS_INLINE_FUNCTION void operator()(const int i,
                                         value_type &update) const {
    kernel(i, update);
  }

  // "Join" intermediate results from different threads.
  // This should normally implement the same reduction
  // operation as operator() above.
  KOKKOS_INLINE_FUNCTION void join(value_type &dst,
                                   const value_type &src) const {
    reduction_op(dst, src);
  }
};

template <typename ResultType, typename Kernel>
KokkosReducer<ResultType, Kernel>
make_kokkos_reducer(Kernel k, ReductionOp<ResultType> reduce_op) {
  return KokkosReducer<ResultType, Kernel>(k, reduce_op);
}
#endif

} // namespace Reduction

/** @brief performs a reductino over all particles
 *
 * @param add_partial is a function that adds a reduction result from a single
 * particle
 * @param reduce_op is a function that joins two reduction results
 *
 * both functions have to implement the same reduction.
 */
template <typename ResultType>
ResultType reduce_over_local_particles(
    const CellStructure &cs,
    Reduction::AddPartialResultKernel<ResultType> add_partial,
    Reduction::ReductionOp<ResultType> reduce_op) {
#ifdef SHARED_MEMORY_PARALLELISM
  ResultType result{};

  auto const &cells = cs.decomposition().local_cells();
  if (cells.size() > 1) { // parallel loop over cells
    auto reducer = Reduction::make_kokkos_reducer<ResultType>(
        [&cells, add_partial](int i, ResultType &res) {
          for (auto &p : cells[i]->particles()) {
            add_partial(p, res);
          }
        },
        reduce_op);
    Kokkos::parallel_reduce( // loop over cells
        "reduce_on_local_particle", cells.size(), reducer, result);
    return result;
  } else { // cells.size()==1
    auto const &cell = cells[0];
    auto reducer = Reduction::make_kokkos_reducer<ResultType>(
        [&cell, add_partial](int i, ResultType &res) {
          add_partial(*(cell->particles().begin() + i), res);
        },
        reduce_op);
    Kokkos::parallel_reduce( // loop over cells
        "reduce_on_local_particle", cell->particles().size(), reducer, result);
    return result;
  }
#endif

  ResultType accumulated{};
  for (const auto &p : cs.local_particles()) {
    add_partial(p, accumulated);
  }
  return accumulated;
}
