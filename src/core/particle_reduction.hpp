#pragma once

#include <functional>
#include <type_traits>

#include "cell_system/CellStructure.hpp"
#include "config/config.hpp"

namespace Reduction {
template <typename ResultType>
using AddPartialResultKernel =
    std::function<void(const Particle &, ResultType &)>;

template <typename ResultType>
using SingleResultKernel = std::function<ResultType(const Particle &)>;

template <typename ResultType>
using ReductionOp = std::function<ResultType(ResultType &, ResultType &)>;

} // namespace Reduction

template <typename ResultType>
ResultType reduce_over_local_particles(
    const CellStructure &cs,
    Reduction::AddPartialResultKernel<ResultType> &kernel,
    Reduction::ReductionOp<ResultType> &reduce_op) {
  ResultType accumulator{};
  for (const auto &p : cs.local_particles()) {
    ResultType temp{};
    kernel(p, temp);
    accumulator = reduce_op(accumulator, temp);
  }
  return accumulator;
}
