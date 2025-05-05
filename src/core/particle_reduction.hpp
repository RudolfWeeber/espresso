#pragma once

#include <functional>
#include <type_traits>

#include "cell_system/CellStructure.hpp"
#include "config/config.hpp"
template <typename ResultType>
using ParticleReduceKernel =
    std::function<void(const Particle &, ResultType &)>;

template <typename ResultType, typename ReductionOp>
ResultType reduce_over_local_particles(const CellStructure &cs,
                                       ParticleReduceKernel<ResultType> &kernel,
                                       ReductionOp reduce_op) {
  ResultType accumulator{};
  for (const auto &p : cs.local_particles()) {
    ResultType temp{};
    kernel(p, temp);
    accumulator = reduce_op(accumulator, temp);
  }
  return accumulator;
}
