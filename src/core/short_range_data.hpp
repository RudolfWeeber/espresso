#pragma once

#include "BoxGeometry.hpp"
#include "cell_system/Cell.hpp"
#include "cells.hpp"
#include "grid.hpp"
#include "utils/Vector.hpp"
#include <functional>
#include <map>
#include <vector>

struct ShortRangeData {
  template <typename T> using VecType = std::vector<T>;
  using FloatType = double;

  using IntType = int;
  using CellRef = Cell *;
  VecType<FloatType> x;
  VecType<FloatType> y;
  VecType<FloatType> z;
  VecType<FloatType> q;
  VecType<IntType> type;
  std::map<CellRef, size_t> start_idx;
  std::map<CellRef, size_t> end_idx;
};

template <typename VecType>
void apply_mi_convention(VecType &d, int coord, BoxGeometry &box_geo) {
  if (!box_geo.periodic(coord))
    return;
  auto const half_length = box_geo.length_half()[coord];
  auto const length = box_geo.length()[coord];
  std::for_each(d.begin(), d.end(), [half_length, length](auto &dx) {
    while (dx >= half_length)
      dx -= length;
    while (dx < -half_length)
      dx += length;
  });
}

template <typename BinaryOp, typename VecType, typename ScalarType>
auto scalar_vec_op_const(int start_idx, int end_idx, VecType &vec,
                         ScalarType &scalar, BinaryOp op) {
  VecType res(end_idx - start_idx);
  auto begin = vec.begin() + start_idx;
  auto end = vec.begin() + end_idx;
  std::transform(begin, end, res.begin(),
                 [scalar, op](auto a) { return op(scalar, a); });
  return res;
}

template <typename ShortRangeData> struct LHS {
  typename ShortRangeData::FloatType x;
  typename ShortRangeData::FloatType y;
  typename ShortRangeData::FloatType z;
  typename ShortRangeData::FloatType q;
  typename ShortRangeData::IntType type;
};
inline auto populate_lhs_from_arrays(const ShortRangeData &sd, int i) {
  LHS<ShortRangeData> res = {.x = sd.x[i],
                             .y = sd.y[i],
                             .z = sd.z[i],
                             .q = sd.q[i],
                             .type = sd.type[i]};
  return res;
}

template <typename VecType>
auto calc_norm(const VecType &x, const VecType &y, const VecType &z) {
  VecType res(x.size());
  for (int i = 0; i < x.size(); i++) {
    res[i] = std::sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
  }
  return res;
}

template <typename PairKernel>
void run_one_on_many(ShortRangeData &sd, int lhs_idx, int rhs_begin,
                     int rhs_end, PairKernel kernel) {
  printf("rom: %d | %d %d\n", lhs_idx, rhs_begin, rhs_end);
  auto left = populate_lhs_from_arrays(sd, lhs_idx);

  auto d_x = scalar_vec_op_const(rhs_begin, rhs_end, sd.x, left.x,
                                 std::minus<double>());
  auto d_y = scalar_vec_op_const(rhs_begin, rhs_end, sd.y, left.y,
                                 std::minus<double>());
  auto d_z = scalar_vec_op_const(rhs_begin, rhs_end, sd.z, left.z,
                                 std::minus<double>());

  apply_mi_convention(d_x, 0, box_geo);
  apply_mi_convention(d_y, 1, box_geo);
  apply_mi_convention(d_z, 2, box_geo);
  // scalar distance
  auto d = calc_norm(d_x, d_y, d_z);

  // Charge product
  auto qq = scalar_vec_op_const(rhs_begin, rhs_end, sd.q, left.q,
                                std::multiplies<double>());

  // Execute kernel on all pairs
  for (int right_idx = rhs_begin; right_idx < rhs_end; right_idx++) {
    const int data_idx = right_idx - rhs_begin;
    kernel(sd, lhs_idx, right_idx,
           Utils::Vector3d{d_x[data_idx], d_y[data_idx], d_z[data_idx]},
           d[data_idx], qq[data_idx]);
  }
}
template <typename PairKernel>
void run_short_range_kernel(CellStructure &cs, ShortRangeData &sd,
                            PairKernel kernel) {

  for (auto const &cell : cs.decomposition().local_cells()) {
    const int start_idx = sd.start_idx[cell];
    const int end_idx = sd.end_idx[cell];
    if (start_idx == end_idx)
      continue;
    for (int i = start_idx; i < end_idx; i++) {
      run_one_on_many(sd, i, i + 1, end_idx, kernel);
      for (auto const &neighbor_cell : cell->neighbors().red()) {
        const int neighbor_start_idx = sd.start_idx[neighbor_cell];
        const int neighbor_end_idx = sd.end_idx[neighbor_cell];
        if (neighbor_start_idx == neighbor_end_idx)
          continue;
        run_one_on_many(sd, i, neighbor_start_idx, neighbor_end_idx, kernel);
      }
    }
  }
}
void calc_lj_forces();
