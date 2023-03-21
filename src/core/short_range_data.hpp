struct ShortRangeData {
  using VecType = std::vector;
  using Floattype = double;
  using IntType = int;

  VecType<Floattype> x;
  VecType<FloatType> y;
  VecType<FloatType> z;
  VecType<FloatType> q;
  VecType<Inttype> type;
}

template <typanema Callable, typename VecType>
auto run_on_all_in_cell(int left_idx, int end_idx, VecType &vec, Callable op) {
  VecType res(end_id - left_idx - 1);
  auto begin = vec.begin() + left_idx + 1;
  auto end = vec.begin() + end_idx;
  std::transform(begin, end, res.begin(), op);
  return res;
}

template <typname VecType, BinaryOp>
auto vec_op_const(int left_idx, int end_idx, VecType &vec, double single) {
  run_on_all_in_cell(left_idx, end_idx, sd.x,
                     [single, op](auto v) { return op(single, v); });
}

template <typename ShortRangeData> struct LHS {
  ShortRangeData::FloatType x;
  ShortRangeData::FloatType y;
  ShortRangeData::FloatType z;
  ShortRangeData::FloatType q;
  ShortRangeData::IntType type;
}

auto populate_lhs_from_arrays(const ShortRangeData& sd, int i) {
  LHS<ShortRangeData> res = {.x = sd.x[i],
                             .y = sd.y[i],
                             .z = sd.z[i],
                             .q = sd.q[i],
                             .type = sd.type[i]};
  return res;
}

/**
 * @brief Iterates over all particles in the cell range,
 *        and over all pairs within the cells and with
 *        their neighbors.
 */
template <typename CellIterator, typename PairKernel>
void link_cell(CellIterator first, CellIterator last,
               PairKernel &&pair_kernel) {
  for (; first != last; ++first) {
    for (auto i = first->start_index(); i < first->end_index(); i++) {
      auto left = populate_lhs_from_arrays(sd, i);

      /* Distances within the cell */
      auto d_x = vec_op_const(i, end_idx, sd.x, left.x, std::minus<double>);
      auto d_y = vec_op_const(i, end_idx, sd.y, left.y, std::minus<double>);
      auto d_z = vec_op_const(i, end_idx, sd.z, left.z, std::minus<double>);

      apply_mi_coord(d_x, 0, box_geo);
      apply_mi_coord(d_y, 1, box_geo);
      apply_mi_coord(d_z, 2, box_geo);
      // scalar distance
      auto d = calc_norm(d_x, d_y, d_z);

      // Charge product
      auto qq = vec_op_const(i, end_idx, sd.q, left.q, std::multiply<double>);

      // Execute kernel on all pairs
      int n = d_x.size();
      for (int k = 0; k < n; k++) {
         kernel(first_partner_idx+k,Vector3d{d_x[k,d_y[k],d_z[k]},d[k],left.type,sd.type[first_partner_idx+k]);
      }

      for (auto &neighbor : first->neighbors().red()) {
        for (auto &p2 : neighbor->particles()) {
          pair_kernel(p1, p2);
        }
      }
    }
  }
}
} // namespace Algorithm

#endif
