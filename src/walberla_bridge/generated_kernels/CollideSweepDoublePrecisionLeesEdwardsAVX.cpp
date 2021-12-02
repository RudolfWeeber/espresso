// kernel generated with pystencils v0.4.3+4.g30da657, lbmpy v0.4.3+2.g0e17e61,
// lbmpy_walberla/pystencils_walberla from commit
// 88f85eb7a979f81d68e76009811aeed53ec3014e

//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \\file CollideSweepDoublePrecisionLeesEdwardsAVX.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepDoublePrecisionLeesEdwardsAVX.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#include <immintrin.h>

#define FUNC_PREFIX

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning push
#pragma warning(disable : 1599)
#endif

using namespace std;

namespace walberla {
namespace pystencils {

namespace internal_f11a519921c681cbc9d0b2f51454c920 {
static FUNC_PREFIX void
collidesweepdoubleprecisionleesedwardsavx_collidesweepdoubleprecisionleesedwardsavx(
    double *RESTRICT const _data_force, double *RESTRICT _data_pdfs,
    double *RESTRICT const _data_velocity, int64_t const _size_force_0,
    int64_t const _size_force_1, int64_t const _size_force_2,
    int64_t const _stride_force_1, int64_t const _stride_force_2,
    int64_t const _stride_force_3, int64_t const _stride_pdfs_1,
    int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3,
    int64_t const _stride_velocity_1, int64_t const _stride_velocity_2,
    int64_t const _stride_velocity_3, double omega_bulk, bool points_down,
    bool points_up) {
  const double xi_1 = omega_bulk * -0.5 + 1.0;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    double *RESTRICT _data_velocity_20_30 =
        _data_velocity + _stride_velocity_2 * ctr_2;
    double *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    double *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    double *RESTRICT _data_velocity_20_32 =
        _data_velocity + _stride_velocity_2 * ctr_2 + 2 * _stride_velocity_3;
    double *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    double *RESTRICT _data_velocity_20_31 =
        _data_velocity + _stride_velocity_2 * ctr_2 + _stride_velocity_3;
    double *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      double *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      double *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      double *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      double *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      double *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      double *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      double *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      double *RESTRICT _data_velocity_20_30_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_30;
      double *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      double *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      double *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      double *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      double *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      double *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      double *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      double *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      double *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      double *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      double *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      double *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      double *RESTRICT _data_velocity_20_32_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_32;
      double *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      double *RESTRICT _data_velocity_20_31_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_31;
      double *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      {
        for (int64_t ctr_0 = 0; ctr_0 < (int64_t)((_size_force_0) / (4)) * (4);
             ctr_0 += 4) {
          const __m256d xi_63 = _mm256_load_pd(&_data_pdfs_20_312_10[ctr_0]);
          const __m256d xi_64 = _mm256_load_pd(&_data_pdfs_20_314_10[ctr_0]);
          const __m256d xi_65 = _mm256_load_pd(&_data_pdfs_20_315_10[ctr_0]);
          const __m256d xi_66 = _mm256_load_pd(&_data_pdfs_20_38_10[ctr_0]);
          const __m256d xi_67 = _mm256_load_pd(&_data_pdfs_20_316_10[ctr_0]);
          const __m256d xi_68 = _mm256_load_pd(&_data_pdfs_20_33_10[ctr_0]);
          const __m256d xi_69 = _mm256_load_pd(&_data_pdfs_20_30_10[ctr_0]);
          const __m256d xi_70 = _mm256_load_pd(&_data_pdfs_20_37_10[ctr_0]);
          const __m256d xi_71 = _mm256_load_pd(&_data_velocity_20_30_10[ctr_0]);
          const __m256d xi_72 = _mm256_load_pd(&_data_force_20_30_10[ctr_0]);
          const __m256d xi_73 = _mm256_load_pd(&_data_pdfs_20_35_10[ctr_0]);
          const __m256d xi_74 = _mm256_load_pd(&_data_pdfs_20_32_10[ctr_0]);
          const __m256d xi_75 = _mm256_load_pd(&_data_pdfs_20_318_10[ctr_0]);
          const __m256d xi_76 = _mm256_load_pd(&_data_force_20_32_10[ctr_0]);
          const __m256d xi_77 = _mm256_load_pd(&_data_pdfs_20_317_10[ctr_0]);
          const __m256d xi_78 = _mm256_load_pd(&_data_pdfs_20_39_10[ctr_0]);
          const __m256d xi_79 = _mm256_load_pd(&_data_pdfs_20_313_10[ctr_0]);
          const __m256d xi_80 = _mm256_load_pd(&_data_force_20_31_10[ctr_0]);
          const __m256d xi_81 = _mm256_load_pd(&_data_pdfs_20_310_10[ctr_0]);
          const __m256d xi_82 = _mm256_load_pd(&_data_pdfs_20_311_10[ctr_0]);
          const __m256d xi_83 = _mm256_load_pd(&_data_pdfs_20_31_10[ctr_0]);
          const __m256d xi_84 = _mm256_load_pd(&_data_velocity_20_32_10[ctr_0]);
          const __m256d xi_85 = _mm256_load_pd(&_data_pdfs_20_36_10[ctr_0]);
          const __m256d xi_86 = _mm256_load_pd(&_data_velocity_20_31_10[ctr_0]);
          const __m256d xi_87 = _mm256_load_pd(&_data_pdfs_20_34_10[ctr_0]);
          const __m256d xi_3 = _mm256_mul_pd(xi_80, xi_86);
          const __m256d xi_4 = _mm256_mul_pd(xi_76, xi_84);
          const __m256d xi_5 =
              _mm256_mul_pd(xi_86, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
          const __m256d xi_6 =
              _mm256_add_pd(xi_5, _mm256_set_pd(1.0, 1.0, 1.0, 1.0));
          const __m256d xi_7 = _mm256_mul_pd(
              xi_80, _mm256_set_pd(0.166666666666667, 0.166666666666667,
                                   0.166666666666667, 0.166666666666667));
          const __m256d xi_9 = _mm256_mul_pd(
              xi_4, _mm256_set_pd(-0.166666666666667, -0.166666666666667,
                                  -0.166666666666667, -0.166666666666667));
          const __m256d xi_11 =
              _mm256_add_pd(xi_5, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_14 = _mm256_mul_pd(
              xi_72, _mm256_set_pd(0.166666666666667, 0.166666666666667,
                                   0.166666666666667, 0.166666666666667));
          const __m256d xi_15 = _mm256_mul_pd(
              xi_3, _mm256_set_pd(-0.166666666666667, -0.166666666666667,
                                  -0.166666666666667, -0.166666666666667));
          const __m256d xi_16 = _mm256_add_pd(xi_15, xi_9);
          const __m256d xi_18 =
              _mm256_mul_pd(xi_84, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
          const __m256d xi_19 =
              _mm256_add_pd(xi_18, _mm256_set_pd(1.0, 1.0, 1.0, 1.0));
          const __m256d xi_20 = _mm256_mul_pd(
              xi_76, _mm256_set_pd(0.166666666666667, 0.166666666666667,
                                   0.166666666666667, 0.166666666666667));
          const __m256d xi_22 =
              _mm256_add_pd(xi_18, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_23 = _mm256_mul_pd(
              xi_4, _mm256_set_pd(-0.0833333333333333, -0.0833333333333333,
                                  -0.0833333333333333, -0.0833333333333333));
          const __m256d xi_24 =
              _mm256_mul_pd(xi_86, _mm256_set_pd(3.0, 3.0, 3.0, 3.0));
          const __m256d xi_26 = _mm256_mul_pd(
              xi_72, _mm256_set_pd(0.0833333333333333, 0.0833333333333333,
                                   0.0833333333333333, 0.0833333333333333));
          const __m256d xi_29 = _mm256_mul_pd(
              xi_80, _mm256_set_pd(0.0833333333333333, 0.0833333333333333,
                                   0.0833333333333333, 0.0833333333333333));
          const __m256d xi_30 =
              _mm256_mul_pd(xi_24, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_31 = _mm256_add_pd(
              _mm256_mul_pd(xi_5, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
              _mm256_set_pd(1.0, 1.0, 1.0, 1.0));
          const __m256d xi_33 =
              _mm256_mul_pd(xi_84, _mm256_set_pd(3.0, 3.0, 3.0, 3.0));
          const __m256d xi_34 = _mm256_mul_pd(
              xi_76, _mm256_set_pd(0.0833333333333333, 0.0833333333333333,
                                   0.0833333333333333, 0.0833333333333333));
          const __m256d xi_35 = _mm256_mul_pd(
              xi_3, _mm256_set_pd(-0.0833333333333333, -0.0833333333333333,
                                  -0.0833333333333333, -0.0833333333333333));
          const __m256d xi_36 =
              _mm256_mul_pd(xi_33, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_37 = _mm256_add_pd(
              _mm256_mul_pd(xi_18, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
              _mm256_set_pd(1.0, 1.0, 1.0, 1.0));
          const __m256d xi_38 =
              _mm256_mul_pd(xi_84, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_39 = (_mm256_mul_pd(xi_86, xi_86));
          const __m256d xi_41 = (_mm256_mul_pd(xi_84, xi_84));
          const __m256d xi_44 = _mm256_mul_pd(
              xi_86, _mm256_set_pd(0.166666666666667, 0.166666666666667,
                                   0.166666666666667, 0.166666666666667));
          const __m256d xi_45 =
              _mm256_mul_pd(xi_39, _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
          const __m256d xi_48 = _mm256_mul_pd(
              xi_84, _mm256_set_pd(0.166666666666667, 0.166666666666667,
                                   0.166666666666667, 0.166666666666667));
          const __m256d xi_49 =
              _mm256_mul_pd(xi_41, _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
          const __m256d rho = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_add_pd(
                                      _mm256_add_pd(
                                          _mm256_add_pd(
                                              _mm256_add_pd(
                                                  _mm256_add_pd(
                                                      _mm256_add_pd(
                                                          _mm256_add_pd(
                                                              _mm256_add_pd(
                                                                  _mm256_add_pd(
                                                                      _mm256_add_pd(
                                                                          _mm256_add_pd(
                                                                              _mm256_add_pd(
                                                                                  xi_63,
                                                                                  xi_64),
                                                                              xi_65),
                                                                          xi_66),
                                                                      xi_67),
                                                                  xi_68),
                                                              xi_69),
                                                          xi_70),
                                                      xi_73),
                                                  xi_74),
                                              xi_75),
                                          xi_77),
                                      xi_78),
                                  xi_79),
                              xi_81),
                          xi_82),
                      xi_83),
                  xi_85),
              xi_87);
          const __m256d xi_40 =
              _mm256_mul_pd(rho, _mm256_set_pd(1.5, 1.5, 1.5, 1.5));
          const __m256d u_0 = _mm256_add_pd(
              xi_71, _mm256_set_pd(
                         ((points_down && ctr_1 <= 0)
                              ? (1.0)
                              : ((points_up && ctr_1 >= 63) ? (-1.0) : (0.0))) *
                             0.050000000000000003,
                         ((points_down && ctr_1 <= 0)
                              ? (1.0)
                              : ((points_up && ctr_1 >= 63) ? (-1.0) : (0.0))) *
                             0.050000000000000003,
                         ((points_down && ctr_1 <= 0)
                              ? (1.0)
                              : ((points_up && ctr_1 >= 63) ? (-1.0) : (0.0))) *
                             0.050000000000000003,
                         ((points_down && ctr_1 <= 0)
                              ? (1.0)
                              : ((points_up && ctr_1 >= 63) ? (-1.0) : (0.0))) *
                             0.050000000000000003));
          const __m256d xi_2 = _mm256_mul_pd(u_0, xi_72);
          const __m256d xi_8 = _mm256_mul_pd(
              xi_2, _mm256_set_pd(-0.166666666666667, -0.166666666666667,
                                  -0.166666666666667, -0.166666666666667));
          const __m256d xi_10 = _mm256_add_pd(xi_8, xi_9);
          const __m256d xi_12 =
              _mm256_mul_pd(u_0, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
          const __m256d xi_13 =
              _mm256_add_pd(xi_12, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_17 =
              _mm256_add_pd(xi_12, _mm256_set_pd(1.0, 1.0, 1.0, 1.0));
          const __m256d xi_21 = _mm256_add_pd(xi_15, xi_8);
          const __m256d xi_25 = _mm256_add_pd(
              _mm256_mul_pd(xi_12, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
              _mm256_set_pd(1.0, 1.0, 1.0, 1.0));
          const __m256d xi_27 =
              _mm256_mul_pd(u_0, _mm256_set_pd(3.0, 3.0, 3.0, 3.0));
          const __m256d xi_28 =
              _mm256_mul_pd(xi_27, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_32 = _mm256_mul_pd(
              xi_2, _mm256_set_pd(-0.0833333333333333, -0.0833333333333333,
                                  -0.0833333333333333, -0.0833333333333333));
          const __m256d xi_42 = (_mm256_mul_pd(u_0, u_0));
          const __m256d xi_46 = _mm256_mul_pd(
              u_0, _mm256_set_pd(0.166666666666667, 0.166666666666667,
                                 0.166666666666667, 0.166666666666667));
          const __m256d xi_47 =
              _mm256_mul_pd(xi_42, _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
          const __m256d forceTerm_0 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(_mm256_mul_pd(xi_2, _mm256_set_pd(-1.0, -1.0,
                                                                  -1.0, -1.0)),
                                _mm256_mul_pd(xi_3, _mm256_set_pd(-1.0, -1.0,
                                                                  -1.0, -1.0))),
                  _mm256_mul_pd(xi_4, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_1 =
              _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(xi_6, xi_7), xi_10),
                            _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_2 =
              _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(xi_11, xi_7), xi_10),
                            _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_3 =
              _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(xi_13, xi_14), xi_16),
                            _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_4 =
              _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(xi_14, xi_17), xi_16),
                            _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_5 =
              _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(xi_19, xi_20), xi_21),
                            _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_6 =
              _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(xi_20, xi_22), xi_21),
                            _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_7 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_29, _mm256_add_pd(xi_28, xi_6)),
                      _mm256_mul_pd(
                          _mm256_mul_pd(xi_26, _mm256_add_pd(xi_24, xi_25)),
                          _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  xi_23),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_8 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_26, _mm256_add_pd(xi_17, xi_24)),
                      _mm256_mul_pd(xi_29, _mm256_add_pd(xi_27, xi_6))),
                  xi_23),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_9 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_29, _mm256_add_pd(xi_11, xi_27)),
                      _mm256_mul_pd(xi_26, _mm256_add_pd(xi_13, xi_24))),
                  xi_23),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_10 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_26, _mm256_add_pd(xi_17, xi_30)),
                      _mm256_mul_pd(
                          _mm256_mul_pd(xi_29, _mm256_add_pd(xi_27, xi_31)),
                          _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  xi_23),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_11 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_34, _mm256_add_pd(xi_19, xi_24)),
                      _mm256_mul_pd(xi_29, _mm256_add_pd(xi_33, xi_6))),
                  xi_32),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_12 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_34, _mm256_add_pd(xi_19, xi_30)),
                      _mm256_mul_pd(
                          _mm256_mul_pd(xi_29, _mm256_add_pd(xi_31, xi_33)),
                          _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  xi_32),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_13 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_34, _mm256_add_pd(xi_19, xi_28)),
                      _mm256_mul_pd(
                          _mm256_mul_pd(xi_26, _mm256_add_pd(xi_25, xi_33)),
                          _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  xi_35),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_14 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_26, _mm256_add_pd(xi_17, xi_33)),
                      _mm256_mul_pd(xi_34, _mm256_add_pd(xi_19, xi_27))),
                  xi_35),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_15 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_29, _mm256_add_pd(xi_36, xi_6)),
                      _mm256_mul_pd(
                          _mm256_mul_pd(xi_34, _mm256_add_pd(xi_24, xi_37)),
                          _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  xi_32),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_16 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_29, _mm256_add_pd(xi_11, xi_33)),
                      _mm256_mul_pd(xi_34, _mm256_add_pd(xi_22, xi_24))),
                  xi_32),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_17 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_26, _mm256_add_pd(xi_13, xi_33)),
                      _mm256_mul_pd(xi_34, _mm256_add_pd(xi_22, xi_27))),
                  xi_35),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d forceTerm_18 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_26, _mm256_add_pd(xi_17, xi_36)),
                      _mm256_mul_pd(
                          _mm256_mul_pd(xi_34, _mm256_add_pd(xi_27, xi_37)),
                          _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  xi_35),
              _mm256_set_pd(xi_1, xi_1, xi_1, xi_1));
          const __m256d u0Mu1 = _mm256_add_pd(
              _mm256_mul_pd(xi_86, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), u_0);
          const __m256d xi_51 = _mm256_mul_pd(
              u0Mu1, _mm256_set_pd(0.0833333333333333, 0.0833333333333333,
                                   0.0833333333333333, 0.0833333333333333));
          const __m256d xi_52 =
              _mm256_mul_pd((_mm256_mul_pd(u0Mu1, u0Mu1)),
                            _mm256_set_pd(0.125, 0.125, 0.125, 0.125));
          const __m256d u0Pu1 = _mm256_add_pd(u_0, xi_86);
          const __m256d xi_53 = _mm256_mul_pd(
              u0Pu1, _mm256_set_pd(0.0833333333333333, 0.0833333333333333,
                                   0.0833333333333333, 0.0833333333333333));
          const __m256d xi_54 =
              _mm256_mul_pd((_mm256_mul_pd(u0Pu1, u0Pu1)),
                            _mm256_set_pd(0.125, 0.125, 0.125, 0.125));
          const __m256d u1Pu2 = _mm256_add_pd(xi_84, xi_86);
          const __m256d xi_55 = _mm256_mul_pd(
              u1Pu2, _mm256_set_pd(0.0833333333333333, 0.0833333333333333,
                                   0.0833333333333333, 0.0833333333333333));
          const __m256d xi_56 =
              _mm256_mul_pd((_mm256_mul_pd(u1Pu2, u1Pu2)),
                            _mm256_set_pd(0.125, 0.125, 0.125, 0.125));
          const __m256d u1Mu2 = _mm256_add_pd(xi_38, xi_86);
          const __m256d xi_57 = _mm256_mul_pd(
              u1Mu2, _mm256_set_pd(0.0833333333333333, 0.0833333333333333,
                                   0.0833333333333333, 0.0833333333333333));
          const __m256d xi_58 =
              _mm256_mul_pd((_mm256_mul_pd(u1Mu2, u1Mu2)),
                            _mm256_set_pd(0.125, 0.125, 0.125, 0.125));
          const __m256d u0Mu2 = _mm256_add_pd(u_0, xi_38);
          const __m256d xi_59 = _mm256_mul_pd(
              u0Mu2, _mm256_set_pd(0.0833333333333333, 0.0833333333333333,
                                   0.0833333333333333, 0.0833333333333333));
          const __m256d xi_60 =
              _mm256_mul_pd((_mm256_mul_pd(u0Mu2, u0Mu2)),
                            _mm256_set_pd(0.125, 0.125, 0.125, 0.125));
          const __m256d u0Pu2 = _mm256_add_pd(u_0, xi_84);
          const __m256d xi_61 = _mm256_mul_pd(
              u0Pu2, _mm256_set_pd(0.0833333333333333, 0.0833333333333333,
                                   0.0833333333333333, 0.0833333333333333));
          const __m256d xi_62 =
              _mm256_mul_pd((_mm256_mul_pd(u0Pu2, u0Pu2)),
                            _mm256_set_pd(0.125, 0.125, 0.125, 0.125));
          const __m256d f_eq_common = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(_mm256_mul_pd(xi_39, xi_40),
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                      _mm256_mul_pd(_mm256_mul_pd(xi_40, xi_41),
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  _mm256_mul_pd(_mm256_mul_pd(xi_40, xi_42),
                                _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
              rho);
          const __m256d xi_43 = _mm256_mul_pd(
              f_eq_common,
              _mm256_set_pd(0.0555555555555556, 0.0555555555555556,
                            0.0555555555555556, 0.0555555555555556));
          const __m256d xi_50 = _mm256_mul_pd(
              f_eq_common,
              _mm256_set_pd(0.0277777777777778, 0.0277777777777778,
                            0.0277777777777778, 0.0277777777777778));
          _mm256_store_pd(
              &_data_pdfs_20_30_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_mul_pd(f_eq_common,
                                            _mm256_set_pd(0.333333333333333,
                                                          0.333333333333333,
                                                          0.333333333333333,
                                                          0.333333333333333)),
                              _mm256_mul_pd(xi_69, _mm256_set_pd(-1.0, -1.0,
                                                                 -1.0, -1.0))),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_0),
                  xi_69));
          _mm256_store_pd(
              &_data_pdfs_20_31_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(rho,
                                                _mm256_add_pd(xi_44, xi_45)),
                                  _mm256_mul_pd(
                                      xi_83,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_43),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_1),
                  xi_83));
          _mm256_store_pd(
              &_data_pdfs_20_32_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      rho,
                                      _mm256_add_pd(
                                          _mm256_mul_pd(
                                              xi_44, _mm256_set_pd(-1.0, -1.0,
                                                                   -1.0, -1.0)),
                                          xi_45)),
                                  _mm256_mul_pd(
                                      xi_74,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_43),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_2),
                  xi_74));
          _mm256_store_pd(
              &_data_pdfs_20_33_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      rho,
                                      _mm256_add_pd(
                                          _mm256_mul_pd(
                                              xi_46, _mm256_set_pd(-1.0, -1.0,
                                                                   -1.0, -1.0)),
                                          xi_47)),
                                  _mm256_mul_pd(
                                      xi_68,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_43),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_3),
                  xi_68));
          _mm256_store_pd(
              &_data_pdfs_20_34_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(rho,
                                                _mm256_add_pd(xi_46, xi_47)),
                                  _mm256_mul_pd(
                                      xi_87,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_43),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_4),
                  xi_87));
          _mm256_store_pd(
              &_data_pdfs_20_35_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(rho,
                                                _mm256_add_pd(xi_48, xi_49)),
                                  _mm256_mul_pd(
                                      xi_73,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_43),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_5),
                  xi_73));
          _mm256_store_pd(
              &_data_pdfs_20_36_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      rho,
                                      _mm256_add_pd(
                                          _mm256_mul_pd(
                                              xi_48, _mm256_set_pd(-1.0, -1.0,
                                                                   -1.0, -1.0)),
                                          xi_49)),
                                  _mm256_mul_pd(
                                      xi_85,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_43),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_6),
                  xi_85));
          _mm256_store_pd(
              &_data_pdfs_20_37_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      rho,
                                      _mm256_add_pd(
                                          _mm256_mul_pd(
                                              xi_51, _mm256_set_pd(-1.0, -1.0,
                                                                   -1.0, -1.0)),
                                          xi_52)),
                                  _mm256_mul_pd(
                                      xi_70,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_50),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_7),
                  xi_70));
          _mm256_store_pd(
              &_data_pdfs_20_38_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(rho,
                                                _mm256_add_pd(xi_53, xi_54)),
                                  _mm256_mul_pd(
                                      xi_66,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_50),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_8),
                  xi_66));
          _mm256_store_pd(
              &_data_pdfs_20_39_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      rho,
                                      _mm256_add_pd(
                                          _mm256_mul_pd(
                                              xi_53, _mm256_set_pd(-1.0, -1.0,
                                                                   -1.0, -1.0)),
                                          xi_54)),
                                  _mm256_mul_pd(
                                      xi_78,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_50),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_9),
                  xi_78));
          _mm256_store_pd(
              &_data_pdfs_20_310_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(rho,
                                                _mm256_add_pd(xi_51, xi_52)),
                                  _mm256_mul_pd(
                                      xi_81,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_50),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_10),
                  xi_81));
          _mm256_store_pd(
              &_data_pdfs_20_311_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(rho,
                                                _mm256_add_pd(xi_55, xi_56)),
                                  _mm256_mul_pd(
                                      xi_82,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_50),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_11),
                  xi_82));
          _mm256_store_pd(
              &_data_pdfs_20_312_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      rho,
                                      _mm256_add_pd(
                                          _mm256_mul_pd(
                                              xi_57, _mm256_set_pd(-1.0, -1.0,
                                                                   -1.0, -1.0)),
                                          xi_58)),
                                  _mm256_mul_pd(
                                      xi_63,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_50),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_12),
                  xi_63));
          _mm256_store_pd(
              &_data_pdfs_20_313_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      rho,
                                      _mm256_add_pd(
                                          _mm256_mul_pd(
                                              xi_59, _mm256_set_pd(-1.0, -1.0,
                                                                   -1.0, -1.0)),
                                          xi_60)),
                                  _mm256_mul_pd(
                                      xi_79,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_50),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_13),
                  xi_79));
          _mm256_store_pd(
              &_data_pdfs_20_314_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(rho,
                                                _mm256_add_pd(xi_61, xi_62)),
                                  _mm256_mul_pd(
                                      xi_64,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_50),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_14),
                  xi_64));
          _mm256_store_pd(
              &_data_pdfs_20_315_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(rho,
                                                _mm256_add_pd(xi_57, xi_58)),
                                  _mm256_mul_pd(
                                      xi_65,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_50),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_15),
                  xi_65));
          _mm256_store_pd(
              &_data_pdfs_20_316_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      rho,
                                      _mm256_add_pd(
                                          _mm256_mul_pd(
                                              xi_55, _mm256_set_pd(-1.0, -1.0,
                                                                   -1.0, -1.0)),
                                          xi_56)),
                                  _mm256_mul_pd(
                                      xi_67,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_50),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_16),
                  xi_67));
          _mm256_store_pd(
              &_data_pdfs_20_317_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      rho,
                                      _mm256_add_pd(
                                          _mm256_mul_pd(
                                              xi_61, _mm256_set_pd(-1.0, -1.0,
                                                                   -1.0, -1.0)),
                                          xi_62)),
                                  _mm256_mul_pd(
                                      xi_77,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_50),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_17),
                  xi_77));
          _mm256_store_pd(
              &_data_pdfs_20_318_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(rho,
                                                _mm256_add_pd(xi_59, xi_60)),
                                  _mm256_mul_pd(
                                      xi_75,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_50),
                          _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk)),
                      forceTerm_18),
                  xi_75));
        }
        for (int64_t ctr_0 = (int64_t)((_size_force_0) / (4)) * (4);
             ctr_0 < _size_force_0; ctr_0 += 1) {
          const double xi_63 = _data_pdfs_20_312_10[ctr_0];
          const double xi_64 = _data_pdfs_20_314_10[ctr_0];
          const double xi_65 = _data_pdfs_20_315_10[ctr_0];
          const double xi_66 = _data_pdfs_20_38_10[ctr_0];
          const double xi_67 = _data_pdfs_20_316_10[ctr_0];
          const double xi_68 = _data_pdfs_20_33_10[ctr_0];
          const double xi_69 = _data_pdfs_20_30_10[ctr_0];
          const double xi_70 = _data_pdfs_20_37_10[ctr_0];
          const double xi_71 = _data_velocity_20_30_10[ctr_0];
          const double xi_72 = _data_force_20_30_10[ctr_0];
          const double xi_73 = _data_pdfs_20_35_10[ctr_0];
          const double xi_74 = _data_pdfs_20_32_10[ctr_0];
          const double xi_75 = _data_pdfs_20_318_10[ctr_0];
          const double xi_76 = _data_force_20_32_10[ctr_0];
          const double xi_77 = _data_pdfs_20_317_10[ctr_0];
          const double xi_78 = _data_pdfs_20_39_10[ctr_0];
          const double xi_79 = _data_pdfs_20_313_10[ctr_0];
          const double xi_80 = _data_force_20_31_10[ctr_0];
          const double xi_81 = _data_pdfs_20_310_10[ctr_0];
          const double xi_82 = _data_pdfs_20_311_10[ctr_0];
          const double xi_83 = _data_pdfs_20_31_10[ctr_0];
          const double xi_84 = _data_velocity_20_32_10[ctr_0];
          const double xi_85 = _data_pdfs_20_36_10[ctr_0];
          const double xi_86 = _data_velocity_20_31_10[ctr_0];
          const double xi_87 = _data_pdfs_20_34_10[ctr_0];
          const double xi_3 = xi_80 * xi_86;
          const double xi_4 = xi_76 * xi_84;
          const double xi_5 = xi_86 * 2.0;
          const double xi_6 = xi_5 + 1.0;
          const double xi_7 = xi_80 * 0.166666666666667;
          const double xi_9 = xi_4 * -0.166666666666667;
          const double xi_11 = xi_5 - 1.0;
          const double xi_14 = xi_72 * 0.166666666666667;
          const double xi_15 = xi_3 * -0.166666666666667;
          const double xi_16 = xi_15 + xi_9;
          const double xi_18 = xi_84 * 2.0;
          const double xi_19 = xi_18 + 1.0;
          const double xi_20 = xi_76 * 0.166666666666667;
          const double xi_22 = xi_18 - 1.0;
          const double xi_23 = xi_4 * -0.0833333333333333;
          const double xi_24 = xi_86 * 3.0;
          const double xi_26 = xi_72 * 0.0833333333333333;
          const double xi_29 = xi_80 * 0.0833333333333333;
          const double xi_30 = -xi_24;
          const double xi_31 = -xi_5 + 1.0;
          const double xi_33 = xi_84 * 3.0;
          const double xi_34 = xi_76 * 0.0833333333333333;
          const double xi_35 = xi_3 * -0.0833333333333333;
          const double xi_36 = -xi_33;
          const double xi_37 = -xi_18 + 1.0;
          const double xi_38 = -xi_84;
          const double xi_39 = (xi_86 * xi_86);
          const double xi_41 = (xi_84 * xi_84);
          const double xi_44 = xi_86 * 0.166666666666667;
          const double xi_45 = xi_39 * 0.25;
          const double xi_48 = xi_84 * 0.166666666666667;
          const double xi_49 = xi_41 * 0.25;
          const double rho = xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68 +
                             xi_69 + xi_70 + xi_73 + xi_74 + xi_75 + xi_77 +
                             xi_78 + xi_79 + xi_81 + xi_82 + xi_83 + xi_85 +
                             xi_87;
          const double xi_40 = rho * 1.5;
          const double u_0 =
              xi_71 + ((points_down && ctr_1 <= 0)
                           ? (1.0)
                           : ((points_up && ctr_1 >= 63) ? (-1.0) : (0.0))) *
                          0.050000000000000003;
          const double xi_2 = u_0 * xi_72;
          const double xi_8 = xi_2 * -0.166666666666667;
          const double xi_10 = xi_8 + xi_9;
          const double xi_12 = u_0 * 2.0;
          const double xi_13 = xi_12 - 1.0;
          const double xi_17 = xi_12 + 1.0;
          const double xi_21 = xi_15 + xi_8;
          const double xi_25 = -xi_12 + 1.0;
          const double xi_27 = u_0 * 3.0;
          const double xi_28 = -xi_27;
          const double xi_32 = xi_2 * -0.0833333333333333;
          const double xi_42 = (u_0 * u_0);
          const double xi_46 = u_0 * 0.166666666666667;
          const double xi_47 = xi_42 * 0.25;
          const double forceTerm_0 = xi_1 * (-xi_2 - xi_3 - xi_4);
          const double forceTerm_1 = xi_1 * (xi_10 + xi_6 * xi_7);
          const double forceTerm_2 = xi_1 * (xi_10 + xi_11 * xi_7);
          const double forceTerm_3 = xi_1 * (xi_13 * xi_14 + xi_16);
          const double forceTerm_4 = xi_1 * (xi_14 * xi_17 + xi_16);
          const double forceTerm_5 = xi_1 * (xi_19 * xi_20 + xi_21);
          const double forceTerm_6 = xi_1 * (xi_20 * xi_22 + xi_21);
          const double forceTerm_7 =
              xi_1 * (xi_23 - xi_26 * (xi_24 + xi_25) + xi_29 * (xi_28 + xi_6));
          const double forceTerm_8 =
              xi_1 * (xi_23 + xi_26 * (xi_17 + xi_24) + xi_29 * (xi_27 + xi_6));
          const double forceTerm_9 = xi_1 * (xi_23 + xi_26 * (xi_13 + xi_24) +
                                             xi_29 * (xi_11 + xi_27));
          const double forceTerm_10 = xi_1 * (xi_23 + xi_26 * (xi_17 + xi_30) -
                                              xi_29 * (xi_27 + xi_31));
          const double forceTerm_11 =
              xi_1 * (xi_29 * (xi_33 + xi_6) + xi_32 + xi_34 * (xi_19 + xi_24));
          const double forceTerm_12 = xi_1 * (-xi_29 * (xi_31 + xi_33) + xi_32 +
                                              xi_34 * (xi_19 + xi_30));
          const double forceTerm_13 = xi_1 * (-xi_26 * (xi_25 + xi_33) +
                                              xi_34 * (xi_19 + xi_28) + xi_35);
          const double forceTerm_14 = xi_1 * (xi_26 * (xi_17 + xi_33) +
                                              xi_34 * (xi_19 + xi_27) + xi_35);
          const double forceTerm_15 =
              xi_1 * (xi_29 * (xi_36 + xi_6) + xi_32 - xi_34 * (xi_24 + xi_37));
          const double forceTerm_16 = xi_1 * (xi_29 * (xi_11 + xi_33) + xi_32 +
                                              xi_34 * (xi_22 + xi_24));
          const double forceTerm_17 = xi_1 * (xi_26 * (xi_13 + xi_33) +
                                              xi_34 * (xi_22 + xi_27) + xi_35);
          const double forceTerm_18 = xi_1 * (xi_26 * (xi_17 + xi_36) -
                                              xi_34 * (xi_27 + xi_37) + xi_35);
          const double u0Mu1 = u_0 - xi_86;
          const double xi_51 = u0Mu1 * 0.0833333333333333;
          const double xi_52 = (u0Mu1 * u0Mu1) * 0.125;
          const double u0Pu1 = u_0 + xi_86;
          const double xi_53 = u0Pu1 * 0.0833333333333333;
          const double xi_54 = (u0Pu1 * u0Pu1) * 0.125;
          const double u1Pu2 = xi_84 + xi_86;
          const double xi_55 = u1Pu2 * 0.0833333333333333;
          const double xi_56 = (u1Pu2 * u1Pu2) * 0.125;
          const double u1Mu2 = xi_38 + xi_86;
          const double xi_57 = u1Mu2 * 0.0833333333333333;
          const double xi_58 = (u1Mu2 * u1Mu2) * 0.125;
          const double u0Mu2 = u_0 + xi_38;
          const double xi_59 = u0Mu2 * 0.0833333333333333;
          const double xi_60 = (u0Mu2 * u0Mu2) * 0.125;
          const double u0Pu2 = u_0 + xi_84;
          const double xi_61 = u0Pu2 * 0.0833333333333333;
          const double xi_62 = (u0Pu2 * u0Pu2) * 0.125;
          const double f_eq_common =
              rho - xi_39 * xi_40 - xi_40 * xi_41 - xi_40 * xi_42;
          const double xi_43 = f_eq_common * 0.0555555555555556;
          const double xi_50 = f_eq_common * 0.0277777777777778;
          _data_pdfs_20_30_10[ctr_0] =
              forceTerm_0 +
              omega_bulk * (f_eq_common * 0.333333333333333 - xi_69) + xi_69;
          _data_pdfs_20_31_10[ctr_0] =
              forceTerm_1 +
              omega_bulk * (rho * (xi_44 + xi_45) + xi_43 - xi_83) + xi_83;
          _data_pdfs_20_32_10[ctr_0] =
              forceTerm_2 +
              omega_bulk * (rho * (-xi_44 + xi_45) + xi_43 - xi_74) + xi_74;
          _data_pdfs_20_33_10[ctr_0] =
              forceTerm_3 +
              omega_bulk * (rho * (-xi_46 + xi_47) + xi_43 - xi_68) + xi_68;
          _data_pdfs_20_34_10[ctr_0] =
              forceTerm_4 +
              omega_bulk * (rho * (xi_46 + xi_47) + xi_43 - xi_87) + xi_87;
          _data_pdfs_20_35_10[ctr_0] =
              forceTerm_5 +
              omega_bulk * (rho * (xi_48 + xi_49) + xi_43 - xi_73) + xi_73;
          _data_pdfs_20_36_10[ctr_0] =
              forceTerm_6 +
              omega_bulk * (rho * (-xi_48 + xi_49) + xi_43 - xi_85) + xi_85;
          _data_pdfs_20_37_10[ctr_0] =
              forceTerm_7 +
              omega_bulk * (rho * (-xi_51 + xi_52) + xi_50 - xi_70) + xi_70;
          _data_pdfs_20_38_10[ctr_0] =
              forceTerm_8 +
              omega_bulk * (rho * (xi_53 + xi_54) + xi_50 - xi_66) + xi_66;
          _data_pdfs_20_39_10[ctr_0] =
              forceTerm_9 +
              omega_bulk * (rho * (-xi_53 + xi_54) + xi_50 - xi_78) + xi_78;
          _data_pdfs_20_310_10[ctr_0] =
              forceTerm_10 +
              omega_bulk * (rho * (xi_51 + xi_52) + xi_50 - xi_81) + xi_81;
          _data_pdfs_20_311_10[ctr_0] =
              forceTerm_11 +
              omega_bulk * (rho * (xi_55 + xi_56) + xi_50 - xi_82) + xi_82;
          _data_pdfs_20_312_10[ctr_0] =
              forceTerm_12 +
              omega_bulk * (rho * (-xi_57 + xi_58) + xi_50 - xi_63) + xi_63;
          _data_pdfs_20_313_10[ctr_0] =
              forceTerm_13 +
              omega_bulk * (rho * (-xi_59 + xi_60) + xi_50 - xi_79) + xi_79;
          _data_pdfs_20_314_10[ctr_0] =
              forceTerm_14 +
              omega_bulk * (rho * (xi_61 + xi_62) + xi_50 - xi_64) + xi_64;
          _data_pdfs_20_315_10[ctr_0] =
              forceTerm_15 +
              omega_bulk * (rho * (xi_57 + xi_58) + xi_50 - xi_65) + xi_65;
          _data_pdfs_20_316_10[ctr_0] =
              forceTerm_16 +
              omega_bulk * (rho * (-xi_55 + xi_56) + xi_50 - xi_67) + xi_67;
          _data_pdfs_20_317_10[ctr_0] =
              forceTerm_17 +
              omega_bulk * (rho * (-xi_61 + xi_62) + xi_50 - xi_77) + xi_77;
          _data_pdfs_20_318_10[ctr_0] =
              forceTerm_18 +
              omega_bulk * (rho * (xi_59 + xi_60) + xi_50 - xi_75) + xi_75;
        }
      }
    }
  }
}
} // namespace internal_f11a519921c681cbc9d0b2f51454c920

void CollideSweepDoublePrecisionLeesEdwardsAVX::run(IBlock *block) {
  auto velocity = block->getData<field::GhostLayerField<double, 3>>(velocityID);
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  auto &omega_bulk = this->omega_bulk_;
  auto &points_up = this->points_up_;
  auto &points_down = this->points_down_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()));
  double *RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(force->xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(force->ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->zSize()) + 0));
  const int64_t _size_force_2 = int64_t(cell_idx_c(force->zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  internal_f11a519921c681cbc9d0b2f51454c920::
      collidesweepdoubleprecisionleesedwardsavx_collidesweepdoubleprecisionleesedwardsavx(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, omega_bulk, points_down,
          points_up);
}

void CollideSweepDoublePrecisionLeesEdwardsAVX::runOnCellInterval(
    const shared_ptr<StructuredBlockStorage> &blocks,
    const CellInterval &globalCellInterval, cell_idx_t ghostLayers,
    IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto velocity = block->getData<field::GhostLayerField<double, 3>>(velocityID);
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  auto &omega_bulk = this->omega_bulk_;
  auto &points_up = this->points_up_;
  auto &points_down = this->points_down_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force =
      force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs =
      pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()));
  double *RESTRICT const _data_velocity =
      velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.zSize()) + 0));
  const int64_t _size_force_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  internal_f11a519921c681cbc9d0b2f51454c920::
      collidesweepdoubleprecisionleesedwardsavx_collidesweepdoubleprecisionleesedwardsavx(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, omega_bulk, points_down,
          points_up);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif