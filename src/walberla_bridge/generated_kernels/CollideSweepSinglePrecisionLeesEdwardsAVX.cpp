// kernel generated with pystencils v0.4.3+12.g29e0e84, lbmpy v0.4.3+2.g0e17e61,
// lbmpy_walberla/pystencils_walberla from commit
// 55e6cf598e7e55f496ffaecd40bde632de76930e

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
//! \\file CollideSweepSinglePrecisionLeesEdwardsAVX.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepSinglePrecisionLeesEdwardsAVX.h"
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

namespace internal_9a18f2f4073cdcc5365cdfddb752069e {
static FUNC_PREFIX void
collidesweepsingleprecisionleesedwardsavx_collidesweepsingleprecisionleesedwardsavx(
    float *RESTRICT const _data_force, float *RESTRICT _data_pdfs,
    float *RESTRICT const _data_velocity, int64_t const _size_force_0,
    int64_t const _size_force_1, int64_t const _size_force_2,
    int64_t const _stride_force_1, int64_t const _stride_force_2,
    int64_t const _stride_force_3, int64_t const _stride_pdfs_1,
    int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3,
    int64_t const _stride_velocity_1, int64_t const _stride_velocity_2,
    int64_t const _stride_velocity_3, float omega_bulk, bool points_down,
    bool points_up) {
  const float xi_1 = omega_bulk * -0.5f + 1.0f;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_velocity_20_30 =
        _data_velocity + _stride_velocity_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    float *RESTRICT _data_velocity_20_32 =
        _data_velocity + _stride_velocity_2 * ctr_2 + 2 * _stride_velocity_3;
    float *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_velocity_20_31 =
        _data_velocity + _stride_velocity_2 * ctr_2 + _stride_velocity_3;
    float *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      float *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      float *RESTRICT _data_velocity_20_30_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_30;
      float *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      float *RESTRICT _data_velocity_20_32_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_32;
      float *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      float *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      float *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      float *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_velocity_20_31_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_31;
      float *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      {
        for (int64_t ctr_0 = 0; ctr_0 < (int64_t)((_size_force_0) / (8)) * (8);
             ctr_0 += 8) {
          const __m256 xi_63 = _mm256_load_ps(&_data_pdfs_20_313_10[ctr_0]);
          const __m256 xi_64 = _mm256_load_ps(&_data_pdfs_20_311_10[ctr_0]);
          const __m256 xi_65 = _mm256_load_ps(&_data_pdfs_20_310_10[ctr_0]);
          const __m256 xi_66 = _mm256_load_ps(&_data_pdfs_20_312_10[ctr_0]);
          const __m256 xi_67 = _mm256_load_ps(&_data_force_20_31_10[ctr_0]);
          const __m256 xi_68 = _mm256_load_ps(&_data_pdfs_20_32_10[ctr_0]);
          const __m256 xi_69 = _mm256_load_ps(&_data_pdfs_20_38_10[ctr_0]);
          const __m256 xi_70 = _mm256_load_ps(&_data_velocity_20_30_10[ctr_0]);
          const __m256 xi_71 = _mm256_load_ps(&_data_pdfs_20_316_10[ctr_0]);
          const __m256 xi_72 = _mm256_load_ps(&_data_pdfs_20_37_10[ctr_0]);
          const __m256 xi_73 = _mm256_load_ps(&_data_velocity_20_32_10[ctr_0]);
          const __m256 xi_74 = _mm256_load_ps(&_data_pdfs_20_33_10[ctr_0]);
          const __m256 xi_75 = _mm256_load_ps(&_data_pdfs_20_39_10[ctr_0]);
          const __m256 xi_76 = _mm256_load_ps(&_data_pdfs_20_314_10[ctr_0]);
          const __m256 xi_77 = _mm256_load_ps(&_data_pdfs_20_36_10[ctr_0]);
          const __m256 xi_78 = _mm256_load_ps(&_data_pdfs_20_34_10[ctr_0]);
          const __m256 xi_79 = _mm256_load_ps(&_data_pdfs_20_318_10[ctr_0]);
          const __m256 xi_80 = _mm256_load_ps(&_data_pdfs_20_35_10[ctr_0]);
          const __m256 xi_81 = _mm256_load_ps(&_data_pdfs_20_31_10[ctr_0]);
          const __m256 xi_82 = _mm256_load_ps(&_data_pdfs_20_317_10[ctr_0]);
          const __m256 xi_83 = _mm256_load_ps(&_data_velocity_20_31_10[ctr_0]);
          const __m256 xi_84 = _mm256_load_ps(&_data_force_20_32_10[ctr_0]);
          const __m256 xi_85 = _mm256_load_ps(&_data_pdfs_20_315_10[ctr_0]);
          const __m256 xi_86 = _mm256_load_ps(&_data_pdfs_20_30_10[ctr_0]);
          const __m256 xi_87 = _mm256_load_ps(&_data_force_20_30_10[ctr_0]);
          const __m256 xi_3 = _mm256_mul_ps(xi_67, xi_83);
          const __m256 xi_4 = _mm256_mul_ps(xi_73, xi_84);
          const __m256 xi_5 =
              _mm256_mul_ps(xi_83, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                 2.0f, 2.0f, 2.0f));
          const __m256 xi_6 =
              _mm256_add_ps(xi_5, _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                                                1.0f, 1.0f, 1.0f));
          const __m256 xi_7 = _mm256_mul_ps(
              xi_67, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f));
          const __m256 xi_9 = _mm256_mul_ps(
              xi_4, _mm256_set_ps(-0.166666666666667f, -0.166666666666667f,
                                  -0.166666666666667f, -0.166666666666667f,
                                  -0.166666666666667f, -0.166666666666667f,
                                  -0.166666666666667f, -0.166666666666667f));
          const __m256 xi_11 =
              _mm256_add_ps(xi_5, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f,
                                                -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_14 = _mm256_mul_ps(
              xi_87, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f));
          const __m256 xi_15 = _mm256_mul_ps(
              xi_3, _mm256_set_ps(-0.166666666666667f, -0.166666666666667f,
                                  -0.166666666666667f, -0.166666666666667f,
                                  -0.166666666666667f, -0.166666666666667f,
                                  -0.166666666666667f, -0.166666666666667f));
          const __m256 xi_16 = _mm256_add_ps(xi_15, xi_9);
          const __m256 xi_18 =
              _mm256_mul_ps(xi_73, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                 2.0f, 2.0f, 2.0f));
          const __m256 xi_19 =
              _mm256_add_ps(xi_18, _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                                                 1.0f, 1.0f, 1.0f));
          const __m256 xi_20 = _mm256_mul_ps(
              xi_84, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f));
          const __m256 xi_22 =
              _mm256_add_ps(xi_18, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f,
                                                 -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_23 = _mm256_mul_ps(
              xi_4, _mm256_set_ps(-0.0833333333333333f, -0.0833333333333333f,
                                  -0.0833333333333333f, -0.0833333333333333f,
                                  -0.0833333333333333f, -0.0833333333333333f,
                                  -0.0833333333333333f, -0.0833333333333333f));
          const __m256 xi_24 =
              _mm256_mul_ps(xi_83, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f,
                                                 3.0f, 3.0f, 3.0f));
          const __m256 xi_26 = _mm256_mul_ps(
              xi_87, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f));
          const __m256 xi_29 = _mm256_mul_ps(
              xi_67, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f));
          const __m256 xi_30 =
              _mm256_mul_ps(xi_24, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0));
          const __m256 xi_31 = _mm256_add_ps(
              _mm256_mul_ps(xi_5, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
              _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f));
          const __m256 xi_33 =
              _mm256_mul_ps(xi_73, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f,
                                                 3.0f, 3.0f, 3.0f));
          const __m256 xi_34 = _mm256_mul_ps(
              xi_84, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f));
          const __m256 xi_35 = _mm256_mul_ps(
              xi_3, _mm256_set_ps(-0.0833333333333333f, -0.0833333333333333f,
                                  -0.0833333333333333f, -0.0833333333333333f,
                                  -0.0833333333333333f, -0.0833333333333333f,
                                  -0.0833333333333333f, -0.0833333333333333f));
          const __m256 xi_36 =
              _mm256_mul_ps(xi_33, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0));
          const __m256 xi_37 = _mm256_add_ps(
              _mm256_mul_ps(xi_18, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0)),
              _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f));
          const __m256 xi_38 =
              _mm256_mul_ps(xi_73, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0));
          const __m256 xi_39 = (_mm256_mul_ps(xi_83, xi_83));
          const __m256 xi_41 = (_mm256_mul_ps(xi_73, xi_73));
          const __m256 xi_44 = _mm256_mul_ps(
              xi_83, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f));
          const __m256 xi_45 =
              _mm256_mul_ps(xi_39, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f,
                                                 0.25f, 0.25f, 0.25f, 0.25f));
          const __m256 xi_48 = _mm256_mul_ps(
              xi_73, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f,
                                   0.166666666666667f, 0.166666666666667f));
          const __m256 xi_49 =
              _mm256_mul_ps(xi_41, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f,
                                                 0.25f, 0.25f, 0.25f, 0.25f));
          const __m256 rho = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_add_ps(
                                                                  _mm256_add_ps(
                                                                      _mm256_add_ps(
                                                                          _mm256_add_ps(
                                                                              _mm256_add_ps(
                                                                                  xi_63,
                                                                                  xi_64),
                                                                              xi_65),
                                                                          xi_66),
                                                                      xi_68),
                                                                  xi_69),
                                                              xi_71),
                                                          xi_72),
                                                      xi_74),
                                                  xi_75),
                                              xi_76),
                                          xi_77),
                                      xi_78),
                                  xi_79),
                              xi_80),
                          xi_81),
                      xi_82),
                  xi_85),
              xi_86);
          const __m256 xi_40 =
              _mm256_mul_ps(rho, _mm256_set_ps(1.5f, 1.5f, 1.5f, 1.5f, 1.5f,
                                               1.5f, 1.5f, 1.5f));
          const __m256 u_0 = _mm256_add_ps(
              xi_70,
              _mm256_set_ps(
                  ((points_down && ctr_1 <= 0)
                       ? (1.0f)
                       : ((points_up && ctr_1 >= 63) ? (-1.0f) : (0.0f))) *
                      0.050000000000000003f,
                  ((points_down && ctr_1 <= 0)
                       ? (1.0f)
                       : ((points_up && ctr_1 >= 63) ? (-1.0f) : (0.0f))) *
                      0.050000000000000003f,
                  ((points_down && ctr_1 <= 0)
                       ? (1.0f)
                       : ((points_up && ctr_1 >= 63) ? (-1.0f) : (0.0f))) *
                      0.050000000000000003f,
                  ((points_down && ctr_1 <= 0)
                       ? (1.0f)
                       : ((points_up && ctr_1 >= 63) ? (-1.0f) : (0.0f))) *
                      0.050000000000000003f,
                  ((points_down && ctr_1 <= 0)
                       ? (1.0f)
                       : ((points_up && ctr_1 >= 63) ? (-1.0f) : (0.0f))) *
                      0.050000000000000003f,
                  ((points_down && ctr_1 <= 0)
                       ? (1.0f)
                       : ((points_up && ctr_1 >= 63) ? (-1.0f) : (0.0f))) *
                      0.050000000000000003f,
                  ((points_down && ctr_1 <= 0)
                       ? (1.0f)
                       : ((points_up && ctr_1 >= 63) ? (-1.0f) : (0.0f))) *
                      0.050000000000000003f,
                  ((points_down && ctr_1 <= 0)
                       ? (1.0f)
                       : ((points_up && ctr_1 >= 63) ? (-1.0f) : (0.0f))) *
                      0.050000000000000003f));
          const __m256 xi_2 = _mm256_mul_ps(u_0, xi_87);
          const __m256 xi_8 = _mm256_mul_ps(
              xi_2, _mm256_set_ps(-0.166666666666667f, -0.166666666666667f,
                                  -0.166666666666667f, -0.166666666666667f,
                                  -0.166666666666667f, -0.166666666666667f,
                                  -0.166666666666667f, -0.166666666666667f));
          const __m256 xi_10 = _mm256_add_ps(xi_8, xi_9);
          const __m256 xi_12 =
              _mm256_mul_ps(u_0, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                               2.0f, 2.0f, 2.0f));
          const __m256 xi_13 =
              _mm256_add_ps(xi_12, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f,
                                                 -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_17 =
              _mm256_add_ps(xi_12, _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                                                 1.0f, 1.0f, 1.0f));
          const __m256 xi_21 = _mm256_add_ps(xi_15, xi_8);
          const __m256 xi_25 = _mm256_add_ps(
              _mm256_mul_ps(xi_12, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0)),
              _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f));
          const __m256 xi_27 =
              _mm256_mul_ps(u_0, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f,
                                               3.0f, 3.0f, 3.0f));
          const __m256 xi_28 =
              _mm256_mul_ps(xi_27, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0));
          const __m256 xi_32 = _mm256_mul_ps(
              xi_2, _mm256_set_ps(-0.0833333333333333f, -0.0833333333333333f,
                                  -0.0833333333333333f, -0.0833333333333333f,
                                  -0.0833333333333333f, -0.0833333333333333f,
                                  -0.0833333333333333f, -0.0833333333333333f));
          const __m256 xi_42 = (_mm256_mul_ps(u_0, u_0));
          const __m256 xi_46 = _mm256_mul_ps(
              u_0, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                 0.166666666666667f, 0.166666666666667f,
                                 0.166666666666667f, 0.166666666666667f,
                                 0.166666666666667f, 0.166666666666667f));
          const __m256 xi_47 =
              _mm256_mul_ps(xi_42, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f,
                                                 0.25f, 0.25f, 0.25f, 0.25f));
          const __m256 forceTerm_0 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_2,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0)),
                      _mm256_mul_ps(xi_3,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0))),
                  _mm256_mul_ps(xi_4, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_1 = _mm256_mul_ps(
              _mm256_add_ps(_mm256_mul_ps(xi_6, xi_7), xi_10),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_2 = _mm256_mul_ps(
              _mm256_add_ps(_mm256_mul_ps(xi_11, xi_7), xi_10),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_3 = _mm256_mul_ps(
              _mm256_add_ps(_mm256_mul_ps(xi_13, xi_14), xi_16),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_4 = _mm256_mul_ps(
              _mm256_add_ps(_mm256_mul_ps(xi_14, xi_17), xi_16),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_5 = _mm256_mul_ps(
              _mm256_add_ps(_mm256_mul_ps(xi_19, xi_20), xi_21),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_6 = _mm256_mul_ps(
              _mm256_add_ps(_mm256_mul_ps(xi_20, xi_22), xi_21),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_7 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_29, _mm256_add_ps(xi_28, xi_6)),
                      _mm256_mul_ps(
                          _mm256_mul_ps(xi_26, _mm256_add_ps(xi_24, xi_25)),
                          _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                                        -1.0, -1.0))),
                  xi_23),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_8 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_26, _mm256_add_ps(xi_17, xi_24)),
                      _mm256_mul_ps(xi_29, _mm256_add_ps(xi_27, xi_6))),
                  xi_23),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_9 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_29, _mm256_add_ps(xi_11, xi_27)),
                      _mm256_mul_ps(xi_26, _mm256_add_ps(xi_13, xi_24))),
                  xi_23),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_10 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_26, _mm256_add_ps(xi_17, xi_30)),
                      _mm256_mul_ps(
                          _mm256_mul_ps(xi_29, _mm256_add_ps(xi_27, xi_31)),
                          _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                                        -1.0, -1.0))),
                  xi_23),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_11 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_34, _mm256_add_ps(xi_19, xi_24)),
                      _mm256_mul_ps(xi_29, _mm256_add_ps(xi_33, xi_6))),
                  xi_32),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_12 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_34, _mm256_add_ps(xi_19, xi_30)),
                      _mm256_mul_ps(
                          _mm256_mul_ps(xi_29, _mm256_add_ps(xi_31, xi_33)),
                          _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                                        -1.0, -1.0))),
                  xi_32),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_13 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_34, _mm256_add_ps(xi_19, xi_28)),
                      _mm256_mul_ps(
                          _mm256_mul_ps(xi_26, _mm256_add_ps(xi_25, xi_33)),
                          _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                                        -1.0, -1.0))),
                  xi_35),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_14 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_26, _mm256_add_ps(xi_17, xi_33)),
                      _mm256_mul_ps(xi_34, _mm256_add_ps(xi_19, xi_27))),
                  xi_35),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_15 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_29, _mm256_add_ps(xi_36, xi_6)),
                      _mm256_mul_ps(
                          _mm256_mul_ps(xi_34, _mm256_add_ps(xi_24, xi_37)),
                          _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                                        -1.0, -1.0))),
                  xi_32),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_16 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_29, _mm256_add_ps(xi_11, xi_33)),
                      _mm256_mul_ps(xi_34, _mm256_add_ps(xi_22, xi_24))),
                  xi_32),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_17 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_26, _mm256_add_ps(xi_13, xi_33)),
                      _mm256_mul_ps(xi_34, _mm256_add_ps(xi_22, xi_27))),
                  xi_35),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 forceTerm_18 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_26, _mm256_add_ps(xi_17, xi_36)),
                      _mm256_mul_ps(
                          _mm256_mul_ps(xi_34, _mm256_add_ps(xi_27, xi_37)),
                          _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                                        -1.0, -1.0))),
                  xi_35),
              _mm256_set_ps(xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1, xi_1));
          const __m256 u0Mu1 = _mm256_add_ps(
              _mm256_mul_ps(xi_83, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0)),
              u_0);
          const __m256 xi_51 = _mm256_mul_ps(
              u0Mu1, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f));
          const __m256 xi_52 =
              _mm256_mul_ps((_mm256_mul_ps(u0Mu1, u0Mu1)),
                            _mm256_set_ps(0.125f, 0.125f, 0.125f, 0.125f,
                                          0.125f, 0.125f, 0.125f, 0.125f));
          const __m256 u0Pu1 = _mm256_add_ps(u_0, xi_83);
          const __m256 xi_53 = _mm256_mul_ps(
              u0Pu1, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f));
          const __m256 xi_54 =
              _mm256_mul_ps((_mm256_mul_ps(u0Pu1, u0Pu1)),
                            _mm256_set_ps(0.125f, 0.125f, 0.125f, 0.125f,
                                          0.125f, 0.125f, 0.125f, 0.125f));
          const __m256 u1Pu2 = _mm256_add_ps(xi_73, xi_83);
          const __m256 xi_55 = _mm256_mul_ps(
              u1Pu2, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f));
          const __m256 xi_56 =
              _mm256_mul_ps((_mm256_mul_ps(u1Pu2, u1Pu2)),
                            _mm256_set_ps(0.125f, 0.125f, 0.125f, 0.125f,
                                          0.125f, 0.125f, 0.125f, 0.125f));
          const __m256 u1Mu2 = _mm256_add_ps(xi_38, xi_83);
          const __m256 xi_57 = _mm256_mul_ps(
              u1Mu2, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f));
          const __m256 xi_58 =
              _mm256_mul_ps((_mm256_mul_ps(u1Mu2, u1Mu2)),
                            _mm256_set_ps(0.125f, 0.125f, 0.125f, 0.125f,
                                          0.125f, 0.125f, 0.125f, 0.125f));
          const __m256 u0Mu2 = _mm256_add_ps(u_0, xi_38);
          const __m256 xi_59 = _mm256_mul_ps(
              u0Mu2, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f));
          const __m256 xi_60 =
              _mm256_mul_ps((_mm256_mul_ps(u0Mu2, u0Mu2)),
                            _mm256_set_ps(0.125f, 0.125f, 0.125f, 0.125f,
                                          0.125f, 0.125f, 0.125f, 0.125f));
          const __m256 u0Pu2 = _mm256_add_ps(u_0, xi_73);
          const __m256 xi_61 = _mm256_mul_ps(
              u0Pu2, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f,
                                   0.0833333333333333f, 0.0833333333333333f));
          const __m256 xi_62 =
              _mm256_mul_ps((_mm256_mul_ps(u0Pu2, u0Pu2)),
                            _mm256_set_ps(0.125f, 0.125f, 0.125f, 0.125f,
                                          0.125f, 0.125f, 0.125f, 0.125f));
          const __m256 f_eq_common = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(_mm256_mul_ps(xi_39, xi_40),
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0)),
                      _mm256_mul_ps(_mm256_mul_ps(xi_40, xi_41),
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0))),
                  _mm256_mul_ps(_mm256_mul_ps(xi_40, xi_42),
                                _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0))),
              rho);
          const __m256 xi_43 = _mm256_mul_ps(
              f_eq_common,
              _mm256_set_ps(0.0555555555555556f, 0.0555555555555556f,
                            0.0555555555555556f, 0.0555555555555556f,
                            0.0555555555555556f, 0.0555555555555556f,
                            0.0555555555555556f, 0.0555555555555556f));
          const __m256 xi_50 = _mm256_mul_ps(
              f_eq_common,
              _mm256_set_ps(0.0277777777777778f, 0.0277777777777778f,
                            0.0277777777777778f, 0.0277777777777778f,
                            0.0277777777777778f, 0.0277777777777778f,
                            0.0277777777777778f, 0.0277777777777778f));
          _mm256_store_ps(
              &_data_pdfs_20_30_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_mul_ps(
                                  f_eq_common,
                                  _mm256_set_ps(
                                      0.333333333333333f, 0.333333333333333f,
                                      0.333333333333333f, 0.333333333333333f,
                                      0.333333333333333f, 0.333333333333333f,
                                      0.333333333333333f, 0.333333333333333f)),
                              _mm256_mul_ps(xi_86,
                                            _mm256_set_ps(-1.0, -1.0, -1.0,
                                                          -1.0, -1.0, -1.0,
                                                          -1.0, -1.0))),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_0),
                  xi_86));
          _mm256_store_ps(
              &_data_pdfs_20_31_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(rho,
                                                _mm256_add_ps(xi_44, xi_45)),
                                  _mm256_mul_ps(xi_81,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_43),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_1),
                  xi_81));
          _mm256_store_ps(
              &_data_pdfs_20_32_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   xi_44,
                                                   _mm256_set_ps(
                                                       -1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                                               xi_45)),
                                  _mm256_mul_ps(xi_68,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_43),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_2),
                  xi_68));
          _mm256_store_ps(
              &_data_pdfs_20_33_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   xi_46,
                                                   _mm256_set_ps(
                                                       -1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                                               xi_47)),
                                  _mm256_mul_ps(xi_74,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_43),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_3),
                  xi_74));
          _mm256_store_ps(
              &_data_pdfs_20_34_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(rho,
                                                _mm256_add_ps(xi_46, xi_47)),
                                  _mm256_mul_ps(xi_78,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_43),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_4),
                  xi_78));
          _mm256_store_ps(
              &_data_pdfs_20_35_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(rho,
                                                _mm256_add_ps(xi_48, xi_49)),
                                  _mm256_mul_ps(xi_80,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_43),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_5),
                  xi_80));
          _mm256_store_ps(
              &_data_pdfs_20_36_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   xi_48,
                                                   _mm256_set_ps(
                                                       -1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                                               xi_49)),
                                  _mm256_mul_ps(xi_77,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_43),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_6),
                  xi_77));
          _mm256_store_ps(
              &_data_pdfs_20_37_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   xi_51,
                                                   _mm256_set_ps(
                                                       -1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                                               xi_52)),
                                  _mm256_mul_ps(xi_72,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_50),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_7),
                  xi_72));
          _mm256_store_ps(
              &_data_pdfs_20_38_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(rho,
                                                _mm256_add_ps(xi_53, xi_54)),
                                  _mm256_mul_ps(xi_69,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_50),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_8),
                  xi_69));
          _mm256_store_ps(
              &_data_pdfs_20_39_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   xi_53,
                                                   _mm256_set_ps(
                                                       -1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                                               xi_54)),
                                  _mm256_mul_ps(xi_75,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_50),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_9),
                  xi_75));
          _mm256_store_ps(
              &_data_pdfs_20_310_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(rho,
                                                _mm256_add_ps(xi_51, xi_52)),
                                  _mm256_mul_ps(xi_65,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_50),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_10),
                  xi_65));
          _mm256_store_ps(
              &_data_pdfs_20_311_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(rho,
                                                _mm256_add_ps(xi_55, xi_56)),
                                  _mm256_mul_ps(xi_64,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_50),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_11),
                  xi_64));
          _mm256_store_ps(
              &_data_pdfs_20_312_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   xi_57,
                                                   _mm256_set_ps(
                                                       -1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                                               xi_58)),
                                  _mm256_mul_ps(xi_66,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_50),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_12),
                  xi_66));
          _mm256_store_ps(
              &_data_pdfs_20_313_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   xi_59,
                                                   _mm256_set_ps(
                                                       -1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                                               xi_60)),
                                  _mm256_mul_ps(xi_63,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_50),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_13),
                  xi_63));
          _mm256_store_ps(
              &_data_pdfs_20_314_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(rho,
                                                _mm256_add_ps(xi_61, xi_62)),
                                  _mm256_mul_ps(xi_76,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_50),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_14),
                  xi_76));
          _mm256_store_ps(
              &_data_pdfs_20_315_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(rho,
                                                _mm256_add_ps(xi_57, xi_58)),
                                  _mm256_mul_ps(xi_85,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_50),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_15),
                  xi_85));
          _mm256_store_ps(
              &_data_pdfs_20_316_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   xi_55,
                                                   _mm256_set_ps(
                                                       -1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                                               xi_56)),
                                  _mm256_mul_ps(xi_71,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_50),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_16),
                  xi_71));
          _mm256_store_ps(
              &_data_pdfs_20_317_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   xi_61,
                                                   _mm256_set_ps(
                                                       -1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                                               xi_62)),
                                  _mm256_mul_ps(xi_82,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_50),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_17),
                  xi_82));
          _mm256_store_ps(
              &_data_pdfs_20_318_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(rho,
                                                _mm256_add_ps(xi_59, xi_60)),
                                  _mm256_mul_ps(xi_79,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_50),
                          _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk, omega_bulk,
                                        omega_bulk, omega_bulk)),
                      forceTerm_18),
                  xi_79));
        }
        for (int64_t ctr_0 = (int64_t)((_size_force_0) / (8)) * (8);
             ctr_0 < _size_force_0; ctr_0 += 1) {
          const float xi_63 = _data_pdfs_20_313_10[ctr_0];
          const float xi_64 = _data_pdfs_20_311_10[ctr_0];
          const float xi_65 = _data_pdfs_20_310_10[ctr_0];
          const float xi_66 = _data_pdfs_20_312_10[ctr_0];
          const float xi_67 = _data_force_20_31_10[ctr_0];
          const float xi_68 = _data_pdfs_20_32_10[ctr_0];
          const float xi_69 = _data_pdfs_20_38_10[ctr_0];
          const float xi_70 = _data_velocity_20_30_10[ctr_0];
          const float xi_71 = _data_pdfs_20_316_10[ctr_0];
          const float xi_72 = _data_pdfs_20_37_10[ctr_0];
          const float xi_73 = _data_velocity_20_32_10[ctr_0];
          const float xi_74 = _data_pdfs_20_33_10[ctr_0];
          const float xi_75 = _data_pdfs_20_39_10[ctr_0];
          const float xi_76 = _data_pdfs_20_314_10[ctr_0];
          const float xi_77 = _data_pdfs_20_36_10[ctr_0];
          const float xi_78 = _data_pdfs_20_34_10[ctr_0];
          const float xi_79 = _data_pdfs_20_318_10[ctr_0];
          const float xi_80 = _data_pdfs_20_35_10[ctr_0];
          const float xi_81 = _data_pdfs_20_31_10[ctr_0];
          const float xi_82 = _data_pdfs_20_317_10[ctr_0];
          const float xi_83 = _data_velocity_20_31_10[ctr_0];
          const float xi_84 = _data_force_20_32_10[ctr_0];
          const float xi_85 = _data_pdfs_20_315_10[ctr_0];
          const float xi_86 = _data_pdfs_20_30_10[ctr_0];
          const float xi_87 = _data_force_20_30_10[ctr_0];
          const float xi_3 = xi_67 * xi_83;
          const float xi_4 = xi_73 * xi_84;
          const float xi_5 = xi_83 * 2.0f;
          const float xi_6 = xi_5 + 1.0f;
          const float xi_7 = xi_67 * 0.166666666666667f;
          const float xi_9 = xi_4 * -0.166666666666667f;
          const float xi_11 = xi_5 - 1.0f;
          const float xi_14 = xi_87 * 0.166666666666667f;
          const float xi_15 = xi_3 * -0.166666666666667f;
          const float xi_16 = xi_15 + xi_9;
          const float xi_18 = xi_73 * 2.0f;
          const float xi_19 = xi_18 + 1.0f;
          const float xi_20 = xi_84 * 0.166666666666667f;
          const float xi_22 = xi_18 - 1.0f;
          const float xi_23 = xi_4 * -0.0833333333333333f;
          const float xi_24 = xi_83 * 3.0f;
          const float xi_26 = xi_87 * 0.0833333333333333f;
          const float xi_29 = xi_67 * 0.0833333333333333f;
          const float xi_30 = -xi_24;
          const float xi_31 = -xi_5 + 1.0f;
          const float xi_33 = xi_73 * 3.0f;
          const float xi_34 = xi_84 * 0.0833333333333333f;
          const float xi_35 = xi_3 * -0.0833333333333333f;
          const float xi_36 = -xi_33;
          const float xi_37 = -xi_18 + 1.0f;
          const float xi_38 = -xi_73;
          const float xi_39 = (xi_83 * xi_83);
          const float xi_41 = (xi_73 * xi_73);
          const float xi_44 = xi_83 * 0.166666666666667f;
          const float xi_45 = xi_39 * 0.25f;
          const float xi_48 = xi_73 * 0.166666666666667f;
          const float xi_49 = xi_41 * 0.25f;
          const float rho = xi_63 + xi_64 + xi_65 + xi_66 + xi_68 + xi_69 +
                            xi_71 + xi_72 + xi_74 + xi_75 + xi_76 + xi_77 +
                            xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_85 +
                            xi_86;
          const float xi_40 = rho * 1.5f;
          const float u_0 =
              xi_70 + ((points_down && ctr_1 <= 0)
                           ? (1.0f)
                           : ((points_up && ctr_1 >= 63) ? (-1.0f) : (0.0f))) *
                          0.050000000000000003f;
          const float xi_2 = u_0 * xi_87;
          const float xi_8 = xi_2 * -0.166666666666667f;
          const float xi_10 = xi_8 + xi_9;
          const float xi_12 = u_0 * 2.0f;
          const float xi_13 = xi_12 - 1.0f;
          const float xi_17 = xi_12 + 1.0f;
          const float xi_21 = xi_15 + xi_8;
          const float xi_25 = -xi_12 + 1.0f;
          const float xi_27 = u_0 * 3.0f;
          const float xi_28 = -xi_27;
          const float xi_32 = xi_2 * -0.0833333333333333f;
          const float xi_42 = (u_0 * u_0);
          const float xi_46 = u_0 * 0.166666666666667f;
          const float xi_47 = xi_42 * 0.25f;
          const float forceTerm_0 = xi_1 * (-xi_2 - xi_3 - xi_4);
          const float forceTerm_1 = xi_1 * (xi_10 + xi_6 * xi_7);
          const float forceTerm_2 = xi_1 * (xi_10 + xi_11 * xi_7);
          const float forceTerm_3 = xi_1 * (xi_13 * xi_14 + xi_16);
          const float forceTerm_4 = xi_1 * (xi_14 * xi_17 + xi_16);
          const float forceTerm_5 = xi_1 * (xi_19 * xi_20 + xi_21);
          const float forceTerm_6 = xi_1 * (xi_20 * xi_22 + xi_21);
          const float forceTerm_7 =
              xi_1 * (xi_23 - xi_26 * (xi_24 + xi_25) + xi_29 * (xi_28 + xi_6));
          const float forceTerm_8 =
              xi_1 * (xi_23 + xi_26 * (xi_17 + xi_24) + xi_29 * (xi_27 + xi_6));
          const float forceTerm_9 = xi_1 * (xi_23 + xi_26 * (xi_13 + xi_24) +
                                            xi_29 * (xi_11 + xi_27));
          const float forceTerm_10 = xi_1 * (xi_23 + xi_26 * (xi_17 + xi_30) -
                                             xi_29 * (xi_27 + xi_31));
          const float forceTerm_11 =
              xi_1 * (xi_29 * (xi_33 + xi_6) + xi_32 + xi_34 * (xi_19 + xi_24));
          const float forceTerm_12 = xi_1 * (-xi_29 * (xi_31 + xi_33) + xi_32 +
                                             xi_34 * (xi_19 + xi_30));
          const float forceTerm_13 = xi_1 * (-xi_26 * (xi_25 + xi_33) +
                                             xi_34 * (xi_19 + xi_28) + xi_35);
          const float forceTerm_14 = xi_1 * (xi_26 * (xi_17 + xi_33) +
                                             xi_34 * (xi_19 + xi_27) + xi_35);
          const float forceTerm_15 =
              xi_1 * (xi_29 * (xi_36 + xi_6) + xi_32 - xi_34 * (xi_24 + xi_37));
          const float forceTerm_16 = xi_1 * (xi_29 * (xi_11 + xi_33) + xi_32 +
                                             xi_34 * (xi_22 + xi_24));
          const float forceTerm_17 = xi_1 * (xi_26 * (xi_13 + xi_33) +
                                             xi_34 * (xi_22 + xi_27) + xi_35);
          const float forceTerm_18 = xi_1 * (xi_26 * (xi_17 + xi_36) -
                                             xi_34 * (xi_27 + xi_37) + xi_35);
          const float u0Mu1 = u_0 - xi_83;
          const float xi_51 = u0Mu1 * 0.0833333333333333f;
          const float xi_52 = (u0Mu1 * u0Mu1) * 0.125f;
          const float u0Pu1 = u_0 + xi_83;
          const float xi_53 = u0Pu1 * 0.0833333333333333f;
          const float xi_54 = (u0Pu1 * u0Pu1) * 0.125f;
          const float u1Pu2 = xi_73 + xi_83;
          const float xi_55 = u1Pu2 * 0.0833333333333333f;
          const float xi_56 = (u1Pu2 * u1Pu2) * 0.125f;
          const float u1Mu2 = xi_38 + xi_83;
          const float xi_57 = u1Mu2 * 0.0833333333333333f;
          const float xi_58 = (u1Mu2 * u1Mu2) * 0.125f;
          const float u0Mu2 = u_0 + xi_38;
          const float xi_59 = u0Mu2 * 0.0833333333333333f;
          const float xi_60 = (u0Mu2 * u0Mu2) * 0.125f;
          const float u0Pu2 = u_0 + xi_73;
          const float xi_61 = u0Pu2 * 0.0833333333333333f;
          const float xi_62 = (u0Pu2 * u0Pu2) * 0.125f;
          const float f_eq_common =
              rho - xi_39 * xi_40 - xi_40 * xi_41 - xi_40 * xi_42;
          const float xi_43 = f_eq_common * 0.0555555555555556f;
          const float xi_50 = f_eq_common * 0.0277777777777778f;
          _data_pdfs_20_30_10[ctr_0] =
              forceTerm_0 +
              omega_bulk * (f_eq_common * 0.333333333333333f - xi_86) + xi_86;
          _data_pdfs_20_31_10[ctr_0] =
              forceTerm_1 +
              omega_bulk * (rho * (xi_44 + xi_45) + xi_43 - xi_81) + xi_81;
          _data_pdfs_20_32_10[ctr_0] =
              forceTerm_2 +
              omega_bulk * (rho * (-xi_44 + xi_45) + xi_43 - xi_68) + xi_68;
          _data_pdfs_20_33_10[ctr_0] =
              forceTerm_3 +
              omega_bulk * (rho * (-xi_46 + xi_47) + xi_43 - xi_74) + xi_74;
          _data_pdfs_20_34_10[ctr_0] =
              forceTerm_4 +
              omega_bulk * (rho * (xi_46 + xi_47) + xi_43 - xi_78) + xi_78;
          _data_pdfs_20_35_10[ctr_0] =
              forceTerm_5 +
              omega_bulk * (rho * (xi_48 + xi_49) + xi_43 - xi_80) + xi_80;
          _data_pdfs_20_36_10[ctr_0] =
              forceTerm_6 +
              omega_bulk * (rho * (-xi_48 + xi_49) + xi_43 - xi_77) + xi_77;
          _data_pdfs_20_37_10[ctr_0] =
              forceTerm_7 +
              omega_bulk * (rho * (-xi_51 + xi_52) + xi_50 - xi_72) + xi_72;
          _data_pdfs_20_38_10[ctr_0] =
              forceTerm_8 +
              omega_bulk * (rho * (xi_53 + xi_54) + xi_50 - xi_69) + xi_69;
          _data_pdfs_20_39_10[ctr_0] =
              forceTerm_9 +
              omega_bulk * (rho * (-xi_53 + xi_54) + xi_50 - xi_75) + xi_75;
          _data_pdfs_20_310_10[ctr_0] =
              forceTerm_10 +
              omega_bulk * (rho * (xi_51 + xi_52) + xi_50 - xi_65) + xi_65;
          _data_pdfs_20_311_10[ctr_0] =
              forceTerm_11 +
              omega_bulk * (rho * (xi_55 + xi_56) + xi_50 - xi_64) + xi_64;
          _data_pdfs_20_312_10[ctr_0] =
              forceTerm_12 +
              omega_bulk * (rho * (-xi_57 + xi_58) + xi_50 - xi_66) + xi_66;
          _data_pdfs_20_313_10[ctr_0] =
              forceTerm_13 +
              omega_bulk * (rho * (-xi_59 + xi_60) + xi_50 - xi_63) + xi_63;
          _data_pdfs_20_314_10[ctr_0] =
              forceTerm_14 +
              omega_bulk * (rho * (xi_61 + xi_62) + xi_50 - xi_76) + xi_76;
          _data_pdfs_20_315_10[ctr_0] =
              forceTerm_15 +
              omega_bulk * (rho * (xi_57 + xi_58) + xi_50 - xi_85) + xi_85;
          _data_pdfs_20_316_10[ctr_0] =
              forceTerm_16 +
              omega_bulk * (rho * (-xi_55 + xi_56) + xi_50 - xi_71) + xi_71;
          _data_pdfs_20_317_10[ctr_0] =
              forceTerm_17 +
              omega_bulk * (rho * (-xi_61 + xi_62) + xi_50 - xi_82) + xi_82;
          _data_pdfs_20_318_10[ctr_0] =
              forceTerm_18 +
              omega_bulk * (rho * (xi_59 + xi_60) + xi_50 - xi_79) + xi_79;
        }
      }
    }
  }
}
} // namespace internal_9a18f2f4073cdcc5365cdfddb752069e

void CollideSweepSinglePrecisionLeesEdwardsAVX::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  auto velocity = block->getData<field::GhostLayerField<float, 3>>(velocityID);

  auto &omega_bulk = this->omega_bulk_;
  auto &points_up = this->points_up_;
  auto &points_down = this->points_down_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()));
  float *RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
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
  internal_9a18f2f4073cdcc5365cdfddb752069e::
      collidesweepsingleprecisionleesedwardsavx_collidesweepsingleprecisionleesedwardsavx(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, omega_bulk, points_down,
          points_up);
}

void CollideSweepSinglePrecisionLeesEdwardsAVX::runOnCellInterval(
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

  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  auto velocity = block->getData<field::GhostLayerField<float, 3>>(velocityID);

  auto &omega_bulk = this->omega_bulk_;
  auto &points_up = this->points_up_;
  auto &points_down = this->points_down_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force =
      force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()));
  float *RESTRICT const _data_velocity =
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
  internal_9a18f2f4073cdcc5365cdfddb752069e::
      collidesweepsingleprecisionleesedwardsavx_collidesweepsingleprecisionleesedwardsavx(
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