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
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.3, lbmpy v1.3.3,
// lbmpy_walberla/pystencils_walberla from waLBerla commit
// 04f4adbdfc0af983e2d9b72e244d775f37d77034

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
    int64_t const _size_force_0, int64_t const _size_force_1,
    int64_t const _size_force_2, int64_t const _stride_force_1,
    int64_t const _stride_force_2, int64_t const _stride_force_3,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, float grid_size, float omega_shear,
    float v_s) {
  const float xi_0 = ((1.0f) / (omega_shear * -0.25f + 2.0f));
  const float rr_0 = xi_0 * (omega_shear * -2.0f + 4.0f);
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      {
        for (int64_t ctr_0 = 0; ctr_0 < (int64_t)((_size_force_0) / (8)) * (8);
             ctr_0 += 8) {
          const __m256 xi_25 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          14 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_26 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          13 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_27 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          11 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_28 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          5 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_29 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          10 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_30 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          15 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_31 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          18 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_32 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          6 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_33 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          _stride_pdfs_3 + ctr_0]);
          const __m256 xi_34 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          17 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_35 = _mm256_load_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          16 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_36 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          9 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_37 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          2 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_38 = _mm256_load_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          8 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_39 =
              _mm256_load_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 +
                                         _stride_pdfs_2 * ctr_2 + ctr_0]);
          const __m256 xi_40 =
              _mm256_load_ps(&_data_force[_stride_force_1 * ctr_1 +
                                          _stride_force_2 * ctr_2 + ctr_0]);
          const __m256 xi_41 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          7 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_42 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          3 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_43 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          4 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_44 = _mm256_loadu_ps(
              &_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 +
                           2 * _stride_force_3 + ctr_0]);
          const __m256 xi_45 = _mm256_loadu_ps(
              &_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 +
                           _stride_force_3 + ctr_0]);
          const __m256 xi_46 = _mm256_loadu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          12 * _stride_pdfs_3 + ctr_0]);
          const __m256 xi_3 = xi_33;
          const __m256 xi_4 = xi_38;
          const __m256 xi_5 = xi_42;
          const __m256 xi_6 = xi_29;
          const __m256 xi_7 = xi_45;
          const __m256 xi_8 = xi_31;
          const __m256 xi_9 = xi_25;
          const __m256 xi_10 = xi_46;
          const __m256 xi_11 = xi_36;
          const __m256 xi_12 = xi_27;
          const __m256 xi_13 = xi_39;
          const __m256 xi_14 = xi_30;
          const __m256 xi_15 = xi_34;
          const __m256 xi_16 = xi_32;
          const __m256 xi_17 = xi_35;
          const __m256 xi_18 = xi_37;
          const __m256 xi_19 = xi_26;
          const __m256 xi_20 = xi_44;
          const __m256 xi_21 = xi_41;
          const __m256 xi_22 = xi_28;
          const __m256 xi_23 = xi_40;
          const __m256 xi_24 = xi_43;
          const __m256 vel0Term = _mm256_add_ps(
              _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_24, xi_4), xi_6),
                            xi_8),
              xi_9);
          const __m256 vel1Term = _mm256_add_ps(
              _mm256_add_ps(_mm256_add_ps(xi_12, xi_14), xi_21), xi_3);
          const __m256 vel2Term =
              _mm256_add_ps(_mm256_add_ps(xi_10, xi_19), xi_22);
          const __m256 rho = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(vel0Term, vel1Term),
                                          vel2Term),
                                      xi_11),
                                  xi_13),
                              xi_15),
                          xi_16),
                      xi_17),
                  xi_18),
              xi_5);
          const __m256 xi_1 = _mm256_div_ps(
              _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f),
              rho);
          const __m256 u_0 = _mm256_add_ps(
              _mm256_mul_ps(
                  xi_1,
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_11,
                                          _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f)),
                                      _mm256_mul_ps(
                                          xi_15,
                                          _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f))),
                                  _mm256_mul_ps(
                                      xi_19, _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                           -1.0f, -1.0f, -1.0f,
                                                           -1.0f, -1.0f))),
                              _mm256_mul_ps(xi_21,
                                            _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                          -1.0f, -1.0f, -1.0f,
                                                          -1.0f, -1.0f))),
                          _mm256_mul_ps(xi_5, _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f))),
                      vel0Term)),
              _mm256_mul_ps(_mm256_mul_ps(xi_1, xi_23),
                            _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                                          0.5f, 0.5f)));
          const __m256 u_1 = _mm256_add_ps(
              _mm256_mul_ps(
                  xi_1,
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              xi_10,
                                              _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f)),
                                          _mm256_mul_ps(
                                              xi_11,
                                              _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f))),
                                      _mm256_mul_ps(
                                          xi_17,
                                          _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f))),
                                  _mm256_mul_ps(
                                      xi_18, _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                           -1.0f, -1.0f, -1.0f,
                                                           -1.0f, -1.0f))),
                              _mm256_mul_ps(xi_6,
                                            _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                          -1.0f, -1.0f, -1.0f,
                                                          -1.0f, -1.0f))),
                          vel1Term),
                      xi_4)),
              _mm256_mul_ps(_mm256_mul_ps(xi_1, xi_7),
                            _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                                          0.5f, 0.5f)));
          const __m256 u_2 = _mm256_add_ps(
              _mm256_mul_ps(
                  xi_1,
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  xi_14,
                                                  _mm256_set_ps(-1.0f, -1.0f,
                                                                -1.0f, -1.0f,
                                                                -1.0f, -1.0f,
                                                                -1.0f, -1.0f)),
                                              _mm256_mul_ps(
                                                  xi_15,
                                                  _mm256_set_ps(-1.0f, -1.0f,
                                                                -1.0f, -1.0f,
                                                                -1.0f, -1.0f,
                                                                -1.0f, -1.0f))),
                                          _mm256_mul_ps(
                                              xi_16,
                                              _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f))),
                                      _mm256_mul_ps(
                                          xi_17,
                                          _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f))),
                                  _mm256_mul_ps(
                                      xi_8, _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                          -1.0f, -1.0f, -1.0f,
                                                          -1.0f, -1.0f))),
                              vel2Term),
                          xi_12),
                      xi_9)),
              _mm256_mul_ps(_mm256_mul_ps(xi_1, xi_20),
                            _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                                          0.5f, 0.5f)));
          const __m256 forceTerm_0 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_0, xi_23),
                                            _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                          -1.0f, -1.0f, -1.0f,
                                                          -1.0f, -1.0f)),
                              _mm256_mul_ps(_mm256_mul_ps(u_1, xi_7),
                                            _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                          -1.0f, -1.0f, -1.0f,
                                                          -1.0f, -1.0f))),
                          _mm256_mul_ps(_mm256_mul_ps(u_2, xi_20),
                                        _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                      -1.0f, -1.0f, -1.0f,
                                                      -1.0f, -1.0f))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(_mm256_mul_ps(u_0, xi_23),
                                        _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f,
                                                      0.5f, 0.5f, 0.5f, 0.5f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(_mm256_mul_ps(u_1, xi_7),
                                    _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                                                  0.5f, 0.5f, 0.5f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(_mm256_mul_ps(u_2, xi_20),
                                _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                                              0.5f, 0.5f, 0.5f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_1 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_7,
                                          _mm256_set_ps(0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f)),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(u_1, xi_7),
                                          _mm256_set_ps(0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_7, _mm256_set_ps(
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f)),
                                      _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                    rr_0, rr_0, rr_0, rr_0))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(u_0, xi_23),
                                  _mm256_set_ps(-0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_2, xi_20),
                              _mm256_set_ps(
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f,
                                  -0.16666666666666666f))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_0, xi_23),
                              _mm256_set_ps(
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f,
                                  0.083333333333333329f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_2, xi_20),
                          _mm256_set_ps(
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_1, xi_7),
                      _mm256_set_ps(
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_2 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_7,
                                          _mm256_set_ps(-0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f)),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              xi_7, _mm256_set_ps(
                                                        0.083333333333333329f,
                                                        0.083333333333333329f,
                                                        0.083333333333333329f,
                                                        0.083333333333333329f,
                                                        0.083333333333333329f,
                                                        0.083333333333333329f,
                                                        0.083333333333333329f,
                                                        0.083333333333333329f)),
                                          _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                        rr_0, rr_0, rr_0,
                                                        rr_0))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_1, xi_7),
                                      _mm256_set_ps(0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(u_0, xi_23),
                                  _mm256_set_ps(-0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_2, xi_20),
                              _mm256_set_ps(
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f,
                                  -0.16666666666666666f))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_0, xi_23),
                              _mm256_set_ps(
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f,
                                  0.083333333333333329f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_2, xi_20),
                          _mm256_set_ps(
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_1, xi_7),
                      _mm256_set_ps(
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_3 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_23,
                                          _mm256_set_ps(-0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f)),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              xi_23,
                                              _mm256_set_ps(
                                                  0.083333333333333329f,
                                                  0.083333333333333329f,
                                                  0.083333333333333329f,
                                                  0.083333333333333329f,
                                                  0.083333333333333329f,
                                                  0.083333333333333329f,
                                                  0.083333333333333329f,
                                                  0.083333333333333329f)),
                                          _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                        rr_0, rr_0, rr_0,
                                                        rr_0))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_0, xi_23),
                                      _mm256_set_ps(0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(u_1, xi_7),
                                  _mm256_set_ps(-0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_2, xi_20),
                              _mm256_set_ps(
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f,
                                  -0.16666666666666666f))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_1, xi_7),
                              _mm256_set_ps(
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f,
                                  0.083333333333333329f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_2, xi_20),
                          _mm256_set_ps(
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_0, xi_23),
                      _mm256_set_ps(
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_4 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_23,
                                          _mm256_set_ps(0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f)),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(u_0, xi_23),
                                          _mm256_set_ps(0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_23, _mm256_set_ps(
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f)),
                                      _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                    rr_0, rr_0, rr_0, rr_0))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(u_1, xi_7),
                                  _mm256_set_ps(-0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_2, xi_20),
                              _mm256_set_ps(
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f,
                                  -0.16666666666666666f))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_1, xi_7),
                              _mm256_set_ps(
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f,
                                  0.083333333333333329f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_2, xi_20),
                          _mm256_set_ps(
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_0, xi_23),
                      _mm256_set_ps(
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_5 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_20,
                                          _mm256_set_ps(0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f)),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(u_2, xi_20),
                                          _mm256_set_ps(0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f,
                                                        0.33333333333333331f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_20, _mm256_set_ps(
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f,
                                                     -0.083333333333333329f)),
                                      _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                    rr_0, rr_0, rr_0, rr_0))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(u_0, xi_23),
                                  _mm256_set_ps(-0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_1, xi_7),
                              _mm256_set_ps(
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f,
                                  -0.16666666666666666f))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_0, xi_23),
                              _mm256_set_ps(
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f,
                                  0.083333333333333329f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_1, xi_7),
                          _mm256_set_ps(
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_2, xi_20),
                      _mm256_set_ps(
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_6 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_20,
                                          _mm256_set_ps(-0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f,
                                                        -0.16666666666666666f)),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              xi_20,
                                              _mm256_set_ps(
                                                  0.083333333333333329f,
                                                  0.083333333333333329f,
                                                  0.083333333333333329f,
                                                  0.083333333333333329f,
                                                  0.083333333333333329f,
                                                  0.083333333333333329f,
                                                  0.083333333333333329f,
                                                  0.083333333333333329f)),
                                          _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                        rr_0, rr_0, rr_0,
                                                        rr_0))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_2, xi_20),
                                      _mm256_set_ps(0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f,
                                                    0.33333333333333331f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(u_0, xi_23),
                                  _mm256_set_ps(-0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f,
                                                -0.16666666666666666f))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_1, xi_7),
                              _mm256_set_ps(
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f, -0.16666666666666666f,
                                  -0.16666666666666666f,
                                  -0.16666666666666666f))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_0, xi_23),
                              _mm256_set_ps(
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f, 0.083333333333333329f,
                                  0.083333333333333329f,
                                  0.083333333333333329f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_1, xi_7),
                          _mm256_set_ps(
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f,
                              0.083333333333333329f, 0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_2, xi_20),
                      _mm256_set_ps(
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f,
                          -0.16666666666666666f, -0.16666666666666666f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_7 = _mm256_add_ps(
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
                                                              _mm256_mul_ps(
                                                                  xi_7,
                                                                  _mm256_set_ps(
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f)),
                                                              _mm256_mul_ps(
                                                                  xi_23,
                                                                  _mm256_set_ps(
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f))),
                                                          _mm256_mul_ps(
                                                              _mm256_mul_ps(
                                                                  xi_23,
                                                                  _mm256_set_ps(
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f)),
                                                              _mm256_set_ps(
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0))),
                                                      _mm256_mul_ps(
                                                          _mm256_mul_ps(u_0,
                                                                        xi_23),
                                                          _mm256_set_ps(
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f))),
                                                  _mm256_mul_ps(
                                                      _mm256_mul_ps(u_1, xi_7),
                                                      _mm256_set_ps(
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f))),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(u_0, xi_7),
                                                  _mm256_set_ps(
                                                      -0.25f, -0.25f, -0.25f,
                                                      -0.25f, -0.25f, -0.25f,
                                                      -0.25f, -0.25f))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(u_1, xi_23),
                                              _mm256_set_ps(-0.25f, -0.25f,
                                                            -0.25f, -0.25f,
                                                            -0.25f, -0.25f,
                                                            -0.25f, -0.25f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              xi_7,
                                              _mm256_set_ps(
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f)),
                                          _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                        rr_0, rr_0, rr_0,
                                                        rr_0))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_2, xi_20),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(_mm256_mul_ps(u_0, xi_7),
                                                _mm256_set_ps(0.125f, 0.125f,
                                                              0.125f, 0.125f,
                                                              0.125f, 0.125f,
                                                              0.125f, 0.125f)),
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_1, xi_23),
                                            _mm256_set_ps(0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f)),
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_2, xi_20),
                              _mm256_set_ps(
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f,
                                  0.041666666666666664f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_0, xi_23),
                          _mm256_set_ps(
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_1, xi_7),
                      _mm256_set_ps(
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_8 = _mm256_add_ps(
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
                                                              _mm256_mul_ps(
                                                                  xi_23,
                                                                  _mm256_set_ps(
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f)),
                                                              _mm256_mul_ps(
                                                                  xi_7,
                                                                  _mm256_set_ps(
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f))),
                                                          _mm256_mul_ps(
                                                              _mm256_mul_ps(
                                                                  u_0, xi_7),
                                                              _mm256_set_ps(
                                                                  0.25f, 0.25f,
                                                                  0.25f, 0.25f,
                                                                  0.25f, 0.25f,
                                                                  0.25f,
                                                                  0.25f))),
                                                      _mm256_mul_ps(
                                                          _mm256_mul_ps(u_1,
                                                                        xi_23),
                                                          _mm256_set_ps(
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f))),
                                                  _mm256_mul_ps(
                                                      _mm256_mul_ps(u_0, xi_23),
                                                      _mm256_set_ps(
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f))),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(u_1, xi_7),
                                                  _mm256_set_ps(
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(
                                                  xi_23,
                                                  _mm256_set_ps(
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f)),
                                              _mm256_set_ps(rr_0, rr_0, rr_0,
                                                            rr_0, rr_0, rr_0,
                                                            rr_0, rr_0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              xi_7,
                                              _mm256_set_ps(
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f)),
                                          _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                        rr_0, rr_0, rr_0,
                                                        rr_0))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_2, xi_20),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_2, xi_20),
                                      _mm256_set_ps(0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f)),
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_0, xi_7),
                                            _mm256_set_ps(-0.125f, -0.125f,
                                                          -0.125f, -0.125f,
                                                          -0.125f, -0.125f,
                                                          -0.125f, -0.125f)),
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(_mm256_mul_ps(u_1, xi_23),
                                        _mm256_set_ps(-0.125f, -0.125f, -0.125f,
                                                      -0.125f, -0.125f, -0.125f,
                                                      -0.125f, -0.125f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_0, xi_23),
                          _mm256_set_ps(
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_1, xi_7),
                      _mm256_set_ps(
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_9 = _mm256_add_ps(
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
                                                              _mm256_mul_ps(
                                                                  xi_23,
                                                                  _mm256_set_ps(
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f)),
                                                              _mm256_mul_ps(
                                                                  xi_7,
                                                                  _mm256_set_ps(
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f))),
                                                          _mm256_mul_ps(
                                                              _mm256_mul_ps(
                                                                  u_0, xi_7),
                                                              _mm256_set_ps(
                                                                  0.25f, 0.25f,
                                                                  0.25f, 0.25f,
                                                                  0.25f, 0.25f,
                                                                  0.25f,
                                                                  0.25f))),
                                                      _mm256_mul_ps(
                                                          _mm256_mul_ps(u_1,
                                                                        xi_23),
                                                          _mm256_set_ps(
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f))),
                                                  _mm256_mul_ps(
                                                      _mm256_mul_ps(
                                                          xi_23,
                                                          _mm256_set_ps(
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f)),
                                                      _mm256_set_ps(
                                                          rr_0, rr_0, rr_0,
                                                          rr_0, rr_0, rr_0,
                                                          rr_0, rr_0))),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f)),
                                                  _mm256_set_ps(
                                                      rr_0, rr_0, rr_0, rr_0,
                                                      rr_0, rr_0, rr_0, rr_0))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(u_0, xi_23),
                                              _mm256_set_ps(
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(u_1, xi_7),
                                          _mm256_set_ps(0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_2, xi_20),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_2, xi_20),
                                      _mm256_set_ps(0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f)),
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_0, xi_7),
                                            _mm256_set_ps(-0.125f, -0.125f,
                                                          -0.125f, -0.125f,
                                                          -0.125f, -0.125f,
                                                          -0.125f, -0.125f)),
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(_mm256_mul_ps(u_1, xi_23),
                                        _mm256_set_ps(-0.125f, -0.125f, -0.125f,
                                                      -0.125f, -0.125f, -0.125f,
                                                      -0.125f, -0.125f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_0, xi_23),
                          _mm256_set_ps(
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_1, xi_7),
                      _mm256_set_ps(
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_10 = _mm256_add_ps(
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
                                                              _mm256_mul_ps(
                                                                  xi_23,
                                                                  _mm256_set_ps(
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f)),
                                                              _mm256_mul_ps(
                                                                  xi_7,
                                                                  _mm256_set_ps(
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f))),
                                                          _mm256_mul_ps(
                                                              _mm256_mul_ps(
                                                                  xi_7,
                                                                  _mm256_set_ps(
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f)),
                                                              _mm256_set_ps(
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0))),
                                                      _mm256_mul_ps(
                                                          _mm256_mul_ps(u_0,
                                                                        xi_23),
                                                          _mm256_set_ps(
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f))),
                                                  _mm256_mul_ps(
                                                      _mm256_mul_ps(u_1, xi_7),
                                                      _mm256_set_ps(
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f))),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(u_0, xi_7),
                                                  _mm256_set_ps(
                                                      -0.25f, -0.25f, -0.25f,
                                                      -0.25f, -0.25f, -0.25f,
                                                      -0.25f, -0.25f))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(u_1, xi_23),
                                              _mm256_set_ps(-0.25f, -0.25f,
                                                            -0.25f, -0.25f,
                                                            -0.25f, -0.25f,
                                                            -0.25f, -0.25f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              xi_23,
                                              _mm256_set_ps(
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f)),
                                          _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                        rr_0, rr_0, rr_0,
                                                        rr_0))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_2, xi_20),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(_mm256_mul_ps(u_0, xi_7),
                                                _mm256_set_ps(0.125f, 0.125f,
                                                              0.125f, 0.125f,
                                                              0.125f, 0.125f,
                                                              0.125f, 0.125f)),
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_1, xi_23),
                                            _mm256_set_ps(0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f)),
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_2, xi_20),
                              _mm256_set_ps(
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f,
                                  0.041666666666666664f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_0, xi_23),
                          _mm256_set_ps(
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_1, xi_7),
                      _mm256_set_ps(
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_11 = _mm256_add_ps(
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
                                                              _mm256_mul_ps(
                                                                  xi_20,
                                                                  _mm256_set_ps(
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f)),
                                                              _mm256_mul_ps(
                                                                  xi_7,
                                                                  _mm256_set_ps(
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f))),
                                                          _mm256_mul_ps(
                                                              _mm256_mul_ps(
                                                                  u_1, xi_20),
                                                              _mm256_set_ps(
                                                                  0.25f, 0.25f,
                                                                  0.25f, 0.25f,
                                                                  0.25f, 0.25f,
                                                                  0.25f,
                                                                  0.25f))),
                                                      _mm256_mul_ps(
                                                          _mm256_mul_ps(u_2,
                                                                        xi_7),
                                                          _mm256_set_ps(
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f))),
                                                  _mm256_mul_ps(
                                                      _mm256_mul_ps(u_1, xi_7),
                                                      _mm256_set_ps(
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f))),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(u_2, xi_20),
                                                  _mm256_set_ps(
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(
                                                  xi_20,
                                                  _mm256_set_ps(
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f)),
                                              _mm256_set_ps(rr_0, rr_0, rr_0,
                                                            rr_0, rr_0, rr_0,
                                                            rr_0, rr_0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              xi_7,
                                              _mm256_set_ps(
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f)),
                                          _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                        rr_0, rr_0, rr_0,
                                                        rr_0))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_0, xi_23),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_0, xi_23),
                                      _mm256_set_ps(0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f)),
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_1, xi_20),
                                            _mm256_set_ps(-0.125f, -0.125f,
                                                          -0.125f, -0.125f,
                                                          -0.125f, -0.125f,
                                                          -0.125f, -0.125f)),
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(_mm256_mul_ps(u_2, xi_7),
                                        _mm256_set_ps(-0.125f, -0.125f, -0.125f,
                                                      -0.125f, -0.125f, -0.125f,
                                                      -0.125f, -0.125f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_1, xi_7),
                          _mm256_set_ps(
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_2, xi_20),
                      _mm256_set_ps(
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_12 = _mm256_add_ps(
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
                                                              _mm256_mul_ps(
                                                                  xi_20,
                                                                  _mm256_set_ps(
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f)),
                                                              _mm256_mul_ps(
                                                                  xi_7,
                                                                  _mm256_set_ps(
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f))),
                                                          _mm256_mul_ps(
                                                              _mm256_mul_ps(
                                                                  xi_7,
                                                                  _mm256_set_ps(
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f)),
                                                              _mm256_set_ps(
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0))),
                                                      _mm256_mul_ps(
                                                          _mm256_mul_ps(u_1,
                                                                        xi_7),
                                                          _mm256_set_ps(
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f))),
                                                  _mm256_mul_ps(
                                                      _mm256_mul_ps(u_2, xi_20),
                                                      _mm256_set_ps(
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f))),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(u_1, xi_20),
                                                  _mm256_set_ps(
                                                      -0.25f, -0.25f, -0.25f,
                                                      -0.25f, -0.25f, -0.25f,
                                                      -0.25f, -0.25f))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(u_2, xi_7),
                                              _mm256_set_ps(-0.25f, -0.25f,
                                                            -0.25f, -0.25f,
                                                            -0.25f, -0.25f,
                                                            -0.25f, -0.25f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              xi_20,
                                              _mm256_set_ps(
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f)),
                                          _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                        rr_0, rr_0, rr_0,
                                                        rr_0))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_0, xi_23),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(_mm256_mul_ps(u_1, xi_20),
                                                _mm256_set_ps(0.125f, 0.125f,
                                                              0.125f, 0.125f,
                                                              0.125f, 0.125f,
                                                              0.125f, 0.125f)),
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_2, xi_7),
                                            _mm256_set_ps(0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f)),
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_0, xi_23),
                              _mm256_set_ps(
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f,
                                  0.041666666666666664f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_1, xi_7),
                          _mm256_set_ps(
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_2, xi_20),
                      _mm256_set_ps(
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_13 = _mm256_add_ps(
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
                                                              _mm256_mul_ps(
                                                                  xi_20,
                                                                  _mm256_set_ps(
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f)),
                                                              _mm256_mul_ps(
                                                                  xi_23,
                                                                  _mm256_set_ps(
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f))),
                                                          _mm256_mul_ps(
                                                              _mm256_mul_ps(
                                                                  xi_23,
                                                                  _mm256_set_ps(
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f)),
                                                              _mm256_set_ps(
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0))),
                                                      _mm256_mul_ps(
                                                          _mm256_mul_ps(u_0,
                                                                        xi_23),
                                                          _mm256_set_ps(
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f))),
                                                  _mm256_mul_ps(
                                                      _mm256_mul_ps(u_2, xi_20),
                                                      _mm256_set_ps(
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f))),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(u_0, xi_20),
                                                  _mm256_set_ps(
                                                      -0.25f, -0.25f, -0.25f,
                                                      -0.25f, -0.25f, -0.25f,
                                                      -0.25f, -0.25f))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(u_2, xi_23),
                                              _mm256_set_ps(-0.25f, -0.25f,
                                                            -0.25f, -0.25f,
                                                            -0.25f, -0.25f,
                                                            -0.25f, -0.25f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              xi_20,
                                              _mm256_set_ps(
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f)),
                                          _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                        rr_0, rr_0, rr_0,
                                                        rr_0))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_1, xi_7),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(_mm256_mul_ps(u_0, xi_20),
                                                _mm256_set_ps(0.125f, 0.125f,
                                                              0.125f, 0.125f,
                                                              0.125f, 0.125f,
                                                              0.125f, 0.125f)),
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_2, xi_23),
                                            _mm256_set_ps(0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f)),
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_1, xi_7),
                              _mm256_set_ps(
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f,
                                  0.041666666666666664f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_0, xi_23),
                          _mm256_set_ps(
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_2, xi_20),
                      _mm256_set_ps(
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_14 = _mm256_add_ps(
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
                                                              _mm256_mul_ps(
                                                                  xi_20,
                                                                  _mm256_set_ps(
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f)),
                                                              _mm256_mul_ps(
                                                                  xi_23,
                                                                  _mm256_set_ps(
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f))),
                                                          _mm256_mul_ps(
                                                              _mm256_mul_ps(
                                                                  u_0, xi_20),
                                                              _mm256_set_ps(
                                                                  0.25f, 0.25f,
                                                                  0.25f, 0.25f,
                                                                  0.25f, 0.25f,
                                                                  0.25f,
                                                                  0.25f))),
                                                      _mm256_mul_ps(
                                                          _mm256_mul_ps(u_2,
                                                                        xi_23),
                                                          _mm256_set_ps(
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f))),
                                                  _mm256_mul_ps(
                                                      _mm256_mul_ps(u_0, xi_23),
                                                      _mm256_set_ps(
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f))),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(u_2, xi_20),
                                                  _mm256_set_ps(
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(
                                                  xi_20,
                                                  _mm256_set_ps(
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f,
                                                      -0.041666666666666664f)),
                                              _mm256_set_ps(rr_0, rr_0, rr_0,
                                                            rr_0, rr_0, rr_0,
                                                            rr_0, rr_0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              xi_23,
                                              _mm256_set_ps(
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f)),
                                          _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                        rr_0, rr_0, rr_0,
                                                        rr_0))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_1, xi_7),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_1, xi_7),
                                      _mm256_set_ps(0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f)),
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_0, xi_20),
                                            _mm256_set_ps(-0.125f, -0.125f,
                                                          -0.125f, -0.125f,
                                                          -0.125f, -0.125f,
                                                          -0.125f, -0.125f)),
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(_mm256_mul_ps(u_2, xi_23),
                                        _mm256_set_ps(-0.125f, -0.125f, -0.125f,
                                                      -0.125f, -0.125f, -0.125f,
                                                      -0.125f, -0.125f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_0, xi_23),
                          _mm256_set_ps(
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_2, xi_20),
                      _mm256_set_ps(
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_15 = _mm256_add_ps(
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
                                                              _mm256_mul_ps(
                                                                  xi_7,
                                                                  _mm256_set_ps(
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f)),
                                                              _mm256_mul_ps(
                                                                  xi_20,
                                                                  _mm256_set_ps(
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f))),
                                                          _mm256_mul_ps(
                                                              _mm256_mul_ps(
                                                                  xi_20,
                                                                  _mm256_set_ps(
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f)),
                                                              _mm256_set_ps(
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0))),
                                                      _mm256_mul_ps(
                                                          _mm256_mul_ps(u_1,
                                                                        xi_7),
                                                          _mm256_set_ps(
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f))),
                                                  _mm256_mul_ps(
                                                      _mm256_mul_ps(u_2, xi_20),
                                                      _mm256_set_ps(
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f))),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(u_1, xi_20),
                                                  _mm256_set_ps(
                                                      -0.25f, -0.25f, -0.25f,
                                                      -0.25f, -0.25f, -0.25f,
                                                      -0.25f, -0.25f))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(u_2, xi_7),
                                              _mm256_set_ps(-0.25f, -0.25f,
                                                            -0.25f, -0.25f,
                                                            -0.25f, -0.25f,
                                                            -0.25f, -0.25f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              xi_7,
                                              _mm256_set_ps(
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f)),
                                          _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                        rr_0, rr_0, rr_0,
                                                        rr_0))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_0, xi_23),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(_mm256_mul_ps(u_1, xi_20),
                                                _mm256_set_ps(0.125f, 0.125f,
                                                              0.125f, 0.125f,
                                                              0.125f, 0.125f,
                                                              0.125f, 0.125f)),
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_2, xi_7),
                                            _mm256_set_ps(0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f)),
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_0, xi_23),
                              _mm256_set_ps(
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f,
                                  0.041666666666666664f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_1, xi_7),
                          _mm256_set_ps(
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_2, xi_20),
                      _mm256_set_ps(
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_16 = _mm256_add_ps(
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
                                                              _mm256_mul_ps(
                                                                  xi_20,
                                                                  _mm256_set_ps(
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f)),
                                                              _mm256_mul_ps(
                                                                  xi_7,
                                                                  _mm256_set_ps(
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f))),
                                                          _mm256_mul_ps(
                                                              _mm256_mul_ps(
                                                                  u_1, xi_20),
                                                              _mm256_set_ps(
                                                                  0.25f, 0.25f,
                                                                  0.25f, 0.25f,
                                                                  0.25f, 0.25f,
                                                                  0.25f,
                                                                  0.25f))),
                                                      _mm256_mul_ps(
                                                          _mm256_mul_ps(u_2,
                                                                        xi_7),
                                                          _mm256_set_ps(
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f))),
                                                  _mm256_mul_ps(
                                                      _mm256_mul_ps(
                                                          xi_20,
                                                          _mm256_set_ps(
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f)),
                                                      _mm256_set_ps(
                                                          rr_0, rr_0, rr_0,
                                                          rr_0, rr_0, rr_0,
                                                          rr_0, rr_0))),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f)),
                                                  _mm256_set_ps(
                                                      rr_0, rr_0, rr_0, rr_0,
                                                      rr_0, rr_0, rr_0, rr_0))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(u_1, xi_7),
                                              _mm256_set_ps(
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(u_2, xi_20),
                                          _mm256_set_ps(0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_0, xi_23),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_0, xi_23),
                                      _mm256_set_ps(0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f)),
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_1, xi_20),
                                            _mm256_set_ps(-0.125f, -0.125f,
                                                          -0.125f, -0.125f,
                                                          -0.125f, -0.125f,
                                                          -0.125f, -0.125f)),
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(_mm256_mul_ps(u_2, xi_7),
                                        _mm256_set_ps(-0.125f, -0.125f, -0.125f,
                                                      -0.125f, -0.125f, -0.125f,
                                                      -0.125f, -0.125f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_1, xi_7),
                          _mm256_set_ps(
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_2, xi_20),
                      _mm256_set_ps(
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_17 = _mm256_add_ps(
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
                                                              _mm256_mul_ps(
                                                                  xi_20,
                                                                  _mm256_set_ps(
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f)),
                                                              _mm256_mul_ps(
                                                                  xi_23,
                                                                  _mm256_set_ps(
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f))),
                                                          _mm256_mul_ps(
                                                              _mm256_mul_ps(
                                                                  u_0, xi_20),
                                                              _mm256_set_ps(
                                                                  0.25f, 0.25f,
                                                                  0.25f, 0.25f,
                                                                  0.25f, 0.25f,
                                                                  0.25f,
                                                                  0.25f))),
                                                      _mm256_mul_ps(
                                                          _mm256_mul_ps(u_2,
                                                                        xi_23),
                                                          _mm256_set_ps(
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f,
                                                              0.25f, 0.25f))),
                                                  _mm256_mul_ps(
                                                      _mm256_mul_ps(
                                                          xi_20,
                                                          _mm256_set_ps(
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f,
                                                              0.041666666666666664f)),
                                                      _mm256_set_ps(
                                                          rr_0, rr_0, rr_0,
                                                          rr_0, rr_0, rr_0,
                                                          rr_0, rr_0))),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(
                                                      xi_23,
                                                      _mm256_set_ps(
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f)),
                                                  _mm256_set_ps(
                                                      rr_0, rr_0, rr_0, rr_0,
                                                      rr_0, rr_0, rr_0, rr_0))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(u_0, xi_23),
                                              _mm256_set_ps(
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(u_2, xi_20),
                                          _mm256_set_ps(0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_1, xi_7),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_1, xi_7),
                                      _mm256_set_ps(0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f,
                                                    0.041666666666666664f)),
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_0, xi_20),
                                            _mm256_set_ps(-0.125f, -0.125f,
                                                          -0.125f, -0.125f,
                                                          -0.125f, -0.125f,
                                                          -0.125f, -0.125f)),
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(_mm256_mul_ps(u_2, xi_23),
                                        _mm256_set_ps(-0.125f, -0.125f, -0.125f,
                                                      -0.125f, -0.125f, -0.125f,
                                                      -0.125f, -0.125f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_0, xi_23),
                          _mm256_set_ps(
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_2, xi_20),
                      _mm256_set_ps(
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 forceTerm_18 = _mm256_add_ps(
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
                                                              _mm256_mul_ps(
                                                                  xi_23,
                                                                  _mm256_set_ps(
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f,
                                                                      0.083333333333333329f)),
                                                              _mm256_mul_ps(
                                                                  xi_20,
                                                                  _mm256_set_ps(
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f,
                                                                      -0.083333333333333329f))),
                                                          _mm256_mul_ps(
                                                              _mm256_mul_ps(
                                                                  xi_20,
                                                                  _mm256_set_ps(
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f,
                                                                      0.041666666666666664f)),
                                                              _mm256_set_ps(
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0,
                                                                  rr_0, rr_0))),
                                                      _mm256_mul_ps(
                                                          _mm256_mul_ps(u_0,
                                                                        xi_23),
                                                          _mm256_set_ps(
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f,
                                                              0.16666666666666666f))),
                                                  _mm256_mul_ps(
                                                      _mm256_mul_ps(u_2, xi_20),
                                                      _mm256_set_ps(
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f,
                                                          0.16666666666666666f))),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(u_0, xi_20),
                                                  _mm256_set_ps(
                                                      -0.25f, -0.25f, -0.25f,
                                                      -0.25f, -0.25f, -0.25f,
                                                      -0.25f, -0.25f))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(u_2, xi_23),
                                              _mm256_set_ps(-0.25f, -0.25f,
                                                            -0.25f, -0.25f,
                                                            -0.25f, -0.25f,
                                                            -0.25f, -0.25f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              xi_23,
                                              _mm256_set_ps(
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f,
                                                  -0.041666666666666664f)),
                                          _mm256_set_ps(rr_0, rr_0, rr_0, rr_0,
                                                        rr_0, rr_0, rr_0,
                                                        rr_0))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(u_1, xi_7),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(_mm256_mul_ps(u_0, xi_20),
                                                _mm256_set_ps(0.125f, 0.125f,
                                                              0.125f, 0.125f,
                                                              0.125f, 0.125f,
                                                              0.125f, 0.125f)),
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear))),
                          _mm256_mul_ps(
                              _mm256_mul_ps(_mm256_mul_ps(u_2, xi_23),
                                            _mm256_set_ps(0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f)),
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear))),
                      _mm256_mul_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(u_1, xi_7),
                              _mm256_set_ps(
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f, 0.041666666666666664f,
                                  0.041666666666666664f,
                                  0.041666666666666664f)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(u_0, xi_23),
                          _mm256_set_ps(
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f,
                              -0.083333333333333329f, -0.083333333333333329f)),
                      _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear, omega_shear,
                                    omega_shear, omega_shear))),
              _mm256_mul_ps(
                  _mm256_mul_ps(
                      _mm256_mul_ps(u_2, xi_20),
                      _mm256_set_ps(
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f,
                          -0.083333333333333329f, -0.083333333333333329f)),
                  _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear, omega_shear,
                                omega_shear, omega_shear)));
          const __m256 u0Mu1 = _mm256_add_ps(
              _mm256_mul_ps(u_1, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f,
                                               -1.0f, -1.0f, -1.0f, -1.0f)),
              u_0);
          const __m256 u0Pu1 = _mm256_add_ps(u_0, u_1);
          const __m256 u1Pu2 = _mm256_add_ps(u_1, u_2);
          const __m256 u1Mu2 = _mm256_add_ps(
              _mm256_mul_ps(u_2, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f,
                                               -1.0f, -1.0f, -1.0f, -1.0f)),
              u_1);
          const __m256 u0Mu2 = _mm256_add_ps(
              _mm256_mul_ps(u_2, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f,
                                               -1.0f, -1.0f, -1.0f, -1.0f)),
              u_0);
          const __m256 u0Pu2 = _mm256_add_ps(u_0, u_2);
          const __m256 f_eq_common = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_mul_ps(rho, _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                           -1.0f, -1.0f, -1.0f,
                                                           -1.0f, -1.0f)),
                          _mm256_mul_ps(u_0, u_0)),
                      _mm256_mul_ps(
                          _mm256_mul_ps(rho, _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                           -1.0f, -1.0f, -1.0f,
                                                           -1.0f, -1.0f)),
                          _mm256_mul_ps(u_1, u_1))),
                  _mm256_mul_ps(
                      _mm256_mul_ps(rho,
                                    _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f,
                                                  -1.0f, -1.0f, -1.0f, -1.0f)),
                      _mm256_mul_ps(u_2, u_2))),
              rho);
          _mm256_store_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_mul_ps(
                                  f_eq_common,
                                  _mm256_set_ps(0.33333333333333331f,
                                                0.33333333333333331f,
                                                0.33333333333333331f,
                                                0.33333333333333331f,
                                                0.33333333333333331f,
                                                0.33333333333333331f,
                                                0.33333333333333331f,
                                                0.33333333333333331f)),
                              _mm256_mul_ps(xi_13,
                                            _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                          -1.0f, -1.0f, -1.0f,
                                                          -1.0f, -1.0f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear)),
                      forceTerm_0),
                  xi_13));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_mul_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              xi_18,
                                              _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                            0.5f, 0.5f, 0.5f,
                                                            0.5f, 0.5f)),
                                          _mm256_mul_ps(
                                              xi_3,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(rho, u_1),
                                          _mm256_set_ps(0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f,
                                                        0.16666666666666666f))),
                                  _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0,
                                                rr_0, rr_0, rr_0)),
                              _mm256_mul_ps(
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear),
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  f_eq_common,
                                                  _mm256_set_ps(
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f)),
                                              _mm256_mul_ps(
                                                  xi_18,
                                                  _mm256_set_ps(-0.5f, -0.5f,
                                                                -0.5f, -0.5f,
                                                                -0.5f, -0.5f,
                                                                -0.5f, -0.5f))),
                                          _mm256_mul_ps(
                                              xi_3,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          rho,
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(
                                                      0.33333333333333331f,
                                                      0.33333333333333331f,
                                                      0.33333333333333331f,
                                                      0.33333333333333331f,
                                                      0.33333333333333331f,
                                                      0.33333333333333331f,
                                                      0.33333333333333331f,
                                                      0.33333333333333331f),
                                                  _mm256_mul_ps(u_1, u_1)),
                                              _mm256_set_ps(
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f)))))),
                          _mm256_blendv_ps(
                              _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                            0.0f, 0.0f),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho,
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  u_0,
                                                  _mm256_set_ps(
                                                      2.0f, 2.0f, 2.0f, 2.0f,
                                                      2.0f, 2.0f, 2.0f, 2.0f)),
                                              _mm256_set_ps(v_s, v_s, v_s, v_s,
                                                            v_s, v_s, v_s,
                                                            v_s))),
                                      _mm256_set_ps(0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f)),
                                  _mm256_set_ps(v_s, v_s, v_s, v_s, v_s, v_s,
                                                v_s, v_s)),
                              _mm256_cmp_ps(
                                  _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f,
                                                -1.0f, -1.0f, -1.0f, -1.0f),
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f),
                                          _mm256_set_ps(grid_size, grid_size,
                                                        grid_size, grid_size,
                                                        grid_size, grid_size,
                                                        grid_size, grid_size)),
                                      _mm256_set_ps(
                                          ((float)(ctr_1)), ((float)(ctr_1)),
                                          ((float)(ctr_1)), ((float)(ctr_1)),
                                          ((float)(ctr_1)), ((float)(ctr_1)),
                                          ((float)(ctr_1)), ((float)(ctr_1)))),
                                  _CMP_LE_OQ))),
                      forceTerm_1),
                  xi_3));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          2 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_mul_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              xi_3,
                                              _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                            0.5f, 0.5f, 0.5f,
                                                            0.5f, 0.5f)),
                                          _mm256_mul_ps(
                                              xi_18,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(rho, u_1),
                                          _mm256_set_ps(
                                              -0.16666666666666666f,
                                              -0.16666666666666666f,
                                              -0.16666666666666666f,
                                              -0.16666666666666666f,
                                              -0.16666666666666666f,
                                              -0.16666666666666666f,
                                              -0.16666666666666666f,
                                              -0.16666666666666666f))),
                                  _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0,
                                                rr_0, rr_0, rr_0)),
                              _mm256_mul_ps(
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear),
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  f_eq_common,
                                                  _mm256_set_ps(
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f,
                                                      0.16666666666666666f)),
                                              _mm256_mul_ps(
                                                  xi_18,
                                                  _mm256_set_ps(-0.5f, -0.5f,
                                                                -0.5f, -0.5f,
                                                                -0.5f, -0.5f,
                                                                -0.5f, -0.5f))),
                                          _mm256_mul_ps(
                                              xi_3,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          rho,
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(
                                                      0.33333333333333331f,
                                                      0.33333333333333331f,
                                                      0.33333333333333331f,
                                                      0.33333333333333331f,
                                                      0.33333333333333331f,
                                                      0.33333333333333331f,
                                                      0.33333333333333331f,
                                                      0.33333333333333331f),
                                                  _mm256_mul_ps(u_1, u_1)),
                                              _mm256_set_ps(
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f,
                                                  -0.1111111111111111f)))))),
                          _mm256_blendv_ps(
                              _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                            0.0f, 0.0f),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho,
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  u_0, _mm256_set_ps(
                                                           -2.0f, -2.0f, -2.0f,
                                                           -2.0f, -2.0f, -2.0f,
                                                           -2.0f, -2.0f)),
                                              _mm256_set_ps(v_s, v_s, v_s, v_s,
                                                            v_s, v_s, v_s,
                                                            v_s))),
                                      _mm256_set_ps(0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f)),
                                  _mm256_set_ps(v_s, v_s, v_s, v_s, v_s, v_s,
                                                v_s, v_s)),
                              _mm256_cmp_ps(
                                  _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                                0.0f, 0.0f, 0.0f),
                                  _mm256_set_ps(
                                      ((float)(ctr_1)), ((float)(ctr_1)),
                                      ((float)(ctr_1)), ((float)(ctr_1)),
                                      ((float)(ctr_1)), ((float)(ctr_1)),
                                      ((float)(ctr_1)), ((float)(ctr_1))),
                                  _CMP_GE_OQ))),
                      forceTerm_2),
                  xi_18));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          3 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_24, _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f)),
                                      _mm256_mul_ps(
                                          xi_5,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(rho, u_0),
                                      _mm256_set_ps(-0.16666666666666666f,
                                                    -0.16666666666666666f,
                                                    -0.16666666666666666f,
                                                    -0.16666666666666666f,
                                                    -0.16666666666666666f,
                                                    -0.16666666666666666f,
                                                    -0.16666666666666666f,
                                                    -0.16666666666666666f))),
                              _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0,
                                            rr_0, rr_0)),
                          _mm256_mul_ps(
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear),
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              f_eq_common,
                                              _mm256_set_ps(
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f)),
                                          _mm256_mul_ps(
                                              xi_24,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          xi_5,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   _mm256_set_ps(
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f),
                                                   _mm256_mul_ps(u_0, u_0)),
                                               _mm256_set_ps(
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f)))))),
                      forceTerm_3),
                  xi_5));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          4 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_5, _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                              0.5f, 0.5f, 0.5f,
                                                              0.5f, 0.5f)),
                                      _mm256_mul_ps(
                                          xi_24,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(rho, u_0),
                                      _mm256_set_ps(0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f))),
                              _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0,
                                            rr_0, rr_0)),
                          _mm256_mul_ps(
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear),
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              f_eq_common,
                                              _mm256_set_ps(
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f)),
                                          _mm256_mul_ps(
                                              xi_24,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          xi_5,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   _mm256_set_ps(
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f),
                                                   _mm256_mul_ps(u_0, u_0)),
                                               _mm256_set_ps(
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f)))))),
                      forceTerm_4),
                  xi_24));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          5 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_16, _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f)),
                                      _mm256_mul_ps(
                                          xi_22,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(rho, u_2),
                                      _mm256_set_ps(0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f,
                                                    0.16666666666666666f))),
                              _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0,
                                            rr_0, rr_0)),
                          _mm256_mul_ps(
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear),
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              f_eq_common,
                                              _mm256_set_ps(
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f)),
                                          _mm256_mul_ps(
                                              xi_16,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          xi_22,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   _mm256_set_ps(
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f),
                                                   _mm256_mul_ps(u_2, u_2)),
                                               _mm256_set_ps(
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f)))))),
                      forceTerm_5),
                  xi_22));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          6 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_22, _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f)),
                                      _mm256_mul_ps(
                                          xi_16,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(rho, u_2),
                                      _mm256_set_ps(-0.16666666666666666f,
                                                    -0.16666666666666666f,
                                                    -0.16666666666666666f,
                                                    -0.16666666666666666f,
                                                    -0.16666666666666666f,
                                                    -0.16666666666666666f,
                                                    -0.16666666666666666f,
                                                    -0.16666666666666666f))),
                              _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0,
                                            rr_0, rr_0)),
                          _mm256_mul_ps(
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear),
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              f_eq_common,
                                              _mm256_set_ps(
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f,
                                                  0.16666666666666666f)),
                                          _mm256_mul_ps(
                                              xi_16,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          xi_22,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      rho, _mm256_add_ps(
                                               _mm256_mul_ps(
                                                   _mm256_set_ps(
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f,
                                                       0.33333333333333331f),
                                                   _mm256_mul_ps(u_2, u_2)),
                                               _mm256_set_ps(
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f,
                                                   -0.1111111111111111f)))))),
                      forceTerm_6),
                  xi_16));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          7 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_mul_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              xi_6,
                                              _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                            0.5f, 0.5f, 0.5f,
                                                            0.5f, 0.5f)),
                                          _mm256_mul_ps(
                                              xi_21,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(rho, u0Mu1),
                                          _mm256_set_ps(
                                              -0.083333333333333329f,
                                              -0.083333333333333329f,
                                              -0.083333333333333329f,
                                              -0.083333333333333329f,
                                              -0.083333333333333329f,
                                              -0.083333333333333329f,
                                              -0.083333333333333329f,
                                              -0.083333333333333329f))),
                                  _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0,
                                                rr_0, rr_0, rr_0)),
                              _mm256_mul_ps(
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear),
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  f_eq_common,
                                                  _mm256_set_ps(
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f)),
                                              _mm256_mul_ps(
                                                  xi_21,
                                                  _mm256_set_ps(-0.5f, -0.5f,
                                                                -0.5f, -0.5f,
                                                                -0.5f, -0.5f,
                                                                -0.5f, -0.5f))),
                                          _mm256_mul_ps(
                                              xi_6,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          rho,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      _mm256_set_ps(
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f),
                                                      _mm256_mul_ps(u0Mu1,
                                                                    u0Mu1)),
                                                  _mm256_mul_ps(
                                                      _mm256_set_ps(
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f),
                                                      _mm256_mul_ps(u_2, u_2))),
                                              _mm256_set_ps(
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f)))))),
                          _mm256_blendv_ps(
                              _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                            0.0f, 0.0f),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_add_ps(
                                                      _mm256_mul_ps(
                                                          _mm256_set_ps(
                                                              -1.0f, -1.0f,
                                                              -1.0f, -1.0f,
                                                              -1.0f, -1.0f,
                                                              -1.0f, -1.0f),
                                                          _mm256_set_ps(
                                                              v_s, v_s, v_s,
                                                              v_s, v_s, v_s,
                                                              v_s, v_s)),
                                                      _mm256_mul_ps(
                                                          u_1,
                                                          _mm256_set_ps(
                                                              3.0f, 3.0f, 3.0f,
                                                              3.0f, 3.0f, 3.0f,
                                                              3.0f, 3.0f))),
                                                  _mm256_mul_ps(
                                                      u_0,
                                                      _mm256_set_ps(
                                                          -2.0f, -2.0f, -2.0f,
                                                          -2.0f, -2.0f, -2.0f,
                                                          -2.0f, -2.0f))),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f)),
                                  _mm256_set_ps(v_s, v_s, v_s, v_s, v_s, v_s,
                                                v_s, v_s)),
                              _mm256_cmp_ps(
                                  _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f,
                                                -1.0f, -1.0f, -1.0f, -1.0f),
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f),
                                          _mm256_set_ps(grid_size, grid_size,
                                                        grid_size, grid_size,
                                                        grid_size, grid_size,
                                                        grid_size, grid_size)),
                                      _mm256_set_ps(
                                          ((float)(ctr_1)), ((float)(ctr_1)),
                                          ((float)(ctr_1)), ((float)(ctr_1)),
                                          ((float)(ctr_1)), ((float)(ctr_1)),
                                          ((float)(ctr_1)), ((float)(ctr_1)))),
                                  _CMP_LE_OQ))),
                      forceTerm_7),
                  xi_21));
          _mm256_store_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          8 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_mul_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              xi_11,
                                              _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                            0.5f, 0.5f, 0.5f,
                                                            0.5f, 0.5f)),
                                          _mm256_mul_ps(
                                              xi_4,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(rho, u0Pu1),
                                          _mm256_set_ps(
                                              0.083333333333333329f,
                                              0.083333333333333329f,
                                              0.083333333333333329f,
                                              0.083333333333333329f,
                                              0.083333333333333329f,
                                              0.083333333333333329f,
                                              0.083333333333333329f,
                                              0.083333333333333329f))),
                                  _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0,
                                                rr_0, rr_0, rr_0)),
                              _mm256_mul_ps(
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear),
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  f_eq_common,
                                                  _mm256_set_ps(
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f)),
                                              _mm256_mul_ps(
                                                  xi_11,
                                                  _mm256_set_ps(-0.5f, -0.5f,
                                                                -0.5f, -0.5f,
                                                                -0.5f, -0.5f,
                                                                -0.5f, -0.5f))),
                                          _mm256_mul_ps(
                                              xi_4,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          rho,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      _mm256_set_ps(
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f),
                                                      _mm256_mul_ps(u0Pu1,
                                                                    u0Pu1)),
                                                  _mm256_mul_ps(
                                                      _mm256_set_ps(
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f),
                                                      _mm256_mul_ps(u_2, u_2))),
                                              _mm256_set_ps(
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f)))))),
                          _mm256_blendv_ps(
                              _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                            0.0f, 0.0f),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_add_ps(
                                                      _mm256_mul_ps(
                                                          u_0,
                                                          _mm256_set_ps(
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f)),
                                                      _mm256_mul_ps(
                                                          u_1,
                                                          _mm256_set_ps(
                                                              3.0f, 3.0f, 3.0f,
                                                              3.0f, 3.0f, 3.0f,
                                                              3.0f, 3.0f))),
                                                  _mm256_set_ps(
                                                      1.0f, 1.0f, 1.0f, 1.0f,
                                                      1.0f, 1.0f, 1.0f, 1.0f)),
                                              _mm256_set_ps(v_s, v_s, v_s, v_s,
                                                            v_s, v_s, v_s,
                                                            v_s))),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f)),
                                  _mm256_set_ps(v_s, v_s, v_s, v_s, v_s, v_s,
                                                v_s, v_s)),
                              _mm256_cmp_ps(
                                  _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f,
                                                -1.0f, -1.0f, -1.0f, -1.0f),
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f),
                                          _mm256_set_ps(grid_size, grid_size,
                                                        grid_size, grid_size,
                                                        grid_size, grid_size,
                                                        grid_size, grid_size)),
                                      _mm256_set_ps(
                                          ((float)(ctr_1)), ((float)(ctr_1)),
                                          ((float)(ctr_1)), ((float)(ctr_1)),
                                          ((float)(ctr_1)), ((float)(ctr_1)),
                                          ((float)(ctr_1)), ((float)(ctr_1)))),
                                  _CMP_LE_OQ))),
                      forceTerm_8),
                  xi_4));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          9 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_mul_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              xi_4,
                                              _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                            0.5f, 0.5f, 0.5f,
                                                            0.5f, 0.5f)),
                                          _mm256_mul_ps(
                                              xi_11,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(rho, u0Pu1),
                                          _mm256_set_ps(
                                              -0.083333333333333329f,
                                              -0.083333333333333329f,
                                              -0.083333333333333329f,
                                              -0.083333333333333329f,
                                              -0.083333333333333329f,
                                              -0.083333333333333329f,
                                              -0.083333333333333329f,
                                              -0.083333333333333329f))),
                                  _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0,
                                                rr_0, rr_0, rr_0)),
                              _mm256_mul_ps(
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear),
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  f_eq_common,
                                                  _mm256_set_ps(
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f)),
                                              _mm256_mul_ps(
                                                  xi_11,
                                                  _mm256_set_ps(-0.5f, -0.5f,
                                                                -0.5f, -0.5f,
                                                                -0.5f, -0.5f,
                                                                -0.5f, -0.5f))),
                                          _mm256_mul_ps(
                                              xi_4,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          rho,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      _mm256_set_ps(
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f),
                                                      _mm256_mul_ps(u0Pu1,
                                                                    u0Pu1)),
                                                  _mm256_mul_ps(
                                                      _mm256_set_ps(
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f),
                                                      _mm256_mul_ps(u_2, u_2))),
                                              _mm256_set_ps(
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f)))))),
                          _mm256_blendv_ps(
                              _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                            0.0f, 0.0f),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_add_ps(
                                                      _mm256_mul_ps(
                                                          _mm256_set_ps(
                                                              -1.0f, -1.0f,
                                                              -1.0f, -1.0f,
                                                              -1.0f, -1.0f,
                                                              -1.0f, -1.0f),
                                                          _mm256_set_ps(
                                                              v_s, v_s, v_s,
                                                              v_s, v_s, v_s,
                                                              v_s, v_s)),
                                                      _mm256_mul_ps(
                                                          u_0,
                                                          _mm256_set_ps(
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f))),
                                                  _mm256_mul_ps(
                                                      u_1, _mm256_set_ps(
                                                               3.0f, 3.0f, 3.0f,
                                                               3.0f, 3.0f, 3.0f,
                                                               3.0f, 3.0f))),
                                              _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f))),
                                      _mm256_set_ps(0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f)),
                                  _mm256_set_ps(v_s, v_s, v_s, v_s, v_s, v_s,
                                                v_s, v_s)),
                              _mm256_cmp_ps(
                                  _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                                0.0f, 0.0f, 0.0f),
                                  _mm256_set_ps(
                                      ((float)(ctr_1)), ((float)(ctr_1)),
                                      ((float)(ctr_1)), ((float)(ctr_1)),
                                      ((float)(ctr_1)), ((float)(ctr_1)),
                                      ((float)(ctr_1)), ((float)(ctr_1))),
                                  _CMP_GE_OQ))),
                      forceTerm_9),
                  xi_11));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          10 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_mul_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              xi_21,
                                              _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                            0.5f, 0.5f, 0.5f,
                                                            0.5f, 0.5f)),
                                          _mm256_mul_ps(
                                              xi_6,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(rho, u0Mu1),
                                          _mm256_set_ps(
                                              0.083333333333333329f,
                                              0.083333333333333329f,
                                              0.083333333333333329f,
                                              0.083333333333333329f,
                                              0.083333333333333329f,
                                              0.083333333333333329f,
                                              0.083333333333333329f,
                                              0.083333333333333329f))),
                                  _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0,
                                                rr_0, rr_0, rr_0)),
                              _mm256_mul_ps(
                                  _mm256_set_ps(omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear,
                                                omega_shear, omega_shear),
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  f_eq_common,
                                                  _mm256_set_ps(
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f)),
                                              _mm256_mul_ps(
                                                  xi_21,
                                                  _mm256_set_ps(-0.5f, -0.5f,
                                                                -0.5f, -0.5f,
                                                                -0.5f, -0.5f,
                                                                -0.5f, -0.5f))),
                                          _mm256_mul_ps(
                                              xi_6,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          rho,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      _mm256_set_ps(
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f,
                                                          0.125f, 0.125f),
                                                      _mm256_mul_ps(u0Mu1,
                                                                    u0Mu1)),
                                                  _mm256_mul_ps(
                                                      _mm256_set_ps(
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f,
                                                          0.041666666666666664f),
                                                      _mm256_mul_ps(u_2, u_2))),
                                              _mm256_set_ps(
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f,
                                                  -0.013888888888888888f)))))),
                          _mm256_blendv_ps(
                              _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                            0.0f, 0.0f),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_add_ps(
                                                      _mm256_mul_ps(
                                                          _mm256_set_ps(
                                                              -1.0f, -1.0f,
                                                              -1.0f, -1.0f,
                                                              -1.0f, -1.0f,
                                                              -1.0f, -1.0f),
                                                          _mm256_set_ps(
                                                              v_s, v_s, v_s,
                                                              v_s, v_s, v_s,
                                                              v_s, v_s)),
                                                      _mm256_mul_ps(
                                                          u_0,
                                                          _mm256_set_ps(
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f))),
                                                  _mm256_mul_ps(
                                                      u_1,
                                                      _mm256_set_ps(
                                                          -3.0f, -3.0f, -3.0f,
                                                          -3.0f, -3.0f, -3.0f,
                                                          -3.0f, -3.0f))),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f)),
                                  _mm256_set_ps(v_s, v_s, v_s, v_s, v_s, v_s,
                                                v_s, v_s)),
                              _mm256_cmp_ps(
                                  _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                                0.0f, 0.0f, 0.0f),
                                  _mm256_set_ps(
                                      ((float)(ctr_1)), ((float)(ctr_1)),
                                      ((float)(ctr_1)), ((float)(ctr_1)),
                                      ((float)(ctr_1)), ((float)(ctr_1)),
                                      ((float)(ctr_1)), ((float)(ctr_1))),
                                  _CMP_GE_OQ))),
                      forceTerm_10),
                  xi_6));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          11 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_17, _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f)),
                                      _mm256_mul_ps(
                                          xi_12,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(rho, u1Pu2),
                                      _mm256_set_ps(0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f))),
                              _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0,
                                            rr_0, rr_0)),
                          _mm256_mul_ps(
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear),
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              f_eq_common,
                                              _mm256_set_ps(
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f)),
                                          _mm256_mul_ps(
                                              xi_12,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          xi_17,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      rho,
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f),
                                                  _mm256_mul_ps(u1Pu2, u1Pu2)),
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f),
                                                  _mm256_mul_ps(u_0, u_0))),
                                          _mm256_set_ps(
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f)))))),
                      forceTerm_11),
                  xi_12));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          12 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_14, _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f)),
                                      _mm256_mul_ps(
                                          xi_10,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(rho, u1Mu2),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0,
                                            rr_0, rr_0)),
                          _mm256_mul_ps(
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear),
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              f_eq_common,
                                              _mm256_set_ps(
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f)),
                                          _mm256_mul_ps(
                                              xi_10,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          xi_14,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      rho,
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f),
                                                  _mm256_mul_ps(u1Mu2, u1Mu2)),
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f),
                                                  _mm256_mul_ps(u_0, u_0))),
                                          _mm256_set_ps(
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f)))))),
                      forceTerm_12),
                  xi_10));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          13 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_8, _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                              0.5f, 0.5f, 0.5f,
                                                              0.5f, 0.5f)),
                                      _mm256_mul_ps(
                                          xi_19,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(rho, u0Mu2),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0,
                                            rr_0, rr_0)),
                          _mm256_mul_ps(
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear),
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              f_eq_common,
                                              _mm256_set_ps(
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f)),
                                          _mm256_mul_ps(
                                              xi_19,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          xi_8,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      rho,
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f),
                                                  _mm256_mul_ps(u0Mu2, u0Mu2)),
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f),
                                                  _mm256_mul_ps(u_1, u_1))),
                                          _mm256_set_ps(
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f)))))),
                      forceTerm_13),
                  xi_19));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          14 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_15, _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f)),
                                      _mm256_mul_ps(
                                          xi_9,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(rho, u0Pu2),
                                      _mm256_set_ps(0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f))),
                              _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0,
                                            rr_0, rr_0)),
                          _mm256_mul_ps(
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear),
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              f_eq_common,
                                              _mm256_set_ps(
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f)),
                                          _mm256_mul_ps(
                                              xi_15,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          xi_9,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      rho,
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f),
                                                  _mm256_mul_ps(u0Pu2, u0Pu2)),
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f),
                                                  _mm256_mul_ps(u_1, u_1))),
                                          _mm256_set_ps(
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f)))))),
                      forceTerm_14),
                  xi_9));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          15 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_10, _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f)),
                                      _mm256_mul_ps(
                                          xi_14,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(rho, u1Mu2),
                                      _mm256_set_ps(0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f))),
                              _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0,
                                            rr_0, rr_0)),
                          _mm256_mul_ps(
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear),
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              f_eq_common,
                                              _mm256_set_ps(
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f)),
                                          _mm256_mul_ps(
                                              xi_10,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          xi_14,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      rho,
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f),
                                                  _mm256_mul_ps(u1Mu2, u1Mu2)),
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f),
                                                  _mm256_mul_ps(u_0, u_0))),
                                          _mm256_set_ps(
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f)))))),
                      forceTerm_15),
                  xi_14));
          _mm256_store_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          16 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_12, _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f)),
                                      _mm256_mul_ps(
                                          xi_17,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(rho, u1Pu2),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0,
                                            rr_0, rr_0)),
                          _mm256_mul_ps(
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear),
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              f_eq_common,
                                              _mm256_set_ps(
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f)),
                                          _mm256_mul_ps(
                                              xi_12,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          xi_17,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      rho,
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f),
                                                  _mm256_mul_ps(u1Pu2, u1Pu2)),
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f),
                                                  _mm256_mul_ps(u_0, u_0))),
                                          _mm256_set_ps(
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f)))))),
                      forceTerm_16),
                  xi_17));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          17 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_9, _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                              0.5f, 0.5f, 0.5f,
                                                              0.5f, 0.5f)),
                                      _mm256_mul_ps(
                                          xi_15,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(rho, u0Pu2),
                                      _mm256_set_ps(-0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f,
                                                    -0.083333333333333329f))),
                              _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0,
                                            rr_0, rr_0)),
                          _mm256_mul_ps(
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear),
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              f_eq_common,
                                              _mm256_set_ps(
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f)),
                                          _mm256_mul_ps(
                                              xi_15,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          xi_9,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      rho,
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f),
                                                  _mm256_mul_ps(u0Pu2, u0Pu2)),
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f),
                                                  _mm256_mul_ps(u_1, u_1))),
                                          _mm256_set_ps(
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f)))))),
                      forceTerm_17),
                  xi_15));
          _mm256_storeu_ps(
              &_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                          18 * _stride_pdfs_3 + ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_19, _mm256_set_ps(0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f, 0.5f,
                                                               0.5f, 0.5f)),
                                      _mm256_mul_ps(
                                          xi_8,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(rho, u0Mu2),
                                      _mm256_set_ps(0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f,
                                                    0.083333333333333329f))),
                              _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0,
                                            rr_0, rr_0)),
                          _mm256_mul_ps(
                              _mm256_set_ps(omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear,
                                            omega_shear, omega_shear),
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              f_eq_common,
                                              _mm256_set_ps(
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f,
                                                  0.041666666666666664f)),
                                          _mm256_mul_ps(
                                              xi_19,
                                              _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f, -0.5f,
                                                            -0.5f, -0.5f))),
                                      _mm256_mul_ps(
                                          xi_8,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                                  _mm256_mul_ps(
                                      rho,
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f,
                                                                0.125f, 0.125f),
                                                  _mm256_mul_ps(u0Mu2, u0Mu2)),
                                              _mm256_mul_ps(
                                                  _mm256_set_ps(
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f,
                                                      0.041666666666666664f),
                                                  _mm256_mul_ps(u_1, u_1))),
                                          _mm256_set_ps(
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f,
                                              -0.013888888888888888f)))))),
                      forceTerm_18),
                  xi_8));
        }
        for (int64_t ctr_0 = (int64_t)((_size_force_0) / (8)) * (8);
             ctr_0 < _size_force_0; ctr_0 += 1) {
          const float xi_25 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         14 * _stride_pdfs_3 + ctr_0];
          const float xi_26 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         13 * _stride_pdfs_3 + ctr_0];
          const float xi_27 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         11 * _stride_pdfs_3 + ctr_0];
          const float xi_28 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         5 * _stride_pdfs_3 + ctr_0];
          const float xi_29 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         10 * _stride_pdfs_3 + ctr_0];
          const float xi_30 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         15 * _stride_pdfs_3 + ctr_0];
          const float xi_31 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         18 * _stride_pdfs_3 + ctr_0];
          const float xi_32 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         6 * _stride_pdfs_3 + ctr_0];
          const float xi_33 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         _stride_pdfs_3 + ctr_0];
          const float xi_34 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         17 * _stride_pdfs_3 + ctr_0];
          const float xi_35 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         16 * _stride_pdfs_3 + ctr_0];
          const float xi_36 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         9 * _stride_pdfs_3 + ctr_0];
          const float xi_37 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         2 * _stride_pdfs_3 + ctr_0];
          const float xi_38 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         8 * _stride_pdfs_3 + ctr_0];
          const float xi_39 = _data_pdfs[_stride_pdfs_1 * ctr_1 +
                                         _stride_pdfs_2 * ctr_2 + ctr_0];
          const float xi_40 = _data_force[_stride_force_1 * ctr_1 +
                                          _stride_force_2 * ctr_2 + ctr_0];
          const float xi_41 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         7 * _stride_pdfs_3 + ctr_0];
          const float xi_42 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         3 * _stride_pdfs_3 + ctr_0];
          const float xi_43 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         4 * _stride_pdfs_3 + ctr_0];
          const float xi_44 =
              _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 +
                          2 * _stride_force_3 + ctr_0];
          const float xi_45 =
              _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 +
                          _stride_force_3 + ctr_0];
          const float xi_46 =
              _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                         12 * _stride_pdfs_3 + ctr_0];
          const float xi_3 = xi_33;
          const float xi_4 = xi_38;
          const float xi_5 = xi_42;
          const float xi_6 = xi_29;
          const float xi_7 = xi_45;
          const float xi_8 = xi_31;
          const float xi_9 = xi_25;
          const float xi_10 = xi_46;
          const float xi_11 = xi_36;
          const float xi_12 = xi_27;
          const float xi_13 = xi_39;
          const float xi_14 = xi_30;
          const float xi_15 = xi_34;
          const float xi_16 = xi_32;
          const float xi_17 = xi_35;
          const float xi_18 = xi_37;
          const float xi_19 = xi_26;
          const float xi_20 = xi_44;
          const float xi_21 = xi_41;
          const float xi_22 = xi_28;
          const float xi_23 = xi_40;
          const float xi_24 = xi_43;
          const float vel0Term = xi_24 + xi_4 + xi_6 + xi_8 + xi_9;
          const float vel1Term = xi_12 + xi_14 + xi_21 + xi_3;
          const float vel2Term = xi_10 + xi_19 + xi_22;
          const float rho = vel0Term + vel1Term + vel2Term + xi_11 + xi_13 +
                            xi_15 + xi_16 + xi_17 + xi_18 + xi_5;
          const float xi_1 = ((1.0f) / (rho));
          const float u_0 =
              xi_1 * xi_23 * 0.5f +
              xi_1 * (vel0Term - xi_11 - xi_15 - xi_19 - xi_21 - xi_5);
          const float u_1 =
              xi_1 * xi_7 * 0.5f +
              xi_1 * (vel1Term - xi_10 - xi_11 - xi_17 - xi_18 + xi_4 - xi_6);
          const float u_2 =
              xi_1 * xi_20 * 0.5f + xi_1 * (vel2Term + xi_12 - xi_14 - xi_15 -
                                            xi_16 - xi_17 - xi_8 + xi_9);
          const float forceTerm_0 = omega_shear * u_0 * xi_23 * 0.5f +
                                    omega_shear * u_1 * xi_7 * 0.5f +
                                    omega_shear * u_2 * xi_20 * 0.5f -
                                    u_0 * xi_23 - u_1 * xi_7 - u_2 * xi_20;
          const float forceTerm_1 =
              omega_shear * u_0 * xi_23 * 0.083333333333333329f +
              omega_shear * u_1 * xi_7 * -0.16666666666666666f +
              omega_shear * u_2 * xi_20 * 0.083333333333333329f +
              rr_0 * xi_7 * -0.083333333333333329f +
              u_0 * xi_23 * -0.16666666666666666f +
              u_1 * xi_7 * 0.33333333333333331f +
              u_2 * xi_20 * -0.16666666666666666f + xi_7 * 0.16666666666666666f;
          const float forceTerm_2 =
              omega_shear * u_0 * xi_23 * 0.083333333333333329f +
              omega_shear * u_1 * xi_7 * -0.16666666666666666f +
              omega_shear * u_2 * xi_20 * 0.083333333333333329f +
              rr_0 * xi_7 * 0.083333333333333329f +
              u_0 * xi_23 * -0.16666666666666666f +
              u_1 * xi_7 * 0.33333333333333331f +
              u_2 * xi_20 * -0.16666666666666666f +
              xi_7 * -0.16666666666666666f;
          const float forceTerm_3 =
              omega_shear * u_0 * xi_23 * -0.16666666666666666f +
              omega_shear * u_1 * xi_7 * 0.083333333333333329f +
              omega_shear * u_2 * xi_20 * 0.083333333333333329f +
              rr_0 * xi_23 * 0.083333333333333329f +
              u_0 * xi_23 * 0.33333333333333331f +
              u_1 * xi_7 * -0.16666666666666666f +
              u_2 * xi_20 * -0.16666666666666666f +
              xi_23 * -0.16666666666666666f;
          const float forceTerm_4 =
              omega_shear * u_0 * xi_23 * -0.16666666666666666f +
              omega_shear * u_1 * xi_7 * 0.083333333333333329f +
              omega_shear * u_2 * xi_20 * 0.083333333333333329f +
              rr_0 * xi_23 * -0.083333333333333329f +
              u_0 * xi_23 * 0.33333333333333331f +
              u_1 * xi_7 * -0.16666666666666666f +
              u_2 * xi_20 * -0.16666666666666666f +
              xi_23 * 0.16666666666666666f;
          const float forceTerm_5 =
              omega_shear * u_0 * xi_23 * 0.083333333333333329f +
              omega_shear * u_1 * xi_7 * 0.083333333333333329f +
              omega_shear * u_2 * xi_20 * -0.16666666666666666f +
              rr_0 * xi_20 * -0.083333333333333329f +
              u_0 * xi_23 * -0.16666666666666666f +
              u_1 * xi_7 * -0.16666666666666666f +
              u_2 * xi_20 * 0.33333333333333331f + xi_20 * 0.16666666666666666f;
          const float forceTerm_6 =
              omega_shear * u_0 * xi_23 * 0.083333333333333329f +
              omega_shear * u_1 * xi_7 * 0.083333333333333329f +
              omega_shear * u_2 * xi_20 * -0.16666666666666666f +
              rr_0 * xi_20 * 0.083333333333333329f +
              u_0 * xi_23 * -0.16666666666666666f +
              u_1 * xi_7 * -0.16666666666666666f +
              u_2 * xi_20 * 0.33333333333333331f +
              xi_20 * -0.16666666666666666f;
          const float forceTerm_7 =
              omega_shear * u_0 * xi_23 * -0.083333333333333329f +
              omega_shear * u_0 * xi_7 * 0.125f +
              omega_shear * u_1 * xi_23 * 0.125f +
              omega_shear * u_1 * xi_7 * -0.083333333333333329f +
              omega_shear * u_2 * xi_20 * 0.041666666666666664f +
              rr_0 * xi_23 * 0.041666666666666664f +
              rr_0 * xi_7 * -0.041666666666666664f +
              u_0 * xi_23 * 0.16666666666666666f + u_0 * xi_7 * -0.25f +
              u_1 * xi_23 * -0.25f + u_1 * xi_7 * 0.16666666666666666f +
              u_2 * xi_20 * -0.083333333333333329f +
              xi_23 * -0.083333333333333329f + xi_7 * 0.083333333333333329f;
          const float forceTerm_8 =
              omega_shear * u_0 * xi_23 * -0.083333333333333329f +
              omega_shear * u_0 * xi_7 * -0.125f +
              omega_shear * u_1 * xi_23 * -0.125f +
              omega_shear * u_1 * xi_7 * -0.083333333333333329f +
              omega_shear * u_2 * xi_20 * 0.041666666666666664f +
              rr_0 * xi_23 * -0.041666666666666664f +
              rr_0 * xi_7 * -0.041666666666666664f +
              u_0 * xi_23 * 0.16666666666666666f + u_0 * xi_7 * 0.25f +
              u_1 * xi_23 * 0.25f + u_1 * xi_7 * 0.16666666666666666f +
              u_2 * xi_20 * -0.083333333333333329f +
              xi_23 * 0.083333333333333329f + xi_7 * 0.083333333333333329f;
          const float forceTerm_9 =
              omega_shear * u_0 * xi_23 * -0.083333333333333329f +
              omega_shear * u_0 * xi_7 * -0.125f +
              omega_shear * u_1 * xi_23 * -0.125f +
              omega_shear * u_1 * xi_7 * -0.083333333333333329f +
              omega_shear * u_2 * xi_20 * 0.041666666666666664f +
              rr_0 * xi_23 * 0.041666666666666664f +
              rr_0 * xi_7 * 0.041666666666666664f +
              u_0 * xi_23 * 0.16666666666666666f + u_0 * xi_7 * 0.25f +
              u_1 * xi_23 * 0.25f + u_1 * xi_7 * 0.16666666666666666f +
              u_2 * xi_20 * -0.083333333333333329f +
              xi_23 * -0.083333333333333329f + xi_7 * -0.083333333333333329f;
          const float forceTerm_10 =
              omega_shear * u_0 * xi_23 * -0.083333333333333329f +
              omega_shear * u_0 * xi_7 * 0.125f +
              omega_shear * u_1 * xi_23 * 0.125f +
              omega_shear * u_1 * xi_7 * -0.083333333333333329f +
              omega_shear * u_2 * xi_20 * 0.041666666666666664f +
              rr_0 * xi_23 * -0.041666666666666664f +
              rr_0 * xi_7 * 0.041666666666666664f +
              u_0 * xi_23 * 0.16666666666666666f + u_0 * xi_7 * -0.25f +
              u_1 * xi_23 * -0.25f + u_1 * xi_7 * 0.16666666666666666f +
              u_2 * xi_20 * -0.083333333333333329f +
              xi_23 * 0.083333333333333329f + xi_7 * -0.083333333333333329f;
          const float forceTerm_11 =
              omega_shear * u_0 * xi_23 * 0.041666666666666664f +
              omega_shear * u_1 * xi_20 * -0.125f +
              omega_shear * u_1 * xi_7 * -0.083333333333333329f +
              omega_shear * u_2 * xi_20 * -0.083333333333333329f +
              omega_shear * u_2 * xi_7 * -0.125f +
              rr_0 * xi_20 * -0.041666666666666664f +
              rr_0 * xi_7 * -0.041666666666666664f +
              u_0 * xi_23 * -0.083333333333333329f + u_1 * xi_20 * 0.25f +
              u_1 * xi_7 * 0.16666666666666666f +
              u_2 * xi_20 * 0.16666666666666666f + u_2 * xi_7 * 0.25f +
              xi_20 * 0.083333333333333329f + xi_7 * 0.083333333333333329f;
          const float forceTerm_12 =
              omega_shear * u_0 * xi_23 * 0.041666666666666664f +
              omega_shear * u_1 * xi_20 * 0.125f +
              omega_shear * u_1 * xi_7 * -0.083333333333333329f +
              omega_shear * u_2 * xi_20 * -0.083333333333333329f +
              omega_shear * u_2 * xi_7 * 0.125f +
              rr_0 * xi_20 * -0.041666666666666664f +
              rr_0 * xi_7 * 0.041666666666666664f +
              u_0 * xi_23 * -0.083333333333333329f + u_1 * xi_20 * -0.25f +
              u_1 * xi_7 * 0.16666666666666666f +
              u_2 * xi_20 * 0.16666666666666666f + u_2 * xi_7 * -0.25f +
              xi_20 * 0.083333333333333329f + xi_7 * -0.083333333333333329f;
          const float forceTerm_13 =
              omega_shear * u_0 * xi_20 * 0.125f +
              omega_shear * u_0 * xi_23 * -0.083333333333333329f +
              omega_shear * u_1 * xi_7 * 0.041666666666666664f +
              omega_shear * u_2 * xi_20 * -0.083333333333333329f +
              omega_shear * u_2 * xi_23 * 0.125f +
              rr_0 * xi_20 * -0.041666666666666664f +
              rr_0 * xi_23 * 0.041666666666666664f + u_0 * xi_20 * -0.25f +
              u_0 * xi_23 * 0.16666666666666666f +
              u_1 * xi_7 * -0.083333333333333329f +
              u_2 * xi_20 * 0.16666666666666666f + u_2 * xi_23 * -0.25f +
              xi_20 * 0.083333333333333329f + xi_23 * -0.083333333333333329f;
          const float forceTerm_14 =
              omega_shear * u_0 * xi_20 * -0.125f +
              omega_shear * u_0 * xi_23 * -0.083333333333333329f +
              omega_shear * u_1 * xi_7 * 0.041666666666666664f +
              omega_shear * u_2 * xi_20 * -0.083333333333333329f +
              omega_shear * u_2 * xi_23 * -0.125f +
              rr_0 * xi_20 * -0.041666666666666664f +
              rr_0 * xi_23 * -0.041666666666666664f + u_0 * xi_20 * 0.25f +
              u_0 * xi_23 * 0.16666666666666666f +
              u_1 * xi_7 * -0.083333333333333329f +
              u_2 * xi_20 * 0.16666666666666666f + u_2 * xi_23 * 0.25f +
              xi_20 * 0.083333333333333329f + xi_23 * 0.083333333333333329f;
          const float forceTerm_15 =
              omega_shear * u_0 * xi_23 * 0.041666666666666664f +
              omega_shear * u_1 * xi_20 * 0.125f +
              omega_shear * u_1 * xi_7 * -0.083333333333333329f +
              omega_shear * u_2 * xi_20 * -0.083333333333333329f +
              omega_shear * u_2 * xi_7 * 0.125f +
              rr_0 * xi_20 * 0.041666666666666664f +
              rr_0 * xi_7 * -0.041666666666666664f +
              u_0 * xi_23 * -0.083333333333333329f + u_1 * xi_20 * -0.25f +
              u_1 * xi_7 * 0.16666666666666666f +
              u_2 * xi_20 * 0.16666666666666666f + u_2 * xi_7 * -0.25f +
              xi_20 * -0.083333333333333329f + xi_7 * 0.083333333333333329f;
          const float forceTerm_16 =
              omega_shear * u_0 * xi_23 * 0.041666666666666664f +
              omega_shear * u_1 * xi_20 * -0.125f +
              omega_shear * u_1 * xi_7 * -0.083333333333333329f +
              omega_shear * u_2 * xi_20 * -0.083333333333333329f +
              omega_shear * u_2 * xi_7 * -0.125f +
              rr_0 * xi_20 * 0.041666666666666664f +
              rr_0 * xi_7 * 0.041666666666666664f +
              u_0 * xi_23 * -0.083333333333333329f + u_1 * xi_20 * 0.25f +
              u_1 * xi_7 * 0.16666666666666666f +
              u_2 * xi_20 * 0.16666666666666666f + u_2 * xi_7 * 0.25f +
              xi_20 * -0.083333333333333329f + xi_7 * -0.083333333333333329f;
          const float forceTerm_17 =
              omega_shear * u_0 * xi_20 * -0.125f +
              omega_shear * u_0 * xi_23 * -0.083333333333333329f +
              omega_shear * u_1 * xi_7 * 0.041666666666666664f +
              omega_shear * u_2 * xi_20 * -0.083333333333333329f +
              omega_shear * u_2 * xi_23 * -0.125f +
              rr_0 * xi_20 * 0.041666666666666664f +
              rr_0 * xi_23 * 0.041666666666666664f + u_0 * xi_20 * 0.25f +
              u_0 * xi_23 * 0.16666666666666666f +
              u_1 * xi_7 * -0.083333333333333329f +
              u_2 * xi_20 * 0.16666666666666666f + u_2 * xi_23 * 0.25f +
              xi_20 * -0.083333333333333329f + xi_23 * -0.083333333333333329f;
          const float forceTerm_18 =
              omega_shear * u_0 * xi_20 * 0.125f +
              omega_shear * u_0 * xi_23 * -0.083333333333333329f +
              omega_shear * u_1 * xi_7 * 0.041666666666666664f +
              omega_shear * u_2 * xi_20 * -0.083333333333333329f +
              omega_shear * u_2 * xi_23 * 0.125f +
              rr_0 * xi_20 * 0.041666666666666664f +
              rr_0 * xi_23 * -0.041666666666666664f + u_0 * xi_20 * -0.25f +
              u_0 * xi_23 * 0.16666666666666666f +
              u_1 * xi_7 * -0.083333333333333329f +
              u_2 * xi_20 * 0.16666666666666666f + u_2 * xi_23 * -0.25f +
              xi_20 * -0.083333333333333329f + xi_23 * 0.083333333333333329f;
          const float u0Mu1 = u_0 - u_1;
          const float u0Pu1 = u_0 + u_1;
          const float u1Pu2 = u_1 + u_2;
          const float u1Mu2 = u_1 - u_2;
          const float u0Mu2 = u_0 - u_2;
          const float u0Pu2 = u_0 + u_2;
          const float f_eq_common =
              rho - rho * (u_0 * u_0) - rho * (u_1 * u_1) - rho * (u_2 * u_2);
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0] =
              forceTerm_0 +
              omega_shear * (f_eq_common * 0.33333333333333331f - xi_13) +
              xi_13;
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     _stride_pdfs_3 + ctr_0] =
              forceTerm_1 +
              omega_shear * (f_eq_common * 0.16666666666666666f +
                             rho * (-0.1111111111111111f +
                                    0.33333333333333331f * (u_1 * u_1)) +
                             xi_18 * -0.5f + xi_3 * -0.5f) +
              rr_0 * (rho * u_1 * 0.16666666666666666f + xi_18 * 0.5f +
                      xi_3 * -0.5f) +
              xi_3 +
              ((-1.0f <= -grid_size + ((float)(ctr_1)))
                   ? (rho * v_s * (u_0 * 2.0f + v_s) * 0.16666666666666666f)
                   : (0.0f));
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     2 * _stride_pdfs_3 + ctr_0] =
              forceTerm_2 +
              omega_shear * (f_eq_common * 0.16666666666666666f +
                             rho * (-0.1111111111111111f +
                                    0.33333333333333331f * (u_1 * u_1)) +
                             xi_18 * -0.5f + xi_3 * -0.5f) +
              rr_0 * (rho * u_1 * -0.16666666666666666f + xi_18 * -0.5f +
                      xi_3 * 0.5f) +
              xi_18 +
              ((0.0f >= ((float)(ctr_1)))
                   ? (rho * v_s * (u_0 * -2.0f + v_s) * 0.16666666666666666f)
                   : (0.0f));
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     3 * _stride_pdfs_3 + ctr_0] =
              forceTerm_3 +
              omega_shear * (f_eq_common * 0.16666666666666666f +
                             rho * (-0.1111111111111111f +
                                    0.33333333333333331f * (u_0 * u_0)) +
                             xi_24 * -0.5f + xi_5 * -0.5f) +
              rr_0 * (rho * u_0 * -0.16666666666666666f + xi_24 * 0.5f +
                      xi_5 * -0.5f) +
              xi_5;
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     4 * _stride_pdfs_3 + ctr_0] =
              forceTerm_4 +
              omega_shear * (f_eq_common * 0.16666666666666666f +
                             rho * (-0.1111111111111111f +
                                    0.33333333333333331f * (u_0 * u_0)) +
                             xi_24 * -0.5f + xi_5 * -0.5f) +
              rr_0 * (rho * u_0 * 0.16666666666666666f + xi_24 * -0.5f +
                      xi_5 * 0.5f) +
              xi_24;
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     5 * _stride_pdfs_3 + ctr_0] =
              forceTerm_5 +
              omega_shear * (f_eq_common * 0.16666666666666666f +
                             rho * (-0.1111111111111111f +
                                    0.33333333333333331f * (u_2 * u_2)) +
                             xi_16 * -0.5f + xi_22 * -0.5f) +
              rr_0 * (rho * u_2 * 0.16666666666666666f + xi_16 * 0.5f +
                      xi_22 * -0.5f) +
              xi_22;
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     6 * _stride_pdfs_3 + ctr_0] =
              forceTerm_6 +
              omega_shear * (f_eq_common * 0.16666666666666666f +
                             rho * (-0.1111111111111111f +
                                    0.33333333333333331f * (u_2 * u_2)) +
                             xi_16 * -0.5f + xi_22 * -0.5f) +
              rr_0 * (rho * u_2 * -0.16666666666666666f + xi_16 * -0.5f +
                      xi_22 * 0.5f) +
              xi_16;
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     7 * _stride_pdfs_3 + ctr_0] =
              forceTerm_7 +
              omega_shear * (f_eq_common * 0.041666666666666664f +
                             rho * (-0.013888888888888888f +
                                    0.041666666666666664f * (u_2 * u_2) +
                                    0.125f * (u0Mu1 * u0Mu1)) +
                             xi_21 * -0.5f + xi_6 * -0.5f) +
              rr_0 * (rho * u0Mu1 * -0.083333333333333329f + xi_21 * -0.5f +
                      xi_6 * 0.5f) +
              xi_21 +
              ((-1.0f <= -grid_size + ((float)(ctr_1)))
                   ? (rho * v_s * (u_0 * -2.0f + u_1 * 3.0f - v_s + 1.0f) *
                      0.083333333333333329f)
                   : (0.0f));
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     8 * _stride_pdfs_3 + ctr_0] =
              forceTerm_8 +
              omega_shear * (f_eq_common * 0.041666666666666664f +
                             rho * (-0.013888888888888888f +
                                    0.041666666666666664f * (u_2 * u_2) +
                                    0.125f * (u0Pu1 * u0Pu1)) +
                             xi_11 * -0.5f + xi_4 * -0.5f) +
              rr_0 * (rho * u0Pu1 * 0.083333333333333329f + xi_11 * 0.5f +
                      xi_4 * -0.5f) +
              xi_4 +
              ((-1.0f <= -grid_size + ((float)(ctr_1)))
                   ? (rho * v_s * (u_0 * 2.0f + u_1 * 3.0f + v_s + 1.0f) *
                      -0.083333333333333329f)
                   : (0.0f));
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     9 * _stride_pdfs_3 + ctr_0] =
              forceTerm_9 +
              omega_shear * (f_eq_common * 0.041666666666666664f +
                             rho * (-0.013888888888888888f +
                                    0.041666666666666664f * (u_2 * u_2) +
                                    0.125f * (u0Pu1 * u0Pu1)) +
                             xi_11 * -0.5f + xi_4 * -0.5f) +
              rr_0 * (rho * u0Pu1 * -0.083333333333333329f + xi_11 * -0.5f +
                      xi_4 * 0.5f) +
              xi_11 +
              ((0.0f >= ((float)(ctr_1)))
                   ? (rho * v_s * (u_0 * 2.0f + u_1 * 3.0f - v_s - 1.0f) *
                      0.083333333333333329f)
                   : (0.0f));
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     10 * _stride_pdfs_3 + ctr_0] =
              forceTerm_10 +
              omega_shear * (f_eq_common * 0.041666666666666664f +
                             rho * (-0.013888888888888888f +
                                    0.041666666666666664f * (u_2 * u_2) +
                                    0.125f * (u0Mu1 * u0Mu1)) +
                             xi_21 * -0.5f + xi_6 * -0.5f) +
              rr_0 * (rho * u0Mu1 * 0.083333333333333329f + xi_21 * 0.5f +
                      xi_6 * -0.5f) +
              xi_6 +
              ((0.0f >= ((float)(ctr_1)))
                   ? (rho * v_s * (u_0 * 2.0f + u_1 * -3.0f - v_s + 1.0f) *
                      0.083333333333333329f)
                   : (0.0f));
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     11 * _stride_pdfs_3 + ctr_0] =
              forceTerm_11 +
              omega_shear * (f_eq_common * 0.041666666666666664f +
                             rho * (-0.013888888888888888f +
                                    0.041666666666666664f * (u_0 * u_0) +
                                    0.125f * (u1Pu2 * u1Pu2)) +
                             xi_12 * -0.5f + xi_17 * -0.5f) +
              rr_0 * (rho * u1Pu2 * 0.083333333333333329f + xi_12 * -0.5f +
                      xi_17 * 0.5f) +
              xi_12;
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     12 * _stride_pdfs_3 + ctr_0] =
              forceTerm_12 +
              omega_shear * (f_eq_common * 0.041666666666666664f +
                             rho * (-0.013888888888888888f +
                                    0.041666666666666664f * (u_0 * u_0) +
                                    0.125f * (u1Mu2 * u1Mu2)) +
                             xi_10 * -0.5f + xi_14 * -0.5f) +
              rr_0 * (rho * u1Mu2 * -0.083333333333333329f + xi_10 * -0.5f +
                      xi_14 * 0.5f) +
              xi_10;
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     13 * _stride_pdfs_3 + ctr_0] =
              forceTerm_13 +
              omega_shear * (f_eq_common * 0.041666666666666664f +
                             rho * (-0.013888888888888888f +
                                    0.041666666666666664f * (u_1 * u_1) +
                                    0.125f * (u0Mu2 * u0Mu2)) +
                             xi_19 * -0.5f + xi_8 * -0.5f) +
              rr_0 * (rho * u0Mu2 * -0.083333333333333329f + xi_19 * -0.5f +
                      xi_8 * 0.5f) +
              xi_19;
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     14 * _stride_pdfs_3 + ctr_0] =
              forceTerm_14 +
              omega_shear * (f_eq_common * 0.041666666666666664f +
                             rho * (-0.013888888888888888f +
                                    0.041666666666666664f * (u_1 * u_1) +
                                    0.125f * (u0Pu2 * u0Pu2)) +
                             xi_15 * -0.5f + xi_9 * -0.5f) +
              rr_0 * (rho * u0Pu2 * 0.083333333333333329f + xi_15 * 0.5f +
                      xi_9 * -0.5f) +
              xi_9;
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     15 * _stride_pdfs_3 + ctr_0] =
              forceTerm_15 +
              omega_shear * (f_eq_common * 0.041666666666666664f +
                             rho * (-0.013888888888888888f +
                                    0.041666666666666664f * (u_0 * u_0) +
                                    0.125f * (u1Mu2 * u1Mu2)) +
                             xi_10 * -0.5f + xi_14 * -0.5f) +
              rr_0 * (rho * u1Mu2 * 0.083333333333333329f + xi_10 * 0.5f +
                      xi_14 * -0.5f) +
              xi_14;
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     16 * _stride_pdfs_3 + ctr_0] =
              forceTerm_16 +
              omega_shear * (f_eq_common * 0.041666666666666664f +
                             rho * (-0.013888888888888888f +
                                    0.041666666666666664f * (u_0 * u_0) +
                                    0.125f * (u1Pu2 * u1Pu2)) +
                             xi_12 * -0.5f + xi_17 * -0.5f) +
              rr_0 * (rho * u1Pu2 * -0.083333333333333329f + xi_12 * 0.5f +
                      xi_17 * -0.5f) +
              xi_17;
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     17 * _stride_pdfs_3 + ctr_0] =
              forceTerm_17 +
              omega_shear * (f_eq_common * 0.041666666666666664f +
                             rho * (-0.013888888888888888f +
                                    0.041666666666666664f * (u_1 * u_1) +
                                    0.125f * (u0Pu2 * u0Pu2)) +
                             xi_15 * -0.5f + xi_9 * -0.5f) +
              rr_0 * (rho * u0Pu2 * -0.083333333333333329f + xi_15 * -0.5f +
                      xi_9 * 0.5f) +
              xi_15;
          _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 +
                     18 * _stride_pdfs_3 + ctr_0] =
              forceTerm_18 +
              omega_shear * (f_eq_common * 0.041666666666666664f +
                             rho * (-0.013888888888888888f +
                                    0.041666666666666664f * (u_1 * u_1) +
                                    0.125f * (u0Mu2 * u0Mu2)) +
                             xi_19 * -0.5f + xi_8 * -0.5f) +
              rr_0 * (rho * u0Mu2 * 0.083333333333333329f + xi_19 * 0.5f +
                      xi_8 * -0.5f) +
              xi_8;
        }
      }
    }
  }
}
} // namespace internal_9a18f2f4073cdcc5365cdfddb752069e

void CollideSweepSinglePrecisionLeesEdwardsAVX::run(IBlock *block) {

  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto &omega_shear = this->omega_shear_;
  auto &grid_size = this->grid_size_;
  auto &v_s = this->v_s_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()))
  float *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()))
  float *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(int64_c(force->xSize()) + 0))
  const int64_t _size_force_0 = int64_t(int64_c(force->xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(int64_c(force->ySize()) + 0))
  const int64_t _size_force_1 = int64_t(int64_c(force->ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(int64_c(force->zSize()) + 0))
  const int64_t _size_force_2 = int64_t(int64_c(force->zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_9a18f2f4073cdcc5365cdfddb752069e::
      collidesweepsingleprecisionleesedwardsavx_collidesweepsingleprecisionleesedwardsavx(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1,
          _stride_pdfs_2, _stride_pdfs_3, grid_size, omega_shear, v_s);
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

  auto &omega_shear = this->omega_shear_;
  auto &grid_size = this->grid_size_;
  auto &v_s = this->v_s_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()))
  float *RESTRICT const _data_force =
      force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
  float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(int64_c(ci.xSize()) + 0))
  const int64_t _size_force_0 = int64_t(int64_c(ci.xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(int64_c(ci.ySize()) + 0))
  const int64_t _size_force_1 = int64_t(int64_c(ci.ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(int64_c(ci.zSize()) + 0))
  const int64_t _size_force_2 = int64_t(int64_c(ci.zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_9a18f2f4073cdcc5365cdfddb752069e::
      collidesweepsingleprecisionleesedwardsavx_collidesweepsingleprecisionleesedwardsavx(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1,
          _stride_pdfs_2, _stride_pdfs_3, grid_size, omega_shear, v_s);
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
