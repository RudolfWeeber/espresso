// kernel generated with pystencils v1.0+0.g354fede.dirty, lbmpy v1.0,
// lbmpy_walberla/pystencils_walberla from commit
// e1fe2ad1dcbe8f31ea79d95e8a5a5cc0ee3691f3

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
//! \\file CollideSweepSinglePrecision.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepSinglePrecision.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

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

namespace internal_ac3754ab06497f05a231a858fd02d610 {
static FUNC_PREFIX void collidesweepsingleprecision_collidesweepsingleprecision(
    float *RESTRICT const _data_force, float *RESTRICT _data_pdfs,
    int64_t const _size_force_0, int64_t const _size_force_1,
    int64_t const _size_force_2, int64_t const _stride_force_0,
    int64_t const _stride_force_1, int64_t const _stride_force_2,
    int64_t const _stride_force_3, int64_t const _stride_pdfs_0,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, float omega_bulk, float omega_even,
    float omega_odd, float omega_shear) {
  const float xi_22 = omega_shear * -0.5f + 1.0f;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    float *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      float *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      float *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      float *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      float *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      float *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      float *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const float xi_59 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const float xi_60 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const float xi_61 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const float xi_62 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const float xi_63 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const float xi_64 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const float xi_65 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const float xi_66 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const float xi_67 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const float xi_68 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const float xi_69 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const float xi_70 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const float xi_71 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const float xi_72 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const float xi_73 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const float xi_74 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const float xi_75 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const float xi_76 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const float xi_77 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const float xi_78 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const float xi_79 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const float xi_80 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const float xi_0 = xi_60 + xi_72;
        const float xi_1 = xi_66 + xi_69;
        const float xi_2 = xi_65 + xi_76;
        const float xi_3 = xi_68 + xi_73;
        const float xi_4 = xi_67 + xi_75;
        const float xi_6 = xi_59 + xi_61;
        const float xi_7 = xi_79 + xi_80;
        const float partial_m_m1_0_e_0 = xi_0 + xi_70;
        const float partial_m_m1_e_0_0 = partial_m_m1_0_e_0 + xi_1;
        const float partial_m_0_m1_e_0 = xi_2 + xi_74;
        const float partial_m_0_0_e_0 = xi_3 + xi_78;
        const float partial_m_0_1_e_0 = xi_4 + xi_62;
        const float xi_5 = partial_m_0_1_e_0 + partial_m_0_m1_e_0;
        const float partial_m_0_e_0_0 = partial_m_0_0_e_0 + xi_5;
        const float partial_m_1_0_e_0 = xi_6 + xi_71;
        const float partial_m_1_e_0_0 = partial_m_1_0_e_0 + xi_7;
        const float xi_10 = partial_m_1_e_0_0 + partial_m_m1_e_0_0;
        const float partial_m_m1_e_1_0 = xi_66 + xi_69 * -1.0f;
        const float partial_m_0_e_1_0 =
            partial_m_0_1_e_0 + partial_m_0_m1_e_0 * -1.0f;
        const float partial_m_1_e_1_0 = xi_79 * -1.0f + xi_80;
        const float xi_11 = partial_m_1_e_1_0 + partial_m_m1_e_1_0;
        const float partial_m_m1_0_e_1 = xi_60 * -1.0f + xi_72;
        const float partial_m_0_m1_e_1 = xi_65 + xi_76 * -1.0f;
        const float partial_m_0_0_e_1 = xi_68 + xi_73 * -1.0f;
        const float partial_m_0_1_e_1 = xi_67 + xi_75 * -1.0f;
        const float xi_8 = partial_m_0_1_e_1 + partial_m_0_m1_e_1;
        const float partial_m_0_e_0_1 = partial_m_0_0_e_1 + xi_8;
        const float partial_m_1_0_e_1 = xi_59 + xi_61 * -1.0f;
        const float xi_12 = partial_m_1_0_e_1 + partial_m_m1_0_e_1;
        const float partial_m_m1_e_2_0 = xi_1;
        const float partial_m_0_e_2_0 = xi_5;
        const float partial_m_1_e_2_0 = xi_7;
        const float xi_13 = partial_m_1_e_2_0 + partial_m_m1_e_2_0;
        const float partial_m_m1_0_e_2 = xi_0;
        const float partial_m_0_m1_e_2 = xi_2;
        const float partial_m_0_0_e_2 = xi_3;
        const float partial_m_0_1_e_2 = xi_4;
        const float xi_9 = partial_m_0_1_e_2 + partial_m_0_m1_e_2;
        const float partial_m_0_e_0_2 = partial_m_0_0_e_2 + xi_9;
        const float partial_m_1_0_e_2 = xi_6;
        const float xi_14 = partial_m_1_0_e_2 + partial_m_m1_0_e_2;
        const float partial_m_0_e_1_1 =
            partial_m_0_1_e_1 + partial_m_0_m1_e_1 * -1.0f;
        const float partial_m_0_e_2_1 = xi_8;
        const float partial_m_0_e_1_2 =
            partial_m_0_1_e_2 + partial_m_0_m1_e_2 * -1.0f;
        const float partial_m_0_e_2_2 = xi_9;
        const float m_000 = partial_m_0_e_0_0 + xi_10;
        const float m_100 = partial_m_1_e_0_0 + partial_m_m1_e_0_0 * -1.0f;
        const float xi_19 = m_100 * -1.0f;
        const float m_010 = partial_m_0_e_1_0 + xi_11;
        const float xi_17 = m_010 * -1.0f;
        const float m_001 = partial_m_0_e_0_1 + xi_12;
        const float xi_18 = m_001 * -1.0f;
        const float m_200 = xi_10;
        const float xi_20 = m_000 + m_200 * -6.0f;
        const float m_020 = partial_m_0_e_2_0 + xi_13;
        const float m_002 = partial_m_0_e_0_2 + xi_14;
        const float xi_16 = m_002 * -1.0f;
        const float m_110 = partial_m_1_e_1_0 + partial_m_m1_e_1_0 * -1.0f;
        const float m_101 = partial_m_1_0_e_1 + partial_m_m1_0_e_1 * -1.0f;
        const float m_210 = xi_11;
        const float m_201 = xi_12;
        const float m_120 = partial_m_1_e_2_0 + partial_m_m1_e_2_0 * -1.0f;
        const float m_102 = partial_m_1_0_e_2 + partial_m_m1_0_e_2 * -1.0f;
        const float m_220 = xi_13;
        const float xi_21 = m_002 * -4.0f + m_220 * 4.0f;
        const float m_202 = xi_14;
        const float sub_f_to_m_0 = ((1.0f) / (m_000));
        const float xi_15 = sub_f_to_m_0 * 0.5f;
        const float u_0 = m_100 * sub_f_to_m_0 + xi_15 * xi_63;
        const float xi_23 = u_0 * xi_63;
        const float xi_27 = m_000 * (u_0 * u_0);
        const float xi_31 = m_000 * u_0;
        const float u_1 = m_010 * sub_f_to_m_0 + xi_15 * xi_77;
        const float xi_24 = u_1 * xi_77 * 2.0f;
        const float xi_28 = m_000 * (u_1 * u_1);
        const float u_2 = m_001 * sub_f_to_m_0 + xi_15 * xi_64;
        const float xi_25 = u_2 * xi_64 * 2.0f;
        const float xi_26 = xi_25 * -1.0f;
        const float xi_29 = m_000 * (u_2 * u_2);
        const float xi_30 = xi_29 * -1.0f;
        const float xi_32 = xi_29 * 0.66666666666666663f;
        const float M_4 = m_020 * -1.0f + m_200 * 2.0f + xi_16;
        const float M_5 = m_020 + xi_16;
        const float M_9 = m_000 * -1.0f + m_002 + m_020 + m_200;
        const float M_10 = m_210 * 3.0f + xi_17;
        const float M_11 = m_201 * 3.0f + xi_18;
        const float M_12 = m_120 * 3.0f + xi_19;
        const float M_13 = m_102 * 2.0f + m_120 + xi_19;
        const float M_14 = m_201 + partial_m_0_e_2_1 * 2.0f + xi_18;
        const float M_15 = m_210 + partial_m_0_e_1_2 * 2.0f + xi_17;
        const float M_16 = m_002 * 3.0f + m_020 * -6.0f + m_220 * 18.0f + xi_20;
        const float M_17 = m_020 + m_202 * 14.0f + xi_20 + xi_21;
        const float M_18 = m_000 + m_020 * -4.0f + m_200 * -1.0f +
                           m_202 * 4.0f + partial_m_0_e_2_2 * 10.0f + xi_21;
        const float M_post_1 = m_100 + xi_63;
        const float xi_39 = M_post_1 * 0.33333333333333331f;
        const float M_post_2 = m_010 + xi_77;
        const float xi_37 = M_post_2 * 0.33333333333333331f;
        const float M_post_3 = m_001 + xi_64;
        const float xi_38 = M_post_3 * 0.33333333333333331f;
        const float M_post_4 =
            M_4 +
            omega_shear * (M_4 * -1.0f + xi_27 * 2.0f + xi_28 * -1.0f + xi_30) +
            xi_22 * (xi_23 * 4.0f + xi_24 * -1.0f + xi_26);
        const float xi_35 = M_post_4 * -0.16666666666666666f;
        const float M_post_5 = M_5 +
                               omega_shear * (M_5 * -1.0f + xi_28 + xi_30) +
                               xi_22 * (xi_24 + xi_26);
        const float xi_34 = M_post_5 * 0.5f;
        const float xi_40 = M_post_5 * 0.25f;
        const float M_post_6 = m_110 +
                               omega_shear * (m_110 * -1.0f + u_1 * xi_31) +
                               xi_22 * (u_0 * xi_77 + u_1 * xi_63);
        const float xi_47 = M_post_6 * 0.25f;
        const float M_post_7 = m_101 +
                               omega_shear * (m_101 * -1.0f + u_2 * xi_31) +
                               xi_22 * (u_0 * xi_64 + u_2 * xi_63);
        const float xi_55 = M_post_7 * 0.25f;
        const float M_post_8 =
            omega_shear * (m_000 * u_1 * u_2 + partial_m_0_e_1_1 * -1.0f) +
            partial_m_0_e_1_1 + xi_22 * (u_1 * xi_64 + u_2 * xi_77);
        const float xi_51 = M_post_8 * 0.25f;
        const float M_post_9 =
            M_9 + omega_bulk * (M_9 * -1.0f + xi_27 + xi_28 + xi_29) +
            (omega_bulk * -0.5f + 1.0f) * (xi_23 * 2.0f + xi_24 + xi_25);
        const float xi_33 =
            M_post_9 * 0.33333333333333331f + m_000 * 0.33333333333333331f;
        const float xi_36 = xi_33 + xi_35;
        const float xi_41 =
            M_post_9 * 0.16666666666666666f + m_000 * 0.1111111111111111f;
        const float xi_42 = M_post_4 * 0.083333333333333329f + xi_41;
        const float M_post_10 = M_10 * omega_odd * -1.0f + M_10;
        const float M_post_11 = M_11 * omega_odd * -1.0f + M_11;
        const float M_post_12 = M_12 * omega_odd * -1.0f + M_12;
        const float M_post_13 = M_13 * omega_odd * -1.0f + M_13;
        const float M_post_14 = M_14 * omega_odd * -1.0f + M_14;
        const float M_post_15 = M_15 * omega_odd * -1.0f + M_15;
        const float M_post_16 =
            M_16 + omega_even * (M_16 * -1.0f + xi_29 * 3.0f);
        const float xi_43 = M_post_16 * -0.015873015873015872f;
        const float M_post_17 =
            M_17 +
            omega_even * (M_17 * -1.0f + xi_28 * 2.3333333333333335f + xi_32);
        const float M_post_18 =
            M_18 + omega_even * (M_18 * -1.0f + xi_27 * 1.6666666666666667f +
                                 xi_28 * 0.66666666666666663f + xi_32);
        const float m_post_200 = M_post_4 * 0.33333333333333331f + xi_33;
        const float m_post_020 = xi_34 + xi_36;
        const float m_post_002 = xi_34 * -1.0f + xi_36;
        const float m_post_210 = M_post_10 * 0.33333333333333331f + xi_37;
        const float xi_50 = m_post_210 * 0.25f;
        const float m_post_201 = M_post_11 * 0.33333333333333331f + xi_38;
        const float xi_58 = m_post_201 * 0.25f;
        const float m_post_120 = M_post_12 * 0.33333333333333331f + xi_39;
        const float xi_49 = m_post_120 * 0.25f;
        const float m_post_102 =
            M_post_12 * -0.16666666666666666f + M_post_13 * 0.5f + xi_39;
        const float xi_57 = m_post_102 * 0.25f;
        const float m_post_021 =
            M_post_11 * -0.16666666666666666f + M_post_14 * 0.5f + xi_38;
        const float xi_54 = m_post_021 * 0.25f;
        const float m_post_012 =
            M_post_10 * -0.16666666666666666f + M_post_15 * 0.5f + xi_37;
        const float xi_53 = m_post_012 * 0.25f;
        const float m_post_220 =
            M_post_16 * 0.055555555555555552f + xi_40 + xi_42;
        const float xi_45 = m_post_220 * -0.5f;
        const float xi_48 = m_post_220 * 0.25f;
        const float m_post_202 =
            M_post_17 * 0.071428571428571425f + xi_40 * -1.0f + xi_42 + xi_43;
        const float xi_46 = m_post_202 * -0.5f;
        const float xi_56 = m_post_202 * 0.25f;
        const float m_post_022 = M_post_17 * -0.028571428571428571f +
                                 M_post_18 * 0.10000000000000001f + xi_35 +
                                 xi_41 + xi_43;
        const float xi_44 = m_post_022 * -0.5f;
        const float xi_52 = m_post_022 * 0.25f;
        const float sub_k_to_f_20 = m_post_020 * 0.5f + xi_44 + xi_45;
        const float sub_k_to_f_21 =
            M_post_2 * 0.5f + m_post_012 * -0.5f + m_post_210 * -0.5f;
        const float sub_k_to_f_22 = m_post_200 * 0.5f + xi_45 + xi_46;
        const float sub_k_to_f_23 =
            M_post_1 * -0.5f + m_post_102 * 0.5f + m_post_120 * 0.5f;
        const float sub_k_to_f_24 = m_post_002 * 0.5f + xi_44 + xi_46;
        const float sub_k_to_f_25 =
            M_post_3 * 0.5f + m_post_021 * -0.5f + m_post_201 * -0.5f;
        const float sub_k_to_f_26 = xi_47 * -1.0f + xi_48;
        const float sub_k_to_f_27 = xi_49 * -1.0f + xi_50;
        const float sub_k_to_f_28 = xi_47 + xi_48;
        const float sub_k_to_f_29 = xi_49 + xi_50;
        const float sub_k_to_f_30 = xi_51 + xi_52;
        const float sub_k_to_f_31 = xi_53 + xi_54;
        const float sub_k_to_f_32 = xi_51 * -1.0f + xi_52;
        const float sub_k_to_f_33 = xi_53 * -1.0f + xi_54;
        const float sub_k_to_f_34 = xi_55 * -1.0f + xi_56;
        const float sub_k_to_f_35 = xi_57 * -1.0f + xi_58;
        const float sub_k_to_f_36 = xi_55 + xi_56;
        const float sub_k_to_f_37 = xi_57 + xi_58;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] =
            m_000 + m_post_002 * -1.0f + m_post_020 * -1.0f + m_post_022 +
            m_post_200 * -1.0f + m_post_202 + m_post_220;
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_20 + sub_k_to_f_21;
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_20 + sub_k_to_f_21 * -1.0f;
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_22 + sub_k_to_f_23;
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_22 + sub_k_to_f_23 * -1.0f;
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_24 + sub_k_to_f_25;
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_24 + sub_k_to_f_25 * -1.0f;
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_26 + sub_k_to_f_27;
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_28 + sub_k_to_f_29;
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_28 + sub_k_to_f_29 * -1.0f;
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_26 + sub_k_to_f_27 * -1.0f;
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_30 + sub_k_to_f_31;
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_32 + sub_k_to_f_33;
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_34 + sub_k_to_f_35;
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_36 + sub_k_to_f_37;
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_32 + sub_k_to_f_33 * -1.0f;
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_30 + sub_k_to_f_31 * -1.0f;
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_36 + sub_k_to_f_37 * -1.0f;
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_34 + sub_k_to_f_35 * -1.0f;
      }
    }
  }
}
} // namespace internal_ac3754ab06497f05a231a858fd02d610

void CollideSweepSinglePrecision::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto &omega_bulk = this->omega_bulk_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_odd = this->omega_odd_;
  auto &omega_even = this->omega_even_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
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
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_ac3754ab06497f05a231a858fd02d610::
      collidesweepsingleprecision_collidesweepsingleprecision(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
          omega_bulk, omega_even, omega_odd, omega_shear);
}

void CollideSweepSinglePrecision::runOnCellInterval(
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

  auto &omega_bulk = this->omega_bulk_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_odd = this->omega_odd_;
  auto &omega_even = this->omega_even_;
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
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_ac3754ab06497f05a231a858fd02d610::
      collidesweepsingleprecision_collidesweepsingleprecision(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
          omega_bulk, omega_even, omega_odd, omega_shear);
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