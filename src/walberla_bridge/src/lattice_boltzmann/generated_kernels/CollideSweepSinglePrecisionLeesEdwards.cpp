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
//! \\file CollideSweepSinglePrecisionLeesEdwards.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7, lbmpy v1.3.7, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 59c9b8b185782eba184e0fdfb2144793343213f0


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "CollideSweepSinglePrecisionLeesEdwards.h"




#define FUNC_PREFIX

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning push
#pragma warning( disable :  1599 )
#endif

using namespace std;

namespace walberla {
namespace pystencils {


namespace internal_ab1f3bc3368574afb482da84ccb58898 {
static FUNC_PREFIX void collidesweepsingleprecisionleesedwards_collidesweepsingleprecisionleesedwards(float * RESTRICT const _data_force, float * RESTRICT  _data_pdfs, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_0, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, float grid_size, float omega_bulk, float omega_even, float omega_odd, float omega_shear, float rho, float v_s)
{
#ifdef _OPENMP
   #pragma omp parallel
#endif
   {
      const float xi_2 = 0.50000000000000000f;
      const int64_t xi_3 = 2;
      const int64_t xi_4 = 3;
      const int64_t xi_5 = 6;
      const int64_t xi_6 = 18;
      const int64_t xi_7 = 4;
      const int64_t xi_8 = 14;
      const int64_t xi_9 = 10;
      const float xi_10 = 0.66666666666666663f;
      const float xi_11 = 2.3333333333333335f;
      const float xi_12 = 1.6666666666666667f;
      const float xi_13 = 0.33333333333333331f;
      const float xi_14 = 0.16666666666666666f;
      const float xi_15 = 0.25000000000000000f;
      const float xi_16 = 0.1111111111111111f;
      const float xi_17 = 0.083333333333333329f;
      const float xi_18 = 0.055555555555555552f;
      const float xi_19 = 0.015873015873015872f;
      const float xi_20 = 0.071428571428571425f;
      const float xi_21 = 0.028571428571428571f;
      const float xi_22 = 0.10000000000000001f;
#ifdef _OPENMP
      #pragma omp for schedule(static)
#endif
      for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1)
      {
         for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1)
         {
            for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1)
            {
               const float xi_45 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3];
               const float xi_46 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3];
               const float xi_47 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3];
               const float xi_48 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3];
               const float xi_49 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3];
               const float xi_50 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3];
               const float xi_51 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3];
               const float xi_52 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3];
               const float xi_53 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3];
               const float xi_54 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3];
               const float xi_55 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3];
               const float xi_56 = _data_force[_stride_force_0*ctr_0 + _stride_force_1*ctr_1 + _stride_force_2*ctr_2];
               const float xi_57 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3];
               const float xi_58 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3];
               const float xi_59 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3];
               const float xi_60 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3];
               const float xi_61 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3];
               const float xi_62 = _data_force[_stride_force_0*ctr_0 + _stride_force_1*ctr_1 + _stride_force_2*ctr_2 + _stride_force_3];
               const float xi_63 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3];
               const float xi_64 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3];
               const float xi_65 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2];
               const float xi_66 = _data_force[_stride_force_0*ctr_0 + _stride_force_1*ctr_1 + _stride_force_2*ctr_2 + 2*_stride_force_3];
               const float xi_23 = xi_45;
               const float xi_24 = xi_46;
               const float xi_25 = xi_47;
               const float xi_26 = xi_48;
               const float xi_27 = xi_49;
               const float xi_28 = xi_50;
               const float xi_29 = xi_51;
               const float xi_30 = xi_52;
               const float xi_31 = xi_53;
               const float xi_32 = xi_54;
               const float xi_33 = xi_55;
               const float xi_34 = xi_56;
               const float xi_35 = xi_57;
               const float xi_36 = xi_58;
               const float xi_37 = xi_59;
               const float xi_38 = xi_60;
               const float xi_39 = xi_61;
               const float xi_40 = xi_62;
               const float xi_41 = xi_63;
               const float xi_42 = xi_64;
               const float xi_43 = xi_65;
               const float xi_44 = xi_66;
               const float partial_m_m1_0_e_0 = xi_25 + xi_36 + xi_41;
               const float partial_m_m1_e_0_0 = partial_m_m1_0_e_0 + xi_28 + xi_38;
               const float partial_m_0_m1_e_0 = xi_23 + xi_33 + xi_39;
               const float partial_m_0_0_e_0 = xi_26 + xi_42 + xi_43;
               const float partial_m_0_1_e_0 = xi_29 + xi_32 + xi_35;
               const float partial_m_0_e_0_0 = partial_m_0_0_e_0 + partial_m_0_1_e_0 + partial_m_0_m1_e_0;
               const float partial_m_1_0_e_0 = xi_24 + xi_27 + xi_31;
               const float partial_m_1_e_0_0 = partial_m_1_0_e_0 + xi_30 + xi_37;
               const float partial_m_m1_e_1_0 = -xi_28 + xi_38;
               const float partial_m_0_e_1_0 = partial_m_0_1_e_0 - partial_m_0_m1_e_0;
               const float partial_m_1_e_1_0 = xi_30 - xi_37;
               const float partial_m_m1_0_e_1 = -xi_25 + xi_36;
               const float partial_m_0_m1_e_1 = -xi_33 + xi_39;
               const float partial_m_0_0_e_1 = -xi_26 + xi_42;
               const float partial_m_0_1_e_1 = xi_29 - xi_35;
               const float partial_m_0_e_0_1 = partial_m_0_0_e_1 + partial_m_0_1_e_1 + partial_m_0_m1_e_1;
               const float partial_m_1_0_e_1 = -xi_24 + xi_27;
               const float partial_m_m1_e_2_0 = xi_28 + xi_38;
               const float partial_m_0_e_2_0 = partial_m_0_1_e_0 + partial_m_0_m1_e_0;
               const float partial_m_1_e_2_0 = xi_30 + xi_37;
               const float partial_m_m1_0_e_2 = xi_25 + xi_36;
               const float partial_m_0_m1_e_2 = xi_33 + xi_39;
               const float partial_m_0_0_e_2 = xi_26 + xi_42;
               const float partial_m_0_1_e_2 = xi_29 + xi_35;
               const float partial_m_0_e_0_2 = partial_m_0_0_e_2 + partial_m_0_1_e_2 + partial_m_0_m1_e_2;
               const float partial_m_1_0_e_2 = xi_24 + xi_27;
               const float partial_m_0_e_1_1 = partial_m_0_1_e_1 - partial_m_0_m1_e_1;
               const float partial_m_0_e_2_1 = partial_m_0_1_e_1 + partial_m_0_m1_e_1;
               const float partial_m_0_e_1_2 = partial_m_0_1_e_2 - partial_m_0_m1_e_2;
               const float partial_m_0_e_2_2 = partial_m_0_1_e_2 + partial_m_0_m1_e_2;
               const float m_000 = partial_m_0_e_0_0 + partial_m_1_e_0_0 + partial_m_m1_e_0_0;
               const float xi_0 = ((1.0f) / (m_000));
               const float m_100 = partial_m_1_e_0_0 - partial_m_m1_e_0_0;
               const float m_010 = partial_m_0_e_1_0 + partial_m_1_e_1_0 + partial_m_m1_e_1_0;
               const float m_001 = partial_m_0_e_0_1 + partial_m_1_0_e_1 + partial_m_m1_0_e_1;
               const float m_200 = partial_m_1_e_0_0 + partial_m_m1_e_0_0;
               const float m_020 = partial_m_0_e_2_0 + partial_m_1_e_2_0 + partial_m_m1_e_2_0;
               const float m_002 = partial_m_0_e_0_2 + partial_m_1_0_e_2 + partial_m_m1_0_e_2;
               const float m_110 = partial_m_1_e_1_0 - partial_m_m1_e_1_0;
               const float m_101 = partial_m_1_0_e_1 - partial_m_m1_0_e_1;
               const float m_210 = partial_m_1_e_1_0 + partial_m_m1_e_1_0;
               const float m_201 = partial_m_1_0_e_1 + partial_m_m1_0_e_1;
               const float m_120 = partial_m_1_e_2_0 - partial_m_m1_e_2_0;
               const float m_102 = partial_m_1_0_e_2 - partial_m_m1_0_e_2;
               const float m_220 = partial_m_1_e_2_0 + partial_m_m1_e_2_0;
               const float m_202 = partial_m_1_0_e_2 + partial_m_m1_0_e_2;
               const float u_0 = m_100*xi_0 + xi_0*xi_2*xi_34;
               const float u_1 = m_010*xi_0 + xi_0*xi_2*xi_40;
               const float u_2 = m_001*xi_0 + xi_0*xi_2*xi_44;
               const float M_4 = -m_002 - m_020 + m_200*((float)(xi_3));
               const float M_5 = -m_002 + m_020;
               const float M_9 = -m_000 + m_002 + m_020 + m_200;
               const float M_10 = -m_010 + m_210*((float)(xi_4));
               const float M_11 = -m_001 + m_201*((float)(xi_4));
               const float M_12 = -m_100 + m_120*((float)(xi_4));
               const float M_13 = -m_100 + m_102*((float)(xi_3)) + m_120;
               const float M_14 = -m_001 + m_201 + partial_m_0_e_2_1*((float)(xi_3));
               const float M_15 = -m_010 + m_210 + partial_m_0_e_1_2*((float)(xi_3));
               const float M_16 = m_000 + m_002*((float)(xi_4)) + m_220*((float)(xi_6)) + (-m_020 - m_200)*((float)(xi_5));
               const float M_17 = m_000 + m_020 - m_200*((float)(xi_5)) + m_202*((float)(xi_8)) + (-m_002 + m_220)*((float)(xi_7));
               const float M_18 = m_000 - m_200 + partial_m_0_e_2_2*((float)(xi_9)) + (-m_002 - m_020 + m_202 + m_220)*((float)(xi_7));
               const float M_post_1 = m_100 + xi_34;
               const float M_post_2 = m_010 + xi_40;
               const float M_post_3 = m_001 + xi_44;
               const float M_post_4 = -M_4*omega_shear + M_4 - m_000*omega_shear*u_1*u_1 - m_000*omega_shear*u_2*u_2 + (-omega_shear*u_0*xi_2*xi_34 + u_0*xi_34)*((float)(xi_7)) + ((float)(xi_3))*(m_000*omega_shear*(u_0*u_0) - u_1*xi_40 - u_2*xi_44 + xi_2*(omega_shear*u_1*xi_40 + omega_shear*u_2*xi_44));
               const float M_post_5 = -M_5*omega_shear + M_5 + m_000*omega_shear*(u_1*u_1) - m_000*omega_shear*u_2*u_2 + (u_1*xi_40 - u_2*xi_44 + xi_2*(-omega_shear*u_1*xi_40 + omega_shear*u_2*xi_44))*((float)(xi_3));
               const float M_post_6 = m_000*omega_shear*u_0*u_1 - m_110*omega_shear + m_110 + u_0*xi_40 + u_1*xi_34 + xi_2*(-omega_shear*u_0*xi_40 - omega_shear*u_1*xi_34);
               const float M_post_7 = m_000*omega_shear*u_0*u_2 - m_101*omega_shear + m_101 + u_0*xi_44 + u_2*xi_34 + xi_2*(-omega_shear*u_0*xi_44 - omega_shear*u_2*xi_34);
               const float M_post_8 = m_000*omega_shear*u_1*u_2 - omega_shear*partial_m_0_e_1_1 + partial_m_0_e_1_1 + u_1*xi_44 + u_2*xi_40 + xi_2*(-omega_shear*u_1*xi_44 - omega_shear*u_2*xi_40);
               const float M_post_9 = -M_9*omega_bulk + M_9 + m_000*omega_bulk*(u_0*u_0) + m_000*omega_bulk*(u_1*u_1) + m_000*omega_bulk*(u_2*u_2) + (u_0*xi_34 + u_1*xi_40 + u_2*xi_44 + xi_2*(-omega_bulk*u_0*xi_34 - omega_bulk*u_1*xi_40 - omega_bulk*u_2*xi_44))*((float)(xi_3));
               const float M_post_10 = -M_10*omega_odd + M_10;
               const float M_post_11 = -M_11*omega_odd + M_11;
               const float M_post_12 = -M_12*omega_odd + M_12;
               const float M_post_13 = -M_13*omega_odd + M_13;
               const float M_post_14 = -M_14*omega_odd + M_14;
               const float M_post_15 = -M_15*omega_odd + M_15;
               const float M_post_16 = -M_16*omega_even + M_16 + m_000*omega_even*((float)(xi_4))*(u_2*u_2);
               const float M_post_17 = -M_17*omega_even + M_17 + m_000*omega_even*xi_10*(u_2*u_2) + m_000*omega_even*xi_11*(u_1*u_1);
               const float M_post_18 = -M_18*omega_even + M_18 + m_000*omega_even*xi_12*(u_0*u_0) + xi_10*(m_000*omega_even*(u_1*u_1) + m_000*omega_even*(u_2*u_2));
               const float m_post_200 = xi_13*(M_post_4 + M_post_9 + m_000);
               const float m_post_020 = -M_post_4*xi_14 + M_post_5*xi_2 + xi_13*(M_post_9 + m_000);
               const float m_post_002 = -M_post_4*xi_14 - M_post_5*xi_2 + xi_13*(M_post_9 + m_000);
               const float m_post_210 = xi_13*(M_post_10 + M_post_2);
               const float m_post_201 = xi_13*(M_post_11 + M_post_3);
               const float m_post_120 = xi_13*(M_post_1 + M_post_12);
               const float m_post_102 = M_post_1*xi_13 - M_post_12*xi_14 + M_post_13*xi_2;
               const float m_post_021 = -M_post_11*xi_14 + M_post_14*xi_2 + M_post_3*xi_13;
               const float m_post_012 = -M_post_10*xi_14 + M_post_15*xi_2 + M_post_2*xi_13;
               const float m_post_220 = M_post_16*xi_18 + M_post_4*xi_17 + M_post_5*xi_15 + M_post_9*xi_14 + m_000*xi_16;
               const float m_post_202 = -M_post_16*xi_19 + M_post_17*xi_20 + M_post_4*xi_17 - M_post_5*xi_15 + M_post_9*xi_14 + m_000*xi_16;
               const float m_post_022 = -M_post_16*xi_19 - M_post_17*xi_21 + M_post_18*xi_22 + m_000*xi_16 + xi_14*(-M_post_4 + M_post_9);
               const float sub_k_to_f_20 = xi_2*(m_post_020 - m_post_022 - m_post_220);
               const float sub_k_to_f_21 = xi_2*(M_post_2 - m_post_012 - m_post_210);
               const float sub_k_to_f_22 = xi_2*(m_post_200 - m_post_202 - m_post_220);
               const float sub_k_to_f_23 = xi_2*(-M_post_1 + m_post_102 + m_post_120);
               const float sub_k_to_f_24 = xi_2*(m_post_002 - m_post_022 - m_post_202);
               const float sub_k_to_f_25 = xi_2*(M_post_3 - m_post_021 - m_post_201);
               const float sub_k_to_f_26 = xi_15*(-M_post_6 + m_post_220);
               const float sub_k_to_f_27 = xi_15*(-m_post_120 + m_post_210);
               const float sub_k_to_f_28 = xi_15*(M_post_6 + m_post_220);
               const float sub_k_to_f_29 = xi_15*(m_post_120 + m_post_210);
               const float sub_k_to_f_30 = xi_15*(M_post_8 + m_post_022);
               const float sub_k_to_f_31 = xi_15*(m_post_012 + m_post_021);
               const float sub_k_to_f_32 = xi_15*(-M_post_8 + m_post_022);
               const float sub_k_to_f_33 = xi_15*(-m_post_012 + m_post_021);
               const float sub_k_to_f_34 = xi_15*(-M_post_7 + m_post_202);
               const float sub_k_to_f_35 = xi_15*(-m_post_102 + m_post_201);
               const float sub_k_to_f_36 = xi_15*(M_post_7 + m_post_202);
               const float sub_k_to_f_37 = xi_15*(m_post_102 + m_post_201);
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2] = m_000 - m_post_002 - m_post_020 + m_post_022 - m_post_200 + m_post_202 + m_post_220;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3] = sub_k_to_f_20 + sub_k_to_f_21 + ((-1.0f <= -grid_size + ((float)(ctr_1))) ? (rho*v_s*(u_0*2.0f + v_s)*0.16666666666666666f): (0.0f));
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3] = sub_k_to_f_20 - sub_k_to_f_21 + ((0.0f >= ((float)(ctr_1))) ? (rho*v_s*(u_0*-2.0f + v_s)*0.16666666666666666f): (0.0f));
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3] = sub_k_to_f_22 + sub_k_to_f_23;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3] = sub_k_to_f_22 - sub_k_to_f_23;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3] = sub_k_to_f_24 + sub_k_to_f_25;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3] = sub_k_to_f_24 - sub_k_to_f_25;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3] = sub_k_to_f_26 + sub_k_to_f_27 + ((-1.0f <= -grid_size + ((float)(ctr_1))) ? (rho*v_s*(u_0*-2.0f + u_1*3.0f - v_s + 1.0f)*0.083333333333333329f): (0.0f));
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3] = sub_k_to_f_28 + sub_k_to_f_29 + ((-1.0f <= -grid_size + ((float)(ctr_1))) ? (rho*v_s*(u_0*-2.0f + u_1*-3.0f - v_s - 1.0f)*0.083333333333333329f): (0.0f));
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3] = sub_k_to_f_28 - sub_k_to_f_29 + ((0.0f >= ((float)(ctr_1))) ? (rho*v_s*(u_0*2.0f + u_1*3.0f - v_s - 1.0f)*0.083333333333333329f): (0.0f));
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3] = sub_k_to_f_26 - sub_k_to_f_27 + ((0.0f >= ((float)(ctr_1))) ? (rho*v_s*(u_0*2.0f + u_1*-3.0f - v_s + 1.0f)*0.083333333333333329f): (0.0f));
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3] = sub_k_to_f_30 + sub_k_to_f_31;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3] = sub_k_to_f_32 + sub_k_to_f_33;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3] = sub_k_to_f_34 + sub_k_to_f_35;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3] = sub_k_to_f_36 + sub_k_to_f_37;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3] = sub_k_to_f_32 - sub_k_to_f_33;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3] = sub_k_to_f_30 - sub_k_to_f_31;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3] = sub_k_to_f_36 - sub_k_to_f_37;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3] = sub_k_to_f_34 - sub_k_to_f_35;
            }
         }
      }
   }
}
}


void CollideSweepSinglePrecisionLeesEdwards::run(IBlock * block)
{
   
    auto force = block->getData< field::GhostLayerField<float, 3> >(forceID);
    auto pdfs = block->getData< field::GhostLayerField<float, 19> >(pdfsID);

    auto & grid_size = this->grid_size_;
    auto & omega_bulk = this->omega_bulk_;
    auto & omega_shear = this->omega_shear_;
    auto & omega_odd = this->omega_odd_;
    auto & omega_even = this->omega_even_;
    auto & v_s = this->v_s_;
    auto & rho = this->rho_;
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()))
    float * RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(int64_c(force->xSize()) + 0))
    const int64_t _size_force_0 = int64_t(int64_c(force->xSize()) + 0);
    WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(int64_c(force->ySize()) + 0))
    const int64_t _size_force_1 = int64_t(int64_c(force->ySize()) + 0);
    WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(int64_c(force->zSize()) + 0))
    const int64_t _size_force_2 = int64_t(int64_c(force->zSize()) + 0);
    WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
    const int64_t _stride_force_0 = int64_t(force->xStride());
    const int64_t _stride_force_1 = int64_t(force->yStride());
    const int64_t _stride_force_2 = int64_t(force->zStride());
    const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_ab1f3bc3368574afb482da84ccb58898::collidesweepsingleprecisionleesedwards_collidesweepsingleprecisionleesedwards(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, grid_size, omega_bulk, omega_even, omega_odd, omega_shear, rho, v_s);
    
}


void CollideSweepSinglePrecisionLeesEdwards::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block)
{
   
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    auto force = block->getData< field::GhostLayerField<float, 3> >(forceID);
    auto pdfs = block->getData< field::GhostLayerField<float, 19> >(pdfsID);

    auto & grid_size = this->grid_size_;
    auto & omega_bulk = this->omega_bulk_;
    auto & omega_shear = this->omega_shear_;
    auto & omega_odd = this->omega_odd_;
    auto & omega_even = this->omega_even_;
    auto & v_s = this->v_s_;
    auto & rho = this->rho_;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()))
    float * RESTRICT const _data_force = force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_force_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_force_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_force_2 = int64_t(int64_c(ci.zSize()) + 0);
    WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
    const int64_t _stride_force_0 = int64_t(force->xStride());
    const int64_t _stride_force_1 = int64_t(force->yStride());
    const int64_t _stride_force_2 = int64_t(force->zStride());
    const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_ab1f3bc3368574afb482da84ccb58898::collidesweepsingleprecisionleesedwards_collidesweepsingleprecisionleesedwards(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, grid_size, omega_bulk, omega_even, omega_odd, omega_shear, rho, v_s);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
