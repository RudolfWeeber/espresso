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
//! \\file CollideSweepDoublePrecisionLeesEdwards.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7, lbmpy v1.3.7, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 59c9b8b185782eba184e0fdfb2144793343213f0


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "CollideSweepDoublePrecisionLeesEdwards.h"




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


namespace internal_607d8a5c7ac58c25acf09ad94bb82cf4 {
static FUNC_PREFIX void collidesweepdoubleprecisionleesedwards_collidesweepdoubleprecisionleesedwards(double * RESTRICT const _data_force, double * RESTRICT  _data_pdfs, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_0, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, double grid_size, double omega_bulk, double omega_even, double omega_odd, double omega_shear, double rho, double v_s)
{
#ifdef _OPENMP
   #pragma omp parallel
#endif
   {
      const double xi_2 = 0.50000000000000000;
      const int64_t xi_3 = 2;
      const int64_t xi_4 = 3;
      const int64_t xi_5 = 6;
      const int64_t xi_6 = 18;
      const int64_t xi_7 = 4;
      const int64_t xi_8 = 14;
      const int64_t xi_9 = 10;
      const double xi_10 = 0.66666666666666663;
      const double xi_11 = 2.3333333333333335;
      const double xi_12 = 1.6666666666666667;
      const double xi_13 = 0.33333333333333331;
      const double xi_14 = 0.16666666666666666;
      const double xi_15 = 0.25000000000000000;
      const double xi_16 = 0.1111111111111111;
      const double xi_17 = 0.083333333333333329;
      const double xi_18 = 0.055555555555555552;
      const double xi_19 = 0.015873015873015872;
      const double xi_20 = 0.071428571428571425;
      const double xi_21 = 0.028571428571428571;
      const double xi_22 = 0.10000000000000001;
#ifdef _OPENMP
      #pragma omp for schedule(static)
#endif
      for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1)
      {
         for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1)
         {
            for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1)
            {
               const double xi_45 = _data_force[_stride_force_0*ctr_0 + _stride_force_1*ctr_1 + _stride_force_2*ctr_2 + 2*_stride_force_3];
               const double xi_46 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3];
               const double xi_47 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3];
               const double xi_48 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3];
               const double xi_49 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3];
               const double xi_50 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3];
               const double xi_51 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3];
               const double xi_52 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3];
               const double xi_53 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3];
               const double xi_54 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2];
               const double xi_55 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3];
               const double xi_56 = _data_force[_stride_force_0*ctr_0 + _stride_force_1*ctr_1 + _stride_force_2*ctr_2 + _stride_force_3];
               const double xi_57 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3];
               const double xi_58 = _data_force[_stride_force_0*ctr_0 + _stride_force_1*ctr_1 + _stride_force_2*ctr_2];
               const double xi_59 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3];
               const double xi_60 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3];
               const double xi_61 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3];
               const double xi_62 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3];
               const double xi_63 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3];
               const double xi_64 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3];
               const double xi_65 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3];
               const double xi_66 = _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3];
               const double xi_23 = xi_56;
               const double xi_24 = xi_66;
               const double xi_25 = xi_65;
               const double xi_26 = xi_53;
               const double xi_27 = xi_45;
               const double xi_28 = xi_52;
               const double xi_29 = xi_64;
               const double xi_30 = xi_47;
               const double xi_31 = xi_59;
               const double xi_32 = xi_51;
               const double xi_33 = xi_63;
               const double xi_34 = xi_54;
               const double xi_35 = xi_61;
               const double xi_36 = xi_49;
               const double xi_37 = xi_57;
               const double xi_38 = xi_48;
               const double xi_39 = xi_58;
               const double xi_40 = xi_62;
               const double xi_41 = xi_55;
               const double xi_42 = xi_50;
               const double xi_43 = xi_46;
               const double xi_44 = xi_60;
               const double partial_m_m1_0_e_0 = xi_29 + xi_37 + xi_42;
               const double partial_m_m1_e_0_0 = partial_m_m1_0_e_0 + xi_30 + xi_38;
               const double partial_m_0_m1_e_0 = xi_24 + xi_25 + xi_40;
               const double partial_m_0_0_e_0 = xi_34 + xi_35 + xi_36;
               const double partial_m_0_1_e_0 = xi_26 + xi_28 + xi_43;
               const double partial_m_0_e_0_0 = partial_m_0_0_e_0 + partial_m_0_1_e_0 + partial_m_0_m1_e_0;
               const double partial_m_1_0_e_0 = xi_31 + xi_41 + xi_44;
               const double partial_m_1_e_0_0 = partial_m_1_0_e_0 + xi_32 + xi_33;
               const double partial_m_m1_e_1_0 = -xi_30 + xi_38;
               const double partial_m_0_e_1_0 = partial_m_0_1_e_0 - partial_m_0_m1_e_0;
               const double partial_m_1_e_1_0 = xi_32 - xi_33;
               const double partial_m_m1_0_e_1 = -xi_37 + xi_42;
               const double partial_m_0_m1_e_1 = xi_25 - xi_40;
               const double partial_m_0_0_e_1 = xi_35 - xi_36;
               const double partial_m_0_1_e_1 = xi_26 - xi_43;
               const double partial_m_0_e_0_1 = partial_m_0_0_e_1 + partial_m_0_1_e_1 + partial_m_0_m1_e_1;
               const double partial_m_1_0_e_1 = -xi_31 + xi_44;
               const double partial_m_m1_e_2_0 = xi_30 + xi_38;
               const double partial_m_0_e_2_0 = partial_m_0_1_e_0 + partial_m_0_m1_e_0;
               const double partial_m_1_e_2_0 = xi_32 + xi_33;
               const double partial_m_m1_0_e_2 = xi_37 + xi_42;
               const double partial_m_0_m1_e_2 = xi_25 + xi_40;
               const double partial_m_0_0_e_2 = xi_35 + xi_36;
               const double partial_m_0_1_e_2 = xi_26 + xi_43;
               const double partial_m_0_e_0_2 = partial_m_0_0_e_2 + partial_m_0_1_e_2 + partial_m_0_m1_e_2;
               const double partial_m_1_0_e_2 = xi_31 + xi_44;
               const double partial_m_0_e_1_1 = partial_m_0_1_e_1 - partial_m_0_m1_e_1;
               const double partial_m_0_e_2_1 = partial_m_0_1_e_1 + partial_m_0_m1_e_1;
               const double partial_m_0_e_1_2 = partial_m_0_1_e_2 - partial_m_0_m1_e_2;
               const double partial_m_0_e_2_2 = partial_m_0_1_e_2 + partial_m_0_m1_e_2;
               const double m_000 = partial_m_0_e_0_0 + partial_m_1_e_0_0 + partial_m_m1_e_0_0;
               const double xi_0 = ((1.0) / (m_000));
               const double m_100 = partial_m_1_e_0_0 - partial_m_m1_e_0_0;
               const double m_010 = partial_m_0_e_1_0 + partial_m_1_e_1_0 + partial_m_m1_e_1_0;
               const double m_001 = partial_m_0_e_0_1 + partial_m_1_0_e_1 + partial_m_m1_0_e_1;
               const double m_200 = partial_m_1_e_0_0 + partial_m_m1_e_0_0;
               const double m_020 = partial_m_0_e_2_0 + partial_m_1_e_2_0 + partial_m_m1_e_2_0;
               const double m_002 = partial_m_0_e_0_2 + partial_m_1_0_e_2 + partial_m_m1_0_e_2;
               const double m_110 = partial_m_1_e_1_0 - partial_m_m1_e_1_0;
               const double m_101 = partial_m_1_0_e_1 - partial_m_m1_0_e_1;
               const double m_210 = partial_m_1_e_1_0 + partial_m_m1_e_1_0;
               const double m_201 = partial_m_1_0_e_1 + partial_m_m1_0_e_1;
               const double m_120 = partial_m_1_e_2_0 - partial_m_m1_e_2_0;
               const double m_102 = partial_m_1_0_e_2 - partial_m_m1_0_e_2;
               const double m_220 = partial_m_1_e_2_0 + partial_m_m1_e_2_0;
               const double m_202 = partial_m_1_0_e_2 + partial_m_m1_0_e_2;
               const double u_0 = m_100*xi_0 + xi_0*xi_2*xi_39;
               const double u_1 = m_010*xi_0 + xi_0*xi_2*xi_23;
               const double u_2 = m_001*xi_0 + xi_0*xi_2*xi_27;
               const double M_4 = -m_002 - m_020 + m_200*((double)(xi_3));
               const double M_5 = -m_002 + m_020;
               const double M_9 = -m_000 + m_002 + m_020 + m_200;
               const double M_10 = -m_010 + m_210*((double)(xi_4));
               const double M_11 = -m_001 + m_201*((double)(xi_4));
               const double M_12 = -m_100 + m_120*((double)(xi_4));
               const double M_13 = -m_100 + m_102*((double)(xi_3)) + m_120;
               const double M_14 = -m_001 + m_201 + partial_m_0_e_2_1*((double)(xi_3));
               const double M_15 = -m_010 + m_210 + partial_m_0_e_1_2*((double)(xi_3));
               const double M_16 = m_000 + m_002*((double)(xi_4)) + m_220*((double)(xi_6)) + (-m_020 - m_200)*((double)(xi_5));
               const double M_17 = m_000 + m_020 - m_200*((double)(xi_5)) + m_202*((double)(xi_8)) + (-m_002 + m_220)*((double)(xi_7));
               const double M_18 = m_000 - m_200 + partial_m_0_e_2_2*((double)(xi_9)) + (-m_002 - m_020 + m_202 + m_220)*((double)(xi_7));
               const double M_post_1 = m_100 + xi_39;
               const double M_post_2 = m_010 + xi_23;
               const double M_post_3 = m_001 + xi_27;
               const double M_post_4 = -M_4*omega_shear + M_4 - m_000*omega_shear*u_1*u_1 - m_000*omega_shear*u_2*u_2 + (-omega_shear*u_0*xi_2*xi_39 + u_0*xi_39)*((double)(xi_7)) + ((double)(xi_3))*(m_000*omega_shear*(u_0*u_0) - u_1*xi_23 - u_2*xi_27 + xi_2*(omega_shear*u_1*xi_23 + omega_shear*u_2*xi_27));
               const double M_post_5 = -M_5*omega_shear + M_5 + m_000*omega_shear*(u_1*u_1) - m_000*omega_shear*u_2*u_2 + (u_1*xi_23 - u_2*xi_27 + xi_2*(-omega_shear*u_1*xi_23 + omega_shear*u_2*xi_27))*((double)(xi_3));
               const double M_post_6 = m_000*omega_shear*u_0*u_1 - m_110*omega_shear + m_110 + u_0*xi_23 + u_1*xi_39 + xi_2*(-omega_shear*u_0*xi_23 - omega_shear*u_1*xi_39);
               const double M_post_7 = m_000*omega_shear*u_0*u_2 - m_101*omega_shear + m_101 + u_0*xi_27 + u_2*xi_39 + xi_2*(-omega_shear*u_0*xi_27 - omega_shear*u_2*xi_39);
               const double M_post_8 = m_000*omega_shear*u_1*u_2 - omega_shear*partial_m_0_e_1_1 + partial_m_0_e_1_1 + u_1*xi_27 + u_2*xi_23 + xi_2*(-omega_shear*u_1*xi_27 - omega_shear*u_2*xi_23);
               const double M_post_9 = -M_9*omega_bulk + M_9 + m_000*omega_bulk*(u_0*u_0) + m_000*omega_bulk*(u_1*u_1) + m_000*omega_bulk*(u_2*u_2) + (u_0*xi_39 + u_1*xi_23 + u_2*xi_27 + xi_2*(-omega_bulk*u_0*xi_39 - omega_bulk*u_1*xi_23 - omega_bulk*u_2*xi_27))*((double)(xi_3));
               const double M_post_10 = -M_10*omega_odd + M_10;
               const double M_post_11 = -M_11*omega_odd + M_11;
               const double M_post_12 = -M_12*omega_odd + M_12;
               const double M_post_13 = -M_13*omega_odd + M_13;
               const double M_post_14 = -M_14*omega_odd + M_14;
               const double M_post_15 = -M_15*omega_odd + M_15;
               const double M_post_16 = -M_16*omega_even + M_16 + m_000*omega_even*((double)(xi_4))*(u_2*u_2);
               const double M_post_17 = -M_17*omega_even + M_17 + m_000*omega_even*xi_10*(u_2*u_2) + m_000*omega_even*xi_11*(u_1*u_1);
               const double M_post_18 = -M_18*omega_even + M_18 + m_000*omega_even*xi_12*(u_0*u_0) + xi_10*(m_000*omega_even*(u_1*u_1) + m_000*omega_even*(u_2*u_2));
               const double m_post_200 = xi_13*(M_post_4 + M_post_9 + m_000);
               const double m_post_020 = -M_post_4*xi_14 + M_post_5*xi_2 + xi_13*(M_post_9 + m_000);
               const double m_post_002 = -M_post_4*xi_14 - M_post_5*xi_2 + xi_13*(M_post_9 + m_000);
               const double m_post_210 = xi_13*(M_post_10 + M_post_2);
               const double m_post_201 = xi_13*(M_post_11 + M_post_3);
               const double m_post_120 = xi_13*(M_post_1 + M_post_12);
               const double m_post_102 = M_post_1*xi_13 - M_post_12*xi_14 + M_post_13*xi_2;
               const double m_post_021 = -M_post_11*xi_14 + M_post_14*xi_2 + M_post_3*xi_13;
               const double m_post_012 = -M_post_10*xi_14 + M_post_15*xi_2 + M_post_2*xi_13;
               const double m_post_220 = M_post_16*xi_18 + M_post_4*xi_17 + M_post_5*xi_15 + M_post_9*xi_14 + m_000*xi_16;
               const double m_post_202 = -M_post_16*xi_19 + M_post_17*xi_20 + M_post_4*xi_17 - M_post_5*xi_15 + M_post_9*xi_14 + m_000*xi_16;
               const double m_post_022 = -M_post_16*xi_19 - M_post_17*xi_21 + M_post_18*xi_22 + m_000*xi_16 + xi_14*(-M_post_4 + M_post_9);
               const double sub_k_to_f_20 = xi_2*(m_post_020 - m_post_022 - m_post_220);
               const double sub_k_to_f_21 = xi_2*(M_post_2 - m_post_012 - m_post_210);
               const double sub_k_to_f_22 = xi_2*(m_post_200 - m_post_202 - m_post_220);
               const double sub_k_to_f_23 = xi_2*(-M_post_1 + m_post_102 + m_post_120);
               const double sub_k_to_f_24 = xi_2*(m_post_002 - m_post_022 - m_post_202);
               const double sub_k_to_f_25 = xi_2*(M_post_3 - m_post_021 - m_post_201);
               const double sub_k_to_f_26 = xi_15*(-M_post_6 + m_post_220);
               const double sub_k_to_f_27 = xi_15*(-m_post_120 + m_post_210);
               const double sub_k_to_f_28 = xi_15*(M_post_6 + m_post_220);
               const double sub_k_to_f_29 = xi_15*(m_post_120 + m_post_210);
               const double sub_k_to_f_30 = xi_15*(M_post_8 + m_post_022);
               const double sub_k_to_f_31 = xi_15*(m_post_012 + m_post_021);
               const double sub_k_to_f_32 = xi_15*(-M_post_8 + m_post_022);
               const double sub_k_to_f_33 = xi_15*(-m_post_012 + m_post_021);
               const double sub_k_to_f_34 = xi_15*(-M_post_7 + m_post_202);
               const double sub_k_to_f_35 = xi_15*(-m_post_102 + m_post_201);
               const double sub_k_to_f_36 = xi_15*(M_post_7 + m_post_202);
               const double sub_k_to_f_37 = xi_15*(m_post_102 + m_post_201);
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2] = m_000 - m_post_002 - m_post_020 + m_post_022 - m_post_200 + m_post_202 + m_post_220;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + _stride_pdfs_3] = sub_k_to_f_20 + sub_k_to_f_21 + ((-1.0 <= -grid_size + ((double)(ctr_1))) ? (rho*v_s*(u_0*2.0 + v_s)*0.16666666666666666): (0.0));
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3] = sub_k_to_f_20 - sub_k_to_f_21 + ((0.0 >= ((double)(ctr_1))) ? (rho*v_s*(u_0*-2.0 + v_s)*0.16666666666666666): (0.0));
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3] = sub_k_to_f_22 + sub_k_to_f_23;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3] = sub_k_to_f_22 - sub_k_to_f_23;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3] = sub_k_to_f_24 + sub_k_to_f_25;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3] = sub_k_to_f_24 - sub_k_to_f_25;
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3] = sub_k_to_f_26 + sub_k_to_f_27 + ((-1.0 <= -grid_size + ((double)(ctr_1))) ? (rho*v_s*(u_0*-2.0 + u_1*3.0 - v_s + 1.0)*0.083333333333333329): (0.0));
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3] = sub_k_to_f_28 + sub_k_to_f_29 + ((-1.0 <= -grid_size + ((double)(ctr_1))) ? (rho*v_s*(u_0*-2.0 + u_1*-3.0 - v_s - 1.0)*0.083333333333333329): (0.0));
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3] = sub_k_to_f_28 - sub_k_to_f_29 + ((0.0 >= ((double)(ctr_1))) ? (rho*v_s*(u_0*2.0 + u_1*3.0 - v_s - 1.0)*0.083333333333333329): (0.0));
               _data_pdfs[_stride_pdfs_0*ctr_0 + _stride_pdfs_1*ctr_1 + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3] = sub_k_to_f_26 - sub_k_to_f_27 + ((0.0 >= ((double)(ctr_1))) ? (rho*v_s*(u_0*2.0 + u_1*-3.0 - v_s + 1.0)*0.083333333333333329): (0.0));
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


void CollideSweepDoublePrecisionLeesEdwards::run(IBlock * block)
{
   
    auto force = block->getData< field::GhostLayerField<double, 3> >(forceID);
    auto pdfs = block->getData< field::GhostLayerField<double, 19> >(pdfsID);

    auto & omega_shear = this->omega_shear_;
    auto & rho = this->rho_;
    auto & v_s = this->v_s_;
    auto & omega_odd = this->omega_odd_;
    auto & grid_size = this->grid_size_;
    auto & omega_bulk = this->omega_bulk_;
    auto & omega_even = this->omega_even_;
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()))
    double * RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
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
    internal_607d8a5c7ac58c25acf09ad94bb82cf4::collidesweepdoubleprecisionleesedwards_collidesweepdoubleprecisionleesedwards(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, grid_size, omega_bulk, omega_even, omega_odd, omega_shear, rho, v_s);
    
}


void CollideSweepDoublePrecisionLeesEdwards::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block)
{
   
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    auto force = block->getData< field::GhostLayerField<double, 3> >(forceID);
    auto pdfs = block->getData< field::GhostLayerField<double, 19> >(pdfsID);

    auto & omega_shear = this->omega_shear_;
    auto & rho = this->rho_;
    auto & v_s = this->v_s_;
    auto & omega_odd = this->omega_odd_;
    auto & grid_size = this->grid_size_;
    auto & omega_bulk = this->omega_bulk_;
    auto & omega_even = this->omega_even_;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()))
    double * RESTRICT const _data_force = force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
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
    internal_607d8a5c7ac58c25acf09ad94bb82cf4::collidesweepdoubleprecisionleesedwards_collidesweepdoubleprecisionleesedwards(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, grid_size, omega_bulk, omega_even, omega_odd, omega_shear, rho, v_s);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
