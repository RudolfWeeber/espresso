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
//! \\author Martin Bauer <martin.bauer@fau.de>
//======================================================================================================================

#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "lbm/field/PdfField.h"
#include "lbm/sweeps/Streaming.h"
#include "MRTLatticeModelAvx.h"

#ifdef _MSC_VER
#  pragma warning( disable : 4458 )
#endif

#define FUNC_PREFIX

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wunused-variable"
#   pragma GCC diagnostic ignored "-Wunused-parameter"
#endif


#include <immintrin.h>



using namespace std;

namespace walberla {
namespace lbm {

namespace internal_kernel_streamCollide {
static FUNC_PREFIX void kernel_streamCollide(double * RESTRICT const _data_force, double * RESTRICT const _data_pdfs, double * RESTRICT _data_pdfs_tmp, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3, double omega_bulk, double omega_even, double omega_odd, double omega_shear)
{
   const double xi_50 = -omega_shear + 2.0;
   const double xi_51 = xi_50*0.5;
   const double xi_56 = xi_50*0.0833333333333333;
   const double xi_61 = xi_50*0.166666666666667;
   const double xi_71 = xi_50*0.25;
   const double xi_76 = xi_50*0.0416666666666667;
   const double xi_120 = omega_odd*0.25;
   const double xi_126 = omega_odd*0.0833333333333333;
   const double xi_159 = omega_shear*0.25;
   const double xi_181 = omega_odd*0.0416666666666667;
   const double xi_183 = omega_odd*0.125;
   const int64_t rr_0 = 0.0;
   const double xi_131 = rr_0*0.166666666666667;
   const double xi_164 = rr_0*0.0833333333333333;
   for (int64_t ctr_2 = 1; ctr_2 < _size_force_2 - 1; ctr_2 += 1)
   {
      double * RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3;
      double * RESTRICT _data_force_20_31 = _data_force + _stride_force_2*ctr_2 + _stride_force_3;
      double * RESTRICT _data_force_20_30 = _data_force + _stride_force_2*ctr_2;
      double * RESTRICT _data_force_20_32 = _data_force + _stride_force_2*ctr_2 + 2*_stride_force_3;
      double * RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * RESTRICT _data_pdfs_tmp_20_30 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2;
      double * RESTRICT _data_pdfs_tmp_20_31 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + _stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_32 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 2*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_33 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 3*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_34 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 4*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_35 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 5*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_36 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 6*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_37 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 7*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_38 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 8*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_39 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 9*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_310 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 10*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_311 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 11*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_312 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 12*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_313 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 13*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_314 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 14*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_315 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 15*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_316 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 16*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_317 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 17*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_318 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 18*_stride_pdfs_tmp_3;
      for (int64_t ctr_1 = 1; ctr_1 < _size_force_1 - 1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_20_310_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
         double * RESTRICT _data_pdfs_20_38_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
         double * RESTRICT _data_pdfs_2m1_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_314;
         double * RESTRICT _data_pdfs_21_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_318;
         double * RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * RESTRICT _data_pdfs_2m1_311_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
         double * RESTRICT _data_pdfs_21_315_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
         double * RESTRICT _data_pdfs_20_31_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
         double * RESTRICT _data_pdfs_20_37_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
         double * RESTRICT _data_pdfs_2m1_312_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
         double * RESTRICT _data_pdfs_2m1_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_35;
         double * RESTRICT _data_pdfs_2m1_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_313;
         double * RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * RESTRICT _data_pdfs_20_39_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
         double * RESTRICT _data_pdfs_20_32_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
         double * RESTRICT _data_pdfs_21_316_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
         double * RESTRICT _data_pdfs_21_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_317;
         double * RESTRICT _data_pdfs_21_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_36;
         double * RESTRICT _data_force_20_31_10 = _stride_force_1*ctr_1 + _data_force_20_31;
         double * RESTRICT _data_force_20_30_10 = _stride_force_1*ctr_1 + _data_force_20_30;
         double * RESTRICT _data_force_20_32_10 = _stride_force_1*ctr_1 + _data_force_20_32;
         double * RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * RESTRICT _data_pdfs_tmp_20_30_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_30;
         double * RESTRICT _data_pdfs_tmp_20_31_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_31;
         double * RESTRICT _data_pdfs_tmp_20_32_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_32;
         double * RESTRICT _data_pdfs_tmp_20_33_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_33;
         double * RESTRICT _data_pdfs_tmp_20_34_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_34;
         double * RESTRICT _data_pdfs_tmp_20_35_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_35;
         double * RESTRICT _data_pdfs_tmp_20_36_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_36;
         double * RESTRICT _data_pdfs_tmp_20_37_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_37;
         double * RESTRICT _data_pdfs_tmp_20_38_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_38;
         double * RESTRICT _data_pdfs_tmp_20_39_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_39;
         double * RESTRICT _data_pdfs_tmp_20_310_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_310;
         double * RESTRICT _data_pdfs_tmp_20_311_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_311;
         double * RESTRICT _data_pdfs_tmp_20_312_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_312;
         double * RESTRICT _data_pdfs_tmp_20_313_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_313;
         double * RESTRICT _data_pdfs_tmp_20_314_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_314;
         double * RESTRICT _data_pdfs_tmp_20_315_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_315;
         double * RESTRICT _data_pdfs_tmp_20_316_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_316;
         double * RESTRICT _data_pdfs_tmp_20_317_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_317;
         double * RESTRICT _data_pdfs_tmp_20_318_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_318;
         for (int64_t ctr_0 = 1; ctr_0 < ((_size_force_0 - 2) % (4) == 0 ? _size_force_0 - 2 : ((int64_t)((_size_force_0 - 2) / (4))+1) * (4)) + 1; ctr_0 += 4)
         {
            const __m256d xi_0 = _mm256_add_pd(_mm256_loadu_pd(& _data_pdfs_20_310_11[ctr_0 - 1]),_mm256_loadu_pd(& _data_pdfs_20_38_1m1[ctr_0 - 1]));
            const __m256d xi_1 = _mm256_add_pd(_mm256_loadu_pd(& _data_pdfs_21_318_10[ctr_0 - 1]),_mm256_loadu_pd(& _data_pdfs_2m1_314_10[ctr_0 - 1]));
            const __m256d xi_2 = _mm256_add_pd(_mm256_add_pd(xi_0,xi_1),_mm256_loadu_pd(& _data_pdfs_20_34_10[ctr_0 - 1]));
            const __m256d xi_3 = _mm256_add_pd(_mm256_load_pd(& _data_pdfs_21_315_1m1[ctr_0]),_mm256_load_pd(& _data_pdfs_2m1_311_1m1[ctr_0]));
            const __m256d xi_4 = _mm256_add_pd(xi_3,_mm256_load_pd(& _data_pdfs_20_31_1m1[ctr_0]));
            const __m256d xi_5 = _mm256_add_pd(xi_4,_mm256_loadu_pd(& _data_pdfs_20_37_1m1[ctr_0 + 1]));
            const __m256d xi_6 = _mm256_add_pd(_mm256_load_pd(& _data_pdfs_2m1_312_11[ctr_0]),_mm256_load_pd(& _data_pdfs_2m1_35_10[ctr_0]));
            const __m256d xi_7 = _mm256_add_pd(xi_6,_mm256_loadu_pd(& _data_pdfs_2m1_313_10[ctr_0 + 1]));
            const __m256d xi_8 = _mm256_add_pd(_mm256_loadu_pd(& _data_pdfs_20_33_10[ctr_0 + 1]),_mm256_loadu_pd(& _data_pdfs_20_39_11[ctr_0 + 1]));
            const __m256d xi_9 = _mm256_add_pd(_mm256_load_pd(& _data_pdfs_20_32_11[ctr_0]),_mm256_load_pd(& _data_pdfs_21_316_11[ctr_0]));
            const __m256d xi_10 = _mm256_add_pd(_mm256_loadu_pd(& _data_pdfs_21_317_10[ctr_0 + 1]),_mm256_load_pd(& _data_pdfs_21_36_10[ctr_0]));
            const __m256d xi_12 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_loadu_pd(& _data_pdfs_21_317_10[ctr_0 + 1]));
            const __m256d xi_13 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_loadu_pd(& _data_pdfs_20_33_10[ctr_0 + 1]));
            const __m256d xi_14 = _mm256_add_pd(xi_12,xi_13);
            const __m256d xi_15 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_loadu_pd(& _data_pdfs_2m1_313_10[ctr_0 + 1]));
            const __m256d xi_16 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_loadu_pd(& _data_pdfs_20_37_1m1[ctr_0 + 1]));
            const __m256d xi_17 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_loadu_pd(& _data_pdfs_20_39_11[ctr_0 + 1]));
            const __m256d xi_18 = _mm256_add_pd(xi_16,xi_17);
            const __m256d xi_19 = _mm256_add_pd(xi_15,xi_18);
            const __m256d xi_20 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_loadu_pd(& _data_pdfs_20_310_11[ctr_0 - 1]));
            const __m256d xi_21 = _mm256_add_pd(xi_20,_mm256_loadu_pd(& _data_pdfs_20_38_1m1[ctr_0 - 1]));
            const __m256d xi_22 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_load_pd(& _data_pdfs_2m1_312_11[ctr_0]));
            const __m256d xi_23 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_load_pd(& _data_pdfs_20_32_11[ctr_0]));
            const __m256d xi_24 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_load_pd(& _data_pdfs_21_316_11[ctr_0]));
            const __m256d xi_25 = _mm256_add_pd(_mm256_add_pd(xi_22,xi_23),xi_24);
            const __m256d xi_26 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_loadu_pd(& _data_pdfs_21_318_10[ctr_0 - 1]));
            const __m256d xi_27 = _mm256_add_pd(xi_12,xi_26);
            const __m256d xi_28 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_load_pd(& _data_pdfs_21_315_1m1[ctr_0]));
            const __m256d xi_29 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_load_pd(& _data_pdfs_21_36_10[ctr_0]));
            const __m256d xi_30 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_24,xi_28),xi_29),_mm256_load_pd(& _data_pdfs_2m1_311_1m1[ctr_0]));
            const __m256d xi_31 = _mm256_add_pd(_mm256_loadu_pd(& _data_pdfs_21_317_10[ctr_0 + 1]),_mm256_loadu_pd(& _data_pdfs_2m1_313_10[ctr_0 + 1]));
            const __m256d xi_32 = _mm256_add_pd(_mm256_add_pd(xi_31,xi_8),_mm256_loadu_pd(& _data_pdfs_20_37_1m1[ctr_0 + 1]));
            const __m256d xi_33 = _mm256_add_pd(_mm256_add_pd(xi_16,xi_21),_mm256_loadu_pd(& _data_pdfs_20_39_11[ctr_0 + 1]));
            const __m256d xi_34 = _mm256_add_pd(xi_26,_mm256_loadu_pd(& _data_pdfs_21_317_10[ctr_0 + 1]));
            const __m256d xi_35 = _mm256_add_pd(xi_15,_mm256_loadu_pd(& _data_pdfs_2m1_314_10[ctr_0 - 1]));
            const __m256d xi_36 = _mm256_add_pd(xi_34,xi_35);
            const __m256d xi_37 = _mm256_add_pd(_mm256_add_pd(xi_9,_mm256_loadu_pd(& _data_pdfs_20_39_11[ctr_0 + 1])),_mm256_load_pd(& _data_pdfs_2m1_312_11[ctr_0]));
            const __m256d xi_38 = _mm256_add_pd(xi_28,_mm256_load_pd(& _data_pdfs_21_316_11[ctr_0]));
            const __m256d xi_39 = _mm256_add_pd(_mm256_add_pd(xi_22,xi_38),_mm256_load_pd(& _data_pdfs_2m1_311_1m1[ctr_0]));
            const __m256d xi_55 = _mm256_mul_pd(_mm256_set_pd(0.166666666666667,0.166666666666667,0.166666666666667,0.166666666666667),_mm256_load_pd(& _data_force_20_31_10[ctr_0]));
            const __m256d xi_63 = _mm256_mul_pd(_mm256_set_pd(0.166666666666667,0.166666666666667,0.166666666666667,0.166666666666667),_mm256_load_pd(& _data_force_20_30_10[ctr_0]));
            const __m256d xi_67 = _mm256_mul_pd(_mm256_set_pd(0.166666666666667,0.166666666666667,0.166666666666667,0.166666666666667),_mm256_load_pd(& _data_force_20_32_10[ctr_0]));
            const __m256d xi_70 = _mm256_mul_pd(_mm256_set_pd(0.5,0.5,0.5,0.5),_mm256_load_pd(& _data_force_20_31_10[ctr_0]));
            const __m256d xi_74 = _mm256_mul_pd(_mm256_set_pd(0.0833333333333333,0.0833333333333333,0.0833333333333333,0.0833333333333333),_mm256_load_pd(& _data_force_20_30_10[ctr_0]));
            const __m256d xi_78 = _mm256_mul_pd(_mm256_set_pd(0.0833333333333333,0.0833333333333333,0.0833333333333333,0.0833333333333333),_mm256_load_pd(& _data_force_20_31_10[ctr_0]));
            const __m256d xi_88 = _mm256_mul_pd(_mm256_set_pd(0.0833333333333333,0.0833333333333333,0.0833333333333333,0.0833333333333333),_mm256_load_pd(& _data_force_20_32_10[ctr_0]));
            const __m256d xi_99 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_load_pd(& _data_pdfs_20_30_10[ctr_0]));
            const __m256d xi_100 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(3.0,3.0,3.0,3.0),_mm256_load_pd(& _data_pdfs_21_36_10[ctr_0])),_mm256_mul_pd(_mm256_set_pd(3.0,3.0,3.0,3.0),_mm256_load_pd(& _data_pdfs_2m1_35_10[ctr_0]))),xi_99);
            const __m256d xi_101 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-3.0,-3.0,-3.0,-3.0),_mm256_load_pd(& _data_pdfs_21_315_1m1[ctr_0])),_mm256_mul_pd(_mm256_set_pd(-3.0,-3.0,-3.0,-3.0),_mm256_load_pd(& _data_pdfs_21_316_11[ctr_0]))),_mm256_mul_pd(_mm256_set_pd(-3.0,-3.0,-3.0,-3.0),_mm256_load_pd(& _data_pdfs_2m1_311_1m1[ctr_0]))),_mm256_mul_pd(_mm256_set_pd(-3.0,-3.0,-3.0,-3.0),_mm256_load_pd(& _data_pdfs_2m1_312_11[ctr_0]))),_mm256_mul_pd(_mm256_set_pd(3.0,3.0,3.0,3.0),_mm256_load_pd(& _data_pdfs_20_31_1m1[ctr_0]))),_mm256_mul_pd(_mm256_set_pd(3.0,3.0,3.0,3.0),_mm256_load_pd(& _data_pdfs_20_32_11[ctr_0]))),xi_100),_mm256_set_pd(omega_even,omega_even,omega_even,omega_even));
            const __m256d xi_102 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(2.0,2.0,2.0,2.0),_mm256_load_pd(& _data_pdfs_21_315_1m1[ctr_0])),_mm256_mul_pd(_mm256_set_pd(2.0,2.0,2.0,2.0),_mm256_load_pd(& _data_pdfs_21_316_11[ctr_0]))),_mm256_mul_pd(_mm256_set_pd(2.0,2.0,2.0,2.0),_mm256_load_pd(& _data_pdfs_2m1_311_1m1[ctr_0]))),_mm256_mul_pd(_mm256_set_pd(2.0,2.0,2.0,2.0),_mm256_load_pd(& _data_pdfs_2m1_312_11[ctr_0])));
            const __m256d xi_103 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(5.0,5.0,5.0,5.0),_mm256_loadu_pd(& _data_pdfs_20_33_10[ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(5.0,5.0,5.0,5.0),_mm256_loadu_pd(& _data_pdfs_20_34_10[ctr_0 - 1]))),xi_102);
            const __m256d xi_104 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-5.0,-5.0,-5.0,-5.0),_mm256_loadu_pd(& _data_pdfs_21_317_10[ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(-5.0,-5.0,-5.0,-5.0),_mm256_loadu_pd(& _data_pdfs_21_318_10[ctr_0 - 1]))),_mm256_mul_pd(_mm256_set_pd(-5.0,-5.0,-5.0,-5.0),_mm256_loadu_pd(& _data_pdfs_2m1_313_10[ctr_0 + 1]))),_mm256_mul_pd(_mm256_set_pd(-5.0,-5.0,-5.0,-5.0),_mm256_loadu_pd(& _data_pdfs_2m1_314_10[ctr_0 - 1]))),_mm256_mul_pd(_mm256_set_pd(-2.0,-2.0,-2.0,-2.0),_mm256_load_pd(& _data_pdfs_20_31_1m1[ctr_0]))),_mm256_mul_pd(_mm256_set_pd(-2.0,-2.0,-2.0,-2.0),_mm256_load_pd(& _data_pdfs_20_32_11[ctr_0]))),xi_100),xi_103),_mm256_set_pd(omega_even,omega_even,omega_even,omega_even));
            const __m256d xi_107 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_load_pd(& _data_pdfs_2m1_311_1m1[ctr_0]));
            const __m256d xi_108 = _mm256_add_pd(xi_107,xi_22);
            const __m256d xi_109 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_loadu_pd(& _data_pdfs_20_38_1m1[ctr_0 - 1]));
            const __m256d xi_112 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_loadu_pd(& _data_pdfs_2m1_314_10[ctr_0 - 1]));
            const __m256d xi_113 = _mm256_add_pd(_mm256_add_pd(xi_112,xi_19),xi_27);
            const __m256d xi_115 = _mm256_mul_pd(_mm256_set_pd(2.0,2.0,2.0,2.0),_mm256_loadu_pd(& _data_pdfs_2m1_313_10[ctr_0 + 1]));
            const __m256d xi_116 = _mm256_mul_pd(_mm256_set_pd(2.0,2.0,2.0,2.0),_mm256_loadu_pd(& _data_pdfs_2m1_314_10[ctr_0 - 1]));
            const __m256d xi_117 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(2.0,2.0,2.0,2.0),_mm256_loadu_pd(& _data_pdfs_21_317_10[ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(2.0,2.0,2.0,2.0),_mm256_loadu_pd(& _data_pdfs_21_318_10[ctr_0 - 1])));
            const __m256d xi_118 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-7.0,-7.0,-7.0,-7.0),_mm256_loadu_pd(& _data_pdfs_20_310_11[ctr_0 - 1])),_mm256_mul_pd(_mm256_set_pd(-7.0,-7.0,-7.0,-7.0),_mm256_loadu_pd(& _data_pdfs_20_37_1m1[ctr_0 + 1]))),_mm256_mul_pd(_mm256_set_pd(-7.0,-7.0,-7.0,-7.0),_mm256_loadu_pd(& _data_pdfs_20_38_1m1[ctr_0 - 1]))),_mm256_mul_pd(_mm256_set_pd(-7.0,-7.0,-7.0,-7.0),_mm256_loadu_pd(& _data_pdfs_20_39_11[ctr_0 + 1]))),_mm256_mul_pd(_mm256_set_pd(-4.0,-4.0,-4.0,-4.0),_mm256_load_pd(& _data_pdfs_21_36_10[ctr_0]))),_mm256_mul_pd(_mm256_set_pd(-4.0,-4.0,-4.0,-4.0),_mm256_load_pd(& _data_pdfs_2m1_35_10[ctr_0]))),_mm256_mul_pd(_mm256_set_pd(5.0,5.0,5.0,5.0),_mm256_load_pd(& _data_pdfs_20_31_1m1[ctr_0]))),_mm256_mul_pd(_mm256_set_pd(5.0,5.0,5.0,5.0),_mm256_load_pd(& _data_pdfs_20_32_11[ctr_0]))),xi_103),xi_115),xi_116),xi_117),xi_99),_mm256_set_pd(omega_even,omega_even,omega_even,omega_even));
            const __m256d xi_119 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_107,xi_23),xi_38),_mm256_load_pd(& _data_pdfs_20_31_1m1[ctr_0])),_mm256_load_pd(& _data_pdfs_2m1_312_11[ctr_0]));
            const __m256d xi_121 = _mm256_mul_pd(xi_119,_mm256_set_pd(xi_120,xi_120,xi_120,xi_120));
            const __m256d xi_122 = _mm256_mul_pd(_mm256_set_pd(2.0,2.0,2.0,2.0),_mm256_loadu_pd(& _data_pdfs_20_37_1m1[ctr_0 + 1]));
            const __m256d xi_123 = _mm256_mul_pd(_mm256_set_pd(2.0,2.0,2.0,2.0),_mm256_loadu_pd(& _data_pdfs_20_310_11[ctr_0 - 1]));
            const __m256d xi_124 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-2.0,-2.0,-2.0,-2.0),_mm256_loadu_pd(& _data_pdfs_20_38_1m1[ctr_0 - 1])),_mm256_mul_pd(_mm256_set_pd(2.0,2.0,2.0,2.0),_mm256_loadu_pd(& _data_pdfs_20_39_11[ctr_0 + 1])));
            const __m256d xi_125 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_122,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_123),xi_124),xi_25),xi_4);
            const __m256d xi_127 = _mm256_mul_pd(xi_125,_mm256_set_pd(xi_126,xi_126,xi_126,xi_126));
            const __m256d xi_128 = _mm256_mul_pd(xi_127,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_133 = _mm256_mul_pd(xi_118,_mm256_set_pd(-0.0198412698412698,-0.0198412698412698,-0.0198412698412698,-0.0198412698412698));
            const __m256d xi_144 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_112,xi_13),xi_34),_mm256_loadu_pd(& _data_pdfs_20_34_10[ctr_0 - 1])),_mm256_loadu_pd(& _data_pdfs_2m1_313_10[ctr_0 + 1]));
            const __m256d xi_145 = _mm256_mul_pd(xi_144,_mm256_set_pd(xi_120,xi_120,xi_120,xi_120));
            const __m256d xi_146 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_123,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_122),xi_124),xi_14),xi_35),_mm256_loadu_pd(& _data_pdfs_20_34_10[ctr_0 - 1])),_mm256_loadu_pd(& _data_pdfs_21_318_10[ctr_0 - 1]));
            const __m256d xi_147 = _mm256_mul_pd(xi_146,_mm256_set_pd(xi_126,xi_126,xi_126,xi_126));
            const __m256d xi_149 = _mm256_mul_pd(xi_147,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_150 = _mm256_add_pd(_mm256_load_pd(& _data_pdfs_21_315_1m1[ctr_0]),_mm256_load_pd(& _data_pdfs_21_316_11[ctr_0]));
            const __m256d xi_151 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_108,xi_150),xi_29),_mm256_load_pd(& _data_pdfs_2m1_35_10[ctr_0]));
            const __m256d xi_152 = _mm256_mul_pd(xi_151,_mm256_set_pd(xi_120,xi_120,xi_120,xi_120));
            const __m256d xi_153 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_115,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_116,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_117),xi_30),xi_6);
            const __m256d xi_154 = _mm256_mul_pd(xi_153,_mm256_set_pd(xi_126,xi_126,xi_126,xi_126));
            const __m256d xi_155 = _mm256_mul_pd(xi_154,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_162 = _mm256_mul_pd(xi_118,_mm256_set_pd(0.0138888888888889,0.0138888888888889,0.0138888888888889,0.0138888888888889));
            const __m256d xi_177 = _mm256_mul_pd(xi_104,_mm256_set_pd(-0.00714285714285714,-0.00714285714285714,-0.00714285714285714,-0.00714285714285714));
            const __m256d xi_179 = _mm256_mul_pd(xi_101,_mm256_set_pd(0.025,0.025,0.025,0.025));
            const __m256d xi_182 = _mm256_mul_pd(xi_153,_mm256_set_pd(xi_181,xi_181,xi_181,xi_181));
            const __m256d xi_184 = _mm256_mul_pd(xi_151,_mm256_set_pd(xi_183,xi_183,xi_183,xi_183));
            const __m256d xi_185 = _mm256_mul_pd(xi_118,_mm256_set_pd(-0.00396825396825397,-0.00396825396825397,-0.00396825396825397,-0.00396825396825397));
            const __m256d xi_189 = _mm256_mul_pd(xi_125,_mm256_set_pd(xi_181,xi_181,xi_181,xi_181));
            const __m256d xi_190 = _mm256_mul_pd(xi_119,_mm256_set_pd(xi_183,xi_183,xi_183,xi_183));
            const __m256d xi_196 = _mm256_mul_pd(xi_104,_mm256_set_pd(0.0178571428571429,0.0178571428571429,0.0178571428571429,0.0178571428571429));
            const __m256d xi_199 = _mm256_mul_pd(xi_144,_mm256_set_pd(xi_183,xi_183,xi_183,xi_183));
            const __m256d xi_200 = _mm256_mul_pd(xi_146,_mm256_set_pd(xi_181,xi_181,xi_181,xi_181));
            const __m256d vel0Term = xi_2;
            const __m256d vel1Term = xi_5;
            const __m256d vel2Term = xi_7;
            const __m256d rho = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(vel0Term,vel1Term),vel2Term),xi_10),xi_8),xi_9),_mm256_load_pd(& _data_pdfs_20_30_10[ctr_0]));
            const __m256d xi_11 = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),rho);
            const __m256d u_0 = _mm256_mul_pd(xi_11,_mm256_add_pd(_mm256_add_pd(vel0Term,xi_14),xi_19));
            const __m256d xi_40 = _mm256_mul_pd(u_0,_mm256_load_pd(& _data_force_20_30_10[ctr_0]));
            const __m256d xi_41 = _mm256_mul_pd(xi_40,_mm256_set_pd(0.333333333333333,0.333333333333333,0.333333333333333,0.333333333333333));
            const __m256d xi_47 = _mm256_mul_pd(xi_41,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_105 = _mm256_mul_pd(rho,(_mm256_mul_pd(u_0,u_0)));
            const __m256d xi_141 = _mm256_mul_pd(rho,u_0);
            const __m256d xi_142 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(vel0Term,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_141),xi_32);
            const __m256d xi_143 = _mm256_mul_pd(xi_142,_mm256_set_pd(xi_131,xi_131,xi_131,xi_131));
            const __m256d xi_167 = _mm256_mul_pd(xi_142,_mm256_set_pd(xi_164,xi_164,xi_164,xi_164));
            const __m256d u_1 = _mm256_mul_pd(xi_11,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(vel1Term,xi_17),xi_21),xi_25));
            const __m256d xi_42 = _mm256_mul_pd(u_1,_mm256_load_pd(& _data_force_20_31_10[ctr_0]));
            const __m256d xi_43 = _mm256_mul_pd(xi_42,_mm256_set_pd(0.333333333333333,0.333333333333333,0.333333333333333,0.333333333333333));
            const __m256d xi_48 = _mm256_mul_pd(xi_43,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_69 = _mm256_mul_pd(u_1,_mm256_set_pd(0.5,0.5,0.5,0.5));
            const __m256d xi_72 = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(u_0,xi_70),_mm256_mul_pd(xi_69,_mm256_load_pd(& _data_force_20_30_10[ctr_0]))),_mm256_set_pd(xi_71,xi_71,xi_71,xi_71));
            const __m256d xi_73 = _mm256_mul_pd(xi_72,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_110 = _mm256_mul_pd(rho,(_mm256_mul_pd(u_1,u_1)));
            const __m256d xi_111 = _mm256_add_pd(_mm256_add_pd(xi_109,xi_110),xi_20);
            const __m256d xi_129 = _mm256_mul_pd(rho,u_1);
            const __m256d xi_130 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(vel1Term,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_109),xi_129),xi_37),_mm256_loadu_pd(& _data_pdfs_20_310_11[ctr_0 - 1]));
            const __m256d xi_132 = _mm256_mul_pd(xi_130,_mm256_set_pd(xi_131,xi_131,xi_131,xi_131));
            const __m256d xi_165 = _mm256_mul_pd(xi_130,_mm256_set_pd(xi_164,xi_164,xi_164,xi_164));
            const __m256d xi_166 = _mm256_add_pd(xi_127,xi_165);
            const __m256d xi_175 = _mm256_mul_pd(xi_165,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_176 = _mm256_add_pd(xi_128,xi_175);
            const __m256d xi_191 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_189,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_165),xi_190);
            const __m256d xi_192 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_190,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_175),xi_189);
            const __m256d u_2 = _mm256_mul_pd(xi_11,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(vel2Term,xi_27),xi_30),_mm256_loadu_pd(& _data_pdfs_2m1_314_10[ctr_0 - 1])));
            const __m256d xi_44 = _mm256_mul_pd(u_2,_mm256_load_pd(& _data_force_20_32_10[ctr_0]));
            const __m256d xi_45 = _mm256_mul_pd(xi_44,_mm256_set_pd(0.333333333333333,0.333333333333333,0.333333333333333,0.333333333333333));
            const __m256d xi_46 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(xi_41,xi_43),xi_45),_mm256_set_pd(-omega_bulk + 2.0,-omega_bulk + 2.0,-omega_bulk + 2.0,-omega_bulk + 2.0));
            const __m256d xi_49 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_44,_mm256_set_pd(0.666666666666667,0.666666666666667,0.666666666666667,0.666666666666667)),xi_47),xi_48);
            const __m256d xi_52 = _mm256_mul_pd(xi_45,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_53 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_42,_mm256_set_pd(0.666666666666667,0.666666666666667,0.666666666666667,0.666666666666667)),xi_47),xi_52);
            const __m256d xi_54 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_40,_mm256_set_pd(0.666666666666667,0.666666666666667,0.666666666666667,0.666666666666667)),xi_48),xi_52);
            const __m256d xi_57 = _mm256_mul_pd(xi_49,_mm256_set_pd(xi_56,xi_56,xi_56,xi_56));
            const __m256d xi_58 = _mm256_mul_pd(xi_57,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_59 = _mm256_mul_pd(xi_54,_mm256_set_pd(xi_56,xi_56,xi_56,xi_56));
            const __m256d xi_60 = _mm256_mul_pd(xi_59,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_62 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_53,_mm256_set_pd(xi_61,xi_61,xi_61,xi_61)),xi_58),xi_60);
            const __m256d xi_64 = _mm256_mul_pd(xi_53,_mm256_set_pd(xi_56,xi_56,xi_56,xi_56));
            const __m256d xi_65 = _mm256_mul_pd(xi_64,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_66 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_54,_mm256_set_pd(xi_61,xi_61,xi_61,xi_61)),xi_58),xi_65);
            const __m256d xi_68 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_49,_mm256_set_pd(xi_61,xi_61,xi_61,xi_61)),xi_60),xi_65);
            const __m256d xi_75 = _mm256_add_pd(_mm256_mul_pd(xi_74,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_59);
            const __m256d xi_77 = _mm256_mul_pd(_mm256_mul_pd(xi_49,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_76,xi_76,xi_76,xi_76));
            const __m256d xi_79 = _mm256_mul_pd(xi_46,_mm256_set_pd(0.125,0.125,0.125,0.125));
            const __m256d xi_80 = _mm256_add_pd(xi_64,xi_79);
            const __m256d xi_81 = _mm256_add_pd(xi_78,xi_80);
            const __m256d xi_82 = _mm256_add_pd(xi_77,xi_81);
            const __m256d xi_83 = _mm256_add_pd(xi_59,xi_74);
            const __m256d xi_84 = _mm256_add_pd(_mm256_mul_pd(xi_78,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_80);
            const __m256d xi_85 = _mm256_add_pd(xi_77,xi_84);
            const __m256d xi_86 = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(u_2,xi_70),_mm256_mul_pd(xi_69,_mm256_load_pd(& _data_force_20_32_10[ctr_0]))),_mm256_set_pd(xi_71,xi_71,xi_71,xi_71));
            const __m256d xi_87 = _mm256_mul_pd(_mm256_mul_pd(xi_54,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_76,xi_76,xi_76,xi_76));
            const __m256d xi_89 = _mm256_add_pd(xi_57,xi_88);
            const __m256d xi_90 = _mm256_add_pd(xi_87,xi_89);
            const __m256d xi_91 = _mm256_mul_pd(xi_86,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_92 = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(u_0,_mm256_set_pd(0.5,0.5,0.5,0.5)),_mm256_load_pd(& _data_force_20_32_10[ctr_0])),_mm256_mul_pd(_mm256_mul_pd(u_2,_mm256_set_pd(0.5,0.5,0.5,0.5)),_mm256_load_pd(& _data_force_20_30_10[ctr_0]))),_mm256_set_pd(xi_71,xi_71,xi_71,xi_71));
            const __m256d xi_93 = _mm256_mul_pd(xi_92,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_94 = _mm256_mul_pd(_mm256_mul_pd(xi_53,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_76,xi_76,xi_76,xi_76));
            const __m256d xi_95 = _mm256_add_pd(_mm256_add_pd(xi_79,xi_89),xi_94);
            const __m256d xi_96 = _mm256_add_pd(_mm256_mul_pd(xi_88,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_57);
            const __m256d xi_97 = _mm256_add_pd(xi_87,xi_96);
            const __m256d xi_98 = _mm256_add_pd(_mm256_add_pd(xi_79,xi_94),xi_96);
            const __m256d xi_106 = _mm256_mul_pd(rho,(_mm256_mul_pd(u_2,u_2)));
            const __m256d xi_114 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_105,xi_106),xi_108),xi_111),xi_113),xi_24),xi_28),_mm256_load_pd(& _data_pdfs_20_30_10[ctr_0])),_mm256_set_pd(omega_bulk,omega_bulk,omega_bulk,omega_bulk));
            const __m256d xi_134 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_106,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_load_pd(& _data_pdfs_21_36_10[ctr_0])),_mm256_load_pd(& _data_pdfs_2m1_35_10[ctr_0]));
            const __m256d xi_135 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_load_pd(& _data_pdfs_20_31_1m1[ctr_0])),xi_1),xi_111),xi_134),xi_18),xi_23),xi_31),_mm256_set_pd(omega_shear,omega_shear,omega_shear,omega_shear));
            const __m256d xi_136 = _mm256_mul_pd(xi_135,_mm256_set_pd(0.125,0.125,0.125,0.125));
            const __m256d xi_137 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_105,_mm256_set_pd(2.0,2.0,2.0,2.0)),_mm256_mul_pd(xi_110,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_set_pd(-2.0,-2.0,-2.0,-2.0),_mm256_loadu_pd(& _data_pdfs_20_33_10[ctr_0 + 1]))),_mm256_mul_pd(_mm256_set_pd(-2.0,-2.0,-2.0,-2.0),_mm256_loadu_pd(& _data_pdfs_20_34_10[ctr_0 - 1]))),xi_102),xi_109),xi_113),xi_134),xi_20),_mm256_load_pd(& _data_pdfs_20_31_1m1[ctr_0])),_mm256_load_pd(& _data_pdfs_20_32_11[ctr_0])),_mm256_set_pd(omega_shear,omega_shear,omega_shear,omega_shear));
            const __m256d xi_138 = _mm256_mul_pd(xi_137,_mm256_set_pd(-0.0416666666666667,-0.0416666666666667,-0.0416666666666667,-0.0416666666666667));
            const __m256d xi_139 = _mm256_add_pd(_mm256_mul_pd(xi_101,_mm256_set_pd(-0.05,-0.05,-0.05,-0.05)),xi_138);
            const __m256d xi_140 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_104,_mm256_set_pd(0.0142857142857143,0.0142857142857143,0.0142857142857143,0.0142857142857143)),xi_133),xi_136),xi_139);
            const __m256d xi_148 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_104,_mm256_set_pd(-0.0357142857142857,-0.0357142857142857,-0.0357142857142857,-0.0357142857142857)),_mm256_mul_pd(xi_137,_mm256_set_pd(0.0833333333333333,0.0833333333333333,0.0833333333333333,0.0833333333333333))),xi_133);
            const __m256d xi_156 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(rho,u_2),_mm256_mul_pd(vel2Term,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_10),xi_107),xi_112),xi_150),_mm256_loadu_pd(& _data_pdfs_21_318_10[ctr_0 - 1]));
            const __m256d xi_157 = _mm256_mul_pd(xi_156,_mm256_set_pd(xi_131,xi_131,xi_131,xi_131));
            const __m256d xi_158 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_104,_mm256_set_pd(-0.0214285714285714,-0.0214285714285714,-0.0214285714285714,-0.0214285714285714)),_mm256_mul_pd(xi_118,_mm256_set_pd(0.0158730158730159,0.0158730158730159,0.0158730158730159,0.0158730158730159))),_mm256_mul_pd(xi_136,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_139);
            const __m256d xi_161 = _mm256_mul_pd(xi_135,_mm256_set_pd(0.0625,0.0625,0.0625,0.0625));
            const __m256d xi_168 = _mm256_mul_pd(xi_114,_mm256_set_pd(0.0416666666666667,0.0416666666666667,0.0416666666666667,0.0416666666666667));
            const __m256d xi_169 = _mm256_add_pd(_mm256_mul_pd(xi_137,_mm256_set_pd(0.0208333333333333,0.0208333333333333,0.0208333333333333,0.0208333333333333)),xi_168);
            const __m256d xi_170 = _mm256_add_pd(_mm256_mul_pd(xi_167,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_169);
            const __m256d xi_171 = _mm256_add_pd(xi_149,xi_170);
            const __m256d xi_173 = _mm256_add_pd(xi_167,xi_169);
            const __m256d xi_174 = _mm256_add_pd(xi_147,xi_173);
            const __m256d xi_186 = _mm256_mul_pd(xi_156,_mm256_set_pd(xi_164,xi_164,xi_164,xi_164));
            const __m256d xi_187 = _mm256_add_pd(xi_185,xi_186);
            const __m256d xi_188 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_182,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_184),xi_187);
            const __m256d xi_195 = _mm256_mul_pd(xi_161,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_198 = _mm256_add_pd(xi_154,xi_187);
            const __m256d xi_201 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_199,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_170),xi_200);
            const __m256d xi_203 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_200,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_173),xi_199);
            const __m256d xi_204 = _mm256_add_pd(_mm256_mul_pd(xi_186,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_185);
            const __m256d xi_205 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_184,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_182),xi_204);
            const __m256d xi_206 = _mm256_add_pd(xi_155,xi_204);
            const __m256d p_0 = _mm256_add_pd(xi_2,xi_32);
            const __m256d p_1 = xi_33;
            const __m256d xi_160 = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(p_1,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_0,xi_129)),_mm256_set_pd(xi_159,xi_159,xi_159,xi_159));
            const __m256d xi_163 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_160,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_161),xi_162);
            const __m256d xi_172 = _mm256_add_pd(_mm256_add_pd(xi_160,xi_161),xi_162);
            const __m256d p_2 = xi_36;
            const __m256d xi_194 = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(p_2,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_2,xi_141)),_mm256_set_pd(xi_159,xi_159,xi_159,xi_159));
            const __m256d xi_197 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_194,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_195),xi_196);
            const __m256d xi_202 = _mm256_add_pd(_mm256_add_pd(xi_194,xi_195),xi_196);
            const __m256d p_3 = xi_33;
            const __m256d p_4 = _mm256_add_pd(_mm256_add_pd(xi_0,xi_37),xi_5);
            const __m256d p_5 = xi_39;
            const __m256d xi_178 = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(p_5,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_2,xi_129)),_mm256_set_pd(xi_159,xi_159,xi_159,xi_159));
            const __m256d xi_180 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_138,xi_168),xi_177),xi_178),xi_179);
            const __m256d xi_193 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_178,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_138),xi_168),xi_177),xi_179);
            const __m256d p_6 = xi_36;
            const __m256d p_7 = xi_39;
            const __m256d p_8 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_1,xi_10),xi_3),xi_7),_mm256_load_pd(& _data_pdfs_21_316_11[ctr_0]));
            const __m256d forceTerm_0 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_46,_mm256_set_pd(-1.5,-1.5,-1.5,-1.5)),_mm256_mul_pd(_mm256_mul_pd(xi_49,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_51,xi_51,xi_51,xi_51))),_mm256_mul_pd(_mm256_mul_pd(xi_53,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_51,xi_51,xi_51,xi_51))),_mm256_mul_pd(_mm256_mul_pd(xi_54,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_51,xi_51,xi_51,xi_51)));
            const __m256d forceTerm_1 = _mm256_add_pd(xi_55,xi_62);
            const __m256d forceTerm_2 = _mm256_add_pd(_mm256_mul_pd(xi_55,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_62);
            const __m256d forceTerm_3 = _mm256_add_pd(_mm256_mul_pd(xi_63,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_66);
            const __m256d forceTerm_4 = _mm256_add_pd(xi_63,xi_66);
            const __m256d forceTerm_5 = _mm256_add_pd(xi_67,xi_68);
            const __m256d forceTerm_6 = _mm256_add_pd(_mm256_mul_pd(xi_67,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_68);
            const __m256d forceTerm_7 = _mm256_add_pd(_mm256_add_pd(xi_73,xi_75),xi_82);
            const __m256d forceTerm_8 = _mm256_add_pd(_mm256_add_pd(xi_72,xi_82),xi_83);
            const __m256d forceTerm_9 = _mm256_add_pd(_mm256_add_pd(xi_72,xi_75),xi_85);
            const __m256d forceTerm_10 = _mm256_add_pd(_mm256_add_pd(xi_73,xi_83),xi_85);
            const __m256d forceTerm_11 = _mm256_add_pd(_mm256_add_pd(xi_81,xi_86),xi_90);
            const __m256d forceTerm_12 = _mm256_add_pd(_mm256_add_pd(xi_84,xi_90),xi_91);
            const __m256d forceTerm_13 = _mm256_add_pd(_mm256_add_pd(xi_75,xi_93),xi_95);
            const __m256d forceTerm_14 = _mm256_add_pd(_mm256_add_pd(xi_83,xi_92),xi_95);
            const __m256d forceTerm_15 = _mm256_add_pd(_mm256_add_pd(xi_81,xi_91),xi_97);
            const __m256d forceTerm_16 = _mm256_add_pd(_mm256_add_pd(xi_84,xi_86),xi_97);
            const __m256d forceTerm_17 = _mm256_add_pd(_mm256_add_pd(xi_75,xi_92),xi_98);
            const __m256d forceTerm_18 = _mm256_add_pd(_mm256_add_pd(xi_83,xi_93),xi_98);
            _mm256_store_pd(&_data_pdfs_tmp_20_30_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_101,_mm256_set_pd(0.1,0.1,0.1,0.1)),_mm256_mul_pd(xi_104,_mm256_set_pd(0.0428571428571429,0.0428571428571429,0.0428571428571429,0.0428571428571429))),_mm256_mul_pd(xi_114,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5))),_mm256_mul_pd(xi_118,_mm256_set_pd(0.0238095238095238,0.0238095238095238,0.0238095238095238,0.0238095238095238))),forceTerm_0),_mm256_load_pd(& _data_pdfs_20_30_10[ctr_0])));
            _mm256_store_pd(&_data_pdfs_tmp_20_31_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_121,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),forceTerm_1),xi_128),xi_132),xi_140),_mm256_load_pd(& _data_pdfs_20_31_1m1[ctr_0])));
            _mm256_store_pd(&_data_pdfs_tmp_20_32_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_132,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),forceTerm_2),xi_121),xi_127),xi_140),_mm256_load_pd(& _data_pdfs_20_32_11[ctr_0])));
            _mm256_store_pd(&_data_pdfs_tmp_20_33_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_143,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),forceTerm_3),xi_145),xi_147),xi_148),_mm256_loadu_pd(& _data_pdfs_20_33_10[ctr_0 + 1])));
            _mm256_store_pd(&_data_pdfs_tmp_20_34_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_145,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),forceTerm_4),xi_143),xi_148),xi_149),_mm256_loadu_pd(& _data_pdfs_20_34_10[ctr_0 - 1])));
            _mm256_store_pd(&_data_pdfs_tmp_20_35_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_152,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),forceTerm_5),xi_155),xi_157),xi_158),_mm256_load_pd(& _data_pdfs_2m1_35_10[ctr_0])));
            _mm256_store_pd(&_data_pdfs_tmp_20_36_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_157,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),forceTerm_6),xi_152),xi_154),xi_158),_mm256_load_pd(& _data_pdfs_21_36_10[ctr_0])));
            _mm256_store_pd(&_data_pdfs_tmp_20_37_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_7,xi_163),xi_166),xi_171),_mm256_loadu_pd(& _data_pdfs_20_37_1m1[ctr_0 + 1])));
            _mm256_store_pd(&_data_pdfs_tmp_20_38_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_8,xi_166),xi_172),xi_174),_mm256_loadu_pd(& _data_pdfs_20_38_1m1[ctr_0 - 1])));
            _mm256_store_pd(&_data_pdfs_tmp_20_39_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_9,xi_171),xi_172),xi_176),_mm256_loadu_pd(& _data_pdfs_20_39_11[ctr_0 + 1])));
            _mm256_store_pd(&_data_pdfs_tmp_20_310_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_10,xi_163),xi_174),xi_176),_mm256_loadu_pd(& _data_pdfs_20_310_11[ctr_0 - 1])));
            _mm256_store_pd(&_data_pdfs_tmp_20_311_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_11,xi_180),xi_188),xi_191),_mm256_load_pd(& _data_pdfs_2m1_311_1m1[ctr_0])));
            _mm256_store_pd(&_data_pdfs_tmp_20_312_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_12,xi_188),xi_192),xi_193),_mm256_load_pd(& _data_pdfs_2m1_312_11[ctr_0])));
            _mm256_store_pd(&_data_pdfs_tmp_20_313_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_13,xi_197),xi_198),xi_201),_mm256_loadu_pd(& _data_pdfs_2m1_313_10[ctr_0 + 1])));
            _mm256_store_pd(&_data_pdfs_tmp_20_314_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_14,xi_198),xi_202),xi_203),_mm256_loadu_pd(& _data_pdfs_2m1_314_10[ctr_0 - 1])));
            _mm256_store_pd(&_data_pdfs_tmp_20_315_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_15,xi_191),xi_193),xi_205),_mm256_load_pd(& _data_pdfs_21_315_1m1[ctr_0])));
            _mm256_store_pd(&_data_pdfs_tmp_20_316_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_16,xi_180),xi_192),xi_205),_mm256_load_pd(& _data_pdfs_21_316_11[ctr_0])));
            _mm256_store_pd(&_data_pdfs_tmp_20_317_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_17,xi_201),xi_202),xi_206),_mm256_loadu_pd(& _data_pdfs_21_317_10[ctr_0 + 1])));
            _mm256_store_pd(&_data_pdfs_tmp_20_318_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_18,xi_197),xi_203),xi_206),_mm256_loadu_pd(& _data_pdfs_21_318_10[ctr_0 - 1])));
         }
      }
   }
}
}
namespace internal_kernel_collide {
static FUNC_PREFIX void kernel_collide(double * RESTRICT const _data_force, double * RESTRICT _data_pdfs, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, double omega_bulk, double omega_even, double omega_odd, double omega_shear)
{
   const double xi_50 = -omega_shear + 2.0;
   const double xi_51 = xi_50*0.5;
   const double xi_56 = xi_50*0.0833333333333333;
   const double xi_61 = xi_50*0.166666666666667;
   const double xi_71 = xi_50*0.25;
   const double xi_76 = xi_50*0.0416666666666667;
   const double xi_120 = omega_odd*0.25;
   const double xi_126 = omega_odd*0.0833333333333333;
   const double xi_159 = omega_shear*0.25;
   const double xi_181 = omega_odd*0.0416666666666667;
   const double xi_183 = omega_odd*0.125;
   const int64_t rr_0 = 0.0;
   const double xi_131 = rr_0*0.166666666666667;
   const double xi_164 = rr_0*0.0833333333333333;
   for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1)
   {
      double * RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_35 = _data_pdfs + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * RESTRICT _data_pdfs_20_311 = _data_pdfs + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      double * RESTRICT _data_force_20_31 = _data_force + _stride_force_2*ctr_2 + _stride_force_3;
      double * RESTRICT _data_pdfs_20_314 = _data_pdfs + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      double * RESTRICT _data_force_20_32 = _data_force + _stride_force_2*ctr_2 + 2*_stride_force_3;
      double * RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_313 = _data_pdfs + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      double * RESTRICT _data_force_20_30 = _data_force + _stride_force_2*ctr_2;
      double * RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_312 = _data_pdfs + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_20_38_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_38;
         double * RESTRICT _data_pdfs_20_31_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_31;
         double * RESTRICT _data_pdfs_20_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_35;
         double * RESTRICT _data_pdfs_20_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_318;
         double * RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * RESTRICT _data_pdfs_20_311_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_311;
         double * RESTRICT _data_pdfs_20_310_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_310;
         double * RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * RESTRICT _data_pdfs_20_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_36;
         double * RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * RESTRICT _data_pdfs_20_37_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_37;
         double * RESTRICT _data_pdfs_20_315_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_315;
         double * RESTRICT _data_force_20_31_10 = _stride_force_1*ctr_1 + _data_force_20_31;
         double * RESTRICT _data_pdfs_20_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_314;
         double * RESTRICT _data_force_20_32_10 = _stride_force_1*ctr_1 + _data_force_20_32;
         double * RESTRICT _data_pdfs_20_39_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_39;
         double * RESTRICT _data_pdfs_20_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_313;
         double * RESTRICT _data_pdfs_20_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_317;
         double * RESTRICT _data_force_20_30_10 = _stride_force_1*ctr_1 + _data_force_20_30;
         double * RESTRICT _data_pdfs_20_32_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_32;
         double * RESTRICT _data_pdfs_20_312_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_312;
         double * RESTRICT _data_pdfs_20_316_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_316;
         for (int64_t ctr_0 = 0; ctr_0 < ((_size_force_0) % (4) == 0 ? _size_force_0 : ((int64_t)((_size_force_0) / (4))+1) * (4)); ctr_0 += 4)
         {
            const __m256d xi_207 = _mm256_load_pd(& _data_pdfs_20_38_10[ctr_0]);
            const __m256d xi_208 = _mm256_load_pd(& _data_pdfs_20_31_10[ctr_0]);
            const __m256d xi_209 = _mm256_load_pd(& _data_pdfs_20_35_10[ctr_0]);
            const __m256d xi_210 = _mm256_load_pd(& _data_pdfs_20_318_10[ctr_0]);
            const __m256d xi_211 = _mm256_load_pd(& _data_pdfs_20_30_10[ctr_0]);
            const __m256d xi_212 = _mm256_load_pd(& _data_pdfs_20_311_10[ctr_0]);
            const __m256d xi_213 = _mm256_load_pd(& _data_pdfs_20_310_10[ctr_0]);
            const __m256d xi_214 = _mm256_load_pd(& _data_pdfs_20_34_10[ctr_0]);
            const __m256d xi_215 = _mm256_load_pd(& _data_pdfs_20_36_10[ctr_0]);
            const __m256d xi_216 = _mm256_load_pd(& _data_pdfs_20_33_10[ctr_0]);
            const __m256d xi_217 = _mm256_load_pd(& _data_pdfs_20_37_10[ctr_0]);
            const __m256d xi_218 = _mm256_load_pd(& _data_pdfs_20_315_10[ctr_0]);
            const __m256d xi_219 = _mm256_load_pd(& _data_force_20_31_10[ctr_0]);
            const __m256d xi_220 = _mm256_load_pd(& _data_pdfs_20_314_10[ctr_0]);
            const __m256d xi_221 = _mm256_load_pd(& _data_force_20_32_10[ctr_0]);
            const __m256d xi_222 = _mm256_load_pd(& _data_pdfs_20_39_10[ctr_0]);
            const __m256d xi_223 = _mm256_load_pd(& _data_pdfs_20_313_10[ctr_0]);
            const __m256d xi_224 = _mm256_load_pd(& _data_pdfs_20_317_10[ctr_0]);
            const __m256d xi_225 = _mm256_load_pd(& _data_force_20_30_10[ctr_0]);
            const __m256d xi_226 = _mm256_load_pd(& _data_pdfs_20_32_10[ctr_0]);
            const __m256d xi_227 = _mm256_load_pd(& _data_pdfs_20_312_10[ctr_0]);
            const __m256d xi_228 = _mm256_load_pd(& _data_pdfs_20_316_10[ctr_0]);
            const __m256d xi_0 = _mm256_add_pd(xi_207,xi_213);
            const __m256d xi_1 = _mm256_add_pd(xi_210,xi_220);
            const __m256d xi_2 = _mm256_add_pd(_mm256_add_pd(xi_0,xi_1),xi_214);
            const __m256d xi_3 = _mm256_add_pd(xi_212,xi_218);
            const __m256d xi_4 = _mm256_add_pd(xi_208,xi_3);
            const __m256d xi_5 = _mm256_add_pd(xi_217,xi_4);
            const __m256d xi_6 = _mm256_add_pd(xi_209,xi_227);
            const __m256d xi_7 = _mm256_add_pd(xi_223,xi_6);
            const __m256d xi_8 = _mm256_add_pd(xi_216,xi_222);
            const __m256d xi_9 = _mm256_add_pd(xi_226,xi_228);
            const __m256d xi_10 = _mm256_add_pd(xi_215,xi_224);
            const __m256d xi_12 = _mm256_mul_pd(xi_224,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_13 = _mm256_mul_pd(xi_216,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_14 = _mm256_add_pd(xi_12,xi_13);
            const __m256d xi_15 = _mm256_mul_pd(xi_223,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_16 = _mm256_mul_pd(xi_217,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_17 = _mm256_mul_pd(xi_222,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_18 = _mm256_add_pd(xi_16,xi_17);
            const __m256d xi_19 = _mm256_add_pd(xi_15,xi_18);
            const __m256d xi_20 = _mm256_mul_pd(xi_213,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_21 = _mm256_add_pd(xi_20,xi_207);
            const __m256d xi_22 = _mm256_mul_pd(xi_227,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_23 = _mm256_mul_pd(xi_226,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_24 = _mm256_mul_pd(xi_228,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_25 = _mm256_add_pd(_mm256_add_pd(xi_22,xi_23),xi_24);
            const __m256d xi_26 = _mm256_mul_pd(xi_210,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_27 = _mm256_add_pd(xi_12,xi_26);
            const __m256d xi_28 = _mm256_mul_pd(xi_218,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_29 = _mm256_mul_pd(xi_215,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_30 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_212,xi_24),xi_28),xi_29);
            const __m256d xi_31 = _mm256_add_pd(xi_223,xi_224);
            const __m256d xi_32 = _mm256_add_pd(_mm256_add_pd(xi_217,xi_31),xi_8);
            const __m256d xi_33 = _mm256_add_pd(_mm256_add_pd(xi_16,xi_21),xi_222);
            const __m256d xi_34 = _mm256_add_pd(xi_224,xi_26);
            const __m256d xi_35 = _mm256_add_pd(xi_15,xi_220);
            const __m256d xi_36 = _mm256_add_pd(xi_34,xi_35);
            const __m256d xi_37 = _mm256_add_pd(_mm256_add_pd(xi_222,xi_227),xi_9);
            const __m256d xi_38 = _mm256_add_pd(xi_228,xi_28);
            const __m256d xi_39 = _mm256_add_pd(_mm256_add_pd(xi_212,xi_22),xi_38);
            const __m256d xi_55 = _mm256_mul_pd(xi_219,_mm256_set_pd(0.166666666666667,0.166666666666667,0.166666666666667,0.166666666666667));
            const __m256d xi_63 = _mm256_mul_pd(xi_225,_mm256_set_pd(0.166666666666667,0.166666666666667,0.166666666666667,0.166666666666667));
            const __m256d xi_67 = _mm256_mul_pd(xi_221,_mm256_set_pd(0.166666666666667,0.166666666666667,0.166666666666667,0.166666666666667));
            const __m256d xi_70 = _mm256_mul_pd(xi_219,_mm256_set_pd(0.5,0.5,0.5,0.5));
            const __m256d xi_74 = _mm256_mul_pd(xi_225,_mm256_set_pd(0.0833333333333333,0.0833333333333333,0.0833333333333333,0.0833333333333333));
            const __m256d xi_78 = _mm256_mul_pd(xi_219,_mm256_set_pd(0.0833333333333333,0.0833333333333333,0.0833333333333333,0.0833333333333333));
            const __m256d xi_88 = _mm256_mul_pd(xi_221,_mm256_set_pd(0.0833333333333333,0.0833333333333333,0.0833333333333333,0.0833333333333333));
            const __m256d xi_99 = _mm256_mul_pd(xi_211,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_100 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_209,_mm256_set_pd(3.0,3.0,3.0,3.0)),_mm256_mul_pd(xi_215,_mm256_set_pd(3.0,3.0,3.0,3.0))),xi_99);
            const __m256d xi_101 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_208,_mm256_set_pd(3.0,3.0,3.0,3.0)),_mm256_mul_pd(xi_212,_mm256_set_pd(-3.0,-3.0,-3.0,-3.0))),_mm256_mul_pd(xi_218,_mm256_set_pd(-3.0,-3.0,-3.0,-3.0))),_mm256_mul_pd(xi_226,_mm256_set_pd(3.0,3.0,3.0,3.0))),_mm256_mul_pd(xi_227,_mm256_set_pd(-3.0,-3.0,-3.0,-3.0))),_mm256_mul_pd(xi_228,_mm256_set_pd(-3.0,-3.0,-3.0,-3.0))),xi_100),_mm256_set_pd(omega_even,omega_even,omega_even,omega_even));
            const __m256d xi_102 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_212,_mm256_set_pd(2.0,2.0,2.0,2.0)),_mm256_mul_pd(xi_218,_mm256_set_pd(2.0,2.0,2.0,2.0))),_mm256_mul_pd(xi_227,_mm256_set_pd(2.0,2.0,2.0,2.0))),_mm256_mul_pd(xi_228,_mm256_set_pd(2.0,2.0,2.0,2.0)));
            const __m256d xi_103 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_214,_mm256_set_pd(5.0,5.0,5.0,5.0)),_mm256_mul_pd(xi_216,_mm256_set_pd(5.0,5.0,5.0,5.0))),xi_102);
            const __m256d xi_104 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_208,_mm256_set_pd(-2.0,-2.0,-2.0,-2.0)),_mm256_mul_pd(xi_210,_mm256_set_pd(-5.0,-5.0,-5.0,-5.0))),_mm256_mul_pd(xi_220,_mm256_set_pd(-5.0,-5.0,-5.0,-5.0))),_mm256_mul_pd(xi_223,_mm256_set_pd(-5.0,-5.0,-5.0,-5.0))),_mm256_mul_pd(xi_224,_mm256_set_pd(-5.0,-5.0,-5.0,-5.0))),_mm256_mul_pd(xi_226,_mm256_set_pd(-2.0,-2.0,-2.0,-2.0))),xi_100),xi_103),_mm256_set_pd(omega_even,omega_even,omega_even,omega_even));
            const __m256d xi_107 = _mm256_mul_pd(xi_212,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_108 = _mm256_add_pd(xi_107,xi_22);
            const __m256d xi_109 = _mm256_mul_pd(xi_207,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_112 = _mm256_mul_pd(xi_220,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_113 = _mm256_add_pd(_mm256_add_pd(xi_112,xi_19),xi_27);
            const __m256d xi_115 = _mm256_mul_pd(xi_223,_mm256_set_pd(2.0,2.0,2.0,2.0));
            const __m256d xi_116 = _mm256_mul_pd(xi_220,_mm256_set_pd(2.0,2.0,2.0,2.0));
            const __m256d xi_117 = _mm256_add_pd(_mm256_mul_pd(xi_210,_mm256_set_pd(2.0,2.0,2.0,2.0)),_mm256_mul_pd(xi_224,_mm256_set_pd(2.0,2.0,2.0,2.0)));
            const __m256d xi_118 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_207,_mm256_set_pd(-7.0,-7.0,-7.0,-7.0)),_mm256_mul_pd(xi_208,_mm256_set_pd(5.0,5.0,5.0,5.0))),_mm256_mul_pd(xi_209,_mm256_set_pd(-4.0,-4.0,-4.0,-4.0))),_mm256_mul_pd(xi_213,_mm256_set_pd(-7.0,-7.0,-7.0,-7.0))),_mm256_mul_pd(xi_215,_mm256_set_pd(-4.0,-4.0,-4.0,-4.0))),_mm256_mul_pd(xi_217,_mm256_set_pd(-7.0,-7.0,-7.0,-7.0))),_mm256_mul_pd(xi_222,_mm256_set_pd(-7.0,-7.0,-7.0,-7.0))),_mm256_mul_pd(xi_226,_mm256_set_pd(5.0,5.0,5.0,5.0))),xi_103),xi_115),xi_116),xi_117),xi_99),_mm256_set_pd(omega_even,omega_even,omega_even,omega_even));
            const __m256d xi_119 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_107,xi_208),xi_227),xi_23),xi_38);
            const __m256d xi_121 = _mm256_mul_pd(xi_119,_mm256_set_pd(xi_120,xi_120,xi_120,xi_120));
            const __m256d xi_122 = _mm256_mul_pd(xi_217,_mm256_set_pd(2.0,2.0,2.0,2.0));
            const __m256d xi_123 = _mm256_mul_pd(xi_213,_mm256_set_pd(2.0,2.0,2.0,2.0));
            const __m256d xi_124 = _mm256_add_pd(_mm256_mul_pd(xi_207,_mm256_set_pd(-2.0,-2.0,-2.0,-2.0)),_mm256_mul_pd(xi_222,_mm256_set_pd(2.0,2.0,2.0,2.0)));
            const __m256d xi_125 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_122,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_123),xi_124),xi_25),xi_4);
            const __m256d xi_127 = _mm256_mul_pd(xi_125,_mm256_set_pd(xi_126,xi_126,xi_126,xi_126));
            const __m256d xi_128 = _mm256_mul_pd(xi_127,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_133 = _mm256_mul_pd(xi_118,_mm256_set_pd(-0.0198412698412698,-0.0198412698412698,-0.0198412698412698,-0.0198412698412698));
            const __m256d xi_144 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_112,xi_13),xi_214),xi_223),xi_34);
            const __m256d xi_145 = _mm256_mul_pd(xi_144,_mm256_set_pd(xi_120,xi_120,xi_120,xi_120));
            const __m256d xi_146 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_123,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_122),xi_124),xi_14),xi_210),xi_214),xi_35);
            const __m256d xi_147 = _mm256_mul_pd(xi_146,_mm256_set_pd(xi_126,xi_126,xi_126,xi_126));
            const __m256d xi_149 = _mm256_mul_pd(xi_147,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_150 = _mm256_add_pd(xi_218,xi_228);
            const __m256d xi_151 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_108,xi_150),xi_209),xi_29);
            const __m256d xi_152 = _mm256_mul_pd(xi_151,_mm256_set_pd(xi_120,xi_120,xi_120,xi_120));
            const __m256d xi_153 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_115,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_116,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_117),xi_30),xi_6);
            const __m256d xi_154 = _mm256_mul_pd(xi_153,_mm256_set_pd(xi_126,xi_126,xi_126,xi_126));
            const __m256d xi_155 = _mm256_mul_pd(xi_154,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_162 = _mm256_mul_pd(xi_118,_mm256_set_pd(0.0138888888888889,0.0138888888888889,0.0138888888888889,0.0138888888888889));
            const __m256d xi_177 = _mm256_mul_pd(xi_104,_mm256_set_pd(-0.00714285714285714,-0.00714285714285714,-0.00714285714285714,-0.00714285714285714));
            const __m256d xi_179 = _mm256_mul_pd(xi_101,_mm256_set_pd(0.025,0.025,0.025,0.025));
            const __m256d xi_182 = _mm256_mul_pd(xi_153,_mm256_set_pd(xi_181,xi_181,xi_181,xi_181));
            const __m256d xi_184 = _mm256_mul_pd(xi_151,_mm256_set_pd(xi_183,xi_183,xi_183,xi_183));
            const __m256d xi_185 = _mm256_mul_pd(xi_118,_mm256_set_pd(-0.00396825396825397,-0.00396825396825397,-0.00396825396825397,-0.00396825396825397));
            const __m256d xi_189 = _mm256_mul_pd(xi_125,_mm256_set_pd(xi_181,xi_181,xi_181,xi_181));
            const __m256d xi_190 = _mm256_mul_pd(xi_119,_mm256_set_pd(xi_183,xi_183,xi_183,xi_183));
            const __m256d xi_196 = _mm256_mul_pd(xi_104,_mm256_set_pd(0.0178571428571429,0.0178571428571429,0.0178571428571429,0.0178571428571429));
            const __m256d xi_199 = _mm256_mul_pd(xi_144,_mm256_set_pd(xi_183,xi_183,xi_183,xi_183));
            const __m256d xi_200 = _mm256_mul_pd(xi_146,_mm256_set_pd(xi_181,xi_181,xi_181,xi_181));
            const __m256d vel0Term = xi_2;
            const __m256d vel1Term = xi_5;
            const __m256d vel2Term = xi_7;
            const __m256d rho = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(vel0Term,vel1Term),vel2Term),xi_10),xi_211),xi_8),xi_9);
            const __m256d xi_11 = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),rho);
            const __m256d u_0 = _mm256_mul_pd(xi_11,_mm256_add_pd(_mm256_add_pd(vel0Term,xi_14),xi_19));
            const __m256d xi_40 = _mm256_mul_pd(u_0,xi_225);
            const __m256d xi_41 = _mm256_mul_pd(xi_40,_mm256_set_pd(0.333333333333333,0.333333333333333,0.333333333333333,0.333333333333333));
            const __m256d xi_47 = _mm256_mul_pd(xi_41,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_105 = _mm256_mul_pd(rho,(_mm256_mul_pd(u_0,u_0)));
            const __m256d xi_141 = _mm256_mul_pd(rho,u_0);
            const __m256d xi_142 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(vel0Term,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_141),xi_32);
            const __m256d xi_143 = _mm256_mul_pd(xi_142,_mm256_set_pd(xi_131,xi_131,xi_131,xi_131));
            const __m256d xi_167 = _mm256_mul_pd(xi_142,_mm256_set_pd(xi_164,xi_164,xi_164,xi_164));
            const __m256d u_1 = _mm256_mul_pd(xi_11,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(vel1Term,xi_17),xi_21),xi_25));
            const __m256d xi_42 = _mm256_mul_pd(u_1,xi_219);
            const __m256d xi_43 = _mm256_mul_pd(xi_42,_mm256_set_pd(0.333333333333333,0.333333333333333,0.333333333333333,0.333333333333333));
            const __m256d xi_48 = _mm256_mul_pd(xi_43,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_69 = _mm256_mul_pd(u_1,_mm256_set_pd(0.5,0.5,0.5,0.5));
            const __m256d xi_72 = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(u_0,xi_70),_mm256_mul_pd(xi_225,xi_69)),_mm256_set_pd(xi_71,xi_71,xi_71,xi_71));
            const __m256d xi_73 = _mm256_mul_pd(xi_72,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_110 = _mm256_mul_pd(rho,(_mm256_mul_pd(u_1,u_1)));
            const __m256d xi_111 = _mm256_add_pd(_mm256_add_pd(xi_109,xi_110),xi_20);
            const __m256d xi_129 = _mm256_mul_pd(rho,u_1);
            const __m256d xi_130 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(vel1Term,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_109),xi_129),xi_213),xi_37);
            const __m256d xi_132 = _mm256_mul_pd(xi_130,_mm256_set_pd(xi_131,xi_131,xi_131,xi_131));
            const __m256d xi_165 = _mm256_mul_pd(xi_130,_mm256_set_pd(xi_164,xi_164,xi_164,xi_164));
            const __m256d xi_166 = _mm256_add_pd(xi_127,xi_165);
            const __m256d xi_175 = _mm256_mul_pd(xi_165,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_176 = _mm256_add_pd(xi_128,xi_175);
            const __m256d xi_191 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_189,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_165),xi_190);
            const __m256d xi_192 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_190,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_175),xi_189);
            const __m256d u_2 = _mm256_mul_pd(xi_11,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(vel2Term,xi_220),xi_27),xi_30));
            const __m256d xi_44 = _mm256_mul_pd(u_2,xi_221);
            const __m256d xi_45 = _mm256_mul_pd(xi_44,_mm256_set_pd(0.333333333333333,0.333333333333333,0.333333333333333,0.333333333333333));
            const __m256d xi_46 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(xi_41,xi_43),xi_45),_mm256_set_pd(-omega_bulk + 2.0,-omega_bulk + 2.0,-omega_bulk + 2.0,-omega_bulk + 2.0));
            const __m256d xi_49 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_44,_mm256_set_pd(0.666666666666667,0.666666666666667,0.666666666666667,0.666666666666667)),xi_47),xi_48);
            const __m256d xi_52 = _mm256_mul_pd(xi_45,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_53 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_42,_mm256_set_pd(0.666666666666667,0.666666666666667,0.666666666666667,0.666666666666667)),xi_47),xi_52);
            const __m256d xi_54 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_40,_mm256_set_pd(0.666666666666667,0.666666666666667,0.666666666666667,0.666666666666667)),xi_48),xi_52);
            const __m256d xi_57 = _mm256_mul_pd(xi_49,_mm256_set_pd(xi_56,xi_56,xi_56,xi_56));
            const __m256d xi_58 = _mm256_mul_pd(xi_57,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_59 = _mm256_mul_pd(xi_54,_mm256_set_pd(xi_56,xi_56,xi_56,xi_56));
            const __m256d xi_60 = _mm256_mul_pd(xi_59,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_62 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_53,_mm256_set_pd(xi_61,xi_61,xi_61,xi_61)),xi_58),xi_60);
            const __m256d xi_64 = _mm256_mul_pd(xi_53,_mm256_set_pd(xi_56,xi_56,xi_56,xi_56));
            const __m256d xi_65 = _mm256_mul_pd(xi_64,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_66 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_54,_mm256_set_pd(xi_61,xi_61,xi_61,xi_61)),xi_58),xi_65);
            const __m256d xi_68 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_49,_mm256_set_pd(xi_61,xi_61,xi_61,xi_61)),xi_60),xi_65);
            const __m256d xi_75 = _mm256_add_pd(_mm256_mul_pd(xi_74,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_59);
            const __m256d xi_77 = _mm256_mul_pd(_mm256_mul_pd(xi_49,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_76,xi_76,xi_76,xi_76));
            const __m256d xi_79 = _mm256_mul_pd(xi_46,_mm256_set_pd(0.125,0.125,0.125,0.125));
            const __m256d xi_80 = _mm256_add_pd(xi_64,xi_79);
            const __m256d xi_81 = _mm256_add_pd(xi_78,xi_80);
            const __m256d xi_82 = _mm256_add_pd(xi_77,xi_81);
            const __m256d xi_83 = _mm256_add_pd(xi_59,xi_74);
            const __m256d xi_84 = _mm256_add_pd(_mm256_mul_pd(xi_78,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_80);
            const __m256d xi_85 = _mm256_add_pd(xi_77,xi_84);
            const __m256d xi_86 = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(u_2,xi_70),_mm256_mul_pd(xi_221,xi_69)),_mm256_set_pd(xi_71,xi_71,xi_71,xi_71));
            const __m256d xi_87 = _mm256_mul_pd(_mm256_mul_pd(xi_54,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_76,xi_76,xi_76,xi_76));
            const __m256d xi_89 = _mm256_add_pd(xi_57,xi_88);
            const __m256d xi_90 = _mm256_add_pd(xi_87,xi_89);
            const __m256d xi_91 = _mm256_mul_pd(xi_86,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_92 = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(u_0,xi_221),_mm256_set_pd(0.5,0.5,0.5,0.5)),_mm256_mul_pd(_mm256_mul_pd(u_2,xi_225),_mm256_set_pd(0.5,0.5,0.5,0.5))),_mm256_set_pd(xi_71,xi_71,xi_71,xi_71));
            const __m256d xi_93 = _mm256_mul_pd(xi_92,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_94 = _mm256_mul_pd(_mm256_mul_pd(xi_53,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_76,xi_76,xi_76,xi_76));
            const __m256d xi_95 = _mm256_add_pd(_mm256_add_pd(xi_79,xi_89),xi_94);
            const __m256d xi_96 = _mm256_add_pd(_mm256_mul_pd(xi_88,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_57);
            const __m256d xi_97 = _mm256_add_pd(xi_87,xi_96);
            const __m256d xi_98 = _mm256_add_pd(_mm256_add_pd(xi_79,xi_94),xi_96);
            const __m256d xi_106 = _mm256_mul_pd(rho,(_mm256_mul_pd(u_2,u_2)));
            const __m256d xi_114 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_105,xi_106),xi_108),xi_111),xi_113),xi_211),xi_24),xi_28),_mm256_set_pd(omega_bulk,omega_bulk,omega_bulk,omega_bulk));
            const __m256d xi_134 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_106,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_209),xi_215);
            const __m256d xi_135 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_208,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_1),xi_111),xi_134),xi_18),xi_23),xi_31),_mm256_set_pd(omega_shear,omega_shear,omega_shear,omega_shear));
            const __m256d xi_136 = _mm256_mul_pd(xi_135,_mm256_set_pd(0.125,0.125,0.125,0.125));
            const __m256d xi_137 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_105,_mm256_set_pd(2.0,2.0,2.0,2.0)),_mm256_mul_pd(xi_110,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_214,_mm256_set_pd(-2.0,-2.0,-2.0,-2.0))),_mm256_mul_pd(xi_216,_mm256_set_pd(-2.0,-2.0,-2.0,-2.0))),xi_102),xi_109),xi_113),xi_134),xi_20),xi_208),xi_226),_mm256_set_pd(omega_shear,omega_shear,omega_shear,omega_shear));
            const __m256d xi_138 = _mm256_mul_pd(xi_137,_mm256_set_pd(-0.0416666666666667,-0.0416666666666667,-0.0416666666666667,-0.0416666666666667));
            const __m256d xi_139 = _mm256_add_pd(_mm256_mul_pd(xi_101,_mm256_set_pd(-0.05,-0.05,-0.05,-0.05)),xi_138);
            const __m256d xi_140 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_104,_mm256_set_pd(0.0142857142857143,0.0142857142857143,0.0142857142857143,0.0142857142857143)),xi_133),xi_136),xi_139);
            const __m256d xi_148 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_104,_mm256_set_pd(-0.0357142857142857,-0.0357142857142857,-0.0357142857142857,-0.0357142857142857)),_mm256_mul_pd(xi_137,_mm256_set_pd(0.0833333333333333,0.0833333333333333,0.0833333333333333,0.0833333333333333))),xi_133);
            const __m256d xi_156 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(rho,u_2),_mm256_mul_pd(vel2Term,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_10),xi_107),xi_112),xi_150),xi_210);
            const __m256d xi_157 = _mm256_mul_pd(xi_156,_mm256_set_pd(xi_131,xi_131,xi_131,xi_131));
            const __m256d xi_158 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_104,_mm256_set_pd(-0.0214285714285714,-0.0214285714285714,-0.0214285714285714,-0.0214285714285714)),_mm256_mul_pd(xi_118,_mm256_set_pd(0.0158730158730159,0.0158730158730159,0.0158730158730159,0.0158730158730159))),_mm256_mul_pd(xi_136,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_139);
            const __m256d xi_161 = _mm256_mul_pd(xi_135,_mm256_set_pd(0.0625,0.0625,0.0625,0.0625));
            const __m256d xi_168 = _mm256_mul_pd(xi_114,_mm256_set_pd(0.0416666666666667,0.0416666666666667,0.0416666666666667,0.0416666666666667));
            const __m256d xi_169 = _mm256_add_pd(_mm256_mul_pd(xi_137,_mm256_set_pd(0.0208333333333333,0.0208333333333333,0.0208333333333333,0.0208333333333333)),xi_168);
            const __m256d xi_170 = _mm256_add_pd(_mm256_mul_pd(xi_167,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_169);
            const __m256d xi_171 = _mm256_add_pd(xi_149,xi_170);
            const __m256d xi_173 = _mm256_add_pd(xi_167,xi_169);
            const __m256d xi_174 = _mm256_add_pd(xi_147,xi_173);
            const __m256d xi_186 = _mm256_mul_pd(xi_156,_mm256_set_pd(xi_164,xi_164,xi_164,xi_164));
            const __m256d xi_187 = _mm256_add_pd(xi_185,xi_186);
            const __m256d xi_188 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_182,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_184),xi_187);
            const __m256d xi_195 = _mm256_mul_pd(xi_161,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
            const __m256d xi_198 = _mm256_add_pd(xi_154,xi_187);
            const __m256d xi_201 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_199,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_170),xi_200);
            const __m256d xi_203 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_200,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_173),xi_199);
            const __m256d xi_204 = _mm256_add_pd(_mm256_mul_pd(xi_186,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_185);
            const __m256d xi_205 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_184,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_182),xi_204);
            const __m256d xi_206 = _mm256_add_pd(xi_155,xi_204);
            const __m256d p_0 = _mm256_add_pd(xi_2,xi_32);
            const __m256d p_1 = xi_33;
            const __m256d xi_160 = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(p_1,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_0,xi_129)),_mm256_set_pd(xi_159,xi_159,xi_159,xi_159));
            const __m256d xi_163 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_160,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_161),xi_162);
            const __m256d xi_172 = _mm256_add_pd(_mm256_add_pd(xi_160,xi_161),xi_162);
            const __m256d p_2 = xi_36;
            const __m256d xi_194 = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(p_2,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_2,xi_141)),_mm256_set_pd(xi_159,xi_159,xi_159,xi_159));
            const __m256d xi_197 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_194,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_195),xi_196);
            const __m256d xi_202 = _mm256_add_pd(_mm256_add_pd(xi_194,xi_195),xi_196);
            const __m256d p_3 = xi_33;
            const __m256d p_4 = _mm256_add_pd(_mm256_add_pd(xi_0,xi_37),xi_5);
            const __m256d p_5 = xi_39;
            const __m256d xi_178 = _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(p_5,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_2,xi_129)),_mm256_set_pd(xi_159,xi_159,xi_159,xi_159));
            const __m256d xi_180 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_138,xi_168),xi_177),xi_178),xi_179);
            const __m256d xi_193 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_178,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_138),xi_168),xi_177),xi_179);
            const __m256d p_6 = xi_36;
            const __m256d p_7 = xi_39;
            const __m256d p_8 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_1,xi_10),xi_228),xi_3),xi_7);
            const __m256d forceTerm_0 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_46,_mm256_set_pd(-1.5,-1.5,-1.5,-1.5)),_mm256_mul_pd(_mm256_mul_pd(xi_49,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_51,xi_51,xi_51,xi_51))),_mm256_mul_pd(_mm256_mul_pd(xi_53,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_51,xi_51,xi_51,xi_51))),_mm256_mul_pd(_mm256_mul_pd(xi_54,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_51,xi_51,xi_51,xi_51)));
            const __m256d forceTerm_1 = _mm256_add_pd(xi_55,xi_62);
            const __m256d forceTerm_2 = _mm256_add_pd(_mm256_mul_pd(xi_55,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_62);
            const __m256d forceTerm_3 = _mm256_add_pd(_mm256_mul_pd(xi_63,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_66);
            const __m256d forceTerm_4 = _mm256_add_pd(xi_63,xi_66);
            const __m256d forceTerm_5 = _mm256_add_pd(xi_67,xi_68);
            const __m256d forceTerm_6 = _mm256_add_pd(_mm256_mul_pd(xi_67,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_68);
            const __m256d forceTerm_7 = _mm256_add_pd(_mm256_add_pd(xi_73,xi_75),xi_82);
            const __m256d forceTerm_8 = _mm256_add_pd(_mm256_add_pd(xi_72,xi_82),xi_83);
            const __m256d forceTerm_9 = _mm256_add_pd(_mm256_add_pd(xi_72,xi_75),xi_85);
            const __m256d forceTerm_10 = _mm256_add_pd(_mm256_add_pd(xi_73,xi_83),xi_85);
            const __m256d forceTerm_11 = _mm256_add_pd(_mm256_add_pd(xi_81,xi_86),xi_90);
            const __m256d forceTerm_12 = _mm256_add_pd(_mm256_add_pd(xi_84,xi_90),xi_91);
            const __m256d forceTerm_13 = _mm256_add_pd(_mm256_add_pd(xi_75,xi_93),xi_95);
            const __m256d forceTerm_14 = _mm256_add_pd(_mm256_add_pd(xi_83,xi_92),xi_95);
            const __m256d forceTerm_15 = _mm256_add_pd(_mm256_add_pd(xi_81,xi_91),xi_97);
            const __m256d forceTerm_16 = _mm256_add_pd(_mm256_add_pd(xi_84,xi_86),xi_97);
            const __m256d forceTerm_17 = _mm256_add_pd(_mm256_add_pd(xi_75,xi_92),xi_98);
            const __m256d forceTerm_18 = _mm256_add_pd(_mm256_add_pd(xi_83,xi_93),xi_98);
            _mm256_store_pd(&_data_pdfs_20_30_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_101,_mm256_set_pd(0.1,0.1,0.1,0.1)),_mm256_mul_pd(xi_104,_mm256_set_pd(0.0428571428571429,0.0428571428571429,0.0428571428571429,0.0428571428571429))),_mm256_mul_pd(xi_114,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5))),_mm256_mul_pd(xi_118,_mm256_set_pd(0.0238095238095238,0.0238095238095238,0.0238095238095238,0.0238095238095238))),forceTerm_0),xi_211));
            _mm256_store_pd(&_data_pdfs_20_31_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_121,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),forceTerm_1),xi_128),xi_132),xi_140),xi_208));
            _mm256_store_pd(&_data_pdfs_20_32_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_132,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),forceTerm_2),xi_121),xi_127),xi_140),xi_226));
            _mm256_store_pd(&_data_pdfs_20_33_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_143,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),forceTerm_3),xi_145),xi_147),xi_148),xi_216));
            _mm256_store_pd(&_data_pdfs_20_34_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_145,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),forceTerm_4),xi_143),xi_148),xi_149),xi_214));
            _mm256_store_pd(&_data_pdfs_20_35_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_152,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),forceTerm_5),xi_155),xi_157),xi_158),xi_209));
            _mm256_store_pd(&_data_pdfs_20_36_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_157,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),forceTerm_6),xi_152),xi_154),xi_158),xi_215));
            _mm256_store_pd(&_data_pdfs_20_37_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_7,xi_163),xi_166),xi_171),xi_217));
            _mm256_store_pd(&_data_pdfs_20_38_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_8,xi_166),xi_172),xi_174),xi_207));
            _mm256_store_pd(&_data_pdfs_20_39_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_9,xi_171),xi_172),xi_176),xi_222));
            _mm256_store_pd(&_data_pdfs_20_310_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_10,xi_163),xi_174),xi_176),xi_213));
            _mm256_store_pd(&_data_pdfs_20_311_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_11,xi_180),xi_188),xi_191),xi_212));
            _mm256_store_pd(&_data_pdfs_20_312_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_12,xi_188),xi_192),xi_193),xi_227));
            _mm256_store_pd(&_data_pdfs_20_313_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_13,xi_197),xi_198),xi_201),xi_223));
            _mm256_store_pd(&_data_pdfs_20_314_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_14,xi_198),xi_202),xi_203),xi_220));
            _mm256_store_pd(&_data_pdfs_20_315_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_15,xi_191),xi_193),xi_205),xi_218));
            _mm256_store_pd(&_data_pdfs_20_316_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_16,xi_180),xi_192),xi_205),xi_228));
            _mm256_store_pd(&_data_pdfs_20_317_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_17,xi_201),xi_202),xi_206),xi_224));
            _mm256_store_pd(&_data_pdfs_20_318_10[ctr_0],_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_18,xi_197),xi_203),xi_206),xi_210));
         }
      }
   }
}
}
namespace internal_kernel_stream {
static FUNC_PREFIX void kernel_stream(double * RESTRICT const _data_pdfs, double * RESTRICT _data_pdfs_tmp, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3)
{
   for (int64_t ctr_2 = 1; ctr_2 < _size_pdfs_2 - 1; ctr_2 += 1)
   {
      double * RESTRICT _data_pdfs_tmp_20_30 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2;
      double * RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * RESTRICT _data_pdfs_tmp_20_31 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + _stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_32 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 2*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_33 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 3*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_34 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 4*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_35 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 5*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_36 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 6*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_37 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 7*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_38 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 8*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_39 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 9*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_310 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 10*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_311 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 11*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_312 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 12*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_313 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 13*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_314 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 14*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_315 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 15*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_316 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 16*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_317 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 17*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_318 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 18*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3;
      for (int64_t ctr_1 = 1; ctr_1 < _size_pdfs_1 - 1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_tmp_20_30_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_30;
         double * RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * RESTRICT _data_pdfs_tmp_20_31_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_31;
         double * RESTRICT _data_pdfs_20_31_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
         double * RESTRICT _data_pdfs_tmp_20_32_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_32;
         double * RESTRICT _data_pdfs_20_32_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
         double * RESTRICT _data_pdfs_tmp_20_33_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_33;
         double * RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * RESTRICT _data_pdfs_tmp_20_34_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_34;
         double * RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * RESTRICT _data_pdfs_tmp_20_35_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_35;
         double * RESTRICT _data_pdfs_2m1_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_35;
         double * RESTRICT _data_pdfs_tmp_20_36_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_36;
         double * RESTRICT _data_pdfs_21_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_36;
         double * RESTRICT _data_pdfs_tmp_20_37_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_37;
         double * RESTRICT _data_pdfs_20_37_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
         double * RESTRICT _data_pdfs_tmp_20_38_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_38;
         double * RESTRICT _data_pdfs_20_38_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
         double * RESTRICT _data_pdfs_tmp_20_39_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_39;
         double * RESTRICT _data_pdfs_20_39_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
         double * RESTRICT _data_pdfs_tmp_20_310_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_310;
         double * RESTRICT _data_pdfs_20_310_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
         double * RESTRICT _data_pdfs_tmp_20_311_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_311;
         double * RESTRICT _data_pdfs_2m1_311_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
         double * RESTRICT _data_pdfs_tmp_20_312_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_312;
         double * RESTRICT _data_pdfs_2m1_312_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
         double * RESTRICT _data_pdfs_tmp_20_313_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_313;
         double * RESTRICT _data_pdfs_2m1_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_313;
         double * RESTRICT _data_pdfs_tmp_20_314_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_314;
         double * RESTRICT _data_pdfs_2m1_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_314;
         double * RESTRICT _data_pdfs_tmp_20_315_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_315;
         double * RESTRICT _data_pdfs_21_315_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
         double * RESTRICT _data_pdfs_tmp_20_316_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_316;
         double * RESTRICT _data_pdfs_21_316_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
         double * RESTRICT _data_pdfs_tmp_20_317_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_317;
         double * RESTRICT _data_pdfs_21_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_317;
         double * RESTRICT _data_pdfs_tmp_20_318_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_318;
         double * RESTRICT _data_pdfs_21_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_318;
         for (int64_t ctr_0 = 1; ctr_0 < ((_size_pdfs_0 - 2) % (4) == 0 ? _size_pdfs_0 - 2 : ((int64_t)((_size_pdfs_0 - 2) / (4))+1) * (4)) + 1; ctr_0 += 4)
         {
            _mm256_store_pd(&_data_pdfs_tmp_20_30_10[ctr_0],_mm256_load_pd(& _data_pdfs_20_30_10[ctr_0]));
            _mm256_store_pd(&_data_pdfs_tmp_20_31_10[ctr_0],_mm256_load_pd(& _data_pdfs_20_31_1m1[ctr_0]));
            _mm256_store_pd(&_data_pdfs_tmp_20_32_10[ctr_0],_mm256_load_pd(& _data_pdfs_20_32_11[ctr_0]));
            _mm256_store_pd(&_data_pdfs_tmp_20_33_10[ctr_0],_mm256_loadu_pd(& _data_pdfs_20_33_10[ctr_0 + 1]));
            _mm256_store_pd(&_data_pdfs_tmp_20_34_10[ctr_0],_mm256_loadu_pd(& _data_pdfs_20_34_10[ctr_0 - 1]));
            _mm256_store_pd(&_data_pdfs_tmp_20_35_10[ctr_0],_mm256_load_pd(& _data_pdfs_2m1_35_10[ctr_0]));
            _mm256_store_pd(&_data_pdfs_tmp_20_36_10[ctr_0],_mm256_load_pd(& _data_pdfs_21_36_10[ctr_0]));
            _mm256_store_pd(&_data_pdfs_tmp_20_37_10[ctr_0],_mm256_loadu_pd(& _data_pdfs_20_37_1m1[ctr_0 + 1]));
            _mm256_store_pd(&_data_pdfs_tmp_20_38_10[ctr_0],_mm256_loadu_pd(& _data_pdfs_20_38_1m1[ctr_0 - 1]));
            _mm256_store_pd(&_data_pdfs_tmp_20_39_10[ctr_0],_mm256_loadu_pd(& _data_pdfs_20_39_11[ctr_0 + 1]));
            _mm256_store_pd(&_data_pdfs_tmp_20_310_10[ctr_0],_mm256_loadu_pd(& _data_pdfs_20_310_11[ctr_0 - 1]));
            _mm256_store_pd(&_data_pdfs_tmp_20_311_10[ctr_0],_mm256_load_pd(& _data_pdfs_2m1_311_1m1[ctr_0]));
            _mm256_store_pd(&_data_pdfs_tmp_20_312_10[ctr_0],_mm256_load_pd(& _data_pdfs_2m1_312_11[ctr_0]));
            _mm256_store_pd(&_data_pdfs_tmp_20_313_10[ctr_0],_mm256_loadu_pd(& _data_pdfs_2m1_313_10[ctr_0 + 1]));
            _mm256_store_pd(&_data_pdfs_tmp_20_314_10[ctr_0],_mm256_loadu_pd(& _data_pdfs_2m1_314_10[ctr_0 - 1]));
            _mm256_store_pd(&_data_pdfs_tmp_20_315_10[ctr_0],_mm256_load_pd(& _data_pdfs_21_315_1m1[ctr_0]));
            _mm256_store_pd(&_data_pdfs_tmp_20_316_10[ctr_0],_mm256_load_pd(& _data_pdfs_21_316_11[ctr_0]));
            _mm256_store_pd(&_data_pdfs_tmp_20_317_10[ctr_0],_mm256_loadu_pd(& _data_pdfs_21_317_10[ctr_0 + 1]));
            _mm256_store_pd(&_data_pdfs_tmp_20_318_10[ctr_0],_mm256_loadu_pd(& _data_pdfs_21_318_10[ctr_0 - 1]));
         }
      }
   }
}
}


const real_t MRTLatticeModelAvx::w[19] = { 0.333333333333333,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778 };
const real_t MRTLatticeModelAvx::wInv[19] = { 3.00000000000000,18.0000000000000,18.0000000000000,18.0000000000000,18.0000000000000,18.0000000000000,18.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000 };

void MRTLatticeModelAvx::Sweep::streamCollide( IBlock * block, const uint_t numberOfGhostLayersToInclude )
{
    auto pdfs = block->getData< field::GhostLayerField<double, 19> >(pdfsID);
    field::GhostLayerField<double, 19> * pdfs_tmp;
    {
        // Getting temporary field pdfs_tmp
        auto it = cache_pdfs_.find( pdfs );
        if( it != cache_pdfs_.end() )
        {
            pdfs_tmp = *it;
        }
        else
        {
            pdfs_tmp = pdfs->cloneUninitialized();
            cache_pdfs_.insert(pdfs_tmp);
        }
    }


    auto & lm = dynamic_cast< lbm::PdfField<MRTLatticeModelAvx> * > (pdfs)->latticeModel();
    WALBERLA_ASSERT_EQUAL( *(lm.blockId_), block->getId() );

    auto & force = lm.force_;
    auto & omega_even = lm.omega_even_;
    auto & omega_odd = lm.omega_odd_;
    auto & omega_shear = lm.omega_shear_;
    auto & omega_bulk = lm.omega_bulk_;
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(force->nrOfGhostLayers()));
    double * RESTRICT const _data_force = force->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs->nrOfGhostLayers()));
    double * RESTRICT const _data_pdfs = pdfs->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs_tmp->nrOfGhostLayers()));
    double * RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(cell_idx_c(force->xSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_force_0 = int64_t(cell_idx_c(force->xSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(cell_idx_c(force->ySize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_force_1 = int64_t(cell_idx_c(force->ySize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(cell_idx_c(force->zSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_force_2 = int64_t(cell_idx_c(force->zSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    const int64_t _stride_force_1 = int64_t(force->yStride());
    const int64_t _stride_force_2 = int64_t(force->zStride());
    const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
    const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
    const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
    internal_kernel_streamCollide::kernel_streamCollide(_data_force, _data_pdfs, _data_pdfs_tmp, _size_force_0, _size_force_1, _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, omega_bulk, omega_even, omega_odd, omega_shear);
    pdfs->swapDataPointers(pdfs_tmp);

}

void MRTLatticeModelAvx::Sweep::collide( IBlock * block, const uint_t numberOfGhostLayersToInclude )
{
   auto pdfs = block->getData< field::GhostLayerField<double, 19> >(pdfsID);


    auto & lm = dynamic_cast< lbm::PdfField<MRTLatticeModelAvx> * > (pdfs)->latticeModel();
    WALBERLA_ASSERT_EQUAL( *(lm.blockId_), block->getId() );

    auto & force = lm.force_;
    auto & omega_even = lm.omega_even_;
    auto & omega_odd = lm.omega_odd_;
    auto & omega_shear = lm.omega_shear_;
    auto & omega_bulk = lm.omega_bulk_;
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude), -int_c(force->nrOfGhostLayers()));
    double * RESTRICT const _data_force = force->dataAt(-cell_idx_c(numberOfGhostLayersToInclude), -cell_idx_c(numberOfGhostLayersToInclude), -cell_idx_c(numberOfGhostLayersToInclude), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude), -int_c(pdfs->nrOfGhostLayers()));
    double * RESTRICT _data_pdfs = pdfs->dataAt(-cell_idx_c(numberOfGhostLayersToInclude), -cell_idx_c(numberOfGhostLayersToInclude), -cell_idx_c(numberOfGhostLayersToInclude), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(cell_idx_c(force->xSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude)));
    const int64_t _size_force_0 = int64_t(cell_idx_c(force->xSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude));
    WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(cell_idx_c(force->ySize()) + 2*cell_idx_c(numberOfGhostLayersToInclude)));
    const int64_t _size_force_1 = int64_t(cell_idx_c(force->ySize()) + 2*cell_idx_c(numberOfGhostLayersToInclude));
    WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(cell_idx_c(force->zSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude)));
    const int64_t _size_force_2 = int64_t(cell_idx_c(force->zSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude));
    const int64_t _stride_force_1 = int64_t(force->yStride());
    const int64_t _stride_force_2 = int64_t(force->zStride());
    const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_kernel_collide::kernel_collide(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, omega_bulk, omega_even, omega_odd, omega_shear);
}


void MRTLatticeModelAvx::Sweep::stream( IBlock * block, const uint_t numberOfGhostLayersToInclude )
{
    auto pdfs = block->getData< field::GhostLayerField<double, 19> >(pdfsID);
    field::GhostLayerField<double, 19> * pdfs_tmp;
    {
        // Getting temporary field pdfs_tmp
        auto it = cache_pdfs_.find( pdfs );
        if( it != cache_pdfs_.end() )
        {
            pdfs_tmp = *it;
        }
        else
        {
            pdfs_tmp = pdfs->cloneUninitialized();
            cache_pdfs_.insert(pdfs_tmp);
        }
    }


    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs->nrOfGhostLayers()));
    double * RESTRICT const _data_pdfs = pdfs->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs_tmp->nrOfGhostLayers()));
    double * RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(pdfs->xSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(pdfs->xSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(pdfs->ySize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(pdfs->ySize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(pdfs->zSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(pdfs->zSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
    const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
    const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
    internal_kernel_stream::kernel_stream(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3);
    
    pdfs->swapDataPointers(pdfs_tmp);

}


} // namespace lbm
} // namespace walberla




// Buffer Packing

namespace walberla {
namespace mpi {

mpi::SendBuffer & operator<< (mpi::SendBuffer & buf, const ::walberla::lbm::MRTLatticeModelAvx & lm)
{
    buf << lm.currentLevel;
    return buf;
}

mpi::RecvBuffer & operator>> (mpi::RecvBuffer & buf, ::walberla::lbm::MRTLatticeModelAvx & lm)
{
    buf >> lm.currentLevel;
    return buf;
}


} // namespace mpi
} // namespace walberla

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif