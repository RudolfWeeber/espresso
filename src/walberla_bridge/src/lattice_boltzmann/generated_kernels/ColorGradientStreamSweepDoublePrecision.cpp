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
//! \\file ColorGradientStreamSweepDoublePrecision.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 17fc54c872bd8ceabf271a7e9e636c7c583f55af


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "ColorGradientStreamSweepDoublePrecision.h"




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


namespace internal_e4397140a5a3f258ae7331867e48f30b {
static FUNC_PREFIX void colorgradientstreamsweepdoubleprecision_colorgradientstreamsweepdoubleprecision(double * RESTRICT const _data_force_a, double * RESTRICT const _data_force_b, double * RESTRICT const _data_pdfs_a, double * RESTRICT  _data_pdfs_a_tmp, double * RESTRICT const _data_pdfs_b, double * RESTRICT  _data_pdfs_b_tmp, double * RESTRICT  _data_phasefield, double * RESTRICT  _data_rho_a, double * RESTRICT  _data_rho_b, double * RESTRICT  _data_velocity, int64_t const _size_force_a_0, int64_t const _size_force_a_1, int64_t const _size_force_a_2, int64_t const _stride_force_a_0, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_0, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_0, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_a_tmp_0, int64_t const _stride_pdfs_a_tmp_1, int64_t const _stride_pdfs_a_tmp_2, int64_t const _stride_pdfs_a_tmp_3, int64_t const _stride_pdfs_b_0, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_pdfs_b_tmp_0, int64_t const _stride_pdfs_b_tmp_1, int64_t const _stride_pdfs_b_tmp_2, int64_t const _stride_pdfs_b_tmp_3, int64_t const _stride_phasefield_0, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_0, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_0, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3)
{
#ifdef _OPENMP
   #pragma omp parallel
#endif
   {
#ifdef _OPENMP
      #pragma omp for schedule(static)
#endif
      for (int64_t ctr_2 = 1; ctr_2 < _size_force_a_2 - 1; ctr_2 += 1)
      {
         for (int64_t ctr_1 = 1; ctr_1 < _size_force_a_1 - 1; ctr_1 += 1)
         {
            for (int64_t ctr_0 = 1; ctr_0 < _size_force_a_0 - 1; ctr_0 += 1)
            {
               const double streamed_0_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2];
               const double streamed_1_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 - _stride_pdfs_a_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3];
               const double streamed_2_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3];
               const double streamed_3_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3];
               const double streamed_4_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 - _stride_pdfs_a_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3];
               const double streamed_5_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 - _stride_pdfs_a_2 + 5*_stride_pdfs_a_3];
               const double streamed_6_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_2 + 6*_stride_pdfs_a_3];
               const double streamed_7_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_0 + _stride_pdfs_a_1*ctr_1 - _stride_pdfs_a_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3];
               const double streamed_8_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 - _stride_pdfs_a_0 + _stride_pdfs_a_1*ctr_1 - _stride_pdfs_a_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3];
               const double streamed_9_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3];
               const double streamed_10_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 - _stride_pdfs_a_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3];
               const double streamed_11_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 - _stride_pdfs_a_1 + _stride_pdfs_a_2*ctr_2 - _stride_pdfs_a_2 + 11*_stride_pdfs_a_3];
               const double streamed_12_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_1 + _stride_pdfs_a_2*ctr_2 - _stride_pdfs_a_2 + 12*_stride_pdfs_a_3];
               const double streamed_13_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 - _stride_pdfs_a_2 + 13*_stride_pdfs_a_3];
               const double streamed_14_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 - _stride_pdfs_a_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 - _stride_pdfs_a_2 + 14*_stride_pdfs_a_3];
               const double streamed_15_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 - _stride_pdfs_a_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_2 + 15*_stride_pdfs_a_3];
               const double streamed_16_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_2 + 16*_stride_pdfs_a_3];
               const double streamed_17_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_2 + 17*_stride_pdfs_a_3];
               const double streamed_18_a = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 - _stride_pdfs_a_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_2 + 18*_stride_pdfs_a_3];
               const double vel0Term_a = streamed_10_a + streamed_14_a + streamed_18_a + streamed_4_a + streamed_8_a;
               const double momdensity_0_a = -streamed_13_a - streamed_17_a - streamed_3_a - streamed_7_a - streamed_9_a + vel0Term_a;
               const double vel1Term_a = streamed_11_a + streamed_15_a + streamed_1_a + streamed_7_a;
               const double momdensity_1_a = -streamed_10_a - streamed_12_a - streamed_16_a - streamed_2_a + streamed_8_a - streamed_9_a + vel1Term_a;
               const double vel2Term_a = streamed_12_a + streamed_13_a + streamed_5_a;
               const double rho_a = streamed_0_a + streamed_16_a + streamed_17_a + streamed_2_a + streamed_3_a + streamed_6_a + streamed_9_a + vel0Term_a + vel1Term_a + vel2Term_a;
               const double momdensity_2_a = streamed_11_a + streamed_14_a - streamed_15_a - streamed_16_a - streamed_17_a - streamed_18_a - streamed_6_a + vel2Term_a;
               const double u_0_a = momdensity_0_a*((1.0) / (rho_a)) + 0.5*((1.0) / (rho_a))*_data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2];
               const double u_1_a = momdensity_1_a*((1.0) / (rho_a)) + 0.5*((1.0) / (rho_a))*_data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3];
               const double u_2_a = momdensity_2_a*((1.0) / (rho_a)) + 0.5*((1.0) / (rho_a))*_data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3];
               const double rho_tmp_a = rho_a;
               const double vel_0_a = u_0_a;
               const double vel_1_a = u_1_a;
               const double vel_2_a = u_2_a;
               const double streamed_0_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2];
               const double streamed_1_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 - _stride_pdfs_b_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3];
               const double streamed_2_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3];
               const double streamed_3_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3];
               const double streamed_4_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 - _stride_pdfs_b_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3];
               const double streamed_5_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 - _stride_pdfs_b_2 + 5*_stride_pdfs_b_3];
               const double streamed_6_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_2 + 6*_stride_pdfs_b_3];
               const double streamed_7_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_0 + _stride_pdfs_b_1*ctr_1 - _stride_pdfs_b_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3];
               const double streamed_8_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 - _stride_pdfs_b_0 + _stride_pdfs_b_1*ctr_1 - _stride_pdfs_b_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3];
               const double streamed_9_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3];
               const double streamed_10_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 - _stride_pdfs_b_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3];
               const double streamed_11_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 - _stride_pdfs_b_1 + _stride_pdfs_b_2*ctr_2 - _stride_pdfs_b_2 + 11*_stride_pdfs_b_3];
               const double streamed_12_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_1 + _stride_pdfs_b_2*ctr_2 - _stride_pdfs_b_2 + 12*_stride_pdfs_b_3];
               const double streamed_13_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 - _stride_pdfs_b_2 + 13*_stride_pdfs_b_3];
               const double streamed_14_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 - _stride_pdfs_b_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 - _stride_pdfs_b_2 + 14*_stride_pdfs_b_3];
               const double streamed_15_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 - _stride_pdfs_b_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_2 + 15*_stride_pdfs_b_3];
               const double streamed_16_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_2 + 16*_stride_pdfs_b_3];
               const double streamed_17_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_2 + 17*_stride_pdfs_b_3];
               const double streamed_18_b = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 - _stride_pdfs_b_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_2 + 18*_stride_pdfs_b_3];
               const double vel0Term_b = streamed_10_b + streamed_14_b + streamed_18_b + streamed_4_b + streamed_8_b;
               const double momdensity_0_b = -streamed_13_b - streamed_17_b - streamed_3_b - streamed_7_b - streamed_9_b + vel0Term_b;
               const double vel1Term_b = streamed_11_b + streamed_15_b + streamed_1_b + streamed_7_b;
               const double momdensity_1_b = -streamed_10_b - streamed_12_b - streamed_16_b - streamed_2_b + streamed_8_b - streamed_9_b + vel1Term_b;
               const double vel2Term_b = streamed_12_b + streamed_13_b + streamed_5_b;
               const double rho_b = streamed_0_b + streamed_16_b + streamed_17_b + streamed_2_b + streamed_3_b + streamed_6_b + streamed_9_b + vel0Term_b + vel1Term_b + vel2Term_b;
               const double momdensity_2_b = streamed_11_b + streamed_14_b - streamed_15_b - streamed_16_b - streamed_17_b - streamed_18_b - streamed_6_b + vel2Term_b;
               const double u_0_b = momdensity_0_b*((1.0) / (rho_b)) + 0.5*((1.0) / (rho_b))*_data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2];
               const double u_1_b = momdensity_1_b*((1.0) / (rho_b)) + 0.5*((1.0) / (rho_b))*_data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3];
               const double u_2_b = momdensity_2_b*((1.0) / (rho_b)) + 0.5*((1.0) / (rho_b))*_data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3];
               const double rho_tmp_b = rho_b;
               const double vel_0_b = u_0_b;
               const double vel_1_b = u_1_b;
               const double vel_2_b = u_2_b;
               const double rho_total_inv = ((1.0) / (rho_tmp_a + rho_tmp_b));
               const double phi = rho_total_inv*(rho_tmp_a - rho_tmp_b);
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2] = streamed_0_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + _stride_pdfs_a_tmp_3] = streamed_1_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 2*_stride_pdfs_a_tmp_3] = streamed_2_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 3*_stride_pdfs_a_tmp_3] = streamed_3_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 4*_stride_pdfs_a_tmp_3] = streamed_4_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 5*_stride_pdfs_a_tmp_3] = streamed_5_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 6*_stride_pdfs_a_tmp_3] = streamed_6_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 7*_stride_pdfs_a_tmp_3] = streamed_7_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 8*_stride_pdfs_a_tmp_3] = streamed_8_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 9*_stride_pdfs_a_tmp_3] = streamed_9_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 10*_stride_pdfs_a_tmp_3] = streamed_10_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 11*_stride_pdfs_a_tmp_3] = streamed_11_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 12*_stride_pdfs_a_tmp_3] = streamed_12_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 13*_stride_pdfs_a_tmp_3] = streamed_13_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 14*_stride_pdfs_a_tmp_3] = streamed_14_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 15*_stride_pdfs_a_tmp_3] = streamed_15_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 16*_stride_pdfs_a_tmp_3] = streamed_16_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 17*_stride_pdfs_a_tmp_3] = streamed_17_a;
               _data_pdfs_a_tmp[_stride_pdfs_a_tmp_0*ctr_0 + _stride_pdfs_a_tmp_1*ctr_1 + _stride_pdfs_a_tmp_2*ctr_2 + 18*_stride_pdfs_a_tmp_3] = streamed_18_a;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2] = streamed_0_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + _stride_pdfs_b_tmp_3] = streamed_1_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 2*_stride_pdfs_b_tmp_3] = streamed_2_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 3*_stride_pdfs_b_tmp_3] = streamed_3_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 4*_stride_pdfs_b_tmp_3] = streamed_4_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 5*_stride_pdfs_b_tmp_3] = streamed_5_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 6*_stride_pdfs_b_tmp_3] = streamed_6_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 7*_stride_pdfs_b_tmp_3] = streamed_7_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 8*_stride_pdfs_b_tmp_3] = streamed_8_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 9*_stride_pdfs_b_tmp_3] = streamed_9_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 10*_stride_pdfs_b_tmp_3] = streamed_10_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 11*_stride_pdfs_b_tmp_3] = streamed_11_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 12*_stride_pdfs_b_tmp_3] = streamed_12_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 13*_stride_pdfs_b_tmp_3] = streamed_13_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 14*_stride_pdfs_b_tmp_3] = streamed_14_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 15*_stride_pdfs_b_tmp_3] = streamed_15_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 16*_stride_pdfs_b_tmp_3] = streamed_16_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 17*_stride_pdfs_b_tmp_3] = streamed_17_b;
               _data_pdfs_b_tmp[_stride_pdfs_b_tmp_0*ctr_0 + _stride_pdfs_b_tmp_1*ctr_1 + _stride_pdfs_b_tmp_2*ctr_2 + 18*_stride_pdfs_b_tmp_3] = streamed_18_b;
               _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] = rho_tmp_a;
               _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2] = rho_tmp_b;
               _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2] = rho_total_inv*(rho_tmp_a*vel_0_a + rho_tmp_b*vel_0_b);
               _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3] = rho_total_inv*(rho_tmp_a*vel_1_a + rho_tmp_b*vel_1_b);
               _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3] = rho_total_inv*(rho_tmp_a*vel_2_a + rho_tmp_b*vel_2_b);
               _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] = phi;
            }
         }
      }
   }
}
}


void ColorGradientStreamSweepDoublePrecision::run(IBlock * block)
{
   
    auto rho_a = block->getData< field::GhostLayerField<double, 1> >(rho_aID);
    auto phasefield = block->getData< field::GhostLayerField<double, 1> >(phasefieldID);
    auto force_a = block->getData< field::GhostLayerField<double, 3> >(force_aID);
    auto pdfs_a = block->getData< field::GhostLayerField<double, 19> >(pdfs_aID);
    auto rho_b = block->getData< field::GhostLayerField<double, 1> >(rho_bID);
    auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
    auto pdfs_b = block->getData< field::GhostLayerField<double, 19> >(pdfs_bID);
    auto force_b = block->getData< field::GhostLayerField<double, 3> >(force_bID);
    field::GhostLayerField<double, 19> * pdfs_b_tmp;
    {
        if (cache_pdfs_b_.find(block) == cache_pdfs_b_.end())
        {
            pdfs_b_tmp = pdfs_b->cloneUninitialized();
            cache_pdfs_b_[block] = pdfs_b_tmp;
        }
        else
        {
            pdfs_b_tmp = cache_pdfs_b_[block];
        }
    }

    field::GhostLayerField<double, 19> * pdfs_a_tmp;
    {
        if (cache_pdfs_a_.find(block) == cache_pdfs_a_.end())
        {
            pdfs_a_tmp = pdfs_a->cloneUninitialized();
            cache_pdfs_a_[block] = pdfs_a_tmp;
        }
        else
        {
            pdfs_a_tmp = cache_pdfs_a_[block];
        }
    }

    
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_a->nrOfGhostLayers()))
    double * RESTRICT const _data_force_a = force_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT const _data_pdfs_a = pdfs_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_a_tmp->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a_tmp = pdfs_a_tmp->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a_tmp->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT const _data_pdfs_b = pdfs_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_b_tmp->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b_tmp = pdfs_b_tmp->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b_tmp->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT  _data_phasefield = phasefield->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT  _data_rho_a = rho_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT  _data_rho_b = rho_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT  _data_velocity = velocity->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(force_a->xSize()) + 2))
    const int64_t _size_force_a_0 = int64_t(int64_c(force_a->xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(force_a->ySize()) + 2))
    const int64_t _size_force_a_1 = int64_t(int64_c(force_a->ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(force_a->zSize()) + 2))
    const int64_t _size_force_a_2 = int64_t(int64_c(force_a->zSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    const int64_t _stride_force_a_0 = int64_t(force_a->xStride());
    const int64_t _stride_force_a_1 = int64_t(force_a->yStride());
    const int64_t _stride_force_a_2 = int64_t(force_a->zStride());
    const int64_t _stride_force_a_3 = int64_t(1 * int64_t(force_a->fStride()));
    const int64_t _stride_force_b_0 = int64_t(force_b->xStride());
    const int64_t _stride_force_b_1 = int64_t(force_b->yStride());
    const int64_t _stride_force_b_2 = int64_t(force_b->zStride());
    const int64_t _stride_force_b_3 = int64_t(1 * int64_t(force_b->fStride()));
    const int64_t _stride_pdfs_a_0 = int64_t(pdfs_a->xStride());
    const int64_t _stride_pdfs_a_1 = int64_t(pdfs_a->yStride());
    const int64_t _stride_pdfs_a_2 = int64_t(pdfs_a->zStride());
    const int64_t _stride_pdfs_a_3 = int64_t(1 * int64_t(pdfs_a->fStride()));
    const int64_t _stride_pdfs_a_tmp_0 = int64_t(pdfs_a_tmp->xStride());
    const int64_t _stride_pdfs_a_tmp_1 = int64_t(pdfs_a_tmp->yStride());
    const int64_t _stride_pdfs_a_tmp_2 = int64_t(pdfs_a_tmp->zStride());
    const int64_t _stride_pdfs_a_tmp_3 = int64_t(1 * int64_t(pdfs_a_tmp->fStride()));
    const int64_t _stride_pdfs_b_0 = int64_t(pdfs_b->xStride());
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
    const int64_t _stride_pdfs_b_tmp_0 = int64_t(pdfs_b_tmp->xStride());
    const int64_t _stride_pdfs_b_tmp_1 = int64_t(pdfs_b_tmp->yStride());
    const int64_t _stride_pdfs_b_tmp_2 = int64_t(pdfs_b_tmp->zStride());
    const int64_t _stride_pdfs_b_tmp_3 = int64_t(1 * int64_t(pdfs_b_tmp->fStride()));
    const int64_t _stride_phasefield_0 = int64_t(phasefield->xStride());
    const int64_t _stride_phasefield_1 = int64_t(phasefield->yStride());
    const int64_t _stride_phasefield_2 = int64_t(phasefield->zStride());
    const int64_t _stride_rho_a_0 = int64_t(rho_a->xStride());
    const int64_t _stride_rho_a_1 = int64_t(rho_a->yStride());
    const int64_t _stride_rho_a_2 = int64_t(rho_a->zStride());
    const int64_t _stride_rho_b_0 = int64_t(rho_b->xStride());
    const int64_t _stride_rho_b_1 = int64_t(rho_b->yStride());
    const int64_t _stride_rho_b_2 = int64_t(rho_b->zStride());
    const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
    const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
    internal_e4397140a5a3f258ae7331867e48f30b::colorgradientstreamsweepdoubleprecision_colorgradientstreamsweepdoubleprecision(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_a_tmp, _data_pdfs_b, _data_pdfs_b_tmp, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_a_tmp_0, _stride_pdfs_a_tmp_1, _stride_pdfs_a_tmp_2, _stride_pdfs_a_tmp_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_pdfs_b_tmp_0, _stride_pdfs_b_tmp_1, _stride_pdfs_b_tmp_2, _stride_pdfs_b_tmp_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3);
    pdfs_a->swapDataPointers(pdfs_a_tmp);
    pdfs_b->swapDataPointers(pdfs_b_tmp);

}


void ColorGradientStreamSweepDoublePrecision::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block)
{
   
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    auto rho_a = block->getData< field::GhostLayerField<double, 1> >(rho_aID);
    auto phasefield = block->getData< field::GhostLayerField<double, 1> >(phasefieldID);
    auto force_a = block->getData< field::GhostLayerField<double, 3> >(force_aID);
    auto pdfs_a = block->getData< field::GhostLayerField<double, 19> >(pdfs_aID);
    auto rho_b = block->getData< field::GhostLayerField<double, 1> >(rho_bID);
    auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
    auto pdfs_b = block->getData< field::GhostLayerField<double, 19> >(pdfs_bID);
    auto force_b = block->getData< field::GhostLayerField<double, 3> >(force_bID);
    field::GhostLayerField<double, 19> * pdfs_b_tmp;
    {
        if (cache_pdfs_b_.find(block) == cache_pdfs_b_.end())
        {
            pdfs_b_tmp = pdfs_b->cloneUninitialized();
            cache_pdfs_b_[block] = pdfs_b_tmp;
        }
        else
        {
            pdfs_b_tmp = cache_pdfs_b_[block];
        }
    }

    field::GhostLayerField<double, 19> * pdfs_a_tmp;
    {
        if (cache_pdfs_a_.find(block) == cache_pdfs_a_.end())
        {
            pdfs_a_tmp = pdfs_a->cloneUninitialized();
            cache_pdfs_a_[block] = pdfs_a_tmp;
        }
        else
        {
            pdfs_a_tmp = cache_pdfs_a_[block];
        }
    }

    
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    double * RESTRICT const _data_force_a = force_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT const _data_pdfs_a = pdfs_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_a_tmp->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_a_tmp->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_a_tmp->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a_tmp = pdfs_a_tmp->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a_tmp->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT const _data_pdfs_b = pdfs_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_b_tmp->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_b_tmp->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_b_tmp->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b_tmp = pdfs_b_tmp->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b_tmp->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT  _data_phasefield = phasefield->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT  _data_rho_a = rho_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT  _data_rho_b = rho_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT  _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
    const int64_t _size_force_a_0 = int64_t(int64_c(ci.xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
    const int64_t _size_force_a_1 = int64_t(int64_c(ci.ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
    const int64_t _size_force_a_2 = int64_t(int64_c(ci.zSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    const int64_t _stride_force_a_0 = int64_t(force_a->xStride());
    const int64_t _stride_force_a_1 = int64_t(force_a->yStride());
    const int64_t _stride_force_a_2 = int64_t(force_a->zStride());
    const int64_t _stride_force_a_3 = int64_t(1 * int64_t(force_a->fStride()));
    const int64_t _stride_force_b_0 = int64_t(force_b->xStride());
    const int64_t _stride_force_b_1 = int64_t(force_b->yStride());
    const int64_t _stride_force_b_2 = int64_t(force_b->zStride());
    const int64_t _stride_force_b_3 = int64_t(1 * int64_t(force_b->fStride()));
    const int64_t _stride_pdfs_a_0 = int64_t(pdfs_a->xStride());
    const int64_t _stride_pdfs_a_1 = int64_t(pdfs_a->yStride());
    const int64_t _stride_pdfs_a_2 = int64_t(pdfs_a->zStride());
    const int64_t _stride_pdfs_a_3 = int64_t(1 * int64_t(pdfs_a->fStride()));
    const int64_t _stride_pdfs_a_tmp_0 = int64_t(pdfs_a_tmp->xStride());
    const int64_t _stride_pdfs_a_tmp_1 = int64_t(pdfs_a_tmp->yStride());
    const int64_t _stride_pdfs_a_tmp_2 = int64_t(pdfs_a_tmp->zStride());
    const int64_t _stride_pdfs_a_tmp_3 = int64_t(1 * int64_t(pdfs_a_tmp->fStride()));
    const int64_t _stride_pdfs_b_0 = int64_t(pdfs_b->xStride());
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
    const int64_t _stride_pdfs_b_tmp_0 = int64_t(pdfs_b_tmp->xStride());
    const int64_t _stride_pdfs_b_tmp_1 = int64_t(pdfs_b_tmp->yStride());
    const int64_t _stride_pdfs_b_tmp_2 = int64_t(pdfs_b_tmp->zStride());
    const int64_t _stride_pdfs_b_tmp_3 = int64_t(1 * int64_t(pdfs_b_tmp->fStride()));
    const int64_t _stride_phasefield_0 = int64_t(phasefield->xStride());
    const int64_t _stride_phasefield_1 = int64_t(phasefield->yStride());
    const int64_t _stride_phasefield_2 = int64_t(phasefield->zStride());
    const int64_t _stride_rho_a_0 = int64_t(rho_a->xStride());
    const int64_t _stride_rho_a_1 = int64_t(rho_a->yStride());
    const int64_t _stride_rho_a_2 = int64_t(rho_a->zStride());
    const int64_t _stride_rho_b_0 = int64_t(rho_b->xStride());
    const int64_t _stride_rho_b_1 = int64_t(rho_b->yStride());
    const int64_t _stride_rho_b_2 = int64_t(rho_b->zStride());
    const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
    const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
    internal_e4397140a5a3f258ae7331867e48f30b::colorgradientstreamsweepdoubleprecision_colorgradientstreamsweepdoubleprecision(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_a_tmp, _data_pdfs_b, _data_pdfs_b_tmp, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_a_tmp_0, _stride_pdfs_a_tmp_1, _stride_pdfs_a_tmp_2, _stride_pdfs_a_tmp_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_pdfs_b_tmp_0, _stride_pdfs_b_tmp_1, _stride_pdfs_b_tmp_2, _stride_pdfs_b_tmp_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3);
    pdfs_a->swapDataPointers(pdfs_a_tmp);
    pdfs_b->swapDataPointers(pdfs_b_tmp);

}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
