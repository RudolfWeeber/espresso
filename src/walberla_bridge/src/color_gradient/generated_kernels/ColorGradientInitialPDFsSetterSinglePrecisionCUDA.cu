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
//! \\file ColorGradientInitialPDFsSetterSinglePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 17fc54c872bd8ceabf271a7e9e636c7c583f55af


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "ColorGradientInitialPDFsSetterSinglePrecisionCUDA.h"




#define FUNC_PREFIX __global__

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


namespace internal_colorgradientinitialpdfssettersingleprecisioncuda_colorgradientinitialpdfssettersingleprecisioncuda {
static FUNC_PREFIX __launch_bounds__(256) void colorgradientinitialpdfssettersingleprecisioncuda_colorgradientinitialpdfssettersingleprecisioncuda(float * RESTRICT const _data_force_a, float * RESTRICT const _data_force_b, float * RESTRICT  _data_pdfs_a, float * RESTRICT  _data_pdfs_b, float * RESTRICT  _data_phasefield, float * RESTRICT const _data_rho_a, float * RESTRICT const _data_rho_b, float * RESTRICT const _data_velocity, int64_t const _size_force_a_0, int64_t const _size_force_a_1, int64_t const _size_force_a_2, int64_t const _stride_force_a_0, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_0, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_0, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_b_0, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_phasefield_0, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_0, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_0, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3)
{
   if (blockDim.x*blockIdx.x + threadIdx.x < _size_force_a_0 && blockDim.y*blockIdx.y + threadIdx.y < _size_force_a_1 && blockDim.z*blockIdx.z + threadIdx.z < _size_force_a_2)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z;
      const float rho_a = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
      const float u_0_a = -0.5f*((1.0f) / (rho_a))*_data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2] + _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
      const float u_1_a = -0.5f*((1.0f) / (rho_a))*_data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3] + _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
      const float u_2_a = -0.5f*((1.0f) / (rho_a))*_data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3] + _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
      const float rho_b = _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float u_0_b = -0.5f*((1.0f) / (rho_b))*_data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2] + _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
      const float u_1_b = -0.5f*((1.0f) / (rho_b))*_data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3] + _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
      const float u_2_b = -0.5f*((1.0f) / (rho_b))*_data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3] + _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
      const float phi = (rho_a - rho_b)*((1.0f) / (rho_a + rho_b));
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2] = rho_a*-0.5f*(u_0_a*u_0_a) + rho_a*-0.5f*(u_1_a*u_1_a) + rho_a*-0.5f*(u_2_a*u_2_a) + rho_a*0.44444444444444459f;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3] = rho_a*u_1_a*0.16666666666666666f + rho_a*-0.083333333333333329f*(u_0_a*u_0_a) + rho_a*-0.083333333333333329f*(u_2_a*u_2_a) + rho_a*0.046296296296296259f + rho_a*0.16666666666666666f*(u_1_a*u_1_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3] = rho_a*u_1_a*-0.16666666666666666f + rho_a*-0.083333333333333329f*(u_0_a*u_0_a) + rho_a*-0.083333333333333329f*(u_2_a*u_2_a) + rho_a*0.046296296296296266f + rho_a*0.16666666666666666f*(u_1_a*u_1_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3] = rho_a*u_0_a*-0.16666666666666666f + rho_a*-0.083333333333333329f*(u_1_a*u_1_a) + rho_a*-0.083333333333333329f*(u_2_a*u_2_a) + rho_a*0.046296296296296266f + rho_a*0.16666666666666666f*(u_0_a*u_0_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3] = rho_a*u_0_a*0.16666666666666666f + rho_a*-0.083333333333333329f*(u_1_a*u_1_a) + rho_a*-0.083333333333333329f*(u_2_a*u_2_a) + rho_a*0.046296296296296259f + rho_a*0.16666666666666666f*(u_0_a*u_0_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3] = rho_a*u_2_a*0.16666666666666666f + rho_a*-0.083333333333333329f*(u_0_a*u_0_a) + rho_a*-0.083333333333333329f*(u_1_a*u_1_a) + rho_a*0.046296296296296259f + rho_a*0.16666666666666666f*(u_2_a*u_2_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3] = rho_a*u_2_a*-0.16666666666666666f + rho_a*-0.083333333333333329f*(u_0_a*u_0_a) + rho_a*-0.083333333333333329f*(u_1_a*u_1_a) + rho_a*0.046296296296296266f + rho_a*0.16666666666666666f*(u_2_a*u_2_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3] = rho_a*u_0_a*u_1_a*-0.25f + rho_a*u_0_a*-0.083333333333333329f + rho_a*u_1_a*0.083333333333333329f + rho_a*-0.041666666666666664f*(u_2_a*u_2_a) + rho_a*0.02314814814814815f + rho_a*0.083333333333333329f*(u_0_a*u_0_a) + rho_a*0.083333333333333329f*(u_1_a*u_1_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3] = rho_a*u_0_a*u_1_a*0.25f + rho_a*u_0_a*0.083333333333333329f + rho_a*u_1_a*0.083333333333333329f + rho_a*-0.041666666666666664f*(u_2_a*u_2_a) + rho_a*0.02314814814814815f + rho_a*0.083333333333333329f*(u_0_a*u_0_a) + rho_a*0.083333333333333329f*(u_1_a*u_1_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3] = rho_a*u_0_a*u_1_a*0.25f + rho_a*u_0_a*-0.083333333333333329f + rho_a*u_1_a*-0.083333333333333329f + rho_a*-0.041666666666666664f*(u_2_a*u_2_a) + rho_a*0.02314814814814815f + rho_a*0.083333333333333329f*(u_0_a*u_0_a) + rho_a*0.083333333333333329f*(u_1_a*u_1_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3] = rho_a*u_0_a*u_1_a*-0.25f + rho_a*u_0_a*0.083333333333333329f + rho_a*u_1_a*-0.083333333333333329f + rho_a*-0.041666666666666664f*(u_2_a*u_2_a) + rho_a*0.02314814814814815f + rho_a*0.083333333333333329f*(u_0_a*u_0_a) + rho_a*0.083333333333333329f*(u_1_a*u_1_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3] = rho_a*u_1_a*u_2_a*0.25f + rho_a*u_1_a*0.083333333333333329f + rho_a*u_2_a*0.083333333333333329f + rho_a*-0.041666666666666664f*(u_0_a*u_0_a) + rho_a*0.02314814814814815f + rho_a*0.083333333333333329f*(u_1_a*u_1_a) + rho_a*0.083333333333333329f*(u_2_a*u_2_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3] = rho_a*u_1_a*u_2_a*-0.25f + rho_a*u_1_a*-0.083333333333333329f + rho_a*u_2_a*0.083333333333333329f + rho_a*-0.041666666666666664f*(u_0_a*u_0_a) + rho_a*0.02314814814814815f + rho_a*0.083333333333333329f*(u_1_a*u_1_a) + rho_a*0.083333333333333329f*(u_2_a*u_2_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3] = rho_a*u_0_a*u_2_a*-0.25f + rho_a*u_0_a*-0.083333333333333329f + rho_a*u_2_a*0.083333333333333329f + rho_a*-0.041666666666666664f*(u_1_a*u_1_a) + rho_a*0.02314814814814815f + rho_a*0.083333333333333329f*(u_0_a*u_0_a) + rho_a*0.083333333333333329f*(u_2_a*u_2_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3] = rho_a*u_0_a*u_2_a*0.25f + rho_a*u_0_a*0.083333333333333329f + rho_a*u_2_a*0.083333333333333329f + rho_a*-0.041666666666666664f*(u_1_a*u_1_a) + rho_a*0.02314814814814815f + rho_a*0.083333333333333329f*(u_0_a*u_0_a) + rho_a*0.083333333333333329f*(u_2_a*u_2_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3] = rho_a*u_1_a*u_2_a*-0.25f + rho_a*u_1_a*0.083333333333333329f + rho_a*u_2_a*-0.083333333333333329f + rho_a*-0.041666666666666664f*(u_0_a*u_0_a) + rho_a*0.02314814814814815f + rho_a*0.083333333333333329f*(u_1_a*u_1_a) + rho_a*0.083333333333333329f*(u_2_a*u_2_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3] = rho_a*u_1_a*u_2_a*0.25f + rho_a*u_1_a*-0.083333333333333329f + rho_a*u_2_a*-0.083333333333333329f + rho_a*-0.041666666666666664f*(u_0_a*u_0_a) + rho_a*0.02314814814814815f + rho_a*0.083333333333333329f*(u_1_a*u_1_a) + rho_a*0.083333333333333329f*(u_2_a*u_2_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3] = rho_a*u_0_a*u_2_a*0.25f + rho_a*u_0_a*-0.083333333333333329f + rho_a*u_2_a*-0.083333333333333329f + rho_a*-0.041666666666666664f*(u_1_a*u_1_a) + rho_a*0.02314814814814815f + rho_a*0.083333333333333329f*(u_0_a*u_0_a) + rho_a*0.083333333333333329f*(u_2_a*u_2_a);
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3] = rho_a*u_0_a*u_2_a*-0.25f + rho_a*u_0_a*0.083333333333333329f + rho_a*u_2_a*-0.083333333333333329f + rho_a*-0.041666666666666664f*(u_1_a*u_1_a) + rho_a*0.02314814814814815f + rho_a*0.083333333333333329f*(u_0_a*u_0_a) + rho_a*0.083333333333333329f*(u_2_a*u_2_a);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2] = rho_b*-0.5f*(u_0_b*u_0_b) + rho_b*-0.5f*(u_1_b*u_1_b) + rho_b*-0.5f*(u_2_b*u_2_b) + rho_b*0.44444444444444459f;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3] = rho_b*u_1_b*0.16666666666666666f + rho_b*-0.083333333333333329f*(u_0_b*u_0_b) + rho_b*-0.083333333333333329f*(u_2_b*u_2_b) + rho_b*0.046296296296296259f + rho_b*0.16666666666666666f*(u_1_b*u_1_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3] = rho_b*u_1_b*-0.16666666666666666f + rho_b*-0.083333333333333329f*(u_0_b*u_0_b) + rho_b*-0.083333333333333329f*(u_2_b*u_2_b) + rho_b*0.046296296296296266f + rho_b*0.16666666666666666f*(u_1_b*u_1_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3] = rho_b*u_0_b*-0.16666666666666666f + rho_b*-0.083333333333333329f*(u_1_b*u_1_b) + rho_b*-0.083333333333333329f*(u_2_b*u_2_b) + rho_b*0.046296296296296266f + rho_b*0.16666666666666666f*(u_0_b*u_0_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3] = rho_b*u_0_b*0.16666666666666666f + rho_b*-0.083333333333333329f*(u_1_b*u_1_b) + rho_b*-0.083333333333333329f*(u_2_b*u_2_b) + rho_b*0.046296296296296259f + rho_b*0.16666666666666666f*(u_0_b*u_0_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3] = rho_b*u_2_b*0.16666666666666666f + rho_b*-0.083333333333333329f*(u_0_b*u_0_b) + rho_b*-0.083333333333333329f*(u_1_b*u_1_b) + rho_b*0.046296296296296259f + rho_b*0.16666666666666666f*(u_2_b*u_2_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3] = rho_b*u_2_b*-0.16666666666666666f + rho_b*-0.083333333333333329f*(u_0_b*u_0_b) + rho_b*-0.083333333333333329f*(u_1_b*u_1_b) + rho_b*0.046296296296296266f + rho_b*0.16666666666666666f*(u_2_b*u_2_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3] = rho_b*u_0_b*u_1_b*-0.25f + rho_b*u_0_b*-0.083333333333333329f + rho_b*u_1_b*0.083333333333333329f + rho_b*-0.041666666666666664f*(u_2_b*u_2_b) + rho_b*0.02314814814814815f + rho_b*0.083333333333333329f*(u_0_b*u_0_b) + rho_b*0.083333333333333329f*(u_1_b*u_1_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3] = rho_b*u_0_b*u_1_b*0.25f + rho_b*u_0_b*0.083333333333333329f + rho_b*u_1_b*0.083333333333333329f + rho_b*-0.041666666666666664f*(u_2_b*u_2_b) + rho_b*0.02314814814814815f + rho_b*0.083333333333333329f*(u_0_b*u_0_b) + rho_b*0.083333333333333329f*(u_1_b*u_1_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3] = rho_b*u_0_b*u_1_b*0.25f + rho_b*u_0_b*-0.083333333333333329f + rho_b*u_1_b*-0.083333333333333329f + rho_b*-0.041666666666666664f*(u_2_b*u_2_b) + rho_b*0.02314814814814815f + rho_b*0.083333333333333329f*(u_0_b*u_0_b) + rho_b*0.083333333333333329f*(u_1_b*u_1_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3] = rho_b*u_0_b*u_1_b*-0.25f + rho_b*u_0_b*0.083333333333333329f + rho_b*u_1_b*-0.083333333333333329f + rho_b*-0.041666666666666664f*(u_2_b*u_2_b) + rho_b*0.02314814814814815f + rho_b*0.083333333333333329f*(u_0_b*u_0_b) + rho_b*0.083333333333333329f*(u_1_b*u_1_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3] = rho_b*u_1_b*u_2_b*0.25f + rho_b*u_1_b*0.083333333333333329f + rho_b*u_2_b*0.083333333333333329f + rho_b*-0.041666666666666664f*(u_0_b*u_0_b) + rho_b*0.02314814814814815f + rho_b*0.083333333333333329f*(u_1_b*u_1_b) + rho_b*0.083333333333333329f*(u_2_b*u_2_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3] = rho_b*u_1_b*u_2_b*-0.25f + rho_b*u_1_b*-0.083333333333333329f + rho_b*u_2_b*0.083333333333333329f + rho_b*-0.041666666666666664f*(u_0_b*u_0_b) + rho_b*0.02314814814814815f + rho_b*0.083333333333333329f*(u_1_b*u_1_b) + rho_b*0.083333333333333329f*(u_2_b*u_2_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3] = rho_b*u_0_b*u_2_b*-0.25f + rho_b*u_0_b*-0.083333333333333329f + rho_b*u_2_b*0.083333333333333329f + rho_b*-0.041666666666666664f*(u_1_b*u_1_b) + rho_b*0.02314814814814815f + rho_b*0.083333333333333329f*(u_0_b*u_0_b) + rho_b*0.083333333333333329f*(u_2_b*u_2_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3] = rho_b*u_0_b*u_2_b*0.25f + rho_b*u_0_b*0.083333333333333329f + rho_b*u_2_b*0.083333333333333329f + rho_b*-0.041666666666666664f*(u_1_b*u_1_b) + rho_b*0.02314814814814815f + rho_b*0.083333333333333329f*(u_0_b*u_0_b) + rho_b*0.083333333333333329f*(u_2_b*u_2_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3] = rho_b*u_1_b*u_2_b*-0.25f + rho_b*u_1_b*0.083333333333333329f + rho_b*u_2_b*-0.083333333333333329f + rho_b*-0.041666666666666664f*(u_0_b*u_0_b) + rho_b*0.02314814814814815f + rho_b*0.083333333333333329f*(u_1_b*u_1_b) + rho_b*0.083333333333333329f*(u_2_b*u_2_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3] = rho_b*u_1_b*u_2_b*0.25f + rho_b*u_1_b*-0.083333333333333329f + rho_b*u_2_b*-0.083333333333333329f + rho_b*-0.041666666666666664f*(u_0_b*u_0_b) + rho_b*0.02314814814814815f + rho_b*0.083333333333333329f*(u_1_b*u_1_b) + rho_b*0.083333333333333329f*(u_2_b*u_2_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3] = rho_b*u_0_b*u_2_b*0.25f + rho_b*u_0_b*-0.083333333333333329f + rho_b*u_2_b*-0.083333333333333329f + rho_b*-0.041666666666666664f*(u_1_b*u_1_b) + rho_b*0.02314814814814815f + rho_b*0.083333333333333329f*(u_0_b*u_0_b) + rho_b*0.083333333333333329f*(u_2_b*u_2_b);
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3] = rho_b*u_0_b*u_2_b*-0.25f + rho_b*u_0_b*0.083333333333333329f + rho_b*u_2_b*-0.083333333333333329f + rho_b*-0.041666666666666664f*(u_1_b*u_1_b) + rho_b*0.02314814814814815f + rho_b*0.083333333333333329f*(u_0_b*u_0_b) + rho_b*0.083333333333333329f*(u_2_b*u_2_b);
      _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] = phi;
   } 
}
}


void ColorGradientInitialPDFsSetterSinglePrecisionCUDA::run(IBlock * block, gpuStream_t stream)
{
   
    auto pdfs_a = block->getData< gpu::GPUField<float> >(pdfs_aID);
    auto velocity = block->getData< gpu::GPUField<float> >(velocityID);
    auto phasefield = block->getData< gpu::GPUField<float> >(phasefieldID);
    auto force_b = block->getData< gpu::GPUField<float> >(force_bID);
    auto pdfs_b = block->getData< gpu::GPUField<float> >(pdfs_bID);
    auto rho_a = block->getData< gpu::GPUField<float> >(rho_aID);
    auto force_a = block->getData< gpu::GPUField<float> >(force_aID);
    auto rho_b = block->getData< gpu::GPUField<float> >(rho_bID);

    
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force_a->nrOfGhostLayers()))
    float * RESTRICT const _data_force_a = force_a->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force_b->nrOfGhostLayers()))
    float * RESTRICT const _data_force_b = force_b->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs_a->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs_b->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(phasefield->nrOfGhostLayers()))
    float * RESTRICT  _data_phasefield = phasefield->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_a->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_a = rho_a->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_b->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_b = rho_b->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()))
    float * RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(force_a->xSize()) + 0))
    const int64_t _size_force_a_0 = int64_t(int64_c(force_a->xSize()) + 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(force_a->ySize()) + 0))
    const int64_t _size_force_a_1 = int64_t(int64_c(force_a->ySize()) + 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(force_a->zSize()) + 0))
    const int64_t _size_force_a_2 = int64_t(int64_c(force_a->zSize()) + 0);
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
    const int64_t _stride_pdfs_b_0 = int64_t(pdfs_b->xStride());
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
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
    dim3 _block(uint32_c(((128 < _size_force_a_0) ? 128 : _size_force_a_0)), uint32_c(((1024 < ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))) ? 1024 : ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))), uint32_c(((64 < ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))))) ? 64 : ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))))));
    dim3 _grid(uint32_c(( (_size_force_a_0) % (((128 < _size_force_a_0) ? 128 : _size_force_a_0)) == 0 ? (int64_t)(_size_force_a_0) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)) : ( (int64_t)(_size_force_a_0) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)) ) +1 )), uint32_c(( (_size_force_a_1) % (((1024 < ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))) ? 1024 : ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))) == 0 ? (int64_t)(_size_force_a_1) / (int64_t)(((1024 < ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))) ? 1024 : ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))) : ( (int64_t)(_size_force_a_1) / (int64_t)(((1024 < ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))) ? 1024 : ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))) ) +1 )), uint32_c(( (_size_force_a_2) % (((64 < ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))))) ? 64 : ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))))) == 0 ? (int64_t)(_size_force_a_2) / (int64_t)(((64 < ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))))) ? 64 : ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))))) : ( (int64_t)(_size_force_a_2) / (int64_t)(((64 < ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))))) ? 64 : ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))))) ) +1 )));
    internal_colorgradientinitialpdfssettersingleprecisioncuda_colorgradientinitialpdfssettersingleprecisioncuda::colorgradientinitialpdfssettersingleprecisioncuda_colorgradientinitialpdfssettersingleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3);
    
}


void ColorGradientInitialPDFsSetterSinglePrecisionCUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block, gpuStream_t stream)
{
   
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    auto pdfs_a = block->getData< gpu::GPUField<float> >(pdfs_aID);
    auto velocity = block->getData< gpu::GPUField<float> >(velocityID);
    auto phasefield = block->getData< gpu::GPUField<float> >(phasefieldID);
    auto force_b = block->getData< gpu::GPUField<float> >(force_bID);
    auto pdfs_b = block->getData< gpu::GPUField<float> >(pdfs_bID);
    auto rho_a = block->getData< gpu::GPUField<float> >(rho_aID);
    auto force_a = block->getData< gpu::GPUField<float> >(force_aID);
    auto rho_b = block->getData< gpu::GPUField<float> >(rho_bID);

    
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force_a->nrOfGhostLayers()))
    float * RESTRICT const _data_force_a = force_a->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force_b->nrOfGhostLayers()))
    float * RESTRICT const _data_force_b = force_b->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_a->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_b->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(phasefield->nrOfGhostLayers()))
    float * RESTRICT  _data_phasefield = phasefield->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_a->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_a = rho_a->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_b->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_b = rho_b->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()))
    float * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_force_a_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_force_a_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_force_a_2 = int64_t(int64_c(ci.zSize()) + 0);
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
    const int64_t _stride_pdfs_b_0 = int64_t(pdfs_b->xStride());
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
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
    dim3 _block(uint32_c(((128 < _size_force_a_0) ? 128 : _size_force_a_0)), uint32_c(((1024 < ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))) ? 1024 : ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))), uint32_c(((64 < ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))))) ? 64 : ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))))));
    dim3 _grid(uint32_c(( (_size_force_a_0) % (((128 < _size_force_a_0) ? 128 : _size_force_a_0)) == 0 ? (int64_t)(_size_force_a_0) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)) : ( (int64_t)(_size_force_a_0) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)) ) +1 )), uint32_c(( (_size_force_a_1) % (((1024 < ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))) ? 1024 : ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))) == 0 ? (int64_t)(_size_force_a_1) / (int64_t)(((1024 < ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))) ? 1024 : ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))) : ( (int64_t)(_size_force_a_1) / (int64_t)(((1024 < ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))) ? 1024 : ((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))) ) +1 )), uint32_c(( (_size_force_a_2) % (((64 < ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))))) ? 64 : ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))))) == 0 ? (int64_t)(_size_force_a_2) / (int64_t)(((64 < ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))))) ? 64 : ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))))) : ( (int64_t)(_size_force_a_2) / (int64_t)(((64 < ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))))))) ? 64 : ((_size_force_a_2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))) ? _size_force_a_2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)*((_size_force_a_1 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0)))) ? _size_force_a_1 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0) ? 128 : _size_force_a_0))))))))) ) +1 )));
    internal_colorgradientinitialpdfssettersingleprecisioncuda_colorgradientinitialpdfssettersingleprecisioncuda::colorgradientinitialpdfssettersingleprecisioncuda_colorgradientinitialpdfssettersingleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
