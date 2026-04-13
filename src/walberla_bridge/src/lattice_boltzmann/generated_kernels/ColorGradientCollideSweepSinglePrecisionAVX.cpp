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
//! \\file ColorGradientCollideSweepSinglePrecisionAVX.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 3247aa7395049ca5bfb69d34d55e45db19fa439c


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "ColorGradientCollideSweepSinglePrecisionAVX.h"


#include <immintrin.h>



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


namespace internal_196de4e60a60b4c0f3c29ed6234bf8bc {
static FUNC_PREFIX void colorgradientcollidesweepsingleprecisionavx_colorgradientcollidesweepsingleprecisionavx(float * RESTRICT const _data_color_gradient, float * RESTRICT const _data_force_a, float * RESTRICT const _data_force_b, float * RESTRICT  _data_pdfs_a, float * RESTRICT  _data_pdfs_b, float * RESTRICT const _data_phasefield, float * RESTRICT const _data_rho_a, float * RESTRICT const _data_rho_b, float * RESTRICT const _data_velocity, int64_t const _size_color_gradient_0, int64_t const _size_color_gradient_1, int64_t const _size_color_gradient_2, int64_t const _stride_color_gradient_1, int64_t const _stride_color_gradient_2, int64_t const _stride_color_gradient_3, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3, float beta, float omega_even_a, float omega_even_b, float omega_odd_a, float omega_odd_b, float omega_shear_a, float omega_shear_b, float sigma)
{
#ifdef _OPENMP
   #pragma omp parallel
#endif
   {
      const float xi_3 = omega_even_a*0.5f;
      const float xi_33 = omega_odd_a*0.041666666666666664f;
      const float xi_45 = omega_shear_a*0.125f;
      const float xi_85 = omega_even_b*0.5f;
      const float xi_115 = omega_odd_b*0.041666666666666664f;
      const float xi_127 = omega_shear_b*0.125f;
      const float xi_190 = omega_shear_a*omega_shear_b*((1.0f) / (omega_shear_a + omega_shear_b));
      const float xi_191 = xi_190*2.0f;
      const float xi_192 = xi_190*8.0f;
      const float xi_193 = omega_shear_a*-4.0f + xi_192;
      const float xi_195 = omega_shear_b*-4.0f + xi_192;
      const float xi_204 = omega_odd_a*0.5f;
      const float xi_234 = omega_odd_a*0.25f;
      const float xi_238 = omega_shear_a*0.25f;
      const float xi_307 = omega_odd_b*0.5f;
      const float xi_331 = omega_odd_b*0.25f;
      const float xi_335 = omega_shear_b*0.25f;
      const float rr_0_a_collide = 0.0f;
      const float xi_5 = rr_0_a_collide*0.25f;
      const float rr_0_b_collide = 0.0f;
      const float xi_87 = rr_0_b_collide*0.25f;
#ifdef _OPENMP
      #pragma omp for schedule(static)
#endif
      for (int64_t ctr_2 = 0; ctr_2 < _size_color_gradient_2; ctr_2 += 1)
      {
         for (int64_t ctr_1 = 0; ctr_1 < _size_color_gradient_1; ctr_1 += 1)
         {
            {
               for (int64_t ctr_0 = 0; ctr_0 < (int64_t)((_size_color_gradient_0) / (8)) * (8); ctr_0 += 8)
               {
                  const __m256 xi_185 = _mm256_mul_ps(_mm256_load_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]),_mm256_load_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]));
                  const __m256 xi_186 = _mm256_mul_ps(_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]),_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]));
                  const __m256 xi_187 = _mm256_mul_ps(_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]),_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]));
                  const __m256 xi_188 = _mm256_add_ps(_mm256_add_ps(xi_185,xi_186),xi_187);
                  const __m256 xi_189 = _mm256_sqrt_ps(xi_188);
                  const __m256 xi_194 = _mm256_mul_ps(_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]));
                  const __m256 xi_196 = _mm256_mul_ps(_mm256_mul_ps(xi_189,_mm256_set_ps(sigma,sigma,sigma,sigma,sigma,sigma,sigma,sigma)),_mm256_blendv_ps(_mm256_blendv_ps(_mm256_blendv_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_194,_mm256_set_ps(xi_195,xi_195,xi_195,xi_195,xi_195,xi_195,xi_195,xi_195)),_mm256_mul_ps(_mm256_set_ps(xi_195,xi_195,xi_195,xi_195,xi_195,xi_195,xi_195,xi_195),_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))),_mm256_set_ps(xi_191,xi_191,xi_191,xi_191,xi_191,xi_191,xi_191,xi_191)),_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_194,_mm256_set_ps(xi_193,xi_193,xi_193,xi_193,xi_193,xi_193,xi_193,xi_193)),_mm256_mul_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f),_mm256_set_ps(xi_193,xi_193,xi_193,xi_193,xi_193,xi_193,xi_193,xi_193)),_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))),_mm256_set_ps(xi_191,xi_191,xi_191,xi_191,xi_191,xi_191,xi_191,xi_191)),_mm256_cmp_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NGE_UQ)),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),_mm256_cmp_ps(_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f),_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NLE_UQ)),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),_mm256_cmp_ps(_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f),_mm256_load_ps(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NGE_UQ)));
                  const __m256 xi_197 = _mm256_cmp_ps(xi_189,_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_CMP_NLE_UQ);
                  const __m256 xi_198 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_196,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f)),xi_197);
                  const __m256 xi_212 = _mm256_div_ps(_mm256_set_ps(1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f),xi_188);
                  const __m256 xi_213 = _mm256_mul_ps(xi_212,_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f));
                  const __m256 xi_214 = _mm256_mul_ps(xi_196,_mm256_set_ps(1.125f,1.125f,1.125f,1.125f,1.125f,1.125f,1.125f,1.125f));
                  const __m256 xi_215 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_186,xi_213),_mm256_set_ps(-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f))),xi_197);
                  const __m256 xi_225 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_185,xi_213),_mm256_set_ps(-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f))),xi_197);
                  const __m256 xi_232 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_187,xi_213),_mm256_set_ps(-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f,-0.018518518518518517f))),xi_197);
                  const __m256 xi_241 = _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f),_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0])),_mm256_load_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]));
                  const __m256 xi_242 = _mm256_mul_ps(xi_212,_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f));
                  const __m256 xi_243 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_mul_ps(xi_241,xi_241)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_197);
                  const __m256 xi_247 = _mm256_add_ps(_mm256_load_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]),_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]));
                  const __m256 xi_248 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_mul_ps(xi_247,xi_247)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_197);
                  const __m256 xi_249 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_mul_ps(xi_247,xi_247)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_197);
                  const __m256 xi_251 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_mul_ps(xi_241,xi_241)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_197);
                  const __m256 xi_254 = _mm256_add_ps(_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]),_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]));
                  const __m256 xi_255 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_mul_ps(xi_254,xi_254)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_197);
                  const __m256 xi_261 = _mm256_mul_ps(_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f),_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]));
                  const __m256 xi_262 = _mm256_add_ps(xi_261,_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]));
                  const __m256 xi_263 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_mul_ps(xi_262,xi_262)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_197);
                  const __m256 xi_268 = _mm256_add_ps(xi_261,_mm256_load_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]));
                  const __m256 xi_269 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_mul_ps(xi_268,xi_268)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_197);
                  const __m256 xi_273 = _mm256_add_ps(_mm256_load_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]),_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]));
                  const __m256 xi_274 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_mul_ps(xi_273,xi_273)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_197);
                  const __m256 xi_275 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_mul_ps(xi_262,xi_262)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_197);
                  const __m256 xi_277 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_mul_ps(xi_254,xi_254)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_197);
                  const __m256 xi_278 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_mul_ps(xi_273,xi_273)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_197);
                  const __m256 xi_280 = _mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_214,_mm256_add_ps(_mm256_mul_ps(xi_242,_mm256_mul_ps(xi_268,xi_268)),_mm256_set_ps(-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f,-0.037037037037037035f))),xi_197);
                  const __m256 xia_1_collide = _mm256_load_ps(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0]);
                  const __m256 xia_2_collide = _mm256_loadu_ps(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3 + ctr_0]);
                  const __m256 xi_4 = _mm256_mul_ps(xia_2_collide,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_7 = _mm256_mul_ps(xi_4,_mm256_set_ps(omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a));
                  const __m256 xi_37 = _mm256_mul_ps(xia_2_collide,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_38 = _mm256_mul_ps(xia_2_collide,_mm256_set_ps(xi_33,xi_33,xi_33,xi_33,xi_33,xi_33,xi_33,xi_33));
                  const __m256 xi_39 = _mm256_add_ps(_mm256_mul_ps(xi_37,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_38);
                  const __m256 xi_41 = _mm256_mul_ps(xia_2_collide,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_46 = _mm256_mul_ps(xia_2_collide,_mm256_set_ps(xi_45,xi_45,xi_45,xi_45,xi_45,xi_45,xi_45,xi_45));
                  const __m256 xi_53 = _mm256_add_ps(_mm256_mul_ps(xi_38,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_37);
                  const __m256 xia_3_collide = _mm256_loadu_ps(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0]);
                  const __m256 xia_4_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_220 = _mm256_mul_ps(xia_4_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xia_5_collide = _mm256_load_ps(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + ctr_0]);
                  const __m256 xi_17 = _mm256_mul_ps(xia_5_collide,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_18 = _mm256_mul_ps(xi_17,_mm256_set_ps(omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a));
                  const __m256 xi_32 = _mm256_mul_ps(xia_5_collide,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_34 = _mm256_mul_ps(xia_5_collide,_mm256_set_ps(xi_33,xi_33,xi_33,xi_33,xi_33,xi_33,xi_33,xi_33));
                  const __m256 xi_35 = _mm256_add_ps(_mm256_mul_ps(xi_34,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_32);
                  const __m256 xi_51 = _mm256_add_ps(_mm256_mul_ps(xi_32,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_34);
                  const __m256 xia_6_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_7_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_8_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_206 = _mm256_mul_ps(xia_8_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xia_9_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_221 = _mm256_add_ps(xi_220,xia_9_collide);
                  const __m256 xia_10_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_11_collide = _mm256_load_ps(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]);
                  const __m256 xia_12_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_13_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_177 = _mm256_add_ps(xia_13_collide,xia_7_collide);
                  const __m256 xia_14_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_15_collide = _mm256_load_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_178 = _mm256_add_ps(_mm256_add_ps(xi_177,xia_12_collide),xia_15_collide);
                  const __m256 xi_201 = _mm256_mul_ps(xia_15_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_202 = _mm256_add_ps(xi_201,xia_12_collide);
                  const __m256 xia_16_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_17_collide = _mm256_loadu_ps(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0]);
                  const __m256 xia_18_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_19_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_170 = _mm256_add_ps(xia_18_collide,xia_19_collide);
                  const __m256 xia_20_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_21_collide = _mm256_load_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_171 = _mm256_add_ps(_mm256_add_ps(xi_170,xia_21_collide),xia_8_collide);
                  const __m256 xi_207 = _mm256_add_ps(xi_206,xia_21_collide);
                  const __m256 xia_22_collide = _mm256_loadu_ps(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3 + ctr_0]);
                  const __m256 xi_26 = _mm256_mul_ps(xia_22_collide,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_28 = _mm256_mul_ps(xi_26,_mm256_set_ps(omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a));
                  const __m256 xi_55 = _mm256_mul_ps(xia_22_collide,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_56 = _mm256_mul_ps(xia_22_collide,_mm256_set_ps(xi_33,xi_33,xi_33,xi_33,xi_33,xi_33,xi_33,xi_33));
                  const __m256 xi_57 = _mm256_add_ps(_mm256_mul_ps(xi_56,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_55);
                  const __m256 xi_70 = _mm256_add_ps(_mm256_mul_ps(xi_55,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_56);
                  const __m256 xia_23_collide = _mm256_load_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0]);
                  const __m256 xia_24_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xia_25_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_174 = _mm256_add_ps(xia_25_collide,xia_6_collide);
                  const __m256 xi_175 = _mm256_add_ps(_mm256_add_ps(xi_174,xia_4_collide),xia_9_collide);
                  const __m256 xia_26_collide = _mm256_loadu_ps(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0]);
                  const __m256 xi_6 = _mm256_mul_ps(xia_2_collide,_mm256_set_ps(xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5));
                  const __m256 xi_19 = _mm256_mul_ps(xia_5_collide,_mm256_set_ps(xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5));
                  const __m256 xi_27 = _mm256_mul_ps(xia_22_collide,_mm256_set_ps(xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5,xi_5));
                  const __m256 rho_a_collide = xia_11_collide;
                  const __m256 xi_167 = _mm256_mul_ps(rho_a_collide,_mm256_set_ps(-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f));
                  const __m256 xi_180 = _mm256_mul_ps(rho_a_collide,_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f));
                  const __m256 xi_181 = _mm256_add_ps(xi_178,xi_180);
                  const __m256 xi_199 = _mm256_mul_ps(rho_a_collide,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f));
                  const __m256 u_0_a_collide = xia_1_collide;
                  const __m256 xi_0 = _mm256_mul_ps(u_0_a_collide,xia_5_collide);
                  const __m256 xi_13 = _mm256_mul_ps(xi_0,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_14 = _mm256_add_ps(_mm256_mul_ps(xi_0,_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f)),_mm256_mul_ps(xi_13,_mm256_set_ps(omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a)));
                  const __m256 xi_20 = _mm256_mul_ps(xi_0,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f));
                  const __m256 xi_42 = _mm256_mul_ps(u_0_a_collide,xi_41);
                  const __m256 xi_47 = _mm256_mul_ps(u_0_a_collide,xi_46);
                  const __m256 xi_58 = _mm256_mul_ps(xi_13,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_61 = _mm256_mul_ps(_mm256_mul_ps(xi_0,_mm256_set_ps(0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f)),_mm256_set_ps(omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256 xi_73 = _mm256_mul_ps(u_0_a_collide,xia_22_collide);
                  const __m256 xi_74 = _mm256_mul_ps(xi_73,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_77 = _mm256_mul_ps(xi_73,_mm256_set_ps(xi_45,xi_45,xi_45,xi_45,xi_45,xi_45,xi_45,xi_45));
                  const __m256 xi_166 = _mm256_mul_ps(u_0_a_collide,u_0_a_collide);
                  const __m256 xi_173 = _mm256_mul_ps(_mm256_mul_ps(rho_a_collide,xi_166),_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f));
                  const __m256 xi_182 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_175,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_181,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xia_10_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xia_26_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(rho_a_collide,xi_166)),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a));
                  const __m256 xi_217 = _mm256_mul_ps(u_0_a_collide,xi_199);
                  const __m256 xi_218 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_7_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_202),xi_217),xia_13_collide);
                  const __m256 xi_219 = _mm256_mul_ps(xi_218,_mm256_set_ps(xi_204,xi_204,xi_204,xi_204,xi_204,xi_204,xi_204,xi_204));
                  const __m256 xi_222 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_25_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_217),xi_221),xia_6_collide);
                  const __m256 xi_223 = _mm256_mul_ps(xi_222,_mm256_set_ps(xi_204,xi_204,xi_204,xi_204,xi_204,xi_204,xi_204,xi_204));
                  const __m256 xi_235 = _mm256_mul_ps(xi_218,_mm256_set_ps(xi_234,xi_234,xi_234,xi_234,xi_234,xi_234,xi_234,xi_234));
                  const __m256 xi_236 = _mm256_mul_ps(xi_235,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_264 = _mm256_mul_ps(xi_222,_mm256_set_ps(xi_234,xi_234,xi_234,xi_234,xi_234,xi_234,xi_234,xi_234));
                  const __m256 xi_265 = _mm256_mul_ps(xi_264,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 u_1_a_collide = xia_3_collide;
                  const __m256 xi_1 = _mm256_mul_ps(u_1_a_collide,xia_2_collide);
                  const __m256 xi_8 = _mm256_mul_ps(xi_1,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f));
                  const __m256 xi_21 = _mm256_mul_ps(xi_1,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_22 = _mm256_mul_ps(xi_1,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_23 = _mm256_mul_ps(xi_22,_mm256_set_ps(omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256 xi_24 = _mm256_add_ps(_mm256_mul_ps(xi_21,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_23);
                  const __m256 xi_30 = _mm256_add_ps(xi_14,xi_24);
                  const __m256 xi_43 = _mm256_mul_ps(u_1_a_collide,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_44 = _mm256_mul_ps(xi_43,xia_5_collide);
                  const __m256 xi_48 = _mm256_mul_ps(u_1_a_collide,_mm256_set_ps(xi_45,xi_45,xi_45,xi_45,xi_45,xi_45,xi_45,xi_45));
                  const __m256 xi_49 = _mm256_mul_ps(xi_48,xia_5_collide);
                  const __m256 xi_50 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_47,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_49,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_42),xi_44);
                  const __m256 xi_52 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_42,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_44,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_47),xi_49);
                  const __m256 xi_59 = _mm256_mul_ps(xi_23,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_63 = _mm256_mul_ps(xi_43,xia_22_collide);
                  const __m256 xi_65 = _mm256_mul_ps(xi_48,xia_22_collide);
                  const __m256 xi_71 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(u_1_a_collide,xia_2_collide),_mm256_set_ps(-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f)),_mm256_set_ps(omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256 xi_168 = _mm256_mul_ps(u_1_a_collide,u_1_a_collide);
                  const __m256 xi_169 = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(rho_a_collide,xi_168),_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f)),xi_167);
                  const __m256 xi_183 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_171,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_181,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xia_14_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xia_16_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(rho_a_collide,xi_168)),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a));
                  const __m256 xi_200 = _mm256_mul_ps(u_1_a_collide,xi_199);
                  const __m256 xi_203 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_13_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_200),xi_202),xia_7_collide);
                  const __m256 xi_205 = _mm256_mul_ps(xi_203,_mm256_set_ps(xi_204,xi_204,xi_204,xi_204,xi_204,xi_204,xi_204,xi_204));
                  const __m256 xi_208 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_19_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_200),xi_207),xia_18_collide);
                  const __m256 xi_209 = _mm256_mul_ps(xi_208,_mm256_set_ps(xi_204,xi_204,xi_204,xi_204,xi_204,xi_204,xi_204,xi_204));
                  const __m256 xi_237 = _mm256_mul_ps(rho_a_collide,u_1_a_collide);
                  const __m256 xi_239 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_12_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(u_0_a_collide,xi_237)),xi_177),xi_201),_mm256_set_ps(xi_238,xi_238,xi_238,xi_238,xi_238,xi_238,xi_238,xi_238));
                  const __m256 xi_240 = _mm256_mul_ps(xi_239,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_245 = _mm256_mul_ps(xi_203,_mm256_set_ps(xi_234,xi_234,xi_234,xi_234,xi_234,xi_234,xi_234,xi_234));
                  const __m256 xi_252 = _mm256_mul_ps(xi_208,_mm256_set_ps(xi_234,xi_234,xi_234,xi_234,xi_234,xi_234,xi_234,xi_234));
                  const __m256 xi_259 = _mm256_mul_ps(xi_252,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 u_2_a_collide = xia_17_collide;
                  const __m256 xi_2 = _mm256_mul_ps(u_2_a_collide,xia_22_collide);
                  const __m256 xi_9 = _mm256_mul_ps(xi_2,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_10 = _mm256_mul_ps(xi_2,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_11 = _mm256_mul_ps(xi_10,_mm256_set_ps(omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256 xi_12 = _mm256_add_ps(_mm256_mul_ps(xi_9,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_11);
                  const __m256 xi_15 = _mm256_add_ps(xi_12,xi_14);
                  const __m256 xi_16 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_8,_mm256_set_ps(omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a)),_mm256_mul_ps(_mm256_mul_ps(xi_1,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f)),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),xi_15),xi_8);
                  const __m256 xi_25 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_20,_mm256_set_ps(omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a)),_mm256_mul_ps(_mm256_mul_ps(xi_0,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f)),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),xi_12),xi_20),xi_24);
                  const __m256 xi_29 = _mm256_mul_ps(xi_2,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f));
                  const __m256 xi_31 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_29,_mm256_set_ps(omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a)),_mm256_mul_ps(_mm256_mul_ps(xi_2,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f)),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),xi_29),xi_30);
                  const __m256 xi_36 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(u_2_a_collide,xia_22_collide),_mm256_set_ps(-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f)),_mm256_set_ps(omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256 xi_40 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_10,xi_30),xi_36),xi_39);
                  const __m256 xi_54 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_10,xi_30),xi_36),xi_53);
                  const __m256 xi_60 = _mm256_mul_ps(xi_11,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_62 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_21,xi_53),xi_58),xi_59),xi_60),xi_61),xi_9);
                  const __m256 xi_64 = _mm256_mul_ps(u_2_a_collide,xi_41);
                  const __m256 xi_66 = _mm256_mul_ps(u_2_a_collide,xi_46);
                  const __m256 xi_67 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_65,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_66,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_63),xi_64);
                  const __m256 xi_68 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_63,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_64,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_65),xi_66);
                  const __m256 xi_69 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_21,xi_39),xi_58),xi_59),xi_60),xi_61),xi_9);
                  const __m256 xi_72 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_15,xi_22),xi_35),xi_71);
                  const __m256 xi_75 = _mm256_mul_ps(u_2_a_collide,xia_5_collide);
                  const __m256 xi_76 = _mm256_mul_ps(xi_75,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_78 = _mm256_mul_ps(xi_75,_mm256_set_ps(xi_45,xi_45,xi_45,xi_45,xi_45,xi_45,xi_45,xi_45));
                  const __m256 xi_79 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_77,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_78,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_74),xi_76);
                  const __m256 xi_80 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_74,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_76,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_77),xi_78);
                  const __m256 xi_81 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_15,xi_22),xi_51),xi_71);
                  const __m256 xi_164 = _mm256_mul_ps(u_2_a_collide,u_2_a_collide);
                  const __m256 xi_165 = _mm256_mul_ps(_mm256_mul_ps(rho_a_collide,xi_164),_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f));
                  const __m256 xi_172 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_165,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_169,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_171,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(rho_a_collide,xi_166),_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f))),_mm256_set_ps(omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256 xi_176 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_165,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_167,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_173,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_175,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(rho_a_collide,xi_168),_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f))),_mm256_set_ps(omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256 xi_179 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_169,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_173,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_178,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(rho_a_collide,xi_164),_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f))),_mm256_set_ps(omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256 xi_184 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_171,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_175,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_180,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xia_20_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xia_24_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(rho_a_collide,xi_164)),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a));
                  const __m256 xi_210 = _mm256_mul_ps(xi_179,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f));
                  const __m256 xi_211 = _mm256_mul_ps(xi_172,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f));
                  const __m256 xi_216 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_183,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f)),xi_210),xi_211),xi_215);
                  const __m256 xi_224 = _mm256_mul_ps(xi_176,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f));
                  const __m256 xi_226 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_182,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f)),xi_210),xi_224),xi_225);
                  const __m256 xi_227 = _mm256_mul_ps(u_2_a_collide,xi_199);
                  const __m256 xi_228 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_18_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_207),xi_227),xia_19_collide);
                  const __m256 xi_229 = _mm256_mul_ps(xi_228,_mm256_set_ps(xi_204,xi_204,xi_204,xi_204,xi_204,xi_204,xi_204,xi_204));
                  const __m256 xi_230 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_6_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_221),xi_227),xia_25_collide);
                  const __m256 xi_231 = _mm256_mul_ps(xi_230,_mm256_set_ps(xi_204,xi_204,xi_204,xi_204,xi_204,xi_204,xi_204,xi_204));
                  const __m256 xi_233 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_184,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f)),xi_211),xi_224),xi_232);
                  const __m256 xi_244 = _mm256_mul_ps(xi_179,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_246 = _mm256_add_ps(xi_244,xi_245);
                  const __m256 xi_250 = _mm256_add_ps(_mm256_mul_ps(xi_245,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_244);
                  const __m256 xi_253 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_21_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(u_2_a_collide,xi_237)),xi_170),xi_206),_mm256_set_ps(xi_238,xi_238,xi_238,xi_238,xi_238,xi_238,xi_238,xi_238));
                  const __m256 xi_256 = _mm256_mul_ps(xi_172,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_257 = _mm256_mul_ps(xi_228,_mm256_set_ps(xi_234,xi_234,xi_234,xi_234,xi_234,xi_234,xi_234,xi_234));
                  const __m256 xi_258 = _mm256_add_ps(xi_256,xi_257);
                  const __m256 xi_260 = _mm256_mul_ps(xi_253,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_266 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xia_9_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(_mm256_mul_ps(rho_a_collide,u_0_a_collide),u_2_a_collide)),xi_174),xi_220),_mm256_set_ps(xi_238,xi_238,xi_238,xi_238,xi_238,xi_238,xi_238,xi_238));
                  const __m256 xi_267 = _mm256_mul_ps(xi_266,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_270 = _mm256_mul_ps(xi_176,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_271 = _mm256_mul_ps(xi_230,_mm256_set_ps(xi_234,xi_234,xi_234,xi_234,xi_234,xi_234,xi_234,xi_234));
                  const __m256 xi_272 = _mm256_add_ps(xi_270,xi_271);
                  const __m256 xi_276 = _mm256_add_ps(_mm256_mul_ps(xi_257,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_256);
                  const __m256 xi_279 = _mm256_add_ps(_mm256_mul_ps(xi_271,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_270);
                  const __m256 forceTerm_0_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_0,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_1,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_2,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(xi_0,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_3,xi_3,xi_3,xi_3,xi_3,xi_3,xi_3,xi_3))),_mm256_mul_ps(_mm256_mul_ps(xi_1,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_3,xi_3,xi_3,xi_3,xi_3,xi_3,xi_3,xi_3))),_mm256_mul_ps(_mm256_mul_ps(xi_2,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_3,xi_3,xi_3,xi_3,xi_3,xi_3,xi_3,xi_3))),_mm256_mul_ps(_mm256_mul_ps(u_0_a_collide,xia_5_collide),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),_mm256_mul_ps(_mm256_mul_ps(u_1_a_collide,xia_2_collide),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),_mm256_mul_ps(_mm256_mul_ps(u_2_a_collide,xia_22_collide),_mm256_set_ps(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a)));
                  const __m256 forceTerm_1_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_6,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_16),xi_4),xi_7);
                  const __m256 forceTerm_2_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_4,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_7,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_16),xi_6);
                  const __m256 forceTerm_3_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_17,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_18,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_19),xi_25);
                  const __m256 forceTerm_4_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_19,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_17),xi_18),xi_25);
                  const __m256 forceTerm_5_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_27,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_26),xi_28),xi_31);
                  const __m256 forceTerm_6_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_26,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_28,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_27),xi_31);
                  const __m256 forceTerm_7_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_35,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_40,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_50,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_8_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_40,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_51,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_52,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_9_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_35,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_52,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_54,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_10_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_50,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_51,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_54,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_11_a_collide = _mm256_add_ps(_mm256_add_ps(xi_57,xi_62),xi_67);
                  const __m256 forceTerm_12_a_collide = _mm256_add_ps(_mm256_add_ps(xi_57,xi_68),xi_69);
                  const __m256 forceTerm_13_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_70,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_72,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_79,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_14_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_70,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_80,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_81,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_15_a_collide = _mm256_add_ps(_mm256_add_ps(xi_62,xi_68),xi_70);
                  const __m256 forceTerm_16_a_collide = _mm256_add_ps(_mm256_add_ps(xi_67,xi_69),xi_70);
                  const __m256 forceTerm_17_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_57,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_72,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_80,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_18_a_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_57,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_79,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_81,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 xib_1_collide = _mm256_loadu_ps(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3 + ctr_0]);
                  const __m256 xi_86 = _mm256_mul_ps(xib_1_collide,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_89 = _mm256_mul_ps(xi_86,_mm256_set_ps(omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b));
                  const __m256 xi_119 = _mm256_mul_ps(xib_1_collide,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_120 = _mm256_mul_ps(xib_1_collide,_mm256_set_ps(xi_115,xi_115,xi_115,xi_115,xi_115,xi_115,xi_115,xi_115));
                  const __m256 xi_121 = _mm256_add_ps(_mm256_mul_ps(xi_119,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_120);
                  const __m256 xi_123 = _mm256_mul_ps(xib_1_collide,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_128 = _mm256_mul_ps(xib_1_collide,_mm256_set_ps(xi_127,xi_127,xi_127,xi_127,xi_127,xi_127,xi_127,xi_127));
                  const __m256 xi_135 = _mm256_add_ps(_mm256_mul_ps(xi_120,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_119);
                  const __m256 xib_2_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_3_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_4_collide = _mm256_load_ps(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0]);
                  const __m256 xib_5_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_6_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_7_collide = _mm256_loadu_ps(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0]);
                  const __m256 xib_8_collide = _mm256_loadu_ps(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3 + ctr_0]);
                  const __m256 xi_108 = _mm256_mul_ps(xib_8_collide,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_110 = _mm256_mul_ps(xi_108,_mm256_set_ps(omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b));
                  const __m256 xi_137 = _mm256_mul_ps(xib_8_collide,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_138 = _mm256_mul_ps(xib_8_collide,_mm256_set_ps(xi_115,xi_115,xi_115,xi_115,xi_115,xi_115,xi_115,xi_115));
                  const __m256 xi_139 = _mm256_add_ps(_mm256_mul_ps(xi_138,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_137);
                  const __m256 xi_152 = _mm256_add_ps(_mm256_mul_ps(xi_137,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_138);
                  const __m256 xib_9_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_10_collide = _mm256_load_ps(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const __m256 xib_11_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_12_collide = _mm256_load_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_309 = _mm256_mul_ps(xib_12_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_310 = _mm256_add_ps(xi_309,xib_9_collide);
                  const __m256 xib_13_collide = _mm256_load_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0]);
                  const __m256 xib_14_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_15_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_16_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_304 = _mm256_mul_ps(xib_16_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xib_17_collide = _mm256_load_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_305 = _mm256_add_ps(xi_304,xib_17_collide);
                  const __m256 xib_18_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_19_collide = _mm256_loadu_ps(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0]);
                  const __m256 xib_20_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_21_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0]);
                  const __m256 xib_22_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_294 = _mm256_add_ps(xib_22_collide,xib_5_collide);
                  const __m256 xi_295 = _mm256_add_ps(_mm256_add_ps(xi_294,xib_12_collide),xib_9_collide);
                  const __m256 xib_23_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_319 = _mm256_mul_ps(xib_23_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_320 = _mm256_add_ps(xi_319,xib_2_collide);
                  const __m256 xib_24_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_287 = _mm256_add_ps(xib_24_collide,xib_6_collide);
                  const __m256 xi_288 = _mm256_add_ps(_mm256_add_ps(xi_287,xib_16_collide),xib_17_collide);
                  const __m256 xib_25_collide = _mm256_loadu_ps(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0]);
                  const __m256 xi_291 = _mm256_add_ps(xib_25_collide,xib_3_collide);
                  const __m256 xi_292 = _mm256_add_ps(_mm256_add_ps(xi_291,xib_23_collide),xib_2_collide);
                  const __m256 xib_26_collide = _mm256_load_ps(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + ctr_0]);
                  const __m256 xi_99 = _mm256_mul_ps(xib_26_collide,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_100 = _mm256_mul_ps(xi_99,_mm256_set_ps(omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b));
                  const __m256 xi_114 = _mm256_mul_ps(xib_26_collide,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_116 = _mm256_mul_ps(xib_26_collide,_mm256_set_ps(xi_115,xi_115,xi_115,xi_115,xi_115,xi_115,xi_115,xi_115));
                  const __m256 xi_117 = _mm256_add_ps(_mm256_mul_ps(xi_116,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_114);
                  const __m256 xi_133 = _mm256_add_ps(_mm256_mul_ps(xi_114,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_116);
                  const __m256 xi_88 = _mm256_mul_ps(xib_1_collide,_mm256_set_ps(xi_87,xi_87,xi_87,xi_87,xi_87,xi_87,xi_87,xi_87));
                  const __m256 xi_101 = _mm256_mul_ps(xib_26_collide,_mm256_set_ps(xi_87,xi_87,xi_87,xi_87,xi_87,xi_87,xi_87,xi_87));
                  const __m256 xi_109 = _mm256_mul_ps(xib_8_collide,_mm256_set_ps(xi_87,xi_87,xi_87,xi_87,xi_87,xi_87,xi_87,xi_87));
                  const __m256 rho_b_collide = xib_10_collide;
                  const __m256 xi_284 = _mm256_mul_ps(rho_b_collide,_mm256_set_ps(-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f,-0.1111111111111111f));
                  const __m256 xi_297 = _mm256_mul_ps(rho_b_collide,_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f));
                  const __m256 xi_298 = _mm256_add_ps(xi_295,xi_297);
                  const __m256 xi_302 = _mm256_mul_ps(rho_b_collide,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f));
                  const __m256 u_0_b_collide = xib_4_collide;
                  const __m256 xi_82 = _mm256_mul_ps(u_0_b_collide,xib_26_collide);
                  const __m256 xi_95 = _mm256_mul_ps(xi_82,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_96 = _mm256_add_ps(_mm256_mul_ps(xi_82,_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f)),_mm256_mul_ps(xi_95,_mm256_set_ps(omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b)));
                  const __m256 xi_102 = _mm256_mul_ps(xi_82,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f));
                  const __m256 xi_124 = _mm256_mul_ps(u_0_b_collide,xi_123);
                  const __m256 xi_129 = _mm256_mul_ps(u_0_b_collide,xi_128);
                  const __m256 xi_140 = _mm256_mul_ps(xi_95,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_143 = _mm256_mul_ps(_mm256_mul_ps(xi_82,_mm256_set_ps(0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f,0.041666666666666664f)),_mm256_set_ps(omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256 xi_155 = _mm256_mul_ps(u_0_b_collide,xib_8_collide);
                  const __m256 xi_156 = _mm256_mul_ps(xi_155,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_159 = _mm256_mul_ps(xi_155,_mm256_set_ps(xi_127,xi_127,xi_127,xi_127,xi_127,xi_127,xi_127,xi_127));
                  const __m256 xi_283 = _mm256_mul_ps(u_0_b_collide,u_0_b_collide);
                  const __m256 xi_290 = _mm256_mul_ps(_mm256_mul_ps(rho_b_collide,xi_283),_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f));
                  const __m256 xi_299 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_292,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_298,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xib_14_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xib_15_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(rho_b_collide,xi_283)),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b));
                  const __m256 xi_316 = _mm256_mul_ps(u_0_b_collide,xi_302);
                  const __m256 xi_317 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_5_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_310),xi_316),xib_22_collide);
                  const __m256 xi_318 = _mm256_mul_ps(xi_317,_mm256_set_ps(xi_307,xi_307,xi_307,xi_307,xi_307,xi_307,xi_307,xi_307));
                  const __m256 xi_321 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_25_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_316),xi_320),xib_3_collide);
                  const __m256 xi_322 = _mm256_mul_ps(xi_321,_mm256_set_ps(xi_307,xi_307,xi_307,xi_307,xi_307,xi_307,xi_307,xi_307));
                  const __m256 xi_332 = _mm256_mul_ps(xi_317,_mm256_set_ps(xi_331,xi_331,xi_331,xi_331,xi_331,xi_331,xi_331,xi_331));
                  const __m256 xi_333 = _mm256_mul_ps(xi_332,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_349 = _mm256_mul_ps(xi_321,_mm256_set_ps(xi_331,xi_331,xi_331,xi_331,xi_331,xi_331,xi_331,xi_331));
                  const __m256 xi_350 = _mm256_mul_ps(xi_349,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 u_1_b_collide = xib_7_collide;
                  const __m256 xi_83 = _mm256_mul_ps(u_1_b_collide,xib_1_collide);
                  const __m256 xi_90 = _mm256_mul_ps(xi_83,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f));
                  const __m256 xi_103 = _mm256_mul_ps(xi_83,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_104 = _mm256_mul_ps(xi_83,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_105 = _mm256_mul_ps(xi_104,_mm256_set_ps(omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256 xi_106 = _mm256_add_ps(_mm256_mul_ps(xi_103,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_105);
                  const __m256 xi_112 = _mm256_add_ps(xi_106,xi_96);
                  const __m256 xi_125 = _mm256_mul_ps(u_1_b_collide,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_126 = _mm256_mul_ps(xi_125,xib_26_collide);
                  const __m256 xi_130 = _mm256_mul_ps(u_1_b_collide,_mm256_set_ps(xi_127,xi_127,xi_127,xi_127,xi_127,xi_127,xi_127,xi_127));
                  const __m256 xi_131 = _mm256_mul_ps(xi_130,xib_26_collide);
                  const __m256 xi_132 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_129,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_131,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_124),xi_126);
                  const __m256 xi_134 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_124,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_126,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_129),xi_131);
                  const __m256 xi_141 = _mm256_mul_ps(xi_105,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_145 = _mm256_mul_ps(xi_125,xib_8_collide);
                  const __m256 xi_147 = _mm256_mul_ps(xi_130,xib_8_collide);
                  const __m256 xi_153 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(u_1_b_collide,xib_1_collide),_mm256_set_ps(-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f)),_mm256_set_ps(omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256 xi_285 = _mm256_mul_ps(u_1_b_collide,u_1_b_collide);
                  const __m256 xi_286 = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(rho_b_collide,xi_285),_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f)),xi_284);
                  const __m256 xi_300 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_288,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_298,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xib_11_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xib_21_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(rho_b_collide,xi_285)),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b));
                  const __m256 xi_303 = _mm256_mul_ps(u_1_b_collide,xi_302);
                  const __m256 xi_306 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_24_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_303),xi_305),xib_6_collide);
                  const __m256 xi_308 = _mm256_mul_ps(xi_306,_mm256_set_ps(xi_307,xi_307,xi_307,xi_307,xi_307,xi_307,xi_307,xi_307));
                  const __m256 xi_311 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_22_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_303),xi_310),xib_5_collide);
                  const __m256 xi_312 = _mm256_mul_ps(xi_311,_mm256_set_ps(xi_307,xi_307,xi_307,xi_307,xi_307,xi_307,xi_307,xi_307));
                  const __m256 xi_334 = _mm256_mul_ps(rho_b_collide,u_1_b_collide);
                  const __m256 xi_336 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_9_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(u_0_b_collide,xi_334)),xi_294),xi_309),_mm256_set_ps(xi_335,xi_335,xi_335,xi_335,xi_335,xi_335,xi_335,xi_335));
                  const __m256 xi_337 = _mm256_mul_ps(xi_336,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_339 = _mm256_mul_ps(xi_311,_mm256_set_ps(xi_331,xi_331,xi_331,xi_331,xi_331,xi_331,xi_331,xi_331));
                  const __m256 xi_342 = _mm256_mul_ps(xi_306,_mm256_set_ps(xi_331,xi_331,xi_331,xi_331,xi_331,xi_331,xi_331,xi_331));
                  const __m256 xi_347 = _mm256_mul_ps(xi_342,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 u_2_b_collide = xib_19_collide;
                  const __m256 xi_84 = _mm256_mul_ps(u_2_b_collide,xib_8_collide);
                  const __m256 xi_91 = _mm256_mul_ps(xi_84,_mm256_set_ps(0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f,0.16666666666666666f));
                  const __m256 xi_92 = _mm256_mul_ps(xi_84,_mm256_set_ps(0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f,0.083333333333333329f));
                  const __m256 xi_93 = _mm256_mul_ps(xi_92,_mm256_set_ps(omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256 xi_94 = _mm256_add_ps(_mm256_mul_ps(xi_91,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_93);
                  const __m256 xi_97 = _mm256_add_ps(xi_94,xi_96);
                  const __m256 xi_98 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_90,_mm256_set_ps(omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b)),_mm256_mul_ps(_mm256_mul_ps(xi_83,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f)),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),xi_90),xi_97);
                  const __m256 xi_107 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_102,_mm256_set_ps(omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b)),_mm256_mul_ps(_mm256_mul_ps(xi_82,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f)),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),xi_102),xi_106),xi_94);
                  const __m256 xi_111 = _mm256_mul_ps(xi_84,_mm256_set_ps(0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f,0.33333333333333331f));
                  const __m256 xi_113 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_111,_mm256_set_ps(omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b)),_mm256_mul_ps(_mm256_mul_ps(xi_84,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f)),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),xi_111),xi_112);
                  const __m256 xi_118 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(u_2_b_collide,xib_8_collide),_mm256_set_ps(-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f,-0.041666666666666664f)),_mm256_set_ps(omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256 xi_122 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_112,xi_118),xi_121),xi_92);
                  const __m256 xi_136 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_112,xi_118),xi_135),xi_92);
                  const __m256 xi_142 = _mm256_mul_ps(xi_93,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_144 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_103,xi_135),xi_140),xi_141),xi_142),xi_143),xi_91);
                  const __m256 xi_146 = _mm256_mul_ps(u_2_b_collide,xi_123);
                  const __m256 xi_148 = _mm256_mul_ps(u_2_b_collide,xi_128);
                  const __m256 xi_149 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_147,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_148,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_145),xi_146);
                  const __m256 xi_150 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_145,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_146,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_147),xi_148);
                  const __m256 xi_151 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_103,xi_121),xi_140),xi_141),xi_142),xi_143),xi_91);
                  const __m256 xi_154 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_104,xi_117),xi_153),xi_97);
                  const __m256 xi_157 = _mm256_mul_ps(u_2_b_collide,xib_26_collide);
                  const __m256 xi_158 = _mm256_mul_ps(xi_157,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_160 = _mm256_mul_ps(xi_157,_mm256_set_ps(xi_127,xi_127,xi_127,xi_127,xi_127,xi_127,xi_127,xi_127));
                  const __m256 xi_161 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_159,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_160,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_156),xi_158);
                  const __m256 xi_162 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_156,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_158,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_159),xi_160);
                  const __m256 xi_163 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_104,xi_133),xi_153),xi_97);
                  const __m256 xi_281 = _mm256_mul_ps(u_2_b_collide,u_2_b_collide);
                  const __m256 xi_282 = _mm256_mul_ps(_mm256_mul_ps(rho_b_collide,xi_281),_mm256_set_ps(-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f,-0.33333333333333331f));
                  const __m256 xi_289 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_282,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_286,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_288,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(rho_b_collide,xi_283),_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f))),_mm256_set_ps(omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256 xi_293 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_282,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_284,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_290,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_292,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(rho_b_collide,xi_285),_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f))),_mm256_set_ps(omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256 xi_296 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_286,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_290,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_295,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(rho_b_collide,xi_281),_mm256_set_ps(-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f,-0.16666666666666666f))),_mm256_set_ps(omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256 xi_301 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_288,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_292,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_297,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xib_18_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xib_20_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(rho_b_collide,xi_281)),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b));
                  const __m256 xi_313 = _mm256_mul_ps(xi_296,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f));
                  const __m256 xi_314 = _mm256_mul_ps(xi_289,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f));
                  const __m256 xi_315 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_300,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f)),xi_215),xi_313),xi_314);
                  const __m256 xi_323 = _mm256_mul_ps(xi_293,_mm256_set_ps(-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f,-0.5f));
                  const __m256 xi_324 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_299,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f)),xi_225),xi_313),xi_323);
                  const __m256 xi_325 = _mm256_mul_ps(u_2_b_collide,xi_302);
                  const __m256 xi_326 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_6_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_305),xi_325),xib_24_collide);
                  const __m256 xi_327 = _mm256_mul_ps(xi_326,_mm256_set_ps(xi_307,xi_307,xi_307,xi_307,xi_307,xi_307,xi_307,xi_307));
                  const __m256 xi_328 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_3_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_320),xi_325),xib_25_collide);
                  const __m256 xi_329 = _mm256_mul_ps(xi_328,_mm256_set_ps(xi_307,xi_307,xi_307,xi_307,xi_307,xi_307,xi_307,xi_307));
                  const __m256 xi_330 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_301,_mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f)),xi_232),xi_314),xi_323);
                  const __m256 xi_338 = _mm256_mul_ps(xi_296,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_340 = _mm256_add_ps(xi_338,xi_339);
                  const __m256 xi_341 = _mm256_add_ps(_mm256_mul_ps(xi_339,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_338);
                  const __m256 xi_343 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_17_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(u_2_b_collide,xi_334)),xi_287),xi_304),_mm256_set_ps(xi_335,xi_335,xi_335,xi_335,xi_335,xi_335,xi_335,xi_335));
                  const __m256 xi_344 = _mm256_mul_ps(xi_289,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_345 = _mm256_mul_ps(xi_326,_mm256_set_ps(xi_331,xi_331,xi_331,xi_331,xi_331,xi_331,xi_331,xi_331));
                  const __m256 xi_346 = _mm256_add_ps(xi_344,xi_345);
                  const __m256 xi_348 = _mm256_mul_ps(xi_343,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_351 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xib_2_collide,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(_mm256_mul_ps(rho_b_collide,u_0_b_collide),u_2_b_collide)),xi_291),xi_319),_mm256_set_ps(xi_335,xi_335,xi_335,xi_335,xi_335,xi_335,xi_335,xi_335));
                  const __m256 xi_352 = _mm256_mul_ps(xi_351,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f));
                  const __m256 xi_353 = _mm256_mul_ps(xi_293,_mm256_set_ps(0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f,0.25f));
                  const __m256 xi_354 = _mm256_mul_ps(xi_328,_mm256_set_ps(xi_331,xi_331,xi_331,xi_331,xi_331,xi_331,xi_331,xi_331));
                  const __m256 xi_355 = _mm256_add_ps(xi_353,xi_354);
                  const __m256 xi_356 = _mm256_add_ps(_mm256_mul_ps(xi_345,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_344);
                  const __m256 xi_357 = _mm256_add_ps(_mm256_mul_ps(xi_354,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_353);
                  const __m256 forceTerm_0_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_82,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_83,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_84,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(_mm256_mul_ps(xi_82,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_85,xi_85,xi_85,xi_85,xi_85,xi_85,xi_85,xi_85))),_mm256_mul_ps(_mm256_mul_ps(xi_83,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_85,xi_85,xi_85,xi_85,xi_85,xi_85,xi_85,xi_85))),_mm256_mul_ps(_mm256_mul_ps(xi_84,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_set_ps(xi_85,xi_85,xi_85,xi_85,xi_85,xi_85,xi_85,xi_85))),_mm256_mul_ps(_mm256_mul_ps(u_0_b_collide,xib_26_collide),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),_mm256_mul_ps(_mm256_mul_ps(u_1_b_collide,xib_1_collide),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),_mm256_mul_ps(_mm256_mul_ps(u_2_b_collide,xib_8_collide),_mm256_set_ps(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b)));
                  const __m256 forceTerm_1_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_88,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_86),xi_89),xi_98);
                  const __m256 forceTerm_2_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_86,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_89,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_88),xi_98);
                  const __m256 forceTerm_3_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_100,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_99,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_101),xi_107);
                  const __m256 forceTerm_4_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_101,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_100),xi_107),xi_99);
                  const __m256 forceTerm_5_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_109,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xi_108),xi_110),xi_113);
                  const __m256 forceTerm_6_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_108,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_110,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),xi_109),xi_113);
                  const __m256 forceTerm_7_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_117,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_122,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_132,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_8_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_122,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_133,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_134,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_9_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_117,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_134,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_136,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_10_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_132,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_133,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_136,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_11_b_collide = _mm256_add_ps(_mm256_add_ps(xi_139,xi_144),xi_149);
                  const __m256 forceTerm_12_b_collide = _mm256_add_ps(_mm256_add_ps(xi_139,xi_150),xi_151);
                  const __m256 forceTerm_13_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_152,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_154,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_161,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_14_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_152,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_162,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_163,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_15_b_collide = _mm256_add_ps(_mm256_add_ps(xi_144,xi_150),xi_152);
                  const __m256 forceTerm_16_b_collide = _mm256_add_ps(_mm256_add_ps(xi_149,xi_151),xi_152);
                  const __m256 forceTerm_17_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_139,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_154,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_162,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 forceTerm_18_b_collide = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_139,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_161,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_163,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)));
                  const __m256 tmp_a0 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_182,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_183,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_184,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_0_a_collide),xi_172),xi_176),xi_179),xi_198),xia_23_collide);
                  const __m256 tmp_a1 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_205,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_209,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_1_a_collide),xi_216),xia_14_collide);
                  const __m256 tmp_a2 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_2_a_collide,xi_205),xi_209),xi_216),xia_16_collide);
                  const __m256 tmp_a3 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_3_a_collide,xi_219),xi_223),xi_226),xia_10_collide);
                  const __m256 tmp_a4 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_219,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_223,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_4_a_collide),xi_226),xia_26_collide);
                  const __m256 tmp_a5 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_229,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_231,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_5_a_collide),xi_233),xia_20_collide);
                  const __m256 tmp_a6 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_6_a_collide,xi_229),xi_231),xi_233),xia_24_collide);
                  const __m256 tmp_a7 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_7_a_collide,xi_236),xi_240),xi_243),xi_246),xia_13_collide);
                  const __m256 tmp_a8 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_8_a_collide,xi_235),xi_239),xi_246),xi_248),xia_15_collide);
                  const __m256 tmp_a9 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_9_a_collide,xi_236),xi_239),xi_249),xi_250),xia_12_collide);
                  const __m256 tmp_a10 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_10_a_collide,xi_235),xi_240),xi_250),xi_251),xia_7_collide);
                  const __m256 tmp_a11 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_11_a_collide,xi_252),xi_253),xi_255),xi_258),xia_8_collide);
                  const __m256 tmp_a12 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_12_a_collide,xi_258),xi_259),xi_260),xi_263),xia_18_collide);
                  const __m256 tmp_a13 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_13_a_collide,xi_265),xi_267),xi_269),xi_272),xia_6_collide);
                  const __m256 tmp_a14 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_14_a_collide,xi_264),xi_266),xi_272),xi_274),xia_4_collide);
                  const __m256 tmp_a15 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_15_a_collide,xi_252),xi_260),xi_275),xi_276),xia_19_collide);
                  const __m256 tmp_a16 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_16_a_collide,xi_253),xi_259),xi_276),xi_277),xia_21_collide);
                  const __m256 tmp_a17 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_17_a_collide,xi_265),xi_266),xi_278),xi_279),xia_9_collide);
                  const __m256 tmp_a18 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_18_a_collide,xi_264),xi_267),xi_279),xi_280),xia_25_collide);
                  const __m256 tmp_b0 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_299,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_300,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),_mm256_mul_ps(xi_301,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_0_b_collide),xi_198),xi_289),xi_293),xi_296),xib_13_collide);
                  const __m256 tmp_b1 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_308,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_312,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_1_b_collide),xi_315),xib_21_collide);
                  const __m256 tmp_b2 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_2_b_collide,xi_308),xi_312),xi_315),xib_11_collide);
                  const __m256 tmp_b3 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_3_b_collide,xi_318),xi_322),xi_324),xib_15_collide);
                  const __m256 tmp_b4 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_318,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_322,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_4_b_collide),xi_324),xib_14_collide);
                  const __m256 tmp_b5 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_327,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_329,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f))),forceTerm_5_b_collide),xi_330),xib_18_collide);
                  const __m256 tmp_b6 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_6_b_collide,xi_327),xi_329),xi_330),xib_20_collide);
                  const __m256 tmp_b7 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_7_b_collide,xi_243),xi_333),xi_337),xi_340),xib_22_collide);
                  const __m256 tmp_b8 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_8_b_collide,xi_248),xi_332),xi_336),xi_340),xib_12_collide);
                  const __m256 tmp_b9 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_9_b_collide,xi_249),xi_333),xi_336),xi_341),xib_9_collide);
                  const __m256 tmp_b10 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_10_b_collide,xi_251),xi_332),xi_337),xi_341),xib_5_collide);
                  const __m256 tmp_b11 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_11_b_collide,xi_255),xi_342),xi_343),xi_346),xib_16_collide);
                  const __m256 tmp_b12 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_12_b_collide,xi_263),xi_346),xi_347),xi_348),xib_6_collide);
                  const __m256 tmp_b13 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_13_b_collide,xi_269),xi_350),xi_352),xi_355),xib_3_collide);
                  const __m256 tmp_b14 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_14_b_collide,xi_274),xi_349),xi_351),xi_355),xib_23_collide);
                  const __m256 tmp_b15 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_15_b_collide,xi_275),xi_342),xi_348),xi_356),xib_24_collide);
                  const __m256 tmp_b16 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_16_b_collide,xi_277),xi_343),xi_347),xi_356),xib_17_collide);
                  const __m256 tmp_b17 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_17_b_collide,xi_278),xi_350),xi_351),xi_357),xib_2_collide);
                  const __m256 tmp_b18 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_18_b_collide,xi_280),xi_349),xi_352),xi_357),xib_25_collide);
                  const __m256 xirecolor_0 = _mm256_add_ps(tmp_a0,tmp_b0);
                  const __m256 xirecolor_1 = _mm256_add_ps(_mm256_load_ps(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]),_mm256_load_ps(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256 xirecolor_2 = _mm256_div_ps(_mm256_set_ps(1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f),xirecolor_1);
                  const __m256 xi_364 = _mm256_mul_ps(xirecolor_2,_mm256_load_ps(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256 xirecolor_3 = _mm256_mul_ps(xirecolor_2,_mm256_load_ps(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]));
                  const __m256 xirecolor_4 = _mm256_add_ps(tmp_a1,tmp_b1);
                  const __m256 xirecolor_5 = xi_189;
                  const __m256 xirecolor_6 = _mm256_div_ps(_mm256_set_ps(1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f),xirecolor_5);
                  const __m256 xirecolor_7 = _mm256_mul_ps(xirecolor_6,_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]));
                  const __m256 xirecolor_8 = _mm256_cmp_ps(xirecolor_5,_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_CMP_NLE_UQ);
                  const __m256 xirecolor_9 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set_ps(beta,beta,beta,beta,beta,beta,beta,beta),_mm256_div_ps(_mm256_set_ps(1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f),_mm256_mul_ps(xirecolor_1,xirecolor_1))),_mm256_load_ps(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_load_ps(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256 xirecolor_10 = _mm256_mul_ps(xirecolor_9,_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_load_ps(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_mul_ps(_mm256_set_ps(0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f,0.055555555555555552f),_mm256_load_ps(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]))));
                  const __m256 xirecolor_11 = _mm256_mul_ps(xirecolor_10,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xirecolor_7,xirecolor_8));
                  const __m256 xirecolor_12 = _mm256_add_ps(tmp_a2,tmp_b2);
                  const __m256 xirecolor_13 = _mm256_mul_ps(xirecolor_10,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xirecolor_7,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_8));
                  const __m256 xirecolor_14 = _mm256_add_ps(tmp_a3,tmp_b3);
                  const __m256 xirecolor_15 = _mm256_mul_ps(xirecolor_6,_mm256_load_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]));
                  const __m256 xirecolor_16 = _mm256_mul_ps(xirecolor_10,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xirecolor_15,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_8));
                  const __m256 xirecolor_17 = _mm256_add_ps(tmp_a4,tmp_b4);
                  const __m256 xirecolor_18 = _mm256_mul_ps(xirecolor_10,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xirecolor_15,xirecolor_8));
                  const __m256 xirecolor_19 = _mm256_add_ps(tmp_a5,tmp_b5);
                  const __m256 xirecolor_20 = _mm256_mul_ps(xirecolor_6,_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]));
                  const __m256 xirecolor_21 = _mm256_mul_ps(xirecolor_10,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xirecolor_20,xirecolor_8));
                  const __m256 xirecolor_22 = _mm256_add_ps(tmp_a6,tmp_b6);
                  const __m256 xirecolor_23 = _mm256_mul_ps(xirecolor_10,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xirecolor_20,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_8));
                  const __m256 xirecolor_24 = _mm256_add_ps(tmp_a7,tmp_b7);
                  const __m256 xirecolor_25 = xi_241;
                  const __m256 xirecolor_26 = _mm256_mul_ps(xirecolor_6,_mm256_set_ps(0.70710678118654757f,0.70710678118654757f,0.70710678118654757f,0.70710678118654757f,0.70710678118654757f,0.70710678118654757f,0.70710678118654757f,0.70710678118654757f));
                  const __m256 xi_358 = _mm256_mul_ps(xirecolor_25,xirecolor_26);
                  const __m256 xirecolor_27 = _mm256_mul_ps(xirecolor_9,_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f),_mm256_load_ps(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_mul_ps(_mm256_set_ps(0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f,0.027777777777777776f),_mm256_load_ps(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]))));
                  const __m256 xirecolor_28 = _mm256_mul_ps(xirecolor_27,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_358,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_8));
                  const __m256 xirecolor_29 = _mm256_add_ps(tmp_a8,tmp_b8);
                  const __m256 xirecolor_30 = xi_247;
                  const __m256 xi_359 = _mm256_mul_ps(xirecolor_26,xirecolor_30);
                  const __m256 xirecolor_31 = _mm256_mul_ps(xirecolor_27,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_359,xirecolor_8));
                  const __m256 xirecolor_32 = _mm256_add_ps(tmp_a9,tmp_b9);
                  const __m256 xirecolor_33 = _mm256_mul_ps(xirecolor_27,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_359,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_8));
                  const __m256 xirecolor_34 = _mm256_add_ps(tmp_a10,tmp_b10);
                  const __m256 xirecolor_35 = _mm256_mul_ps(xirecolor_27,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_358,xirecolor_8));
                  const __m256 xirecolor_36 = _mm256_add_ps(tmp_a11,tmp_b11);
                  const __m256 xirecolor_37 = xi_254;
                  const __m256 xi_360 = _mm256_mul_ps(xirecolor_26,xirecolor_37);
                  const __m256 xirecolor_38 = _mm256_mul_ps(xirecolor_27,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_360,xirecolor_8));
                  const __m256 xirecolor_39 = _mm256_add_ps(tmp_a12,tmp_b12);
                  const __m256 xirecolor_40 = xi_261;
                  const __m256 xirecolor_41 = _mm256_add_ps(xirecolor_40,_mm256_loadu_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]));
                  const __m256 xi_361 = _mm256_mul_ps(xirecolor_26,xirecolor_41);
                  const __m256 xirecolor_42 = _mm256_mul_ps(xirecolor_27,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_361,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_8));
                  const __m256 xirecolor_43 = _mm256_add_ps(tmp_a13,tmp_b13);
                  const __m256 xirecolor_44 = _mm256_add_ps(xirecolor_40,_mm256_load_ps(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]));
                  const __m256 xi_362 = _mm256_mul_ps(xirecolor_26,xirecolor_44);
                  const __m256 xirecolor_45 = _mm256_mul_ps(xirecolor_27,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_362,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_8));
                  const __m256 xirecolor_46 = _mm256_add_ps(tmp_a14,tmp_b14);
                  const __m256 xirecolor_47 = xi_273;
                  const __m256 xi_363 = _mm256_mul_ps(xirecolor_26,xirecolor_47);
                  const __m256 xirecolor_48 = _mm256_mul_ps(xirecolor_27,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_363,xirecolor_8));
                  const __m256 xirecolor_49 = _mm256_add_ps(tmp_a15,tmp_b15);
                  const __m256 xirecolor_50 = _mm256_mul_ps(xirecolor_27,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_361,xirecolor_8));
                  const __m256 xirecolor_51 = _mm256_add_ps(tmp_a16,tmp_b16);
                  const __m256 xirecolor_52 = _mm256_mul_ps(xirecolor_27,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_360,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_8));
                  const __m256 xirecolor_53 = _mm256_add_ps(tmp_a17,tmp_b17);
                  const __m256 xirecolor_54 = _mm256_mul_ps(xirecolor_27,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),_mm256_mul_ps(xi_363,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),xirecolor_8));
                  const __m256 xirecolor_55 = _mm256_add_ps(tmp_a18,tmp_b18);
                  const __m256 xirecolor_56 = _mm256_mul_ps(xirecolor_27,_mm256_blendv_ps(_mm256_set_ps(0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f),xi_362,xirecolor_8));
                  _mm256_store_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0],_mm256_mul_ps(xirecolor_0,xirecolor_3));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_4),xirecolor_11));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_12,xirecolor_3),xirecolor_13));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_14,xirecolor_3),xirecolor_16));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_17,xirecolor_3),xirecolor_18));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_19,xirecolor_3),xirecolor_21));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_22,xirecolor_3),xirecolor_23));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_24,xirecolor_3),xirecolor_28));
                  _mm256_store_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_29,xirecolor_3),xirecolor_31));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_32),xirecolor_33));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_34),xirecolor_35));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_36),xirecolor_38));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_39),xirecolor_42));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_43),xirecolor_45));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_46),xirecolor_48));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_49),xirecolor_50));
                  _mm256_store_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_51),xirecolor_52));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_53),xirecolor_54));
                  _mm256_storeu_ps(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_3,xirecolor_55),xirecolor_56));
                  _mm256_store_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0],_mm256_mul_ps(xi_364,xirecolor_0));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_11,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_4)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_13,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_12)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_16,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_14)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_18,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_17)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_21,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_19)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_23,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_22)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_28,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_24)));
                  _mm256_store_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_31,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_29)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_33,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_32)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_35,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_34)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_38,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_36)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_42,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_39)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_45,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_43)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_48,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_46)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_50,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_49)));
                  _mm256_store_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_52,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_51)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_54,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_53)));
                  _mm256_storeu_ps(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0],_mm256_add_ps(_mm256_mul_ps(xirecolor_56,_mm256_set_ps(-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f,-1.0f)),_mm256_mul_ps(xi_364,xirecolor_55)));
               }
               for (int64_t ctr_0 = (int64_t)((_size_color_gradient_0) / (8)) * (8); ctr_0 < _size_color_gradient_0; ctr_0 += 1)
               {
                  const float xi_185 = (_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]*_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]);
                  const float xi_186 = (_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]*_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]);
                  const float xi_187 = (_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]*_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]);
                  const float xi_188 = xi_185 + xi_186 + xi_187;
                  const float xi_189 = powf(xi_188, 0.5f);
                  const float xi_194 = (_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]);
                  const float xi_196 = sigma*xi_189*((0.5f < _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (omega_shear_a): ((-0.5f > _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (omega_shear_b): ((0.0f < _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (xi_191 + xi_193*xi_194 - xi_193*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]): (xi_191 + xi_194*xi_195 + xi_195*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))));
                  const bool xi_197 = xi_189 > 0.0f;
                  const float xi_198 = ((xi_197) ? (xi_196*0.25f): (0.0f));
                  const float xi_212 = ((1.0f) / (xi_188));
                  const float xi_213 = xi_212*0.055555555555555552f;
                  const float xi_214 = xi_196*1.125f;
                  const float xi_215 = ((xi_197) ? (xi_214*(xi_186*xi_213 - 0.018518518518518517f)): (0.0f));
                  const float xi_225 = ((xi_197) ? (xi_214*(xi_185*xi_213 - 0.018518518518518517f)): (0.0f));
                  const float xi_232 = ((xi_197) ? (xi_214*(xi_187*xi_213 - 0.018518518518518517f)): (0.0f));
                  const float xi_241 = -_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0] + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0];
                  const float xi_242 = xi_212*0.027777777777777776f;
                  const float xi_243 = ((xi_197) ? (xi_214*(xi_242*(xi_241*xi_241) - 0.037037037037037035f)): (0.0f));
                  const float xi_247 = _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0] + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0];
                  const float xi_248 = ((xi_197) ? (xi_214*(xi_242*(xi_247*xi_247) - 0.037037037037037035f)): (0.0f));
                  const float xi_249 = ((xi_197) ? (xi_214*(xi_242*(xi_247*xi_247) - 0.037037037037037035f)): (0.0f));
                  const float xi_251 = ((xi_197) ? (xi_214*(xi_242*(xi_241*xi_241) - 0.037037037037037035f)): (0.0f));
                  const float xi_254 = _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0] + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0];
                  const float xi_255 = ((xi_197) ? (xi_214*(xi_242*(xi_254*xi_254) - 0.037037037037037035f)): (0.0f));
                  const float xi_261 = -_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0];
                  const float xi_262 = xi_261 + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0];
                  const float xi_263 = ((xi_197) ? (xi_214*(xi_242*(xi_262*xi_262) - 0.037037037037037035f)): (0.0f));
                  const float xi_268 = xi_261 + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0];
                  const float xi_269 = ((xi_197) ? (xi_214*(xi_242*(xi_268*xi_268) - 0.037037037037037035f)): (0.0f));
                  const float xi_273 = _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0] + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0];
                  const float xi_274 = ((xi_197) ? (xi_214*(xi_242*(xi_273*xi_273) - 0.037037037037037035f)): (0.0f));
                  const float xi_275 = ((xi_197) ? (xi_214*(xi_242*(xi_262*xi_262) - 0.037037037037037035f)): (0.0f));
                  const float xi_277 = ((xi_197) ? (xi_214*(xi_242*(xi_254*xi_254) - 0.037037037037037035f)): (0.0f));
                  const float xi_278 = ((xi_197) ? (xi_214*(xi_242*(xi_273*xi_273) - 0.037037037037037035f)): (0.0f));
                  const float xi_280 = ((xi_197) ? (xi_214*(xi_242*(xi_268*xi_268) - 0.037037037037037035f)): (0.0f));
                  const float xia_1_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0];
                  const float xia_2_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3 + ctr_0];
                  const float xi_4 = xia_2_collide*0.16666666666666666f;
                  const float xi_7 = omega_odd_a*xi_4;
                  const float xi_37 = xia_2_collide*0.083333333333333329f;
                  const float xi_38 = xi_33*xia_2_collide;
                  const float xi_39 = -xi_37 + xi_38;
                  const float xi_41 = xia_2_collide*0.25f;
                  const float xi_46 = xi_45*xia_2_collide;
                  const float xi_53 = xi_37 - xi_38;
                  const float xia_3_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0];
                  const float xia_4_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0];
                  const float xi_220 = -xia_4_collide;
                  const float xia_5_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + ctr_0];
                  const float xi_17 = xia_5_collide*0.16666666666666666f;
                  const float xi_18 = omega_odd_a*xi_17;
                  const float xi_32 = xia_5_collide*0.083333333333333329f;
                  const float xi_34 = xi_33*xia_5_collide;
                  const float xi_35 = xi_32 - xi_34;
                  const float xi_51 = -xi_32 + xi_34;
                  const float xia_6_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0];
                  const float xia_7_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0];
                  const float xia_8_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0];
                  const float xi_206 = -xia_8_collide;
                  const float xia_9_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0];
                  const float xi_221 = xi_220 + xia_9_collide;
                  const float xia_10_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0];
                  const float xia_11_collide = _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0];
                  const float xia_12_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0];
                  const float xia_13_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0];
                  const float xi_177 = xia_13_collide + xia_7_collide;
                  const float xia_14_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0];
                  const float xia_15_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0];
                  const float xi_178 = xi_177 + xia_12_collide + xia_15_collide;
                  const float xi_201 = -xia_15_collide;
                  const float xi_202 = xi_201 + xia_12_collide;
                  const float xia_16_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0];
                  const float xia_17_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0];
                  const float xia_18_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0];
                  const float xia_19_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0];
                  const float xi_170 = xia_18_collide + xia_19_collide;
                  const float xia_20_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0];
                  const float xia_21_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0];
                  const float xi_171 = xi_170 + xia_21_collide + xia_8_collide;
                  const float xi_207 = xi_206 + xia_21_collide;
                  const float xia_22_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3 + ctr_0];
                  const float xi_26 = xia_22_collide*0.16666666666666666f;
                  const float xi_28 = omega_odd_a*xi_26;
                  const float xi_55 = xia_22_collide*0.083333333333333329f;
                  const float xi_56 = xi_33*xia_22_collide;
                  const float xi_57 = xi_55 - xi_56;
                  const float xi_70 = -xi_55 + xi_56;
                  const float xia_23_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0];
                  const float xia_24_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0];
                  const float xia_25_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0];
                  const float xi_174 = xia_25_collide + xia_6_collide;
                  const float xi_175 = xi_174 + xia_4_collide + xia_9_collide;
                  const float xia_26_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0];
                  const float xi_6 = xi_5*xia_2_collide;
                  const float xi_19 = xi_5*xia_5_collide;
                  const float xi_27 = xi_5*xia_22_collide;
                  const float rho_a_collide = xia_11_collide;
                  const float xi_167 = rho_a_collide*-0.1111111111111111f;
                  const float xi_180 = rho_a_collide*-0.33333333333333331f;
                  const float xi_181 = xi_178 + xi_180;
                  const float xi_199 = rho_a_collide*0.33333333333333331f;
                  const float u_0_a_collide = xia_1_collide;
                  const float xi_0 = u_0_a_collide*xia_5_collide;
                  const float xi_13 = xi_0*0.083333333333333329f;
                  const float xi_14 = omega_even_a*xi_13 + xi_0*-0.16666666666666666f;
                  const float xi_20 = xi_0*0.33333333333333331f;
                  const float xi_42 = u_0_a_collide*xi_41;
                  const float xi_47 = u_0_a_collide*xi_46;
                  const float xi_58 = -xi_13;
                  const float xi_61 = omega_even_a*xi_0*0.041666666666666664f;
                  const float xi_73 = u_0_a_collide*xia_22_collide;
                  const float xi_74 = xi_73*0.25f;
                  const float xi_77 = xi_45*xi_73;
                  const float xi_166 = (u_0_a_collide*u_0_a_collide);
                  const float xi_173 = rho_a_collide*xi_166*-0.33333333333333331f;
                  const float xi_182 = omega_shear_a*(rho_a_collide*xi_166 - xi_175 - xi_181 - xia_10_collide - xia_26_collide);
                  const float xi_217 = u_0_a_collide*xi_199;
                  const float xi_218 = xi_202 + xi_217 + xia_13_collide - xia_7_collide;
                  const float xi_219 = xi_204*xi_218;
                  const float xi_222 = xi_217 + xi_221 - xia_25_collide + xia_6_collide;
                  const float xi_223 = xi_204*xi_222;
                  const float xi_235 = xi_218*xi_234;
                  const float xi_236 = -xi_235;
                  const float xi_264 = xi_222*xi_234;
                  const float xi_265 = -xi_264;
                  const float u_1_a_collide = xia_3_collide;
                  const float xi_1 = u_1_a_collide*xia_2_collide;
                  const float xi_8 = xi_1*0.33333333333333331f;
                  const float xi_21 = xi_1*0.16666666666666666f;
                  const float xi_22 = xi_1*0.083333333333333329f;
                  const float xi_23 = omega_even_a*xi_22;
                  const float xi_24 = -xi_21 + xi_23;
                  const float xi_30 = xi_14 + xi_24;
                  const float xi_43 = u_1_a_collide*0.25f;
                  const float xi_44 = xi_43*xia_5_collide;
                  const float xi_48 = u_1_a_collide*xi_45;
                  const float xi_49 = xi_48*xia_5_collide;
                  const float xi_50 = xi_42 + xi_44 - xi_47 - xi_49;
                  const float xi_52 = -xi_42 - xi_44 + xi_47 + xi_49;
                  const float xi_59 = -xi_23;
                  const float xi_63 = xi_43*xia_22_collide;
                  const float xi_65 = xi_48*xia_22_collide;
                  const float xi_71 = omega_even_a*u_1_a_collide*xia_2_collide*-0.041666666666666664f;
                  const float xi_168 = (u_1_a_collide*u_1_a_collide);
                  const float xi_169 = rho_a_collide*xi_168*-0.33333333333333331f + xi_167;
                  const float xi_183 = omega_shear_a*(rho_a_collide*xi_168 - xi_171 - xi_181 - xia_14_collide - xia_16_collide);
                  const float xi_200 = u_1_a_collide*xi_199;
                  const float xi_203 = xi_200 + xi_202 - xia_13_collide + xia_7_collide;
                  const float xi_205 = xi_203*xi_204;
                  const float xi_208 = xi_200 + xi_207 + xia_18_collide - xia_19_collide;
                  const float xi_209 = xi_204*xi_208;
                  const float xi_237 = rho_a_collide*u_1_a_collide;
                  const float xi_239 = xi_238*(u_0_a_collide*xi_237 + xi_177 + xi_201 - xia_12_collide);
                  const float xi_240 = -xi_239;
                  const float xi_245 = xi_203*xi_234;
                  const float xi_252 = xi_208*xi_234;
                  const float xi_259 = -xi_252;
                  const float u_2_a_collide = xia_17_collide;
                  const float xi_2 = u_2_a_collide*xia_22_collide;
                  const float xi_9 = xi_2*0.16666666666666666f;
                  const float xi_10 = xi_2*0.083333333333333329f;
                  const float xi_11 = omega_even_a*xi_10;
                  const float xi_12 = xi_11 - xi_9;
                  const float xi_15 = xi_12 + xi_14;
                  const float xi_16 = omega_even_a*xi_8 + omega_shear_a*xi_1*-0.5f + xi_15 + xi_8;
                  const float xi_25 = omega_even_a*xi_20 + omega_shear_a*xi_0*-0.5f + xi_12 + xi_20 + xi_24;
                  const float xi_29 = xi_2*0.33333333333333331f;
                  const float xi_31 = omega_even_a*xi_29 + omega_shear_a*xi_2*-0.5f + xi_29 + xi_30;
                  const float xi_36 = omega_even_a*u_2_a_collide*xia_22_collide*-0.041666666666666664f;
                  const float xi_40 = xi_10 + xi_30 + xi_36 + xi_39;
                  const float xi_54 = xi_10 + xi_30 + xi_36 + xi_53;
                  const float xi_60 = -xi_11;
                  const float xi_62 = xi_21 + xi_53 + xi_58 + xi_59 + xi_60 + xi_61 + xi_9;
                  const float xi_64 = u_2_a_collide*xi_41;
                  const float xi_66 = u_2_a_collide*xi_46;
                  const float xi_67 = xi_63 + xi_64 - xi_65 - xi_66;
                  const float xi_68 = -xi_63 - xi_64 + xi_65 + xi_66;
                  const float xi_69 = xi_21 + xi_39 + xi_58 + xi_59 + xi_60 + xi_61 + xi_9;
                  const float xi_72 = xi_15 + xi_22 + xi_35 + xi_71;
                  const float xi_75 = u_2_a_collide*xia_5_collide;
                  const float xi_76 = xi_75*0.25f;
                  const float xi_78 = xi_45*xi_75;
                  const float xi_79 = xi_74 + xi_76 - xi_77 - xi_78;
                  const float xi_80 = -xi_74 - xi_76 + xi_77 + xi_78;
                  const float xi_81 = xi_15 + xi_22 + xi_51 + xi_71;
                  const float xi_164 = (u_2_a_collide*u_2_a_collide);
                  const float xi_165 = rho_a_collide*xi_164*-0.33333333333333331f;
                  const float xi_172 = omega_even_a*(rho_a_collide*xi_166*-0.16666666666666666f - xi_165 - xi_169 - xi_171);
                  const float xi_176 = omega_even_a*(rho_a_collide*xi_168*-0.16666666666666666f - xi_165 - xi_167 - xi_173 - xi_175);
                  const float xi_179 = omega_even_a*(rho_a_collide*xi_164*-0.16666666666666666f - xi_169 - xi_173 - xi_178);
                  const float xi_184 = omega_shear_a*(rho_a_collide*xi_164 - xi_171 - xi_175 - xi_180 - xia_20_collide - xia_24_collide);
                  const float xi_210 = xi_179*-0.5f;
                  const float xi_211 = xi_172*-0.5f;
                  const float xi_216 = xi_183*0.5f + xi_210 + xi_211 + xi_215;
                  const float xi_224 = xi_176*-0.5f;
                  const float xi_226 = xi_182*0.5f + xi_210 + xi_224 + xi_225;
                  const float xi_227 = u_2_a_collide*xi_199;
                  const float xi_228 = xi_207 + xi_227 - xia_18_collide + xia_19_collide;
                  const float xi_229 = xi_204*xi_228;
                  const float xi_230 = xi_221 + xi_227 + xia_25_collide - xia_6_collide;
                  const float xi_231 = xi_204*xi_230;
                  const float xi_233 = xi_184*0.5f + xi_211 + xi_224 + xi_232;
                  const float xi_244 = xi_179*0.25f;
                  const float xi_246 = xi_244 + xi_245;
                  const float xi_250 = xi_244 - xi_245;
                  const float xi_253 = xi_238*(u_2_a_collide*xi_237 + xi_170 + xi_206 - xia_21_collide);
                  const float xi_256 = xi_172*0.25f;
                  const float xi_257 = xi_228*xi_234;
                  const float xi_258 = xi_256 + xi_257;
                  const float xi_260 = -xi_253;
                  const float xi_266 = xi_238*(rho_a_collide*u_0_a_collide*u_2_a_collide + xi_174 + xi_220 - xia_9_collide);
                  const float xi_267 = -xi_266;
                  const float xi_270 = xi_176*0.25f;
                  const float xi_271 = xi_230*xi_234;
                  const float xi_272 = xi_270 + xi_271;
                  const float xi_276 = xi_256 - xi_257;
                  const float xi_279 = xi_270 - xi_271;
                  const float forceTerm_0_a_collide = omega_shear_a*u_0_a_collide*xia_5_collide + omega_shear_a*u_1_a_collide*xia_2_collide + omega_shear_a*u_2_a_collide*xia_22_collide - xi_0*xi_3 - xi_0 - xi_1*xi_3 - xi_1 - xi_2*xi_3 - xi_2;
                  const float forceTerm_1_a_collide = xi_16 + xi_4 - xi_6 + xi_7;
                  const float forceTerm_2_a_collide = xi_16 - xi_4 + xi_6 - xi_7;
                  const float forceTerm_3_a_collide = -xi_17 - xi_18 + xi_19 + xi_25;
                  const float forceTerm_4_a_collide = xi_17 + xi_18 - xi_19 + xi_25;
                  const float forceTerm_5_a_collide = xi_26 - xi_27 + xi_28 + xi_31;
                  const float forceTerm_6_a_collide = -xi_26 + xi_27 - xi_28 + xi_31;
                  const float forceTerm_7_a_collide = -xi_35 - xi_40 - xi_50;
                  const float forceTerm_8_a_collide = -xi_40 - xi_51 - xi_52;
                  const float forceTerm_9_a_collide = -xi_35 - xi_52 - xi_54;
                  const float forceTerm_10_a_collide = -xi_50 - xi_51 - xi_54;
                  const float forceTerm_11_a_collide = xi_57 + xi_62 + xi_67;
                  const float forceTerm_12_a_collide = xi_57 + xi_68 + xi_69;
                  const float forceTerm_13_a_collide = -xi_70 - xi_72 - xi_79;
                  const float forceTerm_14_a_collide = -xi_70 - xi_80 - xi_81;
                  const float forceTerm_15_a_collide = xi_62 + xi_68 + xi_70;
                  const float forceTerm_16_a_collide = xi_67 + xi_69 + xi_70;
                  const float forceTerm_17_a_collide = -xi_57 - xi_72 - xi_80;
                  const float forceTerm_18_a_collide = -xi_57 - xi_79 - xi_81;
                  const float xib_1_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3 + ctr_0];
                  const float xi_86 = xib_1_collide*0.16666666666666666f;
                  const float xi_89 = omega_odd_b*xi_86;
                  const float xi_119 = xib_1_collide*0.083333333333333329f;
                  const float xi_120 = xi_115*xib_1_collide;
                  const float xi_121 = -xi_119 + xi_120;
                  const float xi_123 = xib_1_collide*0.25f;
                  const float xi_128 = xi_127*xib_1_collide;
                  const float xi_135 = xi_119 - xi_120;
                  const float xib_2_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0];
                  const float xib_3_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0];
                  const float xib_4_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0];
                  const float xib_5_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0];
                  const float xib_6_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0];
                  const float xib_7_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0];
                  const float xib_8_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3 + ctr_0];
                  const float xi_108 = xib_8_collide*0.16666666666666666f;
                  const float xi_110 = omega_odd_b*xi_108;
                  const float xi_137 = xib_8_collide*0.083333333333333329f;
                  const float xi_138 = xi_115*xib_8_collide;
                  const float xi_139 = xi_137 - xi_138;
                  const float xi_152 = -xi_137 + xi_138;
                  const float xib_9_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0];
                  const float xib_10_collide = _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const float xib_11_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0];
                  const float xib_12_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0];
                  const float xi_309 = -xib_12_collide;
                  const float xi_310 = xi_309 + xib_9_collide;
                  const float xib_13_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0];
                  const float xib_14_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0];
                  const float xib_15_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0];
                  const float xib_16_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0];
                  const float xi_304 = -xib_16_collide;
                  const float xib_17_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0];
                  const float xi_305 = xi_304 + xib_17_collide;
                  const float xib_18_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0];
                  const float xib_19_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0];
                  const float xib_20_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0];
                  const float xib_21_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0];
                  const float xib_22_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0];
                  const float xi_294 = xib_22_collide + xib_5_collide;
                  const float xi_295 = xi_294 + xib_12_collide + xib_9_collide;
                  const float xib_23_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0];
                  const float xi_319 = -xib_23_collide;
                  const float xi_320 = xi_319 + xib_2_collide;
                  const float xib_24_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0];
                  const float xi_287 = xib_24_collide + xib_6_collide;
                  const float xi_288 = xi_287 + xib_16_collide + xib_17_collide;
                  const float xib_25_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0];
                  const float xi_291 = xib_25_collide + xib_3_collide;
                  const float xi_292 = xi_291 + xib_23_collide + xib_2_collide;
                  const float xib_26_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + ctr_0];
                  const float xi_99 = xib_26_collide*0.16666666666666666f;
                  const float xi_100 = omega_odd_b*xi_99;
                  const float xi_114 = xib_26_collide*0.083333333333333329f;
                  const float xi_116 = xi_115*xib_26_collide;
                  const float xi_117 = xi_114 - xi_116;
                  const float xi_133 = -xi_114 + xi_116;
                  const float xi_88 = xi_87*xib_1_collide;
                  const float xi_101 = xi_87*xib_26_collide;
                  const float xi_109 = xi_87*xib_8_collide;
                  const float rho_b_collide = xib_10_collide;
                  const float xi_284 = rho_b_collide*-0.1111111111111111f;
                  const float xi_297 = rho_b_collide*-0.33333333333333331f;
                  const float xi_298 = xi_295 + xi_297;
                  const float xi_302 = rho_b_collide*0.33333333333333331f;
                  const float u_0_b_collide = xib_4_collide;
                  const float xi_82 = u_0_b_collide*xib_26_collide;
                  const float xi_95 = xi_82*0.083333333333333329f;
                  const float xi_96 = omega_even_b*xi_95 + xi_82*-0.16666666666666666f;
                  const float xi_102 = xi_82*0.33333333333333331f;
                  const float xi_124 = u_0_b_collide*xi_123;
                  const float xi_129 = u_0_b_collide*xi_128;
                  const float xi_140 = -xi_95;
                  const float xi_143 = omega_even_b*xi_82*0.041666666666666664f;
                  const float xi_155 = u_0_b_collide*xib_8_collide;
                  const float xi_156 = xi_155*0.25f;
                  const float xi_159 = xi_127*xi_155;
                  const float xi_283 = (u_0_b_collide*u_0_b_collide);
                  const float xi_290 = rho_b_collide*xi_283*-0.33333333333333331f;
                  const float xi_299 = omega_shear_b*(rho_b_collide*xi_283 - xi_292 - xi_298 - xib_14_collide - xib_15_collide);
                  const float xi_316 = u_0_b_collide*xi_302;
                  const float xi_317 = xi_310 + xi_316 + xib_22_collide - xib_5_collide;
                  const float xi_318 = xi_307*xi_317;
                  const float xi_321 = xi_316 + xi_320 - xib_25_collide + xib_3_collide;
                  const float xi_322 = xi_307*xi_321;
                  const float xi_332 = xi_317*xi_331;
                  const float xi_333 = -xi_332;
                  const float xi_349 = xi_321*xi_331;
                  const float xi_350 = -xi_349;
                  const float u_1_b_collide = xib_7_collide;
                  const float xi_83 = u_1_b_collide*xib_1_collide;
                  const float xi_90 = xi_83*0.33333333333333331f;
                  const float xi_103 = xi_83*0.16666666666666666f;
                  const float xi_104 = xi_83*0.083333333333333329f;
                  const float xi_105 = omega_even_b*xi_104;
                  const float xi_106 = -xi_103 + xi_105;
                  const float xi_112 = xi_106 + xi_96;
                  const float xi_125 = u_1_b_collide*0.25f;
                  const float xi_126 = xi_125*xib_26_collide;
                  const float xi_130 = u_1_b_collide*xi_127;
                  const float xi_131 = xi_130*xib_26_collide;
                  const float xi_132 = xi_124 + xi_126 - xi_129 - xi_131;
                  const float xi_134 = -xi_124 - xi_126 + xi_129 + xi_131;
                  const float xi_141 = -xi_105;
                  const float xi_145 = xi_125*xib_8_collide;
                  const float xi_147 = xi_130*xib_8_collide;
                  const float xi_153 = omega_even_b*u_1_b_collide*xib_1_collide*-0.041666666666666664f;
                  const float xi_285 = (u_1_b_collide*u_1_b_collide);
                  const float xi_286 = rho_b_collide*xi_285*-0.33333333333333331f + xi_284;
                  const float xi_300 = omega_shear_b*(rho_b_collide*xi_285 - xi_288 - xi_298 - xib_11_collide - xib_21_collide);
                  const float xi_303 = u_1_b_collide*xi_302;
                  const float xi_306 = xi_303 + xi_305 - xib_24_collide + xib_6_collide;
                  const float xi_308 = xi_306*xi_307;
                  const float xi_311 = xi_303 + xi_310 - xib_22_collide + xib_5_collide;
                  const float xi_312 = xi_307*xi_311;
                  const float xi_334 = rho_b_collide*u_1_b_collide;
                  const float xi_336 = xi_335*(u_0_b_collide*xi_334 + xi_294 + xi_309 - xib_9_collide);
                  const float xi_337 = -xi_336;
                  const float xi_339 = xi_311*xi_331;
                  const float xi_342 = xi_306*xi_331;
                  const float xi_347 = -xi_342;
                  const float u_2_b_collide = xib_19_collide;
                  const float xi_84 = u_2_b_collide*xib_8_collide;
                  const float xi_91 = xi_84*0.16666666666666666f;
                  const float xi_92 = xi_84*0.083333333333333329f;
                  const float xi_93 = omega_even_b*xi_92;
                  const float xi_94 = -xi_91 + xi_93;
                  const float xi_97 = xi_94 + xi_96;
                  const float xi_98 = omega_even_b*xi_90 + omega_shear_b*xi_83*-0.5f + xi_90 + xi_97;
                  const float xi_107 = omega_even_b*xi_102 + omega_shear_b*xi_82*-0.5f + xi_102 + xi_106 + xi_94;
                  const float xi_111 = xi_84*0.33333333333333331f;
                  const float xi_113 = omega_even_b*xi_111 + omega_shear_b*xi_84*-0.5f + xi_111 + xi_112;
                  const float xi_118 = omega_even_b*u_2_b_collide*xib_8_collide*-0.041666666666666664f;
                  const float xi_122 = xi_112 + xi_118 + xi_121 + xi_92;
                  const float xi_136 = xi_112 + xi_118 + xi_135 + xi_92;
                  const float xi_142 = -xi_93;
                  const float xi_144 = xi_103 + xi_135 + xi_140 + xi_141 + xi_142 + xi_143 + xi_91;
                  const float xi_146 = u_2_b_collide*xi_123;
                  const float xi_148 = u_2_b_collide*xi_128;
                  const float xi_149 = xi_145 + xi_146 - xi_147 - xi_148;
                  const float xi_150 = -xi_145 - xi_146 + xi_147 + xi_148;
                  const float xi_151 = xi_103 + xi_121 + xi_140 + xi_141 + xi_142 + xi_143 + xi_91;
                  const float xi_154 = xi_104 + xi_117 + xi_153 + xi_97;
                  const float xi_157 = u_2_b_collide*xib_26_collide;
                  const float xi_158 = xi_157*0.25f;
                  const float xi_160 = xi_127*xi_157;
                  const float xi_161 = xi_156 + xi_158 - xi_159 - xi_160;
                  const float xi_162 = -xi_156 - xi_158 + xi_159 + xi_160;
                  const float xi_163 = xi_104 + xi_133 + xi_153 + xi_97;
                  const float xi_281 = (u_2_b_collide*u_2_b_collide);
                  const float xi_282 = rho_b_collide*xi_281*-0.33333333333333331f;
                  const float xi_289 = omega_even_b*(rho_b_collide*xi_283*-0.16666666666666666f - xi_282 - xi_286 - xi_288);
                  const float xi_293 = omega_even_b*(rho_b_collide*xi_285*-0.16666666666666666f - xi_282 - xi_284 - xi_290 - xi_292);
                  const float xi_296 = omega_even_b*(rho_b_collide*xi_281*-0.16666666666666666f - xi_286 - xi_290 - xi_295);
                  const float xi_301 = omega_shear_b*(rho_b_collide*xi_281 - xi_288 - xi_292 - xi_297 - xib_18_collide - xib_20_collide);
                  const float xi_313 = xi_296*-0.5f;
                  const float xi_314 = xi_289*-0.5f;
                  const float xi_315 = xi_215 + xi_300*0.5f + xi_313 + xi_314;
                  const float xi_323 = xi_293*-0.5f;
                  const float xi_324 = xi_225 + xi_299*0.5f + xi_313 + xi_323;
                  const float xi_325 = u_2_b_collide*xi_302;
                  const float xi_326 = xi_305 + xi_325 + xib_24_collide - xib_6_collide;
                  const float xi_327 = xi_307*xi_326;
                  const float xi_328 = xi_320 + xi_325 + xib_25_collide - xib_3_collide;
                  const float xi_329 = xi_307*xi_328;
                  const float xi_330 = xi_232 + xi_301*0.5f + xi_314 + xi_323;
                  const float xi_338 = xi_296*0.25f;
                  const float xi_340 = xi_338 + xi_339;
                  const float xi_341 = xi_338 - xi_339;
                  const float xi_343 = xi_335*(u_2_b_collide*xi_334 + xi_287 + xi_304 - xib_17_collide);
                  const float xi_344 = xi_289*0.25f;
                  const float xi_345 = xi_326*xi_331;
                  const float xi_346 = xi_344 + xi_345;
                  const float xi_348 = -xi_343;
                  const float xi_351 = xi_335*(rho_b_collide*u_0_b_collide*u_2_b_collide + xi_291 + xi_319 - xib_2_collide);
                  const float xi_352 = -xi_351;
                  const float xi_353 = xi_293*0.25f;
                  const float xi_354 = xi_328*xi_331;
                  const float xi_355 = xi_353 + xi_354;
                  const float xi_356 = xi_344 - xi_345;
                  const float xi_357 = xi_353 - xi_354;
                  const float forceTerm_0_b_collide = omega_shear_b*u_0_b_collide*xib_26_collide + omega_shear_b*u_1_b_collide*xib_1_collide + omega_shear_b*u_2_b_collide*xib_8_collide - xi_82*xi_85 - xi_82 - xi_83*xi_85 - xi_83 - xi_84*xi_85 - xi_84;
                  const float forceTerm_1_b_collide = xi_86 - xi_88 + xi_89 + xi_98;
                  const float forceTerm_2_b_collide = -xi_86 + xi_88 - xi_89 + xi_98;
                  const float forceTerm_3_b_collide = -xi_100 + xi_101 + xi_107 - xi_99;
                  const float forceTerm_4_b_collide = xi_100 - xi_101 + xi_107 + xi_99;
                  const float forceTerm_5_b_collide = xi_108 - xi_109 + xi_110 + xi_113;
                  const float forceTerm_6_b_collide = -xi_108 + xi_109 - xi_110 + xi_113;
                  const float forceTerm_7_b_collide = -xi_117 - xi_122 - xi_132;
                  const float forceTerm_8_b_collide = -xi_122 - xi_133 - xi_134;
                  const float forceTerm_9_b_collide = -xi_117 - xi_134 - xi_136;
                  const float forceTerm_10_b_collide = -xi_132 - xi_133 - xi_136;
                  const float forceTerm_11_b_collide = xi_139 + xi_144 + xi_149;
                  const float forceTerm_12_b_collide = xi_139 + xi_150 + xi_151;
                  const float forceTerm_13_b_collide = -xi_152 - xi_154 - xi_161;
                  const float forceTerm_14_b_collide = -xi_152 - xi_162 - xi_163;
                  const float forceTerm_15_b_collide = xi_144 + xi_150 + xi_152;
                  const float forceTerm_16_b_collide = xi_149 + xi_151 + xi_152;
                  const float forceTerm_17_b_collide = -xi_139 - xi_154 - xi_162;
                  const float forceTerm_18_b_collide = -xi_139 - xi_161 - xi_163;
                  const float tmp_a0 = forceTerm_0_a_collide + xi_172 + xi_176 + xi_179 - xi_182 - xi_183 - xi_184 + xi_198 + xia_23_collide;
                  const float tmp_a1 = forceTerm_1_a_collide - xi_205 - xi_209 + xi_216 + xia_14_collide;
                  const float tmp_a2 = forceTerm_2_a_collide + xi_205 + xi_209 + xi_216 + xia_16_collide;
                  const float tmp_a3 = forceTerm_3_a_collide + xi_219 + xi_223 + xi_226 + xia_10_collide;
                  const float tmp_a4 = forceTerm_4_a_collide - xi_219 - xi_223 + xi_226 + xia_26_collide;
                  const float tmp_a5 = forceTerm_5_a_collide - xi_229 - xi_231 + xi_233 + xia_20_collide;
                  const float tmp_a6 = forceTerm_6_a_collide + xi_229 + xi_231 + xi_233 + xia_24_collide;
                  const float tmp_a7 = forceTerm_7_a_collide + xi_236 + xi_240 + xi_243 + xi_246 + xia_13_collide;
                  const float tmp_a8 = forceTerm_8_a_collide + xi_235 + xi_239 + xi_246 + xi_248 + xia_15_collide;
                  const float tmp_a9 = forceTerm_9_a_collide + xi_236 + xi_239 + xi_249 + xi_250 + xia_12_collide;
                  const float tmp_a10 = forceTerm_10_a_collide + xi_235 + xi_240 + xi_250 + xi_251 + xia_7_collide;
                  const float tmp_a11 = forceTerm_11_a_collide + xi_252 + xi_253 + xi_255 + xi_258 + xia_8_collide;
                  const float tmp_a12 = forceTerm_12_a_collide + xi_258 + xi_259 + xi_260 + xi_263 + xia_18_collide;
                  const float tmp_a13 = forceTerm_13_a_collide + xi_265 + xi_267 + xi_269 + xi_272 + xia_6_collide;
                  const float tmp_a14 = forceTerm_14_a_collide + xi_264 + xi_266 + xi_272 + xi_274 + xia_4_collide;
                  const float tmp_a15 = forceTerm_15_a_collide + xi_252 + xi_260 + xi_275 + xi_276 + xia_19_collide;
                  const float tmp_a16 = forceTerm_16_a_collide + xi_253 + xi_259 + xi_276 + xi_277 + xia_21_collide;
                  const float tmp_a17 = forceTerm_17_a_collide + xi_265 + xi_266 + xi_278 + xi_279 + xia_9_collide;
                  const float tmp_a18 = forceTerm_18_a_collide + xi_264 + xi_267 + xi_279 + xi_280 + xia_25_collide;
                  const float tmp_b0 = forceTerm_0_b_collide + xi_198 + xi_289 + xi_293 + xi_296 - xi_299 - xi_300 - xi_301 + xib_13_collide;
                  const float tmp_b1 = forceTerm_1_b_collide - xi_308 - xi_312 + xi_315 + xib_21_collide;
                  const float tmp_b2 = forceTerm_2_b_collide + xi_308 + xi_312 + xi_315 + xib_11_collide;
                  const float tmp_b3 = forceTerm_3_b_collide + xi_318 + xi_322 + xi_324 + xib_15_collide;
                  const float tmp_b4 = forceTerm_4_b_collide - xi_318 - xi_322 + xi_324 + xib_14_collide;
                  const float tmp_b5 = forceTerm_5_b_collide - xi_327 - xi_329 + xi_330 + xib_18_collide;
                  const float tmp_b6 = forceTerm_6_b_collide + xi_327 + xi_329 + xi_330 + xib_20_collide;
                  const float tmp_b7 = forceTerm_7_b_collide + xi_243 + xi_333 + xi_337 + xi_340 + xib_22_collide;
                  const float tmp_b8 = forceTerm_8_b_collide + xi_248 + xi_332 + xi_336 + xi_340 + xib_12_collide;
                  const float tmp_b9 = forceTerm_9_b_collide + xi_249 + xi_333 + xi_336 + xi_341 + xib_9_collide;
                  const float tmp_b10 = forceTerm_10_b_collide + xi_251 + xi_332 + xi_337 + xi_341 + xib_5_collide;
                  const float tmp_b11 = forceTerm_11_b_collide + xi_255 + xi_342 + xi_343 + xi_346 + xib_16_collide;
                  const float tmp_b12 = forceTerm_12_b_collide + xi_263 + xi_346 + xi_347 + xi_348 + xib_6_collide;
                  const float tmp_b13 = forceTerm_13_b_collide + xi_269 + xi_350 + xi_352 + xi_355 + xib_3_collide;
                  const float tmp_b14 = forceTerm_14_b_collide + xi_274 + xi_349 + xi_351 + xi_355 + xib_23_collide;
                  const float tmp_b15 = forceTerm_15_b_collide + xi_275 + xi_342 + xi_348 + xi_356 + xib_24_collide;
                  const float tmp_b16 = forceTerm_16_b_collide + xi_277 + xi_343 + xi_347 + xi_356 + xib_17_collide;
                  const float tmp_b17 = forceTerm_17_b_collide + xi_278 + xi_350 + xi_351 + xi_357 + xib_2_collide;
                  const float tmp_b18 = forceTerm_18_b_collide + xi_280 + xi_349 + xi_352 + xi_357 + xib_25_collide;
                  const float xirecolor_0 = tmp_a0 + tmp_b0;
                  const float xirecolor_1 = _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const float xirecolor_2 = ((1.0f) / (xirecolor_1));
                  const float xi_364 = xirecolor_2*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const float xirecolor_3 = xirecolor_2*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0];
                  const float xirecolor_4 = tmp_a1 + tmp_b1;
                  const float xirecolor_5 = xi_189;
                  const float xirecolor_6 = ((1.0f) / (xirecolor_5));
                  const float xirecolor_7 = xirecolor_6*_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0];
                  const bool xirecolor_8 = xirecolor_5 > 0.0f;
                  const float xirecolor_9 = beta*((1.0f) / ((xirecolor_1*xirecolor_1)))*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const float xirecolor_10 = xirecolor_9*(0.055555555555555552f*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + 0.055555555555555552f*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const float xirecolor_11 = xirecolor_10*((xirecolor_8) ? (xirecolor_7): (0.0f));
                  const float xirecolor_12 = tmp_a2 + tmp_b2;
                  const float xirecolor_13 = xirecolor_10*((xirecolor_8) ? (-xirecolor_7): (0.0f));
                  const float xirecolor_14 = tmp_a3 + tmp_b3;
                  const float xirecolor_15 = xirecolor_6*_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0];
                  const float xirecolor_16 = xirecolor_10*((xirecolor_8) ? (-xirecolor_15): (0.0f));
                  const float xirecolor_17 = tmp_a4 + tmp_b4;
                  const float xirecolor_18 = xirecolor_10*((xirecolor_8) ? (xirecolor_15): (0.0f));
                  const float xirecolor_19 = tmp_a5 + tmp_b5;
                  const float xirecolor_20 = xirecolor_6*_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0];
                  const float xirecolor_21 = xirecolor_10*((xirecolor_8) ? (xirecolor_20): (0.0f));
                  const float xirecolor_22 = tmp_a6 + tmp_b6;
                  const float xirecolor_23 = xirecolor_10*((xirecolor_8) ? (-xirecolor_20): (0.0f));
                  const float xirecolor_24 = tmp_a7 + tmp_b7;
                  const float xirecolor_25 = xi_241;
                  const float xirecolor_26 = xirecolor_6*0.70710678118654757f;
                  const float xi_358 = xirecolor_25*xirecolor_26;
                  const float xirecolor_27 = xirecolor_9*(0.027777777777777776f*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + 0.027777777777777776f*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const float xirecolor_28 = xirecolor_27*((xirecolor_8) ? (-xi_358): (0.0f));
                  const float xirecolor_29 = tmp_a8 + tmp_b8;
                  const float xirecolor_30 = xi_247;
                  const float xi_359 = xirecolor_26*xirecolor_30;
                  const float xirecolor_31 = xirecolor_27*((xirecolor_8) ? (xi_359): (0.0f));
                  const float xirecolor_32 = tmp_a9 + tmp_b9;
                  const float xirecolor_33 = xirecolor_27*((xirecolor_8) ? (-xi_359): (0.0f));
                  const float xirecolor_34 = tmp_a10 + tmp_b10;
                  const float xirecolor_35 = xirecolor_27*((xirecolor_8) ? (xi_358): (0.0f));
                  const float xirecolor_36 = tmp_a11 + tmp_b11;
                  const float xirecolor_37 = xi_254;
                  const float xi_360 = xirecolor_26*xirecolor_37;
                  const float xirecolor_38 = xirecolor_27*((xirecolor_8) ? (xi_360): (0.0f));
                  const float xirecolor_39 = tmp_a12 + tmp_b12;
                  const float xirecolor_40 = xi_261;
                  const float xirecolor_41 = xirecolor_40 + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0];
                  const float xi_361 = xirecolor_26*xirecolor_41;
                  const float xirecolor_42 = xirecolor_27*((xirecolor_8) ? (-xi_361): (0.0f));
                  const float xirecolor_43 = tmp_a13 + tmp_b13;
                  const float xirecolor_44 = xirecolor_40 + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0];
                  const float xi_362 = xirecolor_26*xirecolor_44;
                  const float xirecolor_45 = xirecolor_27*((xirecolor_8) ? (-xi_362): (0.0f));
                  const float xirecolor_46 = tmp_a14 + tmp_b14;
                  const float xirecolor_47 = xi_273;
                  const float xi_363 = xirecolor_26*xirecolor_47;
                  const float xirecolor_48 = xirecolor_27*((xirecolor_8) ? (xi_363): (0.0f));
                  const float xirecolor_49 = tmp_a15 + tmp_b15;
                  const float xirecolor_50 = xirecolor_27*((xirecolor_8) ? (xi_361): (0.0f));
                  const float xirecolor_51 = tmp_a16 + tmp_b16;
                  const float xirecolor_52 = xirecolor_27*((xirecolor_8) ? (-xi_360): (0.0f));
                  const float xirecolor_53 = tmp_a17 + tmp_b17;
                  const float xirecolor_54 = xirecolor_27*((xirecolor_8) ? (-xi_363): (0.0f));
                  const float xirecolor_55 = tmp_a18 + tmp_b18;
                  const float xirecolor_56 = xirecolor_27*((xirecolor_8) ? (xi_362): (0.0f));
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0] = xirecolor_0*xirecolor_3;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0] = xirecolor_11 + xirecolor_3*xirecolor_4;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0] = xirecolor_12*xirecolor_3 + xirecolor_13;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0] = xirecolor_14*xirecolor_3 + xirecolor_16;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0] = xirecolor_17*xirecolor_3 + xirecolor_18;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0] = xirecolor_19*xirecolor_3 + xirecolor_21;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0] = xirecolor_22*xirecolor_3 + xirecolor_23;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0] = xirecolor_24*xirecolor_3 + xirecolor_28;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0] = xirecolor_29*xirecolor_3 + xirecolor_31;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_32 + xirecolor_33;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_34 + xirecolor_35;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_36 + xirecolor_38;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_39 + xirecolor_42;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_43 + xirecolor_45;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_46 + xirecolor_48;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_49 + xirecolor_50;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_51 + xirecolor_52;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_53 + xirecolor_54;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_55 + xirecolor_56;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0] = xi_364*xirecolor_0;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_4 - xirecolor_11;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_12 - xirecolor_13;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_14 - xirecolor_16;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_17 - xirecolor_18;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_19 - xirecolor_21;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_22 - xirecolor_23;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_24 - xirecolor_28;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_29 - xirecolor_31;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_32 - xirecolor_33;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_34 - xirecolor_35;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_36 - xirecolor_38;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_39 - xirecolor_42;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_43 - xirecolor_45;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_46 - xirecolor_48;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_49 - xirecolor_50;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_51 - xirecolor_52;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_53 - xirecolor_54;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0] = xi_364*xirecolor_55 - xirecolor_56;
               }
            }
         }
      }
   }
}
}


void ColorGradientCollideSweepSinglePrecisionAVX::run(IBlock * block)
{
   
    auto rho_a = block->getData< field::GhostLayerField<float, 1> >(rho_aID);
    auto pdfs_b = block->getData< field::GhostLayerField<float, 19> >(pdfs_bID);
    auto phasefield = block->getData< field::GhostLayerField<float, 1> >(phasefieldID);
    auto force_b = block->getData< field::GhostLayerField<float, 3> >(force_bID);
    auto color_gradient = block->getData< field::GhostLayerField<float, 3> >(color_gradientID);
    auto pdfs_a = block->getData< field::GhostLayerField<float, 19> >(pdfs_aID);
    auto rho_b = block->getData< field::GhostLayerField<float, 1> >(rho_bID);
    auto velocity = block->getData< field::GhostLayerField<float, 3> >(velocityID);
    auto force_a = block->getData< field::GhostLayerField<float, 3> >(force_aID);

    auto & omega_shear_a = this->omega_shear_a_;
    auto & omega_even_a = this->omega_even_a_;
    auto & beta = this->beta_;
    auto & omega_odd_a = this->omega_odd_a_;
    auto & omega_odd_b = this->omega_odd_b_;
    auto & sigma = this->sigma_;
    auto & omega_even_b = this->omega_even_b_;
    auto & omega_shear_b = this->omega_shear_b_;
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(color_gradient->nrOfGhostLayers()))
    float * RESTRICT const _data_color_gradient = color_gradient->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(color_gradient->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) color_gradient->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force_a->nrOfGhostLayers()))
    float * RESTRICT const _data_force_a = force_a->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force_b->nrOfGhostLayers()))
    float * RESTRICT const _data_force_b = force_b->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs_a->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs_b->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(phasefield->nrOfGhostLayers()))
    float * RESTRICT const _data_phasefield = phasefield->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) phasefield->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_a->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_a = rho_a->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_b->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_b = rho_b->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()))
    float * RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) velocity->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(color_gradient->xSizeWithGhostLayer(), int64_t(int64_c(color_gradient->xSize()) + 0))
    const int64_t _size_color_gradient_0 = int64_t(int64_c(color_gradient->xSize()) + 0);
    WALBERLA_ASSERT_EQUAL(color_gradient->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) color_gradient->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(color_gradient->ySizeWithGhostLayer(), int64_t(int64_c(color_gradient->ySize()) + 0))
    const int64_t _size_color_gradient_1 = int64_t(int64_c(color_gradient->ySize()) + 0);
    WALBERLA_ASSERT_EQUAL(color_gradient->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) color_gradient->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(color_gradient->zSizeWithGhostLayer(), int64_t(int64_c(color_gradient->zSize()) + 0))
    const int64_t _size_color_gradient_2 = int64_t(int64_c(color_gradient->zSize()) + 0);
    WALBERLA_ASSERT_EQUAL(color_gradient->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) color_gradient->dataAt(0, 0, 0, 0) %32, 0)
    const int64_t _stride_color_gradient_1 = int64_t(color_gradient->yStride());
    const int64_t _stride_color_gradient_2 = int64_t(color_gradient->zStride());
    const int64_t _stride_color_gradient_3 = int64_t(1 * int64_t(color_gradient->fStride()));
    const int64_t _stride_force_a_1 = int64_t(force_a->yStride());
    const int64_t _stride_force_a_2 = int64_t(force_a->zStride());
    const int64_t _stride_force_a_3 = int64_t(1 * int64_t(force_a->fStride()));
    const int64_t _stride_force_b_1 = int64_t(force_b->yStride());
    const int64_t _stride_force_b_2 = int64_t(force_b->zStride());
    const int64_t _stride_force_b_3 = int64_t(1 * int64_t(force_b->fStride()));
    const int64_t _stride_pdfs_a_1 = int64_t(pdfs_a->yStride());
    const int64_t _stride_pdfs_a_2 = int64_t(pdfs_a->zStride());
    const int64_t _stride_pdfs_a_3 = int64_t(1 * int64_t(pdfs_a->fStride()));
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
    const int64_t _stride_phasefield_1 = int64_t(phasefield->yStride());
    const int64_t _stride_phasefield_2 = int64_t(phasefield->zStride());
    const int64_t _stride_rho_a_1 = int64_t(rho_a->yStride());
    const int64_t _stride_rho_a_2 = int64_t(rho_a->zStride());
    const int64_t _stride_rho_b_1 = int64_t(rho_b->yStride());
    const int64_t _stride_rho_b_2 = int64_t(rho_b->zStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
    const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
    internal_196de4e60a60b4c0f3c29ed6234bf8bc::colorgradientcollidesweepsingleprecisionavx_colorgradientcollidesweepsingleprecisionavx(_data_color_gradient, _data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_color_gradient_0, _size_color_gradient_1, _size_color_gradient_2, _stride_color_gradient_1, _stride_color_gradient_2, _stride_color_gradient_3, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_even_a, omega_even_b, omega_odd_a, omega_odd_b, omega_shear_a, omega_shear_b, sigma);
    
}


void ColorGradientCollideSweepSinglePrecisionAVX::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block)
{
   
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    auto rho_a = block->getData< field::GhostLayerField<float, 1> >(rho_aID);
    auto pdfs_b = block->getData< field::GhostLayerField<float, 19> >(pdfs_bID);
    auto phasefield = block->getData< field::GhostLayerField<float, 1> >(phasefieldID);
    auto force_b = block->getData< field::GhostLayerField<float, 3> >(force_bID);
    auto color_gradient = block->getData< field::GhostLayerField<float, 3> >(color_gradientID);
    auto pdfs_a = block->getData< field::GhostLayerField<float, 19> >(pdfs_aID);
    auto rho_b = block->getData< field::GhostLayerField<float, 1> >(rho_bID);
    auto velocity = block->getData< field::GhostLayerField<float, 3> >(velocityID);
    auto force_a = block->getData< field::GhostLayerField<float, 3> >(force_aID);

    auto & omega_shear_a = this->omega_shear_a_;
    auto & omega_even_a = this->omega_even_a_;
    auto & beta = this->beta_;
    auto & omega_odd_a = this->omega_odd_a_;
    auto & omega_odd_b = this->omega_odd_b_;
    auto & sigma = this->sigma_;
    auto & omega_even_b = this->omega_even_b_;
    auto & omega_shear_b = this->omega_shear_b_;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(color_gradient->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(color_gradient->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(color_gradient->nrOfGhostLayers()))
    float * RESTRICT const _data_color_gradient = color_gradient->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(color_gradient->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) color_gradient->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force_a->nrOfGhostLayers()))
    float * RESTRICT const _data_force_a = force_a->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force_b->nrOfGhostLayers()))
    float * RESTRICT const _data_force_b = force_b->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_a->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_b->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(phasefield->nrOfGhostLayers()))
    float * RESTRICT const _data_phasefield = phasefield->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) phasefield->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_a->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_a = rho_a->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_b->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_b = rho_b->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()))
    float * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) velocity->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(color_gradient->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_color_gradient_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_EQUAL(color_gradient->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) color_gradient->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(color_gradient->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_color_gradient_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_EQUAL(color_gradient->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) color_gradient->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(color_gradient->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_color_gradient_2 = int64_t(int64_c(ci.zSize()) + 0);
    WALBERLA_ASSERT_EQUAL(color_gradient->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) color_gradient->dataAt(0, 0, 0, 0) %32, 0)
    const int64_t _stride_color_gradient_1 = int64_t(color_gradient->yStride());
    const int64_t _stride_color_gradient_2 = int64_t(color_gradient->zStride());
    const int64_t _stride_color_gradient_3 = int64_t(1 * int64_t(color_gradient->fStride()));
    const int64_t _stride_force_a_1 = int64_t(force_a->yStride());
    const int64_t _stride_force_a_2 = int64_t(force_a->zStride());
    const int64_t _stride_force_a_3 = int64_t(1 * int64_t(force_a->fStride()));
    const int64_t _stride_force_b_1 = int64_t(force_b->yStride());
    const int64_t _stride_force_b_2 = int64_t(force_b->zStride());
    const int64_t _stride_force_b_3 = int64_t(1 * int64_t(force_b->fStride()));
    const int64_t _stride_pdfs_a_1 = int64_t(pdfs_a->yStride());
    const int64_t _stride_pdfs_a_2 = int64_t(pdfs_a->zStride());
    const int64_t _stride_pdfs_a_3 = int64_t(1 * int64_t(pdfs_a->fStride()));
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
    const int64_t _stride_phasefield_1 = int64_t(phasefield->yStride());
    const int64_t _stride_phasefield_2 = int64_t(phasefield->zStride());
    const int64_t _stride_rho_a_1 = int64_t(rho_a->yStride());
    const int64_t _stride_rho_a_2 = int64_t(rho_a->zStride());
    const int64_t _stride_rho_b_1 = int64_t(rho_b->yStride());
    const int64_t _stride_rho_b_2 = int64_t(rho_b->zStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
    const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
    internal_196de4e60a60b4c0f3c29ed6234bf8bc::colorgradientcollidesweepsingleprecisionavx_colorgradientcollidesweepsingleprecisionavx(_data_color_gradient, _data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_color_gradient_0, _size_color_gradient_1, _size_color_gradient_2, _stride_color_gradient_1, _stride_color_gradient_2, _stride_color_gradient_3, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_even_a, omega_even_b, omega_odd_a, omega_odd_b, omega_shear_a, omega_shear_b, sigma);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
