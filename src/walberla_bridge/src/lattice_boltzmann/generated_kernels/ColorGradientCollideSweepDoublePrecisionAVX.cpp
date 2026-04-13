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
//! \\file ColorGradientCollideSweepDoublePrecisionAVX.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 3247aa7395049ca5bfb69d34d55e45db19fa439c


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "ColorGradientCollideSweepDoublePrecisionAVX.h"


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


namespace internal_6ed8e43fbb9656b49ba1d2a7ec0e653c {
static FUNC_PREFIX void colorgradientcollidesweepdoubleprecisionavx_colorgradientcollidesweepdoubleprecisionavx(double * RESTRICT const _data_color_gradient, double * RESTRICT const _data_force_a, double * RESTRICT const _data_force_b, double * RESTRICT  _data_pdfs_a, double * RESTRICT  _data_pdfs_b, double * RESTRICT const _data_phasefield, double * RESTRICT const _data_rho_a, double * RESTRICT const _data_rho_b, double * RESTRICT const _data_velocity, int64_t const _size_color_gradient_0, int64_t const _size_color_gradient_1, int64_t const _size_color_gradient_2, int64_t const _stride_color_gradient_1, int64_t const _stride_color_gradient_2, int64_t const _stride_color_gradient_3, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3, double beta, double omega_even_a, double omega_even_b, double omega_odd_a, double omega_odd_b, double omega_shear_a, double omega_shear_b, double sigma)
{
#ifdef _OPENMP
   #pragma omp parallel
#endif
   {
      const double xi_3 = omega_even_a*0.5;
      const double xi_33 = omega_odd_a*0.041666666666666664;
      const double xi_45 = omega_shear_a*0.125;
      const double xi_85 = omega_even_b*0.5;
      const double xi_115 = omega_odd_b*0.041666666666666664;
      const double xi_127 = omega_shear_b*0.125;
      const double xi_190 = omega_shear_a*omega_shear_b*((1.0) / (omega_shear_a + omega_shear_b));
      const double xi_191 = xi_190*2.0;
      const double xi_192 = xi_190*8.0;
      const double xi_193 = omega_shear_a*-4.0 + xi_192;
      const double xi_195 = omega_shear_b*-4.0 + xi_192;
      const double xi_204 = omega_odd_a*0.5;
      const double xi_234 = omega_odd_a*0.25;
      const double xi_238 = omega_shear_a*0.25;
      const double xi_307 = omega_odd_b*0.5;
      const double xi_331 = omega_odd_b*0.25;
      const double xi_335 = omega_shear_b*0.25;
      const double rr_0_a_collide = 0.0;
      const double xi_5 = rr_0_a_collide*0.25;
      const double rr_0_b_collide = 0.0;
      const double xi_87 = rr_0_b_collide*0.25;
#ifdef _OPENMP
      #pragma omp for schedule(static)
#endif
      for (int64_t ctr_2 = 0; ctr_2 < _size_color_gradient_2; ctr_2 += 1)
      {
         for (int64_t ctr_1 = 0; ctr_1 < _size_color_gradient_1; ctr_1 += 1)
         {
            {
               for (int64_t ctr_0 = 0; ctr_0 < (int64_t)((_size_color_gradient_0) / (4)) * (4); ctr_0 += 4)
               {
                  const __m256d xi_185 = _mm256_mul_pd(_mm256_load_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]),_mm256_load_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]));
                  const __m256d xi_186 = _mm256_mul_pd(_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]),_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]));
                  const __m256d xi_187 = _mm256_mul_pd(_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]),_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]));
                  const __m256d xi_188 = _mm256_add_pd(_mm256_add_pd(xi_185,xi_186),xi_187);
                  const __m256d xi_189 = _mm256_sqrt_pd(xi_188);
                  const __m256d xi_194 = _mm256_mul_pd(_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]));
                  const __m256d xi_196 = _mm256_mul_pd(_mm256_mul_pd(xi_189,_mm256_set_pd(sigma,sigma,sigma,sigma)),_mm256_blendv_pd(_mm256_blendv_pd(_mm256_blendv_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_194,_mm256_set_pd(xi_195,xi_195,xi_195,xi_195)),_mm256_mul_pd(_mm256_set_pd(xi_195,xi_195,xi_195,xi_195),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))),_mm256_set_pd(xi_191,xi_191,xi_191,xi_191)),_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_194,_mm256_set_pd(xi_193,xi_193,xi_193,xi_193)),_mm256_mul_pd(_mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_set_pd(xi_193,xi_193,xi_193,xi_193)),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))),_mm256_set_pd(xi_191,xi_191,xi_191,xi_191)),_mm256_cmp_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NGE_UQ)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),_mm256_cmp_pd(_mm256_set_pd(-0.5,-0.5,-0.5,-0.5),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NLE_UQ)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),_mm256_cmp_pd(_mm256_set_pd(0.5,0.5,0.5,0.5),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NGE_UQ)));
                  const __m256d xi_197 = _mm256_cmp_pd(xi_189,_mm256_set_pd(0.0,0.0,0.0,0.0),_CMP_NLE_UQ);
                  const __m256d xi_198 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_196,_mm256_set_pd(0.25,0.25,0.25,0.25)),xi_197);
                  const __m256d xi_212 = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),xi_188);
                  const __m256d xi_213 = _mm256_mul_pd(xi_212,_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552));
                  const __m256d xi_214 = _mm256_mul_pd(xi_196,_mm256_set_pd(1.125,1.125,1.125,1.125));
                  const __m256d xi_215 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_186,xi_213),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_197);
                  const __m256d xi_225 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_185,xi_213),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_197);
                  const __m256d xi_232 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_187,xi_213),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_197);
                  const __m256d xi_241 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0])),_mm256_load_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]));
                  const __m256d xi_242 = _mm256_mul_pd(xi_212,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776));
                  const __m256d xi_243 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_242,_mm256_mul_pd(xi_241,xi_241)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_197);
                  const __m256d xi_247 = _mm256_add_pd(_mm256_load_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]),_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]));
                  const __m256d xi_248 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_242,_mm256_mul_pd(xi_247,xi_247)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_197);
                  const __m256d xi_249 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_242,_mm256_mul_pd(xi_247,xi_247)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_197);
                  const __m256d xi_251 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_242,_mm256_mul_pd(xi_241,xi_241)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_197);
                  const __m256d xi_254 = _mm256_add_pd(_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]),_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]));
                  const __m256d xi_255 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_242,_mm256_mul_pd(xi_254,xi_254)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_197);
                  const __m256d xi_261 = _mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]));
                  const __m256d xi_262 = _mm256_add_pd(xi_261,_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]));
                  const __m256d xi_263 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_242,_mm256_mul_pd(xi_262,xi_262)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_197);
                  const __m256d xi_268 = _mm256_add_pd(xi_261,_mm256_load_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]));
                  const __m256d xi_269 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_242,_mm256_mul_pd(xi_268,xi_268)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_197);
                  const __m256d xi_273 = _mm256_add_pd(_mm256_load_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]),_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]));
                  const __m256d xi_274 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_242,_mm256_mul_pd(xi_273,xi_273)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_197);
                  const __m256d xi_275 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_242,_mm256_mul_pd(xi_262,xi_262)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_197);
                  const __m256d xi_277 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_242,_mm256_mul_pd(xi_254,xi_254)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_197);
                  const __m256d xi_278 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_242,_mm256_mul_pd(xi_273,xi_273)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_197);
                  const __m256d xi_280 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_214,_mm256_add_pd(_mm256_mul_pd(xi_242,_mm256_mul_pd(xi_268,xi_268)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_197);
                  const __m256d xia_1_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_2_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_3_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_4_collide = _mm256_loadu_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0]);
                  const __m256d xia_5_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_6_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_177 = _mm256_add_pd(xia_2_collide,xia_6_collide);
                  const __m256d xia_7_collide = _mm256_load_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0]);
                  const __m256d xia_8_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_9_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_10_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_218 = _mm256_mul_pd(xia_10_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xia_11_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_178 = _mm256_add_pd(_mm256_add_pd(xi_177,xia_11_collide),xia_1_collide);
                  const __m256d xi_201 = _mm256_mul_pd(xia_11_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_202 = _mm256_add_pd(xi_201,xia_1_collide);
                  const __m256d xia_12_collide = _mm256_loadu_pd(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3 + ctr_0]);
                  const __m256d xi_4 = _mm256_mul_pd(xia_12_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_7 = _mm256_mul_pd(xi_4,_mm256_set_pd(omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a));
                  const __m256d xi_37 = _mm256_mul_pd(xia_12_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_38 = _mm256_mul_pd(xia_12_collide,_mm256_set_pd(xi_33,xi_33,xi_33,xi_33));
                  const __m256d xi_39 = _mm256_add_pd(_mm256_mul_pd(xi_37,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_38);
                  const __m256d xi_41 = _mm256_mul_pd(xia_12_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_46 = _mm256_mul_pd(xia_12_collide,_mm256_set_pd(xi_45,xi_45,xi_45,xi_45));
                  const __m256d xi_53 = _mm256_add_pd(_mm256_mul_pd(xi_38,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_37);
                  const __m256d xia_13_collide = _mm256_load_pd(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + ctr_0]);
                  const __m256d xi_17 = _mm256_mul_pd(xia_13_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_18 = _mm256_mul_pd(xi_17,_mm256_set_pd(omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a));
                  const __m256d xi_32 = _mm256_mul_pd(xia_13_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_34 = _mm256_mul_pd(xia_13_collide,_mm256_set_pd(xi_33,xi_33,xi_33,xi_33));
                  const __m256d xi_35 = _mm256_add_pd(_mm256_mul_pd(xi_34,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_32);
                  const __m256d xi_51 = _mm256_add_pd(_mm256_mul_pd(xi_32,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_34);
                  const __m256d xia_14_collide = _mm256_loadu_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0]);
                  const __m256d xia_15_collide = _mm256_loadu_pd(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3 + ctr_0]);
                  const __m256d xi_26 = _mm256_mul_pd(xia_15_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_28 = _mm256_mul_pd(xi_26,_mm256_set_pd(omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a));
                  const __m256d xi_55 = _mm256_mul_pd(xia_15_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_56 = _mm256_mul_pd(xia_15_collide,_mm256_set_pd(xi_33,xi_33,xi_33,xi_33));
                  const __m256d xi_57 = _mm256_add_pd(_mm256_mul_pd(xi_56,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_55);
                  const __m256d xi_70 = _mm256_add_pd(_mm256_mul_pd(xi_55,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_56);
                  const __m256d xia_16_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_219 = _mm256_add_pd(xi_218,xia_16_collide);
                  const __m256d xia_17_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_174 = _mm256_add_pd(xia_17_collide,xia_9_collide);
                  const __m256d xi_175 = _mm256_add_pd(_mm256_add_pd(xi_174,xia_10_collide),xia_16_collide);
                  const __m256d xia_18_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_19_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_20_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_21_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_206 = _mm256_mul_pd(xia_21_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_207 = _mm256_add_pd(xi_206,xia_5_collide);
                  const __m256d xia_22_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_23_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_24_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0]);
                  const __m256d xia_25_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_170 = _mm256_add_pd(xia_25_collide,xia_3_collide);
                  const __m256d xi_171 = _mm256_add_pd(_mm256_add_pd(xi_170,xia_21_collide),xia_5_collide);
                  const __m256d xia_26_collide = _mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]);
                  const __m256d xi_6 = _mm256_mul_pd(xia_12_collide,_mm256_set_pd(xi_5,xi_5,xi_5,xi_5));
                  const __m256d xi_19 = _mm256_mul_pd(xia_13_collide,_mm256_set_pd(xi_5,xi_5,xi_5,xi_5));
                  const __m256d xi_27 = _mm256_mul_pd(xia_15_collide,_mm256_set_pd(xi_5,xi_5,xi_5,xi_5));
                  const __m256d rho_a_collide = xia_26_collide;
                  const __m256d xi_164 = _mm256_mul_pd(rho_a_collide,_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111));
                  const __m256d xi_180 = _mm256_mul_pd(rho_a_collide,_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_181 = _mm256_add_pd(xi_175,xi_180);
                  const __m256d xi_199 = _mm256_mul_pd(rho_a_collide,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d u_0_a_collide = xia_7_collide;
                  const __m256d xi_0 = _mm256_mul_pd(u_0_a_collide,xia_13_collide);
                  const __m256d xi_13 = _mm256_mul_pd(xi_0,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_14 = _mm256_add_pd(_mm256_mul_pd(xi_0,_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666)),_mm256_mul_pd(xi_13,_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a)));
                  const __m256d xi_20 = _mm256_mul_pd(xi_0,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d xi_42 = _mm256_mul_pd(u_0_a_collide,xi_41);
                  const __m256d xi_47 = _mm256_mul_pd(u_0_a_collide,xi_46);
                  const __m256d xi_58 = _mm256_mul_pd(xi_13,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_61 = _mm256_mul_pd(_mm256_mul_pd(xi_0,_mm256_set_pd(0.041666666666666664,0.041666666666666664,0.041666666666666664,0.041666666666666664)),_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_73 = _mm256_mul_pd(u_0_a_collide,xia_15_collide);
                  const __m256d xi_74 = _mm256_mul_pd(xi_73,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_77 = _mm256_mul_pd(xi_73,_mm256_set_pd(xi_45,xi_45,xi_45,xi_45));
                  const __m256d xi_169 = _mm256_mul_pd(u_0_a_collide,u_0_a_collide);
                  const __m256d xi_173 = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_169),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331)),xi_164);
                  const __m256d xi_182 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_178,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_181,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_19_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_22_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_a_collide,xi_169)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a));
                  const __m256d xi_217 = _mm256_mul_pd(u_0_a_collide,xi_199);
                  const __m256d xi_220 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_17_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_217),xi_219),xia_9_collide);
                  const __m256d xi_221 = _mm256_mul_pd(xi_220,_mm256_set_pd(xi_204,xi_204,xi_204,xi_204));
                  const __m256d xi_222 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_2_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_202),xi_217),xia_6_collide);
                  const __m256d xi_223 = _mm256_mul_pd(xi_222,_mm256_set_pd(xi_204,xi_204,xi_204,xi_204));
                  const __m256d xi_235 = _mm256_mul_pd(xi_222,_mm256_set_pd(xi_234,xi_234,xi_234,xi_234));
                  const __m256d xi_236 = _mm256_mul_pd(xi_235,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_264 = _mm256_mul_pd(xi_220,_mm256_set_pd(xi_234,xi_234,xi_234,xi_234));
                  const __m256d xi_265 = _mm256_mul_pd(xi_264,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d u_1_a_collide = xia_14_collide;
                  const __m256d xi_1 = _mm256_mul_pd(u_1_a_collide,xia_12_collide);
                  const __m256d xi_8 = _mm256_mul_pd(xi_1,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d xi_21 = _mm256_mul_pd(xi_1,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_22 = _mm256_mul_pd(xi_1,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_23 = _mm256_mul_pd(xi_22,_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_24 = _mm256_add_pd(_mm256_mul_pd(xi_21,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_23);
                  const __m256d xi_30 = _mm256_add_pd(xi_14,xi_24);
                  const __m256d xi_43 = _mm256_mul_pd(u_1_a_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_44 = _mm256_mul_pd(xi_43,xia_13_collide);
                  const __m256d xi_48 = _mm256_mul_pd(u_1_a_collide,_mm256_set_pd(xi_45,xi_45,xi_45,xi_45));
                  const __m256d xi_49 = _mm256_mul_pd(xi_48,xia_13_collide);
                  const __m256d xi_50 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_47,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_49,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_42),xi_44);
                  const __m256d xi_52 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_42,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_44,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_47),xi_49);
                  const __m256d xi_59 = _mm256_mul_pd(xi_23,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_63 = _mm256_mul_pd(xi_43,xia_15_collide);
                  const __m256d xi_65 = _mm256_mul_pd(xi_48,xia_15_collide);
                  const __m256d xi_71 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(u_1_a_collide,xia_12_collide),_mm256_set_pd(-0.041666666666666664,-0.041666666666666664,-0.041666666666666664,-0.041666666666666664)),_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_167 = _mm256_mul_pd(u_1_a_collide,u_1_a_collide);
                  const __m256d xi_168 = _mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_167),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_183 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_171,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_178,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_180,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_18_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_23_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_a_collide,xi_167)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a));
                  const __m256d xi_200 = _mm256_mul_pd(u_1_a_collide,xi_199);
                  const __m256d xi_203 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_6_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_200),xi_202),xia_2_collide);
                  const __m256d xi_205 = _mm256_mul_pd(xi_203,_mm256_set_pd(xi_204,xi_204,xi_204,xi_204));
                  const __m256d xi_208 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_25_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_200),xi_207),xia_3_collide);
                  const __m256d xi_209 = _mm256_mul_pd(xi_208,_mm256_set_pd(xi_204,xi_204,xi_204,xi_204));
                  const __m256d xi_237 = _mm256_mul_pd(rho_a_collide,u_1_a_collide);
                  const __m256d xi_239 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_1_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_0_a_collide,xi_237)),xi_177),xi_201),_mm256_set_pd(xi_238,xi_238,xi_238,xi_238));
                  const __m256d xi_240 = _mm256_mul_pd(xi_239,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_245 = _mm256_mul_pd(xi_203,_mm256_set_pd(xi_234,xi_234,xi_234,xi_234));
                  const __m256d xi_252 = _mm256_mul_pd(xi_208,_mm256_set_pd(xi_234,xi_234,xi_234,xi_234));
                  const __m256d xi_259 = _mm256_mul_pd(xi_252,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d u_2_a_collide = xia_4_collide;
                  const __m256d xi_2 = _mm256_mul_pd(u_2_a_collide,xia_15_collide);
                  const __m256d xi_9 = _mm256_mul_pd(xi_2,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_10 = _mm256_mul_pd(xi_2,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_11 = _mm256_mul_pd(xi_10,_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_12 = _mm256_add_pd(_mm256_mul_pd(xi_9,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_11);
                  const __m256d xi_15 = _mm256_add_pd(xi_12,xi_14);
                  const __m256d xi_16 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_8,_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a)),_mm256_mul_pd(_mm256_mul_pd(xi_1,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),xi_15),xi_8);
                  const __m256d xi_25 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_20,_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a)),_mm256_mul_pd(_mm256_mul_pd(xi_0,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),xi_12),xi_20),xi_24);
                  const __m256d xi_29 = _mm256_mul_pd(xi_2,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d xi_31 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_29,_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a)),_mm256_mul_pd(_mm256_mul_pd(xi_2,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),xi_29),xi_30);
                  const __m256d xi_36 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(u_2_a_collide,xia_15_collide),_mm256_set_pd(-0.041666666666666664,-0.041666666666666664,-0.041666666666666664,-0.041666666666666664)),_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_40 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_10,xi_30),xi_36),xi_39);
                  const __m256d xi_54 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_10,xi_30),xi_36),xi_53);
                  const __m256d xi_60 = _mm256_mul_pd(xi_11,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_62 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_21,xi_53),xi_58),xi_59),xi_60),xi_61),xi_9);
                  const __m256d xi_64 = _mm256_mul_pd(u_2_a_collide,xi_41);
                  const __m256d xi_66 = _mm256_mul_pd(u_2_a_collide,xi_46);
                  const __m256d xi_67 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_65,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_66,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_63),xi_64);
                  const __m256d xi_68 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_63,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_64,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_65),xi_66);
                  const __m256d xi_69 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_21,xi_39),xi_58),xi_59),xi_60),xi_61),xi_9);
                  const __m256d xi_72 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_15,xi_22),xi_35),xi_71);
                  const __m256d xi_75 = _mm256_mul_pd(u_2_a_collide,xia_13_collide);
                  const __m256d xi_76 = _mm256_mul_pd(xi_75,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_78 = _mm256_mul_pd(xi_75,_mm256_set_pd(xi_45,xi_45,xi_45,xi_45));
                  const __m256d xi_79 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_77,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_78,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_74),xi_76);
                  const __m256d xi_80 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_74,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_76,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_77),xi_78);
                  const __m256d xi_81 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_15,xi_22),xi_51),xi_71);
                  const __m256d xi_165 = _mm256_mul_pd(u_2_a_collide,u_2_a_collide);
                  const __m256d xi_166 = _mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_165),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_172 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_164,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_166,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_168,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_171,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_169),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))),_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_176 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_166,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_173,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_175,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_167),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))),_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_179 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_168,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_173,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_178,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_165),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))),_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_184 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_171,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_181,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_20_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_8_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_a_collide,xi_165)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a));
                  const __m256d xi_210 = _mm256_mul_pd(xi_179,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_211 = _mm256_mul_pd(xi_172,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_216 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_183,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_210),xi_211),xi_215);
                  const __m256d xi_224 = _mm256_mul_pd(xi_176,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_226 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_182,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_210),xi_224),xi_225);
                  const __m256d xi_227 = _mm256_mul_pd(u_2_a_collide,xi_199);
                  const __m256d xi_228 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_9_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_219),xi_227),xia_17_collide);
                  const __m256d xi_229 = _mm256_mul_pd(xi_228,_mm256_set_pd(xi_204,xi_204,xi_204,xi_204));
                  const __m256d xi_230 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_3_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_207),xi_227),xia_25_collide);
                  const __m256d xi_231 = _mm256_mul_pd(xi_230,_mm256_set_pd(xi_204,xi_204,xi_204,xi_204));
                  const __m256d xi_233 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_184,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_211),xi_224),xi_232);
                  const __m256d xi_244 = _mm256_mul_pd(xi_179,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_246 = _mm256_add_pd(xi_244,xi_245);
                  const __m256d xi_250 = _mm256_add_pd(_mm256_mul_pd(xi_245,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_244);
                  const __m256d xi_253 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_5_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_2_a_collide,xi_237)),xi_170),xi_206),_mm256_set_pd(xi_238,xi_238,xi_238,xi_238));
                  const __m256d xi_256 = _mm256_mul_pd(xi_172,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_257 = _mm256_mul_pd(xi_230,_mm256_set_pd(xi_234,xi_234,xi_234,xi_234));
                  const __m256d xi_258 = _mm256_add_pd(xi_256,xi_257);
                  const __m256d xi_260 = _mm256_mul_pd(xi_253,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_266 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_16_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,u_0_a_collide),u_2_a_collide)),xi_174),xi_218),_mm256_set_pd(xi_238,xi_238,xi_238,xi_238));
                  const __m256d xi_267 = _mm256_mul_pd(xi_266,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_270 = _mm256_mul_pd(xi_176,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_271 = _mm256_mul_pd(xi_228,_mm256_set_pd(xi_234,xi_234,xi_234,xi_234));
                  const __m256d xi_272 = _mm256_add_pd(xi_270,xi_271);
                  const __m256d xi_276 = _mm256_add_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_256);
                  const __m256d xi_279 = _mm256_add_pd(_mm256_mul_pd(xi_271,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_270);
                  const __m256d forceTerm_0_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_0,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_1,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_2,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(xi_0,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_3,xi_3,xi_3,xi_3))),_mm256_mul_pd(_mm256_mul_pd(xi_1,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_3,xi_3,xi_3,xi_3))),_mm256_mul_pd(_mm256_mul_pd(xi_2,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_3,xi_3,xi_3,xi_3))),_mm256_mul_pd(_mm256_mul_pd(u_0_a_collide,xia_13_collide),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),_mm256_mul_pd(_mm256_mul_pd(u_1_a_collide,xia_12_collide),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),_mm256_mul_pd(_mm256_mul_pd(u_2_a_collide,xia_15_collide),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a)));
                  const __m256d forceTerm_1_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_6,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_16),xi_4),xi_7);
                  const __m256d forceTerm_2_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_4,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_7,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_16),xi_6);
                  const __m256d forceTerm_3_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_17,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_18,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_19),xi_25);
                  const __m256d forceTerm_4_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_19,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_17),xi_18),xi_25);
                  const __m256d forceTerm_5_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_27,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_26),xi_28),xi_31);
                  const __m256d forceTerm_6_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_26,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_28,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_27),xi_31);
                  const __m256d forceTerm_7_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_35,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_40,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_50,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_8_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_40,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_51,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_52,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_9_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_35,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_52,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_54,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_10_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_50,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_51,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_54,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_11_a_collide = _mm256_add_pd(_mm256_add_pd(xi_57,xi_62),xi_67);
                  const __m256d forceTerm_12_a_collide = _mm256_add_pd(_mm256_add_pd(xi_57,xi_68),xi_69);
                  const __m256d forceTerm_13_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_70,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_72,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_79,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_14_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_70,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_80,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_81,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_15_a_collide = _mm256_add_pd(_mm256_add_pd(xi_62,xi_68),xi_70);
                  const __m256d forceTerm_16_a_collide = _mm256_add_pd(_mm256_add_pd(xi_67,xi_69),xi_70);
                  const __m256d forceTerm_17_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_57,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_72,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_80,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_18_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_57,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_79,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_81,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d xib_1_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_2_collide = _mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const __m256d xib_3_collide = _mm256_loadu_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0]);
                  const __m256d xib_4_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_304 = _mm256_mul_pd(xib_4_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xib_5_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_309 = _mm256_mul_pd(xib_5_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xib_6_collide = _mm256_load_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0]);
                  const __m256d xib_7_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_287 = _mm256_add_pd(xib_1_collide,xib_7_collide);
                  const __m256d xib_8_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_9_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_10_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0]);
                  const __m256d xib_11_collide = _mm256_load_pd(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + ctr_0]);
                  const __m256d xi_99 = _mm256_mul_pd(xib_11_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_100 = _mm256_mul_pd(xi_99,_mm256_set_pd(omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b));
                  const __m256d xi_114 = _mm256_mul_pd(xib_11_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_116 = _mm256_mul_pd(xib_11_collide,_mm256_set_pd(xi_115,xi_115,xi_115,xi_115));
                  const __m256d xi_117 = _mm256_add_pd(_mm256_mul_pd(xi_116,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_114);
                  const __m256d xi_133 = _mm256_add_pd(_mm256_mul_pd(xi_114,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_116);
                  const __m256d xib_12_collide = _mm256_loadu_pd(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3 + ctr_0]);
                  const __m256d xi_108 = _mm256_mul_pd(xib_12_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_110 = _mm256_mul_pd(xi_108,_mm256_set_pd(omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b));
                  const __m256d xi_137 = _mm256_mul_pd(xib_12_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_138 = _mm256_mul_pd(xib_12_collide,_mm256_set_pd(xi_115,xi_115,xi_115,xi_115));
                  const __m256d xi_139 = _mm256_add_pd(_mm256_mul_pd(xi_138,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_137);
                  const __m256d xi_152 = _mm256_add_pd(_mm256_mul_pd(xi_137,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_138);
                  const __m256d xib_13_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_14_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_15_collide = _mm256_loadu_pd(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3 + ctr_0]);
                  const __m256d xi_86 = _mm256_mul_pd(xib_15_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_89 = _mm256_mul_pd(xi_86,_mm256_set_pd(omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b));
                  const __m256d xi_119 = _mm256_mul_pd(xib_15_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_120 = _mm256_mul_pd(xib_15_collide,_mm256_set_pd(xi_115,xi_115,xi_115,xi_115));
                  const __m256d xi_121 = _mm256_add_pd(_mm256_mul_pd(xi_119,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_120);
                  const __m256d xi_123 = _mm256_mul_pd(xib_15_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_128 = _mm256_mul_pd(xib_15_collide,_mm256_set_pd(xi_127,xi_127,xi_127,xi_127));
                  const __m256d xi_135 = _mm256_add_pd(_mm256_mul_pd(xi_120,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_119);
                  const __m256d xib_16_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_288 = _mm256_add_pd(_mm256_add_pd(xi_287,xib_16_collide),xib_5_collide);
                  const __m256d xi_310 = _mm256_add_pd(xi_309,xib_16_collide);
                  const __m256d xib_17_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_18_collide = _mm256_loadu_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0]);
                  const __m256d xib_19_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_20_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_21_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_305 = _mm256_add_pd(xi_304,xib_21_collide);
                  const __m256d xib_22_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_23_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_291 = _mm256_add_pd(xib_23_collide,xib_8_collide);
                  const __m256d xib_24_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_25_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_292 = _mm256_add_pd(_mm256_add_pd(xi_291,xib_20_collide),xib_25_collide);
                  const __m256d xi_317 = _mm256_mul_pd(xib_25_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_318 = _mm256_add_pd(xi_317,xib_20_collide);
                  const __m256d xib_26_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_294 = _mm256_add_pd(xib_14_collide,xib_26_collide);
                  const __m256d xi_295 = _mm256_add_pd(_mm256_add_pd(xi_294,xib_21_collide),xib_4_collide);
                  const __m256d xi_88 = _mm256_mul_pd(xib_15_collide,_mm256_set_pd(xi_87,xi_87,xi_87,xi_87));
                  const __m256d xi_101 = _mm256_mul_pd(xib_11_collide,_mm256_set_pd(xi_87,xi_87,xi_87,xi_87));
                  const __m256d xi_109 = _mm256_mul_pd(xib_12_collide,_mm256_set_pd(xi_87,xi_87,xi_87,xi_87));
                  const __m256d rho_b_collide = xib_2_collide;
                  const __m256d xi_284 = _mm256_mul_pd(rho_b_collide,_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111));
                  const __m256d xi_297 = _mm256_mul_pd(rho_b_collide,_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_299 = _mm256_add_pd(xi_288,xi_297);
                  const __m256d xi_302 = _mm256_mul_pd(rho_b_collide,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d u_0_b_collide = xib_6_collide;
                  const __m256d xi_82 = _mm256_mul_pd(u_0_b_collide,xib_11_collide);
                  const __m256d xi_95 = _mm256_mul_pd(xi_82,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_96 = _mm256_add_pd(_mm256_mul_pd(xi_82,_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666)),_mm256_mul_pd(xi_95,_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b)));
                  const __m256d xi_102 = _mm256_mul_pd(xi_82,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d xi_124 = _mm256_mul_pd(u_0_b_collide,xi_123);
                  const __m256d xi_129 = _mm256_mul_pd(u_0_b_collide,xi_128);
                  const __m256d xi_140 = _mm256_mul_pd(xi_95,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_143 = _mm256_mul_pd(_mm256_mul_pd(xi_82,_mm256_set_pd(0.041666666666666664,0.041666666666666664,0.041666666666666664,0.041666666666666664)),_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_155 = _mm256_mul_pd(u_0_b_collide,xib_12_collide);
                  const __m256d xi_156 = _mm256_mul_pd(xi_155,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_159 = _mm256_mul_pd(xi_155,_mm256_set_pd(xi_127,xi_127,xi_127,xi_127));
                  const __m256d xi_283 = _mm256_mul_pd(u_0_b_collide,u_0_b_collide);
                  const __m256d xi_290 = _mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_283),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_298 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_292,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_295,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_297,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_17_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_22_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_b_collide,xi_283)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b));
                  const __m256d xi_316 = _mm256_mul_pd(u_0_b_collide,xi_302);
                  const __m256d xi_319 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_8_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_316),xi_318),xib_23_collide);
                  const __m256d xi_320 = _mm256_mul_pd(xi_319,_mm256_set_pd(xi_307,xi_307,xi_307,xi_307));
                  const __m256d xi_321 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_14_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_305),xi_316),xib_26_collide);
                  const __m256d xi_322 = _mm256_mul_pd(xi_321,_mm256_set_pd(xi_307,xi_307,xi_307,xi_307));
                  const __m256d xi_332 = _mm256_mul_pd(xi_321,_mm256_set_pd(xi_331,xi_331,xi_331,xi_331));
                  const __m256d xi_333 = _mm256_mul_pd(xi_332,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_349 = _mm256_mul_pd(xi_319,_mm256_set_pd(xi_331,xi_331,xi_331,xi_331));
                  const __m256d xi_350 = _mm256_mul_pd(xi_349,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d u_1_b_collide = xib_18_collide;
                  const __m256d xi_83 = _mm256_mul_pd(u_1_b_collide,xib_15_collide);
                  const __m256d xi_90 = _mm256_mul_pd(xi_83,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d xi_103 = _mm256_mul_pd(xi_83,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_104 = _mm256_mul_pd(xi_83,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_105 = _mm256_mul_pd(xi_104,_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_106 = _mm256_add_pd(_mm256_mul_pd(xi_103,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_105);
                  const __m256d xi_112 = _mm256_add_pd(xi_106,xi_96);
                  const __m256d xi_125 = _mm256_mul_pd(u_1_b_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_126 = _mm256_mul_pd(xi_125,xib_11_collide);
                  const __m256d xi_130 = _mm256_mul_pd(u_1_b_collide,_mm256_set_pd(xi_127,xi_127,xi_127,xi_127));
                  const __m256d xi_131 = _mm256_mul_pd(xi_130,xib_11_collide);
                  const __m256d xi_132 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_129,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_131,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_124),xi_126);
                  const __m256d xi_134 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_124,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_126,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_129),xi_131);
                  const __m256d xi_141 = _mm256_mul_pd(xi_105,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_145 = _mm256_mul_pd(xi_125,xib_12_collide);
                  const __m256d xi_147 = _mm256_mul_pd(xi_130,xib_12_collide);
                  const __m256d xi_153 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(u_1_b_collide,xib_15_collide),_mm256_set_pd(-0.041666666666666664,-0.041666666666666664,-0.041666666666666664,-0.041666666666666664)),_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_285 = _mm256_mul_pd(u_1_b_collide,u_1_b_collide);
                  const __m256d xi_286 = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_285),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331)),xi_284);
                  const __m256d xi_300 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_295,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_299,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_24_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_9_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_b_collide,xi_285)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b));
                  const __m256d xi_303 = _mm256_mul_pd(u_1_b_collide,xi_302);
                  const __m256d xi_306 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_26_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_303),xi_305),xib_14_collide);
                  const __m256d xi_308 = _mm256_mul_pd(xi_306,_mm256_set_pd(xi_307,xi_307,xi_307,xi_307));
                  const __m256d xi_311 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_7_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_303),xi_310),xib_1_collide);
                  const __m256d xi_312 = _mm256_mul_pd(xi_311,_mm256_set_pd(xi_307,xi_307,xi_307,xi_307));
                  const __m256d xi_334 = _mm256_mul_pd(rho_b_collide,u_1_b_collide);
                  const __m256d xi_336 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_21_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_0_b_collide,xi_334)),xi_294),xi_304),_mm256_set_pd(xi_335,xi_335,xi_335,xi_335));
                  const __m256d xi_337 = _mm256_mul_pd(xi_336,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_339 = _mm256_mul_pd(xi_306,_mm256_set_pd(xi_331,xi_331,xi_331,xi_331));
                  const __m256d xi_342 = _mm256_mul_pd(xi_311,_mm256_set_pd(xi_331,xi_331,xi_331,xi_331));
                  const __m256d xi_347 = _mm256_mul_pd(xi_342,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d u_2_b_collide = xib_3_collide;
                  const __m256d xi_84 = _mm256_mul_pd(u_2_b_collide,xib_12_collide);
                  const __m256d xi_91 = _mm256_mul_pd(xi_84,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_92 = _mm256_mul_pd(xi_84,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_93 = _mm256_mul_pd(xi_92,_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_94 = _mm256_add_pd(_mm256_mul_pd(xi_91,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_93);
                  const __m256d xi_97 = _mm256_add_pd(xi_94,xi_96);
                  const __m256d xi_98 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_90,_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b)),_mm256_mul_pd(_mm256_mul_pd(xi_83,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),xi_90),xi_97);
                  const __m256d xi_107 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_102,_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b)),_mm256_mul_pd(_mm256_mul_pd(xi_82,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),xi_102),xi_106),xi_94);
                  const __m256d xi_111 = _mm256_mul_pd(xi_84,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d xi_113 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_111,_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b)),_mm256_mul_pd(_mm256_mul_pd(xi_84,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),xi_111),xi_112);
                  const __m256d xi_118 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(u_2_b_collide,xib_12_collide),_mm256_set_pd(-0.041666666666666664,-0.041666666666666664,-0.041666666666666664,-0.041666666666666664)),_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_122 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_112,xi_118),xi_121),xi_92);
                  const __m256d xi_136 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_112,xi_118),xi_135),xi_92);
                  const __m256d xi_142 = _mm256_mul_pd(xi_93,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_144 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_103,xi_135),xi_140),xi_141),xi_142),xi_143),xi_91);
                  const __m256d xi_146 = _mm256_mul_pd(u_2_b_collide,xi_123);
                  const __m256d xi_148 = _mm256_mul_pd(u_2_b_collide,xi_128);
                  const __m256d xi_149 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_147,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_148,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_145),xi_146);
                  const __m256d xi_150 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_145,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_146,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_147),xi_148);
                  const __m256d xi_151 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_103,xi_121),xi_140),xi_141),xi_142),xi_143),xi_91);
                  const __m256d xi_154 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_104,xi_117),xi_153),xi_97);
                  const __m256d xi_157 = _mm256_mul_pd(u_2_b_collide,xib_11_collide);
                  const __m256d xi_158 = _mm256_mul_pd(xi_157,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_160 = _mm256_mul_pd(xi_157,_mm256_set_pd(xi_127,xi_127,xi_127,xi_127));
                  const __m256d xi_161 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_159,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_160,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_156),xi_158);
                  const __m256d xi_162 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_156,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_158,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_159),xi_160);
                  const __m256d xi_163 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_104,xi_133),xi_153),xi_97);
                  const __m256d xi_281 = _mm256_mul_pd(u_2_b_collide,u_2_b_collide);
                  const __m256d xi_282 = _mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_281),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_289 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_282,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_286,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_288,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_283),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))),_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_293 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_282,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_284,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_290,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_292,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_285),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))),_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_296 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_286,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_290,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_295,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_281),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))),_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_301 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_292,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_299,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_13_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_19_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_b_collide,xi_281)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b));
                  const __m256d xi_313 = _mm256_mul_pd(xi_296,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_314 = _mm256_mul_pd(xi_289,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_315 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_300,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_215),xi_313),xi_314);
                  const __m256d xi_323 = _mm256_mul_pd(xi_293,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_324 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_298,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_225),xi_313),xi_323);
                  const __m256d xi_325 = _mm256_mul_pd(u_2_b_collide,xi_302);
                  const __m256d xi_326 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_1_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_310),xi_325),xib_7_collide);
                  const __m256d xi_327 = _mm256_mul_pd(xi_326,_mm256_set_pd(xi_307,xi_307,xi_307,xi_307));
                  const __m256d xi_328 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_23_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_318),xi_325),xib_8_collide);
                  const __m256d xi_329 = _mm256_mul_pd(xi_328,_mm256_set_pd(xi_307,xi_307,xi_307,xi_307));
                  const __m256d xi_330 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_301,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_232),xi_314),xi_323);
                  const __m256d xi_338 = _mm256_mul_pd(xi_296,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_340 = _mm256_add_pd(xi_338,xi_339);
                  const __m256d xi_341 = _mm256_add_pd(_mm256_mul_pd(xi_339,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_338);
                  const __m256d xi_343 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_16_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_2_b_collide,xi_334)),xi_287),xi_309),_mm256_set_pd(xi_335,xi_335,xi_335,xi_335));
                  const __m256d xi_344 = _mm256_mul_pd(xi_289,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_345 = _mm256_mul_pd(xi_326,_mm256_set_pd(xi_331,xi_331,xi_331,xi_331));
                  const __m256d xi_346 = _mm256_add_pd(xi_344,xi_345);
                  const __m256d xi_348 = _mm256_mul_pd(xi_343,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_351 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_20_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,u_0_b_collide),u_2_b_collide)),xi_291),xi_317),_mm256_set_pd(xi_335,xi_335,xi_335,xi_335));
                  const __m256d xi_352 = _mm256_mul_pd(xi_351,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_353 = _mm256_mul_pd(xi_293,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_354 = _mm256_mul_pd(xi_328,_mm256_set_pd(xi_331,xi_331,xi_331,xi_331));
                  const __m256d xi_355 = _mm256_add_pd(xi_353,xi_354);
                  const __m256d xi_356 = _mm256_add_pd(_mm256_mul_pd(xi_345,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_344);
                  const __m256d xi_357 = _mm256_add_pd(_mm256_mul_pd(xi_354,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_353);
                  const __m256d forceTerm_0_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_82,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_83,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_84,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(xi_82,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_85,xi_85,xi_85,xi_85))),_mm256_mul_pd(_mm256_mul_pd(xi_83,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_85,xi_85,xi_85,xi_85))),_mm256_mul_pd(_mm256_mul_pd(xi_84,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_85,xi_85,xi_85,xi_85))),_mm256_mul_pd(_mm256_mul_pd(u_0_b_collide,xib_11_collide),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),_mm256_mul_pd(_mm256_mul_pd(u_1_b_collide,xib_15_collide),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),_mm256_mul_pd(_mm256_mul_pd(u_2_b_collide,xib_12_collide),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b)));
                  const __m256d forceTerm_1_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_88,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_86),xi_89),xi_98);
                  const __m256d forceTerm_2_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_86,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_89,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_88),xi_98);
                  const __m256d forceTerm_3_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_100,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_99,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_101),xi_107);
                  const __m256d forceTerm_4_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_101,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_100),xi_107),xi_99);
                  const __m256d forceTerm_5_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_109,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_108),xi_110),xi_113);
                  const __m256d forceTerm_6_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_108,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_110,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_109),xi_113);
                  const __m256d forceTerm_7_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_117,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_122,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_132,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_8_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_122,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_133,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_134,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_9_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_117,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_134,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_136,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_10_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_132,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_133,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_136,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_11_b_collide = _mm256_add_pd(_mm256_add_pd(xi_139,xi_144),xi_149);
                  const __m256d forceTerm_12_b_collide = _mm256_add_pd(_mm256_add_pd(xi_139,xi_150),xi_151);
                  const __m256d forceTerm_13_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_152,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_154,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_161,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_14_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_152,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_162,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_163,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_15_b_collide = _mm256_add_pd(_mm256_add_pd(xi_144,xi_150),xi_152);
                  const __m256d forceTerm_16_b_collide = _mm256_add_pd(_mm256_add_pd(xi_149,xi_151),xi_152);
                  const __m256d forceTerm_17_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_139,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_154,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_162,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_18_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_139,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_161,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_163,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d tmp_a0 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_182,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_183,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_184,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_0_a_collide),xi_172),xi_176),xi_179),xi_198),xia_24_collide);
                  const __m256d tmp_a1 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_205,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_209,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_1_a_collide),xi_216),xia_23_collide);
                  const __m256d tmp_a2 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_2_a_collide,xi_205),xi_209),xi_216),xia_18_collide);
                  const __m256d tmp_a3 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_3_a_collide,xi_221),xi_223),xi_226),xia_19_collide);
                  const __m256d tmp_a4 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_221,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_223,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_4_a_collide),xi_226),xia_22_collide);
                  const __m256d tmp_a5 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_229,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_231,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_5_a_collide),xi_233),xia_8_collide);
                  const __m256d tmp_a6 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_6_a_collide,xi_229),xi_231),xi_233),xia_20_collide);
                  const __m256d tmp_a7 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_7_a_collide,xi_236),xi_240),xi_243),xi_246),xia_6_collide);
                  const __m256d tmp_a8 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_8_a_collide,xi_235),xi_239),xi_246),xi_248),xia_11_collide);
                  const __m256d tmp_a9 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_9_a_collide,xi_236),xi_239),xi_249),xi_250),xia_1_collide);
                  const __m256d tmp_a10 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_10_a_collide,xi_235),xi_240),xi_250),xi_251),xia_2_collide);
                  const __m256d tmp_a11 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_11_a_collide,xi_252),xi_253),xi_255),xi_258),xia_21_collide);
                  const __m256d tmp_a12 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_12_a_collide,xi_258),xi_259),xi_260),xi_263),xia_3_collide);
                  const __m256d tmp_a13 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_13_a_collide,xi_265),xi_267),xi_269),xi_272),xia_9_collide);
                  const __m256d tmp_a14 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_14_a_collide,xi_264),xi_266),xi_272),xi_274),xia_10_collide);
                  const __m256d tmp_a15 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_15_a_collide,xi_252),xi_260),xi_275),xi_276),xia_25_collide);
                  const __m256d tmp_a16 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_16_a_collide,xi_253),xi_259),xi_276),xi_277),xia_5_collide);
                  const __m256d tmp_a17 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_17_a_collide,xi_265),xi_266),xi_278),xi_279),xia_16_collide);
                  const __m256d tmp_a18 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_18_a_collide,xi_264),xi_267),xi_279),xi_280),xia_17_collide);
                  const __m256d tmp_b0 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_298,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_300,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_301,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_0_b_collide),xi_198),xi_289),xi_293),xi_296),xib_10_collide);
                  const __m256d tmp_b1 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_308,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_312,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_1_b_collide),xi_315),xib_9_collide);
                  const __m256d tmp_b2 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_2_b_collide,xi_308),xi_312),xi_315),xib_24_collide);
                  const __m256d tmp_b3 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_3_b_collide,xi_320),xi_322),xi_324),xib_17_collide);
                  const __m256d tmp_b4 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_320,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_322,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_4_b_collide),xi_324),xib_22_collide);
                  const __m256d tmp_b5 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_327,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_329,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_5_b_collide),xi_330),xib_19_collide);
                  const __m256d tmp_b6 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_6_b_collide,xi_327),xi_329),xi_330),xib_13_collide);
                  const __m256d tmp_b7 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_7_b_collide,xi_243),xi_333),xi_337),xi_340),xib_26_collide);
                  const __m256d tmp_b8 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_8_b_collide,xi_248),xi_332),xi_336),xi_340),xib_4_collide);
                  const __m256d tmp_b9 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_9_b_collide,xi_249),xi_333),xi_336),xi_341),xib_21_collide);
                  const __m256d tmp_b10 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_10_b_collide,xi_251),xi_332),xi_337),xi_341),xib_14_collide);
                  const __m256d tmp_b11 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_11_b_collide,xi_255),xi_342),xi_343),xi_346),xib_5_collide);
                  const __m256d tmp_b12 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_12_b_collide,xi_263),xi_346),xi_347),xi_348),xib_1_collide);
                  const __m256d tmp_b13 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_13_b_collide,xi_269),xi_350),xi_352),xi_355),xib_23_collide);
                  const __m256d tmp_b14 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_14_b_collide,xi_274),xi_349),xi_351),xi_355),xib_25_collide);
                  const __m256d tmp_b15 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_15_b_collide,xi_275),xi_342),xi_348),xi_356),xib_7_collide);
                  const __m256d tmp_b16 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_16_b_collide,xi_277),xi_343),xi_347),xi_356),xib_16_collide);
                  const __m256d tmp_b17 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_17_b_collide,xi_278),xi_350),xi_351),xi_357),xib_20_collide);
                  const __m256d tmp_b18 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_18_b_collide,xi_280),xi_349),xi_352),xi_357),xib_8_collide);
                  const __m256d xirecolor_0 = _mm256_add_pd(tmp_a0,tmp_b0);
                  const __m256d xirecolor_1 = _mm256_add_pd(_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]),_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_2 = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),xirecolor_1);
                  const __m256d xi_364 = _mm256_mul_pd(xirecolor_2,_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_3 = _mm256_mul_pd(xirecolor_2,_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_4 = _mm256_add_pd(tmp_a1,tmp_b1);
                  const __m256d xirecolor_5 = xi_189;
                  const __m256d xirecolor_6 = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),xirecolor_5);
                  const __m256d xirecolor_7 = _mm256_mul_pd(xirecolor_6,_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]));
                  const __m256d xirecolor_8 = _mm256_cmp_pd(xirecolor_5,_mm256_set_pd(0.0,0.0,0.0,0.0),_CMP_NLE_UQ);
                  const __m256d xirecolor_9 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set_pd(beta,beta,beta,beta),_mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),_mm256_mul_pd(xirecolor_1,xirecolor_1))),_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_10 = _mm256_mul_pd(xirecolor_9,_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]))));
                  const __m256d xirecolor_11 = _mm256_mul_pd(xirecolor_10,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xirecolor_7,xirecolor_8));
                  const __m256d xirecolor_12 = _mm256_add_pd(tmp_a2,tmp_b2);
                  const __m256d xirecolor_13 = _mm256_mul_pd(xirecolor_10,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xirecolor_7,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_8));
                  const __m256d xirecolor_14 = _mm256_add_pd(tmp_a3,tmp_b3);
                  const __m256d xirecolor_15 = _mm256_mul_pd(xirecolor_6,_mm256_load_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_16 = _mm256_mul_pd(xirecolor_10,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xirecolor_15,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_8));
                  const __m256d xirecolor_17 = _mm256_add_pd(tmp_a4,tmp_b4);
                  const __m256d xirecolor_18 = _mm256_mul_pd(xirecolor_10,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xirecolor_15,xirecolor_8));
                  const __m256d xirecolor_19 = _mm256_add_pd(tmp_a5,tmp_b5);
                  const __m256d xirecolor_20 = _mm256_mul_pd(xirecolor_6,_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]));
                  const __m256d xirecolor_21 = _mm256_mul_pd(xirecolor_10,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xirecolor_20,xirecolor_8));
                  const __m256d xirecolor_22 = _mm256_add_pd(tmp_a6,tmp_b6);
                  const __m256d xirecolor_23 = _mm256_mul_pd(xirecolor_10,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xirecolor_20,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_8));
                  const __m256d xirecolor_24 = _mm256_add_pd(tmp_a7,tmp_b7);
                  const __m256d xirecolor_25 = xi_241;
                  const __m256d xirecolor_26 = _mm256_mul_pd(xirecolor_6,_mm256_set_pd(0.70710678118654757,0.70710678118654757,0.70710678118654757,0.70710678118654757));
                  const __m256d xi_358 = _mm256_mul_pd(xirecolor_25,xirecolor_26);
                  const __m256d xirecolor_27 = _mm256_mul_pd(xirecolor_9,_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]))));
                  const __m256d xirecolor_28 = _mm256_mul_pd(xirecolor_27,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_358,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_8));
                  const __m256d xirecolor_29 = _mm256_add_pd(tmp_a8,tmp_b8);
                  const __m256d xirecolor_30 = xi_247;
                  const __m256d xi_359 = _mm256_mul_pd(xirecolor_26,xirecolor_30);
                  const __m256d xirecolor_31 = _mm256_mul_pd(xirecolor_27,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_359,xirecolor_8));
                  const __m256d xirecolor_32 = _mm256_add_pd(tmp_a9,tmp_b9);
                  const __m256d xirecolor_33 = _mm256_mul_pd(xirecolor_27,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_359,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_8));
                  const __m256d xirecolor_34 = _mm256_add_pd(tmp_a10,tmp_b10);
                  const __m256d xirecolor_35 = _mm256_mul_pd(xirecolor_27,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_358,xirecolor_8));
                  const __m256d xirecolor_36 = _mm256_add_pd(tmp_a11,tmp_b11);
                  const __m256d xirecolor_37 = xi_254;
                  const __m256d xi_360 = _mm256_mul_pd(xirecolor_26,xirecolor_37);
                  const __m256d xirecolor_38 = _mm256_mul_pd(xirecolor_27,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_360,xirecolor_8));
                  const __m256d xirecolor_39 = _mm256_add_pd(tmp_a12,tmp_b12);
                  const __m256d xirecolor_40 = xi_261;
                  const __m256d xirecolor_41 = _mm256_add_pd(xirecolor_40,_mm256_loadu_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]));
                  const __m256d xi_361 = _mm256_mul_pd(xirecolor_26,xirecolor_41);
                  const __m256d xirecolor_42 = _mm256_mul_pd(xirecolor_27,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_361,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_8));
                  const __m256d xirecolor_43 = _mm256_add_pd(tmp_a13,tmp_b13);
                  const __m256d xirecolor_44 = _mm256_add_pd(xirecolor_40,_mm256_load_pd(& _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]));
                  const __m256d xi_362 = _mm256_mul_pd(xirecolor_26,xirecolor_44);
                  const __m256d xirecolor_45 = _mm256_mul_pd(xirecolor_27,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_362,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_8));
                  const __m256d xirecolor_46 = _mm256_add_pd(tmp_a14,tmp_b14);
                  const __m256d xirecolor_47 = xi_273;
                  const __m256d xi_363 = _mm256_mul_pd(xirecolor_26,xirecolor_47);
                  const __m256d xirecolor_48 = _mm256_mul_pd(xirecolor_27,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_363,xirecolor_8));
                  const __m256d xirecolor_49 = _mm256_add_pd(tmp_a15,tmp_b15);
                  const __m256d xirecolor_50 = _mm256_mul_pd(xirecolor_27,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_361,xirecolor_8));
                  const __m256d xirecolor_51 = _mm256_add_pd(tmp_a16,tmp_b16);
                  const __m256d xirecolor_52 = _mm256_mul_pd(xirecolor_27,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_360,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_8));
                  const __m256d xirecolor_53 = _mm256_add_pd(tmp_a17,tmp_b17);
                  const __m256d xirecolor_54 = _mm256_mul_pd(xirecolor_27,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_363,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_8));
                  const __m256d xirecolor_55 = _mm256_add_pd(tmp_a18,tmp_b18);
                  const __m256d xirecolor_56 = _mm256_mul_pd(xirecolor_27,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_362,xirecolor_8));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0],_mm256_mul_pd(xirecolor_0,xirecolor_3));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_4),xirecolor_11));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_12,xirecolor_3),xirecolor_13));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_14,xirecolor_3),xirecolor_16));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_17,xirecolor_3),xirecolor_18));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_19,xirecolor_3),xirecolor_21));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_22,xirecolor_3),xirecolor_23));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_24,xirecolor_3),xirecolor_28));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_29,xirecolor_3),xirecolor_31));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_32),xirecolor_33));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_34),xirecolor_35));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_36),xirecolor_38));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_39),xirecolor_42));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_43),xirecolor_45));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_46),xirecolor_48));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_49),xirecolor_50));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_51),xirecolor_52));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_53),xirecolor_54));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_55),xirecolor_56));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0],_mm256_mul_pd(xi_364,xirecolor_0));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_11,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_4)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_13,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_12)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_16,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_14)));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_18,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_17)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_21,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_19)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_23,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_22)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_28,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_24)));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_31,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_29)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_33,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_32)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_35,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_34)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_38,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_36)));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_42,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_39)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_45,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_43)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_48,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_46)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_50,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_49)));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_52,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_51)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_54,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_53)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_56,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_364,xirecolor_55)));
               }
               for (int64_t ctr_0 = (int64_t)((_size_color_gradient_0) / (4)) * (4); ctr_0 < _size_color_gradient_0; ctr_0 += 1)
               {
                  const double xi_185 = (_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]*_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0]);
                  const double xi_186 = (_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]*_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0]);
                  const double xi_187 = (_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]*_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0]);
                  const double xi_188 = xi_185 + xi_186 + xi_187;
                  const double xi_189 = pow(xi_188, 0.5);
                  const double xi_194 = (_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]);
                  const double xi_196 = sigma*xi_189*((0.5 < _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (omega_shear_a): ((-0.5 > _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (omega_shear_b): ((0.0 < _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (xi_191 + xi_193*xi_194 - xi_193*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]): (xi_191 + xi_194*xi_195 + xi_195*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))));
                  const bool xi_197 = xi_189 > 0.0;
                  const double xi_198 = ((xi_197) ? (xi_196*0.25): (0.0));
                  const double xi_212 = ((1.0) / (xi_188));
                  const double xi_213 = xi_212*0.055555555555555552;
                  const double xi_214 = xi_196*1.125;
                  const double xi_215 = ((xi_197) ? (xi_214*(xi_186*xi_213 - 0.018518518518518517)): (0.0));
                  const double xi_225 = ((xi_197) ? (xi_214*(xi_185*xi_213 - 0.018518518518518517)): (0.0));
                  const double xi_232 = ((xi_197) ? (xi_214*(xi_187*xi_213 - 0.018518518518518517)): (0.0));
                  const double xi_241 = -_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0] + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0];
                  const double xi_242 = xi_212*0.027777777777777776;
                  const double xi_243 = ((xi_197) ? (xi_214*(xi_242*(xi_241*xi_241) - 0.037037037037037035)): (0.0));
                  const double xi_247 = _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0] + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0];
                  const double xi_248 = ((xi_197) ? (xi_214*(xi_242*(xi_247*xi_247) - 0.037037037037037035)): (0.0));
                  const double xi_249 = ((xi_197) ? (xi_214*(xi_242*(xi_247*xi_247) - 0.037037037037037035)): (0.0));
                  const double xi_251 = ((xi_197) ? (xi_214*(xi_242*(xi_241*xi_241) - 0.037037037037037035)): (0.0));
                  const double xi_254 = _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0] + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0];
                  const double xi_255 = ((xi_197) ? (xi_214*(xi_242*(xi_254*xi_254) - 0.037037037037037035)): (0.0));
                  const double xi_261 = -_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0];
                  const double xi_262 = xi_261 + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0];
                  const double xi_263 = ((xi_197) ? (xi_214*(xi_242*(xi_262*xi_262) - 0.037037037037037035)): (0.0));
                  const double xi_268 = xi_261 + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0];
                  const double xi_269 = ((xi_197) ? (xi_214*(xi_242*(xi_268*xi_268) - 0.037037037037037035)): (0.0));
                  const double xi_273 = _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0] + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0];
                  const double xi_274 = ((xi_197) ? (xi_214*(xi_242*(xi_273*xi_273) - 0.037037037037037035)): (0.0));
                  const double xi_275 = ((xi_197) ? (xi_214*(xi_242*(xi_262*xi_262) - 0.037037037037037035)): (0.0));
                  const double xi_277 = ((xi_197) ? (xi_214*(xi_242*(xi_254*xi_254) - 0.037037037037037035)): (0.0));
                  const double xi_278 = ((xi_197) ? (xi_214*(xi_242*(xi_273*xi_273) - 0.037037037037037035)): (0.0));
                  const double xi_280 = ((xi_197) ? (xi_214*(xi_242*(xi_268*xi_268) - 0.037037037037037035)): (0.0));
                  const double xia_1_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0];
                  const double xia_2_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0];
                  const double xia_3_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0];
                  const double xia_4_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0];
                  const double xia_5_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0];
                  const double xia_6_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0];
                  const double xi_177 = xia_2_collide + xia_6_collide;
                  const double xia_7_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0];
                  const double xia_8_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0];
                  const double xia_9_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0];
                  const double xia_10_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0];
                  const double xi_218 = -xia_10_collide;
                  const double xia_11_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0];
                  const double xi_178 = xi_177 + xia_11_collide + xia_1_collide;
                  const double xi_201 = -xia_11_collide;
                  const double xi_202 = xi_201 + xia_1_collide;
                  const double xia_12_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3 + ctr_0];
                  const double xi_4 = xia_12_collide*0.16666666666666666;
                  const double xi_7 = omega_odd_a*xi_4;
                  const double xi_37 = xia_12_collide*0.083333333333333329;
                  const double xi_38 = xi_33*xia_12_collide;
                  const double xi_39 = -xi_37 + xi_38;
                  const double xi_41 = xia_12_collide*0.25;
                  const double xi_46 = xi_45*xia_12_collide;
                  const double xi_53 = xi_37 - xi_38;
                  const double xia_13_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + ctr_0];
                  const double xi_17 = xia_13_collide*0.16666666666666666;
                  const double xi_18 = omega_odd_a*xi_17;
                  const double xi_32 = xia_13_collide*0.083333333333333329;
                  const double xi_34 = xi_33*xia_13_collide;
                  const double xi_35 = xi_32 - xi_34;
                  const double xi_51 = -xi_32 + xi_34;
                  const double xia_14_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0];
                  const double xia_15_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3 + ctr_0];
                  const double xi_26 = xia_15_collide*0.16666666666666666;
                  const double xi_28 = omega_odd_a*xi_26;
                  const double xi_55 = xia_15_collide*0.083333333333333329;
                  const double xi_56 = xi_33*xia_15_collide;
                  const double xi_57 = xi_55 - xi_56;
                  const double xi_70 = -xi_55 + xi_56;
                  const double xia_16_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0];
                  const double xi_219 = xi_218 + xia_16_collide;
                  const double xia_17_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0];
                  const double xi_174 = xia_17_collide + xia_9_collide;
                  const double xi_175 = xi_174 + xia_10_collide + xia_16_collide;
                  const double xia_18_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0];
                  const double xia_19_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0];
                  const double xia_20_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0];
                  const double xia_21_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0];
                  const double xi_206 = -xia_21_collide;
                  const double xi_207 = xi_206 + xia_5_collide;
                  const double xia_22_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0];
                  const double xia_23_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0];
                  const double xia_24_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0];
                  const double xia_25_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0];
                  const double xi_170 = xia_25_collide + xia_3_collide;
                  const double xi_171 = xi_170 + xia_21_collide + xia_5_collide;
                  const double xia_26_collide = _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0];
                  const double xi_6 = xi_5*xia_12_collide;
                  const double xi_19 = xi_5*xia_13_collide;
                  const double xi_27 = xi_5*xia_15_collide;
                  const double rho_a_collide = xia_26_collide;
                  const double xi_164 = rho_a_collide*-0.1111111111111111;
                  const double xi_180 = rho_a_collide*-0.33333333333333331;
                  const double xi_181 = xi_175 + xi_180;
                  const double xi_199 = rho_a_collide*0.33333333333333331;
                  const double u_0_a_collide = xia_7_collide;
                  const double xi_0 = u_0_a_collide*xia_13_collide;
                  const double xi_13 = xi_0*0.083333333333333329;
                  const double xi_14 = omega_even_a*xi_13 + xi_0*-0.16666666666666666;
                  const double xi_20 = xi_0*0.33333333333333331;
                  const double xi_42 = u_0_a_collide*xi_41;
                  const double xi_47 = u_0_a_collide*xi_46;
                  const double xi_58 = -xi_13;
                  const double xi_61 = omega_even_a*xi_0*0.041666666666666664;
                  const double xi_73 = u_0_a_collide*xia_15_collide;
                  const double xi_74 = xi_73*0.25;
                  const double xi_77 = xi_45*xi_73;
                  const double xi_169 = (u_0_a_collide*u_0_a_collide);
                  const double xi_173 = rho_a_collide*xi_169*-0.33333333333333331 + xi_164;
                  const double xi_182 = omega_shear_a*(rho_a_collide*xi_169 - xi_178 - xi_181 - xia_19_collide - xia_22_collide);
                  const double xi_217 = u_0_a_collide*xi_199;
                  const double xi_220 = xi_217 + xi_219 - xia_17_collide + xia_9_collide;
                  const double xi_221 = xi_204*xi_220;
                  const double xi_222 = xi_202 + xi_217 - xia_2_collide + xia_6_collide;
                  const double xi_223 = xi_204*xi_222;
                  const double xi_235 = xi_222*xi_234;
                  const double xi_236 = -xi_235;
                  const double xi_264 = xi_220*xi_234;
                  const double xi_265 = -xi_264;
                  const double u_1_a_collide = xia_14_collide;
                  const double xi_1 = u_1_a_collide*xia_12_collide;
                  const double xi_8 = xi_1*0.33333333333333331;
                  const double xi_21 = xi_1*0.16666666666666666;
                  const double xi_22 = xi_1*0.083333333333333329;
                  const double xi_23 = omega_even_a*xi_22;
                  const double xi_24 = -xi_21 + xi_23;
                  const double xi_30 = xi_14 + xi_24;
                  const double xi_43 = u_1_a_collide*0.25;
                  const double xi_44 = xi_43*xia_13_collide;
                  const double xi_48 = u_1_a_collide*xi_45;
                  const double xi_49 = xi_48*xia_13_collide;
                  const double xi_50 = xi_42 + xi_44 - xi_47 - xi_49;
                  const double xi_52 = -xi_42 - xi_44 + xi_47 + xi_49;
                  const double xi_59 = -xi_23;
                  const double xi_63 = xi_43*xia_15_collide;
                  const double xi_65 = xi_48*xia_15_collide;
                  const double xi_71 = omega_even_a*u_1_a_collide*xia_12_collide*-0.041666666666666664;
                  const double xi_167 = (u_1_a_collide*u_1_a_collide);
                  const double xi_168 = rho_a_collide*xi_167*-0.33333333333333331;
                  const double xi_183 = omega_shear_a*(rho_a_collide*xi_167 - xi_171 - xi_178 - xi_180 - xia_18_collide - xia_23_collide);
                  const double xi_200 = u_1_a_collide*xi_199;
                  const double xi_203 = xi_200 + xi_202 + xia_2_collide - xia_6_collide;
                  const double xi_205 = xi_203*xi_204;
                  const double xi_208 = xi_200 + xi_207 - xia_25_collide + xia_3_collide;
                  const double xi_209 = xi_204*xi_208;
                  const double xi_237 = rho_a_collide*u_1_a_collide;
                  const double xi_239 = xi_238*(u_0_a_collide*xi_237 + xi_177 + xi_201 - xia_1_collide);
                  const double xi_240 = -xi_239;
                  const double xi_245 = xi_203*xi_234;
                  const double xi_252 = xi_208*xi_234;
                  const double xi_259 = -xi_252;
                  const double u_2_a_collide = xia_4_collide;
                  const double xi_2 = u_2_a_collide*xia_15_collide;
                  const double xi_9 = xi_2*0.16666666666666666;
                  const double xi_10 = xi_2*0.083333333333333329;
                  const double xi_11 = omega_even_a*xi_10;
                  const double xi_12 = xi_11 - xi_9;
                  const double xi_15 = xi_12 + xi_14;
                  const double xi_16 = omega_even_a*xi_8 + omega_shear_a*xi_1*-0.5 + xi_15 + xi_8;
                  const double xi_25 = omega_even_a*xi_20 + omega_shear_a*xi_0*-0.5 + xi_12 + xi_20 + xi_24;
                  const double xi_29 = xi_2*0.33333333333333331;
                  const double xi_31 = omega_even_a*xi_29 + omega_shear_a*xi_2*-0.5 + xi_29 + xi_30;
                  const double xi_36 = omega_even_a*u_2_a_collide*xia_15_collide*-0.041666666666666664;
                  const double xi_40 = xi_10 + xi_30 + xi_36 + xi_39;
                  const double xi_54 = xi_10 + xi_30 + xi_36 + xi_53;
                  const double xi_60 = -xi_11;
                  const double xi_62 = xi_21 + xi_53 + xi_58 + xi_59 + xi_60 + xi_61 + xi_9;
                  const double xi_64 = u_2_a_collide*xi_41;
                  const double xi_66 = u_2_a_collide*xi_46;
                  const double xi_67 = xi_63 + xi_64 - xi_65 - xi_66;
                  const double xi_68 = -xi_63 - xi_64 + xi_65 + xi_66;
                  const double xi_69 = xi_21 + xi_39 + xi_58 + xi_59 + xi_60 + xi_61 + xi_9;
                  const double xi_72 = xi_15 + xi_22 + xi_35 + xi_71;
                  const double xi_75 = u_2_a_collide*xia_13_collide;
                  const double xi_76 = xi_75*0.25;
                  const double xi_78 = xi_45*xi_75;
                  const double xi_79 = xi_74 + xi_76 - xi_77 - xi_78;
                  const double xi_80 = -xi_74 - xi_76 + xi_77 + xi_78;
                  const double xi_81 = xi_15 + xi_22 + xi_51 + xi_71;
                  const double xi_165 = (u_2_a_collide*u_2_a_collide);
                  const double xi_166 = rho_a_collide*xi_165*-0.33333333333333331;
                  const double xi_172 = omega_even_a*(rho_a_collide*xi_169*-0.16666666666666666 - xi_164 - xi_166 - xi_168 - xi_171);
                  const double xi_176 = omega_even_a*(rho_a_collide*xi_167*-0.16666666666666666 - xi_166 - xi_173 - xi_175);
                  const double xi_179 = omega_even_a*(rho_a_collide*xi_165*-0.16666666666666666 - xi_168 - xi_173 - xi_178);
                  const double xi_184 = omega_shear_a*(rho_a_collide*xi_165 - xi_171 - xi_181 - xia_20_collide - xia_8_collide);
                  const double xi_210 = xi_179*-0.5;
                  const double xi_211 = xi_172*-0.5;
                  const double xi_216 = xi_183*0.5 + xi_210 + xi_211 + xi_215;
                  const double xi_224 = xi_176*-0.5;
                  const double xi_226 = xi_182*0.5 + xi_210 + xi_224 + xi_225;
                  const double xi_227 = u_2_a_collide*xi_199;
                  const double xi_228 = xi_219 + xi_227 + xia_17_collide - xia_9_collide;
                  const double xi_229 = xi_204*xi_228;
                  const double xi_230 = xi_207 + xi_227 + xia_25_collide - xia_3_collide;
                  const double xi_231 = xi_204*xi_230;
                  const double xi_233 = xi_184*0.5 + xi_211 + xi_224 + xi_232;
                  const double xi_244 = xi_179*0.25;
                  const double xi_246 = xi_244 + xi_245;
                  const double xi_250 = xi_244 - xi_245;
                  const double xi_253 = xi_238*(u_2_a_collide*xi_237 + xi_170 + xi_206 - xia_5_collide);
                  const double xi_256 = xi_172*0.25;
                  const double xi_257 = xi_230*xi_234;
                  const double xi_258 = xi_256 + xi_257;
                  const double xi_260 = -xi_253;
                  const double xi_266 = xi_238*(rho_a_collide*u_0_a_collide*u_2_a_collide + xi_174 + xi_218 - xia_16_collide);
                  const double xi_267 = -xi_266;
                  const double xi_270 = xi_176*0.25;
                  const double xi_271 = xi_228*xi_234;
                  const double xi_272 = xi_270 + xi_271;
                  const double xi_276 = xi_256 - xi_257;
                  const double xi_279 = xi_270 - xi_271;
                  const double forceTerm_0_a_collide = omega_shear_a*u_0_a_collide*xia_13_collide + omega_shear_a*u_1_a_collide*xia_12_collide + omega_shear_a*u_2_a_collide*xia_15_collide - xi_0*xi_3 - xi_0 - xi_1*xi_3 - xi_1 - xi_2*xi_3 - xi_2;
                  const double forceTerm_1_a_collide = xi_16 + xi_4 - xi_6 + xi_7;
                  const double forceTerm_2_a_collide = xi_16 - xi_4 + xi_6 - xi_7;
                  const double forceTerm_3_a_collide = -xi_17 - xi_18 + xi_19 + xi_25;
                  const double forceTerm_4_a_collide = xi_17 + xi_18 - xi_19 + xi_25;
                  const double forceTerm_5_a_collide = xi_26 - xi_27 + xi_28 + xi_31;
                  const double forceTerm_6_a_collide = -xi_26 + xi_27 - xi_28 + xi_31;
                  const double forceTerm_7_a_collide = -xi_35 - xi_40 - xi_50;
                  const double forceTerm_8_a_collide = -xi_40 - xi_51 - xi_52;
                  const double forceTerm_9_a_collide = -xi_35 - xi_52 - xi_54;
                  const double forceTerm_10_a_collide = -xi_50 - xi_51 - xi_54;
                  const double forceTerm_11_a_collide = xi_57 + xi_62 + xi_67;
                  const double forceTerm_12_a_collide = xi_57 + xi_68 + xi_69;
                  const double forceTerm_13_a_collide = -xi_70 - xi_72 - xi_79;
                  const double forceTerm_14_a_collide = -xi_70 - xi_80 - xi_81;
                  const double forceTerm_15_a_collide = xi_62 + xi_68 + xi_70;
                  const double forceTerm_16_a_collide = xi_67 + xi_69 + xi_70;
                  const double forceTerm_17_a_collide = -xi_57 - xi_72 - xi_80;
                  const double forceTerm_18_a_collide = -xi_57 - xi_79 - xi_81;
                  const double xib_1_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0];
                  const double xib_2_collide = _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const double xib_3_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0];
                  const double xib_4_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0];
                  const double xi_304 = -xib_4_collide;
                  const double xib_5_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0];
                  const double xi_309 = -xib_5_collide;
                  const double xib_6_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0];
                  const double xib_7_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0];
                  const double xi_287 = xib_1_collide + xib_7_collide;
                  const double xib_8_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0];
                  const double xib_9_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0];
                  const double xib_10_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0];
                  const double xib_11_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + ctr_0];
                  const double xi_99 = xib_11_collide*0.16666666666666666;
                  const double xi_100 = omega_odd_b*xi_99;
                  const double xi_114 = xib_11_collide*0.083333333333333329;
                  const double xi_116 = xi_115*xib_11_collide;
                  const double xi_117 = xi_114 - xi_116;
                  const double xi_133 = -xi_114 + xi_116;
                  const double xib_12_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3 + ctr_0];
                  const double xi_108 = xib_12_collide*0.16666666666666666;
                  const double xi_110 = omega_odd_b*xi_108;
                  const double xi_137 = xib_12_collide*0.083333333333333329;
                  const double xi_138 = xi_115*xib_12_collide;
                  const double xi_139 = xi_137 - xi_138;
                  const double xi_152 = -xi_137 + xi_138;
                  const double xib_13_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0];
                  const double xib_14_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0];
                  const double xib_15_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3 + ctr_0];
                  const double xi_86 = xib_15_collide*0.16666666666666666;
                  const double xi_89 = omega_odd_b*xi_86;
                  const double xi_119 = xib_15_collide*0.083333333333333329;
                  const double xi_120 = xi_115*xib_15_collide;
                  const double xi_121 = -xi_119 + xi_120;
                  const double xi_123 = xib_15_collide*0.25;
                  const double xi_128 = xi_127*xib_15_collide;
                  const double xi_135 = xi_119 - xi_120;
                  const double xib_16_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0];
                  const double xi_288 = xi_287 + xib_16_collide + xib_5_collide;
                  const double xi_310 = xi_309 + xib_16_collide;
                  const double xib_17_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0];
                  const double xib_18_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0];
                  const double xib_19_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0];
                  const double xib_20_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0];
                  const double xib_21_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0];
                  const double xi_305 = xi_304 + xib_21_collide;
                  const double xib_22_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0];
                  const double xib_23_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0];
                  const double xi_291 = xib_23_collide + xib_8_collide;
                  const double xib_24_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0];
                  const double xib_25_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0];
                  const double xi_292 = xi_291 + xib_20_collide + xib_25_collide;
                  const double xi_317 = -xib_25_collide;
                  const double xi_318 = xi_317 + xib_20_collide;
                  const double xib_26_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0];
                  const double xi_294 = xib_14_collide + xib_26_collide;
                  const double xi_295 = xi_294 + xib_21_collide + xib_4_collide;
                  const double xi_88 = xi_87*xib_15_collide;
                  const double xi_101 = xi_87*xib_11_collide;
                  const double xi_109 = xi_87*xib_12_collide;
                  const double rho_b_collide = xib_2_collide;
                  const double xi_284 = rho_b_collide*-0.1111111111111111;
                  const double xi_297 = rho_b_collide*-0.33333333333333331;
                  const double xi_299 = xi_288 + xi_297;
                  const double xi_302 = rho_b_collide*0.33333333333333331;
                  const double u_0_b_collide = xib_6_collide;
                  const double xi_82 = u_0_b_collide*xib_11_collide;
                  const double xi_95 = xi_82*0.083333333333333329;
                  const double xi_96 = omega_even_b*xi_95 + xi_82*-0.16666666666666666;
                  const double xi_102 = xi_82*0.33333333333333331;
                  const double xi_124 = u_0_b_collide*xi_123;
                  const double xi_129 = u_0_b_collide*xi_128;
                  const double xi_140 = -xi_95;
                  const double xi_143 = omega_even_b*xi_82*0.041666666666666664;
                  const double xi_155 = u_0_b_collide*xib_12_collide;
                  const double xi_156 = xi_155*0.25;
                  const double xi_159 = xi_127*xi_155;
                  const double xi_283 = (u_0_b_collide*u_0_b_collide);
                  const double xi_290 = rho_b_collide*xi_283*-0.33333333333333331;
                  const double xi_298 = omega_shear_b*(rho_b_collide*xi_283 - xi_292 - xi_295 - xi_297 - xib_17_collide - xib_22_collide);
                  const double xi_316 = u_0_b_collide*xi_302;
                  const double xi_319 = xi_316 + xi_318 + xib_23_collide - xib_8_collide;
                  const double xi_320 = xi_307*xi_319;
                  const double xi_321 = xi_305 + xi_316 - xib_14_collide + xib_26_collide;
                  const double xi_322 = xi_307*xi_321;
                  const double xi_332 = xi_321*xi_331;
                  const double xi_333 = -xi_332;
                  const double xi_349 = xi_319*xi_331;
                  const double xi_350 = -xi_349;
                  const double u_1_b_collide = xib_18_collide;
                  const double xi_83 = u_1_b_collide*xib_15_collide;
                  const double xi_90 = xi_83*0.33333333333333331;
                  const double xi_103 = xi_83*0.16666666666666666;
                  const double xi_104 = xi_83*0.083333333333333329;
                  const double xi_105 = omega_even_b*xi_104;
                  const double xi_106 = -xi_103 + xi_105;
                  const double xi_112 = xi_106 + xi_96;
                  const double xi_125 = u_1_b_collide*0.25;
                  const double xi_126 = xi_125*xib_11_collide;
                  const double xi_130 = u_1_b_collide*xi_127;
                  const double xi_131 = xi_130*xib_11_collide;
                  const double xi_132 = xi_124 + xi_126 - xi_129 - xi_131;
                  const double xi_134 = -xi_124 - xi_126 + xi_129 + xi_131;
                  const double xi_141 = -xi_105;
                  const double xi_145 = xi_125*xib_12_collide;
                  const double xi_147 = xi_130*xib_12_collide;
                  const double xi_153 = omega_even_b*u_1_b_collide*xib_15_collide*-0.041666666666666664;
                  const double xi_285 = (u_1_b_collide*u_1_b_collide);
                  const double xi_286 = rho_b_collide*xi_285*-0.33333333333333331 + xi_284;
                  const double xi_300 = omega_shear_b*(rho_b_collide*xi_285 - xi_295 - xi_299 - xib_24_collide - xib_9_collide);
                  const double xi_303 = u_1_b_collide*xi_302;
                  const double xi_306 = xi_303 + xi_305 + xib_14_collide - xib_26_collide;
                  const double xi_308 = xi_306*xi_307;
                  const double xi_311 = xi_303 + xi_310 + xib_1_collide - xib_7_collide;
                  const double xi_312 = xi_307*xi_311;
                  const double xi_334 = rho_b_collide*u_1_b_collide;
                  const double xi_336 = xi_335*(u_0_b_collide*xi_334 + xi_294 + xi_304 - xib_21_collide);
                  const double xi_337 = -xi_336;
                  const double xi_339 = xi_306*xi_331;
                  const double xi_342 = xi_311*xi_331;
                  const double xi_347 = -xi_342;
                  const double u_2_b_collide = xib_3_collide;
                  const double xi_84 = u_2_b_collide*xib_12_collide;
                  const double xi_91 = xi_84*0.16666666666666666;
                  const double xi_92 = xi_84*0.083333333333333329;
                  const double xi_93 = omega_even_b*xi_92;
                  const double xi_94 = -xi_91 + xi_93;
                  const double xi_97 = xi_94 + xi_96;
                  const double xi_98 = omega_even_b*xi_90 + omega_shear_b*xi_83*-0.5 + xi_90 + xi_97;
                  const double xi_107 = omega_even_b*xi_102 + omega_shear_b*xi_82*-0.5 + xi_102 + xi_106 + xi_94;
                  const double xi_111 = xi_84*0.33333333333333331;
                  const double xi_113 = omega_even_b*xi_111 + omega_shear_b*xi_84*-0.5 + xi_111 + xi_112;
                  const double xi_118 = omega_even_b*u_2_b_collide*xib_12_collide*-0.041666666666666664;
                  const double xi_122 = xi_112 + xi_118 + xi_121 + xi_92;
                  const double xi_136 = xi_112 + xi_118 + xi_135 + xi_92;
                  const double xi_142 = -xi_93;
                  const double xi_144 = xi_103 + xi_135 + xi_140 + xi_141 + xi_142 + xi_143 + xi_91;
                  const double xi_146 = u_2_b_collide*xi_123;
                  const double xi_148 = u_2_b_collide*xi_128;
                  const double xi_149 = xi_145 + xi_146 - xi_147 - xi_148;
                  const double xi_150 = -xi_145 - xi_146 + xi_147 + xi_148;
                  const double xi_151 = xi_103 + xi_121 + xi_140 + xi_141 + xi_142 + xi_143 + xi_91;
                  const double xi_154 = xi_104 + xi_117 + xi_153 + xi_97;
                  const double xi_157 = u_2_b_collide*xib_11_collide;
                  const double xi_158 = xi_157*0.25;
                  const double xi_160 = xi_127*xi_157;
                  const double xi_161 = xi_156 + xi_158 - xi_159 - xi_160;
                  const double xi_162 = -xi_156 - xi_158 + xi_159 + xi_160;
                  const double xi_163 = xi_104 + xi_133 + xi_153 + xi_97;
                  const double xi_281 = (u_2_b_collide*u_2_b_collide);
                  const double xi_282 = rho_b_collide*xi_281*-0.33333333333333331;
                  const double xi_289 = omega_even_b*(rho_b_collide*xi_283*-0.16666666666666666 - xi_282 - xi_286 - xi_288);
                  const double xi_293 = omega_even_b*(rho_b_collide*xi_285*-0.16666666666666666 - xi_282 - xi_284 - xi_290 - xi_292);
                  const double xi_296 = omega_even_b*(rho_b_collide*xi_281*-0.16666666666666666 - xi_286 - xi_290 - xi_295);
                  const double xi_301 = omega_shear_b*(rho_b_collide*xi_281 - xi_292 - xi_299 - xib_13_collide - xib_19_collide);
                  const double xi_313 = xi_296*-0.5;
                  const double xi_314 = xi_289*-0.5;
                  const double xi_315 = xi_215 + xi_300*0.5 + xi_313 + xi_314;
                  const double xi_323 = xi_293*-0.5;
                  const double xi_324 = xi_225 + xi_298*0.5 + xi_313 + xi_323;
                  const double xi_325 = u_2_b_collide*xi_302;
                  const double xi_326 = xi_310 + xi_325 - xib_1_collide + xib_7_collide;
                  const double xi_327 = xi_307*xi_326;
                  const double xi_328 = xi_318 + xi_325 - xib_23_collide + xib_8_collide;
                  const double xi_329 = xi_307*xi_328;
                  const double xi_330 = xi_232 + xi_301*0.5 + xi_314 + xi_323;
                  const double xi_338 = xi_296*0.25;
                  const double xi_340 = xi_338 + xi_339;
                  const double xi_341 = xi_338 - xi_339;
                  const double xi_343 = xi_335*(u_2_b_collide*xi_334 + xi_287 + xi_309 - xib_16_collide);
                  const double xi_344 = xi_289*0.25;
                  const double xi_345 = xi_326*xi_331;
                  const double xi_346 = xi_344 + xi_345;
                  const double xi_348 = -xi_343;
                  const double xi_351 = xi_335*(rho_b_collide*u_0_b_collide*u_2_b_collide + xi_291 + xi_317 - xib_20_collide);
                  const double xi_352 = -xi_351;
                  const double xi_353 = xi_293*0.25;
                  const double xi_354 = xi_328*xi_331;
                  const double xi_355 = xi_353 + xi_354;
                  const double xi_356 = xi_344 - xi_345;
                  const double xi_357 = xi_353 - xi_354;
                  const double forceTerm_0_b_collide = omega_shear_b*u_0_b_collide*xib_11_collide + omega_shear_b*u_1_b_collide*xib_15_collide + omega_shear_b*u_2_b_collide*xib_12_collide - xi_82*xi_85 - xi_82 - xi_83*xi_85 - xi_83 - xi_84*xi_85 - xi_84;
                  const double forceTerm_1_b_collide = xi_86 - xi_88 + xi_89 + xi_98;
                  const double forceTerm_2_b_collide = -xi_86 + xi_88 - xi_89 + xi_98;
                  const double forceTerm_3_b_collide = -xi_100 + xi_101 + xi_107 - xi_99;
                  const double forceTerm_4_b_collide = xi_100 - xi_101 + xi_107 + xi_99;
                  const double forceTerm_5_b_collide = xi_108 - xi_109 + xi_110 + xi_113;
                  const double forceTerm_6_b_collide = -xi_108 + xi_109 - xi_110 + xi_113;
                  const double forceTerm_7_b_collide = -xi_117 - xi_122 - xi_132;
                  const double forceTerm_8_b_collide = -xi_122 - xi_133 - xi_134;
                  const double forceTerm_9_b_collide = -xi_117 - xi_134 - xi_136;
                  const double forceTerm_10_b_collide = -xi_132 - xi_133 - xi_136;
                  const double forceTerm_11_b_collide = xi_139 + xi_144 + xi_149;
                  const double forceTerm_12_b_collide = xi_139 + xi_150 + xi_151;
                  const double forceTerm_13_b_collide = -xi_152 - xi_154 - xi_161;
                  const double forceTerm_14_b_collide = -xi_152 - xi_162 - xi_163;
                  const double forceTerm_15_b_collide = xi_144 + xi_150 + xi_152;
                  const double forceTerm_16_b_collide = xi_149 + xi_151 + xi_152;
                  const double forceTerm_17_b_collide = -xi_139 - xi_154 - xi_162;
                  const double forceTerm_18_b_collide = -xi_139 - xi_161 - xi_163;
                  const double tmp_a0 = forceTerm_0_a_collide + xi_172 + xi_176 + xi_179 - xi_182 - xi_183 - xi_184 + xi_198 + xia_24_collide;
                  const double tmp_a1 = forceTerm_1_a_collide - xi_205 - xi_209 + xi_216 + xia_23_collide;
                  const double tmp_a2 = forceTerm_2_a_collide + xi_205 + xi_209 + xi_216 + xia_18_collide;
                  const double tmp_a3 = forceTerm_3_a_collide + xi_221 + xi_223 + xi_226 + xia_19_collide;
                  const double tmp_a4 = forceTerm_4_a_collide - xi_221 - xi_223 + xi_226 + xia_22_collide;
                  const double tmp_a5 = forceTerm_5_a_collide - xi_229 - xi_231 + xi_233 + xia_8_collide;
                  const double tmp_a6 = forceTerm_6_a_collide + xi_229 + xi_231 + xi_233 + xia_20_collide;
                  const double tmp_a7 = forceTerm_7_a_collide + xi_236 + xi_240 + xi_243 + xi_246 + xia_6_collide;
                  const double tmp_a8 = forceTerm_8_a_collide + xi_235 + xi_239 + xi_246 + xi_248 + xia_11_collide;
                  const double tmp_a9 = forceTerm_9_a_collide + xi_236 + xi_239 + xi_249 + xi_250 + xia_1_collide;
                  const double tmp_a10 = forceTerm_10_a_collide + xi_235 + xi_240 + xi_250 + xi_251 + xia_2_collide;
                  const double tmp_a11 = forceTerm_11_a_collide + xi_252 + xi_253 + xi_255 + xi_258 + xia_21_collide;
                  const double tmp_a12 = forceTerm_12_a_collide + xi_258 + xi_259 + xi_260 + xi_263 + xia_3_collide;
                  const double tmp_a13 = forceTerm_13_a_collide + xi_265 + xi_267 + xi_269 + xi_272 + xia_9_collide;
                  const double tmp_a14 = forceTerm_14_a_collide + xi_264 + xi_266 + xi_272 + xi_274 + xia_10_collide;
                  const double tmp_a15 = forceTerm_15_a_collide + xi_252 + xi_260 + xi_275 + xi_276 + xia_25_collide;
                  const double tmp_a16 = forceTerm_16_a_collide + xi_253 + xi_259 + xi_276 + xi_277 + xia_5_collide;
                  const double tmp_a17 = forceTerm_17_a_collide + xi_265 + xi_266 + xi_278 + xi_279 + xia_16_collide;
                  const double tmp_a18 = forceTerm_18_a_collide + xi_264 + xi_267 + xi_279 + xi_280 + xia_17_collide;
                  const double tmp_b0 = forceTerm_0_b_collide + xi_198 + xi_289 + xi_293 + xi_296 - xi_298 - xi_300 - xi_301 + xib_10_collide;
                  const double tmp_b1 = forceTerm_1_b_collide - xi_308 - xi_312 + xi_315 + xib_9_collide;
                  const double tmp_b2 = forceTerm_2_b_collide + xi_308 + xi_312 + xi_315 + xib_24_collide;
                  const double tmp_b3 = forceTerm_3_b_collide + xi_320 + xi_322 + xi_324 + xib_17_collide;
                  const double tmp_b4 = forceTerm_4_b_collide - xi_320 - xi_322 + xi_324 + xib_22_collide;
                  const double tmp_b5 = forceTerm_5_b_collide - xi_327 - xi_329 + xi_330 + xib_19_collide;
                  const double tmp_b6 = forceTerm_6_b_collide + xi_327 + xi_329 + xi_330 + xib_13_collide;
                  const double tmp_b7 = forceTerm_7_b_collide + xi_243 + xi_333 + xi_337 + xi_340 + xib_26_collide;
                  const double tmp_b8 = forceTerm_8_b_collide + xi_248 + xi_332 + xi_336 + xi_340 + xib_4_collide;
                  const double tmp_b9 = forceTerm_9_b_collide + xi_249 + xi_333 + xi_336 + xi_341 + xib_21_collide;
                  const double tmp_b10 = forceTerm_10_b_collide + xi_251 + xi_332 + xi_337 + xi_341 + xib_14_collide;
                  const double tmp_b11 = forceTerm_11_b_collide + xi_255 + xi_342 + xi_343 + xi_346 + xib_5_collide;
                  const double tmp_b12 = forceTerm_12_b_collide + xi_263 + xi_346 + xi_347 + xi_348 + xib_1_collide;
                  const double tmp_b13 = forceTerm_13_b_collide + xi_269 + xi_350 + xi_352 + xi_355 + xib_23_collide;
                  const double tmp_b14 = forceTerm_14_b_collide + xi_274 + xi_349 + xi_351 + xi_355 + xib_25_collide;
                  const double tmp_b15 = forceTerm_15_b_collide + xi_275 + xi_342 + xi_348 + xi_356 + xib_7_collide;
                  const double tmp_b16 = forceTerm_16_b_collide + xi_277 + xi_343 + xi_347 + xi_356 + xib_16_collide;
                  const double tmp_b17 = forceTerm_17_b_collide + xi_278 + xi_350 + xi_351 + xi_357 + xib_20_collide;
                  const double tmp_b18 = forceTerm_18_b_collide + xi_280 + xi_349 + xi_352 + xi_357 + xib_8_collide;
                  const double xirecolor_0 = tmp_a0 + tmp_b0;
                  const double xirecolor_1 = _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const double xirecolor_2 = ((1.0) / (xirecolor_1));
                  const double xi_364 = xirecolor_2*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const double xirecolor_3 = xirecolor_2*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0];
                  const double xirecolor_4 = tmp_a1 + tmp_b1;
                  const double xirecolor_5 = xi_189;
                  const double xirecolor_6 = ((1.0) / (xirecolor_5));
                  const double xirecolor_7 = xirecolor_6*_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0];
                  const bool xirecolor_8 = xirecolor_5 > 0.0;
                  const double xirecolor_9 = beta*((1.0) / ((xirecolor_1*xirecolor_1)))*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const double xirecolor_10 = xirecolor_9*(0.055555555555555552*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + 0.055555555555555552*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const double xirecolor_11 = xirecolor_10*((xirecolor_8) ? (xirecolor_7): (0.0));
                  const double xirecolor_12 = tmp_a2 + tmp_b2;
                  const double xirecolor_13 = xirecolor_10*((xirecolor_8) ? (-xirecolor_7): (0.0));
                  const double xirecolor_14 = tmp_a3 + tmp_b3;
                  const double xirecolor_15 = xirecolor_6*_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0];
                  const double xirecolor_16 = xirecolor_10*((xirecolor_8) ? (-xirecolor_15): (0.0));
                  const double xirecolor_17 = tmp_a4 + tmp_b4;
                  const double xirecolor_18 = xirecolor_10*((xirecolor_8) ? (xirecolor_15): (0.0));
                  const double xirecolor_19 = tmp_a5 + tmp_b5;
                  const double xirecolor_20 = xirecolor_6*_data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + 2*_stride_color_gradient_3 + ctr_0];
                  const double xirecolor_21 = xirecolor_10*((xirecolor_8) ? (xirecolor_20): (0.0));
                  const double xirecolor_22 = tmp_a6 + tmp_b6;
                  const double xirecolor_23 = xirecolor_10*((xirecolor_8) ? (-xirecolor_20): (0.0));
                  const double xirecolor_24 = tmp_a7 + tmp_b7;
                  const double xirecolor_25 = xi_241;
                  const double xirecolor_26 = xirecolor_6*0.70710678118654757;
                  const double xi_358 = xirecolor_25*xirecolor_26;
                  const double xirecolor_27 = xirecolor_9*(0.027777777777777776*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + 0.027777777777777776*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const double xirecolor_28 = xirecolor_27*((xirecolor_8) ? (-xi_358): (0.0));
                  const double xirecolor_29 = tmp_a8 + tmp_b8;
                  const double xirecolor_30 = xi_247;
                  const double xi_359 = xirecolor_26*xirecolor_30;
                  const double xirecolor_31 = xirecolor_27*((xirecolor_8) ? (xi_359): (0.0));
                  const double xirecolor_32 = tmp_a9 + tmp_b9;
                  const double xirecolor_33 = xirecolor_27*((xirecolor_8) ? (-xi_359): (0.0));
                  const double xirecolor_34 = tmp_a10 + tmp_b10;
                  const double xirecolor_35 = xirecolor_27*((xirecolor_8) ? (xi_358): (0.0));
                  const double xirecolor_36 = tmp_a11 + tmp_b11;
                  const double xirecolor_37 = xi_254;
                  const double xi_360 = xirecolor_26*xirecolor_37;
                  const double xirecolor_38 = xirecolor_27*((xirecolor_8) ? (xi_360): (0.0));
                  const double xirecolor_39 = tmp_a12 + tmp_b12;
                  const double xirecolor_40 = xi_261;
                  const double xirecolor_41 = xirecolor_40 + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + _stride_color_gradient_3 + ctr_0];
                  const double xi_361 = xirecolor_26*xirecolor_41;
                  const double xirecolor_42 = xirecolor_27*((xirecolor_8) ? (-xi_361): (0.0));
                  const double xirecolor_43 = tmp_a13 + tmp_b13;
                  const double xirecolor_44 = xirecolor_40 + _data_color_gradient[_stride_color_gradient_1*ctr_1 + _stride_color_gradient_2*ctr_2 + ctr_0];
                  const double xi_362 = xirecolor_26*xirecolor_44;
                  const double xirecolor_45 = xirecolor_27*((xirecolor_8) ? (-xi_362): (0.0));
                  const double xirecolor_46 = tmp_a14 + tmp_b14;
                  const double xirecolor_47 = xi_273;
                  const double xi_363 = xirecolor_26*xirecolor_47;
                  const double xirecolor_48 = xirecolor_27*((xirecolor_8) ? (xi_363): (0.0));
                  const double xirecolor_49 = tmp_a15 + tmp_b15;
                  const double xirecolor_50 = xirecolor_27*((xirecolor_8) ? (xi_361): (0.0));
                  const double xirecolor_51 = tmp_a16 + tmp_b16;
                  const double xirecolor_52 = xirecolor_27*((xirecolor_8) ? (-xi_360): (0.0));
                  const double xirecolor_53 = tmp_a17 + tmp_b17;
                  const double xirecolor_54 = xirecolor_27*((xirecolor_8) ? (-xi_363): (0.0));
                  const double xirecolor_55 = tmp_a18 + tmp_b18;
                  const double xirecolor_56 = xirecolor_27*((xirecolor_8) ? (xi_362): (0.0));
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


void ColorGradientCollideSweepDoublePrecisionAVX::run(IBlock * block)
{
   
    auto phasefield = block->getData< field::GhostLayerField<double, 1> >(phasefieldID);
    auto rho_b = block->getData< field::GhostLayerField<double, 1> >(rho_bID);
    auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
    auto rho_a = block->getData< field::GhostLayerField<double, 1> >(rho_aID);
    auto pdfs_b = block->getData< field::GhostLayerField<double, 19> >(pdfs_bID);
    auto force_a = block->getData< field::GhostLayerField<double, 3> >(force_aID);
    auto pdfs_a = block->getData< field::GhostLayerField<double, 19> >(pdfs_aID);
    auto force_b = block->getData< field::GhostLayerField<double, 3> >(force_bID);
    auto color_gradient = block->getData< field::GhostLayerField<double, 3> >(color_gradientID);

    auto & beta = this->beta_;
    auto & omega_odd_b = this->omega_odd_b_;
    auto & sigma = this->sigma_;
    auto & omega_even_b = this->omega_even_b_;
    auto & omega_even_a = this->omega_even_a_;
    auto & omega_shear_b = this->omega_shear_b_;
    auto & omega_shear_a = this->omega_shear_a_;
    auto & omega_odd_a = this->omega_odd_a_;
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(color_gradient->nrOfGhostLayers()))
    double * RESTRICT const _data_color_gradient = color_gradient->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(color_gradient->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) color_gradient->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force_a->nrOfGhostLayers()))
    double * RESTRICT const _data_force_a = force_a->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT const _data_phasefield = phasefield->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) phasefield->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_a = rho_a->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_b = rho_b->dataAt(0, 0, 0, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
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
    internal_6ed8e43fbb9656b49ba1d2a7ec0e653c::colorgradientcollidesweepdoubleprecisionavx_colorgradientcollidesweepdoubleprecisionavx(_data_color_gradient, _data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_color_gradient_0, _size_color_gradient_1, _size_color_gradient_2, _stride_color_gradient_1, _stride_color_gradient_2, _stride_color_gradient_3, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_even_a, omega_even_b, omega_odd_a, omega_odd_b, omega_shear_a, omega_shear_b, sigma);
    
}


void ColorGradientCollideSweepDoublePrecisionAVX::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block)
{
   
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    auto phasefield = block->getData< field::GhostLayerField<double, 1> >(phasefieldID);
    auto rho_b = block->getData< field::GhostLayerField<double, 1> >(rho_bID);
    auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
    auto rho_a = block->getData< field::GhostLayerField<double, 1> >(rho_aID);
    auto pdfs_b = block->getData< field::GhostLayerField<double, 19> >(pdfs_bID);
    auto force_a = block->getData< field::GhostLayerField<double, 3> >(force_aID);
    auto pdfs_a = block->getData< field::GhostLayerField<double, 19> >(pdfs_aID);
    auto force_b = block->getData< field::GhostLayerField<double, 3> >(force_bID);
    auto color_gradient = block->getData< field::GhostLayerField<double, 3> >(color_gradientID);

    auto & beta = this->beta_;
    auto & omega_odd_b = this->omega_odd_b_;
    auto & sigma = this->sigma_;
    auto & omega_even_b = this->omega_even_b_;
    auto & omega_even_a = this->omega_even_a_;
    auto & omega_shear_b = this->omega_shear_b_;
    auto & omega_shear_a = this->omega_shear_a_;
    auto & omega_odd_a = this->omega_odd_a_;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(color_gradient->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(color_gradient->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(color_gradient->nrOfGhostLayers()))
    double * RESTRICT const _data_color_gradient = color_gradient->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(color_gradient->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) color_gradient->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force_a->nrOfGhostLayers()))
    double * RESTRICT const _data_force_a = force_a->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT const _data_phasefield = phasefield->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) phasefield->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_a = rho_a->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_b = rho_b->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
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
    internal_6ed8e43fbb9656b49ba1d2a7ec0e653c::colorgradientcollidesweepdoubleprecisionavx_colorgradientcollidesweepdoubleprecisionavx(_data_color_gradient, _data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_color_gradient_0, _size_color_gradient_1, _size_color_gradient_2, _stride_color_gradient_1, _stride_color_gradient_2, _stride_color_gradient_3, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_even_a, omega_even_b, omega_odd_a, omega_odd_b, omega_shear_a, omega_shear_b, sigma);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
