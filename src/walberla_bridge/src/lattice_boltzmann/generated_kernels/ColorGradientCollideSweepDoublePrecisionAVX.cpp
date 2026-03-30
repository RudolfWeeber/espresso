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

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 17fc54c872bd8ceabf271a7e9e636c7c583f55af


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
static FUNC_PREFIX void colorgradientcollidesweepdoubleprecisionavx_colorgradientcollidesweepdoubleprecisionavx(double * RESTRICT const _data_force_a, double * RESTRICT const _data_force_b, double * RESTRICT  _data_pdfs_a, double * RESTRICT  _data_pdfs_b, double * RESTRICT const _data_phasefield, double * RESTRICT const _data_rho_a, double * RESTRICT const _data_rho_b, double * RESTRICT const _data_velocity, int64_t const _size_force_a_0, int64_t const _size_force_a_1, int64_t const _size_force_a_2, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3, double beta, double omega_even_a, double omega_even_b, double omega_odd_a, double omega_odd_b, double omega_shear_a, double omega_shear_b, double sigma)
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
      const double xi_237 = omega_shear_a*omega_shear_b*((1.0) / (omega_shear_a + omega_shear_b));
      const double xi_238 = xi_237*2.0;
      const double xi_239 = xi_237*8.0;
      const double xi_240 = omega_shear_a*-4.0 + xi_239;
      const double xi_242 = omega_shear_b*-4.0 + xi_239;
      const double xi_251 = omega_odd_a*0.5;
      const double xi_283 = omega_odd_a*0.25;
      const double xi_287 = omega_shear_a*0.25;
      const double xi_372 = omega_odd_b*0.5;
      const double xi_396 = omega_odd_b*0.25;
      const double xi_400 = omega_shear_b*0.25;
      const double rr_0_a_collide = 0.0;
      const double xi_5 = rr_0_a_collide*0.25;
      const double rr_0_b_collide = 0.0;
      const double xi_87 = rr_0_b_collide*0.25;
#ifdef _OPENMP
      #pragma omp for schedule(static)
#endif
      for (int64_t ctr_2 = 1; ctr_2 < _size_force_a_2 - 1; ctr_2 += 1)
      {
         for (int64_t ctr_1 = 1; ctr_1 < _size_force_a_1 - 1; ctr_1 += 1)
         {
            {
               for (int64_t ctr_0 = 1; ctr_0 < (int64_t)((_size_force_a_0 - 2) / (4)) * (4) + 1; ctr_0 += 4)
               {
                  const __m256d xi_185 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_186 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_187 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_188 = _mm256_mul_pd(_mm256_set_pd(-0.013888888888888888,-0.013888888888888888,-0.013888888888888888,-0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_189 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_190 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_191 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_187,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_188),xi_189),xi_190);
                  const __m256d xi_192 = _mm256_mul_pd(_mm256_set_pd(-0.055555555555555552,-0.055555555555555552,-0.055555555555555552,-0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_193 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_194 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_195 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_196 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_193,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_192),xi_194),xi_195);
                  const __m256d xi_197 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_198 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_199 = _mm256_add_pd(_mm256_mul_pd(xi_197,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_198);
                  const __m256d xi_200 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0]));
                  const __m256d xi_201 = _mm256_mul_pd(xi_200,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_202 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0]));
                  const __m256d xi_203 = _mm256_add_pd(xi_201,xi_202);
                  const __m256d xi_204 = _mm256_mul_pd(_mm256_set_pd(-0.22222222222222221,-0.22222222222222221,-0.22222222222222221,-0.22222222222222221),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0]));
                  const __m256d xi_205 = _mm256_mul_pd(_mm256_set_pd(0.22222222222222221,0.22222222222222221,0.22222222222222221,0.22222222222222221),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0]));
                  const __m256d xi_206 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0]));
                  const __m256d xi_207 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0]));
                  const __m256d xi_208 = _mm256_add_pd(_mm256_mul_pd(xi_206,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_207);
                  const __m256d xi_209 = _mm256_add_pd(_mm256_add_pd(xi_204,xi_205),xi_208);
                  const __m256d xi_210 = _mm256_add_pd(xi_203,xi_209);
                  const __m256d xi_211 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_185,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_186),xi_191),xi_196),xi_199),xi_210);
                  const __m256d xi_212 = _mm256_mul_pd(xi_211,xi_211);
                  const __m256d xi_213 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1]));
                  const __m256d xi_214 = _mm256_mul_pd(xi_213,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_215 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1]));
                  const __m256d xi_216 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_201,xi_202),xi_214),xi_215);
                  const __m256d xi_217 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1]));
                  const __m256d xi_218 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1]));
                  const __m256d xi_219 = _mm256_add_pd(_mm256_mul_pd(xi_217,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_218);
                  const __m256d xi_220 = _mm256_add_pd(_mm256_mul_pd(xi_186,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_185);
                  const __m256d xi_221 = _mm256_mul_pd(_mm256_set_pd(0.22222222222222221,0.22222222222222221,0.22222222222222221,0.22222222222222221),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0]));
                  const __m256d xi_222 = _mm256_mul_pd(xi_221,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_223 = _mm256_mul_pd(_mm256_set_pd(0.22222222222222221,0.22222222222222221,0.22222222222222221,0.22222222222222221),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0]));
                  const __m256d xi_224 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_207,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_206),xi_222),xi_223);
                  const __m256d xi_225 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_195,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_191),xi_193),xi_216),xi_219),xi_220),xi_224);
                  const __m256d xi_226 = _mm256_mul_pd(xi_225,xi_225);
                  const __m256d xi_227 = _mm256_add_pd(_mm256_mul_pd(xi_218,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_217);
                  const __m256d xi_228 = _mm256_add_pd(xi_214,xi_215);
                  const __m256d xi_229 = _mm256_mul_pd(_mm256_set_pd(0.22222222222222221,0.22222222222222221,0.22222222222222221,0.22222222222222221),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1]));
                  const __m256d xi_230 = _mm256_mul_pd(xi_229,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_231 = _mm256_mul_pd(_mm256_set_pd(0.22222222222222221,0.22222222222222221,0.22222222222222221,0.22222222222222221),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1]));
                  const __m256d xi_232 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_198,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_197),xi_230),xi_231);
                  const __m256d xi_233 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_190,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_187),xi_188),xi_189),xi_196),xi_220),xi_227),xi_228),xi_232);
                  const __m256d xi_234 = _mm256_mul_pd(xi_233,xi_233);
                  const __m256d xi_235 = _mm256_add_pd(_mm256_add_pd(xi_212,xi_226),xi_234);
                  const __m256d xi_236 = _mm256_sqrt_pd(xi_235);
                  const __m256d xi_241 = _mm256_mul_pd(_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]));
                  const __m256d xi_243 = _mm256_mul_pd(_mm256_mul_pd(xi_236,_mm256_set_pd(sigma,sigma,sigma,sigma)),_mm256_blendv_pd(_mm256_blendv_pd(_mm256_blendv_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_241,_mm256_set_pd(xi_242,xi_242,xi_242,xi_242)),_mm256_mul_pd(_mm256_set_pd(xi_242,xi_242,xi_242,xi_242),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))),_mm256_set_pd(xi_238,xi_238,xi_238,xi_238)),_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_241,_mm256_set_pd(xi_240,xi_240,xi_240,xi_240)),_mm256_mul_pd(_mm256_mul_pd(_mm256_set_pd(-1.0,-1.0,-1.0,-1.0),_mm256_set_pd(xi_240,xi_240,xi_240,xi_240)),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))),_mm256_set_pd(xi_238,xi_238,xi_238,xi_238)),_mm256_cmp_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NGE_UQ)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),_mm256_cmp_pd(_mm256_set_pd(-0.5,-0.5,-0.5,-0.5),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NLE_UQ)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),_mm256_cmp_pd(_mm256_set_pd(0.5,0.5,0.5,0.5),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NGE_UQ)));
                  const __m256d xi_244 = _mm256_cmp_pd(xi_236,_mm256_set_pd(0.0,0.0,0.0,0.0),_CMP_NLE_UQ);
                  const __m256d xi_245 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_243,_mm256_set_pd(0.25,0.25,0.25,0.25)),xi_244);
                  const __m256d xi_257 = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),xi_235);
                  const __m256d xi_258 = _mm256_mul_pd(xi_243,_mm256_set_pd(1.125,1.125,1.125,1.125));
                  const __m256d xi_259 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_226,xi_257),_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552)),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_244);
                  const __m256d xi_263 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552)),_mm256_mul_pd(xi_225,xi_225)),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_244);
                  const __m256d xi_271 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552)),_mm256_mul_pd(xi_233,xi_233)),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_244);
                  const __m256d xi_274 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_234,xi_257),_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552)),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_244);
                  const __m256d xi_280 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552)),_mm256_mul_pd(xi_211,xi_211)),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_244);
                  const __m256d xi_282 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_212,xi_257),_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552)),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_244);
                  const __m256d xi_290 = _mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_291 = _mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_292 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.1111111111111111,0.1111111111111111,0.1111111111111111,0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1])));
                  const __m256d xi_293 = _mm256_add_pd(xi_192,xi_194);
                  const __m256d xi_294 = _mm256_add_pd(xi_232,xi_293);
                  const __m256d xi_295 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(-0.027777777777777776,-0.027777777777777776,-0.027777777777777776,-0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1])));
                  const __m256d xi_296 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_223,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_221),xi_295);
                  const __m256d xi_297 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_202,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_290,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_200),xi_208),xi_291),xi_292),xi_294),xi_296);
                  const __m256d xi_298 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_297,xi_297)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_244);
                  const __m256d xi_302 = _mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_303 = _mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_304 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(-0.027777777777777776,-0.027777777777777776,-0.027777777777777776,-0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1])));
                  const __m256d xi_305 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.1111111111111111,0.1111111111111111,0.1111111111111111,0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1])));
                  const __m256d xi_306 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_302,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_203),xi_224),xi_294),xi_303),xi_304),xi_305);
                  const __m256d xi_307 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_306,xi_306)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_244);
                  const __m256d xi_308 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_306,xi_306)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_244);
                  const __m256d xi_310 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_297,xi_297)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_244);
                  const __m256d xi_313 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.1111111111111111,0.1111111111111111,0.1111111111111111,0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0])),_mm256_mul_pd(_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0])));
                  const __m256d xi_314 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_199,xi_204),xi_205),xi_293);
                  const __m256d xi_315 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_215,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_303,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_213),xi_302);
                  const __m256d xi_316 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_227,xi_296),xi_313),xi_314),xi_315);
                  const __m256d xi_317 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_316,xi_316)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_244);
                  const __m256d xi_323 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.1111111111111111,0.1111111111111111,0.1111111111111111,0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0])),_mm256_mul_pd(_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0])));
                  const __m256d xi_324 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_291,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_219),xi_290);
                  const __m256d xi_325 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_222,xi_223),xi_228),xi_304),xi_314),xi_323),xi_324);
                  const __m256d xi_326 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_325,xi_325)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_244);
                  const __m256d xi_331 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.1111111111111111,0.1111111111111111,0.1111111111111111,0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1])));
                  const __m256d xi_332 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_209,xi_216),xi_227),xi_230),xi_231),xi_295),xi_304),xi_331);
                  const __m256d xi_333 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_332,xi_332)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_244);
                  const __m256d xi_337 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.1111111111111111,0.1111111111111111,0.1111111111111111,0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1])),_mm256_mul_pd(_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1])));
                  const __m256d xi_338 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_231,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_210),xi_229),xi_315),xi_324),xi_337);
                  const __m256d xi_339 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_338,xi_338)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_244);
                  const __m256d xi_340 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_325,xi_325)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_244);
                  const __m256d xi_342 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_316,xi_316)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_244);
                  const __m256d xi_343 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_338,xi_338)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_244);
                  const __m256d xi_345 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_258,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_257,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_332,xi_332)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_244);
                  const __m256d xia_1_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_253 = _mm256_mul_pd(xia_1_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xia_2_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_3_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_4_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_5_collide = _mm256_loadu_pd(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3 + ctr_0]);
                  const __m256d xi_4 = _mm256_mul_pd(xia_5_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_7 = _mm256_mul_pd(xi_4,_mm256_set_pd(omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a));
                  const __m256d xi_37 = _mm256_mul_pd(xia_5_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_38 = _mm256_mul_pd(xia_5_collide,_mm256_set_pd(xi_33,xi_33,xi_33,xi_33));
                  const __m256d xi_39 = _mm256_add_pd(_mm256_mul_pd(xi_37,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_38);
                  const __m256d xi_41 = _mm256_mul_pd(xia_5_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_46 = _mm256_mul_pd(xia_5_collide,_mm256_set_pd(xi_45,xi_45,xi_45,xi_45));
                  const __m256d xi_53 = _mm256_add_pd(_mm256_mul_pd(xi_38,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_37);
                  const __m256d xia_6_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_7_collide = _mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]);
                  const __m256d xia_8_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_9_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_10_collide = _mm256_loadu_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0]);
                  const __m256d xia_11_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_174 = _mm256_add_pd(xia_11_collide,xia_4_collide);
                  const __m256d xia_12_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_13_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_175 = _mm256_add_pd(_mm256_add_pd(xi_174,xia_13_collide),xia_8_collide);
                  const __m256d xi_267 = _mm256_mul_pd(xia_13_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_268 = _mm256_add_pd(xi_267,xia_8_collide);
                  const __m256d xia_14_collide = _mm256_loadu_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0]);
                  const __m256d xia_15_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_16_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_17_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_177 = _mm256_add_pd(xia_17_collide,xia_3_collide);
                  const __m256d xia_18_collide = _mm256_load_pd(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + ctr_0]);
                  const __m256d xi_17 = _mm256_mul_pd(xia_18_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_18 = _mm256_mul_pd(xi_17,_mm256_set_pd(omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a));
                  const __m256d xi_32 = _mm256_mul_pd(xia_18_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_34 = _mm256_mul_pd(xia_18_collide,_mm256_set_pd(xi_33,xi_33,xi_33,xi_33));
                  const __m256d xi_35 = _mm256_add_pd(_mm256_mul_pd(xi_34,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_32);
                  const __m256d xi_51 = _mm256_add_pd(_mm256_mul_pd(xi_32,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_34);
                  const __m256d xia_19_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_20_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_170 = _mm256_add_pd(xia_15_collide,xia_20_collide);
                  const __m256d xia_21_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_178 = _mm256_add_pd(_mm256_add_pd(xi_177,xia_21_collide),xia_9_collide);
                  const __m256d xi_248 = _mm256_mul_pd(xia_21_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_249 = _mm256_add_pd(xi_248,xia_9_collide);
                  const __m256d xia_22_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_23_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_171 = _mm256_add_pd(_mm256_add_pd(xi_170,xia_1_collide),xia_23_collide);
                  const __m256d xi_254 = _mm256_add_pd(xi_253,xia_23_collide);
                  const __m256d xia_24_collide = _mm256_load_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0]);
                  const __m256d xia_25_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0]);
                  const __m256d xia_26_collide = _mm256_loadu_pd(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3 + ctr_0]);
                  const __m256d xi_26 = _mm256_mul_pd(xia_26_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_28 = _mm256_mul_pd(xi_26,_mm256_set_pd(omega_odd_a,omega_odd_a,omega_odd_a,omega_odd_a));
                  const __m256d xi_55 = _mm256_mul_pd(xia_26_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_56 = _mm256_mul_pd(xia_26_collide,_mm256_set_pd(xi_33,xi_33,xi_33,xi_33));
                  const __m256d xi_57 = _mm256_add_pd(_mm256_mul_pd(xi_56,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_55);
                  const __m256d xi_70 = _mm256_add_pd(_mm256_mul_pd(xi_55,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_56);
                  const __m256d xi_6 = _mm256_mul_pd(xia_5_collide,_mm256_set_pd(xi_5,xi_5,xi_5,xi_5));
                  const __m256d xi_19 = _mm256_mul_pd(xia_18_collide,_mm256_set_pd(xi_5,xi_5,xi_5,xi_5));
                  const __m256d xi_27 = _mm256_mul_pd(xia_26_collide,_mm256_set_pd(xi_5,xi_5,xi_5,xi_5));
                  const __m256d rho_a_collide = xia_7_collide;
                  const __m256d xi_167 = _mm256_mul_pd(rho_a_collide,_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111));
                  const __m256d xi_180 = _mm256_mul_pd(rho_a_collide,_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_181 = _mm256_add_pd(xi_175,xi_180);
                  const __m256d xi_246 = _mm256_mul_pd(rho_a_collide,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d u_0_a_collide = xia_24_collide;
                  const __m256d xi_0 = _mm256_mul_pd(u_0_a_collide,xia_18_collide);
                  const __m256d xi_13 = _mm256_mul_pd(xi_0,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_14 = _mm256_add_pd(_mm256_mul_pd(xi_0,_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666)),_mm256_mul_pd(xi_13,_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a)));
                  const __m256d xi_20 = _mm256_mul_pd(xi_0,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d xi_42 = _mm256_mul_pd(u_0_a_collide,xi_41);
                  const __m256d xi_47 = _mm256_mul_pd(u_0_a_collide,xi_46);
                  const __m256d xi_58 = _mm256_mul_pd(xi_13,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_61 = _mm256_mul_pd(_mm256_mul_pd(xi_0,_mm256_set_pd(0.041666666666666664,0.041666666666666664,0.041666666666666664,0.041666666666666664)),_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_73 = _mm256_mul_pd(u_0_a_collide,xia_26_collide);
                  const __m256d xi_74 = _mm256_mul_pd(xi_73,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_77 = _mm256_mul_pd(xi_73,_mm256_set_pd(xi_45,xi_45,xi_45,xi_45));
                  const __m256d xi_166 = _mm256_mul_pd(u_0_a_collide,u_0_a_collide);
                  const __m256d xi_173 = _mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_166),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_182 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_178,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_181,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_12_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_16_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_a_collide,xi_166)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a));
                  const __m256d xi_264 = _mm256_mul_pd(u_0_a_collide,xi_246);
                  const __m256d xi_265 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_17_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_249),xi_264),xia_3_collide);
                  const __m256d xi_266 = _mm256_mul_pd(xi_265,_mm256_set_pd(xi_251,xi_251,xi_251,xi_251));
                  const __m256d xi_269 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_11_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_264),xi_268),xia_4_collide);
                  const __m256d xi_270 = _mm256_mul_pd(xi_269,_mm256_set_pd(xi_251,xi_251,xi_251,xi_251));
                  const __m256d xi_284 = _mm256_mul_pd(xi_265,_mm256_set_pd(xi_283,xi_283,xi_283,xi_283));
                  const __m256d xi_285 = _mm256_mul_pd(xi_284,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_327 = _mm256_mul_pd(xi_269,_mm256_set_pd(xi_283,xi_283,xi_283,xi_283));
                  const __m256d xi_328 = _mm256_mul_pd(xi_327,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d u_1_a_collide = xia_10_collide;
                  const __m256d xi_1 = _mm256_mul_pd(u_1_a_collide,xia_5_collide);
                  const __m256d xi_8 = _mm256_mul_pd(xi_1,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d xi_21 = _mm256_mul_pd(xi_1,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_22 = _mm256_mul_pd(xi_1,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_23 = _mm256_mul_pd(xi_22,_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_24 = _mm256_add_pd(_mm256_mul_pd(xi_21,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_23);
                  const __m256d xi_30 = _mm256_add_pd(xi_14,xi_24);
                  const __m256d xi_43 = _mm256_mul_pd(u_1_a_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_44 = _mm256_mul_pd(xi_43,xia_18_collide);
                  const __m256d xi_48 = _mm256_mul_pd(u_1_a_collide,_mm256_set_pd(xi_45,xi_45,xi_45,xi_45));
                  const __m256d xi_49 = _mm256_mul_pd(xi_48,xia_18_collide);
                  const __m256d xi_50 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_47,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_49,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_42),xi_44);
                  const __m256d xi_52 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_42,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_44,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_47),xi_49);
                  const __m256d xi_59 = _mm256_mul_pd(xi_23,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_63 = _mm256_mul_pd(xi_43,xia_26_collide);
                  const __m256d xi_65 = _mm256_mul_pd(xi_48,xia_26_collide);
                  const __m256d xi_71 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(u_1_a_collide,xia_5_collide),_mm256_set_pd(-0.041666666666666664,-0.041666666666666664,-0.041666666666666664,-0.041666666666666664)),_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_164 = _mm256_mul_pd(u_1_a_collide,u_1_a_collide);
                  const __m256d xi_165 = _mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_164),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_183 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_171,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_178,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_180,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_2_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_6_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_a_collide,xi_164)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a));
                  const __m256d xi_247 = _mm256_mul_pd(u_1_a_collide,xi_246);
                  const __m256d xi_250 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_3_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_247),xi_249),xia_17_collide);
                  const __m256d xi_252 = _mm256_mul_pd(xi_250,_mm256_set_pd(xi_251,xi_251,xi_251,xi_251));
                  const __m256d xi_255 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_15_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_247),xi_254),xia_20_collide);
                  const __m256d xi_256 = _mm256_mul_pd(xi_255,_mm256_set_pd(xi_251,xi_251,xi_251,xi_251));
                  const __m256d xi_286 = _mm256_mul_pd(rho_a_collide,u_1_a_collide);
                  const __m256d xi_288 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_9_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_0_a_collide,xi_286)),xi_177),xi_248),_mm256_set_pd(xi_287,xi_287,xi_287,xi_287));
                  const __m256d xi_289 = _mm256_mul_pd(xi_288,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_300 = _mm256_mul_pd(xi_250,_mm256_set_pd(xi_283,xi_283,xi_283,xi_283));
                  const __m256d xi_311 = _mm256_mul_pd(xi_255,_mm256_set_pd(xi_283,xi_283,xi_283,xi_283));
                  const __m256d xi_321 = _mm256_mul_pd(xi_311,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d u_2_a_collide = xia_14_collide;
                  const __m256d xi_2 = _mm256_mul_pd(u_2_a_collide,xia_26_collide);
                  const __m256d xi_9 = _mm256_mul_pd(xi_2,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_10 = _mm256_mul_pd(xi_2,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_11 = _mm256_mul_pd(xi_10,_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_12 = _mm256_add_pd(_mm256_mul_pd(xi_9,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_11);
                  const __m256d xi_15 = _mm256_add_pd(xi_12,xi_14);
                  const __m256d xi_16 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_8,_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a)),_mm256_mul_pd(_mm256_mul_pd(xi_1,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),xi_15),xi_8);
                  const __m256d xi_25 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_20,_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a)),_mm256_mul_pd(_mm256_mul_pd(xi_0,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),xi_12),xi_20),xi_24);
                  const __m256d xi_29 = _mm256_mul_pd(xi_2,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d xi_31 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_29,_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a)),_mm256_mul_pd(_mm256_mul_pd(xi_2,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),xi_29),xi_30);
                  const __m256d xi_36 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(u_2_a_collide,xia_26_collide),_mm256_set_pd(-0.041666666666666664,-0.041666666666666664,-0.041666666666666664,-0.041666666666666664)),_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
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
                  const __m256d xi_75 = _mm256_mul_pd(u_2_a_collide,xia_18_collide);
                  const __m256d xi_76 = _mm256_mul_pd(xi_75,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_78 = _mm256_mul_pd(xi_75,_mm256_set_pd(xi_45,xi_45,xi_45,xi_45));
                  const __m256d xi_79 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_77,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_78,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_74),xi_76);
                  const __m256d xi_80 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_74,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_76,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_77),xi_78);
                  const __m256d xi_81 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_15,xi_22),xi_51),xi_71);
                  const __m256d xi_168 = _mm256_mul_pd(u_2_a_collide,u_2_a_collide);
                  const __m256d xi_169 = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_168),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331)),xi_167);
                  const __m256d xi_172 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_165,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_169,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_171,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_166),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))),_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_176 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_169,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_173,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_175,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_164),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))),_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_179 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_165,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_167,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_173,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_178,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_168),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))),_mm256_set_pd(omega_even_a,omega_even_a,omega_even_a,omega_even_a));
                  const __m256d xi_184 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_171,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_181,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_19_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_22_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_a_collide,xi_168)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a));
                  const __m256d xi_260 = _mm256_mul_pd(xi_172,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_261 = _mm256_mul_pd(xi_179,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_262 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_183,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_260),xi_261);
                  const __m256d xi_272 = _mm256_mul_pd(xi_176,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_273 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_182,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_261),xi_272);
                  const __m256d xi_275 = _mm256_mul_pd(u_2_a_collide,xi_246);
                  const __m256d xi_276 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_4_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_268),xi_275),xia_11_collide);
                  const __m256d xi_277 = _mm256_mul_pd(xi_276,_mm256_set_pd(xi_251,xi_251,xi_251,xi_251));
                  const __m256d xi_278 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_20_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_254),xi_275),xia_15_collide);
                  const __m256d xi_279 = _mm256_mul_pd(xi_278,_mm256_set_pd(xi_251,xi_251,xi_251,xi_251));
                  const __m256d xi_281 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_184,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_260),xi_272);
                  const __m256d xi_299 = _mm256_mul_pd(xi_179,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_301 = _mm256_add_pd(xi_299,xi_300);
                  const __m256d xi_309 = _mm256_add_pd(_mm256_mul_pd(xi_300,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_299);
                  const __m256d xi_312 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_23_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_2_a_collide,xi_286)),xi_170),xi_253),_mm256_set_pd(xi_287,xi_287,xi_287,xi_287));
                  const __m256d xi_318 = _mm256_mul_pd(xi_172,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_319 = _mm256_mul_pd(xi_278,_mm256_set_pd(xi_283,xi_283,xi_283,xi_283));
                  const __m256d xi_320 = _mm256_add_pd(xi_318,xi_319);
                  const __m256d xi_322 = _mm256_mul_pd(xi_312,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_329 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_8_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,u_0_a_collide),u_2_a_collide)),xi_174),xi_267),_mm256_set_pd(xi_287,xi_287,xi_287,xi_287));
                  const __m256d xi_330 = _mm256_mul_pd(xi_329,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_334 = _mm256_mul_pd(xi_176,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_335 = _mm256_mul_pd(xi_276,_mm256_set_pd(xi_283,xi_283,xi_283,xi_283));
                  const __m256d xi_336 = _mm256_add_pd(xi_334,xi_335);
                  const __m256d xi_341 = _mm256_add_pd(_mm256_mul_pd(xi_319,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_318);
                  const __m256d xi_344 = _mm256_add_pd(_mm256_mul_pd(xi_335,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_334);
                  const __m256d forceTerm_0_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_0,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_1,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_2,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(xi_0,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_3,xi_3,xi_3,xi_3))),_mm256_mul_pd(_mm256_mul_pd(xi_1,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_3,xi_3,xi_3,xi_3))),_mm256_mul_pd(_mm256_mul_pd(xi_2,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_3,xi_3,xi_3,xi_3))),_mm256_mul_pd(_mm256_mul_pd(u_0_a_collide,xia_18_collide),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),_mm256_mul_pd(_mm256_mul_pd(u_1_a_collide,xia_5_collide),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a))),_mm256_mul_pd(_mm256_mul_pd(u_2_a_collide,xia_26_collide),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a)));
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
                  const __m256d xib_1_collide = _mm256_loadu_pd(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3 + ctr_0]);
                  const __m256d xi_86 = _mm256_mul_pd(xib_1_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_89 = _mm256_mul_pd(xi_86,_mm256_set_pd(omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b));
                  const __m256d xi_119 = _mm256_mul_pd(xib_1_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_120 = _mm256_mul_pd(xib_1_collide,_mm256_set_pd(xi_115,xi_115,xi_115,xi_115));
                  const __m256d xi_121 = _mm256_add_pd(_mm256_mul_pd(xi_119,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_120);
                  const __m256d xi_123 = _mm256_mul_pd(xib_1_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_128 = _mm256_mul_pd(xib_1_collide,_mm256_set_pd(xi_127,xi_127,xi_127,xi_127));
                  const __m256d xi_135 = _mm256_add_pd(_mm256_mul_pd(xi_120,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_119);
                  const __m256d xib_2_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_374 = _mm256_mul_pd(xib_2_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xib_3_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_4_collide = _mm256_loadu_pd(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3 + ctr_0]);
                  const __m256d xi_108 = _mm256_mul_pd(xib_4_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_110 = _mm256_mul_pd(xi_108,_mm256_set_pd(omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b));
                  const __m256d xi_137 = _mm256_mul_pd(xib_4_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_138 = _mm256_mul_pd(xib_4_collide,_mm256_set_pd(xi_115,xi_115,xi_115,xi_115));
                  const __m256d xi_139 = _mm256_add_pd(_mm256_mul_pd(xi_138,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_137);
                  const __m256d xi_152 = _mm256_add_pd(_mm256_mul_pd(xi_137,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_138);
                  const __m256d xib_5_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_382 = _mm256_mul_pd(xib_5_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xib_6_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_7_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0]);
                  const __m256d xib_8_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_9_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_383 = _mm256_add_pd(xi_382,xib_9_collide);
                  const __m256d xib_10_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_11_collide = _mm256_loadu_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0]);
                  const __m256d xib_12_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_13_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_14_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_352 = _mm256_add_pd(xib_14_collide,xib_8_collide);
                  const __m256d xib_15_collide = _mm256_loadu_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0]);
                  const __m256d xib_16_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_17_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_375 = _mm256_add_pd(xi_374,xib_17_collide);
                  const __m256d xib_18_collide = _mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const __m256d xib_19_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_356 = _mm256_add_pd(xib_16_collide,xib_19_collide);
                  const __m256d xi_357 = _mm256_add_pd(_mm256_add_pd(xi_356,xib_5_collide),xib_9_collide);
                  const __m256d xib_20_collide = _mm256_load_pd(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + ctr_0]);
                  const __m256d xi_99 = _mm256_mul_pd(xib_20_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_100 = _mm256_mul_pd(xi_99,_mm256_set_pd(omega_odd_b,omega_odd_b,omega_odd_b,omega_odd_b));
                  const __m256d xi_114 = _mm256_mul_pd(xib_20_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_116 = _mm256_mul_pd(xib_20_collide,_mm256_set_pd(xi_115,xi_115,xi_115,xi_115));
                  const __m256d xi_117 = _mm256_add_pd(_mm256_mul_pd(xi_116,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_114);
                  const __m256d xi_133 = _mm256_add_pd(_mm256_mul_pd(xi_114,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_116);
                  const __m256d xib_21_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_22_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_23_collide = _mm256_load_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0]);
                  const __m256d xib_24_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_25_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_353 = _mm256_add_pd(_mm256_add_pd(xi_352,xib_22_collide),xib_25_collide);
                  const __m256d xi_369 = _mm256_mul_pd(xib_25_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_370 = _mm256_add_pd(xi_369,xib_22_collide);
                  const __m256d xib_26_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_359 = _mm256_add_pd(xib_21_collide,xib_26_collide);
                  const __m256d xi_360 = _mm256_add_pd(_mm256_add_pd(xi_359,xib_17_collide),xib_2_collide);
                  const __m256d xi_88 = _mm256_mul_pd(xib_1_collide,_mm256_set_pd(xi_87,xi_87,xi_87,xi_87));
                  const __m256d xi_101 = _mm256_mul_pd(xib_20_collide,_mm256_set_pd(xi_87,xi_87,xi_87,xi_87));
                  const __m256d xi_109 = _mm256_mul_pd(xib_4_collide,_mm256_set_pd(xi_87,xi_87,xi_87,xi_87));
                  const __m256d rho_b_collide = xib_18_collide;
                  const __m256d xi_349 = _mm256_mul_pd(rho_b_collide,_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111));
                  const __m256d xi_362 = _mm256_mul_pd(rho_b_collide,_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_363 = _mm256_add_pd(xi_360,xi_362);
                  const __m256d xi_367 = _mm256_mul_pd(rho_b_collide,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d u_0_b_collide = xib_23_collide;
                  const __m256d xi_82 = _mm256_mul_pd(u_0_b_collide,xib_20_collide);
                  const __m256d xi_95 = _mm256_mul_pd(xi_82,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_96 = _mm256_add_pd(_mm256_mul_pd(xi_82,_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666)),_mm256_mul_pd(xi_95,_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b)));
                  const __m256d xi_102 = _mm256_mul_pd(xi_82,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d xi_124 = _mm256_mul_pd(u_0_b_collide,xi_123);
                  const __m256d xi_129 = _mm256_mul_pd(u_0_b_collide,xi_128);
                  const __m256d xi_140 = _mm256_mul_pd(xi_95,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_143 = _mm256_mul_pd(_mm256_mul_pd(xi_82,_mm256_set_pd(0.041666666666666664,0.041666666666666664,0.041666666666666664,0.041666666666666664)),_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_155 = _mm256_mul_pd(u_0_b_collide,xib_4_collide);
                  const __m256d xi_156 = _mm256_mul_pd(xi_155,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_159 = _mm256_mul_pd(xi_155,_mm256_set_pd(xi_127,xi_127,xi_127,xi_127));
                  const __m256d xi_348 = _mm256_mul_pd(u_0_b_collide,u_0_b_collide);
                  const __m256d xi_355 = _mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_348),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_364 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_357,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_363,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_12_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_6_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_b_collide,xi_348)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b));
                  const __m256d xi_381 = _mm256_mul_pd(u_0_b_collide,xi_367);
                  const __m256d xi_384 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_19_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_381),xi_383),xib_16_collide);
                  const __m256d xi_385 = _mm256_mul_pd(xi_384,_mm256_set_pd(xi_372,xi_372,xi_372,xi_372));
                  const __m256d xi_386 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_21_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_375),xi_381),xib_26_collide);
                  const __m256d xi_387 = _mm256_mul_pd(xi_386,_mm256_set_pd(xi_372,xi_372,xi_372,xi_372));
                  const __m256d xi_397 = _mm256_mul_pd(xi_386,_mm256_set_pd(xi_396,xi_396,xi_396,xi_396));
                  const __m256d xi_398 = _mm256_mul_pd(xi_397,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_414 = _mm256_mul_pd(xi_384,_mm256_set_pd(xi_396,xi_396,xi_396,xi_396));
                  const __m256d xi_415 = _mm256_mul_pd(xi_414,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d u_1_b_collide = xib_11_collide;
                  const __m256d xi_83 = _mm256_mul_pd(u_1_b_collide,xib_1_collide);
                  const __m256d xi_90 = _mm256_mul_pd(xi_83,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d xi_103 = _mm256_mul_pd(xi_83,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_104 = _mm256_mul_pd(xi_83,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_105 = _mm256_mul_pd(xi_104,_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_106 = _mm256_add_pd(_mm256_mul_pd(xi_103,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_105);
                  const __m256d xi_112 = _mm256_add_pd(xi_106,xi_96);
                  const __m256d xi_125 = _mm256_mul_pd(u_1_b_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_126 = _mm256_mul_pd(xi_125,xib_20_collide);
                  const __m256d xi_130 = _mm256_mul_pd(u_1_b_collide,_mm256_set_pd(xi_127,xi_127,xi_127,xi_127));
                  const __m256d xi_131 = _mm256_mul_pd(xi_130,xib_20_collide);
                  const __m256d xi_132 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_129,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_131,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_124),xi_126);
                  const __m256d xi_134 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_124,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_126,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_129),xi_131);
                  const __m256d xi_141 = _mm256_mul_pd(xi_105,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_145 = _mm256_mul_pd(xi_125,xib_4_collide);
                  const __m256d xi_147 = _mm256_mul_pd(xi_130,xib_4_collide);
                  const __m256d xi_153 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(u_1_b_collide,xib_1_collide),_mm256_set_pd(-0.041666666666666664,-0.041666666666666664,-0.041666666666666664,-0.041666666666666664)),_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_346 = _mm256_mul_pd(u_1_b_collide,u_1_b_collide);
                  const __m256d xi_347 = _mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_346),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_365 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_353,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_363,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_10_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_3_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_b_collide,xi_346)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b));
                  const __m256d xi_368 = _mm256_mul_pd(u_1_b_collide,xi_367);
                  const __m256d xi_371 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_8_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_368),xi_370),xib_14_collide);
                  const __m256d xi_373 = _mm256_mul_pd(xi_371,_mm256_set_pd(xi_372,xi_372,xi_372,xi_372));
                  const __m256d xi_376 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_26_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_368),xi_375),xib_21_collide);
                  const __m256d xi_377 = _mm256_mul_pd(xi_376,_mm256_set_pd(xi_372,xi_372,xi_372,xi_372));
                  const __m256d xi_399 = _mm256_mul_pd(rho_b_collide,u_1_b_collide);
                  const __m256d xi_401 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_17_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_0_b_collide,xi_399)),xi_359),xi_374),_mm256_set_pd(xi_400,xi_400,xi_400,xi_400));
                  const __m256d xi_402 = _mm256_mul_pd(xi_401,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_404 = _mm256_mul_pd(xi_376,_mm256_set_pd(xi_396,xi_396,xi_396,xi_396));
                  const __m256d xi_407 = _mm256_mul_pd(xi_371,_mm256_set_pd(xi_396,xi_396,xi_396,xi_396));
                  const __m256d xi_412 = _mm256_mul_pd(xi_407,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d u_2_b_collide = xib_15_collide;
                  const __m256d xi_84 = _mm256_mul_pd(u_2_b_collide,xib_4_collide);
                  const __m256d xi_91 = _mm256_mul_pd(xi_84,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_92 = _mm256_mul_pd(xi_84,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_93 = _mm256_mul_pd(xi_92,_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_94 = _mm256_add_pd(_mm256_mul_pd(xi_91,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_93);
                  const __m256d xi_97 = _mm256_add_pd(xi_94,xi_96);
                  const __m256d xi_98 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_90,_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b)),_mm256_mul_pd(_mm256_mul_pd(xi_83,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),xi_90),xi_97);
                  const __m256d xi_107 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_102,_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b)),_mm256_mul_pd(_mm256_mul_pd(xi_82,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),xi_102),xi_106),xi_94);
                  const __m256d xi_111 = _mm256_mul_pd(xi_84,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d xi_113 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_111,_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b)),_mm256_mul_pd(_mm256_mul_pd(xi_84,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),xi_111),xi_112);
                  const __m256d xi_118 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(u_2_b_collide,xib_4_collide),_mm256_set_pd(-0.041666666666666664,-0.041666666666666664,-0.041666666666666664,-0.041666666666666664)),_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
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
                  const __m256d xi_157 = _mm256_mul_pd(u_2_b_collide,xib_20_collide);
                  const __m256d xi_158 = _mm256_mul_pd(xi_157,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_160 = _mm256_mul_pd(xi_157,_mm256_set_pd(xi_127,xi_127,xi_127,xi_127));
                  const __m256d xi_161 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_159,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_160,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_156),xi_158);
                  const __m256d xi_162 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_156,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_158,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_159),xi_160);
                  const __m256d xi_163 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_104,xi_133),xi_153),xi_97);
                  const __m256d xi_350 = _mm256_mul_pd(u_2_b_collide,u_2_b_collide);
                  const __m256d xi_351 = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_350),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331)),xi_349);
                  const __m256d xi_354 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_347,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_351,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_353,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_348),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))),_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_358 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_351,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_355,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_357,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_346),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))),_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_361 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_347,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_349,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_355,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_360,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_350),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))),_mm256_set_pd(omega_even_b,omega_even_b,omega_even_b,omega_even_b));
                  const __m256d xi_366 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_353,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_357,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_362,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_13_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_24_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_b_collide,xi_350)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b));
                  const __m256d xi_378 = _mm256_mul_pd(xi_354,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_379 = _mm256_mul_pd(xi_361,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_380 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_365,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_378),xi_379);
                  const __m256d xi_388 = _mm256_mul_pd(xi_358,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_389 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_364,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_379),xi_388);
                  const __m256d xi_390 = _mm256_mul_pd(u_2_b_collide,xi_367);
                  const __m256d xi_391 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_16_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_383),xi_390),xib_19_collide);
                  const __m256d xi_392 = _mm256_mul_pd(xi_391,_mm256_set_pd(xi_372,xi_372,xi_372,xi_372));
                  const __m256d xi_393 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_14_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_370),xi_390),xib_8_collide);
                  const __m256d xi_394 = _mm256_mul_pd(xi_393,_mm256_set_pd(xi_372,xi_372,xi_372,xi_372));
                  const __m256d xi_395 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_366,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_378),xi_388);
                  const __m256d xi_403 = _mm256_mul_pd(xi_361,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_405 = _mm256_add_pd(xi_403,xi_404);
                  const __m256d xi_406 = _mm256_add_pd(_mm256_mul_pd(xi_404,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_403);
                  const __m256d xi_408 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_22_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_2_b_collide,xi_399)),xi_352),xi_369),_mm256_set_pd(xi_400,xi_400,xi_400,xi_400));
                  const __m256d xi_409 = _mm256_mul_pd(xi_354,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_410 = _mm256_mul_pd(xi_393,_mm256_set_pd(xi_396,xi_396,xi_396,xi_396));
                  const __m256d xi_411 = _mm256_add_pd(xi_409,xi_410);
                  const __m256d xi_413 = _mm256_mul_pd(xi_408,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_416 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_9_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,u_0_b_collide),u_2_b_collide)),xi_356),xi_382),_mm256_set_pd(xi_400,xi_400,xi_400,xi_400));
                  const __m256d xi_417 = _mm256_mul_pd(xi_416,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_418 = _mm256_mul_pd(xi_358,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_419 = _mm256_mul_pd(xi_391,_mm256_set_pd(xi_396,xi_396,xi_396,xi_396));
                  const __m256d xi_420 = _mm256_add_pd(xi_418,xi_419);
                  const __m256d xi_421 = _mm256_add_pd(_mm256_mul_pd(xi_410,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_409);
                  const __m256d xi_422 = _mm256_add_pd(_mm256_mul_pd(xi_419,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_418);
                  const __m256d forceTerm_0_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_82,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_83,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_84,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(xi_82,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_85,xi_85,xi_85,xi_85))),_mm256_mul_pd(_mm256_mul_pd(xi_83,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_85,xi_85,xi_85,xi_85))),_mm256_mul_pd(_mm256_mul_pd(xi_84,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_85,xi_85,xi_85,xi_85))),_mm256_mul_pd(_mm256_mul_pd(u_0_b_collide,xib_20_collide),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),_mm256_mul_pd(_mm256_mul_pd(u_1_b_collide,xib_1_collide),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b))),_mm256_mul_pd(_mm256_mul_pd(u_2_b_collide,xib_4_collide),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b)));
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
                  const __m256d tmp_a0 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_182,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_183,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_184,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_0_a_collide),xi_172),xi_176),xi_179),xi_245),xia_25_collide);
                  const __m256d tmp_a1 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_252,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_256,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_1_a_collide),xi_259),xi_262),xia_6_collide);
                  const __m256d tmp_a2 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_2_a_collide,xi_252),xi_256),xi_262),xi_263),xia_2_collide);
                  const __m256d tmp_a3 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_3_a_collide,xi_266),xi_270),xi_271),xi_273),xia_16_collide);
                  const __m256d tmp_a4 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_266,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_270,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_4_a_collide),xi_273),xi_274),xia_12_collide);
                  const __m256d tmp_a5 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_277,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_279,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_5_a_collide),xi_280),xi_281),xia_19_collide);
                  const __m256d tmp_a6 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_6_a_collide,xi_277),xi_279),xi_281),xi_282),xia_22_collide);
                  const __m256d tmp_a7 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_7_a_collide,xi_285),xi_289),xi_298),xi_301),xia_3_collide);
                  const __m256d tmp_a8 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_8_a_collide,xi_284),xi_288),xi_301),xi_307),xia_21_collide);
                  const __m256d tmp_a9 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_9_a_collide,xi_285),xi_288),xi_308),xi_309),xia_9_collide);
                  const __m256d tmp_a10 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_10_a_collide,xi_284),xi_289),xi_309),xi_310),xia_17_collide);
                  const __m256d tmp_a11 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_11_a_collide,xi_311),xi_312),xi_317),xi_320),xia_1_collide);
                  const __m256d tmp_a12 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_12_a_collide,xi_320),xi_321),xi_322),xi_326),xia_20_collide);
                  const __m256d tmp_a13 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_13_a_collide,xi_328),xi_330),xi_333),xi_336),xia_4_collide);
                  const __m256d tmp_a14 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_14_a_collide,xi_327),xi_329),xi_336),xi_339),xia_13_collide);
                  const __m256d tmp_a15 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_15_a_collide,xi_311),xi_322),xi_340),xi_341),xia_15_collide);
                  const __m256d tmp_a16 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_16_a_collide,xi_312),xi_321),xi_341),xi_342),xia_23_collide);
                  const __m256d tmp_a17 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_17_a_collide,xi_328),xi_329),xi_343),xi_344),xia_8_collide);
                  const __m256d tmp_a18 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_18_a_collide,xi_327),xi_330),xi_344),xi_345),xia_11_collide);
                  const __m256d tmp_b0 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_364,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_365,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_366,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_0_b_collide),xi_245),xi_354),xi_358),xi_361),xib_7_collide);
                  const __m256d tmp_b1 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_373,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_377,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_1_b_collide),xi_259),xi_380),xib_10_collide);
                  const __m256d tmp_b2 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_2_b_collide,xi_263),xi_373),xi_377),xi_380),xib_3_collide);
                  const __m256d tmp_b3 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_3_b_collide,xi_271),xi_385),xi_387),xi_389),xib_12_collide);
                  const __m256d tmp_b4 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_385,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_387,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_4_b_collide),xi_274),xi_389),xib_6_collide);
                  const __m256d tmp_b5 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_392,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_394,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_5_b_collide),xi_280),xi_395),xib_13_collide);
                  const __m256d tmp_b6 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_6_b_collide,xi_282),xi_392),xi_394),xi_395),xib_24_collide);
                  const __m256d tmp_b7 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_7_b_collide,xi_298),xi_398),xi_402),xi_405),xib_26_collide);
                  const __m256d tmp_b8 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_8_b_collide,xi_307),xi_397),xi_401),xi_405),xib_2_collide);
                  const __m256d tmp_b9 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_9_b_collide,xi_308),xi_398),xi_401),xi_406),xib_17_collide);
                  const __m256d tmp_b10 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_10_b_collide,xi_310),xi_397),xi_402),xi_406),xib_21_collide);
                  const __m256d tmp_b11 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_11_b_collide,xi_317),xi_407),xi_408),xi_411),xib_25_collide);
                  const __m256d tmp_b12 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_12_b_collide,xi_326),xi_411),xi_412),xi_413),xib_14_collide);
                  const __m256d tmp_b13 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_13_b_collide,xi_333),xi_415),xi_417),xi_420),xib_16_collide);
                  const __m256d tmp_b14 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_14_b_collide,xi_339),xi_414),xi_416),xi_420),xib_5_collide);
                  const __m256d tmp_b15 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_15_b_collide,xi_340),xi_407),xi_413),xi_421),xib_8_collide);
                  const __m256d tmp_b16 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_16_b_collide,xi_342),xi_408),xi_412),xi_421),xib_22_collide);
                  const __m256d tmp_b17 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_17_b_collide,xi_343),xi_415),xi_416),xi_422),xib_9_collide);
                  const __m256d tmp_b18 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_18_b_collide,xi_345),xi_414),xi_417),xi_422),xib_19_collide);
                  const __m256d xirecolor_0 = _mm256_add_pd(tmp_a0,tmp_b0);
                  const __m256d xirecolor_1 = _mm256_add_pd(_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]),_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_2 = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),xirecolor_1);
                  const __m256d xi_444 = _mm256_mul_pd(xirecolor_2,_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_3 = _mm256_mul_pd(xirecolor_2,_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_4 = _mm256_add_pd(tmp_a1,tmp_b1);
                  const __m256d xirecolor_5 = xi_195;
                  const __m256d xirecolor_6 = xi_193;
                  const __m256d xirecolor_7 = xi_213;
                  const __m256d xi_423 = _mm256_mul_pd(xirecolor_7,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_8 = xi_423;
                  const __m256d xirecolor_9 = xi_200;
                  const __m256d xirecolor_10 = _mm256_mul_pd(xirecolor_9,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_11 = xi_202;
                  const __m256d xi_424 = _mm256_add_pd(xirecolor_10,xirecolor_11);
                  const __m256d xirecolor_12 = xi_215;
                  const __m256d xi_425 = _mm256_add_pd(xirecolor_12,xirecolor_8);
                  const __m256d xirecolor_13 = _mm256_add_pd(xi_424,xi_425);
                  const __m256d xirecolor_14 = xi_217;
                  const __m256d xirecolor_15 = xi_218;
                  const __m256d xi_426 = _mm256_add_pd(_mm256_mul_pd(xirecolor_15,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_14);
                  const __m256d xirecolor_16 = _mm256_mul_pd(xi_426,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_17 = xi_187;
                  const __m256d xirecolor_18 = xi_188;
                  const __m256d xirecolor_19 = xi_189;
                  const __m256d xi_427 = _mm256_add_pd(xirecolor_18,xirecolor_19);
                  const __m256d xirecolor_20 = xi_190;
                  const __m256d xirecolor_21 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_17,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_427),xirecolor_20);
                  const __m256d xirecolor_22 = xi_186;
                  const __m256d xirecolor_23 = xi_185;
                  const __m256d xi_428 = _mm256_add_pd(_mm256_mul_pd(xirecolor_23,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_22);
                  const __m256d xirecolor_24 = _mm256_mul_pd(xi_428,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_25 = xi_221;
                  const __m256d xirecolor_26 = _mm256_mul_pd(xirecolor_25,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_27 = xi_207;
                  const __m256d xirecolor_28 = xi_206;
                  const __m256d xirecolor_29 = xi_223;
                  const __m256d xi_429 = _mm256_add_pd(xirecolor_26,xirecolor_29);
                  const __m256d xirecolor_30 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_27,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_429),xirecolor_28);
                  const __m256d xirecolor_31 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_5,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_13),xirecolor_16),xirecolor_21),xirecolor_24),xirecolor_30),xirecolor_6);
                  const __m256d xirecolor_32 = xi_192;
                  const __m256d xirecolor_33 = xi_194;
                  const __m256d xi_430 = _mm256_add_pd(xirecolor_32,xirecolor_33);
                  const __m256d xirecolor_34 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_6,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_430),xirecolor_5);
                  const __m256d xirecolor_35 = xi_197;
                  const __m256d xirecolor_36 = xi_198;
                  const __m256d xi_431 = _mm256_add_pd(_mm256_mul_pd(xirecolor_36,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_35);
                  const __m256d xirecolor_37 = _mm256_mul_pd(xi_431,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_38 = xi_424;
                  const __m256d xirecolor_39 = xi_204;
                  const __m256d xirecolor_40 = xi_205;
                  const __m256d xi_432 = _mm256_add_pd(xirecolor_39,xirecolor_40);
                  const __m256d xirecolor_41 = _mm256_add_pd(_mm256_mul_pd(xirecolor_28,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_27);
                  const __m256d xirecolor_42 = _mm256_add_pd(xi_432,xirecolor_41);
                  const __m256d xirecolor_43 = _mm256_add_pd(xirecolor_38,xirecolor_42);
                  const __m256d xirecolor_44 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_428,xirecolor_21),xirecolor_34),xirecolor_37),xirecolor_43);
                  const __m256d xirecolor_45 = xi_426;
                  const __m256d xirecolor_46 = xi_425;
                  const __m256d xirecolor_47 = xi_229;
                  const __m256d xirecolor_48 = _mm256_mul_pd(xirecolor_47,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_49 = xi_231;
                  const __m256d xi_433 = _mm256_add_pd(xirecolor_48,xirecolor_49);
                  const __m256d xirecolor_50 = _mm256_add_pd(xi_431,xi_433);
                  const __m256d xirecolor_51 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_20,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_427),xirecolor_17),xirecolor_24),xirecolor_34),xirecolor_45),xirecolor_46),xirecolor_50);
                  const __m256d xirecolor_52 = _mm256_sqrt_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_31,xirecolor_31),_mm256_mul_pd(xirecolor_44,xirecolor_44)),_mm256_mul_pd(xirecolor_51,xirecolor_51)));
                  const __m256d xirecolor_53 = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),xirecolor_52);
                  const __m256d xi_434 = _mm256_mul_pd(xirecolor_31,xirecolor_53);
                  const __m256d xi_435 = _mm256_mul_pd(xirecolor_51,xirecolor_53);
                  const __m256d xi_436 = _mm256_mul_pd(xirecolor_44,xirecolor_53);
                  const __m256d xirecolor_54 = _mm256_cmp_pd(xirecolor_52,_mm256_set_pd(0.0,0.0,0.0,0.0),_CMP_NLE_UQ);
                  const __m256d xirecolor_55 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set_pd(beta,beta,beta,beta),_mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),_mm256_mul_pd(xirecolor_1,xirecolor_1))),_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_56 = _mm256_mul_pd(xirecolor_55,_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]))));
                  const __m256d xirecolor_57 = _mm256_mul_pd(xirecolor_56,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_434,xirecolor_54));
                  const __m256d xirecolor_58 = _mm256_add_pd(tmp_a2,tmp_b2);
                  const __m256d xirecolor_59 = _mm256_mul_pd(xirecolor_56,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_434,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_60 = _mm256_add_pd(tmp_a3,tmp_b3);
                  const __m256d xirecolor_61 = _mm256_mul_pd(xirecolor_56,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_435,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_62 = _mm256_add_pd(tmp_a4,tmp_b4);
                  const __m256d xirecolor_63 = _mm256_mul_pd(xirecolor_56,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_435,xirecolor_54));
                  const __m256d xirecolor_64 = _mm256_add_pd(tmp_a5,tmp_b5);
                  const __m256d xirecolor_65 = _mm256_mul_pd(xirecolor_56,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_436,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_66 = _mm256_add_pd(tmp_a6,tmp_b6);
                  const __m256d xirecolor_67 = _mm256_mul_pd(xirecolor_56,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_436,xirecolor_54));
                  const __m256d xirecolor_68 = _mm256_add_pd(tmp_a7,tmp_b7);
                  const __m256d xirecolor_69 = xi_290;
                  const __m256d xirecolor_70 = xi_291;
                  const __m256d xirecolor_71 = xi_430;
                  const __m256d xirecolor_72 = _mm256_add_pd(xirecolor_50,xirecolor_71);
                  const __m256d xirecolor_73 = xi_295;
                  const __m256d xirecolor_74 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_29,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_25),xirecolor_73);
                  const __m256d xirecolor_75 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_11,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xirecolor_69,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_292),xirecolor_41),xirecolor_70),xirecolor_72),xirecolor_74),xirecolor_9);
                  const __m256d xirecolor_76 = _mm256_mul_pd(xirecolor_53,_mm256_set_pd(0.70710678118654757,0.70710678118654757,0.70710678118654757,0.70710678118654757));
                  const __m256d xi_437 = _mm256_mul_pd(xirecolor_75,xirecolor_76);
                  const __m256d xirecolor_77 = _mm256_mul_pd(xirecolor_55,_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]))));
                  const __m256d xirecolor_78 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_437,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_79 = _mm256_add_pd(tmp_a8,tmp_b8);
                  const __m256d xirecolor_80 = xi_302;
                  const __m256d xirecolor_81 = xi_303;
                  const __m256d xi_438 = _mm256_add_pd(_mm256_mul_pd(xirecolor_80,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_81);
                  const __m256d xirecolor_82 = xi_304;
                  const __m256d xirecolor_83 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_305,xi_438),xirecolor_30),xirecolor_38),xirecolor_72),xirecolor_82);
                  const __m256d xi_439 = _mm256_mul_pd(xirecolor_76,xirecolor_83);
                  const __m256d xirecolor_84 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_439,xirecolor_54));
                  const __m256d xirecolor_85 = _mm256_add_pd(tmp_a9,tmp_b9);
                  const __m256d xirecolor_86 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_439,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_87 = _mm256_add_pd(tmp_a10,tmp_b10);
                  const __m256d xirecolor_88 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_437,xirecolor_54));
                  const __m256d xirecolor_89 = _mm256_add_pd(tmp_a11,tmp_b11);
                  const __m256d xirecolor_90 = _mm256_add_pd(_mm256_add_pd(xi_432,xirecolor_37),xirecolor_71);
                  const __m256d xirecolor_91 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_423,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_438,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xirecolor_12,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d xirecolor_92 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_313,xirecolor_45),xirecolor_74),xirecolor_90),xirecolor_91);
                  const __m256d xi_440 = _mm256_mul_pd(xirecolor_76,xirecolor_92);
                  const __m256d xirecolor_93 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_440,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_94 = _mm256_add_pd(tmp_a12,tmp_b12);
                  const __m256d xirecolor_95 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_70,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_16),xirecolor_69);
                  const __m256d xirecolor_96 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_323,xi_429),xirecolor_46),xirecolor_82),xirecolor_90),xirecolor_95);
                  const __m256d xi_441 = _mm256_mul_pd(xirecolor_76,xirecolor_96);
                  const __m256d xirecolor_97 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_441,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_98 = _mm256_add_pd(tmp_a13,tmp_b13);
                  const __m256d xirecolor_99 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_331,xi_433),xirecolor_13),xirecolor_42),xirecolor_45),xirecolor_73),xirecolor_82);
                  const __m256d xi_442 = _mm256_mul_pd(xirecolor_76,xirecolor_99);
                  const __m256d xirecolor_100 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_442,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_101 = _mm256_add_pd(tmp_a14,tmp_b14);
                  const __m256d xirecolor_102 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_49,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_337),xirecolor_43),xirecolor_47),xirecolor_91),xirecolor_95);
                  const __m256d xi_443 = _mm256_mul_pd(xirecolor_102,xirecolor_76);
                  const __m256d xirecolor_103 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_443,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_104 = _mm256_add_pd(tmp_a15,tmp_b15);
                  const __m256d xirecolor_105 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_441,xirecolor_54));
                  const __m256d xirecolor_106 = _mm256_add_pd(tmp_a16,tmp_b16);
                  const __m256d xirecolor_107 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_440,xirecolor_54));
                  const __m256d xirecolor_108 = _mm256_add_pd(tmp_a17,tmp_b17);
                  const __m256d xirecolor_109 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_443,xirecolor_54));
                  const __m256d xirecolor_110 = _mm256_add_pd(tmp_a18,tmp_b18);
                  const __m256d xirecolor_111 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_442,xirecolor_54));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0],_mm256_mul_pd(xirecolor_0,xirecolor_3));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_4),xirecolor_57));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_58),xirecolor_59));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_60),xirecolor_61));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_62),xirecolor_63));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_64),xirecolor_65));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_66),xirecolor_67));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_68),xirecolor_78));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_79),xirecolor_84));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_85),xirecolor_86));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_87),xirecolor_88));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_89),xirecolor_93));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_94),xirecolor_97));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_98),xirecolor_100));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_101,xirecolor_3),xirecolor_103));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_104,xirecolor_3),xirecolor_105));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_106,xirecolor_3),xirecolor_107));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_108,xirecolor_3),xirecolor_109));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_110,xirecolor_3),xirecolor_111));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0],_mm256_mul_pd(xi_444,xirecolor_0));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_57,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_4)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_59,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_58)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_61,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_60)));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_63,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_62)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_65,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_64)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_67,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_66)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_78,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_68)));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_84,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_79)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_86,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_85)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_88,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_87)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_93,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_89)));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_97,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_94)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_100,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_98)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_103,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_101)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_105,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_104)));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_107,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_106)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_109,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_108)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_111,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_444,xirecolor_110)));
               }
               for (int64_t ctr_0 = (int64_t)((_size_force_a_0 - 2) / (4)) * (4) + 1; ctr_0 < _size_force_a_0 - 1; ctr_0 += 1)
               {
                  const double xi_185 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_186 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_187 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_188 = -0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_189 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_190 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_191 = -xi_187 + xi_188 + xi_189 + xi_190;
                  const double xi_192 = -0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_193 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_194 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_195 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_196 = xi_192 - xi_193 + xi_194 + xi_195;
                  const double xi_197 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_198 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_199 = -xi_197 + xi_198;
                  const double xi_200 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0];
                  const double xi_201 = -xi_200;
                  const double xi_202 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const double xi_203 = xi_201 + xi_202;
                  const double xi_204 = -0.22222222222222221*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0];
                  const double xi_205 = 0.22222222222222221*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const double xi_206 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0];
                  const double xi_207 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const double xi_208 = -xi_206 + xi_207;
                  const double xi_209 = xi_204 + xi_205 + xi_208;
                  const double xi_210 = xi_203 + xi_209;
                  const double xi_211 = -xi_185 + xi_186 + xi_191 + xi_196 + xi_199 + xi_210;
                  const double xi_212 = (xi_211*xi_211);
                  const double xi_213 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1];
                  const double xi_214 = -xi_213;
                  const double xi_215 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const double xi_216 = xi_201 + xi_202 + xi_214 + xi_215;
                  const double xi_217 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const double xi_218 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1];
                  const double xi_219 = -xi_217 + xi_218;
                  const double xi_220 = xi_185 - xi_186;
                  const double xi_221 = 0.22222222222222221*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0];
                  const double xi_222 = -xi_221;
                  const double xi_223 = 0.22222222222222221*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0];
                  const double xi_224 = xi_206 - xi_207 + xi_222 + xi_223;
                  const double xi_225 = xi_191 + xi_193 - xi_195 + xi_216 + xi_219 + xi_220 + xi_224;
                  const double xi_226 = (xi_225*xi_225);
                  const double xi_227 = xi_217 - xi_218;
                  const double xi_228 = xi_214 + xi_215;
                  const double xi_229 = 0.22222222222222221*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1];
                  const double xi_230 = -xi_229;
                  const double xi_231 = 0.22222222222222221*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const double xi_232 = xi_197 - xi_198 + xi_230 + xi_231;
                  const double xi_233 = xi_187 + xi_188 + xi_189 - xi_190 + xi_196 + xi_220 + xi_227 + xi_228 + xi_232;
                  const double xi_234 = (xi_233*xi_233);
                  const double xi_235 = xi_212 + xi_226 + xi_234;
                  const double xi_236 = pow(xi_235, 0.5);
                  const double xi_241 = (_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]);
                  const double xi_243 = sigma*xi_236*((0.5 < _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (omega_shear_a): ((-0.5 > _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (omega_shear_b): ((0.0 < _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (xi_238 + xi_240*xi_241 - xi_240*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]): (xi_238 + xi_241*xi_242 + xi_242*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))));
                  const bool xi_244 = xi_236 > 0.0;
                  const double xi_245 = ((xi_244) ? (xi_243*0.25): (0.0));
                  const double xi_257 = ((1.0) / (xi_235));
                  const double xi_258 = xi_243*1.125;
                  const double xi_259 = ((xi_244) ? (xi_258*(xi_226*xi_257*0.055555555555555552 - 0.018518518518518517)): (0.0));
                  const double xi_263 = ((xi_244) ? (xi_258*(xi_257*0.055555555555555552*(xi_225*xi_225) - 0.018518518518518517)): (0.0));
                  const double xi_271 = ((xi_244) ? (xi_258*(xi_257*0.055555555555555552*(xi_233*xi_233) - 0.018518518518518517)): (0.0));
                  const double xi_274 = ((xi_244) ? (xi_258*(xi_234*xi_257*0.055555555555555552 - 0.018518518518518517)): (0.0));
                  const double xi_280 = ((xi_244) ? (xi_258*(xi_257*0.055555555555555552*(xi_211*xi_211) - 0.018518518518518517)): (0.0));
                  const double xi_282 = ((xi_244) ? (xi_258*(xi_212*xi_257*0.055555555555555552 - 0.018518518518518517)): (0.0));
                  const double xi_290 = 0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_291 = 0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_292 = -0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1] + 0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const double xi_293 = xi_192 + xi_194;
                  const double xi_294 = xi_232 + xi_293;
                  const double xi_295 = -0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1] + 0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_296 = xi_221 - xi_223 + xi_295;
                  const double xi_297 = xi_200 - xi_202 + xi_208 - xi_290 + xi_291 + xi_292 + xi_294 + xi_296;
                  const double xi_298 = ((xi_244) ? (xi_258*(xi_257*0.027777777777777776*(xi_297*xi_297) - 0.037037037037037035)): (0.0));
                  const double xi_302 = 0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_303 = 0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_304 = -0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1] + 0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_305 = -0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1] + 0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const double xi_306 = xi_203 + xi_224 + xi_294 - xi_302 + xi_303 + xi_304 + xi_305;
                  const double xi_307 = ((xi_244) ? (xi_258*(xi_257*0.027777777777777776*(xi_306*xi_306) - 0.037037037037037035)): (0.0));
                  const double xi_308 = ((xi_244) ? (xi_258*(xi_257*0.027777777777777776*(xi_306*xi_306) - 0.037037037037037035)): (0.0));
                  const double xi_310 = ((xi_244) ? (xi_258*(xi_257*0.027777777777777776*(xi_297*xi_297) - 0.037037037037037035)): (0.0));
                  const double xi_313 = -0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0] + 0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const double xi_314 = xi_199 + xi_204 + xi_205 + xi_293;
                  const double xi_315 = xi_213 - xi_215 + xi_302 - xi_303;
                  const double xi_316 = xi_227 + xi_296 + xi_313 + xi_314 + xi_315;
                  const double xi_317 = ((xi_244) ? (xi_258*(xi_257*0.027777777777777776*(xi_316*xi_316) - 0.037037037037037035)): (0.0));
                  const double xi_323 = -0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0] + 0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const double xi_324 = xi_219 + xi_290 - xi_291;
                  const double xi_325 = xi_222 + xi_223 + xi_228 + xi_304 + xi_314 + xi_323 + xi_324;
                  const double xi_326 = ((xi_244) ? (xi_258*(xi_257*0.027777777777777776*(xi_325*xi_325) - 0.037037037037037035)): (0.0));
                  const double xi_331 = -0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1] + 0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_332 = xi_209 + xi_216 + xi_227 + xi_230 + xi_231 + xi_295 + xi_304 + xi_331;
                  const double xi_333 = ((xi_244) ? (xi_258*(xi_257*0.027777777777777776*(xi_332*xi_332) - 0.037037037037037035)): (0.0));
                  const double xi_337 = -0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1] + 0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_338 = xi_210 + xi_229 - xi_231 + xi_315 + xi_324 + xi_337;
                  const double xi_339 = ((xi_244) ? (xi_258*(xi_257*0.027777777777777776*(xi_338*xi_338) - 0.037037037037037035)): (0.0));
                  const double xi_340 = ((xi_244) ? (xi_258*(xi_257*0.027777777777777776*(xi_325*xi_325) - 0.037037037037037035)): (0.0));
                  const double xi_342 = ((xi_244) ? (xi_258*(xi_257*0.027777777777777776*(xi_316*xi_316) - 0.037037037037037035)): (0.0));
                  const double xi_343 = ((xi_244) ? (xi_258*(xi_257*0.027777777777777776*(xi_338*xi_338) - 0.037037037037037035)): (0.0));
                  const double xi_345 = ((xi_244) ? (xi_258*(xi_257*0.027777777777777776*(xi_332*xi_332) - 0.037037037037037035)): (0.0));
                  const double xia_1_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0];
                  const double xi_253 = -xia_1_collide;
                  const double xia_2_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0];
                  const double xia_3_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0];
                  const double xia_4_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0];
                  const double xia_5_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3 + ctr_0];
                  const double xi_4 = xia_5_collide*0.16666666666666666;
                  const double xi_7 = omega_odd_a*xi_4;
                  const double xi_37 = xia_5_collide*0.083333333333333329;
                  const double xi_38 = xi_33*xia_5_collide;
                  const double xi_39 = -xi_37 + xi_38;
                  const double xi_41 = xia_5_collide*0.25;
                  const double xi_46 = xi_45*xia_5_collide;
                  const double xi_53 = xi_37 - xi_38;
                  const double xia_6_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0];
                  const double xia_7_collide = _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0];
                  const double xia_8_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0];
                  const double xia_9_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0];
                  const double xia_10_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0];
                  const double xia_11_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0];
                  const double xi_174 = xia_11_collide + xia_4_collide;
                  const double xia_12_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0];
                  const double xia_13_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0];
                  const double xi_175 = xi_174 + xia_13_collide + xia_8_collide;
                  const double xi_267 = -xia_13_collide;
                  const double xi_268 = xi_267 + xia_8_collide;
                  const double xia_14_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0];
                  const double xia_15_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0];
                  const double xia_16_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0];
                  const double xia_17_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0];
                  const double xi_177 = xia_17_collide + xia_3_collide;
                  const double xia_18_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + ctr_0];
                  const double xi_17 = xia_18_collide*0.16666666666666666;
                  const double xi_18 = omega_odd_a*xi_17;
                  const double xi_32 = xia_18_collide*0.083333333333333329;
                  const double xi_34 = xi_33*xia_18_collide;
                  const double xi_35 = xi_32 - xi_34;
                  const double xi_51 = -xi_32 + xi_34;
                  const double xia_19_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0];
                  const double xia_20_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0];
                  const double xi_170 = xia_15_collide + xia_20_collide;
                  const double xia_21_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0];
                  const double xi_178 = xi_177 + xia_21_collide + xia_9_collide;
                  const double xi_248 = -xia_21_collide;
                  const double xi_249 = xi_248 + xia_9_collide;
                  const double xia_22_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0];
                  const double xia_23_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0];
                  const double xi_171 = xi_170 + xia_1_collide + xia_23_collide;
                  const double xi_254 = xi_253 + xia_23_collide;
                  const double xia_24_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0];
                  const double xia_25_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0];
                  const double xia_26_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3 + ctr_0];
                  const double xi_26 = xia_26_collide*0.16666666666666666;
                  const double xi_28 = omega_odd_a*xi_26;
                  const double xi_55 = xia_26_collide*0.083333333333333329;
                  const double xi_56 = xi_33*xia_26_collide;
                  const double xi_57 = xi_55 - xi_56;
                  const double xi_70 = -xi_55 + xi_56;
                  const double xi_6 = xi_5*xia_5_collide;
                  const double xi_19 = xi_5*xia_18_collide;
                  const double xi_27 = xi_5*xia_26_collide;
                  const double rho_a_collide = xia_7_collide;
                  const double xi_167 = rho_a_collide*-0.1111111111111111;
                  const double xi_180 = rho_a_collide*-0.33333333333333331;
                  const double xi_181 = xi_175 + xi_180;
                  const double xi_246 = rho_a_collide*0.33333333333333331;
                  const double u_0_a_collide = xia_24_collide;
                  const double xi_0 = u_0_a_collide*xia_18_collide;
                  const double xi_13 = xi_0*0.083333333333333329;
                  const double xi_14 = omega_even_a*xi_13 + xi_0*-0.16666666666666666;
                  const double xi_20 = xi_0*0.33333333333333331;
                  const double xi_42 = u_0_a_collide*xi_41;
                  const double xi_47 = u_0_a_collide*xi_46;
                  const double xi_58 = -xi_13;
                  const double xi_61 = omega_even_a*xi_0*0.041666666666666664;
                  const double xi_73 = u_0_a_collide*xia_26_collide;
                  const double xi_74 = xi_73*0.25;
                  const double xi_77 = xi_45*xi_73;
                  const double xi_166 = (u_0_a_collide*u_0_a_collide);
                  const double xi_173 = rho_a_collide*xi_166*-0.33333333333333331;
                  const double xi_182 = omega_shear_a*(rho_a_collide*xi_166 - xi_178 - xi_181 - xia_12_collide - xia_16_collide);
                  const double xi_264 = u_0_a_collide*xi_246;
                  const double xi_265 = xi_249 + xi_264 - xia_17_collide + xia_3_collide;
                  const double xi_266 = xi_251*xi_265;
                  const double xi_269 = xi_264 + xi_268 - xia_11_collide + xia_4_collide;
                  const double xi_270 = xi_251*xi_269;
                  const double xi_284 = xi_265*xi_283;
                  const double xi_285 = -xi_284;
                  const double xi_327 = xi_269*xi_283;
                  const double xi_328 = -xi_327;
                  const double u_1_a_collide = xia_10_collide;
                  const double xi_1 = u_1_a_collide*xia_5_collide;
                  const double xi_8 = xi_1*0.33333333333333331;
                  const double xi_21 = xi_1*0.16666666666666666;
                  const double xi_22 = xi_1*0.083333333333333329;
                  const double xi_23 = omega_even_a*xi_22;
                  const double xi_24 = -xi_21 + xi_23;
                  const double xi_30 = xi_14 + xi_24;
                  const double xi_43 = u_1_a_collide*0.25;
                  const double xi_44 = xi_43*xia_18_collide;
                  const double xi_48 = u_1_a_collide*xi_45;
                  const double xi_49 = xi_48*xia_18_collide;
                  const double xi_50 = xi_42 + xi_44 - xi_47 - xi_49;
                  const double xi_52 = -xi_42 - xi_44 + xi_47 + xi_49;
                  const double xi_59 = -xi_23;
                  const double xi_63 = xi_43*xia_26_collide;
                  const double xi_65 = xi_48*xia_26_collide;
                  const double xi_71 = omega_even_a*u_1_a_collide*xia_5_collide*-0.041666666666666664;
                  const double xi_164 = (u_1_a_collide*u_1_a_collide);
                  const double xi_165 = rho_a_collide*xi_164*-0.33333333333333331;
                  const double xi_183 = omega_shear_a*(rho_a_collide*xi_164 - xi_171 - xi_178 - xi_180 - xia_2_collide - xia_6_collide);
                  const double xi_247 = u_1_a_collide*xi_246;
                  const double xi_250 = xi_247 + xi_249 + xia_17_collide - xia_3_collide;
                  const double xi_252 = xi_250*xi_251;
                  const double xi_255 = xi_247 + xi_254 - xia_15_collide + xia_20_collide;
                  const double xi_256 = xi_251*xi_255;
                  const double xi_286 = rho_a_collide*u_1_a_collide;
                  const double xi_288 = xi_287*(u_0_a_collide*xi_286 + xi_177 + xi_248 - xia_9_collide);
                  const double xi_289 = -xi_288;
                  const double xi_300 = xi_250*xi_283;
                  const double xi_311 = xi_255*xi_283;
                  const double xi_321 = -xi_311;
                  const double u_2_a_collide = xia_14_collide;
                  const double xi_2 = u_2_a_collide*xia_26_collide;
                  const double xi_9 = xi_2*0.16666666666666666;
                  const double xi_10 = xi_2*0.083333333333333329;
                  const double xi_11 = omega_even_a*xi_10;
                  const double xi_12 = xi_11 - xi_9;
                  const double xi_15 = xi_12 + xi_14;
                  const double xi_16 = omega_even_a*xi_8 + omega_shear_a*xi_1*-0.5 + xi_15 + xi_8;
                  const double xi_25 = omega_even_a*xi_20 + omega_shear_a*xi_0*-0.5 + xi_12 + xi_20 + xi_24;
                  const double xi_29 = xi_2*0.33333333333333331;
                  const double xi_31 = omega_even_a*xi_29 + omega_shear_a*xi_2*-0.5 + xi_29 + xi_30;
                  const double xi_36 = omega_even_a*u_2_a_collide*xia_26_collide*-0.041666666666666664;
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
                  const double xi_75 = u_2_a_collide*xia_18_collide;
                  const double xi_76 = xi_75*0.25;
                  const double xi_78 = xi_45*xi_75;
                  const double xi_79 = xi_74 + xi_76 - xi_77 - xi_78;
                  const double xi_80 = -xi_74 - xi_76 + xi_77 + xi_78;
                  const double xi_81 = xi_15 + xi_22 + xi_51 + xi_71;
                  const double xi_168 = (u_2_a_collide*u_2_a_collide);
                  const double xi_169 = rho_a_collide*xi_168*-0.33333333333333331 + xi_167;
                  const double xi_172 = omega_even_a*(rho_a_collide*xi_166*-0.16666666666666666 - xi_165 - xi_169 - xi_171);
                  const double xi_176 = omega_even_a*(rho_a_collide*xi_164*-0.16666666666666666 - xi_169 - xi_173 - xi_175);
                  const double xi_179 = omega_even_a*(rho_a_collide*xi_168*-0.16666666666666666 - xi_165 - xi_167 - xi_173 - xi_178);
                  const double xi_184 = omega_shear_a*(rho_a_collide*xi_168 - xi_171 - xi_181 - xia_19_collide - xia_22_collide);
                  const double xi_260 = xi_172*-0.5;
                  const double xi_261 = xi_179*-0.5;
                  const double xi_262 = xi_183*0.5 + xi_260 + xi_261;
                  const double xi_272 = xi_176*-0.5;
                  const double xi_273 = xi_182*0.5 + xi_261 + xi_272;
                  const double xi_275 = u_2_a_collide*xi_246;
                  const double xi_276 = xi_268 + xi_275 + xia_11_collide - xia_4_collide;
                  const double xi_277 = xi_251*xi_276;
                  const double xi_278 = xi_254 + xi_275 + xia_15_collide - xia_20_collide;
                  const double xi_279 = xi_251*xi_278;
                  const double xi_281 = xi_184*0.5 + xi_260 + xi_272;
                  const double xi_299 = xi_179*0.25;
                  const double xi_301 = xi_299 + xi_300;
                  const double xi_309 = xi_299 - xi_300;
                  const double xi_312 = xi_287*(u_2_a_collide*xi_286 + xi_170 + xi_253 - xia_23_collide);
                  const double xi_318 = xi_172*0.25;
                  const double xi_319 = xi_278*xi_283;
                  const double xi_320 = xi_318 + xi_319;
                  const double xi_322 = -xi_312;
                  const double xi_329 = xi_287*(rho_a_collide*u_0_a_collide*u_2_a_collide + xi_174 + xi_267 - xia_8_collide);
                  const double xi_330 = -xi_329;
                  const double xi_334 = xi_176*0.25;
                  const double xi_335 = xi_276*xi_283;
                  const double xi_336 = xi_334 + xi_335;
                  const double xi_341 = xi_318 - xi_319;
                  const double xi_344 = xi_334 - xi_335;
                  const double forceTerm_0_a_collide = omega_shear_a*u_0_a_collide*xia_18_collide + omega_shear_a*u_1_a_collide*xia_5_collide + omega_shear_a*u_2_a_collide*xia_26_collide - xi_0*xi_3 - xi_0 - xi_1*xi_3 - xi_1 - xi_2*xi_3 - xi_2;
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
                  const double xib_1_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3 + ctr_0];
                  const double xi_86 = xib_1_collide*0.16666666666666666;
                  const double xi_89 = omega_odd_b*xi_86;
                  const double xi_119 = xib_1_collide*0.083333333333333329;
                  const double xi_120 = xi_115*xib_1_collide;
                  const double xi_121 = -xi_119 + xi_120;
                  const double xi_123 = xib_1_collide*0.25;
                  const double xi_128 = xi_127*xib_1_collide;
                  const double xi_135 = xi_119 - xi_120;
                  const double xib_2_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0];
                  const double xi_374 = -xib_2_collide;
                  const double xib_3_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0];
                  const double xib_4_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3 + ctr_0];
                  const double xi_108 = xib_4_collide*0.16666666666666666;
                  const double xi_110 = omega_odd_b*xi_108;
                  const double xi_137 = xib_4_collide*0.083333333333333329;
                  const double xi_138 = xi_115*xib_4_collide;
                  const double xi_139 = xi_137 - xi_138;
                  const double xi_152 = -xi_137 + xi_138;
                  const double xib_5_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0];
                  const double xi_382 = -xib_5_collide;
                  const double xib_6_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0];
                  const double xib_7_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0];
                  const double xib_8_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0];
                  const double xib_9_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0];
                  const double xi_383 = xi_382 + xib_9_collide;
                  const double xib_10_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0];
                  const double xib_11_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0];
                  const double xib_12_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0];
                  const double xib_13_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0];
                  const double xib_14_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0];
                  const double xi_352 = xib_14_collide + xib_8_collide;
                  const double xib_15_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0];
                  const double xib_16_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0];
                  const double xib_17_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0];
                  const double xi_375 = xi_374 + xib_17_collide;
                  const double xib_18_collide = _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const double xib_19_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0];
                  const double xi_356 = xib_16_collide + xib_19_collide;
                  const double xi_357 = xi_356 + xib_5_collide + xib_9_collide;
                  const double xib_20_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + ctr_0];
                  const double xi_99 = xib_20_collide*0.16666666666666666;
                  const double xi_100 = omega_odd_b*xi_99;
                  const double xi_114 = xib_20_collide*0.083333333333333329;
                  const double xi_116 = xi_115*xib_20_collide;
                  const double xi_117 = xi_114 - xi_116;
                  const double xi_133 = -xi_114 + xi_116;
                  const double xib_21_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0];
                  const double xib_22_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0];
                  const double xib_23_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0];
                  const double xib_24_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0];
                  const double xib_25_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0];
                  const double xi_353 = xi_352 + xib_22_collide + xib_25_collide;
                  const double xi_369 = -xib_25_collide;
                  const double xi_370 = xi_369 + xib_22_collide;
                  const double xib_26_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0];
                  const double xi_359 = xib_21_collide + xib_26_collide;
                  const double xi_360 = xi_359 + xib_17_collide + xib_2_collide;
                  const double xi_88 = xi_87*xib_1_collide;
                  const double xi_101 = xi_87*xib_20_collide;
                  const double xi_109 = xi_87*xib_4_collide;
                  const double rho_b_collide = xib_18_collide;
                  const double xi_349 = rho_b_collide*-0.1111111111111111;
                  const double xi_362 = rho_b_collide*-0.33333333333333331;
                  const double xi_363 = xi_360 + xi_362;
                  const double xi_367 = rho_b_collide*0.33333333333333331;
                  const double u_0_b_collide = xib_23_collide;
                  const double xi_82 = u_0_b_collide*xib_20_collide;
                  const double xi_95 = xi_82*0.083333333333333329;
                  const double xi_96 = omega_even_b*xi_95 + xi_82*-0.16666666666666666;
                  const double xi_102 = xi_82*0.33333333333333331;
                  const double xi_124 = u_0_b_collide*xi_123;
                  const double xi_129 = u_0_b_collide*xi_128;
                  const double xi_140 = -xi_95;
                  const double xi_143 = omega_even_b*xi_82*0.041666666666666664;
                  const double xi_155 = u_0_b_collide*xib_4_collide;
                  const double xi_156 = xi_155*0.25;
                  const double xi_159 = xi_127*xi_155;
                  const double xi_348 = (u_0_b_collide*u_0_b_collide);
                  const double xi_355 = rho_b_collide*xi_348*-0.33333333333333331;
                  const double xi_364 = omega_shear_b*(rho_b_collide*xi_348 - xi_357 - xi_363 - xib_12_collide - xib_6_collide);
                  const double xi_381 = u_0_b_collide*xi_367;
                  const double xi_384 = xi_381 + xi_383 + xib_16_collide - xib_19_collide;
                  const double xi_385 = xi_372*xi_384;
                  const double xi_386 = xi_375 + xi_381 - xib_21_collide + xib_26_collide;
                  const double xi_387 = xi_372*xi_386;
                  const double xi_397 = xi_386*xi_396;
                  const double xi_398 = -xi_397;
                  const double xi_414 = xi_384*xi_396;
                  const double xi_415 = -xi_414;
                  const double u_1_b_collide = xib_11_collide;
                  const double xi_83 = u_1_b_collide*xib_1_collide;
                  const double xi_90 = xi_83*0.33333333333333331;
                  const double xi_103 = xi_83*0.16666666666666666;
                  const double xi_104 = xi_83*0.083333333333333329;
                  const double xi_105 = omega_even_b*xi_104;
                  const double xi_106 = -xi_103 + xi_105;
                  const double xi_112 = xi_106 + xi_96;
                  const double xi_125 = u_1_b_collide*0.25;
                  const double xi_126 = xi_125*xib_20_collide;
                  const double xi_130 = u_1_b_collide*xi_127;
                  const double xi_131 = xi_130*xib_20_collide;
                  const double xi_132 = xi_124 + xi_126 - xi_129 - xi_131;
                  const double xi_134 = -xi_124 - xi_126 + xi_129 + xi_131;
                  const double xi_141 = -xi_105;
                  const double xi_145 = xi_125*xib_4_collide;
                  const double xi_147 = xi_130*xib_4_collide;
                  const double xi_153 = omega_even_b*u_1_b_collide*xib_1_collide*-0.041666666666666664;
                  const double xi_346 = (u_1_b_collide*u_1_b_collide);
                  const double xi_347 = rho_b_collide*xi_346*-0.33333333333333331;
                  const double xi_365 = omega_shear_b*(rho_b_collide*xi_346 - xi_353 - xi_363 - xib_10_collide - xib_3_collide);
                  const double xi_368 = u_1_b_collide*xi_367;
                  const double xi_371 = xi_368 + xi_370 + xib_14_collide - xib_8_collide;
                  const double xi_373 = xi_371*xi_372;
                  const double xi_376 = xi_368 + xi_375 + xib_21_collide - xib_26_collide;
                  const double xi_377 = xi_372*xi_376;
                  const double xi_399 = rho_b_collide*u_1_b_collide;
                  const double xi_401 = xi_400*(u_0_b_collide*xi_399 + xi_359 + xi_374 - xib_17_collide);
                  const double xi_402 = -xi_401;
                  const double xi_404 = xi_376*xi_396;
                  const double xi_407 = xi_371*xi_396;
                  const double xi_412 = -xi_407;
                  const double u_2_b_collide = xib_15_collide;
                  const double xi_84 = u_2_b_collide*xib_4_collide;
                  const double xi_91 = xi_84*0.16666666666666666;
                  const double xi_92 = xi_84*0.083333333333333329;
                  const double xi_93 = omega_even_b*xi_92;
                  const double xi_94 = -xi_91 + xi_93;
                  const double xi_97 = xi_94 + xi_96;
                  const double xi_98 = omega_even_b*xi_90 + omega_shear_b*xi_83*-0.5 + xi_90 + xi_97;
                  const double xi_107 = omega_even_b*xi_102 + omega_shear_b*xi_82*-0.5 + xi_102 + xi_106 + xi_94;
                  const double xi_111 = xi_84*0.33333333333333331;
                  const double xi_113 = omega_even_b*xi_111 + omega_shear_b*xi_84*-0.5 + xi_111 + xi_112;
                  const double xi_118 = omega_even_b*u_2_b_collide*xib_4_collide*-0.041666666666666664;
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
                  const double xi_157 = u_2_b_collide*xib_20_collide;
                  const double xi_158 = xi_157*0.25;
                  const double xi_160 = xi_127*xi_157;
                  const double xi_161 = xi_156 + xi_158 - xi_159 - xi_160;
                  const double xi_162 = -xi_156 - xi_158 + xi_159 + xi_160;
                  const double xi_163 = xi_104 + xi_133 + xi_153 + xi_97;
                  const double xi_350 = (u_2_b_collide*u_2_b_collide);
                  const double xi_351 = rho_b_collide*xi_350*-0.33333333333333331 + xi_349;
                  const double xi_354 = omega_even_b*(rho_b_collide*xi_348*-0.16666666666666666 - xi_347 - xi_351 - xi_353);
                  const double xi_358 = omega_even_b*(rho_b_collide*xi_346*-0.16666666666666666 - xi_351 - xi_355 - xi_357);
                  const double xi_361 = omega_even_b*(rho_b_collide*xi_350*-0.16666666666666666 - xi_347 - xi_349 - xi_355 - xi_360);
                  const double xi_366 = omega_shear_b*(rho_b_collide*xi_350 - xi_353 - xi_357 - xi_362 - xib_13_collide - xib_24_collide);
                  const double xi_378 = xi_354*-0.5;
                  const double xi_379 = xi_361*-0.5;
                  const double xi_380 = xi_365*0.5 + xi_378 + xi_379;
                  const double xi_388 = xi_358*-0.5;
                  const double xi_389 = xi_364*0.5 + xi_379 + xi_388;
                  const double xi_390 = u_2_b_collide*xi_367;
                  const double xi_391 = xi_383 + xi_390 - xib_16_collide + xib_19_collide;
                  const double xi_392 = xi_372*xi_391;
                  const double xi_393 = xi_370 + xi_390 - xib_14_collide + xib_8_collide;
                  const double xi_394 = xi_372*xi_393;
                  const double xi_395 = xi_366*0.5 + xi_378 + xi_388;
                  const double xi_403 = xi_361*0.25;
                  const double xi_405 = xi_403 + xi_404;
                  const double xi_406 = xi_403 - xi_404;
                  const double xi_408 = xi_400*(u_2_b_collide*xi_399 + xi_352 + xi_369 - xib_22_collide);
                  const double xi_409 = xi_354*0.25;
                  const double xi_410 = xi_393*xi_396;
                  const double xi_411 = xi_409 + xi_410;
                  const double xi_413 = -xi_408;
                  const double xi_416 = xi_400*(rho_b_collide*u_0_b_collide*u_2_b_collide + xi_356 + xi_382 - xib_9_collide);
                  const double xi_417 = -xi_416;
                  const double xi_418 = xi_358*0.25;
                  const double xi_419 = xi_391*xi_396;
                  const double xi_420 = xi_418 + xi_419;
                  const double xi_421 = xi_409 - xi_410;
                  const double xi_422 = xi_418 - xi_419;
                  const double forceTerm_0_b_collide = omega_shear_b*u_0_b_collide*xib_20_collide + omega_shear_b*u_1_b_collide*xib_1_collide + omega_shear_b*u_2_b_collide*xib_4_collide - xi_82*xi_85 - xi_82 - xi_83*xi_85 - xi_83 - xi_84*xi_85 - xi_84;
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
                  const double tmp_a0 = forceTerm_0_a_collide + xi_172 + xi_176 + xi_179 - xi_182 - xi_183 - xi_184 + xi_245 + xia_25_collide;
                  const double tmp_a1 = forceTerm_1_a_collide - xi_252 - xi_256 + xi_259 + xi_262 + xia_6_collide;
                  const double tmp_a2 = forceTerm_2_a_collide + xi_252 + xi_256 + xi_262 + xi_263 + xia_2_collide;
                  const double tmp_a3 = forceTerm_3_a_collide + xi_266 + xi_270 + xi_271 + xi_273 + xia_16_collide;
                  const double tmp_a4 = forceTerm_4_a_collide - xi_266 - xi_270 + xi_273 + xi_274 + xia_12_collide;
                  const double tmp_a5 = forceTerm_5_a_collide - xi_277 - xi_279 + xi_280 + xi_281 + xia_19_collide;
                  const double tmp_a6 = forceTerm_6_a_collide + xi_277 + xi_279 + xi_281 + xi_282 + xia_22_collide;
                  const double tmp_a7 = forceTerm_7_a_collide + xi_285 + xi_289 + xi_298 + xi_301 + xia_3_collide;
                  const double tmp_a8 = forceTerm_8_a_collide + xi_284 + xi_288 + xi_301 + xi_307 + xia_21_collide;
                  const double tmp_a9 = forceTerm_9_a_collide + xi_285 + xi_288 + xi_308 + xi_309 + xia_9_collide;
                  const double tmp_a10 = forceTerm_10_a_collide + xi_284 + xi_289 + xi_309 + xi_310 + xia_17_collide;
                  const double tmp_a11 = forceTerm_11_a_collide + xi_311 + xi_312 + xi_317 + xi_320 + xia_1_collide;
                  const double tmp_a12 = forceTerm_12_a_collide + xi_320 + xi_321 + xi_322 + xi_326 + xia_20_collide;
                  const double tmp_a13 = forceTerm_13_a_collide + xi_328 + xi_330 + xi_333 + xi_336 + xia_4_collide;
                  const double tmp_a14 = forceTerm_14_a_collide + xi_327 + xi_329 + xi_336 + xi_339 + xia_13_collide;
                  const double tmp_a15 = forceTerm_15_a_collide + xi_311 + xi_322 + xi_340 + xi_341 + xia_15_collide;
                  const double tmp_a16 = forceTerm_16_a_collide + xi_312 + xi_321 + xi_341 + xi_342 + xia_23_collide;
                  const double tmp_a17 = forceTerm_17_a_collide + xi_328 + xi_329 + xi_343 + xi_344 + xia_8_collide;
                  const double tmp_a18 = forceTerm_18_a_collide + xi_327 + xi_330 + xi_344 + xi_345 + xia_11_collide;
                  const double tmp_b0 = forceTerm_0_b_collide + xi_245 + xi_354 + xi_358 + xi_361 - xi_364 - xi_365 - xi_366 + xib_7_collide;
                  const double tmp_b1 = forceTerm_1_b_collide + xi_259 - xi_373 - xi_377 + xi_380 + xib_10_collide;
                  const double tmp_b2 = forceTerm_2_b_collide + xi_263 + xi_373 + xi_377 + xi_380 + xib_3_collide;
                  const double tmp_b3 = forceTerm_3_b_collide + xi_271 + xi_385 + xi_387 + xi_389 + xib_12_collide;
                  const double tmp_b4 = forceTerm_4_b_collide + xi_274 - xi_385 - xi_387 + xi_389 + xib_6_collide;
                  const double tmp_b5 = forceTerm_5_b_collide + xi_280 - xi_392 - xi_394 + xi_395 + xib_13_collide;
                  const double tmp_b6 = forceTerm_6_b_collide + xi_282 + xi_392 + xi_394 + xi_395 + xib_24_collide;
                  const double tmp_b7 = forceTerm_7_b_collide + xi_298 + xi_398 + xi_402 + xi_405 + xib_26_collide;
                  const double tmp_b8 = forceTerm_8_b_collide + xi_307 + xi_397 + xi_401 + xi_405 + xib_2_collide;
                  const double tmp_b9 = forceTerm_9_b_collide + xi_308 + xi_398 + xi_401 + xi_406 + xib_17_collide;
                  const double tmp_b10 = forceTerm_10_b_collide + xi_310 + xi_397 + xi_402 + xi_406 + xib_21_collide;
                  const double tmp_b11 = forceTerm_11_b_collide + xi_317 + xi_407 + xi_408 + xi_411 + xib_25_collide;
                  const double tmp_b12 = forceTerm_12_b_collide + xi_326 + xi_411 + xi_412 + xi_413 + xib_14_collide;
                  const double tmp_b13 = forceTerm_13_b_collide + xi_333 + xi_415 + xi_417 + xi_420 + xib_16_collide;
                  const double tmp_b14 = forceTerm_14_b_collide + xi_339 + xi_414 + xi_416 + xi_420 + xib_5_collide;
                  const double tmp_b15 = forceTerm_15_b_collide + xi_340 + xi_407 + xi_413 + xi_421 + xib_8_collide;
                  const double tmp_b16 = forceTerm_16_b_collide + xi_342 + xi_408 + xi_412 + xi_421 + xib_22_collide;
                  const double tmp_b17 = forceTerm_17_b_collide + xi_343 + xi_415 + xi_416 + xi_422 + xib_9_collide;
                  const double tmp_b18 = forceTerm_18_b_collide + xi_345 + xi_414 + xi_417 + xi_422 + xib_19_collide;
                  const double xirecolor_0 = tmp_a0 + tmp_b0;
                  const double xirecolor_1 = _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const double xirecolor_2 = ((1.0) / (xirecolor_1));
                  const double xi_444 = xirecolor_2*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const double xirecolor_3 = xirecolor_2*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0];
                  const double xirecolor_4 = tmp_a1 + tmp_b1;
                  const double xirecolor_5 = xi_195;
                  const double xirecolor_6 = xi_193;
                  const double xirecolor_7 = xi_213;
                  const double xi_423 = -xirecolor_7;
                  const double xirecolor_8 = xi_423;
                  const double xirecolor_9 = xi_200;
                  const double xirecolor_10 = -xirecolor_9;
                  const double xirecolor_11 = xi_202;
                  const double xi_424 = xirecolor_10 + xirecolor_11;
                  const double xirecolor_12 = xi_215;
                  const double xi_425 = xirecolor_12 + xirecolor_8;
                  const double xirecolor_13 = xi_424 + xi_425;
                  const double xirecolor_14 = xi_217;
                  const double xirecolor_15 = xi_218;
                  const double xi_426 = xirecolor_14 - xirecolor_15;
                  const double xirecolor_16 = -xi_426;
                  const double xirecolor_17 = xi_187;
                  const double xirecolor_18 = xi_188;
                  const double xirecolor_19 = xi_189;
                  const double xi_427 = xirecolor_18 + xirecolor_19;
                  const double xirecolor_20 = xi_190;
                  const double xirecolor_21 = xi_427 - xirecolor_17 + xirecolor_20;
                  const double xirecolor_22 = xi_186;
                  const double xirecolor_23 = xi_185;
                  const double xi_428 = xirecolor_22 - xirecolor_23;
                  const double xirecolor_24 = -xi_428;
                  const double xirecolor_25 = xi_221;
                  const double xirecolor_26 = -xirecolor_25;
                  const double xirecolor_27 = xi_207;
                  const double xirecolor_28 = xi_206;
                  const double xirecolor_29 = xi_223;
                  const double xi_429 = xirecolor_26 + xirecolor_29;
                  const double xirecolor_30 = xi_429 - xirecolor_27 + xirecolor_28;
                  const double xirecolor_31 = xirecolor_13 + xirecolor_16 + xirecolor_21 + xirecolor_24 + xirecolor_30 - xirecolor_5 + xirecolor_6;
                  const double xirecolor_32 = xi_192;
                  const double xirecolor_33 = xi_194;
                  const double xi_430 = xirecolor_32 + xirecolor_33;
                  const double xirecolor_34 = xi_430 + xirecolor_5 - xirecolor_6;
                  const double xirecolor_35 = xi_197;
                  const double xirecolor_36 = xi_198;
                  const double xi_431 = xirecolor_35 - xirecolor_36;
                  const double xirecolor_37 = -xi_431;
                  const double xirecolor_38 = xi_424;
                  const double xirecolor_39 = xi_204;
                  const double xirecolor_40 = xi_205;
                  const double xi_432 = xirecolor_39 + xirecolor_40;
                  const double xirecolor_41 = xirecolor_27 - xirecolor_28;
                  const double xirecolor_42 = xi_432 + xirecolor_41;
                  const double xirecolor_43 = xirecolor_38 + xirecolor_42;
                  const double xirecolor_44 = xi_428 + xirecolor_21 + xirecolor_34 + xirecolor_37 + xirecolor_43;
                  const double xirecolor_45 = xi_426;
                  const double xirecolor_46 = xi_425;
                  const double xirecolor_47 = xi_229;
                  const double xirecolor_48 = -xirecolor_47;
                  const double xirecolor_49 = xi_231;
                  const double xi_433 = xirecolor_48 + xirecolor_49;
                  const double xirecolor_50 = xi_431 + xi_433;
                  const double xirecolor_51 = xi_427 + xirecolor_17 - xirecolor_20 + xirecolor_24 + xirecolor_34 + xirecolor_45 + xirecolor_46 + xirecolor_50;
                  const double xirecolor_52 = pow((xirecolor_31*xirecolor_31) + (xirecolor_44*xirecolor_44) + (xirecolor_51*xirecolor_51), 0.5);
                  const double xirecolor_53 = ((1.0) / (xirecolor_52));
                  const double xi_434 = xirecolor_31*xirecolor_53;
                  const double xi_435 = xirecolor_51*xirecolor_53;
                  const double xi_436 = xirecolor_44*xirecolor_53;
                  const bool xirecolor_54 = xirecolor_52 > 0.0;
                  const double xirecolor_55 = beta*((1.0) / ((xirecolor_1*xirecolor_1)))*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const double xirecolor_56 = xirecolor_55*(0.055555555555555552*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + 0.055555555555555552*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const double xirecolor_57 = xirecolor_56*((xirecolor_54) ? (xi_434): (0.0));
                  const double xirecolor_58 = tmp_a2 + tmp_b2;
                  const double xirecolor_59 = xirecolor_56*((xirecolor_54) ? (-xi_434): (0.0));
                  const double xirecolor_60 = tmp_a3 + tmp_b3;
                  const double xirecolor_61 = xirecolor_56*((xirecolor_54) ? (-xi_435): (0.0));
                  const double xirecolor_62 = tmp_a4 + tmp_b4;
                  const double xirecolor_63 = xirecolor_56*((xirecolor_54) ? (xi_435): (0.0));
                  const double xirecolor_64 = tmp_a5 + tmp_b5;
                  const double xirecolor_65 = xirecolor_56*((xirecolor_54) ? (-xi_436): (0.0));
                  const double xirecolor_66 = tmp_a6 + tmp_b6;
                  const double xirecolor_67 = xirecolor_56*((xirecolor_54) ? (xi_436): (0.0));
                  const double xirecolor_68 = tmp_a7 + tmp_b7;
                  const double xirecolor_69 = xi_290;
                  const double xirecolor_70 = xi_291;
                  const double xirecolor_71 = xi_430;
                  const double xirecolor_72 = xirecolor_50 + xirecolor_71;
                  const double xirecolor_73 = xi_295;
                  const double xirecolor_74 = xirecolor_25 - xirecolor_29 + xirecolor_73;
                  const double xirecolor_75 = xi_292 - xirecolor_11 + xirecolor_41 - xirecolor_69 + xirecolor_70 + xirecolor_72 + xirecolor_74 + xirecolor_9;
                  const double xirecolor_76 = xirecolor_53*0.70710678118654757;
                  const double xi_437 = xirecolor_75*xirecolor_76;
                  const double xirecolor_77 = xirecolor_55*(0.027777777777777776*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + 0.027777777777777776*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const double xirecolor_78 = xirecolor_77*((xirecolor_54) ? (-xi_437): (0.0));
                  const double xirecolor_79 = tmp_a8 + tmp_b8;
                  const double xirecolor_80 = xi_302;
                  const double xirecolor_81 = xi_303;
                  const double xi_438 = -xirecolor_80 + xirecolor_81;
                  const double xirecolor_82 = xi_304;
                  const double xirecolor_83 = xi_305 + xi_438 + xirecolor_30 + xirecolor_38 + xirecolor_72 + xirecolor_82;
                  const double xi_439 = xirecolor_76*xirecolor_83;
                  const double xirecolor_84 = xirecolor_77*((xirecolor_54) ? (xi_439): (0.0));
                  const double xirecolor_85 = tmp_a9 + tmp_b9;
                  const double xirecolor_86 = xirecolor_77*((xirecolor_54) ? (-xi_439): (0.0));
                  const double xirecolor_87 = tmp_a10 + tmp_b10;
                  const double xirecolor_88 = xirecolor_77*((xirecolor_54) ? (xi_437): (0.0));
                  const double xirecolor_89 = tmp_a11 + tmp_b11;
                  const double xirecolor_90 = xi_432 + xirecolor_37 + xirecolor_71;
                  const double xirecolor_91 = -xi_423 - xi_438 - xirecolor_12;
                  const double xirecolor_92 = xi_313 + xirecolor_45 + xirecolor_74 + xirecolor_90 + xirecolor_91;
                  const double xi_440 = xirecolor_76*xirecolor_92;
                  const double xirecolor_93 = xirecolor_77*((xirecolor_54) ? (-xi_440): (0.0));
                  const double xirecolor_94 = tmp_a12 + tmp_b12;
                  const double xirecolor_95 = xirecolor_16 + xirecolor_69 - xirecolor_70;
                  const double xirecolor_96 = xi_323 + xi_429 + xirecolor_46 + xirecolor_82 + xirecolor_90 + xirecolor_95;
                  const double xi_441 = xirecolor_76*xirecolor_96;
                  const double xirecolor_97 = xirecolor_77*((xirecolor_54) ? (-xi_441): (0.0));
                  const double xirecolor_98 = tmp_a13 + tmp_b13;
                  const double xirecolor_99 = xi_331 + xi_433 + xirecolor_13 + xirecolor_42 + xirecolor_45 + xirecolor_73 + xirecolor_82;
                  const double xi_442 = xirecolor_76*xirecolor_99;
                  const double xirecolor_100 = xirecolor_77*((xirecolor_54) ? (-xi_442): (0.0));
                  const double xirecolor_101 = tmp_a14 + tmp_b14;
                  const double xirecolor_102 = xi_337 + xirecolor_43 + xirecolor_47 - xirecolor_49 + xirecolor_91 + xirecolor_95;
                  const double xi_443 = xirecolor_102*xirecolor_76;
                  const double xirecolor_103 = xirecolor_77*((xirecolor_54) ? (-xi_443): (0.0));
                  const double xirecolor_104 = tmp_a15 + tmp_b15;
                  const double xirecolor_105 = xirecolor_77*((xirecolor_54) ? (xi_441): (0.0));
                  const double xirecolor_106 = tmp_a16 + tmp_b16;
                  const double xirecolor_107 = xirecolor_77*((xirecolor_54) ? (xi_440): (0.0));
                  const double xirecolor_108 = tmp_a17 + tmp_b17;
                  const double xirecolor_109 = xirecolor_77*((xirecolor_54) ? (xi_443): (0.0));
                  const double xirecolor_110 = tmp_a18 + tmp_b18;
                  const double xirecolor_111 = xirecolor_77*((xirecolor_54) ? (xi_442): (0.0));
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0] = xirecolor_0*xirecolor_3;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_4 + xirecolor_57;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_58 + xirecolor_59;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_60 + xirecolor_61;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_62 + xirecolor_63;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_64 + xirecolor_65;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_66 + xirecolor_67;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_68 + xirecolor_78;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_79 + xirecolor_84;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_85 + xirecolor_86;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_87 + xirecolor_88;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_89 + xirecolor_93;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_94 + xirecolor_97;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0] = xirecolor_100 + xirecolor_3*xirecolor_98;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0] = xirecolor_101*xirecolor_3 + xirecolor_103;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0] = xirecolor_104*xirecolor_3 + xirecolor_105;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0] = xirecolor_106*xirecolor_3 + xirecolor_107;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0] = xirecolor_108*xirecolor_3 + xirecolor_109;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0] = xirecolor_110*xirecolor_3 + xirecolor_111;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0] = xi_444*xirecolor_0;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_4 - xirecolor_57;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_58 - xirecolor_59;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_60 - xirecolor_61;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_62 - xirecolor_63;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_64 - xirecolor_65;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_66 - xirecolor_67;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_68 - xirecolor_78;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_79 - xirecolor_84;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_85 - xirecolor_86;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_87 - xirecolor_88;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_89 - xirecolor_93;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_94 - xirecolor_97;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_98 - xirecolor_100;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_101 - xirecolor_103;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_104 - xirecolor_105;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_106 - xirecolor_107;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_108 - xirecolor_109;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0] = xi_444*xirecolor_110 - xirecolor_111;
               }
            }
         }
      }
   }
}
}


void ColorGradientCollideSweepDoublePrecisionAVX::run(IBlock * block)
{
   
    auto rho_a = block->getData< field::GhostLayerField<double, 1> >(rho_aID);
    auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
    auto rho_b = block->getData< field::GhostLayerField<double, 1> >(rho_bID);
    auto force_b = block->getData< field::GhostLayerField<double, 3> >(force_bID);
    auto pdfs_b = block->getData< field::GhostLayerField<double, 19> >(pdfs_bID);
    auto pdfs_a = block->getData< field::GhostLayerField<double, 19> >(pdfs_aID);
    auto phasefield = block->getData< field::GhostLayerField<double, 1> >(phasefieldID);
    auto force_a = block->getData< field::GhostLayerField<double, 3> >(force_aID);

    auto & beta = this->beta_;
    auto & omega_odd_b = this->omega_odd_b_;
    auto & sigma = this->sigma_;
    auto & omega_odd_a = this->omega_odd_a_;
    auto & omega_even_b = this->omega_even_b_;
    auto & omega_shear_b = this->omega_shear_b_;
    auto & omega_even_a = this->omega_even_a_;
    auto & omega_shear_a = this->omega_shear_a_;
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_a->nrOfGhostLayers()))
    double * RESTRICT const _data_force_a = force_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT const _data_phasefield = phasefield->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) phasefield->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_a = rho_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_b = rho_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) velocity->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(force_a->xSize()) + 2))
    const int64_t _size_force_a_0 = int64_t(int64_c(force_a->xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(force_a->ySize()) + 2))
    const int64_t _size_force_a_1 = int64_t(int64_c(force_a->ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(force_a->zSize()) + 2))
    const int64_t _size_force_a_2 = int64_t(int64_c(force_a->zSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
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
    internal_6ed8e43fbb9656b49ba1d2a7ec0e653c::colorgradientcollidesweepdoubleprecisionavx_colorgradientcollidesweepdoubleprecisionavx(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_even_a, omega_even_b, omega_odd_a, omega_odd_b, omega_shear_a, omega_shear_b, sigma);
    
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

    auto rho_a = block->getData< field::GhostLayerField<double, 1> >(rho_aID);
    auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
    auto rho_b = block->getData< field::GhostLayerField<double, 1> >(rho_bID);
    auto force_b = block->getData< field::GhostLayerField<double, 3> >(force_bID);
    auto pdfs_b = block->getData< field::GhostLayerField<double, 19> >(pdfs_bID);
    auto pdfs_a = block->getData< field::GhostLayerField<double, 19> >(pdfs_aID);
    auto phasefield = block->getData< field::GhostLayerField<double, 1> >(phasefieldID);
    auto force_a = block->getData< field::GhostLayerField<double, 3> >(force_aID);

    auto & beta = this->beta_;
    auto & omega_odd_b = this->omega_odd_b_;
    auto & sigma = this->sigma_;
    auto & omega_odd_a = this->omega_odd_a_;
    auto & omega_even_b = this->omega_even_b_;
    auto & omega_shear_b = this->omega_shear_b_;
    auto & omega_even_a = this->omega_even_a_;
    auto & omega_shear_a = this->omega_shear_a_;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    double * RESTRICT const _data_force_a = force_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT const _data_phasefield = phasefield->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) phasefield->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_a = rho_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_b = rho_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) velocity->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
    const int64_t _size_force_a_0 = int64_t(int64_c(ci.xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
    const int64_t _size_force_a_1 = int64_t(int64_c(ci.ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
    const int64_t _size_force_a_2 = int64_t(int64_c(ci.zSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
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
    internal_6ed8e43fbb9656b49ba1d2a7ec0e653c::colorgradientcollidesweepdoubleprecisionavx_colorgradientcollidesweepdoubleprecisionavx(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_even_a, omega_even_b, omega_odd_a, omega_odd_b, omega_shear_a, omega_shear_b, sigma);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
