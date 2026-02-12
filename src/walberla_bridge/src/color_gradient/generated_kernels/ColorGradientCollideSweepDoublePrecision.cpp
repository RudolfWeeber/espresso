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
//! \\file ColorGradientCollideSweepDoublePrecision.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 17fc54c872bd8ceabf271a7e9e636c7c583f55af


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "ColorGradientCollideSweepDoublePrecision.h"




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


namespace internal_c748f1e44263b1ff8564a6b7259a4d4c {
static FUNC_PREFIX void colorgradientcollidesweepdoubleprecision_colorgradientcollidesweepdoubleprecision(double * RESTRICT const _data_force_a, double * RESTRICT const _data_force_b, double * RESTRICT  _data_pdfs_a, double * RESTRICT  _data_pdfs_b, double * RESTRICT const _data_phasefield, double * RESTRICT const _data_rho_a, double * RESTRICT const _data_rho_b, double * RESTRICT const _data_velocity, int64_t const _size_force_a_0, int64_t const _size_force_a_1, int64_t const _size_force_a_2, int64_t const _stride_force_a_0, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_0, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_0, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_b_0, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_phasefield_0, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_0, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_0, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3, double omega_even_a, double omega_even_b, double omega_odd_a, double omega_odd_b, double omega_shear_a, double omega_shear_b)
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
      const double xi_244 = omega_odd_a*0.5;
      const double xi_276 = omega_odd_a*0.25;
      const double xi_280 = omega_shear_a*0.25;
      const double xi_365 = omega_odd_b*0.5;
      const double xi_389 = omega_odd_b*0.25;
      const double xi_393 = omega_shear_b*0.25;
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
            for (int64_t ctr_0 = 1; ctr_0 < _size_force_a_0 - 1; ctr_0 += 1)
            {
               const double xi_185 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
               const double xi_186 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_187 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
               const double xi_188 = -0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
               const double xi_189 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_190 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_191 = -xi_187 + xi_188 + xi_189 + xi_190;
               const double xi_192 = -0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
               const double xi_193 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
               const double xi_194 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_195 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_196 = xi_192 - xi_193 + xi_194 + xi_195;
               const double xi_197 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
               const double xi_198 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_199 = -xi_197 + xi_198;
               const double xi_200 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
               const double xi_201 = -xi_200;
               const double xi_202 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_203 = xi_201 + xi_202;
               const double xi_204 = -0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
               const double xi_205 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_206 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
               const double xi_207 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_208 = -xi_206 + xi_207;
               const double xi_209 = xi_204 + xi_205 + xi_208;
               const double xi_210 = xi_203 + xi_209;
               const double xi_211 = -xi_185 + xi_186 + xi_191 + xi_196 + xi_199 + xi_210;
               const double xi_212 = xi_211*xi_211;
               const double xi_213 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
               const double xi_214 = -xi_213;
               const double xi_215 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
               const double xi_216 = xi_201 + xi_202 + xi_214 + xi_215;
               const double xi_217 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
               const double xi_218 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
               const double xi_219 = -xi_217 + xi_218;
               const double xi_220 = xi_185 - xi_186;
               const double xi_221 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
               const double xi_222 = -xi_221;
               const double xi_223 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
               const double xi_224 = xi_206 - xi_207 + xi_222 + xi_223;
               const double xi_225 = xi_191 + xi_193 - xi_195 + xi_216 + xi_219 + xi_220 + xi_224;
               const double xi_226 = xi_225*xi_225;
               const double xi_227 = xi_217 - xi_218;
               const double xi_228 = xi_214 + xi_215;
               const double xi_229 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
               const double xi_230 = -xi_229;
               const double xi_231 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
               const double xi_232 = xi_197 - xi_198 + xi_230 + xi_231;
               const double xi_233 = xi_187 + xi_188 + xi_189 - xi_190 + xi_196 + xi_220 + xi_227 + xi_228 + xi_232;
               const double xi_234 = xi_233*xi_233;
               const double xi_235 = xi_212 + xi_226 + xi_234;
               const double xi_236 = pow(xi_235, 0.5);
               const bool xi_237 = xi_236 > 0.0;
               const double xi_238 = ((xi_237) ? (xi_236*0.049382716049382713): (0.0));
               const double xi_250 = ((1.0) / (xi_235));
               const double xi_251 = xi_236*0.22222222222222221;
               const double xi_252 = ((xi_237) ? (xi_251*(xi_226*xi_250*0.055555555555555552 - 0.018518518518518517)): (0.0));
               const double xi_256 = ((xi_237) ? (xi_251*(xi_250*0.055555555555555552*(xi_225*xi_225) - 0.018518518518518517)): (0.0));
               const double xi_264 = ((xi_237) ? (xi_251*(xi_250*0.055555555555555552*(xi_233*xi_233) - 0.018518518518518517)): (0.0));
               const double xi_267 = ((xi_237) ? (xi_251*(xi_234*xi_250*0.055555555555555552 - 0.018518518518518517)): (0.0));
               const double xi_273 = ((xi_237) ? (xi_251*(xi_250*0.055555555555555552*(xi_211*xi_211) - 0.018518518518518517)): (0.0));
               const double xi_275 = ((xi_237) ? (xi_251*(xi_212*xi_250*0.055555555555555552 - 0.018518518518518517)): (0.0));
               const double xi_283 = 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_284 = 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
               const double xi_285 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
               const double xi_286 = xi_192 + xi_194;
               const double xi_287 = xi_232 + xi_286;
               const double xi_288 = -0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_289 = xi_221 - xi_223 + xi_288;
               const double xi_290 = xi_200 - xi_202 + xi_208 - xi_283 + xi_284 + xi_285 + xi_287 + xi_289;
               const double xi_291 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776*(xi_290*xi_290) - 0.037037037037037035)): (0.0));
               const double xi_295 = 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_296 = 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
               const double xi_297 = -0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_298 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
               const double xi_299 = xi_203 + xi_224 + xi_287 - xi_295 + xi_296 + xi_297 + xi_298;
               const double xi_300 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776*(xi_299*xi_299) - 0.037037037037037035)): (0.0));
               const double xi_301 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776*(xi_299*xi_299) - 0.037037037037037035)): (0.0));
               const double xi_303 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776*(xi_290*xi_290) - 0.037037037037037035)): (0.0));
               const double xi_306 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_307 = xi_199 + xi_204 + xi_205 + xi_286;
               const double xi_308 = xi_213 - xi_215 + xi_295 - xi_296;
               const double xi_309 = xi_227 + xi_289 + xi_306 + xi_307 + xi_308;
               const double xi_310 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776*(xi_309*xi_309) - 0.037037037037037035)): (0.0));
               const double xi_316 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_317 = xi_219 + xi_283 - xi_284;
               const double xi_318 = xi_222 + xi_223 + xi_228 + xi_297 + xi_307 + xi_316 + xi_317;
               const double xi_319 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776*(xi_318*xi_318) - 0.037037037037037035)): (0.0));
               const double xi_324 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_325 = xi_209 + xi_216 + xi_227 + xi_230 + xi_231 + xi_288 + xi_297 + xi_324;
               const double xi_326 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776*(xi_325*xi_325) - 0.037037037037037035)): (0.0));
               const double xi_330 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
               const double xi_331 = xi_210 + xi_229 - xi_231 + xi_308 + xi_317 + xi_330;
               const double xi_332 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776*(xi_331*xi_331) - 0.037037037037037035)): (0.0));
               const double xi_333 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776*(xi_318*xi_318) - 0.037037037037037035)): (0.0));
               const double xi_335 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776*(xi_309*xi_309) - 0.037037037037037035)): (0.0));
               const double xi_336 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776*(xi_331*xi_331) - 0.037037037037037035)): (0.0));
               const double xi_338 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776*(xi_325*xi_325) - 0.037037037037037035)): (0.0));
               const double xia_1_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3];
               const double xia_2_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2];
               const double xi_17 = xia_2_collide*0.16666666666666666;
               const double xi_18 = omega_odd_a*xi_17;
               const double xi_32 = xia_2_collide*0.083333333333333329;
               const double xi_34 = xi_33*xia_2_collide;
               const double xi_35 = xi_32 - xi_34;
               const double xi_51 = -xi_32 + xi_34;
               const double xia_3_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3];
               const double xia_4_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3];
               const double xia_5_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3];
               const double xia_6_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3];
               const double xia_7_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
               const double xia_8_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3];
               const double xi_241 = -xia_8_collide;
               const double xi_242 = xi_241 + xia_1_collide;
               const double xia_9_collide = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
               const double xia_10_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3];
               const double xia_11_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
               const double xia_12_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3];
               const double xi_246 = -xia_12_collide;
               const double xia_13_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3];
               const double xi_177 = xia_10_collide + xia_13_collide;
               const double xi_178 = xi_177 + xia_1_collide + xia_8_collide;
               const double xia_14_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3];
               const double xi_4 = xia_14_collide*0.16666666666666666;
               const double xi_7 = omega_odd_a*xi_4;
               const double xi_37 = xia_14_collide*0.083333333333333329;
               const double xi_38 = xi_33*xia_14_collide;
               const double xi_39 = -xi_37 + xi_38;
               const double xi_41 = xia_14_collide*0.25;
               const double xi_46 = xi_45*xia_14_collide;
               const double xi_53 = xi_37 - xi_38;
               const double xia_15_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3];
               const double xi_260 = -xia_15_collide;
               const double xi_261 = xi_260 + xia_6_collide;
               const double xia_16_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
               const double xia_17_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3];
               const double xia_18_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3];
               const double xia_19_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3];
               const double xia_20_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3];
               const double xi_26 = xia_20_collide*0.16666666666666666;
               const double xi_28 = omega_odd_a*xi_26;
               const double xi_55 = xia_20_collide*0.083333333333333329;
               const double xi_56 = xi_33*xia_20_collide;
               const double xi_57 = xi_55 - xi_56;
               const double xi_70 = -xi_55 + xi_56;
               const double xia_21_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3];
               const double xi_170 = xia_21_collide + xia_3_collide;
               const double xia_22_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3];
               const double xia_23_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3];
               const double xi_174 = xia_17_collide + xia_23_collide;
               const double xi_175 = xi_174 + xia_15_collide + xia_6_collide;
               const double xia_24_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3];
               const double xi_171 = xi_170 + xia_12_collide + xia_24_collide;
               const double xi_247 = xi_246 + xia_24_collide;
               const double xia_25_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2];
               const double xia_26_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3];
               const double xi_6 = xi_5*xia_14_collide;
               const double xi_19 = xi_5*xia_2_collide;
               const double xi_27 = xi_5*xia_20_collide;
               const double rho_a_collide = xia_9_collide;
               const double xi_167 = rho_a_collide*-0.092592592592592601;
               const double xi_180 = rho_a_collide*-0.27777777777777773;
               const double xi_181 = xi_178 + xi_180;
               const double xi_239 = rho_a_collide*0.33333333333333331;
               const double u_0_a_collide = xia_7_collide;
               const double xi_0 = u_0_a_collide*xia_2_collide;
               const double xi_13 = xi_0*0.083333333333333329;
               const double xi_14 = omega_even_a*xi_13 + xi_0*-0.16666666666666666;
               const double xi_20 = xi_0*0.33333333333333331;
               const double xi_42 = u_0_a_collide*xi_41;
               const double xi_47 = u_0_a_collide*xi_46;
               const double xi_58 = -xi_13;
               const double xi_61 = omega_even_a*xi_0*0.041666666666666664;
               const double xi_73 = u_0_a_collide*xia_20_collide;
               const double xi_74 = xi_73*0.25;
               const double xi_77 = xi_45*xi_73;
               const double xi_166 = u_0_a_collide*u_0_a_collide;
               const double xi_173 = rho_a_collide*xi_166*-0.33333333333333331;
               const double xi_182 = omega_shear_a*(rho_a_collide*xi_166 - xi_175 - xi_181 - xia_19_collide - xia_26_collide);
               const double xi_257 = u_0_a_collide*xi_239;
               const double xi_258 = xi_242 + xi_257 + xia_10_collide - xia_13_collide;
               const double xi_259 = xi_244*xi_258;
               const double xi_262 = xi_257 + xi_261 - xia_17_collide + xia_23_collide;
               const double xi_263 = xi_244*xi_262;
               const double xi_277 = xi_258*xi_276;
               const double xi_278 = -xi_277;
               const double xi_320 = xi_262*xi_276;
               const double xi_321 = -xi_320;
               const double u_1_a_collide = xia_11_collide;
               const double xi_1 = u_1_a_collide*xia_14_collide;
               const double xi_8 = xi_1*0.33333333333333331;
               const double xi_21 = xi_1*0.16666666666666666;
               const double xi_22 = xi_1*0.083333333333333329;
               const double xi_23 = omega_even_a*xi_22;
               const double xi_24 = -xi_21 + xi_23;
               const double xi_30 = xi_14 + xi_24;
               const double xi_43 = u_1_a_collide*0.25;
               const double xi_44 = xi_43*xia_2_collide;
               const double xi_48 = u_1_a_collide*xi_45;
               const double xi_49 = xi_48*xia_2_collide;
               const double xi_50 = xi_42 + xi_44 - xi_47 - xi_49;
               const double xi_52 = -xi_42 - xi_44 + xi_47 + xi_49;
               const double xi_59 = -xi_23;
               const double xi_63 = xi_43*xia_20_collide;
               const double xi_65 = xi_48*xia_20_collide;
               const double xi_71 = omega_even_a*u_1_a_collide*xia_14_collide*-0.041666666666666664;
               const double xi_168 = u_1_a_collide*u_1_a_collide;
               const double xi_169 = rho_a_collide*xi_168*-0.33333333333333331 + xi_167;
               const double xi_183 = omega_shear_a*(rho_a_collide*xi_168 - xi_171 - xi_181 - xia_4_collide - xia_5_collide);
               const double xi_240 = u_1_a_collide*xi_239;
               const double xi_243 = xi_240 + xi_242 - xia_10_collide + xia_13_collide;
               const double xi_245 = xi_243*xi_244;
               const double xi_248 = xi_240 + xi_247 - xia_21_collide + xia_3_collide;
               const double xi_249 = xi_244*xi_248;
               const double xi_279 = rho_a_collide*u_1_a_collide;
               const double xi_281 = xi_280*(u_0_a_collide*xi_279 + xi_177 + xi_241 - xia_1_collide);
               const double xi_282 = -xi_281;
               const double xi_293 = xi_243*xi_276;
               const double xi_304 = xi_248*xi_276;
               const double xi_314 = -xi_304;
               const double u_2_a_collide = xia_16_collide;
               const double xi_2 = u_2_a_collide*xia_20_collide;
               const double xi_9 = xi_2*0.16666666666666666;
               const double xi_10 = xi_2*0.083333333333333329;
               const double xi_11 = omega_even_a*xi_10;
               const double xi_12 = xi_11 - xi_9;
               const double xi_15 = xi_12 + xi_14;
               const double xi_16 = omega_even_a*xi_8 + omega_shear_a*xi_1*-0.5 + xi_15 + xi_8;
               const double xi_25 = omega_even_a*xi_20 + omega_shear_a*xi_0*-0.5 + xi_12 + xi_20 + xi_24;
               const double xi_29 = xi_2*0.33333333333333331;
               const double xi_31 = omega_even_a*xi_29 + omega_shear_a*xi_2*-0.5 + xi_29 + xi_30;
               const double xi_36 = omega_even_a*u_2_a_collide*xia_20_collide*-0.041666666666666664;
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
               const double xi_75 = u_2_a_collide*xia_2_collide;
               const double xi_76 = xi_75*0.25;
               const double xi_78 = xi_45*xi_75;
               const double xi_79 = xi_74 + xi_76 - xi_77 - xi_78;
               const double xi_80 = -xi_74 - xi_76 + xi_77 + xi_78;
               const double xi_81 = xi_15 + xi_22 + xi_51 + xi_71;
               const double xi_164 = u_2_a_collide*u_2_a_collide;
               const double xi_165 = rho_a_collide*xi_164*-0.33333333333333331;
               const double xi_172 = omega_even_a*(rho_a_collide*xi_166*-0.16666666666666666 - xi_165 - xi_169 - xi_171);
               const double xi_176 = omega_even_a*(rho_a_collide*xi_168*-0.16666666666666666 - xi_165 - xi_167 - xi_173 - xi_175);
               const double xi_179 = omega_even_a*(rho_a_collide*xi_164*-0.16666666666666666 - xi_169 - xi_173 - xi_178);
               const double xi_184 = omega_shear_a*(rho_a_collide*xi_164 - xi_171 - xi_175 - xi_180 - xia_18_collide - xia_22_collide);
               const double xi_253 = xi_179*-0.5;
               const double xi_254 = xi_172*-0.5;
               const double xi_255 = xi_183*0.5 + xi_253 + xi_254;
               const double xi_265 = xi_176*-0.5;
               const double xi_266 = xi_182*0.5 + xi_253 + xi_265;
               const double xi_268 = u_2_a_collide*xi_239;
               const double xi_269 = xi_261 + xi_268 + xia_17_collide - xia_23_collide;
               const double xi_270 = xi_244*xi_269;
               const double xi_271 = xi_247 + xi_268 + xia_21_collide - xia_3_collide;
               const double xi_272 = xi_244*xi_271;
               const double xi_274 = xi_184*0.5 + xi_254 + xi_265;
               const double xi_292 = xi_179*0.25;
               const double xi_294 = xi_292 + xi_293;
               const double xi_302 = xi_292 - xi_293;
               const double xi_305 = xi_280*(u_2_a_collide*xi_279 + xi_170 + xi_246 - xia_24_collide);
               const double xi_311 = xi_172*0.25;
               const double xi_312 = xi_271*xi_276;
               const double xi_313 = xi_311 + xi_312;
               const double xi_315 = -xi_305;
               const double xi_322 = xi_280*(rho_a_collide*u_0_a_collide*u_2_a_collide + xi_174 + xi_260 - xia_6_collide);
               const double xi_323 = -xi_322;
               const double xi_327 = xi_176*0.25;
               const double xi_328 = xi_269*xi_276;
               const double xi_329 = xi_327 + xi_328;
               const double xi_334 = xi_311 - xi_312;
               const double xi_337 = xi_327 - xi_328;
               const double forceTerm_0_a_collide = omega_shear_a*u_0_a_collide*xia_2_collide + omega_shear_a*u_1_a_collide*xia_14_collide + omega_shear_a*u_2_a_collide*xia_20_collide - xi_0*xi_3 - xi_0 - xi_1*xi_3 - xi_1 - xi_2*xi_3 - xi_2;
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
               const double xib_1_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3];
               const double xib_2_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3];
               const double xib_3_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3];
               const double xi_108 = xib_3_collide*0.16666666666666666;
               const double xi_110 = omega_odd_b*xi_108;
               const double xi_137 = xib_3_collide*0.083333333333333329;
               const double xi_138 = xi_115*xib_3_collide;
               const double xi_139 = xi_137 - xi_138;
               const double xi_152 = -xi_137 + xi_138;
               const double xib_4_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3];
               const double xib_5_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
               const double xib_6_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3];
               const double xib_7_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3];
               const double xib_8_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2];
               const double xib_9_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3];
               const double xi_375 = -xib_9_collide;
               const double xi_376 = xi_375 + xib_7_collide;
               const double xib_10_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
               const double xib_11_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3];
               const double xib_12_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3];
               const double xi_367 = -xib_12_collide;
               const double xi_368 = xi_367 + xib_1_collide;
               const double xib_13_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3];
               const double xi_86 = xib_13_collide*0.16666666666666666;
               const double xi_89 = omega_odd_b*xi_86;
               const double xi_119 = xib_13_collide*0.083333333333333329;
               const double xi_120 = xi_115*xib_13_collide;
               const double xi_121 = -xi_119 + xi_120;
               const double xi_123 = xib_13_collide*0.25;
               const double xi_128 = xi_127*xib_13_collide;
               const double xi_135 = xi_119 - xi_120;
               const double xib_14_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3];
               const double xib_15_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3];
               const double xib_16_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
               const double xib_17_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3];
               const double xi_352 = xib_17_collide + xib_2_collide;
               const double xib_18_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3];
               const double xi_349 = xib_15_collide + xib_18_collide;
               const double xi_350 = xi_349 + xib_7_collide + xib_9_collide;
               const double xib_19_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3];
               const double xib_20_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3];
               const double xib_21_collide = _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
               const double xib_22_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3];
               const double xi_345 = xib_11_collide + xib_22_collide;
               const double xi_346 = xi_345 + xib_12_collide + xib_1_collide;
               const double xib_23_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2];
               const double xi_99 = xib_23_collide*0.16666666666666666;
               const double xi_100 = omega_odd_b*xi_99;
               const double xi_114 = xib_23_collide*0.083333333333333329;
               const double xi_116 = xi_115*xib_23_collide;
               const double xi_117 = xi_114 - xi_116;
               const double xi_133 = -xi_114 + xi_116;
               const double xib_24_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3];
               const double xib_25_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3];
               const double xib_26_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3];
               const double xi_353 = xi_352 + xib_25_collide + xib_26_collide;
               const double xi_362 = -xib_26_collide;
               const double xi_363 = xi_362 + xib_25_collide;
               const double xi_88 = xi_87*xib_13_collide;
               const double xi_101 = xi_87*xib_23_collide;
               const double xi_109 = xi_87*xib_3_collide;
               const double rho_b_collide = xib_21_collide;
               const double xi_342 = rho_b_collide*-0.092592592592592601;
               const double xi_355 = rho_b_collide*-0.27777777777777773;
               const double xi_357 = xi_346 + xi_355;
               const double xi_360 = rho_b_collide*0.33333333333333331;
               const double u_0_b_collide = xib_5_collide;
               const double xi_82 = u_0_b_collide*xib_23_collide;
               const double xi_95 = xi_82*0.083333333333333329;
               const double xi_96 = omega_even_b*xi_95 + xi_82*-0.16666666666666666;
               const double xi_102 = xi_82*0.33333333333333331;
               const double xi_124 = u_0_b_collide*xi_123;
               const double xi_129 = u_0_b_collide*xi_128;
               const double xi_140 = -xi_95;
               const double xi_143 = omega_even_b*xi_82*0.041666666666666664;
               const double xi_155 = u_0_b_collide*xib_3_collide;
               const double xi_156 = xi_155*0.25;
               const double xi_159 = xi_127*xi_155;
               const double xi_341 = u_0_b_collide*u_0_b_collide;
               const double xi_348 = rho_b_collide*xi_341*-0.33333333333333331;
               const double xi_356 = omega_shear_b*(rho_b_collide*xi_341 - xi_350 - xi_353 - xi_355 - xib_24_collide - xib_4_collide);
               const double xi_374 = u_0_b_collide*xi_360;
               const double xi_377 = xi_374 + xi_376 + xib_15_collide - xib_18_collide;
               const double xi_378 = xi_365*xi_377;
               const double xi_379 = xi_363 + xi_374 - xib_17_collide + xib_2_collide;
               const double xi_380 = xi_365*xi_379;
               const double xi_390 = xi_379*xi_389;
               const double xi_391 = -xi_390;
               const double xi_407 = xi_377*xi_389;
               const double xi_408 = -xi_407;
               const double u_1_b_collide = xib_10_collide;
               const double xi_83 = u_1_b_collide*xib_13_collide;
               const double xi_90 = xi_83*0.33333333333333331;
               const double xi_103 = xi_83*0.16666666666666666;
               const double xi_104 = xi_83*0.083333333333333329;
               const double xi_105 = omega_even_b*xi_104;
               const double xi_106 = -xi_103 + xi_105;
               const double xi_112 = xi_106 + xi_96;
               const double xi_125 = u_1_b_collide*0.25;
               const double xi_126 = xi_125*xib_23_collide;
               const double xi_130 = u_1_b_collide*xi_127;
               const double xi_131 = xi_130*xib_23_collide;
               const double xi_132 = xi_124 + xi_126 - xi_129 - xi_131;
               const double xi_134 = -xi_124 - xi_126 + xi_129 + xi_131;
               const double xi_141 = -xi_105;
               const double xi_145 = xi_125*xib_3_collide;
               const double xi_147 = xi_130*xib_3_collide;
               const double xi_153 = omega_even_b*u_1_b_collide*xib_13_collide*-0.041666666666666664;
               const double xi_339 = u_1_b_collide*u_1_b_collide;
               const double xi_340 = rho_b_collide*xi_339*-0.33333333333333331;
               const double xi_358 = omega_shear_b*(rho_b_collide*xi_339 - xi_353 - xi_357 - xib_19_collide - xib_6_collide);
               const double xi_361 = u_1_b_collide*xi_360;
               const double xi_364 = xi_361 + xi_363 + xib_17_collide - xib_2_collide;
               const double xi_366 = xi_364*xi_365;
               const double xi_369 = xi_361 + xi_368 - xib_11_collide + xib_22_collide;
               const double xi_370 = xi_365*xi_369;
               const double xi_392 = rho_b_collide*u_1_b_collide;
               const double xi_394 = xi_393*(u_0_b_collide*xi_392 + xi_352 + xi_362 - xib_25_collide);
               const double xi_395 = -xi_394;
               const double xi_397 = xi_364*xi_389;
               const double xi_400 = xi_369*xi_389;
               const double xi_405 = -xi_400;
               const double u_2_b_collide = xib_16_collide;
               const double xi_84 = u_2_b_collide*xib_3_collide;
               const double xi_91 = xi_84*0.16666666666666666;
               const double xi_92 = xi_84*0.083333333333333329;
               const double xi_93 = omega_even_b*xi_92;
               const double xi_94 = -xi_91 + xi_93;
               const double xi_97 = xi_94 + xi_96;
               const double xi_98 = omega_even_b*xi_90 + omega_shear_b*xi_83*-0.5 + xi_90 + xi_97;
               const double xi_107 = omega_even_b*xi_102 + omega_shear_b*xi_82*-0.5 + xi_102 + xi_106 + xi_94;
               const double xi_111 = xi_84*0.33333333333333331;
               const double xi_113 = omega_even_b*xi_111 + omega_shear_b*xi_84*-0.5 + xi_111 + xi_112;
               const double xi_118 = omega_even_b*u_2_b_collide*xib_3_collide*-0.041666666666666664;
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
               const double xi_157 = u_2_b_collide*xib_23_collide;
               const double xi_158 = xi_157*0.25;
               const double xi_160 = xi_127*xi_157;
               const double xi_161 = xi_156 + xi_158 - xi_159 - xi_160;
               const double xi_162 = -xi_156 - xi_158 + xi_159 + xi_160;
               const double xi_163 = xi_104 + xi_133 + xi_153 + xi_97;
               const double xi_343 = u_2_b_collide*u_2_b_collide;
               const double xi_344 = rho_b_collide*xi_343*-0.33333333333333331 + xi_342;
               const double xi_347 = omega_even_b*(rho_b_collide*xi_341*-0.16666666666666666 - xi_340 - xi_344 - xi_346);
               const double xi_351 = omega_even_b*(rho_b_collide*xi_339*-0.16666666666666666 - xi_344 - xi_348 - xi_350);
               const double xi_354 = omega_even_b*(rho_b_collide*xi_343*-0.16666666666666666 - xi_340 - xi_342 - xi_348 - xi_353);
               const double xi_359 = omega_shear_b*(rho_b_collide*xi_343 - xi_350 - xi_357 - xib_14_collide - xib_20_collide);
               const double xi_371 = xi_347*-0.5;
               const double xi_372 = xi_354*-0.5;
               const double xi_373 = xi_358*0.5 + xi_371 + xi_372;
               const double xi_381 = xi_351*-0.5;
               const double xi_382 = xi_356*0.5 + xi_372 + xi_381;
               const double xi_383 = u_2_b_collide*xi_360;
               const double xi_384 = xi_368 + xi_383 + xib_11_collide - xib_22_collide;
               const double xi_385 = xi_365*xi_384;
               const double xi_386 = xi_376 + xi_383 - xib_15_collide + xib_18_collide;
               const double xi_387 = xi_365*xi_386;
               const double xi_388 = xi_359*0.5 + xi_371 + xi_381;
               const double xi_396 = xi_354*0.25;
               const double xi_398 = xi_396 + xi_397;
               const double xi_399 = xi_396 - xi_397;
               const double xi_401 = xi_393*(u_2_b_collide*xi_392 + xi_345 + xi_367 - xib_1_collide);
               const double xi_402 = xi_347*0.25;
               const double xi_403 = xi_384*xi_389;
               const double xi_404 = xi_402 + xi_403;
               const double xi_406 = -xi_401;
               const double xi_409 = xi_393*(rho_b_collide*u_0_b_collide*u_2_b_collide + xi_349 + xi_375 - xib_7_collide);
               const double xi_410 = -xi_409;
               const double xi_411 = xi_351*0.25;
               const double xi_412 = xi_386*xi_389;
               const double xi_413 = xi_411 + xi_412;
               const double xi_414 = xi_402 - xi_403;
               const double xi_415 = xi_411 - xi_412;
               const double forceTerm_0_b_collide = omega_shear_b*u_0_b_collide*xib_23_collide + omega_shear_b*u_1_b_collide*xib_13_collide + omega_shear_b*u_2_b_collide*xib_3_collide - xi_82*xi_85 - xi_82 - xi_83*xi_85 - xi_83 - xi_84*xi_85 - xi_84;
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
               const double tmp_a0 = forceTerm_0_a_collide + xi_172 + xi_176 + xi_179 - xi_182 - xi_183 - xi_184 + xi_238 + xia_25_collide;
               const double tmp_a1 = forceTerm_1_a_collide - xi_245 - xi_249 + xi_252 + xi_255 + xia_5_collide;
               const double tmp_a2 = forceTerm_2_a_collide + xi_245 + xi_249 + xi_255 + xi_256 + xia_4_collide;
               const double tmp_a3 = forceTerm_3_a_collide + xi_259 + xi_263 + xi_264 + xi_266 + xia_26_collide;
               const double tmp_a4 = forceTerm_4_a_collide - xi_259 - xi_263 + xi_266 + xi_267 + xia_19_collide;
               const double tmp_a5 = forceTerm_5_a_collide - xi_270 - xi_272 + xi_273 + xi_274 + xia_22_collide;
               const double tmp_a6 = forceTerm_6_a_collide + xi_270 + xi_272 + xi_274 + xi_275 + xia_18_collide;
               const double tmp_a7 = forceTerm_7_a_collide + xi_278 + xi_282 + xi_291 + xi_294 + xia_10_collide;
               const double tmp_a8 = forceTerm_8_a_collide + xi_277 + xi_281 + xi_294 + xi_300 + xia_8_collide;
               const double tmp_a9 = forceTerm_9_a_collide + xi_278 + xi_281 + xi_301 + xi_302 + xia_1_collide;
               const double tmp_a10 = forceTerm_10_a_collide + xi_277 + xi_282 + xi_302 + xi_303 + xia_13_collide;
               const double tmp_a11 = forceTerm_11_a_collide + xi_304 + xi_305 + xi_310 + xi_313 + xia_12_collide;
               const double tmp_a12 = forceTerm_12_a_collide + xi_313 + xi_314 + xi_315 + xi_319 + xia_3_collide;
               const double tmp_a13 = forceTerm_13_a_collide + xi_321 + xi_323 + xi_326 + xi_329 + xia_23_collide;
               const double tmp_a14 = forceTerm_14_a_collide + xi_320 + xi_322 + xi_329 + xi_332 + xia_15_collide;
               const double tmp_a15 = forceTerm_15_a_collide + xi_304 + xi_315 + xi_333 + xi_334 + xia_21_collide;
               const double tmp_a16 = forceTerm_16_a_collide + xi_305 + xi_314 + xi_334 + xi_335 + xia_24_collide;
               const double tmp_a17 = forceTerm_17_a_collide + xi_321 + xi_322 + xi_336 + xi_337 + xia_6_collide;
               const double tmp_a18 = forceTerm_18_a_collide + xi_320 + xi_323 + xi_337 + xi_338 + xia_17_collide;
               const double tmp_b0 = forceTerm_0_b_collide + xi_238 + xi_347 + xi_351 + xi_354 - xi_356 - xi_358 - xi_359 + xib_8_collide;
               const double tmp_b1 = forceTerm_1_b_collide + xi_252 - xi_366 - xi_370 + xi_373 + xib_6_collide;
               const double tmp_b2 = forceTerm_2_b_collide + xi_256 + xi_366 + xi_370 + xi_373 + xib_19_collide;
               const double tmp_b3 = forceTerm_3_b_collide + xi_264 + xi_378 + xi_380 + xi_382 + xib_24_collide;
               const double tmp_b4 = forceTerm_4_b_collide + xi_267 - xi_378 - xi_380 + xi_382 + xib_4_collide;
               const double tmp_b5 = forceTerm_5_b_collide + xi_273 - xi_385 - xi_387 + xi_388 + xib_14_collide;
               const double tmp_b6 = forceTerm_6_b_collide + xi_275 + xi_385 + xi_387 + xi_388 + xib_20_collide;
               const double tmp_b7 = forceTerm_7_b_collide + xi_291 + xi_391 + xi_395 + xi_398 + xib_2_collide;
               const double tmp_b8 = forceTerm_8_b_collide + xi_300 + xi_390 + xi_394 + xi_398 + xib_26_collide;
               const double tmp_b9 = forceTerm_9_b_collide + xi_301 + xi_391 + xi_394 + xi_399 + xib_25_collide;
               const double tmp_b10 = forceTerm_10_b_collide + xi_303 + xi_390 + xi_395 + xi_399 + xib_17_collide;
               const double tmp_b11 = forceTerm_11_b_collide + xi_310 + xi_400 + xi_401 + xi_404 + xib_12_collide;
               const double tmp_b12 = forceTerm_12_b_collide + xi_319 + xi_404 + xi_405 + xi_406 + xib_22_collide;
               const double tmp_b13 = forceTerm_13_b_collide + xi_326 + xi_408 + xi_410 + xi_413 + xib_15_collide;
               const double tmp_b14 = forceTerm_14_b_collide + xi_332 + xi_407 + xi_409 + xi_413 + xib_9_collide;
               const double tmp_b15 = forceTerm_15_b_collide + xi_333 + xi_400 + xi_406 + xi_414 + xib_11_collide;
               const double tmp_b16 = forceTerm_16_b_collide + xi_335 + xi_401 + xi_405 + xi_414 + xib_1_collide;
               const double tmp_b17 = forceTerm_17_b_collide + xi_336 + xi_408 + xi_409 + xi_415 + xib_7_collide;
               const double tmp_b18 = forceTerm_18_b_collide + xi_338 + xi_407 + xi_410 + xi_415 + xib_18_collide;
               const double xirecolor_0 = tmp_a0 + tmp_b0;
               const double xirecolor_1 = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
               const double xirecolor_2 = ((1.0) / (xirecolor_1));
               const double xi_437 = xirecolor_2*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
               const double xirecolor_3 = xirecolor_2*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
               const double xirecolor_4 = tmp_a1 + tmp_b1;
               const double xirecolor_5 = xi_195;
               const double xirecolor_6 = xi_193;
               const double xirecolor_7 = xi_213;
               const double xi_416 = -xirecolor_7;
               const double xirecolor_8 = xi_416;
               const double xirecolor_9 = xi_200;
               const double xirecolor_10 = -xirecolor_9;
               const double xirecolor_11 = xi_202;
               const double xi_417 = xirecolor_10 + xirecolor_11;
               const double xirecolor_12 = xi_215;
               const double xi_418 = xirecolor_12 + xirecolor_8;
               const double xirecolor_13 = xi_417 + xi_418;
               const double xirecolor_14 = xi_217;
               const double xirecolor_15 = xi_218;
               const double xi_419 = xirecolor_14 - xirecolor_15;
               const double xirecolor_16 = -xi_419;
               const double xirecolor_17 = xi_187;
               const double xirecolor_18 = xi_188;
               const double xirecolor_19 = xi_189;
               const double xi_420 = xirecolor_18 + xirecolor_19;
               const double xirecolor_20 = xi_190;
               const double xirecolor_21 = xi_420 - xirecolor_17 + xirecolor_20;
               const double xirecolor_22 = xi_186;
               const double xirecolor_23 = xi_185;
               const double xi_421 = xirecolor_22 - xirecolor_23;
               const double xirecolor_24 = -xi_421;
               const double xirecolor_25 = xi_221;
               const double xirecolor_26 = -xirecolor_25;
               const double xirecolor_27 = xi_207;
               const double xirecolor_28 = xi_206;
               const double xirecolor_29 = xi_223;
               const double xi_422 = xirecolor_26 + xirecolor_29;
               const double xirecolor_30 = xi_422 - xirecolor_27 + xirecolor_28;
               const double xirecolor_31 = xirecolor_13 + xirecolor_16 + xirecolor_21 + xirecolor_24 + xirecolor_30 - xirecolor_5 + xirecolor_6;
               const double xirecolor_32 = xi_192;
               const double xirecolor_33 = xi_194;
               const double xi_423 = xirecolor_32 + xirecolor_33;
               const double xirecolor_34 = xi_423 + xirecolor_5 - xirecolor_6;
               const double xirecolor_35 = xi_197;
               const double xirecolor_36 = xi_198;
               const double xi_424 = xirecolor_35 - xirecolor_36;
               const double xirecolor_37 = -xi_424;
               const double xirecolor_38 = xi_417;
               const double xirecolor_39 = xi_204;
               const double xirecolor_40 = xi_205;
               const double xi_425 = xirecolor_39 + xirecolor_40;
               const double xirecolor_41 = xirecolor_27 - xirecolor_28;
               const double xirecolor_42 = xi_425 + xirecolor_41;
               const double xirecolor_43 = xirecolor_38 + xirecolor_42;
               const double xirecolor_44 = xi_421 + xirecolor_21 + xirecolor_34 + xirecolor_37 + xirecolor_43;
               const double xirecolor_45 = xi_419;
               const double xirecolor_46 = xi_418;
               const double xirecolor_47 = xi_229;
               const double xirecolor_48 = -xirecolor_47;
               const double xirecolor_49 = xi_231;
               const double xi_426 = xirecolor_48 + xirecolor_49;
               const double xirecolor_50 = xi_424 + xi_426;
               const double xirecolor_51 = xi_420 + xirecolor_17 - xirecolor_20 + xirecolor_24 + xirecolor_34 + xirecolor_45 + xirecolor_46 + xirecolor_50;
               const double xirecolor_52 = pow(xirecolor_31*xirecolor_31 + xirecolor_44*xirecolor_44 + xirecolor_51*xirecolor_51, 0.5);
               const double xirecolor_53 = ((1.0) / (xirecolor_52));
               const double xi_427 = xirecolor_31*xirecolor_53;
               const double xi_428 = xirecolor_51*xirecolor_53;
               const double xi_429 = xirecolor_44*xirecolor_53;
               const bool xirecolor_54 = xirecolor_52 > 0.0;
               const double xirecolor_55 = 0.69999999999999996*((1.0) / (xirecolor_1*xirecolor_1))*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2]*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
               const double xirecolor_56 = xirecolor_55*(0.046296296296296259*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.046296296296296259*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
               const double xirecolor_57 = xirecolor_56*((xirecolor_54) ? (xi_427): (0.0));
               const double xirecolor_58 = tmp_a2 + tmp_b2;
               const double xirecolor_59 = xirecolor_55*(0.046296296296296266*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.046296296296296266*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
               const double xirecolor_60 = xirecolor_59*((xirecolor_54) ? (-xi_427): (0.0));
               const double xirecolor_61 = tmp_a3 + tmp_b3;
               const double xirecolor_62 = xirecolor_59*((xirecolor_54) ? (-xi_428): (0.0));
               const double xirecolor_63 = tmp_a4 + tmp_b4;
               const double xirecolor_64 = xirecolor_56*((xirecolor_54) ? (xi_428): (0.0));
               const double xirecolor_65 = tmp_a5 + tmp_b5;
               const double xirecolor_66 = xirecolor_56*((xirecolor_54) ? (-xi_429): (0.0));
               const double xirecolor_67 = tmp_a6 + tmp_b6;
               const double xirecolor_68 = xirecolor_59*((xirecolor_54) ? (xi_429): (0.0));
               const double xirecolor_69 = tmp_a7 + tmp_b7;
               const double xirecolor_70 = xi_283;
               const double xirecolor_71 = xi_284;
               const double xirecolor_72 = xi_423;
               const double xirecolor_73 = xirecolor_50 + xirecolor_72;
               const double xirecolor_74 = xi_288;
               const double xirecolor_75 = xirecolor_25 - xirecolor_29 + xirecolor_74;
               const double xirecolor_76 = xi_285 - xirecolor_11 + xirecolor_41 - xirecolor_70 + xirecolor_71 + xirecolor_73 + xirecolor_75 + xirecolor_9;
               const double xirecolor_77 = xirecolor_53*0.70710678118654757;
               const double xi_430 = xirecolor_76*xirecolor_77;
               const double xirecolor_78 = xirecolor_55*(0.02314814814814815*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.02314814814814815*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
               const double xirecolor_79 = xirecolor_78*((xirecolor_54) ? (-xi_430): (0.0));
               const double xirecolor_80 = tmp_a8 + tmp_b8;
               const double xirecolor_81 = xi_295;
               const double xirecolor_82 = xi_296;
               const double xi_431 = -xirecolor_81 + xirecolor_82;
               const double xirecolor_83 = xi_297;
               const double xirecolor_84 = xi_298 + xi_431 + xirecolor_30 + xirecolor_38 + xirecolor_73 + xirecolor_83;
               const double xi_432 = xirecolor_77*xirecolor_84;
               const double xirecolor_85 = xirecolor_78*((xirecolor_54) ? (xi_432): (0.0));
               const double xirecolor_86 = tmp_a9 + tmp_b9;
               const double xirecolor_87 = xirecolor_78*((xirecolor_54) ? (-xi_432): (0.0));
               const double xirecolor_88 = tmp_a10 + tmp_b10;
               const double xirecolor_89 = xirecolor_78*((xirecolor_54) ? (xi_430): (0.0));
               const double xirecolor_90 = tmp_a11 + tmp_b11;
               const double xirecolor_91 = xi_425 + xirecolor_37 + xirecolor_72;
               const double xirecolor_92 = -xi_416 - xi_431 - xirecolor_12;
               const double xirecolor_93 = xi_306 + xirecolor_45 + xirecolor_75 + xirecolor_91 + xirecolor_92;
               const double xi_433 = xirecolor_77*xirecolor_93;
               const double xirecolor_94 = xirecolor_78*((xirecolor_54) ? (-xi_433): (0.0));
               const double xirecolor_95 = tmp_a12 + tmp_b12;
               const double xirecolor_96 = xirecolor_16 + xirecolor_70 - xirecolor_71;
               const double xirecolor_97 = xi_316 + xi_422 + xirecolor_46 + xirecolor_83 + xirecolor_91 + xirecolor_96;
               const double xi_434 = xirecolor_77*xirecolor_97;
               const double xirecolor_98 = xirecolor_78*((xirecolor_54) ? (-xi_434): (0.0));
               const double xirecolor_99 = tmp_a13 + tmp_b13;
               const double xirecolor_100 = xi_324 + xi_426 + xirecolor_13 + xirecolor_42 + xirecolor_45 + xirecolor_74 + xirecolor_83;
               const double xi_435 = xirecolor_100*xirecolor_77;
               const double xirecolor_101 = xirecolor_78*((xirecolor_54) ? (-xi_435): (0.0));
               const double xirecolor_102 = tmp_a14 + tmp_b14;
               const double xirecolor_103 = xi_330 + xirecolor_43 + xirecolor_47 - xirecolor_49 + xirecolor_92 + xirecolor_96;
               const double xi_436 = xirecolor_103*xirecolor_77;
               const double xirecolor_104 = xirecolor_78*((xirecolor_54) ? (-xi_436): (0.0));
               const double xirecolor_105 = tmp_a15 + tmp_b15;
               const double xirecolor_106 = xirecolor_78*((xirecolor_54) ? (xi_434): (0.0));
               const double xirecolor_107 = tmp_a16 + tmp_b16;
               const double xirecolor_108 = xirecolor_78*((xirecolor_54) ? (xi_433): (0.0));
               const double xirecolor_109 = tmp_a17 + tmp_b17;
               const double xirecolor_110 = xirecolor_78*((xirecolor_54) ? (xi_436): (0.0));
               const double xirecolor_111 = tmp_a18 + tmp_b18;
               const double xirecolor_112 = xirecolor_78*((xirecolor_54) ? (xi_435): (0.0));
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2] = xirecolor_0*xirecolor_3;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3] = xirecolor_3*xirecolor_4 + xirecolor_57;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3] = xirecolor_3*xirecolor_58 + xirecolor_60;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3] = xirecolor_3*xirecolor_61 + xirecolor_62;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3] = xirecolor_3*xirecolor_63 + xirecolor_64;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3] = xirecolor_3*xirecolor_65 + xirecolor_66;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3] = xirecolor_3*xirecolor_67 + xirecolor_68;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3] = xirecolor_3*xirecolor_69 + xirecolor_79;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3] = xirecolor_3*xirecolor_80 + xirecolor_85;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3] = xirecolor_3*xirecolor_86 + xirecolor_87;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3] = xirecolor_3*xirecolor_88 + xirecolor_89;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3] = xirecolor_3*xirecolor_90 + xirecolor_94;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3] = xirecolor_3*xirecolor_95 + xirecolor_98;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3] = xirecolor_101 + xirecolor_3*xirecolor_99;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3] = xirecolor_102*xirecolor_3 + xirecolor_104;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3] = xirecolor_105*xirecolor_3 + xirecolor_106;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3] = xirecolor_107*xirecolor_3 + xirecolor_108;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3] = xirecolor_109*xirecolor_3 + xirecolor_110;
               _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3] = xirecolor_111*xirecolor_3 + xirecolor_112;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2] = xi_437*xirecolor_0;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3] = xi_437*xirecolor_4 - xirecolor_57;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3] = xi_437*xirecolor_58 - xirecolor_60;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3] = xi_437*xirecolor_61 - xirecolor_62;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3] = xi_437*xirecolor_63 - xirecolor_64;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3] = xi_437*xirecolor_65 - xirecolor_66;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3] = xi_437*xirecolor_67 - xirecolor_68;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3] = xi_437*xirecolor_69 - xirecolor_79;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3] = xi_437*xirecolor_80 - xirecolor_85;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3] = xi_437*xirecolor_86 - xirecolor_87;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3] = xi_437*xirecolor_88 - xirecolor_89;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3] = xi_437*xirecolor_90 - xirecolor_94;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3] = xi_437*xirecolor_95 - xirecolor_98;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3] = xi_437*xirecolor_99 - xirecolor_101;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3] = xi_437*xirecolor_102 - xirecolor_104;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3] = xi_437*xirecolor_105 - xirecolor_106;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3] = xi_437*xirecolor_107 - xirecolor_108;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3] = xi_437*xirecolor_109 - xirecolor_110;
               _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3] = xi_437*xirecolor_111 - xirecolor_112;
            }
         }
      }
   }
}
}


void ColorGradientCollideSweepDoublePrecision::run(IBlock * block)
{
   
    auto rho_a = block->getData< field::GhostLayerField<double, 1> >(rho_aID);
    auto phasefield = block->getData< field::GhostLayerField<double, 1> >(phasefieldID);
    auto force_a = block->getData< field::GhostLayerField<double, 3> >(force_aID);
    auto pdfs_a = block->getData< field::GhostLayerField<double, 19> >(pdfs_aID);
    auto rho_b = block->getData< field::GhostLayerField<double, 1> >(rho_bID);
    auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
    auto pdfs_b = block->getData< field::GhostLayerField<double, 19> >(pdfs_bID);
    auto force_b = block->getData< field::GhostLayerField<double, 3> >(force_bID);

    auto & omega_odd_b = this->omega_odd_b_;
    auto & omega_odd_a = this->omega_odd_a_;
    auto & omega_even_b = this->omega_even_b_;
    auto & omega_shear_a = this->omega_shear_a_;
    auto & omega_even_a = this->omega_even_a_;
    auto & omega_shear_b = this->omega_shear_b_;
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_a->nrOfGhostLayers()))
    double * RESTRICT const _data_force_a = force_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT const _data_phasefield = phasefield->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_a = rho_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_b = rho_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(-1, -1, -1, 0);
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
    internal_c748f1e44263b1ff8564a6b7259a4d4c::colorgradientcollidesweepdoubleprecision_colorgradientcollidesweepdoubleprecision(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, omega_even_a, omega_even_b, omega_odd_a, omega_odd_b, omega_shear_a, omega_shear_b);
    
}


void ColorGradientCollideSweepDoublePrecision::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block)
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

    auto & omega_odd_b = this->omega_odd_b_;
    auto & omega_odd_a = this->omega_odd_a_;
    auto & omega_even_b = this->omega_even_b_;
    auto & omega_shear_a = this->omega_shear_a_;
    auto & omega_even_a = this->omega_even_a_;
    auto & omega_shear_b = this->omega_shear_b_;
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
    double * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT const _data_phasefield = phasefield->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_a = rho_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_b = rho_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
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
    internal_c748f1e44263b1ff8564a6b7259a4d4c::colorgradientcollidesweepdoubleprecision_colorgradientcollidesweepdoubleprecision(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, omega_even_a, omega_even_b, omega_odd_a, omega_odd_b, omega_shear_a, omega_shear_b);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
