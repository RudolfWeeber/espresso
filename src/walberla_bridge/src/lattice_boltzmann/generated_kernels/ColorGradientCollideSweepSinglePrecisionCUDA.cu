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
//! \\file ColorGradientCollideSweepSinglePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 17fc54c872bd8ceabf271a7e9e636c7c583f55af


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "ColorGradientCollideSweepSinglePrecisionCUDA.h"




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


namespace internal_colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda {
static FUNC_PREFIX __launch_bounds__(256) void colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda(float * RESTRICT const _data_force_a, float * RESTRICT const _data_force_b, float * RESTRICT  _data_pdfs_a, float * RESTRICT  _data_pdfs_b, float * RESTRICT const _data_phasefield, float * RESTRICT const _data_rho_a, float * RESTRICT const _data_rho_b, float * RESTRICT const _data_velocity, int64_t const _size_force_a_0, int64_t const _size_force_a_1, int64_t const _size_force_a_2, int64_t const _stride_force_a_0, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_0, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_0, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_b_0, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_phasefield_0, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_0, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_0, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3, float omega_even_a, float omega_even_b, float omega_odd_a, float omega_odd_b, float omega_shear_a, float omega_shear_b)
{
   if (blockDim.x*blockIdx.x + threadIdx.x + 1 < _size_force_a_0 - 1 && blockDim.y*blockIdx.y + threadIdx.y + 1 < _size_force_a_1 - 1 && blockDim.z*blockIdx.z + threadIdx.z + 1 < _size_force_a_2 - 1)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x + 1;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y + 1;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z + 1;
      const float xi_3 = omega_even_a*0.5f;
      const float xi_33 = omega_odd_a*0.041666666666666664f;
      const float xi_45 = omega_shear_a*0.125f;
      const float xi_85 = omega_even_b*0.5f;
      const float xi_115 = omega_odd_b*0.041666666666666664f;
      const float xi_127 = omega_shear_b*0.125f;
      const float xi_185 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_186 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_187 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_188 = -0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_189 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_190 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_191 = -xi_187 + xi_188 + xi_189 + xi_190;
      const float xi_192 = -0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_193 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_194 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_195 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_196 = xi_192 - xi_193 + xi_194 + xi_195;
      const float xi_197 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_198 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_199 = -xi_197 + xi_198;
      const float xi_200 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_201 = -xi_200;
      const float xi_202 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_203 = xi_201 + xi_202;
      const float xi_204 = -0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_205 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_206 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_207 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_208 = -xi_206 + xi_207;
      const float xi_209 = xi_204 + xi_205 + xi_208;
      const float xi_210 = xi_203 + xi_209;
      const float xi_211 = -xi_185 + xi_186 + xi_191 + xi_196 + xi_199 + xi_210;
      const float xi_212 = xi_211*xi_211;
      const float xi_213 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_214 = -xi_213;
      const float xi_215 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_216 = xi_201 + xi_202 + xi_214 + xi_215;
      const float xi_217 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_218 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_219 = -xi_217 + xi_218;
      const float xi_220 = xi_185 - xi_186;
      const float xi_221 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_222 = -xi_221;
      const float xi_223 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_224 = xi_206 - xi_207 + xi_222 + xi_223;
      const float xi_225 = xi_191 + xi_193 - xi_195 + xi_216 + xi_219 + xi_220 + xi_224;
      const float xi_226 = xi_225*xi_225;
      const float xi_227 = xi_217 - xi_218;
      const float xi_228 = xi_214 + xi_215;
      const float xi_229 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const float xi_230 = -xi_229;
      const float xi_231 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const float xi_232 = xi_197 - xi_198 + xi_230 + xi_231;
      const float xi_233 = xi_187 + xi_188 + xi_189 - xi_190 + xi_196 + xi_220 + xi_227 + xi_228 + xi_232;
      const float xi_234 = xi_233*xi_233;
      const float xi_235 = xi_212 + xi_226 + xi_234;
      const float xi_236 = powf(xi_235, 0.5f);
      const bool xi_237 = xi_236 > 0.0f;
      const float xi_238 = ((xi_237) ? (xi_236*0.049382716049382713f): (0.0f));
      const float xi_244 = omega_odd_a*0.5f;
      const float xi_250 = ((1.0f) / (xi_235));
      const float xi_251 = xi_236*0.22222222222222221f;
      const float xi_252 = ((xi_237) ? (xi_251*(xi_226*xi_250*0.055555555555555552f - 0.018518518518518517f)): (0.0f));
      const float xi_256 = ((xi_237) ? (xi_251*(xi_250*0.055555555555555552f*(xi_225*xi_225) - 0.018518518518518517f)): (0.0f));
      const float xi_264 = ((xi_237) ? (xi_251*(xi_250*0.055555555555555552f*(xi_233*xi_233) - 0.018518518518518517f)): (0.0f));
      const float xi_267 = ((xi_237) ? (xi_251*(xi_234*xi_250*0.055555555555555552f - 0.018518518518518517f)): (0.0f));
      const float xi_273 = ((xi_237) ? (xi_251*(xi_250*0.055555555555555552f*(xi_211*xi_211) - 0.018518518518518517f)): (0.0f));
      const float xi_275 = ((xi_237) ? (xi_251*(xi_212*xi_250*0.055555555555555552f - 0.018518518518518517f)): (0.0f));
      const float xi_276 = omega_odd_a*0.25f;
      const float xi_280 = omega_shear_a*0.25f;
      const float xi_283 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_284 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_285 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_286 = xi_192 + xi_194;
      const float xi_287 = xi_232 + xi_286;
      const float xi_288 = -0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_289 = xi_221 - xi_223 + xi_288;
      const float xi_290 = xi_200 - xi_202 + xi_208 - xi_283 + xi_284 + xi_285 + xi_287 + xi_289;
      const float xi_291 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776f*(xi_290*xi_290) - 0.037037037037037035f)): (0.0f));
      const float xi_295 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_296 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_297 = -0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_298 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_299 = xi_203 + xi_224 + xi_287 - xi_295 + xi_296 + xi_297 + xi_298;
      const float xi_300 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776f*(xi_299*xi_299) - 0.037037037037037035f)): (0.0f));
      const float xi_301 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776f*(xi_299*xi_299) - 0.037037037037037035f)): (0.0f));
      const float xi_303 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776f*(xi_290*xi_290) - 0.037037037037037035f)): (0.0f));
      const float xi_306 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_307 = xi_199 + xi_204 + xi_205 + xi_286;
      const float xi_308 = xi_213 - xi_215 + xi_295 - xi_296;
      const float xi_309 = xi_227 + xi_289 + xi_306 + xi_307 + xi_308;
      const float xi_310 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776f*(xi_309*xi_309) - 0.037037037037037035f)): (0.0f));
      const float xi_316 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_317 = xi_219 + xi_283 - xi_284;
      const float xi_318 = xi_222 + xi_223 + xi_228 + xi_297 + xi_307 + xi_316 + xi_317;
      const float xi_319 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776f*(xi_318*xi_318) - 0.037037037037037035f)): (0.0f));
      const float xi_324 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_325 = xi_209 + xi_216 + xi_227 + xi_230 + xi_231 + xi_288 + xi_297 + xi_324;
      const float xi_326 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776f*(xi_325*xi_325) - 0.037037037037037035f)): (0.0f));
      const float xi_330 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_331 = xi_210 + xi_229 - xi_231 + xi_308 + xi_317 + xi_330;
      const float xi_332 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776f*(xi_331*xi_331) - 0.037037037037037035f)): (0.0f));
      const float xi_333 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776f*(xi_318*xi_318) - 0.037037037037037035f)): (0.0f));
      const float xi_335 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776f*(xi_309*xi_309) - 0.037037037037037035f)): (0.0f));
      const float xi_336 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776f*(xi_331*xi_331) - 0.037037037037037035f)): (0.0f));
      const float xi_338 = ((xi_237) ? (xi_251*(xi_250*0.027777777777777776f*(xi_325*xi_325) - 0.037037037037037035f)): (0.0f));
      const float xi_365 = omega_odd_b*0.5f;
      const float xi_389 = omega_odd_b*0.25f;
      const float xi_393 = omega_shear_b*0.25f;
      const float xia_1_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3];
      const float xia_2_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3];
      const float xi_241 = -xia_2_collide;
      const float xia_3_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3];
      const float xia_4_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3];
      const float xi_4 = xia_4_collide*0.16666666666666666f;
      const float xi_7 = omega_odd_a*xi_4;
      const float xi_37 = xia_4_collide*0.083333333333333329f;
      const float xi_38 = xi_33*xia_4_collide;
      const float xi_39 = -xi_37 + xi_38;
      const float xi_41 = xia_4_collide*0.25f;
      const float xi_46 = xi_45*xia_4_collide;
      const float xi_53 = xi_37 - xi_38;
      const float xia_5_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3];
      const float xi_260 = -xia_5_collide;
      const float xi_261 = xi_260 + xia_1_collide;
      const float xia_6_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3];
      const float xia_7_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3];
      const float xi_177 = xia_3_collide + xia_7_collide;
      const float xia_8_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3];
      const float xia_9_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3];
      const float xia_10_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2];
      const float xi_17 = xia_10_collide*0.16666666666666666f;
      const float xi_18 = omega_odd_a*xi_17;
      const float xi_32 = xia_10_collide*0.083333333333333329f;
      const float xi_34 = xi_33*xia_10_collide;
      const float xi_35 = xi_32 - xi_34;
      const float xi_51 = -xi_32 + xi_34;
      const float xia_11_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3];
      const float xi_178 = xi_177 + xia_11_collide + xia_2_collide;
      const float xi_242 = xi_241 + xia_11_collide;
      const float xia_12_collide = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
      const float xia_13_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3];
      const float xi_246 = -xia_13_collide;
      const float xia_14_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3];
      const float xia_15_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2];
      const float xia_16_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3];
      const float xia_17_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3];
      const float xia_18_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3];
      const float xia_19_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3];
      const float xia_20_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3];
      const float xi_170 = xia_14_collide + xia_20_collide;
      const float xia_21_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
      const float xia_22_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
      const float xia_23_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3];
      const float xi_171 = xi_170 + xia_13_collide + xia_23_collide;
      const float xi_247 = xi_246 + xia_23_collide;
      const float xia_24_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
      const float xia_25_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3];
      const float xi_174 = xia_19_collide + xia_25_collide;
      const float xi_175 = xi_174 + xia_1_collide + xia_5_collide;
      const float xia_26_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3];
      const float xi_26 = xia_26_collide*0.16666666666666666f;
      const float xi_28 = omega_odd_a*xi_26;
      const float xi_55 = xia_26_collide*0.083333333333333329f;
      const float xi_56 = xi_33*xia_26_collide;
      const float xi_57 = xi_55 - xi_56;
      const float xi_70 = -xi_55 + xi_56;
      const float rr_0_a_collide = 0.0f;
      const float xi_5 = rr_0_a_collide*0.25f;
      const float xi_6 = xi_5*xia_4_collide;
      const float xi_19 = xi_5*xia_10_collide;
      const float xi_27 = xi_5*xia_26_collide;
      const float rho_a_collide = xia_12_collide;
      const float xi_167 = rho_a_collide*-0.092592592592592601f;
      const float xi_180 = rho_a_collide*-0.27777777777777773f;
      const float xi_181 = xi_178 + xi_180;
      const float xi_239 = rho_a_collide*0.33333333333333331f;
      const float u_0_a_collide = xia_21_collide;
      const float xi_0 = u_0_a_collide*xia_10_collide;
      const float xi_13 = xi_0*0.083333333333333329f;
      const float xi_14 = omega_even_a*xi_13 + xi_0*-0.16666666666666666f;
      const float xi_20 = xi_0*0.33333333333333331f;
      const float xi_42 = u_0_a_collide*xi_41;
      const float xi_47 = u_0_a_collide*xi_46;
      const float xi_58 = -xi_13;
      const float xi_61 = omega_even_a*xi_0*0.041666666666666664f;
      const float xi_73 = u_0_a_collide*xia_26_collide;
      const float xi_74 = xi_73*0.25f;
      const float xi_77 = xi_45*xi_73;
      const float xi_166 = u_0_a_collide*u_0_a_collide;
      const float xi_173 = rho_a_collide*xi_166*-0.33333333333333331f;
      const float xi_182 = omega_shear_a*(rho_a_collide*xi_166 - xi_175 - xi_181 - xia_6_collide - xia_8_collide);
      const float xi_257 = u_0_a_collide*xi_239;
      const float xi_258 = xi_242 + xi_257 - xia_3_collide + xia_7_collide;
      const float xi_259 = xi_244*xi_258;
      const float xi_262 = xi_257 + xi_261 - xia_19_collide + xia_25_collide;
      const float xi_263 = xi_244*xi_262;
      const float xi_277 = xi_258*xi_276;
      const float xi_278 = -xi_277;
      const float xi_320 = xi_262*xi_276;
      const float xi_321 = -xi_320;
      const float u_1_a_collide = xia_22_collide;
      const float xi_1 = u_1_a_collide*xia_4_collide;
      const float xi_8 = xi_1*0.33333333333333331f;
      const float xi_21 = xi_1*0.16666666666666666f;
      const float xi_22 = xi_1*0.083333333333333329f;
      const float xi_23 = omega_even_a*xi_22;
      const float xi_24 = -xi_21 + xi_23;
      const float xi_30 = xi_14 + xi_24;
      const float xi_43 = u_1_a_collide*0.25f;
      const float xi_44 = xi_43*xia_10_collide;
      const float xi_48 = u_1_a_collide*xi_45;
      const float xi_49 = xi_48*xia_10_collide;
      const float xi_50 = xi_42 + xi_44 - xi_47 - xi_49;
      const float xi_52 = -xi_42 - xi_44 + xi_47 + xi_49;
      const float xi_59 = -xi_23;
      const float xi_63 = xi_43*xia_26_collide;
      const float xi_65 = xi_48*xia_26_collide;
      const float xi_71 = omega_even_a*u_1_a_collide*xia_4_collide*-0.041666666666666664f;
      const float xi_168 = u_1_a_collide*u_1_a_collide;
      const float xi_169 = rho_a_collide*xi_168*-0.33333333333333331f + xi_167;
      const float xi_183 = omega_shear_a*(rho_a_collide*xi_168 - xi_171 - xi_181 - xia_16_collide - xia_18_collide);
      const float xi_240 = u_1_a_collide*xi_239;
      const float xi_243 = xi_240 + xi_242 + xia_3_collide - xia_7_collide;
      const float xi_245 = xi_243*xi_244;
      const float xi_248 = xi_240 + xi_247 - xia_14_collide + xia_20_collide;
      const float xi_249 = xi_244*xi_248;
      const float xi_279 = rho_a_collide*u_1_a_collide;
      const float xi_281 = xi_280*(u_0_a_collide*xi_279 + xi_177 + xi_241 - xia_11_collide);
      const float xi_282 = -xi_281;
      const float xi_293 = xi_243*xi_276;
      const float xi_304 = xi_248*xi_276;
      const float xi_314 = -xi_304;
      const float u_2_a_collide = xia_24_collide;
      const float xi_2 = u_2_a_collide*xia_26_collide;
      const float xi_9 = xi_2*0.16666666666666666f;
      const float xi_10 = xi_2*0.083333333333333329f;
      const float xi_11 = omega_even_a*xi_10;
      const float xi_12 = xi_11 - xi_9;
      const float xi_15 = xi_12 + xi_14;
      const float xi_16 = omega_even_a*xi_8 + omega_shear_a*xi_1*-0.5f + xi_15 + xi_8;
      const float xi_25 = omega_even_a*xi_20 + omega_shear_a*xi_0*-0.5f + xi_12 + xi_20 + xi_24;
      const float xi_29 = xi_2*0.33333333333333331f;
      const float xi_31 = omega_even_a*xi_29 + omega_shear_a*xi_2*-0.5f + xi_29 + xi_30;
      const float xi_36 = omega_even_a*u_2_a_collide*xia_26_collide*-0.041666666666666664f;
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
      const float xi_75 = u_2_a_collide*xia_10_collide;
      const float xi_76 = xi_75*0.25f;
      const float xi_78 = xi_45*xi_75;
      const float xi_79 = xi_74 + xi_76 - xi_77 - xi_78;
      const float xi_80 = -xi_74 - xi_76 + xi_77 + xi_78;
      const float xi_81 = xi_15 + xi_22 + xi_51 + xi_71;
      const float xi_164 = u_2_a_collide*u_2_a_collide;
      const float xi_165 = rho_a_collide*xi_164*-0.33333333333333331f;
      const float xi_172 = omega_even_a*(rho_a_collide*xi_166*-0.16666666666666666f - xi_165 - xi_169 - xi_171);
      const float xi_176 = omega_even_a*(rho_a_collide*xi_168*-0.16666666666666666f - xi_165 - xi_167 - xi_173 - xi_175);
      const float xi_179 = omega_even_a*(rho_a_collide*xi_164*-0.16666666666666666f - xi_169 - xi_173 - xi_178);
      const float xi_184 = omega_shear_a*(rho_a_collide*xi_164 - xi_171 - xi_175 - xi_180 - xia_17_collide - xia_9_collide);
      const float xi_253 = xi_179*-0.5f;
      const float xi_254 = xi_172*-0.5f;
      const float xi_255 = xi_183*0.5f + xi_253 + xi_254;
      const float xi_265 = xi_176*-0.5f;
      const float xi_266 = xi_182*0.5f + xi_253 + xi_265;
      const float xi_268 = u_2_a_collide*xi_239;
      const float xi_269 = xi_247 + xi_268 + xia_14_collide - xia_20_collide;
      const float xi_270 = xi_244*xi_269;
      const float xi_271 = xi_261 + xi_268 + xia_19_collide - xia_25_collide;
      const float xi_272 = xi_244*xi_271;
      const float xi_274 = xi_184*0.5f + xi_254 + xi_265;
      const float xi_292 = xi_179*0.25f;
      const float xi_294 = xi_292 + xi_293;
      const float xi_302 = xi_292 - xi_293;
      const float xi_305 = xi_280*(u_2_a_collide*xi_279 + xi_170 + xi_246 - xia_23_collide);
      const float xi_311 = xi_172*0.25f;
      const float xi_312 = xi_269*xi_276;
      const float xi_313 = xi_311 + xi_312;
      const float xi_315 = -xi_305;
      const float xi_322 = xi_280*(rho_a_collide*u_0_a_collide*u_2_a_collide + xi_174 + xi_260 - xia_1_collide);
      const float xi_323 = -xi_322;
      const float xi_327 = xi_176*0.25f;
      const float xi_328 = xi_271*xi_276;
      const float xi_329 = xi_327 + xi_328;
      const float xi_334 = xi_311 - xi_312;
      const float xi_337 = xi_327 - xi_328;
      const float forceTerm_0_a_collide = omega_shear_a*u_0_a_collide*xia_10_collide + omega_shear_a*u_1_a_collide*xia_4_collide + omega_shear_a*u_2_a_collide*xia_26_collide - xi_0*xi_3 - xi_0 - xi_1*xi_3 - xi_1 - xi_2*xi_3 - xi_2;
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
      const float xib_1_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3];
      const float xi_362 = -xib_1_collide;
      const float xib_2_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3];
      const float xi_363 = xi_362 + xib_2_collide;
      const float xib_3_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3];
      const float xib_4_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2];
      const float xib_5_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3];
      const float xib_6_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3];
      const float xi_367 = -xib_6_collide;
      const float xi_368 = xi_367 + xib_5_collide;
      const float xib_7_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3];
      const float xib_8_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3];
      const float xib_9_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3];
      const float xi_86 = xib_9_collide*0.16666666666666666f;
      const float xi_89 = omega_odd_b*xi_86;
      const float xi_119 = xib_9_collide*0.083333333333333329f;
      const float xi_120 = xi_115*xib_9_collide;
      const float xi_121 = -xi_119 + xi_120;
      const float xi_123 = xib_9_collide*0.25f;
      const float xi_128 = xi_127*xib_9_collide;
      const float xi_135 = xi_119 - xi_120;
      const float xib_10_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3];
      const float xib_11_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3];
      const float xi_108 = xib_11_collide*0.16666666666666666f;
      const float xi_110 = omega_odd_b*xi_108;
      const float xi_137 = xib_11_collide*0.083333333333333329f;
      const float xi_138 = xi_115*xib_11_collide;
      const float xi_139 = xi_137 - xi_138;
      const float xi_152 = -xi_137 + xi_138;
      const float xib_12_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3];
      const float xib_13_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3];
      const float xib_14_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
      const float xib_15_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3];
      const float xib_16_collide = _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float xib_17_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2];
      const float xi_99 = xib_17_collide*0.16666666666666666f;
      const float xi_100 = omega_odd_b*xi_99;
      const float xi_114 = xib_17_collide*0.083333333333333329f;
      const float xi_116 = xi_115*xib_17_collide;
      const float xi_117 = xi_114 - xi_116;
      const float xi_133 = -xi_114 + xi_116;
      const float xib_18_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3];
      const float xi_349 = xib_12_collide + xib_18_collide;
      const float xib_19_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
      const float xib_20_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3];
      const float xib_21_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3];
      const float xi_345 = xib_10_collide + xib_21_collide;
      const float xi_346 = xi_345 + xib_1_collide + xib_2_collide;
      const float xib_22_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
      const float xib_23_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3];
      const float xib_24_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3];
      const float xi_352 = xib_24_collide + xib_3_collide;
      const float xi_353 = xi_352 + xib_5_collide + xib_6_collide;
      const float xib_25_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3];
      const float xi_350 = xi_349 + xib_25_collide + xib_8_collide;
      const float xi_375 = -xib_25_collide;
      const float xi_376 = xi_375 + xib_8_collide;
      const float xib_26_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3];
      const float rr_0_b_collide = 0.0f;
      const float xi_87 = rr_0_b_collide*0.25f;
      const float xi_88 = xi_87*xib_9_collide;
      const float xi_101 = xi_87*xib_17_collide;
      const float xi_109 = xi_87*xib_11_collide;
      const float rho_b_collide = xib_16_collide;
      const float xi_342 = rho_b_collide*-0.092592592592592601f;
      const float xi_355 = rho_b_collide*-0.27777777777777773f;
      const float xi_357 = xi_346 + xi_355;
      const float xi_360 = rho_b_collide*0.33333333333333331f;
      const float u_0_b_collide = xib_14_collide;
      const float xi_82 = u_0_b_collide*xib_17_collide;
      const float xi_95 = xi_82*0.083333333333333329f;
      const float xi_96 = omega_even_b*xi_95 + xi_82*-0.16666666666666666f;
      const float xi_102 = xi_82*0.33333333333333331f;
      const float xi_124 = u_0_b_collide*xi_123;
      const float xi_129 = u_0_b_collide*xi_128;
      const float xi_140 = -xi_95;
      const float xi_143 = omega_even_b*xi_82*0.041666666666666664f;
      const float xi_155 = u_0_b_collide*xib_11_collide;
      const float xi_156 = xi_155*0.25f;
      const float xi_159 = xi_127*xi_155;
      const float xi_341 = u_0_b_collide*u_0_b_collide;
      const float xi_348 = rho_b_collide*xi_341*-0.33333333333333331f;
      const float xi_356 = omega_shear_b*(rho_b_collide*xi_341 - xi_350 - xi_353 - xi_355 - xib_13_collide - xib_23_collide);
      const float xi_374 = u_0_b_collide*xi_360;
      const float xi_377 = xi_374 + xi_376 + xib_12_collide - xib_18_collide;
      const float xi_378 = xi_365*xi_377;
      const float xi_379 = xi_368 + xi_374 + xib_24_collide - xib_3_collide;
      const float xi_380 = xi_365*xi_379;
      const float xi_390 = xi_379*xi_389;
      const float xi_391 = -xi_390;
      const float xi_407 = xi_377*xi_389;
      const float xi_408 = -xi_407;
      const float u_1_b_collide = xib_19_collide;
      const float xi_83 = u_1_b_collide*xib_9_collide;
      const float xi_90 = xi_83*0.33333333333333331f;
      const float xi_103 = xi_83*0.16666666666666666f;
      const float xi_104 = xi_83*0.083333333333333329f;
      const float xi_105 = omega_even_b*xi_104;
      const float xi_106 = -xi_103 + xi_105;
      const float xi_112 = xi_106 + xi_96;
      const float xi_125 = u_1_b_collide*0.25f;
      const float xi_126 = xi_125*xib_17_collide;
      const float xi_130 = u_1_b_collide*xi_127;
      const float xi_131 = xi_130*xib_17_collide;
      const float xi_132 = xi_124 + xi_126 - xi_129 - xi_131;
      const float xi_134 = -xi_124 - xi_126 + xi_129 + xi_131;
      const float xi_141 = -xi_105;
      const float xi_145 = xi_125*xib_11_collide;
      const float xi_147 = xi_130*xib_11_collide;
      const float xi_153 = omega_even_b*u_1_b_collide*xib_9_collide*-0.041666666666666664f;
      const float xi_339 = u_1_b_collide*u_1_b_collide;
      const float xi_340 = rho_b_collide*xi_339*-0.33333333333333331f;
      const float xi_358 = omega_shear_b*(rho_b_collide*xi_339 - xi_353 - xi_357 - xib_15_collide - xib_26_collide);
      const float xi_361 = u_1_b_collide*xi_360;
      const float xi_364 = xi_361 + xi_363 - xib_10_collide + xib_21_collide;
      const float xi_366 = xi_364*xi_365;
      const float xi_369 = xi_361 + xi_368 - xib_24_collide + xib_3_collide;
      const float xi_370 = xi_365*xi_369;
      const float xi_392 = rho_b_collide*u_1_b_collide;
      const float xi_394 = xi_393*(u_0_b_collide*xi_392 + xi_352 + xi_367 - xib_5_collide);
      const float xi_395 = -xi_394;
      const float xi_397 = xi_369*xi_389;
      const float xi_400 = xi_364*xi_389;
      const float xi_405 = -xi_400;
      const float u_2_b_collide = xib_22_collide;
      const float xi_84 = u_2_b_collide*xib_11_collide;
      const float xi_91 = xi_84*0.16666666666666666f;
      const float xi_92 = xi_84*0.083333333333333329f;
      const float xi_93 = omega_even_b*xi_92;
      const float xi_94 = -xi_91 + xi_93;
      const float xi_97 = xi_94 + xi_96;
      const float xi_98 = omega_even_b*xi_90 + omega_shear_b*xi_83*-0.5f + xi_90 + xi_97;
      const float xi_107 = omega_even_b*xi_102 + omega_shear_b*xi_82*-0.5f + xi_102 + xi_106 + xi_94;
      const float xi_111 = xi_84*0.33333333333333331f;
      const float xi_113 = omega_even_b*xi_111 + omega_shear_b*xi_84*-0.5f + xi_111 + xi_112;
      const float xi_118 = omega_even_b*u_2_b_collide*xib_11_collide*-0.041666666666666664f;
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
      const float xi_157 = u_2_b_collide*xib_17_collide;
      const float xi_158 = xi_157*0.25f;
      const float xi_160 = xi_127*xi_157;
      const float xi_161 = xi_156 + xi_158 - xi_159 - xi_160;
      const float xi_162 = -xi_156 - xi_158 + xi_159 + xi_160;
      const float xi_163 = xi_104 + xi_133 + xi_153 + xi_97;
      const float xi_343 = u_2_b_collide*u_2_b_collide;
      const float xi_344 = rho_b_collide*xi_343*-0.33333333333333331f + xi_342;
      const float xi_347 = omega_even_b*(rho_b_collide*xi_341*-0.16666666666666666f - xi_340 - xi_344 - xi_346);
      const float xi_351 = omega_even_b*(rho_b_collide*xi_339*-0.16666666666666666f - xi_344 - xi_348 - xi_350);
      const float xi_354 = omega_even_b*(rho_b_collide*xi_343*-0.16666666666666666f - xi_340 - xi_342 - xi_348 - xi_353);
      const float xi_359 = omega_shear_b*(rho_b_collide*xi_343 - xi_350 - xi_357 - xib_20_collide - xib_7_collide);
      const float xi_371 = xi_347*-0.5f;
      const float xi_372 = xi_354*-0.5f;
      const float xi_373 = xi_358*0.5f + xi_371 + xi_372;
      const float xi_381 = xi_351*-0.5f;
      const float xi_382 = xi_356*0.5f + xi_372 + xi_381;
      const float xi_383 = u_2_b_collide*xi_360;
      const float xi_384 = xi_363 + xi_383 + xib_10_collide - xib_21_collide;
      const float xi_385 = xi_365*xi_384;
      const float xi_386 = xi_376 + xi_383 - xib_12_collide + xib_18_collide;
      const float xi_387 = xi_365*xi_386;
      const float xi_388 = xi_359*0.5f + xi_371 + xi_381;
      const float xi_396 = xi_354*0.25f;
      const float xi_398 = xi_396 + xi_397;
      const float xi_399 = xi_396 - xi_397;
      const float xi_401 = xi_393*(u_2_b_collide*xi_392 + xi_345 + xi_362 - xib_2_collide);
      const float xi_402 = xi_347*0.25f;
      const float xi_403 = xi_384*xi_389;
      const float xi_404 = xi_402 + xi_403;
      const float xi_406 = -xi_401;
      const float xi_409 = xi_393*(rho_b_collide*u_0_b_collide*u_2_b_collide + xi_349 + xi_375 - xib_8_collide);
      const float xi_410 = -xi_409;
      const float xi_411 = xi_351*0.25f;
      const float xi_412 = xi_386*xi_389;
      const float xi_413 = xi_411 + xi_412;
      const float xi_414 = xi_402 - xi_403;
      const float xi_415 = xi_411 - xi_412;
      const float forceTerm_0_b_collide = omega_shear_b*u_0_b_collide*xib_17_collide + omega_shear_b*u_1_b_collide*xib_9_collide + omega_shear_b*u_2_b_collide*xib_11_collide - xi_82*xi_85 - xi_82 - xi_83*xi_85 - xi_83 - xi_84*xi_85 - xi_84;
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
      const float tmp_a0 = forceTerm_0_a_collide + xi_172 + xi_176 + xi_179 - xi_182 - xi_183 - xi_184 + xi_238 + xia_15_collide;
      const float tmp_a1 = forceTerm_1_a_collide - xi_245 - xi_249 + xi_252 + xi_255 + xia_16_collide;
      const float tmp_a2 = forceTerm_2_a_collide + xi_245 + xi_249 + xi_255 + xi_256 + xia_18_collide;
      const float tmp_a3 = forceTerm_3_a_collide + xi_259 + xi_263 + xi_264 + xi_266 + xia_8_collide;
      const float tmp_a4 = forceTerm_4_a_collide - xi_259 - xi_263 + xi_266 + xi_267 + xia_6_collide;
      const float tmp_a5 = forceTerm_5_a_collide - xi_270 - xi_272 + xi_273 + xi_274 + xia_9_collide;
      const float tmp_a6 = forceTerm_6_a_collide + xi_270 + xi_272 + xi_274 + xi_275 + xia_17_collide;
      const float tmp_a7 = forceTerm_7_a_collide + xi_278 + xi_282 + xi_291 + xi_294 + xia_7_collide;
      const float tmp_a8 = forceTerm_8_a_collide + xi_277 + xi_281 + xi_294 + xi_300 + xia_2_collide;
      const float tmp_a9 = forceTerm_9_a_collide + xi_278 + xi_281 + xi_301 + xi_302 + xia_11_collide;
      const float tmp_a10 = forceTerm_10_a_collide + xi_277 + xi_282 + xi_302 + xi_303 + xia_3_collide;
      const float tmp_a11 = forceTerm_11_a_collide + xi_304 + xi_305 + xi_310 + xi_313 + xia_13_collide;
      const float tmp_a12 = forceTerm_12_a_collide + xi_313 + xi_314 + xi_315 + xi_319 + xia_20_collide;
      const float tmp_a13 = forceTerm_13_a_collide + xi_321 + xi_323 + xi_326 + xi_329 + xia_25_collide;
      const float tmp_a14 = forceTerm_14_a_collide + xi_320 + xi_322 + xi_329 + xi_332 + xia_5_collide;
      const float tmp_a15 = forceTerm_15_a_collide + xi_304 + xi_315 + xi_333 + xi_334 + xia_14_collide;
      const float tmp_a16 = forceTerm_16_a_collide + xi_305 + xi_314 + xi_334 + xi_335 + xia_23_collide;
      const float tmp_a17 = forceTerm_17_a_collide + xi_321 + xi_322 + xi_336 + xi_337 + xia_1_collide;
      const float tmp_a18 = forceTerm_18_a_collide + xi_320 + xi_323 + xi_337 + xi_338 + xia_19_collide;
      const float tmp_b0 = forceTerm_0_b_collide + xi_238 + xi_347 + xi_351 + xi_354 - xi_356 - xi_358 - xi_359 + xib_4_collide;
      const float tmp_b1 = forceTerm_1_b_collide + xi_252 - xi_366 - xi_370 + xi_373 + xib_26_collide;
      const float tmp_b2 = forceTerm_2_b_collide + xi_256 + xi_366 + xi_370 + xi_373 + xib_15_collide;
      const float tmp_b3 = forceTerm_3_b_collide + xi_264 + xi_378 + xi_380 + xi_382 + xib_23_collide;
      const float tmp_b4 = forceTerm_4_b_collide + xi_267 - xi_378 - xi_380 + xi_382 + xib_13_collide;
      const float tmp_b5 = forceTerm_5_b_collide + xi_273 - xi_385 - xi_387 + xi_388 + xib_20_collide;
      const float tmp_b6 = forceTerm_6_b_collide + xi_275 + xi_385 + xi_387 + xi_388 + xib_7_collide;
      const float tmp_b7 = forceTerm_7_b_collide + xi_291 + xi_391 + xi_395 + xi_398 + xib_24_collide;
      const float tmp_b8 = forceTerm_8_b_collide + xi_300 + xi_390 + xi_394 + xi_398 + xib_6_collide;
      const float tmp_b9 = forceTerm_9_b_collide + xi_301 + xi_391 + xi_394 + xi_399 + xib_5_collide;
      const float tmp_b10 = forceTerm_10_b_collide + xi_303 + xi_390 + xi_395 + xi_399 + xib_3_collide;
      const float tmp_b11 = forceTerm_11_b_collide + xi_310 + xi_400 + xi_401 + xi_404 + xib_1_collide;
      const float tmp_b12 = forceTerm_12_b_collide + xi_319 + xi_404 + xi_405 + xi_406 + xib_21_collide;
      const float tmp_b13 = forceTerm_13_b_collide + xi_326 + xi_408 + xi_410 + xi_413 + xib_12_collide;
      const float tmp_b14 = forceTerm_14_b_collide + xi_332 + xi_407 + xi_409 + xi_413 + xib_25_collide;
      const float tmp_b15 = forceTerm_15_b_collide + xi_333 + xi_400 + xi_406 + xi_414 + xib_10_collide;
      const float tmp_b16 = forceTerm_16_b_collide + xi_335 + xi_401 + xi_405 + xi_414 + xib_2_collide;
      const float tmp_b17 = forceTerm_17_b_collide + xi_336 + xi_408 + xi_409 + xi_415 + xib_8_collide;
      const float tmp_b18 = forceTerm_18_b_collide + xi_338 + xi_407 + xi_410 + xi_415 + xib_18_collide;
      const float xirecolor_0 = tmp_a0 + tmp_b0;
      const float xirecolor_1 = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float xirecolor_2 = ((1.0f) / (xirecolor_1));
      const float xi_437 = xirecolor_2*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float xirecolor_3 = xirecolor_2*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
      const float xirecolor_4 = tmp_a1 + tmp_b1;
      const float xirecolor_5 = xi_195;
      const float xirecolor_6 = xi_193;
      const float xirecolor_7 = xi_213;
      const float xi_416 = -xirecolor_7;
      const float xirecolor_8 = xi_416;
      const float xirecolor_9 = xi_200;
      const float xirecolor_10 = -xirecolor_9;
      const float xirecolor_11 = xi_202;
      const float xi_417 = xirecolor_10 + xirecolor_11;
      const float xirecolor_12 = xi_215;
      const float xi_418 = xirecolor_12 + xirecolor_8;
      const float xirecolor_13 = xi_417 + xi_418;
      const float xirecolor_14 = xi_217;
      const float xirecolor_15 = xi_218;
      const float xi_419 = xirecolor_14 - xirecolor_15;
      const float xirecolor_16 = -xi_419;
      const float xirecolor_17 = xi_187;
      const float xirecolor_18 = xi_188;
      const float xirecolor_19 = xi_189;
      const float xi_420 = xirecolor_18 + xirecolor_19;
      const float xirecolor_20 = xi_190;
      const float xirecolor_21 = xi_420 - xirecolor_17 + xirecolor_20;
      const float xirecolor_22 = xi_186;
      const float xirecolor_23 = xi_185;
      const float xi_421 = xirecolor_22 - xirecolor_23;
      const float xirecolor_24 = -xi_421;
      const float xirecolor_25 = xi_221;
      const float xirecolor_26 = -xirecolor_25;
      const float xirecolor_27 = xi_207;
      const float xirecolor_28 = xi_206;
      const float xirecolor_29 = xi_223;
      const float xi_422 = xirecolor_26 + xirecolor_29;
      const float xirecolor_30 = xi_422 - xirecolor_27 + xirecolor_28;
      const float xirecolor_31 = xirecolor_13 + xirecolor_16 + xirecolor_21 + xirecolor_24 + xirecolor_30 - xirecolor_5 + xirecolor_6;
      const float xirecolor_32 = xi_192;
      const float xirecolor_33 = xi_194;
      const float xi_423 = xirecolor_32 + xirecolor_33;
      const float xirecolor_34 = xi_423 + xirecolor_5 - xirecolor_6;
      const float xirecolor_35 = xi_197;
      const float xirecolor_36 = xi_198;
      const float xi_424 = xirecolor_35 - xirecolor_36;
      const float xirecolor_37 = -xi_424;
      const float xirecolor_38 = xi_417;
      const float xirecolor_39 = xi_204;
      const float xirecolor_40 = xi_205;
      const float xi_425 = xirecolor_39 + xirecolor_40;
      const float xirecolor_41 = xirecolor_27 - xirecolor_28;
      const float xirecolor_42 = xi_425 + xirecolor_41;
      const float xirecolor_43 = xirecolor_38 + xirecolor_42;
      const float xirecolor_44 = xi_421 + xirecolor_21 + xirecolor_34 + xirecolor_37 + xirecolor_43;
      const float xirecolor_45 = xi_419;
      const float xirecolor_46 = xi_418;
      const float xirecolor_47 = xi_229;
      const float xirecolor_48 = -xirecolor_47;
      const float xirecolor_49 = xi_231;
      const float xi_426 = xirecolor_48 + xirecolor_49;
      const float xirecolor_50 = xi_424 + xi_426;
      const float xirecolor_51 = xi_420 + xirecolor_17 - xirecolor_20 + xirecolor_24 + xirecolor_34 + xirecolor_45 + xirecolor_46 + xirecolor_50;
      const float xirecolor_52 = powf(xirecolor_31*xirecolor_31 + xirecolor_44*xirecolor_44 + xirecolor_51*xirecolor_51, 0.5f);
      const float xirecolor_53 = ((1.0f) / (xirecolor_52));
      const float xi_427 = xirecolor_31*xirecolor_53;
      const float xi_428 = xirecolor_51*xirecolor_53;
      const float xi_429 = xirecolor_44*xirecolor_53;
      const bool xirecolor_54 = xirecolor_52 > 0.0f;
      const float xirecolor_55 = 0.69999999999999996f*((1.0f) / (xirecolor_1*xirecolor_1))*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2]*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float xirecolor_56 = xirecolor_55*(0.046296296296296259f*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.046296296296296259f*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
      const float xirecolor_57 = xirecolor_56*((xirecolor_54) ? (xi_427): (0.0f));
      const float xirecolor_58 = tmp_a2 + tmp_b2;
      const float xirecolor_59 = xirecolor_55*(0.046296296296296266f*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.046296296296296266f*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
      const float xirecolor_60 = xirecolor_59*((xirecolor_54) ? (-xi_427): (0.0f));
      const float xirecolor_61 = tmp_a3 + tmp_b3;
      const float xirecolor_62 = xirecolor_59*((xirecolor_54) ? (-xi_428): (0.0f));
      const float xirecolor_63 = tmp_a4 + tmp_b4;
      const float xirecolor_64 = xirecolor_56*((xirecolor_54) ? (xi_428): (0.0f));
      const float xirecolor_65 = tmp_a5 + tmp_b5;
      const float xirecolor_66 = xirecolor_56*((xirecolor_54) ? (-xi_429): (0.0f));
      const float xirecolor_67 = tmp_a6 + tmp_b6;
      const float xirecolor_68 = xirecolor_59*((xirecolor_54) ? (xi_429): (0.0f));
      const float xirecolor_69 = tmp_a7 + tmp_b7;
      const float xirecolor_70 = xi_283;
      const float xirecolor_71 = xi_284;
      const float xirecolor_72 = xi_423;
      const float xirecolor_73 = xirecolor_50 + xirecolor_72;
      const float xirecolor_74 = xi_288;
      const float xirecolor_75 = xirecolor_25 - xirecolor_29 + xirecolor_74;
      const float xirecolor_76 = xi_285 - xirecolor_11 + xirecolor_41 - xirecolor_70 + xirecolor_71 + xirecolor_73 + xirecolor_75 + xirecolor_9;
      const float xirecolor_77 = xirecolor_53*0.70710678118654757f;
      const float xi_430 = xirecolor_76*xirecolor_77;
      const float xirecolor_78 = xirecolor_55*(0.02314814814814815f*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.02314814814814815f*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
      const float xirecolor_79 = xirecolor_78*((xirecolor_54) ? (-xi_430): (0.0f));
      const float xirecolor_80 = tmp_a8 + tmp_b8;
      const float xirecolor_81 = xi_295;
      const float xirecolor_82 = xi_296;
      const float xi_431 = -xirecolor_81 + xirecolor_82;
      const float xirecolor_83 = xi_297;
      const float xirecolor_84 = xi_298 + xi_431 + xirecolor_30 + xirecolor_38 + xirecolor_73 + xirecolor_83;
      const float xi_432 = xirecolor_77*xirecolor_84;
      const float xirecolor_85 = xirecolor_78*((xirecolor_54) ? (xi_432): (0.0f));
      const float xirecolor_86 = tmp_a9 + tmp_b9;
      const float xirecolor_87 = xirecolor_78*((xirecolor_54) ? (-xi_432): (0.0f));
      const float xirecolor_88 = tmp_a10 + tmp_b10;
      const float xirecolor_89 = xirecolor_78*((xirecolor_54) ? (xi_430): (0.0f));
      const float xirecolor_90 = tmp_a11 + tmp_b11;
      const float xirecolor_91 = xi_425 + xirecolor_37 + xirecolor_72;
      const float xirecolor_92 = -xi_416 - xi_431 - xirecolor_12;
      const float xirecolor_93 = xi_306 + xirecolor_45 + xirecolor_75 + xirecolor_91 + xirecolor_92;
      const float xi_433 = xirecolor_77*xirecolor_93;
      const float xirecolor_94 = xirecolor_78*((xirecolor_54) ? (-xi_433): (0.0f));
      const float xirecolor_95 = tmp_a12 + tmp_b12;
      const float xirecolor_96 = xirecolor_16 + xirecolor_70 - xirecolor_71;
      const float xirecolor_97 = xi_316 + xi_422 + xirecolor_46 + xirecolor_83 + xirecolor_91 + xirecolor_96;
      const float xi_434 = xirecolor_77*xirecolor_97;
      const float xirecolor_98 = xirecolor_78*((xirecolor_54) ? (-xi_434): (0.0f));
      const float xirecolor_99 = tmp_a13 + tmp_b13;
      const float xirecolor_100 = xi_324 + xi_426 + xirecolor_13 + xirecolor_42 + xirecolor_45 + xirecolor_74 + xirecolor_83;
      const float xi_435 = xirecolor_100*xirecolor_77;
      const float xirecolor_101 = xirecolor_78*((xirecolor_54) ? (-xi_435): (0.0f));
      const float xirecolor_102 = tmp_a14 + tmp_b14;
      const float xirecolor_103 = xi_330 + xirecolor_43 + xirecolor_47 - xirecolor_49 + xirecolor_92 + xirecolor_96;
      const float xi_436 = xirecolor_103*xirecolor_77;
      const float xirecolor_104 = xirecolor_78*((xirecolor_54) ? (-xi_436): (0.0f));
      const float xirecolor_105 = tmp_a15 + tmp_b15;
      const float xirecolor_106 = xirecolor_78*((xirecolor_54) ? (xi_434): (0.0f));
      const float xirecolor_107 = tmp_a16 + tmp_b16;
      const float xirecolor_108 = xirecolor_78*((xirecolor_54) ? (xi_433): (0.0f));
      const float xirecolor_109 = tmp_a17 + tmp_b17;
      const float xirecolor_110 = xirecolor_78*((xirecolor_54) ? (xi_436): (0.0f));
      const float xirecolor_111 = tmp_a18 + tmp_b18;
      const float xirecolor_112 = xirecolor_78*((xirecolor_54) ? (xi_435): (0.0f));
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


void ColorGradientCollideSweepSinglePrecisionCUDA::run(IBlock * block, gpuStream_t stream)
{
   
    auto pdfs_a = block->getData< gpu::GPUField<float> >(pdfs_aID);
    auto velocity = block->getData< gpu::GPUField<float> >(velocityID);
    auto phasefield = block->getData< gpu::GPUField<float> >(phasefieldID);
    auto force_b = block->getData< gpu::GPUField<float> >(force_bID);
    auto pdfs_b = block->getData< gpu::GPUField<float> >(pdfs_bID);
    auto rho_a = block->getData< gpu::GPUField<float> >(rho_aID);
    auto force_a = block->getData< gpu::GPUField<float> >(force_aID);
    auto rho_b = block->getData< gpu::GPUField<float> >(rho_bID);

    auto & omega_even_b = this->omega_even_b_;
    auto & omega_shear_a = this->omega_shear_a_;
    auto & omega_odd_a = this->omega_odd_a_;
    auto & omega_odd_b = this->omega_odd_b_;
    auto & omega_shear_b = this->omega_shear_b_;
    auto & omega_even_a = this->omega_even_a_;
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_a->nrOfGhostLayers()))
    float * RESTRICT const _data_force_a = force_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_b->nrOfGhostLayers()))
    float * RESTRICT const _data_force_b = force_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_a->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_b->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phasefield->nrOfGhostLayers()))
    float * RESTRICT const _data_phasefield = phasefield->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_a->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_a = rho_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_b->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_b = rho_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
    float * RESTRICT const _data_velocity = velocity->dataAt(-1, -1, -1, 0);
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
    dim3 _block(uint32_c(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)), uint32_c(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))), uint32_c(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))));
    dim3 _grid(uint32_c(( (_size_force_a_0 - 2) % (((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) == 0 ? (int64_t)(_size_force_a_0 - 2) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) : ( (int64_t)(_size_force_a_0 - 2) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) ) +1 )), uint32_c(( (_size_force_a_1 - 2) % (((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) == 0 ? (int64_t)(_size_force_a_1 - 2) / (int64_t)(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) : ( (int64_t)(_size_force_a_1 - 2) / (int64_t)(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) ) +1 )), uint32_c(( (_size_force_a_2 - 2) % (((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) == 0 ? (int64_t)(_size_force_a_2 - 2) / (int64_t)(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) : ( (int64_t)(_size_force_a_2 - 2) / (int64_t)(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) ) +1 )));
    internal_colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda::colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, omega_even_a, omega_even_b, omega_odd_a, omega_odd_b, omega_shear_a, omega_shear_b);
    
}


void ColorGradientCollideSweepSinglePrecisionCUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block, gpuStream_t stream)
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

    auto & omega_even_b = this->omega_even_b_;
    auto & omega_shear_a = this->omega_shear_a_;
    auto & omega_odd_a = this->omega_odd_a_;
    auto & omega_odd_b = this->omega_odd_b_;
    auto & omega_shear_b = this->omega_shear_b_;
    auto & omega_even_a = this->omega_even_a_;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    float * RESTRICT const _data_force_a = force_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    float * RESTRICT const _data_force_b = force_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    float * RESTRICT const _data_phasefield = phasefield->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_a = rho_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_b = rho_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    float * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
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
    dim3 _block(uint32_c(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)), uint32_c(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))), uint32_c(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))));
    dim3 _grid(uint32_c(( (_size_force_a_0 - 2) % (((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) == 0 ? (int64_t)(_size_force_a_0 - 2) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) : ( (int64_t)(_size_force_a_0 - 2) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) ) +1 )), uint32_c(( (_size_force_a_1 - 2) % (((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) == 0 ? (int64_t)(_size_force_a_1 - 2) / (int64_t)(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) : ( (int64_t)(_size_force_a_1 - 2) / (int64_t)(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) ) +1 )), uint32_c(( (_size_force_a_2 - 2) % (((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) == 0 ? (int64_t)(_size_force_a_2 - 2) / (int64_t)(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) : ( (int64_t)(_size_force_a_2 - 2) / (int64_t)(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) ) +1 )));
    internal_colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda::colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, omega_even_a, omega_even_b, omega_odd_a, omega_odd_b, omega_shear_a, omega_shear_b);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
