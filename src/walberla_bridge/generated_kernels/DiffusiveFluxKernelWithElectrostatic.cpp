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
//! \\file DiffusiveFluxKernelWithElectrostatic.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "DiffusiveFluxKernelWithElectrostatic.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#define FUNC_PREFIX

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning push
#pragma warning(disable : 1599)
#endif

using namespace std;

namespace walberla {
namespace pystencils {

namespace internal_diffusivefluxkernelwithelectrostatic_diffusivefluxkernelwithelectrostatic {
static FUNC_PREFIX void
diffusivefluxkernelwithelectrostatic_diffusivefluxkernelwithelectrostatic(
    double D, double *RESTRICT _data_j, double *RESTRICT const _data_phi,
    double *RESTRICT const _data_rho, int64_t const _size_j_0,
    int64_t const _size_j_1, int64_t const _size_j_2, int64_t const _stride_j_0,
    int64_t const _stride_j_1, int64_t const _stride_j_2,
    int64_t const _stride_j_3, int64_t const _stride_phi_0,
    int64_t const _stride_phi_1, int64_t const _stride_phi_2,
    int64_t const _stride_rho_0, int64_t const _stride_rho_1,
    int64_t const _stride_rho_2, double f_ext0, double f_ext1, double f_ext2,
    double kT, double z) {
#pragma omp parallel
  {
    {
      {
        {
          if (_size_j_1 > 1 && _size_j_2 > 1) {
            double *RESTRICT _data_j_20_312 = _data_j + 12 * _stride_j_3;
            double *RESTRICT _data_j_20_312_10 = _data_j_20_312;
            double *RESTRICT _data_rho_20_30 = _data_rho;
            double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
            double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
            double *RESTRICT _data_rho_21_30_11 =
                _stride_rho_1 + _data_rho_21_30;
            double *RESTRICT _data_phi_20_30 = _data_phi;
            double *RESTRICT _data_phi_20_30_11 =
                _stride_phi_1 + _data_phi_20_30;
            double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
            double *RESTRICT _data_phi_21_30_11 =
                _stride_phi_1 + _data_phi_21_30;
            double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
            double *RESTRICT _data_phi_21_30_10 = _data_phi_21_30;
            double *RESTRICT _data_rho_20_30_11 =
                _stride_rho_1 + _data_rho_20_30;
            double *RESTRICT _data_rho_21_30_10 = _data_rho_21_30;
            _data_j_20_312_10[_stride_j_0] =
                D *
                (f_ext0 * z *
                     (_data_rho_20_30_10[_stride_rho_0] +
                      _data_rho_21_30_11[0]) *
                     2.0 +
                 f_ext1 * z *
                     (_data_rho_20_30_10[_stride_rho_0] +
                      _data_rho_21_30_11[0]) *
                     -2.0 +
                 f_ext2 * z *
                     (_data_rho_20_30_10[_stride_rho_0] +
                      _data_rho_21_30_11[0]) *
                     -2.0 +
                 kT *
                     (-3.0 * _data_rho_20_30_10[_stride_rho_0] +
                      3.0 * _data_rho_21_30_11[0] + _data_rho_20_30_10[0] -
                      _data_rho_20_30_11[0] +
                      _data_rho_20_30_11[_stride_rho_0] -
                      _data_rho_21_30_10[0] +
                      _data_rho_21_30_10[_stride_rho_0] -
                      _data_rho_21_30_11[_stride_rho_0]) *
                     -2.0 +
                 z *
                     (_data_rho_20_30_10[_stride_rho_0] +
                      _data_rho_21_30_11[0]) *
                     (-_data_phi_20_30_10[0] +
                      _data_phi_20_30_10[_stride_phi_0] -
                      _data_phi_21_30_11[0] +
                      _data_phi_21_30_11[_stride_phi_0]) +
                 z *
                     (_data_rho_20_30_10[_stride_rho_0] +
                      _data_rho_21_30_11[0]) *
                     (_data_phi_20_30_10[_stride_phi_0] +
                      _data_phi_20_30_11[0] -
                      _data_phi_21_30_10[_stride_phi_0] -
                      _data_phi_21_30_11[0]) +
                 z *
                     (_data_rho_20_30_10[_stride_rho_0] +
                      _data_rho_21_30_11[0]) *
                     (_data_phi_20_30_10[_stride_phi_0] -
                      _data_phi_20_30_11[_stride_phi_0] +
                      _data_phi_21_30_10[0] - _data_phi_21_30_11[0])) *
                0.0833333333333333 * 1.7320508075688772 /
                (kT * (1.0 + 2.8284271247461903 +
                       1.7320508075688772 * 1.33333333333333));
          }
#pragma omp for schedule(static)
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (_size_j_1 > 1 && _size_j_2 > 1) {
              double *RESTRICT _data_j_20_312 = _data_j + 12 * _stride_j_3;
              double *RESTRICT _data_j_20_312_10 = _data_j_20_312;
              double *RESTRICT _data_rho_20_30 = _data_rho;
              double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
              double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
              double *RESTRICT _data_rho_21_30_11 =
                  _stride_rho_1 + _data_rho_21_30;
              double *RESTRICT _data_phi_20_30 = _data_phi;
              double *RESTRICT _data_phi_20_30_11 =
                  _stride_phi_1 + _data_phi_20_30;
              double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
              double *RESTRICT _data_phi_21_30_11 =
                  _stride_phi_1 + _data_phi_21_30;
              double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
              double *RESTRICT _data_phi_21_30_10 = _data_phi_21_30;
              double *RESTRICT _data_rho_20_30_11 =
                  _stride_rho_1 + _data_rho_20_30;
              double *RESTRICT _data_rho_21_30_10 = _data_rho_21_30;
              _data_j_20_312_10[_stride_j_0 * ctr_0] =
                  D *
                  (f_ext0 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0]) *
                       2.0 +
                   f_ext1 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0]) *
                       -2.0 +
                   f_ext2 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0]) *
                       -2.0 +
                   kT *
                       (-3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        3.0 * _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                                 _stride_rho_0] +
                        _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0] -
                        _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0] +
                        _data_rho_20_30_11[_stride_rho_0 * ctr_0] -
                        _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0] +
                        _data_rho_21_30_10[_stride_rho_0 * ctr_0] -
                        _data_rho_21_30_11[_stride_rho_0 * ctr_0]) *
                       -2.0 +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0]) *
                       (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                            _stride_phi_0] +
                        _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                        _data_phi_21_30_11[_stride_phi_0 * ctr_0 -
                                           _stride_phi_0] +
                        _data_phi_21_30_11[_stride_phi_0 * ctr_0]) +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0]) *
                       (_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                        _data_phi_20_30_11[_stride_phi_0 * ctr_0 -
                                           _stride_phi_0] -
                        _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                        _data_phi_21_30_11[_stride_phi_0 * ctr_0 -
                                           _stride_phi_0]) +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0]) *
                       (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                        _data_phi_20_30_11[_stride_phi_0 * ctr_0] +
                        _data_phi_21_30_10[_stride_phi_0 * ctr_0 -
                                           _stride_phi_0] -
                        _data_phi_21_30_11[_stride_phi_0 * ctr_0 -
                                           _stride_phi_0])) *
                  0.0833333333333333 * 1.7320508075688772 /
                  (kT * (1.0 + 2.8284271247461903 +
                         1.7320508075688772 * 1.33333333333333));
            }
          }
          if (_size_j_1 > 1 && _size_j_2 > 1) {
            double *RESTRICT _data_j_20_312 = _data_j + 12 * _stride_j_3;
            double *RESTRICT _data_j_20_312_10 = _data_j_20_312;
            double *RESTRICT _data_rho_20_30 = _data_rho;
            double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
            double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
            double *RESTRICT _data_rho_21_30_11 =
                _stride_rho_1 + _data_rho_21_30;
            double *RESTRICT _data_phi_20_30 = _data_phi;
            double *RESTRICT _data_phi_20_30_11 =
                _stride_phi_1 + _data_phi_20_30;
            double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
            double *RESTRICT _data_phi_21_30_11 =
                _stride_phi_1 + _data_phi_21_30;
            double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
            double *RESTRICT _data_phi_21_30_10 = _data_phi_21_30;
            double *RESTRICT _data_rho_20_30_11 =
                _stride_rho_1 + _data_rho_20_30;
            double *RESTRICT _data_rho_21_30_10 = _data_rho_21_30;
            _data_j_20_312_10[_stride_j_0 * (_size_j_0 - 1)] =
                D *
                (f_ext0 * z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0]) *
                     2.0 +
                 f_ext1 * z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0]) *
                     -2.0 +
                 f_ext2 * z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0]) *
                     -2.0 +
                 kT *
                     (-3.0 *
                          _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      3.0 * _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0] +
                      _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0] -
                      _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0] +
                      _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1)] -
                      _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0] +
                      _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                      _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1)]) *
                     -2.0 +
                 z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0]) *
                     (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                          _stride_phi_0] +
                      _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                      _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                         _stride_phi_0] +
                      _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1)]) +
                 z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0]) *
                     (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                      _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                         _stride_phi_0] -
                      _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                      _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                         _stride_phi_0]) +
                 z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0]) *
                     (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                      _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1)] +
                      _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                         _stride_phi_0] -
                      _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                         _stride_phi_0])) *
                0.0833333333333333 * 1.7320508075688772 /
                (kT * (1.0 + 2.8284271247461903 +
                       1.7320508075688772 * 1.33333333333333));
          }
        }
#pragma omp for schedule(static)
        for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
          {
            {
              if (_size_j_2 > 1 && ctr_1 > 0 && _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_36 = _data_j + 6 * _stride_j_3;
                double *RESTRICT _data_j_20_36_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_36;
                double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_21_30;
                double *RESTRICT _data_rho_20_30 = _data_rho;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30 = _data_phi;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                _data_j_20_36_10[_stride_j_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_10[0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_10[0]) *
                         -2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0] -
                          _data_rho_21_30_10[0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_10[0]) *
                         (-_data_phi_20_30_10[0] +
                          _data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_21_30_10[0] +
                          _data_phi_21_30_10[_stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_10[0]) *
                         (_data_phi_20_30_10[0] +
                          _data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_21_30_10[0] -
                          _data_phi_21_30_10[_stride_phi_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_0 > 2 && _size_j_2 > 1 && ctr_1 > 0) {
                double *RESTRICT _data_j_20_38 = _data_j + 8 * _stride_j_3;
                double *RESTRICT _data_j_20_38_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_38;
                double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21_30;
                double *RESTRICT _data_rho_20_30 = _data_rho;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_phi_20_30 = _data_phi;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_21_30;
                _data_j_20_38_10[_stride_j_0] =
                    D *
                    (f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[_stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[_stride_rho_0]) *
                         -2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0] -
                          _data_rho_21_30_1m1[_stride_rho_0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[_stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_1m1[_stride_phi_0] +
                          _data_phi_21_30_10[_stride_phi_0] -
                          _data_phi_21_30_1m1[_stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[_stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_20_30_1m1[_stride_phi_0] -
                          _data_phi_21_30_10[_stride_phi_0] -
                          _data_phi_21_30_1m1[_stride_phi_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && ctr_1 > 0) {
                double *RESTRICT _data_j_20_310 = _data_j + 10 * _stride_j_3;
                double *RESTRICT _data_j_20_310_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_310;
                double *RESTRICT _data_rho_20_30 = _data_rho;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21_30;
                double *RESTRICT _data_phi_20_30 = _data_phi;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_21_30;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_21_30;
                _data_j_20_310_10[_stride_j_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[0]) *
                         -2.0 +
                     kT *
                         (-3.0 * _data_rho_21_30_1m1[0] +
                          3.0 * _data_rho_20_30_10[_stride_rho_0] -
                          _data_rho_20_30_10[0] + _data_rho_20_30_1m1[0] -
                          _data_rho_20_30_1m1[_stride_rho_0] +
                          _data_rho_21_30_10[0] -
                          _data_rho_21_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[_stride_rho_0]) *
                         2.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[0]) *
                         (-_data_phi_20_30_10[0] +
                          _data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_21_30_1m1[0] +
                          _data_phi_21_30_1m1[_stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[0]) *
                         (_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_20_30_1m1[0] -
                          _data_phi_21_30_10[_stride_phi_0] -
                          _data_phi_21_30_1m1[0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[0]) *
                         (_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_1m1[_stride_phi_0] +
                          _data_phi_21_30_10[0] - _data_phi_21_30_1m1[0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_312 = _data_j + 12 * _stride_j_3;
                double *RESTRICT _data_j_20_312_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_312;
                double *RESTRICT _data_rho_20_30 = _data_rho;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_11 =
                    _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_21_30;
                double *RESTRICT _data_phi_20_30 = _data_phi;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_11 =
                    _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_21_30;
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_21_30;
                _data_j_20_312_10[_stride_j_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_11[0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_11[0]) *
                         -2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_11[0]) *
                         -2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0] +
                          3.0 * _data_rho_21_30_11[0] + _data_rho_20_30_10[0] -
                          _data_rho_20_30_11[0] +
                          _data_rho_20_30_11[_stride_rho_0] -
                          _data_rho_21_30_10[0] +
                          _data_rho_21_30_10[_stride_rho_0] -
                          _data_rho_21_30_11[_stride_rho_0]) *
                         -2.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_11[0]) *
                         (-_data_phi_20_30_10[0] +
                          _data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_21_30_11[0] +
                          _data_phi_21_30_11[_stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_11[0]) *
                         (_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_20_30_11[0] -
                          _data_phi_21_30_10[_stride_phi_0] -
                          _data_phi_21_30_11[0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_11[0]) *
                         (_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_11[_stride_phi_0] +
                          _data_phi_21_30_10[0] - _data_phi_21_30_11[0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
            }
            for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
              if (_size_j_2 > 1 && ctr_1 > 0 && _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_36 = _data_j + 6 * _stride_j_3;
                double *RESTRICT _data_j_20_36_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_36;
                double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_21_30;
                double *RESTRICT _data_rho_20_30 = _data_rho;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30 = _data_phi;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                _data_j_20_36_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         -2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] +
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && ctr_1 > 0 && _size_j_0 - ctr_0 > 1) {
                double *RESTRICT _data_j_20_38 = _data_j + 8 * _stride_j_3;
                double *RESTRICT _data_j_20_38_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_38;
                double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21_30;
                double *RESTRICT _data_rho_20_30 = _data_rho;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_phi_20_30 = _data_phi;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_21_30;
                _data_j_20_38_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                         -2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] +
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_1m1[_stride_phi_0 * ctr_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_1m1[_stride_phi_0 * ctr_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && ctr_1 > 0) {
                double *RESTRICT _data_j_20_310 = _data_j + 10 * _stride_j_3;
                double *RESTRICT _data_j_20_310_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_310;
                double *RESTRICT _data_rho_20_30 = _data_rho;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21_30;
                double *RESTRICT _data_phi_20_30 = _data_phi;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_21_30;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_21_30;
                _data_j_20_310_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         -2.0 +
                     kT *
                         (-3.0 * _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                                     _stride_rho_0] +
                          3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0] -
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] -
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                         2.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_1m1[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_21_30_1m1[_stride_phi_0 * ctr_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] -
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_1m1[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] +
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_21_30_1m1[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_312 = _data_j + 12 * _stride_j_3;
                double *RESTRICT _data_j_20_312_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_312;
                double *RESTRICT _data_rho_20_30 = _data_rho;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_11 =
                    _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_21_30;
                double *RESTRICT _data_phi_20_30 = _data_phi;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_11 =
                    _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_21_30;
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_21_30;
                _data_j_20_312_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         -2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         -2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          3.0 * _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                                   _stride_rho_0] +
                          _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] -
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] +
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0] -
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] +
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0]) *
                         -2.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_11[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] +
                          _data_phi_21_30_11[_stride_phi_0 * ctr_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_20_30_11[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_11[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_11[_stride_phi_0 * ctr_0] +
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_21_30_11[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
            }
            {
              if (_size_j_2 > 1 && ctr_1 > 0 && _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_36 = _data_j + 6 * _stride_j_3;
                double *RESTRICT _data_j_20_36_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_36;
                double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_21_30;
                double *RESTRICT _data_rho_20_30 = _data_rho;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30 = _data_phi;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                _data_j_20_36_10[_stride_j_0 * (_size_j_0 - 1)] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         -2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] +
                          _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1)]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_21_30_10[_stride_phi_0 *
                                             (_size_j_0 - 1)])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && ctr_1 > 0) {
                double *RESTRICT _data_j_20_310 = _data_j + 10 * _stride_j_3;
                double *RESTRICT _data_j_20_310_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_310;
                double *RESTRICT _data_rho_20_30 = _data_rho;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21_30;
                double *RESTRICT _data_phi_20_30 = _data_phi;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_21_30;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_21_30;
                _data_j_20_310_10[_stride_j_0 * (_size_j_0 - 1)] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         -2.0 +
                     kT *
                         (-3.0 * _data_rho_21_30_1m1[_stride_rho_0 *
                                                         (_size_j_0 - 1) -
                                                     _stride_rho_0] +
                          3.0 * _data_rho_20_30_10[_stride_rho_0 *
                                                   (_size_j_0 - 1)] -
                          _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0] -
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] -
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 *
                                              (_size_j_0 - 1)]) *
                         2.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_21_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_21_30_1m1[_stride_phi_0 *
                                              (_size_j_0 - 1)]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] -
                          _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_21_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_21_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_312 = _data_j + 12 * _stride_j_3;
                double *RESTRICT _data_j_20_312_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_312;
                double *RESTRICT _data_rho_20_30 = _data_rho;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_11 =
                    _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_21_30;
                double *RESTRICT _data_phi_20_30 = _data_phi;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_11 =
                    _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_21_30;
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_21_30;
                _data_j_20_312_10[_stride_j_0 * (_size_j_0 - 1)] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         -2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         -2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1)] +
                          3.0 * _data_rho_21_30_11[_stride_rho_0 *
                                                       (_size_j_0 - 1) -
                                                   _stride_rho_0] +
                          _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] -
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] +
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] +
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1)]) *
                         -2.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] +
                          _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1)]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
            }
          }
        }
        {
          {
            if (_size_j_0 > 2 && _size_j_1 > 1 && _size_j_2 > 1) {
              double *RESTRICT _data_j_20_38 = _data_j + 8 * _stride_j_3;
              double *RESTRICT _data_j_20_38_10 =
                  _stride_j_1 * (_size_j_1 - 1) + _data_j_20_38;
              double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
              double *RESTRICT _data_rho_21_30_1m1 =
                  _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                  _data_rho_21_30;
              double *RESTRICT _data_rho_20_30 = _data_rho;
              double *RESTRICT _data_rho_20_30_10 =
                  _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
              double *RESTRICT _data_phi_20_30 = _data_phi;
              double *RESTRICT _data_phi_20_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_20_30;
              double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
              double *RESTRICT _data_phi_21_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_21_30;
              double *RESTRICT _data_phi_20_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
              double *RESTRICT _data_phi_21_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21_30;
              _data_j_20_38_10[_stride_j_0] =
                  D *
                  (f_ext1 * z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_21_30_1m1[_stride_rho_0]) *
                       2.0 +
                   f_ext2 * z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_21_30_1m1[_stride_rho_0]) *
                       -2.0 +
                   kT *
                       (_data_rho_20_30_10[_stride_rho_0] -
                        _data_rho_21_30_1m1[_stride_rho_0]) *
                       4.0 +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_21_30_1m1[_stride_rho_0]) *
                       (_data_phi_20_30_10[_stride_phi_0] -
                        _data_phi_20_30_1m1[_stride_phi_0] +
                        _data_phi_21_30_10[_stride_phi_0] -
                        _data_phi_21_30_1m1[_stride_phi_0]) +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_21_30_1m1[_stride_rho_0]) *
                       (_data_phi_20_30_10[_stride_phi_0] +
                        _data_phi_20_30_1m1[_stride_phi_0] -
                        _data_phi_21_30_10[_stride_phi_0] -
                        _data_phi_21_30_1m1[_stride_phi_0])) *
                  0.125 * 1.4142135623730951 /
                  (kT * (1.0 + 2.8284271247461903 +
                         1.7320508075688772 * 1.33333333333333));
            }
            if (_size_j_1 > 1 && _size_j_2 > 1) {
              double *RESTRICT _data_j_20_310 = _data_j + 10 * _stride_j_3;
              double *RESTRICT _data_j_20_310_10 =
                  _stride_j_1 * (_size_j_1 - 1) + _data_j_20_310;
              double *RESTRICT _data_rho_20_30 = _data_rho;
              double *RESTRICT _data_rho_20_30_10 =
                  _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
              double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
              double *RESTRICT _data_rho_21_30_1m1 =
                  _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                  _data_rho_21_30;
              double *RESTRICT _data_phi_20_30 = _data_phi;
              double *RESTRICT _data_phi_20_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_20_30;
              double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
              double *RESTRICT _data_phi_21_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_21_30;
              double *RESTRICT _data_phi_20_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
              double *RESTRICT _data_phi_21_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21_30;
              double *RESTRICT _data_rho_20_30_1m1 =
                  _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                  _data_rho_20_30;
              double *RESTRICT _data_rho_21_30_10 =
                  _stride_rho_1 * (_size_j_1 - 1) + _data_rho_21_30;
              _data_j_20_310_10[_stride_j_0] =
                  D *
                  (f_ext0 * z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_21_30_1m1[0]) *
                       2.0 +
                   f_ext1 * z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_21_30_1m1[0]) *
                       2.0 +
                   f_ext2 * z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_21_30_1m1[0]) *
                       -2.0 +
                   kT *
                       (-3.0 * _data_rho_21_30_1m1[0] +
                        3.0 * _data_rho_20_30_10[_stride_rho_0] -
                        _data_rho_20_30_10[0] + _data_rho_20_30_1m1[0] -
                        _data_rho_20_30_1m1[_stride_rho_0] +
                        _data_rho_21_30_10[0] -
                        _data_rho_21_30_10[_stride_rho_0] +
                        _data_rho_21_30_1m1[_stride_rho_0]) *
                       2.0 +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_21_30_1m1[0]) *
                       (-_data_phi_20_30_10[0] +
                        _data_phi_20_30_10[_stride_phi_0] -
                        _data_phi_21_30_1m1[0] +
                        _data_phi_21_30_1m1[_stride_phi_0]) +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_21_30_1m1[0]) *
                       (_data_phi_20_30_10[_stride_phi_0] +
                        _data_phi_20_30_1m1[0] -
                        _data_phi_21_30_10[_stride_phi_0] -
                        _data_phi_21_30_1m1[0]) +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_21_30_1m1[0]) *
                       (_data_phi_20_30_10[_stride_phi_0] -
                        _data_phi_20_30_1m1[_stride_phi_0] +
                        _data_phi_21_30_10[0] - _data_phi_21_30_1m1[0])) *
                  0.0833333333333333 * 1.7320508075688772 /
                  (kT * (1.0 + 2.8284271247461903 +
                         1.7320508075688772 * 1.33333333333333));
            }
          }
#pragma omp for schedule(static)
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (_size_j_1 > 1 && _size_j_2 > 1 && _size_j_0 - ctr_0 > 1) {
              double *RESTRICT _data_j_20_38 = _data_j + 8 * _stride_j_3;
              double *RESTRICT _data_j_20_38_10 =
                  _stride_j_1 * (_size_j_1 - 1) + _data_j_20_38;
              double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
              double *RESTRICT _data_rho_21_30_1m1 =
                  _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                  _data_rho_21_30;
              double *RESTRICT _data_rho_20_30 = _data_rho;
              double *RESTRICT _data_rho_20_30_10 =
                  _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
              double *RESTRICT _data_phi_20_30 = _data_phi;
              double *RESTRICT _data_phi_20_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_20_30;
              double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
              double *RESTRICT _data_phi_21_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_21_30;
              double *RESTRICT _data_phi_20_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
              double *RESTRICT _data_phi_21_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21_30;
              _data_j_20_38_10[_stride_j_0 * ctr_0] =
                  D *
                  (f_ext1 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                       2.0 +
                   f_ext2 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                       -2.0 +
                   kT *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                        _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                       4.0 +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                       (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                        _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] +
                        _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                        _data_phi_21_30_1m1[_stride_phi_0 * ctr_0]) +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                       (_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                        _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] -
                        _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                        _data_phi_21_30_1m1[_stride_phi_0 * ctr_0])) *
                  0.125 * 1.4142135623730951 /
                  (kT * (1.0 + 2.8284271247461903 +
                         1.7320508075688772 * 1.33333333333333));
            }
            if (_size_j_1 > 1 && _size_j_2 > 1) {
              double *RESTRICT _data_j_20_310 = _data_j + 10 * _stride_j_3;
              double *RESTRICT _data_j_20_310_10 =
                  _stride_j_1 * (_size_j_1 - 1) + _data_j_20_310;
              double *RESTRICT _data_rho_20_30 = _data_rho;
              double *RESTRICT _data_rho_20_30_10 =
                  _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
              double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
              double *RESTRICT _data_rho_21_30_1m1 =
                  _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                  _data_rho_21_30;
              double *RESTRICT _data_phi_20_30 = _data_phi;
              double *RESTRICT _data_phi_20_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_20_30;
              double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
              double *RESTRICT _data_phi_21_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_21_30;
              double *RESTRICT _data_phi_20_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
              double *RESTRICT _data_phi_21_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21_30;
              double *RESTRICT _data_rho_20_30_1m1 =
                  _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                  _data_rho_20_30;
              double *RESTRICT _data_rho_21_30_10 =
                  _stride_rho_1 * (_size_j_1 - 1) + _data_rho_21_30;
              _data_j_20_310_10[_stride_j_0 * ctr_0] =
                  D *
                  (f_ext0 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0]) *
                       2.0 +
                   f_ext1 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0]) *
                       2.0 +
                   f_ext2 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0]) *
                       -2.0 +
                   kT *
                       (-3.0 * _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                                   _stride_rho_0] +
                        3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                        _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0] +
                        _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0] -
                        _data_rho_20_30_1m1[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0] -
                        _data_rho_21_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                       2.0 +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0]) *
                       (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                            _stride_phi_0] +
                        _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                        _data_phi_21_30_1m1[_stride_phi_0 * ctr_0 -
                                            _stride_phi_0] +
                        _data_phi_21_30_1m1[_stride_phi_0 * ctr_0]) +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0]) *
                       (_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                        _data_phi_20_30_1m1[_stride_phi_0 * ctr_0 -
                                            _stride_phi_0] -
                        _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                        _data_phi_21_30_1m1[_stride_phi_0 * ctr_0 -
                                            _stride_phi_0]) +
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0]) *
                       (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                        _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] +
                        _data_phi_21_30_10[_stride_phi_0 * ctr_0 -
                                           _stride_phi_0] -
                        _data_phi_21_30_1m1[_stride_phi_0 * ctr_0 -
                                            _stride_phi_0])) *
                  0.0833333333333333 * 1.7320508075688772 /
                  (kT * (1.0 + 2.8284271247461903 +
                         1.7320508075688772 * 1.33333333333333));
            }
          }
          if (_size_j_1 > 1 && _size_j_2 > 1) {
            double *RESTRICT _data_j_20_310 = _data_j + 10 * _stride_j_3;
            double *RESTRICT _data_j_20_310_10 =
                _stride_j_1 * (_size_j_1 - 1) + _data_j_20_310;
            double *RESTRICT _data_rho_20_30 = _data_rho;
            double *RESTRICT _data_rho_20_30_10 =
                _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
            double *RESTRICT _data_rho_21_30 = _data_rho + _stride_rho_2;
            double *RESTRICT _data_rho_21_30_1m1 =
                _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                _data_rho_21_30;
            double *RESTRICT _data_phi_20_30 = _data_phi;
            double *RESTRICT _data_phi_20_30_1m1 =
                _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                _data_phi_20_30;
            double *RESTRICT _data_phi_21_30 = _data_phi + _stride_phi_2;
            double *RESTRICT _data_phi_21_30_1m1 =
                _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                _data_phi_21_30;
            double *RESTRICT _data_phi_20_30_10 =
                _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
            double *RESTRICT _data_phi_21_30_10 =
                _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21_30;
            double *RESTRICT _data_rho_20_30_1m1 =
                _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                _data_rho_20_30;
            double *RESTRICT _data_rho_21_30_10 =
                _stride_rho_1 * (_size_j_1 - 1) + _data_rho_21_30;
            _data_j_20_310_10[_stride_j_0 * (_size_j_0 - 1)] =
                D *
                (f_ext0 * z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0]) *
                     2.0 +
                 f_ext1 * z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0]) *
                     2.0 +
                 f_ext2 * z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0]) *
                     -2.0 +
                 kT *
                     (-3.0 *
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0] +
                      3.0 *
                          _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                      _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0] +
                      _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0] -
                      _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0] -
                      _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1)]) *
                     2.0 +
                 z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0]) *
                     (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                          _stride_phi_0] +
                      _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                      _data_phi_21_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                          _stride_phi_0] +
                      _data_phi_21_30_1m1[_stride_phi_0 * (_size_j_0 - 1)]) +
                 z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0]) *
                     (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                      _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                          _stride_phi_0] -
                      _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                      _data_phi_21_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                          _stride_phi_0]) +
                 z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0]) *
                     (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                      _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1)] +
                      _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                         _stride_phi_0] -
                      _data_phi_21_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                          _stride_phi_0])) *
                0.0833333333333333 * 1.7320508075688772 /
                (kT * (1.0 + 2.8284271247461903 +
                       1.7320508075688772 * 1.33333333333333));
          }
        }
      }
#pragma omp for schedule(static)
      for (int64_t ctr_2 = 1; ctr_2 < _size_j_2 - 1; ctr_2 += 1) {
        {
          {
            {
              if (_size_j_1 > 1 && ctr_2 > 0 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_34 =
                    _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
                double *RESTRICT _data_j_20_34_10 = _data_j_20_34;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
                _data_j_20_34_10[_stride_j_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_20_30_11[0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_20_30_11[0]) *
                         -2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0] -
                          _data_rho_20_30_11[0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_20_30_11[0]) *
                         (-_data_phi_20_30_10[0] +
                          _data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_11[0] +
                          _data_phi_20_30_11[_stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_20_30_11[0]) *
                         (_data_phi_20_30_10[0] +
                          _data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_11[0] -
                          _data_phi_20_30_11[_stride_phi_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && ctr_2 > 0) {
                double *RESTRICT _data_j_20_311 =
                    _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
                double *RESTRICT _data_j_20_311_10 = _data_j_20_311;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_11 =
                    _stride_rho_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_2m1_30_10 = _data_rho_2m1_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_11 =
                    _stride_phi_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_10 = _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 + _data_phi_20_30;
                _data_j_20_311_10[_stride_j_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_11[0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_11[0]) *
                         -2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_11[0]) *
                         2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0] +
                          3.0 * _data_rho_2m1_30_11[0] + _data_rho_20_30_10[0] -
                          _data_rho_20_30_11[0] +
                          _data_rho_20_30_11[_stride_rho_0] -
                          _data_rho_2m1_30_10[0] +
                          _data_rho_2m1_30_10[_stride_rho_0] -
                          _data_rho_2m1_30_11[_stride_rho_0]) *
                         -2.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_11[0]) *
                         (_data_phi_20_30_10[0] -
                          _data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_11[0] -
                          _data_phi_2m1_30_11[_stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_11[0]) *
                         (-_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_11[0] +
                          _data_phi_2m1_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_11[0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_11[0]) *
                         (-_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_20_30_11[_stride_phi_0] -
                          _data_phi_2m1_30_10[0] + _data_phi_2m1_30_11[0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_312 =
                    _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
                double *RESTRICT _data_j_20_312_10 = _data_j_20_312;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
                double *RESTRICT _data_rho_21_30 =
                    _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_11 =
                    _stride_rho_1 + _data_rho_21_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 =
                    _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_11 =
                    _stride_phi_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 = _data_phi_21_30;
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30_10 = _data_rho_21_30;
                _data_j_20_312_10[_stride_j_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_11[0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_11[0]) *
                         -2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_11[0]) *
                         -2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0] +
                          3.0 * _data_rho_21_30_11[0] + _data_rho_20_30_10[0] -
                          _data_rho_20_30_11[0] +
                          _data_rho_20_30_11[_stride_rho_0] -
                          _data_rho_21_30_10[0] +
                          _data_rho_21_30_10[_stride_rho_0] -
                          _data_rho_21_30_11[_stride_rho_0]) *
                         -2.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_11[0]) *
                         (-_data_phi_20_30_10[0] +
                          _data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_21_30_11[0] +
                          _data_phi_21_30_11[_stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_11[0]) *
                         (_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_20_30_11[0] -
                          _data_phi_21_30_10[_stride_phi_0] -
                          _data_phi_21_30_11[0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_11[0]) *
                         (_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_11[_stride_phi_0] +
                          _data_phi_21_30_10[0] - _data_phi_21_30_11[0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
            }
            for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
              if (_size_j_1 > 1 && ctr_2 > 0 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_34 =
                    _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
                double *RESTRICT _data_j_20_34_10 = _data_j_20_34;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
                _data_j_20_34_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         -2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_11[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] +
                          _data_phi_20_30_11[_stride_phi_0 * ctr_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_11[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_20_30_11[_stride_phi_0 * ctr_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && ctr_2 > 0) {
                double *RESTRICT _data_j_20_311 =
                    _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
                double *RESTRICT _data_j_20_311_10 = _data_j_20_311;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_11 =
                    _stride_rho_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_2m1_30_10 = _data_rho_2m1_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_11 =
                    _stride_phi_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_10 = _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 + _data_phi_20_30;
                _data_j_20_311_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         -2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          3.0 * _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                                    _stride_rho_0] +
                          _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] -
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] +
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0] -
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0]) *
                         -2.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_11[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] -
                          _data_phi_2m1_30_11[_stride_phi_0 * ctr_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_11[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_11[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_20_30_11[_stride_phi_0 * ctr_0] -
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_11[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_312 =
                    _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
                double *RESTRICT _data_j_20_312_10 = _data_j_20_312;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
                double *RESTRICT _data_rho_21_30 =
                    _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_11 =
                    _stride_rho_1 + _data_rho_21_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 =
                    _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_11 =
                    _stride_phi_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 = _data_phi_21_30;
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30_10 = _data_rho_21_30;
                _data_j_20_312_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         -2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         -2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          3.0 * _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                                   _stride_rho_0] +
                          _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] -
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] +
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0] -
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] +
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0]) *
                         -2.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_11[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] +
                          _data_phi_21_30_11[_stride_phi_0 * ctr_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_20_30_11[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_11[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_11[_stride_phi_0 * ctr_0] +
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_21_30_11[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
            }
            {
              if (_size_j_1 > 1 && ctr_2 > 0 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_34 =
                    _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
                double *RESTRICT _data_j_20_34_10 = _data_j_20_34;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
                _data_j_20_34_10[_stride_j_0 * (_size_j_0 - 1)] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         -2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] +
                          _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1)]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_20_30_11[_stride_phi_0 *
                                             (_size_j_0 - 1)])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && ctr_2 > 0) {
                double *RESTRICT _data_j_20_311 =
                    _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
                double *RESTRICT _data_j_20_311_10 = _data_j_20_311;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_11 =
                    _stride_rho_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_2m1_30_10 = _data_rho_2m1_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_11 =
                    _stride_phi_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_10 = _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 + _data_phi_20_30;
                _data_j_20_311_10[_stride_j_0 * (_size_j_0 - 1)] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         -2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1)] +
                          3.0 * _data_rho_2m1_30_11[_stride_rho_0 *
                                                        (_size_j_0 - 1) -
                                                    _stride_rho_0] +
                          _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] -
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] +
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_2m1_30_11[_stride_rho_0 *
                                              (_size_j_0 - 1)]) *
                         -2.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_2m1_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] -
                          _data_phi_2m1_30_11[_stride_phi_0 *
                                              (_size_j_0 - 1)]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_2m1_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_312 =
                    _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
                double *RESTRICT _data_j_20_312_10 = _data_j_20_312;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
                double *RESTRICT _data_rho_21_30 =
                    _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_11 =
                    _stride_rho_1 + _data_rho_21_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 =
                    _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_11 =
                    _stride_phi_1 + _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 = _data_phi_21_30;
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30_10 = _data_rho_21_30;
                _data_j_20_312_10[_stride_j_0 * (_size_j_0 - 1)] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         -2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         -2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1)] +
                          3.0 * _data_rho_21_30_11[_stride_rho_0 *
                                                       (_size_j_0 - 1) -
                                                   _stride_rho_0] +
                          _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] -
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] +
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] +
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1)]) *
                         -2.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] +
                          _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1)]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
            }
          }
          for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
            {
              {
                if (ctr_1 > 0 && ctr_2 > 0 && _size_j_1 - ctr_1 > 1 &&
                    _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_30 =
                      _data_j + _stride_j_2 * ctr_2;
                  double *RESTRICT _data_j_20_30_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  _data_j_20_30_10[_stride_j_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[0] +
                            _data_rho_20_30_10[_stride_rho_0]) +
                       kT *
                           (-_data_rho_20_30_10[0] +
                            _data_rho_20_30_10[_stride_rho_0]) *
                           2.0 +
                       z *
                           (-_data_phi_20_30_10[0] +
                            _data_phi_20_30_10[_stride_phi_0]) *
                           (_data_rho_20_30_10[0] +
                            _data_rho_20_30_10[_stride_rho_0])) *
                      0.5 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (_size_j_0 > 2 && ctr_1 > 0 && ctr_2 > 0 &&
                    _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_31 =
                      _data_j + _stride_j_2 * ctr_2 + _stride_j_3;
                  double *RESTRICT _data_j_20_31_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_31;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                  _data_j_20_31_10[_stride_j_0] =
                      D *
                      (f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_20_30_1m1[_stride_rho_0]) +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0] -
                            _data_rho_20_30_1m1[_stride_rho_0]) *
                           2.0 +
                       z *
                           (_data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_20_30_1m1[_stride_phi_0]) *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_20_30_1m1[_stride_rho_0])) *
                      0.5 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (_size_j_0 > 2 && ctr_1 > 0 && ctr_2 > 0 &&
                    _size_j_1 - ctr_1 > 1) {
                  double *RESTRICT _data_j_20_32 =
                      _data_j + _stride_j_2 * ctr_2 + 2 * _stride_j_3;
                  double *RESTRICT _data_j_20_32_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_32;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  _data_j_20_32_10[_stride_j_0] =
                      D *
                      (f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_10[_stride_rho_0]) +
                       kT *
                           (-_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_10[_stride_rho_0]) *
                           -2.0 -
                       z *
                           (-_data_phi_20_30_10[_stride_phi_0] +
                            _data_phi_2m1_30_10[_stride_phi_0]) *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_10[_stride_rho_0])) *
                      0.5 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && ctr_2 > 0 && _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_33 =
                      _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
                  double *RESTRICT _data_j_20_33_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_33;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  _data_j_20_33_10[_stride_j_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_20_30_1m1[0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_20_30_1m1[0]) *
                           2.0 +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0] -
                            _data_rho_20_30_1m1[0]) *
                           4.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_20_30_1m1[0]) *
                           (-_data_phi_20_30_10[0] +
                            _data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_20_30_1m1[0] +
                            _data_phi_20_30_1m1[_stride_phi_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_20_30_1m1[0]) *
                           (_data_phi_20_30_10[0] +
                            _data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_20_30_1m1[0] -
                            _data_phi_20_30_1m1[_stride_phi_0])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_2 > 0 && _size_j_1 - ctr_1 > 1 &&
                    _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_34 =
                      _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
                  double *RESTRICT _data_j_20_34_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_34;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  _data_j_20_34_10[_stride_j_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_20_30_11[0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_20_30_11[0]) *
                           -2.0 +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0] -
                            _data_rho_20_30_11[0]) *
                           4.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_20_30_11[0]) *
                           (-_data_phi_20_30_10[0] +
                            _data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_20_30_11[0] +
                            _data_phi_20_30_11[_stride_phi_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_20_30_11[0]) *
                           (_data_phi_20_30_10[0] +
                            _data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_20_30_11[0] -
                            _data_phi_20_30_11[_stride_phi_0])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && ctr_2 > 0 && _size_j_1 - ctr_1 > 1) {
                  double *RESTRICT _data_j_20_35 =
                      _data_j + _stride_j_2 * ctr_2 + 5 * _stride_j_3;
                  double *RESTRICT _data_j_20_35_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_35;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  _data_j_20_35_10[_stride_j_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_10[0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_10[0]) *
                           2.0 +
                       kT *
                           (-_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_10[0]) *
                           -4.0 -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_10[0]) *
                           (-_data_phi_20_30_10[0] -
                            _data_phi_20_30_10[_stride_phi_0] +
                            _data_phi_2m1_30_10[0] +
                            _data_phi_2m1_30_10[_stride_phi_0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_10[0]) *
                           (_data_phi_20_30_10[0] -
                            _data_phi_20_30_10[_stride_phi_0] +
                            _data_phi_2m1_30_10[0] -
                            _data_phi_2m1_30_10[_stride_phi_0])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && _size_j_1 - ctr_1 > 1 &&
                    _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_36 =
                      _data_j + _stride_j_2 * ctr_2 + 6 * _stride_j_3;
                  double *RESTRICT _data_j_20_36_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_36;
                  double *RESTRICT _data_rho_21_30 =
                      _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                  double *RESTRICT _data_rho_21_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_21_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_21_30 =
                      _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                  double *RESTRICT _data_phi_21_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_21_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  _data_j_20_36_10[_stride_j_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_10[0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_10[0]) *
                           -2.0 +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0] -
                            _data_rho_21_30_10[0]) *
                           4.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_10[0]) *
                           (-_data_phi_20_30_10[0] +
                            _data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_21_30_10[0] +
                            _data_phi_21_30_10[_stride_phi_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_10[0]) *
                           (_data_phi_20_30_10[0] +
                            _data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_21_30_10[0] -
                            _data_phi_21_30_10[_stride_phi_0])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (_size_j_0 > 2 && ctr_1 > 0 && ctr_2 > 0) {
                  double *RESTRICT _data_j_20_37 =
                      _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
                  double *RESTRICT _data_j_20_37_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_37;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_2m1_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  _data_j_20_37_10[_stride_j_0] =
                      D *
                      (f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0]) *
                           2.0 +
                       kT *
                           (-_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0]) *
                           -4.0 -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_20_30_1m1[_stride_phi_0] +
                            _data_phi_2m1_30_10[_stride_phi_0] +
                            _data_phi_2m1_30_1m1[_stride_phi_0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0] +
                            _data_phi_20_30_1m1[_stride_phi_0] -
                            _data_phi_2m1_30_10[_stride_phi_0] +
                            _data_phi_2m1_30_1m1[_stride_phi_0])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (_size_j_0 > 2 && ctr_1 > 0 && _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_38 =
                      _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
                  double *RESTRICT _data_j_20_38_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_38;
                  double *RESTRICT _data_rho_21_30 =
                      _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                  double *RESTRICT _data_rho_21_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30 =
                      _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                  double *RESTRICT _data_phi_21_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_21_30;
                  _data_j_20_38_10[_stride_j_0] =
                      D *
                      (f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_1m1[_stride_rho_0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_1m1[_stride_rho_0]) *
                           -2.0 +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0] -
                            _data_rho_21_30_1m1[_stride_rho_0]) *
                           4.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_1m1[_stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_20_30_1m1[_stride_phi_0] +
                            _data_phi_21_30_10[_stride_phi_0] -
                            _data_phi_21_30_1m1[_stride_phi_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_1m1[_stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0] +
                            _data_phi_20_30_1m1[_stride_phi_0] -
                            _data_phi_21_30_10[_stride_phi_0] -
                            _data_phi_21_30_1m1[_stride_phi_0])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && ctr_2 > 0) {
                  double *RESTRICT _data_j_20_39 =
                      _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
                  double *RESTRICT _data_j_20_39_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_39;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_rho_2m1_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  _data_j_20_39_10[_stride_j_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_1m1[0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_1m1[0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_1m1[0]) *
                           2.0 +
                       kT *
                           (-3.0 * _data_rho_20_30_10[_stride_rho_0] +
                            3.0 * _data_rho_2m1_30_1m1[0] +
                            _data_rho_20_30_10[0] - _data_rho_20_30_1m1[0] +
                            _data_rho_20_30_1m1[_stride_rho_0] -
                            _data_rho_2m1_30_10[0] +
                            _data_rho_2m1_30_10[_stride_rho_0] -
                            _data_rho_2m1_30_1m1[_stride_rho_0]) *
                           -2.0 -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_1m1[0]) *
                           (_data_phi_20_30_10[0] -
                            _data_phi_20_30_10[_stride_phi_0] +
                            _data_phi_2m1_30_1m1[0] -
                            _data_phi_2m1_30_1m1[_stride_phi_0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_1m1[0]) *
                           (-_data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_20_30_1m1[0] +
                            _data_phi_2m1_30_10[_stride_phi_0] +
                            _data_phi_2m1_30_1m1[0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_1m1[0]) *
                           (-_data_phi_20_30_10[_stride_phi_0] +
                            _data_phi_20_30_1m1[_stride_phi_0] -
                            _data_phi_2m1_30_10[0] + _data_phi_2m1_30_1m1[0])) *
                      0.0833333333333333 * 1.7320508075688772 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_310 =
                      _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
                  double *RESTRICT _data_j_20_310_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_310;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_21_30 =
                      _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                  double *RESTRICT _data_rho_21_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30 =
                      _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                  double *RESTRICT _data_phi_21_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_21_30;
                  double *RESTRICT _data_rho_20_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_21_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_21_30;
                  _data_j_20_310_10[_stride_j_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_1m1[0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_1m1[0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_1m1[0]) *
                           -2.0 +
                       kT *
                           (-3.0 * _data_rho_21_30_1m1[0] +
                            3.0 * _data_rho_20_30_10[_stride_rho_0] -
                            _data_rho_20_30_10[0] + _data_rho_20_30_1m1[0] -
                            _data_rho_20_30_1m1[_stride_rho_0] +
                            _data_rho_21_30_10[0] -
                            _data_rho_21_30_10[_stride_rho_0] +
                            _data_rho_21_30_1m1[_stride_rho_0]) *
                           2.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_1m1[0]) *
                           (-_data_phi_20_30_10[0] +
                            _data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_21_30_1m1[0] +
                            _data_phi_21_30_1m1[_stride_phi_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_1m1[0]) *
                           (_data_phi_20_30_10[_stride_phi_0] +
                            _data_phi_20_30_1m1[0] -
                            _data_phi_21_30_10[_stride_phi_0] -
                            _data_phi_21_30_1m1[0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_1m1[0]) *
                           (_data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_20_30_1m1[_stride_phi_0] +
                            _data_phi_21_30_10[0] - _data_phi_21_30_1m1[0])) *
                      0.0833333333333333 * 1.7320508075688772 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_2 > 0 && _size_j_1 - ctr_1 > 1) {
                  double *RESTRICT _data_j_20_311 =
                      _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
                  double *RESTRICT _data_j_20_311_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_311;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_rho_2m1_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                  _data_j_20_311_10[_stride_j_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_11[0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_11[0]) *
                           -2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_11[0]) *
                           2.0 +
                       kT *
                           (-3.0 * _data_rho_20_30_10[_stride_rho_0] +
                            3.0 * _data_rho_2m1_30_11[0] +
                            _data_rho_20_30_10[0] - _data_rho_20_30_11[0] +
                            _data_rho_20_30_11[_stride_rho_0] -
                            _data_rho_2m1_30_10[0] +
                            _data_rho_2m1_30_10[_stride_rho_0] -
                            _data_rho_2m1_30_11[_stride_rho_0]) *
                           -2.0 -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_11[0]) *
                           (_data_phi_20_30_10[0] -
                            _data_phi_20_30_10[_stride_phi_0] +
                            _data_phi_2m1_30_11[0] -
                            _data_phi_2m1_30_11[_stride_phi_0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_11[0]) *
                           (-_data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_20_30_11[0] +
                            _data_phi_2m1_30_10[_stride_phi_0] +
                            _data_phi_2m1_30_11[0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_2m1_30_11[0]) *
                           (-_data_phi_20_30_10[_stride_phi_0] +
                            _data_phi_20_30_11[_stride_phi_0] -
                            _data_phi_2m1_30_10[0] + _data_phi_2m1_30_11[0])) *
                      0.0833333333333333 * 1.7320508075688772 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (_size_j_1 - ctr_1 > 1 && _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_312 =
                      _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
                  double *RESTRICT _data_j_20_312_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_312;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_21_30 =
                      _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                  double *RESTRICT _data_rho_21_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_21_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30 =
                      _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                  double *RESTRICT _data_phi_21_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_21_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_21_30;
                  double *RESTRICT _data_rho_20_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_21_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_21_30;
                  _data_j_20_312_10[_stride_j_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_11[0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_11[0]) *
                           -2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_11[0]) *
                           -2.0 +
                       kT *
                           (-3.0 * _data_rho_20_30_10[_stride_rho_0] +
                            3.0 * _data_rho_21_30_11[0] +
                            _data_rho_20_30_10[0] - _data_rho_20_30_11[0] +
                            _data_rho_20_30_11[_stride_rho_0] -
                            _data_rho_21_30_10[0] +
                            _data_rho_21_30_10[_stride_rho_0] -
                            _data_rho_21_30_11[_stride_rho_0]) *
                           -2.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_11[0]) *
                           (-_data_phi_20_30_10[0] +
                            _data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_21_30_11[0] +
                            _data_phi_21_30_11[_stride_phi_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_11[0]) *
                           (_data_phi_20_30_10[_stride_phi_0] +
                            _data_phi_20_30_11[0] -
                            _data_phi_21_30_10[_stride_phi_0] -
                            _data_phi_21_30_11[0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0] +
                            _data_rho_21_30_11[0]) *
                           (_data_phi_20_30_10[_stride_phi_0] -
                            _data_phi_20_30_11[_stride_phi_0] +
                            _data_phi_21_30_10[0] - _data_phi_21_30_11[0])) *
                      0.0833333333333333 * 1.7320508075688772 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
              }
              for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
                if (ctr_1 > 0 && ctr_2 > 0 && _size_j_1 - ctr_1 > 1 &&
                    _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_30 =
                      _data_j + _stride_j_2 * ctr_2;
                  double *RESTRICT _data_j_20_30_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  _data_j_20_30_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0] +
                            _data_rho_20_30_10[_stride_rho_0 * ctr_0]) +
                       kT *
                           (-_data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0] +
                            _data_rho_20_30_10[_stride_rho_0 * ctr_0]) *
                           2.0 +
                       z *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0]) *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0] +
                            _data_rho_20_30_10[_stride_rho_0 * ctr_0])) *
                      0.5 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && ctr_2 > 0 && _size_j_0 - ctr_0 > 1 &&
                    _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_31 =
                      _data_j + _stride_j_2 * ctr_2 + _stride_j_3;
                  double *RESTRICT _data_j_20_31_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_31;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                  _data_j_20_31_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_20_30_1m1[_stride_rho_0 * ctr_0]) +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                            _data_rho_20_30_1m1[_stride_rho_0 * ctr_0]) *
                           2.0 +
                       z *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0]) *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_20_30_1m1[_stride_rho_0 * ctr_0])) *
                      0.5 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && ctr_2 > 0 && _size_j_0 - ctr_0 > 1 &&
                    _size_j_1 - ctr_1 > 1) {
                  double *RESTRICT _data_j_20_32 =
                      _data_j + _stride_j_2 * ctr_2 + 2 * _stride_j_3;
                  double *RESTRICT _data_j_20_32_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_32;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  _data_j_20_32_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_10[_stride_rho_0 * ctr_0]) +
                       kT *
                           (-_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_10[_stride_rho_0 * ctr_0]) *
                           -2.0 -
                       z *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_2m1_30_10[_stride_phi_0 * ctr_0]) *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_10[_stride_rho_0 * ctr_0])) *
                      0.5 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && ctr_2 > 0 && _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_33 =
                      _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
                  double *RESTRICT _data_j_20_33_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_33;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  _data_j_20_33_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           2.0 +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                            _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           4.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] +
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] -
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_2 > 0 && _size_j_1 - ctr_1 > 1 &&
                    _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_34 =
                      _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
                  double *RESTRICT _data_j_20_34_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_34;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  _data_j_20_34_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           -2.0 +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                            _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           4.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_20_30_11[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] +
                            _data_phi_20_30_11[_stride_phi_0 * ctr_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_20_30_11[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] -
                            _data_phi_20_30_11[_stride_phi_0 * ctr_0])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && ctr_2 > 0 && _size_j_1 - ctr_1 > 1) {
                  double *RESTRICT _data_j_20_35 =
                      _data_j + _stride_j_2 * ctr_2 + 5 * _stride_j_3;
                  double *RESTRICT _data_j_20_35_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_35;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  _data_j_20_35_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           2.0 +
                       kT *
                           (-_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           -4.0 -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] -
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_2m1_30_10[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] +
                            _data_phi_2m1_30_10[_stride_phi_0 * ctr_0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] -
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_2m1_30_10[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] -
                            _data_phi_2m1_30_10[_stride_phi_0 * ctr_0])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && _size_j_1 - ctr_1 > 1 &&
                    _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_36 =
                      _data_j + _stride_j_2 * ctr_2 + 6 * _stride_j_3;
                  double *RESTRICT _data_j_20_36_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_36;
                  double *RESTRICT _data_rho_21_30 =
                      _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                  double *RESTRICT _data_rho_21_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_21_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_21_30 =
                      _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                  double *RESTRICT _data_phi_21_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_21_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  _data_j_20_36_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           -2.0 +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                            _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           4.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_21_30_10[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] +
                            _data_phi_21_30_10[_stride_phi_0 * ctr_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_21_30_10[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] -
                            _data_phi_21_30_10[_stride_phi_0 * ctr_0])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && ctr_2 > 0 && _size_j_0 - ctr_0 > 1) {
                  double *RESTRICT _data_j_20_37 =
                      _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
                  double *RESTRICT _data_j_20_37_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_37;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_2m1_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  _data_j_20_37_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                           2.0 +
                       kT *
                           (-_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                           -4.0 -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] +
                            _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] -
                            _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && _size_j_0 - ctr_0 > 1 &&
                    _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_38 =
                      _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
                  double *RESTRICT _data_j_20_38_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_38;
                  double *RESTRICT _data_rho_21_30 =
                      _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                  double *RESTRICT _data_rho_21_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30 =
                      _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                  double *RESTRICT _data_phi_21_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_21_30;
                  _data_j_20_38_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                           -2.0 +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                            _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                           4.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] +
                            _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_21_30_1m1[_stride_phi_0 * ctr_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] -
                            _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_21_30_1m1[_stride_phi_0 * ctr_0])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && ctr_2 > 0) {
                  double *RESTRICT _data_j_20_39 =
                      _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
                  double *RESTRICT _data_j_20_39_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_39;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_rho_2m1_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  _data_j_20_39_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                                 _stride_rho_0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                                 _stride_rho_0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                                 _stride_rho_0]) *
                           2.0 +
                       kT *
                           (-3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            3.0 * _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                                       _stride_rho_0] +
                            _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0] -
                            _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0] +
                            _data_rho_20_30_1m1[_stride_rho_0 * ctr_0] -
                            _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0] +
                            _data_rho_2m1_30_10[_stride_rho_0 * ctr_0] -
                            _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                           -2.0 -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                                 _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] -
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0 -
                                                 _stride_phi_0] -
                            _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                                 _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] +
                            _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0 -
                                                 _stride_phi_0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                                 _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] -
                            _data_phi_2m1_30_10[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] +
                            _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0 -
                                                 _stride_phi_0])) *
                      0.0833333333333333 * 1.7320508075688772 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_310 =
                      _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
                  double *RESTRICT _data_j_20_310_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_310;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_21_30 =
                      _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                  double *RESTRICT _data_rho_21_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30 =
                      _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                  double *RESTRICT _data_phi_21_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_21_30;
                  double *RESTRICT _data_rho_20_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_21_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_21_30;
                  _data_j_20_310_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           -2.0 +
                       kT *
                           (-3.0 * _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                                       _stride_rho_0] +
                            3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                            _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0] +
                            _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0] -
                            _data_rho_20_30_1m1[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0] -
                            _data_rho_21_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                           2.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_21_30_1m1[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] +
                            _data_phi_21_30_1m1[_stride_phi_0 * ctr_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] -
                            _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_21_30_1m1[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] +
                            _data_phi_21_30_10[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] -
                            _data_phi_21_30_1m1[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0])) *
                      0.0833333333333333 * 1.7320508075688772 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_2 > 0 && _size_j_1 - ctr_1 > 1) {
                  double *RESTRICT _data_j_20_311 =
                      _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
                  double *RESTRICT _data_j_20_311_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_311;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_rho_2m1_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                  _data_j_20_311_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           -2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           2.0 +
                       kT *
                           (-3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            3.0 * _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                                      _stride_rho_0] +
                            _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0] -
                            _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0] +
                            _data_rho_20_30_11[_stride_rho_0 * ctr_0] -
                            _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0] +
                            _data_rho_2m1_30_10[_stride_rho_0 * ctr_0] -
                            _data_rho_2m1_30_11[_stride_rho_0 * ctr_0]) *
                           -2.0 -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] -
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_2m1_30_11[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] -
                            _data_phi_2m1_30_11[_stride_phi_0 * ctr_0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_20_30_11[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] +
                            _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_2m1_30_11[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                                _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_20_30_11[_stride_phi_0 * ctr_0] -
                            _data_phi_2m1_30_10[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] +
                            _data_phi_2m1_30_11[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0])) *
                      0.0833333333333333 * 1.7320508075688772 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (_size_j_1 - ctr_1 > 1 && _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_312 =
                      _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
                  double *RESTRICT _data_j_20_312_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_312;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_21_30 =
                      _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                  double *RESTRICT _data_rho_21_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_21_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30 =
                      _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                  double *RESTRICT _data_phi_21_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_21_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_21_30;
                  double *RESTRICT _data_rho_20_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_21_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_21_30;
                  _data_j_20_312_10[_stride_j_0 * ctr_0] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           -2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           -2.0 +
                       kT *
                           (-3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            3.0 * _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                                     _stride_rho_0] +
                            _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0] -
                            _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0] +
                            _data_rho_20_30_11[_stride_rho_0 * ctr_0] -
                            _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0] +
                            _data_rho_21_30_10[_stride_rho_0 * ctr_0] -
                            _data_rho_21_30_11[_stride_rho_0 * ctr_0]) *
                           -2.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                                _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_21_30_11[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] +
                            _data_phi_21_30_11[_stride_phi_0 * ctr_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                            _data_phi_20_30_11[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] -
                            _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_21_30_11[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                            _data_rho_21_30_11[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                            _data_phi_20_30_11[_stride_phi_0 * ctr_0] +
                            _data_phi_21_30_10[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] -
                            _data_phi_21_30_11[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0])) *
                      0.0833333333333333 * 1.7320508075688772 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
              }
              {
                if (ctr_1 > 0 && ctr_2 > 0 && _size_j_1 - ctr_1 > 1 &&
                    _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_30 =
                      _data_j + _stride_j_2 * ctr_2;
                  double *RESTRICT _data_j_20_30_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  _data_j_20_30_10[_stride_j_0 * (_size_j_0 - 1)] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0] +
                            _data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)]) +
                       kT *
                           (-_data_rho_20_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0] +
                            _data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)]) *
                           2.0 +
                       z *
                           (-_data_phi_20_30_10[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)]) *
                           (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0] +
                            _data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)])) *
                      0.5 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && ctr_2 > 0 && _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_33 =
                      _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
                  double *RESTRICT _data_j_20_33_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_33;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  _data_j_20_33_10[_stride_j_0 * (_size_j_0 - 1)] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_20_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_20_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           2.0 +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] -
                            _data_rho_20_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           4.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_20_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_20_30_1m1[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] +
                            _data_phi_20_30_1m1[_stride_phi_0 *
                                                (_size_j_0 - 1)]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_20_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_20_30_1m1[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] -
                            _data_phi_20_30_1m1[_stride_phi_0 *
                                                (_size_j_0 - 1)])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_2 > 0 && _size_j_1 - ctr_1 > 1 &&
                    _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_34 =
                      _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
                  double *RESTRICT _data_j_20_34_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_34;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  _data_j_20_34_10[_stride_j_0 * (_size_j_0 - 1)] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           -2.0 +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] -
                            _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           4.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] +
                            _data_phi_20_30_11[_stride_phi_0 *
                                               (_size_j_0 - 1)]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] -
                            _data_phi_20_30_11[_stride_phi_0 *
                                               (_size_j_0 - 1)])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && ctr_2 > 0 && _size_j_1 - ctr_1 > 1) {
                  double *RESTRICT _data_j_20_35 =
                      _data_j + _stride_j_2 * ctr_2 + 5 * _stride_j_3;
                  double *RESTRICT _data_j_20_35_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_35;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  _data_j_20_35_10[_stride_j_0 * (_size_j_0 - 1)] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           2.0 +
                       kT *
                           (-_data_rho_20_30_10[_stride_rho_0 *
                                                (_size_j_0 - 1)] +
                            _data_rho_2m1_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           -4.0 -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] -
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] +
                            _data_phi_2m1_30_10[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] +
                            _data_phi_2m1_30_10[_stride_phi_0 *
                                                (_size_j_0 - 1)]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] -
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] +
                            _data_phi_2m1_30_10[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] -
                            _data_phi_2m1_30_10[_stride_phi_0 *
                                                (_size_j_0 - 1)])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && _size_j_1 - ctr_1 > 1 &&
                    _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_36 =
                      _data_j + _stride_j_2 * ctr_2 + 6 * _stride_j_3;
                  double *RESTRICT _data_j_20_36_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_36;
                  double *RESTRICT _data_rho_21_30 =
                      _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                  double *RESTRICT _data_rho_21_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_21_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_21_30 =
                      _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                  double *RESTRICT _data_phi_21_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_21_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  _data_j_20_36_10[_stride_j_0 * (_size_j_0 - 1)] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           -2.0 +
                       kT *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] -
                            _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           4.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] +
                            _data_phi_21_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] -
                            _data_phi_21_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)])) *
                      0.125 * 1.4142135623730951 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && ctr_2 > 0) {
                  double *RESTRICT _data_j_20_39 =
                      _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
                  double *RESTRICT _data_j_20_39_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_39;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_rho_2m1_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  _data_j_20_39_10[_stride_j_0 * (_size_j_0 - 1)] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 *
                                                     (_size_j_0 - 1) -
                                                 _stride_rho_0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 *
                                                     (_size_j_0 - 1) -
                                                 _stride_rho_0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 *
                                                     (_size_j_0 - 1) -
                                                 _stride_rho_0]) *
                           2.0 +
                       kT *
                           (-3.0 * _data_rho_20_30_10[_stride_rho_0 *
                                                      (_size_j_0 - 1)] +
                            3.0 * _data_rho_2m1_30_1m1[_stride_rho_0 *
                                                           (_size_j_0 - 1) -
                                                       _stride_rho_0] +
                            _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0] -
                            _data_rho_20_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0] +
                            _data_rho_20_30_1m1[_stride_rho_0 *
                                                (_size_j_0 - 1)] -
                            _data_rho_2m1_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0] +
                            _data_rho_2m1_30_10[_stride_rho_0 *
                                                (_size_j_0 - 1)] -
                            _data_rho_2m1_30_1m1[_stride_rho_0 *
                                                 (_size_j_0 - 1)]) *
                           -2.0 -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 *
                                                     (_size_j_0 - 1) -
                                                 _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] -
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] +
                            _data_phi_2m1_30_1m1[_stride_phi_0 *
                                                     (_size_j_0 - 1) -
                                                 _stride_phi_0] -
                            _data_phi_2m1_30_1m1[_stride_phi_0 *
                                                 (_size_j_0 - 1)]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 *
                                                     (_size_j_0 - 1) -
                                                 _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 *
                                                (_size_j_0 - 1)] -
                            _data_phi_20_30_1m1[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] +
                            _data_phi_2m1_30_10[_stride_phi_0 *
                                                (_size_j_0 - 1)] +
                            _data_phi_2m1_30_1m1[_stride_phi_0 *
                                                     (_size_j_0 - 1) -
                                                 _stride_phi_0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_1m1[_stride_rho_0 *
                                                     (_size_j_0 - 1) -
                                                 _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 *
                                                (_size_j_0 - 1)] +
                            _data_phi_20_30_1m1[_stride_phi_0 *
                                                (_size_j_0 - 1)] -
                            _data_phi_2m1_30_10[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] +
                            _data_phi_2m1_30_1m1[_stride_phi_0 *
                                                     (_size_j_0 - 1) -
                                                 _stride_phi_0])) *
                      0.0833333333333333 * 1.7320508075688772 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_1 > 0 && _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_310 =
                      _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
                  double *RESTRICT _data_j_20_310_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_310;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_21_30 =
                      _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                  double *RESTRICT _data_rho_21_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30 =
                      _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                  double *RESTRICT _data_phi_21_30_1m1 =
                      _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_21_30;
                  double *RESTRICT _data_rho_20_30_1m1 =
                      _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_21_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_21_30;
                  _data_j_20_310_10[_stride_j_0 * (_size_j_0 - 1)] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           -2.0 +
                       kT *
                           (-3.0 * _data_rho_21_30_1m1[_stride_rho_0 *
                                                           (_size_j_0 - 1) -
                                                       _stride_rho_0] +
                            3.0 * _data_rho_20_30_10[_stride_rho_0 *
                                                     (_size_j_0 - 1)] -
                            _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0] +
                            _data_rho_20_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0] -
                            _data_rho_20_30_1m1[_stride_rho_0 *
                                                (_size_j_0 - 1)] +
                            _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0] -
                            _data_rho_21_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_1m1[_stride_rho_0 *
                                                (_size_j_0 - 1)]) *
                           2.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_21_30_1m1[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] +
                            _data_phi_21_30_1m1[_stride_phi_0 *
                                                (_size_j_0 - 1)]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] +
                            _data_phi_20_30_1m1[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] -
                            _data_phi_21_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_21_30_1m1[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_1m1[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_20_30_1m1[_stride_phi_0 *
                                                (_size_j_0 - 1)] +
                            _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] -
                            _data_phi_21_30_1m1[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0])) *
                      0.0833333333333333 * 1.7320508075688772 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (ctr_2 > 0 && _size_j_1 - ctr_1 > 1) {
                  double *RESTRICT _data_j_20_311 =
                      _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
                  double *RESTRICT _data_j_20_311_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_311;
                  double *RESTRICT _data_rho_2m1_30 =
                      _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                  double *RESTRICT _data_rho_2m1_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_rho_2m1_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_phi_2m1_30 =
                      _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                  double *RESTRICT _data_phi_2m1_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_2m1_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                  double *RESTRICT _data_phi_20_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                  _data_j_20_311_10[_stride_j_0 * (_size_j_0 - 1)] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_11[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_11[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           -2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_11[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           2.0 +
                       kT *
                           (-3.0 * _data_rho_20_30_10[_stride_rho_0 *
                                                      (_size_j_0 - 1)] +
                            3.0 * _data_rho_2m1_30_11[_stride_rho_0 *
                                                          (_size_j_0 - 1) -
                                                      _stride_rho_0] +
                            _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0] -
                            _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0] +
                            _data_rho_20_30_11[_stride_rho_0 *
                                               (_size_j_0 - 1)] -
                            _data_rho_2m1_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0] +
                            _data_rho_2m1_30_10[_stride_rho_0 *
                                                (_size_j_0 - 1)] -
                            _data_rho_2m1_30_11[_stride_rho_0 *
                                                (_size_j_0 - 1)]) *
                           -2.0 -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_11[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] -
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] +
                            _data_phi_2m1_30_11[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] -
                            _data_phi_2m1_30_11[_stride_phi_0 *
                                                (_size_j_0 - 1)]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_11[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 *
                                                (_size_j_0 - 1)] -
                            _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] +
                            _data_phi_2m1_30_10[_stride_phi_0 *
                                                (_size_j_0 - 1)] +
                            _data_phi_2m1_30_11[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0]) -
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_2m1_30_11[_stride_rho_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 *
                                                (_size_j_0 - 1)] +
                            _data_phi_20_30_11[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_2m1_30_10[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] +
                            _data_phi_2m1_30_11[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0])) *
                      0.0833333333333333 * 1.7320508075688772 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
                if (_size_j_1 - ctr_1 > 1 && _size_j_2 - ctr_2 > 1) {
                  double *RESTRICT _data_j_20_312 =
                      _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
                  double *RESTRICT _data_j_20_312_10 =
                      _stride_j_1 * ctr_1 + _data_j_20_312;
                  double *RESTRICT _data_rho_20_30 =
                      _data_rho + _stride_rho_2 * ctr_2;
                  double *RESTRICT _data_rho_20_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_21_30 =
                      _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                  double *RESTRICT _data_rho_21_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_21_30;
                  double *RESTRICT _data_phi_20_30 =
                      _data_phi + _stride_phi_2 * ctr_2;
                  double *RESTRICT _data_phi_20_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30 =
                      _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                  double *RESTRICT _data_phi_21_30_11 =
                      _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_21_30;
                  double *RESTRICT _data_phi_20_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_20_30;
                  double *RESTRICT _data_phi_21_30_10 =
                      _stride_phi_1 * ctr_1 + _data_phi_21_30;
                  double *RESTRICT _data_rho_20_30_11 =
                      _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                  double *RESTRICT _data_rho_21_30_10 =
                      _stride_rho_1 * ctr_1 + _data_rho_21_30;
                  _data_j_20_312_10[_stride_j_0 * (_size_j_0 - 1)] =
                      D *
                      (f_ext0 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           2.0 +
                       f_ext1 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           -2.0 +
                       f_ext2 * z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           -2.0 +
                       kT *
                           (-3.0 * _data_rho_20_30_10[_stride_rho_0 *
                                                      (_size_j_0 - 1)] +
                            3.0 * _data_rho_21_30_11[_stride_rho_0 *
                                                         (_size_j_0 - 1) -
                                                     _stride_rho_0] +
                            _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0] -
                            _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0] +
                            _data_rho_20_30_11[_stride_rho_0 *
                                               (_size_j_0 - 1)] -
                            _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0] +
                            _data_rho_21_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] -
                            _data_rho_21_30_11[_stride_rho_0 *
                                               (_size_j_0 - 1)]) *
                           -2.0 +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           (-_data_phi_20_30_10[_stride_phi_0 *
                                                    (_size_j_0 - 1) -
                                                _stride_phi_0] +
                            _data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] +
                            _data_phi_21_30_11[_stride_phi_0 *
                                               (_size_j_0 - 1)]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] +
                            _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] -
                            _data_phi_21_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0]) +
                       z *
                           (_data_rho_20_30_10[_stride_rho_0 *
                                               (_size_j_0 - 1)] +
                            _data_rho_21_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                           (_data_phi_20_30_10[_stride_phi_0 *
                                               (_size_j_0 - 1)] -
                            _data_phi_20_30_11[_stride_phi_0 *
                                               (_size_j_0 - 1)] +
                            _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] -
                            _data_phi_21_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0])) *
                      0.0833333333333333 * 1.7320508075688772 /
                      (kT * (1.0 + 2.8284271247461903 +
                             1.7320508075688772 * 1.33333333333333));
                }
              }
            }
          }
          {
            {
              if (_size_j_0 > 2 && _size_j_1 > 1 && ctr_2 > 0 &&
                  _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_31 =
                    _data_j + _stride_j_2 * ctr_2 + _stride_j_3;
                double *RESTRICT _data_j_20_31_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_31;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_20_30;
                _data_j_20_31_10[_stride_j_0] =
                    D *
                    (f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_20_30_1m1[_stride_rho_0]) +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0] -
                          _data_rho_20_30_1m1[_stride_rho_0]) *
                         2.0 +
                     z *
                         (_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_1m1[_stride_phi_0]) *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_20_30_1m1[_stride_rho_0])) *
                    0.5 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && ctr_2 > 0 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_33 =
                    _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
                double *RESTRICT _data_j_20_33_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_33;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                _data_j_20_33_10[_stride_j_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_20_30_1m1[0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_20_30_1m1[0]) *
                         2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0] -
                          _data_rho_20_30_1m1[0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_20_30_1m1[0]) *
                         (-_data_phi_20_30_10[0] +
                          _data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_1m1[0] +
                          _data_phi_20_30_1m1[_stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_20_30_1m1[0]) *
                         (_data_phi_20_30_10[0] +
                          _data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_1m1[0] -
                          _data_phi_20_30_1m1[_stride_phi_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_0 > 2 && _size_j_1 > 1 && ctr_2 > 0) {
                double *RESTRICT _data_j_20_37 =
                    _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
                double *RESTRICT _data_j_20_37_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_37;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_2m1_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                _data_j_20_37_10[_stride_j_0] =
                    D *
                    (f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0]) *
                         2.0 +
                     kT *
                         (-_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0]) *
                         -4.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_1m1[_stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_20_30_1m1[_stride_phi_0] -
                          _data_phi_2m1_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_0 > 2 && _size_j_1 > 1 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_38 =
                    _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
                double *RESTRICT _data_j_20_38_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_38;
                double *RESTRICT _data_rho_21_30 =
                    _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_21_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 =
                    _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21_30;
                _data_j_20_38_10[_stride_j_0] =
                    D *
                    (f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[_stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[_stride_rho_0]) *
                         -2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0] -
                          _data_rho_21_30_1m1[_stride_rho_0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[_stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_1m1[_stride_phi_0] +
                          _data_phi_21_30_10[_stride_phi_0] -
                          _data_phi_21_30_1m1[_stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[_stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_20_30_1m1[_stride_phi_0] -
                          _data_phi_21_30_10[_stride_phi_0] -
                          _data_phi_21_30_1m1[_stride_phi_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && ctr_2 > 0) {
                double *RESTRICT _data_j_20_39 =
                    _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
                double *RESTRICT _data_j_20_39_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_39;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_2m1_30;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_2m1_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                _data_j_20_39_10[_stride_j_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[0]) *
                         2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0] +
                          3.0 * _data_rho_2m1_30_1m1[0] +
                          _data_rho_20_30_10[0] - _data_rho_20_30_1m1[0] +
                          _data_rho_20_30_1m1[_stride_rho_0] -
                          _data_rho_2m1_30_10[0] +
                          _data_rho_2m1_30_10[_stride_rho_0] -
                          _data_rho_2m1_30_1m1[_stride_rho_0]) *
                         -2.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[0]) *
                         (_data_phi_20_30_10[0] -
                          _data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_1m1[0] -
                          _data_phi_2m1_30_1m1[_stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[0]) *
                         (-_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_1m1[0] +
                          _data_phi_2m1_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_1m1[0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[0]) *
                         (-_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_20_30_1m1[_stride_phi_0] -
                          _data_phi_2m1_30_10[0] + _data_phi_2m1_30_1m1[0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_310 =
                    _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
                double *RESTRICT _data_j_20_310_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_310;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30 =
                    _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_21_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 =
                    _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21_30;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_20_30;
                double *RESTRICT _data_rho_21_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_21_30;
                _data_j_20_310_10[_stride_j_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[0]) *
                         -2.0 +
                     kT *
                         (-3.0 * _data_rho_21_30_1m1[0] +
                          3.0 * _data_rho_20_30_10[_stride_rho_0] -
                          _data_rho_20_30_10[0] + _data_rho_20_30_1m1[0] -
                          _data_rho_20_30_1m1[_stride_rho_0] +
                          _data_rho_21_30_10[0] -
                          _data_rho_21_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[_stride_rho_0]) *
                         2.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[0]) *
                         (-_data_phi_20_30_10[0] +
                          _data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_21_30_1m1[0] +
                          _data_phi_21_30_1m1[_stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[0]) *
                         (_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_20_30_1m1[0] -
                          _data_phi_21_30_10[_stride_phi_0] -
                          _data_phi_21_30_1m1[0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_21_30_1m1[0]) *
                         (_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_1m1[_stride_phi_0] +
                          _data_phi_21_30_10[0] - _data_phi_21_30_1m1[0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
            }
            for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
              if (_size_j_1 > 1 && ctr_2 > 0 && _size_j_0 - ctr_0 > 1 &&
                  _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_31 =
                    _data_j + _stride_j_2 * ctr_2 + _stride_j_3;
                double *RESTRICT _data_j_20_31_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_31;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_20_30;
                _data_j_20_31_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0]) +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0]) *
                         2.0 +
                     z *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0]) *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0])) *
                    0.5 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && ctr_2 > 0 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_33 =
                    _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
                double *RESTRICT _data_j_20_33_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_33;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                _data_j_20_33_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && ctr_2 > 0 && _size_j_0 - ctr_0 > 1) {
                double *RESTRICT _data_j_20_37 =
                    _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
                double *RESTRICT _data_j_20_37_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_37;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_2m1_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                _data_j_20_37_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                         2.0 +
                     kT *
                         (-_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                         -4.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] -
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && _size_j_0 - ctr_0 > 1 &&
                  _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_38 =
                    _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
                double *RESTRICT _data_j_20_38_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_38;
                double *RESTRICT _data_rho_21_30 =
                    _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_21_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 =
                    _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21_30;
                _data_j_20_38_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                         -2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] +
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_1m1[_stride_phi_0 * ctr_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_1m1[_stride_phi_0 * ctr_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && ctr_2 > 0) {
                double *RESTRICT _data_j_20_39 =
                    _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
                double *RESTRICT _data_j_20_39_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_39;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_2m1_30;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_2m1_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                _data_j_20_39_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                         2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          3.0 * _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                                     _stride_rho_0] +
                          _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] -
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0] -
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                         -2.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] -
                          _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] -
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_310 =
                    _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
                double *RESTRICT _data_j_20_310_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_310;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30 =
                    _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_21_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 =
                    _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21_30;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_20_30;
                double *RESTRICT _data_rho_21_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_21_30;
                _data_j_20_310_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         -2.0 +
                     kT *
                         (-3.0 * _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                                     _stride_rho_0] +
                          3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0] -
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] -
                          _data_rho_21_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0]) *
                         2.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_1m1[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_21_30_1m1[_stride_phi_0 * ctr_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] -
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_21_30_1m1[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_21_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] +
                          _data_phi_21_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_21_30_1m1[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
            }
            {
              if (_size_j_1 > 1 && ctr_2 > 0 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_33 =
                    _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
                double *RESTRICT _data_j_20_33_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_33;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                _data_j_20_33_10[_stride_j_0 * (_size_j_0 - 1)] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         2.0 +
                     kT *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         4.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_20_30_1m1[_stride_phi_0 *
                                              (_size_j_0 - 1)]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 *
                                              (_size_j_0 - 1)])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && ctr_2 > 0) {
                double *RESTRICT _data_j_20_39 =
                    _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
                double *RESTRICT _data_j_20_39_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_39;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_2m1_30;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_2m1_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                _data_j_20_39_10[_stride_j_0 * (_size_j_0 - 1)] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                         2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1)] +
                          3.0 * _data_rho_2m1_30_1m1[_stride_rho_0 *
                                                         (_size_j_0 - 1) -
                                                     _stride_rho_0] +
                          _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] -
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_2m1_30_1m1[_stride_rho_0 *
                                               (_size_j_0 - 1)]) *
                         -2.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] -
                          _data_phi_2m1_30_1m1[_stride_phi_0 *
                                               (_size_j_0 - 1)]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_1 > 1 && _size_j_2 - ctr_2 > 1) {
                double *RESTRICT _data_j_20_310 =
                    _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
                double *RESTRICT _data_j_20_310_10 =
                    _stride_j_1 * (_size_j_1 - 1) + _data_j_20_310;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * ctr_2;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
                double *RESTRICT _data_rho_21_30 =
                    _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                double *RESTRICT _data_rho_21_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_21_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * ctr_2;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_20_30;
                double *RESTRICT _data_phi_21_30 =
                    _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                double *RESTRICT _data_phi_21_30_1m1 =
                    _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                    _data_phi_21_30;
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
                double *RESTRICT _data_phi_21_30_10 =
                    _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21_30;
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                    _data_rho_20_30;
                double *RESTRICT _data_rho_21_30_10 =
                    _stride_rho_1 * (_size_j_1 - 1) + _data_rho_21_30;
                _data_j_20_310_10[_stride_j_0 * (_size_j_0 - 1)] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         -2.0 +
                     kT *
                         (-3.0 * _data_rho_21_30_1m1[_stride_rho_0 *
                                                         (_size_j_0 - 1) -
                                                     _stride_rho_0] +
                          3.0 * _data_rho_20_30_10[_stride_rho_0 *
                                                   (_size_j_0 - 1)] -
                          _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0] -
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] -
                          _data_rho_21_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 *
                                              (_size_j_0 - 1)]) *
                         2.0 +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_21_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_21_30_1m1[_stride_phi_0 *
                                              (_size_j_0 - 1)]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] -
                          _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_21_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0]) +
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_21_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_21_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_21_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
            }
          }
        }
      }
      {
        {
          if (_size_j_1 > 1 && _size_j_2 > 1) {
            double *RESTRICT _data_j_20_311 =
                _data_j + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3;
            double *RESTRICT _data_j_20_311_10 = _data_j_20_311;
            double *RESTRICT _data_rho_2m1_30 =
                _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
            double *RESTRICT _data_rho_2m1_30_11 =
                _stride_rho_1 + _data_rho_2m1_30;
            double *RESTRICT _data_rho_2m1_30_10 = _data_rho_2m1_30;
            double *RESTRICT _data_rho_20_30 =
                _data_rho + _stride_rho_2 * (_size_j_2 - 1);
            double *RESTRICT _data_rho_20_30_11 =
                _stride_rho_1 + _data_rho_20_30;
            double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
            double *RESTRICT _data_phi_2m1_30 =
                _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
            double *RESTRICT _data_phi_2m1_30_11 =
                _stride_phi_1 + _data_phi_2m1_30;
            double *RESTRICT _data_phi_20_30 =
                _data_phi + _stride_phi_2 * (_size_j_2 - 1);
            double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
            double *RESTRICT _data_phi_2m1_30_10 = _data_phi_2m1_30;
            double *RESTRICT _data_phi_20_30_11 =
                _stride_phi_1 + _data_phi_20_30;
            _data_j_20_311_10[_stride_j_0] =
                D *
                (f_ext0 * z *
                     (_data_rho_20_30_10[_stride_rho_0] +
                      _data_rho_2m1_30_11[0]) *
                     2.0 +
                 f_ext1 * z *
                     (_data_rho_20_30_10[_stride_rho_0] +
                      _data_rho_2m1_30_11[0]) *
                     -2.0 +
                 f_ext2 * z *
                     (_data_rho_20_30_10[_stride_rho_0] +
                      _data_rho_2m1_30_11[0]) *
                     2.0 +
                 kT *
                     (-3.0 * _data_rho_20_30_10[_stride_rho_0] +
                      3.0 * _data_rho_2m1_30_11[0] + _data_rho_20_30_10[0] -
                      _data_rho_20_30_11[0] +
                      _data_rho_20_30_11[_stride_rho_0] -
                      _data_rho_2m1_30_10[0] +
                      _data_rho_2m1_30_10[_stride_rho_0] -
                      _data_rho_2m1_30_11[_stride_rho_0]) *
                     -2.0 -
                 z *
                     (_data_rho_20_30_10[_stride_rho_0] +
                      _data_rho_2m1_30_11[0]) *
                     (_data_phi_20_30_10[0] -
                      _data_phi_20_30_10[_stride_phi_0] +
                      _data_phi_2m1_30_11[0] -
                      _data_phi_2m1_30_11[_stride_phi_0]) -
                 z *
                     (_data_rho_20_30_10[_stride_rho_0] +
                      _data_rho_2m1_30_11[0]) *
                     (-_data_phi_20_30_10[_stride_phi_0] -
                      _data_phi_20_30_11[0] +
                      _data_phi_2m1_30_10[_stride_phi_0] +
                      _data_phi_2m1_30_11[0]) -
                 z *
                     (_data_rho_20_30_10[_stride_rho_0] +
                      _data_rho_2m1_30_11[0]) *
                     (-_data_phi_20_30_10[_stride_phi_0] +
                      _data_phi_20_30_11[_stride_phi_0] -
                      _data_phi_2m1_30_10[0] + _data_phi_2m1_30_11[0])) *
                0.0833333333333333 * 1.7320508075688772 /
                (kT * (1.0 + 2.8284271247461903 +
                       1.7320508075688772 * 1.33333333333333));
          }
#pragma omp for schedule(static)
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (_size_j_1 > 1 && _size_j_2 > 1) {
              double *RESTRICT _data_j_20_311 =
                  _data_j + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3;
              double *RESTRICT _data_j_20_311_10 = _data_j_20_311;
              double *RESTRICT _data_rho_2m1_30 =
                  _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              double *RESTRICT _data_rho_2m1_30_11 =
                  _stride_rho_1 + _data_rho_2m1_30;
              double *RESTRICT _data_rho_2m1_30_10 = _data_rho_2m1_30;
              double *RESTRICT _data_rho_20_30 =
                  _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              double *RESTRICT _data_rho_20_30_11 =
                  _stride_rho_1 + _data_rho_20_30;
              double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
              double *RESTRICT _data_phi_2m1_30 =
                  _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              double *RESTRICT _data_phi_2m1_30_11 =
                  _stride_phi_1 + _data_phi_2m1_30;
              double *RESTRICT _data_phi_20_30 =
                  _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
              double *RESTRICT _data_phi_2m1_30_10 = _data_phi_2m1_30;
              double *RESTRICT _data_phi_20_30_11 =
                  _stride_phi_1 + _data_phi_20_30;
              _data_j_20_311_10[_stride_j_0 * ctr_0] =
                  D *
                  (f_ext0 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0]) *
                       2.0 +
                   f_ext1 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0]) *
                       -2.0 +
                   f_ext2 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0]) *
                       2.0 +
                   kT *
                       (-3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        3.0 * _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                                  _stride_rho_0] +
                        _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0] -
                        _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0] +
                        _data_rho_20_30_11[_stride_rho_0 * ctr_0] -
                        _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0] +
                        _data_rho_2m1_30_10[_stride_rho_0 * ctr_0] -
                        _data_rho_2m1_30_11[_stride_rho_0 * ctr_0]) *
                       -2.0 -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0]) *
                       (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                           _stride_phi_0] -
                        _data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                        _data_phi_2m1_30_11[_stride_phi_0 * ctr_0 -
                                            _stride_phi_0] -
                        _data_phi_2m1_30_11[_stride_phi_0 * ctr_0]) -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0]) *
                       (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                        _data_phi_20_30_11[_stride_phi_0 * ctr_0 -
                                           _stride_phi_0] +
                        _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                        _data_phi_2m1_30_11[_stride_phi_0 * ctr_0 -
                                            _stride_phi_0]) -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0]) *
                       (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                        _data_phi_20_30_11[_stride_phi_0 * ctr_0] -
                        _data_phi_2m1_30_10[_stride_phi_0 * ctr_0 -
                                            _stride_phi_0] +
                        _data_phi_2m1_30_11[_stride_phi_0 * ctr_0 -
                                            _stride_phi_0])) *
                  0.0833333333333333 * 1.7320508075688772 /
                  (kT * (1.0 + 2.8284271247461903 +
                         1.7320508075688772 * 1.33333333333333));
            }
          }
          if (_size_j_1 > 1 && _size_j_2 > 1) {
            double *RESTRICT _data_j_20_311 =
                _data_j + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3;
            double *RESTRICT _data_j_20_311_10 = _data_j_20_311;
            double *RESTRICT _data_rho_2m1_30 =
                _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
            double *RESTRICT _data_rho_2m1_30_11 =
                _stride_rho_1 + _data_rho_2m1_30;
            double *RESTRICT _data_rho_2m1_30_10 = _data_rho_2m1_30;
            double *RESTRICT _data_rho_20_30 =
                _data_rho + _stride_rho_2 * (_size_j_2 - 1);
            double *RESTRICT _data_rho_20_30_11 =
                _stride_rho_1 + _data_rho_20_30;
            double *RESTRICT _data_rho_20_30_10 = _data_rho_20_30;
            double *RESTRICT _data_phi_2m1_30 =
                _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
            double *RESTRICT _data_phi_2m1_30_11 =
                _stride_phi_1 + _data_phi_2m1_30;
            double *RESTRICT _data_phi_20_30 =
                _data_phi + _stride_phi_2 * (_size_j_2 - 1);
            double *RESTRICT _data_phi_20_30_10 = _data_phi_20_30;
            double *RESTRICT _data_phi_2m1_30_10 = _data_phi_2m1_30;
            double *RESTRICT _data_phi_20_30_11 =
                _stride_phi_1 + _data_phi_20_30;
            _data_j_20_311_10[_stride_j_0 * (_size_j_0 - 1)] =
                D *
                (f_ext0 * z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0]) *
                     2.0 +
                 f_ext1 * z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0]) *
                     -2.0 +
                 f_ext2 * z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0]) *
                     2.0 +
                 kT *
                     (-3.0 *
                          _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      3.0 *
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0] +
                      _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0] -
                      _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0] +
                      _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1)] -
                      _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0] +
                      _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                      _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1)]) *
                     -2.0 -
                 z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0]) *
                     (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                         _stride_phi_0] -
                      _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                      _data_phi_2m1_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                          _stride_phi_0] -
                      _data_phi_2m1_30_11[_stride_phi_0 * (_size_j_0 - 1)]) -
                 z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0]) *
                     (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                      _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                         _stride_phi_0] +
                      _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                      _data_phi_2m1_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                          _stride_phi_0]) -
                 z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0]) *
                     (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                      _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1)] -
                      _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                          _stride_phi_0] +
                      _data_phi_2m1_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                          _stride_phi_0])) *
                0.0833333333333333 * 1.7320508075688772 /
                (kT * (1.0 + 2.8284271247461903 +
                       1.7320508075688772 * 1.33333333333333));
          }
        }
#pragma omp for schedule(static)
        for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
          {
            {
              if (_size_j_0 > 2 && _size_j_2 > 1 && ctr_1 > 0 &&
                  _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_32 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 2 * _stride_j_3;
                double *RESTRICT _data_j_20_32_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_32;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                _data_j_20_32_10[_stride_j_0] =
                    D *
                    (f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_10[_stride_rho_0]) +
                     kT *
                         (-_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_10[_stride_rho_0]) *
                         -2.0 -
                     z *
                         (-_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0]) *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_10[_stride_rho_0])) *
                    0.5 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && ctr_1 > 0 && _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_35 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 5 * _stride_j_3;
                double *RESTRICT _data_j_20_35_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_35;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                _data_j_20_35_10[_stride_j_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_10[0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_10[0]) *
                         2.0 +
                     kT *
                         (-_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_10[0]) *
                         -4.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_10[0]) *
                         (-_data_phi_20_30_10[0] -
                          _data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_10[0] +
                          _data_phi_2m1_30_10[_stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_10[0]) *
                         (_data_phi_20_30_10[0] -
                          _data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_10[0] -
                          _data_phi_2m1_30_10[_stride_phi_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_0 > 2 && _size_j_2 > 1 && ctr_1 > 0) {
                double *RESTRICT _data_j_20_37 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3;
                double *RESTRICT _data_j_20_37_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_37;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                _data_j_20_37_10[_stride_j_0] =
                    D *
                    (f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0]) *
                         2.0 +
                     kT *
                         (-_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0]) *
                         -4.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_1m1[_stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_20_30_1m1[_stride_phi_0] -
                          _data_phi_2m1_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && ctr_1 > 0) {
                double *RESTRICT _data_j_20_39 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3;
                double *RESTRICT _data_j_20_39_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_39;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                _data_j_20_39_10[_stride_j_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[0]) *
                         2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0] +
                          3.0 * _data_rho_2m1_30_1m1[0] +
                          _data_rho_20_30_10[0] - _data_rho_20_30_1m1[0] +
                          _data_rho_20_30_1m1[_stride_rho_0] -
                          _data_rho_2m1_30_10[0] +
                          _data_rho_2m1_30_10[_stride_rho_0] -
                          _data_rho_2m1_30_1m1[_stride_rho_0]) *
                         -2.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[0]) *
                         (_data_phi_20_30_10[0] -
                          _data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_1m1[0] -
                          _data_phi_2m1_30_1m1[_stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[0]) *
                         (-_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_1m1[0] +
                          _data_phi_2m1_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_1m1[0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_1m1[0]) *
                         (-_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_20_30_1m1[_stride_phi_0] -
                          _data_phi_2m1_30_10[0] + _data_phi_2m1_30_1m1[0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_311 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3;
                double *RESTRICT _data_j_20_311_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_311;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_11 =
                    _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_11 =
                    _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                _data_j_20_311_10[_stride_j_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_11[0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_11[0]) *
                         -2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_11[0]) *
                         2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0] +
                          3.0 * _data_rho_2m1_30_11[0] + _data_rho_20_30_10[0] -
                          _data_rho_20_30_11[0] +
                          _data_rho_20_30_11[_stride_rho_0] -
                          _data_rho_2m1_30_10[0] +
                          _data_rho_2m1_30_10[_stride_rho_0] -
                          _data_rho_2m1_30_11[_stride_rho_0]) *
                         -2.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_11[0]) *
                         (_data_phi_20_30_10[0] -
                          _data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_11[0] -
                          _data_phi_2m1_30_11[_stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_11[0]) *
                         (-_data_phi_20_30_10[_stride_phi_0] -
                          _data_phi_20_30_11[0] +
                          _data_phi_2m1_30_10[_stride_phi_0] +
                          _data_phi_2m1_30_11[0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0] +
                          _data_rho_2m1_30_11[0]) *
                         (-_data_phi_20_30_10[_stride_phi_0] +
                          _data_phi_20_30_11[_stride_phi_0] -
                          _data_phi_2m1_30_10[0] + _data_phi_2m1_30_11[0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
            }
            for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
              if (_size_j_2 > 1 && ctr_1 > 0 && _size_j_0 - ctr_0 > 1 &&
                  _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_32 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 2 * _stride_j_3;
                double *RESTRICT _data_j_20_32_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_32;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                _data_j_20_32_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0]) +
                     kT *
                         (-_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0]) *
                         -2.0 -
                     z *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0]) *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0])) *
                    0.5 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && ctr_1 > 0 && _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_35 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 5 * _stride_j_3;
                double *RESTRICT _data_j_20_35_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_35;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                _data_j_20_35_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         2.0 +
                     kT *
                         (-_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         -4.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] -
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] -
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && ctr_1 > 0 && _size_j_0 - ctr_0 > 1) {
                double *RESTRICT _data_j_20_37 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3;
                double *RESTRICT _data_j_20_37_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_37;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                _data_j_20_37_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                         2.0 +
                     kT *
                         (-_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                         -4.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] -
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && ctr_1 > 0) {
                double *RESTRICT _data_j_20_39 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3;
                double *RESTRICT _data_j_20_39_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_39;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                _data_j_20_39_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                         2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          3.0 * _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                                     _stride_rho_0] +
                          _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] -
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * ctr_0] -
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                         -2.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0] -
                          _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                               _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] -
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0 -
                                               _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_311 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3;
                double *RESTRICT _data_j_20_311_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_311;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_11 =
                    _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_11 =
                    _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                _data_j_20_311_10[_stride_j_0 * ctr_0] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         -2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          3.0 * _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                                    _stride_rho_0] +
                          _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] -
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0] +
                          _data_rho_20_30_11[_stride_rho_0 * ctr_0] -
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * ctr_0] -
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0]) *
                         -2.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                          _data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_11[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] -
                          _data_phi_2m1_30_11[_stride_phi_0 * ctr_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                          _data_phi_20_30_11[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_2m1_30_11[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                          _data_rho_2m1_30_11[_stride_rho_0 * ctr_0 -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                          _data_phi_20_30_11[_stride_phi_0 * ctr_0] -
                          _data_phi_2m1_30_10[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_11[_stride_phi_0 * ctr_0 -
                                              _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
            }
            {
              if (_size_j_2 > 1 && ctr_1 > 0 && _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_35 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 5 * _stride_j_3;
                double *RESTRICT _data_j_20_35_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_35;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                _data_j_20_35_10[_stride_j_0 * (_size_j_0 - 1)] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         2.0 +
                     kT *
                         (-_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         -4.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] -
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 *
                                              (_size_j_0 - 1)]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] -
                          _data_phi_2m1_30_10[_stride_phi_0 *
                                              (_size_j_0 - 1)])) *
                    0.125 * 1.4142135623730951 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && ctr_1 > 0) {
                double *RESTRICT _data_j_20_39 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3;
                double *RESTRICT _data_j_20_39_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_39;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_1m1 =
                    _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_1m1 =
                    _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20_30;
                _data_j_20_39_10[_stride_j_0 * (_size_j_0 - 1)] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                         2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                         2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1)] +
                          3.0 * _data_rho_2m1_30_1m1[_stride_rho_0 *
                                                         (_size_j_0 - 1) -
                                                     _stride_rho_0] +
                          _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] -
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0] +
                          _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_2m1_30_1m1[_stride_rho_0 *
                                               (_size_j_0 - 1)]) *
                         -2.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0] -
                          _data_phi_2m1_30_1m1[_stride_phi_0 *
                                               (_size_j_0 - 1)]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                               _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
              if (_size_j_2 > 1 && _size_j_1 - ctr_1 > 1) {
                double *RESTRICT _data_j_20_311 =
                    _data_j + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3;
                double *RESTRICT _data_j_20_311_10 =
                    _stride_j_1 * ctr_1 + _data_j_20_311;
                double *RESTRICT _data_rho_2m1_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
                double *RESTRICT _data_rho_2m1_30_11 =
                    _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_2m1_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_2m1_30;
                double *RESTRICT _data_rho_20_30 =
                    _data_rho + _stride_rho_2 * (_size_j_2 - 1);
                double *RESTRICT _data_rho_20_30_11 =
                    _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20_30;
                double *RESTRICT _data_rho_20_30_10 =
                    _stride_rho_1 * ctr_1 + _data_rho_20_30;
                double *RESTRICT _data_phi_2m1_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
                double *RESTRICT _data_phi_2m1_30_11 =
                    _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30 =
                    _data_phi + _stride_phi_2 * (_size_j_2 - 1);
                double *RESTRICT _data_phi_20_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_20_30;
                double *RESTRICT _data_phi_2m1_30_10 =
                    _stride_phi_1 * ctr_1 + _data_phi_2m1_30;
                double *RESTRICT _data_phi_20_30_11 =
                    _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20_30;
                _data_j_20_311_10[_stride_j_0 * (_size_j_0 - 1)] =
                    D *
                    (f_ext0 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         2.0 +
                     f_ext1 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         -2.0 +
                     f_ext2 * z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         2.0 +
                     kT *
                         (-3.0 * _data_rho_20_30_10[_stride_rho_0 *
                                                    (_size_j_0 - 1)] +
                          3.0 * _data_rho_2m1_30_11[_stride_rho_0 *
                                                        (_size_j_0 - 1) -
                                                    _stride_rho_0] +
                          _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] -
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                             _stride_rho_0] +
                          _data_rho_20_30_11[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0] +
                          _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                          _data_rho_2m1_30_11[_stride_rho_0 *
                                              (_size_j_0 - 1)]) *
                         -2.0 -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] -
                          _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_2m1_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] -
                          _data_phi_2m1_30_11[_stride_phi_0 *
                                              (_size_j_0 - 1)]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                             _stride_phi_0] +
                          _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_2m1_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0]) -
                     z *
                         (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                          _data_rho_2m1_30_11[_stride_rho_0 * (_size_j_0 - 1) -
                                              _stride_rho_0]) *
                         (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                          _data_phi_20_30_11[_stride_phi_0 * (_size_j_0 - 1)] -
                          _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0] +
                          _data_phi_2m1_30_11[_stride_phi_0 * (_size_j_0 - 1) -
                                              _stride_phi_0])) *
                    0.0833333333333333 * 1.7320508075688772 /
                    (kT * (1.0 + 2.8284271247461903 +
                           1.7320508075688772 * 1.33333333333333));
              }
            }
          }
        }
        {
          {
            if (_size_j_0 > 2 && _size_j_1 > 1 && _size_j_2 > 1) {
              double *RESTRICT _data_j_20_37 =
                  _data_j + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3;
              double *RESTRICT _data_j_20_37_10 =
                  _stride_j_1 * (_size_j_1 - 1) + _data_j_20_37;
              double *RESTRICT _data_rho_20_30 =
                  _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              double *RESTRICT _data_rho_20_30_10 =
                  _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
              double *RESTRICT _data_rho_2m1_30 =
                  _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              double *RESTRICT _data_rho_2m1_30_1m1 =
                  _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                  _data_rho_2m1_30;
              double *RESTRICT _data_phi_2m1_30 =
                  _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              double *RESTRICT _data_phi_2m1_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1_30;
              double *RESTRICT _data_phi_20_30 =
                  _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              double *RESTRICT _data_phi_20_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
              double *RESTRICT _data_phi_2m1_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_2m1_30;
              double *RESTRICT _data_phi_20_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_20_30;
              _data_j_20_37_10[_stride_j_0] =
                  D *
                  (f_ext1 * z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0]) *
                       2.0 +
                   f_ext2 * z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0]) *
                       2.0 +
                   kT *
                       (-_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0]) *
                       -4.0 -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0]) *
                       (-_data_phi_20_30_10[_stride_phi_0] -
                        _data_phi_20_30_1m1[_stride_phi_0] +
                        _data_phi_2m1_30_10[_stride_phi_0] +
                        _data_phi_2m1_30_1m1[_stride_phi_0]) -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0]) *
                       (-_data_phi_20_30_10[_stride_phi_0] +
                        _data_phi_20_30_1m1[_stride_phi_0] -
                        _data_phi_2m1_30_10[_stride_phi_0] +
                        _data_phi_2m1_30_1m1[_stride_phi_0])) *
                  0.125 * 1.4142135623730951 /
                  (kT * (1.0 + 2.8284271247461903 +
                         1.7320508075688772 * 1.33333333333333));
            }
            if (_size_j_1 > 1 && _size_j_2 > 1) {
              double *RESTRICT _data_j_20_39 =
                  _data_j + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3;
              double *RESTRICT _data_j_20_39_10 =
                  _stride_j_1 * (_size_j_1 - 1) + _data_j_20_39;
              double *RESTRICT _data_rho_2m1_30 =
                  _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              double *RESTRICT _data_rho_2m1_30_1m1 =
                  _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                  _data_rho_2m1_30;
              double *RESTRICT _data_rho_2m1_30_10 =
                  _stride_rho_1 * (_size_j_1 - 1) + _data_rho_2m1_30;
              double *RESTRICT _data_rho_20_30 =
                  _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              double *RESTRICT _data_rho_20_30_1m1 =
                  _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                  _data_rho_20_30;
              double *RESTRICT _data_rho_20_30_10 =
                  _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
              double *RESTRICT _data_phi_2m1_30 =
                  _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              double *RESTRICT _data_phi_2m1_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_2m1_30;
              double *RESTRICT _data_phi_20_30 =
                  _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              double *RESTRICT _data_phi_20_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
              double *RESTRICT _data_phi_2m1_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1_30;
              double *RESTRICT _data_phi_20_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_20_30;
              _data_j_20_39_10[_stride_j_0] =
                  D *
                  (f_ext0 * z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_2m1_30_1m1[0]) *
                       2.0 +
                   f_ext1 * z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_2m1_30_1m1[0]) *
                       2.0 +
                   f_ext2 * z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_2m1_30_1m1[0]) *
                       2.0 +
                   kT *
                       (-3.0 * _data_rho_20_30_10[_stride_rho_0] +
                        3.0 * _data_rho_2m1_30_1m1[0] + _data_rho_20_30_10[0] -
                        _data_rho_20_30_1m1[0] +
                        _data_rho_20_30_1m1[_stride_rho_0] -
                        _data_rho_2m1_30_10[0] +
                        _data_rho_2m1_30_10[_stride_rho_0] -
                        _data_rho_2m1_30_1m1[_stride_rho_0]) *
                       -2.0 -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_2m1_30_1m1[0]) *
                       (_data_phi_20_30_10[0] -
                        _data_phi_20_30_10[_stride_phi_0] +
                        _data_phi_2m1_30_1m1[0] -
                        _data_phi_2m1_30_1m1[_stride_phi_0]) -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_2m1_30_1m1[0]) *
                       (-_data_phi_20_30_10[_stride_phi_0] -
                        _data_phi_20_30_1m1[0] +
                        _data_phi_2m1_30_10[_stride_phi_0] +
                        _data_phi_2m1_30_1m1[0]) -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0] +
                        _data_rho_2m1_30_1m1[0]) *
                       (-_data_phi_20_30_10[_stride_phi_0] +
                        _data_phi_20_30_1m1[_stride_phi_0] -
                        _data_phi_2m1_30_10[0] + _data_phi_2m1_30_1m1[0])) *
                  0.0833333333333333 * 1.7320508075688772 /
                  (kT * (1.0 + 2.8284271247461903 +
                         1.7320508075688772 * 1.33333333333333));
            }
          }
#pragma omp for schedule(static)
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (_size_j_1 > 1 && _size_j_2 > 1 && _size_j_0 - ctr_0 > 1) {
              double *RESTRICT _data_j_20_37 =
                  _data_j + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3;
              double *RESTRICT _data_j_20_37_10 =
                  _stride_j_1 * (_size_j_1 - 1) + _data_j_20_37;
              double *RESTRICT _data_rho_20_30 =
                  _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              double *RESTRICT _data_rho_20_30_10 =
                  _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
              double *RESTRICT _data_rho_2m1_30 =
                  _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              double *RESTRICT _data_rho_2m1_30_1m1 =
                  _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                  _data_rho_2m1_30;
              double *RESTRICT _data_phi_2m1_30 =
                  _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              double *RESTRICT _data_phi_2m1_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1_30;
              double *RESTRICT _data_phi_20_30 =
                  _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              double *RESTRICT _data_phi_20_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
              double *RESTRICT _data_phi_2m1_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_2m1_30;
              double *RESTRICT _data_phi_20_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_20_30;
              _data_j_20_37_10[_stride_j_0 * ctr_0] =
                  D *
                  (f_ext1 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                       2.0 +
                   f_ext2 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                       2.0 +
                   kT *
                       (-_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                       -4.0 -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                       (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                        _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] +
                        _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                        _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0]) -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                       (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                        _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] -
                        _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                        _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0])) *
                  0.125 * 1.4142135623730951 /
                  (kT * (1.0 + 2.8284271247461903 +
                         1.7320508075688772 * 1.33333333333333));
            }
            if (_size_j_1 > 1 && _size_j_2 > 1) {
              double *RESTRICT _data_j_20_39 =
                  _data_j + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3;
              double *RESTRICT _data_j_20_39_10 =
                  _stride_j_1 * (_size_j_1 - 1) + _data_j_20_39;
              double *RESTRICT _data_rho_2m1_30 =
                  _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              double *RESTRICT _data_rho_2m1_30_1m1 =
                  _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                  _data_rho_2m1_30;
              double *RESTRICT _data_rho_2m1_30_10 =
                  _stride_rho_1 * (_size_j_1 - 1) + _data_rho_2m1_30;
              double *RESTRICT _data_rho_20_30 =
                  _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              double *RESTRICT _data_rho_20_30_1m1 =
                  _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                  _data_rho_20_30;
              double *RESTRICT _data_rho_20_30_10 =
                  _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
              double *RESTRICT _data_phi_2m1_30 =
                  _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              double *RESTRICT _data_phi_2m1_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_2m1_30;
              double *RESTRICT _data_phi_20_30 =
                  _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              double *RESTRICT _data_phi_20_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
              double *RESTRICT _data_phi_2m1_30_10 =
                  _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1_30;
              double *RESTRICT _data_phi_20_30_1m1 =
                  _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                  _data_phi_20_30;
              _data_j_20_39_10[_stride_j_0 * ctr_0] =
                  D *
                  (f_ext0 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                       2.0 +
                   f_ext1 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                       2.0 +
                   f_ext2 * z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                       2.0 +
                   kT *
                       (-3.0 * _data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        3.0 * _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                                   _stride_rho_0] +
                        _data_rho_20_30_10[_stride_rho_0 * ctr_0 -
                                           _stride_rho_0] -
                        _data_rho_20_30_1m1[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0] +
                        _data_rho_20_30_1m1[_stride_rho_0 * ctr_0] -
                        _data_rho_2m1_30_10[_stride_rho_0 * ctr_0 -
                                            _stride_rho_0] +
                        _data_rho_2m1_30_10[_stride_rho_0 * ctr_0] -
                        _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0]) *
                       -2.0 -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                       (_data_phi_20_30_10[_stride_phi_0 * ctr_0 -
                                           _stride_phi_0] -
                        _data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                        _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0] -
                        _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0]) -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                       (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] -
                        _data_phi_20_30_1m1[_stride_phi_0 * ctr_0 -
                                            _stride_phi_0] +
                        _data_phi_2m1_30_10[_stride_phi_0 * ctr_0] +
                        _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0]) -
                   z *
                       (_data_rho_20_30_10[_stride_rho_0 * ctr_0] +
                        _data_rho_2m1_30_1m1[_stride_rho_0 * ctr_0 -
                                             _stride_rho_0]) *
                       (-_data_phi_20_30_10[_stride_phi_0 * ctr_0] +
                        _data_phi_20_30_1m1[_stride_phi_0 * ctr_0] -
                        _data_phi_2m1_30_10[_stride_phi_0 * ctr_0 -
                                            _stride_phi_0] +
                        _data_phi_2m1_30_1m1[_stride_phi_0 * ctr_0 -
                                             _stride_phi_0])) *
                  0.0833333333333333 * 1.7320508075688772 /
                  (kT * (1.0 + 2.8284271247461903 +
                         1.7320508075688772 * 1.33333333333333));
            }
          }
          if (_size_j_1 > 1 && _size_j_2 > 1) {
            double *RESTRICT _data_j_20_39 =
                _data_j + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3;
            double *RESTRICT _data_j_20_39_10 =
                _stride_j_1 * (_size_j_1 - 1) + _data_j_20_39;
            double *RESTRICT _data_rho_2m1_30 =
                _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
            double *RESTRICT _data_rho_2m1_30_1m1 =
                _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                _data_rho_2m1_30;
            double *RESTRICT _data_rho_2m1_30_10 =
                _stride_rho_1 * (_size_j_1 - 1) + _data_rho_2m1_30;
            double *RESTRICT _data_rho_20_30 =
                _data_rho + _stride_rho_2 * (_size_j_2 - 1);
            double *RESTRICT _data_rho_20_30_1m1 =
                _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 +
                _data_rho_20_30;
            double *RESTRICT _data_rho_20_30_10 =
                _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20_30;
            double *RESTRICT _data_phi_2m1_30 =
                _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
            double *RESTRICT _data_phi_2m1_30_1m1 =
                _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                _data_phi_2m1_30;
            double *RESTRICT _data_phi_20_30 =
                _data_phi + _stride_phi_2 * (_size_j_2 - 1);
            double *RESTRICT _data_phi_20_30_10 =
                _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20_30;
            double *RESTRICT _data_phi_2m1_30_10 =
                _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1_30;
            double *RESTRICT _data_phi_20_30_1m1 =
                _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 +
                _data_phi_20_30;
            _data_j_20_39_10[_stride_j_0 * (_size_j_0 - 1)] =
                D *
                (f_ext0 * z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                           _stride_rho_0]) *
                     2.0 +
                 f_ext1 * z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                           _stride_rho_0]) *
                     2.0 +
                 f_ext2 * z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                           _stride_rho_0]) *
                     2.0 +
                 kT *
                     (-3.0 *
                          _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      3.0 *
                          _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                               _stride_rho_0] +
                      _data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                         _stride_rho_0] -
                      _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0] +
                      _data_rho_20_30_1m1[_stride_rho_0 * (_size_j_0 - 1)] -
                      _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1) -
                                          _stride_rho_0] +
                      _data_rho_2m1_30_10[_stride_rho_0 * (_size_j_0 - 1)] -
                      _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1)]) *
                     -2.0 -
                 z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                           _stride_rho_0]) *
                     (_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                         _stride_phi_0] -
                      _data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                      _data_phi_2m1_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                           _stride_phi_0] -
                      _data_phi_2m1_30_1m1[_stride_phi_0 * (_size_j_0 - 1)]) -
                 z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                           _stride_rho_0]) *
                     (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] -
                      _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                          _stride_phi_0] +
                      _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                      _data_phi_2m1_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                           _stride_phi_0]) -
                 z *
                     (_data_rho_20_30_10[_stride_rho_0 * (_size_j_0 - 1)] +
                      _data_rho_2m1_30_1m1[_stride_rho_0 * (_size_j_0 - 1) -
                                           _stride_rho_0]) *
                     (-_data_phi_20_30_10[_stride_phi_0 * (_size_j_0 - 1)] +
                      _data_phi_20_30_1m1[_stride_phi_0 * (_size_j_0 - 1)] -
                      _data_phi_2m1_30_10[_stride_phi_0 * (_size_j_0 - 1) -
                                          _stride_phi_0] +
                      _data_phi_2m1_30_1m1[_stride_phi_0 * (_size_j_0 - 1) -
                                           _stride_phi_0])) *
                0.0833333333333333 * 1.7320508075688772 /
                (kT * (1.0 + 2.8284271247461903 +
                       1.7320508075688772 * 1.33333333333333));
          }
        }
      }
    }
  }
}
} // namespace
  // internal_diffusivefluxkernelwithelectrostatic_diffusivefluxkernelwithelectrostatic

void DiffusiveFluxKernelWithElectrostatic::run(IBlock *block) {
  auto rho = block->getData<field::GhostLayerField<double, 1>>(rhoID);
  auto j = block->getData<field::GhostLayerField<double, 13>>(jID);
  auto phi = block->getData<field::GhostLayerField<double, 1>>(phiID);

  auto &kT = this->kT_;
  auto &f_ext0 = this->f_ext0_;
  auto &f_ext2 = this->f_ext2_;
  auto &z = this->z_;
  auto &D = this->D_;
  auto &f_ext1 = this->f_ext1_;
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(j->nrOfGhostLayers()));
  double *RESTRICT _data_j = j->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phi->nrOfGhostLayers()));
  double *RESTRICT const _data_phi = phi->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho->nrOfGhostLayers()));
  double *RESTRICT const _data_rho = rho->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(j->xSize()) + 2));
  const int64_t _size_j_0 = int64_t(cell_idx_c(j->xSize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(j->ySize()) + 2));
  const int64_t _size_j_1 = int64_t(cell_idx_c(j->ySize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(j->zSize()) + 2));
  const int64_t _size_j_2 = int64_t(cell_idx_c(j->zSize()) + 2);
  const int64_t _stride_j_0 = int64_t(j->xStride());
  const int64_t _stride_j_1 = int64_t(j->yStride());
  const int64_t _stride_j_2 = int64_t(j->zStride());
  const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
  const int64_t _stride_phi_0 = int64_t(phi->xStride());
  const int64_t _stride_phi_1 = int64_t(phi->yStride());
  const int64_t _stride_phi_2 = int64_t(phi->zStride());
  const int64_t _stride_rho_0 = int64_t(rho->xStride());
  const int64_t _stride_rho_1 = int64_t(rho->yStride());
  const int64_t _stride_rho_2 = int64_t(rho->zStride());
  internal_diffusivefluxkernelwithelectrostatic_diffusivefluxkernelwithelectrostatic::
      diffusivefluxkernelwithelectrostatic_diffusivefluxkernelwithelectrostatic(
          D, _data_j, _data_phi, _data_rho, _size_j_0, _size_j_1, _size_j_2,
          _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_phi_0,
          _stride_phi_1, _stride_phi_2, _stride_rho_0, _stride_rho_1,
          _stride_rho_2, f_ext0, f_ext1, f_ext2, kT, z);
}

void DiffusiveFluxKernelWithElectrostatic::runOnCellInterval(
    const shared_ptr<StructuredBlockStorage> &blocks,
    const CellInterval &globalCellInterval, cell_idx_t ghostLayers,
    IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto rho = block->getData<field::GhostLayerField<double, 1>>(rhoID);
  auto j = block->getData<field::GhostLayerField<double, 13>>(jID);
  auto phi = block->getData<field::GhostLayerField<double, 1>>(phiID);

  auto &kT = this->kT_;
  auto &f_ext0 = this->f_ext0_;
  auto &f_ext2 = this->f_ext2_;
  auto &z = this->z_;
  auto &D = this->D_;
  auto &f_ext1 = this->f_ext1_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(j->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(j->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(j->nrOfGhostLayers()));
  double *RESTRICT _data_j =
      j->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phi->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phi->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phi->nrOfGhostLayers()));
  double *RESTRICT const _data_phi =
      phi->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho->nrOfGhostLayers()));
  double *RESTRICT const _data_rho =
      rho->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.xSize()) + 2));
  const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.ySize()) + 2));
  const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.zSize()) + 2));
  const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 2);
  const int64_t _stride_j_0 = int64_t(j->xStride());
  const int64_t _stride_j_1 = int64_t(j->yStride());
  const int64_t _stride_j_2 = int64_t(j->zStride());
  const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
  const int64_t _stride_phi_0 = int64_t(phi->xStride());
  const int64_t _stride_phi_1 = int64_t(phi->yStride());
  const int64_t _stride_phi_2 = int64_t(phi->zStride());
  const int64_t _stride_rho_0 = int64_t(rho->xStride());
  const int64_t _stride_rho_1 = int64_t(rho->yStride());
  const int64_t _stride_rho_2 = int64_t(rho->zStride());
  internal_diffusivefluxkernelwithelectrostatic_diffusivefluxkernelwithelectrostatic::
      diffusivefluxkernelwithelectrostatic_diffusivefluxkernelwithelectrostatic(
          D, _data_j, _data_phi, _data_rho, _size_j_0, _size_j_1, _size_j_2,
          _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_phi_0,
          _stride_phi_1, _stride_phi_2, _stride_rho_0, _stride_rho_1,
          _stride_rho_2, f_ext0, f_ext1, f_ext2, kT, z);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif