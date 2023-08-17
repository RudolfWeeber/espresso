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
//! \\file ContinuityKernel_single_precision.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

// kernel generated with pystencils v1.2, lbmpy v1.3.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 065ce5f311850371a97ac4766f47dbb5ca8424ba

#include <cmath>

#include "ContinuityKernel_single_precision.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#define FUNC_PREFIX

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
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

namespace internal_990034b4e4dd57d2802b4bcb5f716e46 {
static FUNC_PREFIX void continuitykernel_single_precision_continuitykernel_single_precision(float *RESTRICT const _data_j, float *RESTRICT _data_rho, int64_t const _size_j_0, int64_t const _size_j_1, int64_t const _size_j_2, int64_t const _stride_j_0, int64_t const _stride_j_1, int64_t const _stride_j_2, int64_t const _stride_j_3, int64_t const _stride_rho_0, int64_t const _stride_rho_1, int64_t const _stride_rho_2) {
  for (int64_t ctr_2 = 1; ctr_2 < _size_j_2 - 1; ctr_2 += 1) {
    float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
    float *RESTRICT _data_j_20_30 = _data_j + _stride_j_2 * ctr_2;
    float *RESTRICT _data_j_20_31 = _data_j + _stride_j_2 * ctr_2 + _stride_j_3;
    float *RESTRICT _data_j_20_310 = _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    float *RESTRICT _data_j_20_311 = _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    float *RESTRICT _data_j_20_312 = _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    float *RESTRICT _data_j_20_32 = _data_j + _stride_j_2 * ctr_2 + 2 * _stride_j_3;
    float *RESTRICT _data_j_20_33 = _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
    float *RESTRICT _data_j_20_34 = _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
    float *RESTRICT _data_j_20_35 = _data_j + _stride_j_2 * ctr_2 + 5 * _stride_j_3;
    float *RESTRICT _data_j_20_36 = _data_j + _stride_j_2 * ctr_2 + 6 * _stride_j_3;
    float *RESTRICT _data_j_20_37 = _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
    float *RESTRICT _data_j_20_38 = _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
    float *RESTRICT _data_j_20_39 = _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_2m1_36 = _data_j + _stride_j_2 * ctr_2 - _stride_j_2 + 6 * _stride_j_3;
    float *RESTRICT _data_j_2m1_310 = _data_j + _stride_j_2 * ctr_2 - _stride_j_2 + 10 * _stride_j_3;
    float *RESTRICT _data_j_2m1_38 = _data_j + _stride_j_2 * ctr_2 - _stride_j_2 + 8 * _stride_j_3;
    float *RESTRICT _data_j_2m1_312 = _data_j + _stride_j_2 * ctr_2 - _stride_j_2 + 12 * _stride_j_3;
    float *RESTRICT _data_j_21_35 = _data_j + _stride_j_2 * ctr_2 + _stride_j_2 + 5 * _stride_j_3;
    float *RESTRICT _data_j_21_39 = _data_j + _stride_j_2 * ctr_2 + _stride_j_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_21_37 = _data_j + _stride_j_2 * ctr_2 + _stride_j_2 + 7 * _stride_j_3;
    float *RESTRICT _data_j_21_311 = _data_j + _stride_j_2 * ctr_2 + _stride_j_2 + 11 * _stride_j_3;
    float *RESTRICT _data_j_21_32 = _data_j + _stride_j_2 * ctr_2 + _stride_j_2 + 2 * _stride_j_3;
    for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
      float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
      float *RESTRICT _data_j_20_30_10 = _stride_j_1 * ctr_1 + _data_j_20_30;
      float *RESTRICT _data_j_20_31_10 = _stride_j_1 * ctr_1 + _data_j_20_31;
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      float *RESTRICT _data_j_20_32_10 = _stride_j_1 * ctr_1 + _data_j_20_32;
      float *RESTRICT _data_j_20_33_10 = _stride_j_1 * ctr_1 + _data_j_20_33;
      float *RESTRICT _data_j_20_34_10 = _stride_j_1 * ctr_1 + _data_j_20_34;
      float *RESTRICT _data_j_20_35_10 = _stride_j_1 * ctr_1 + _data_j_20_35;
      float *RESTRICT _data_j_20_36_10 = _stride_j_1 * ctr_1 + _data_j_20_36;
      float *RESTRICT _data_j_20_37_10 = _stride_j_1 * ctr_1 + _data_j_20_37;
      float *RESTRICT _data_j_20_38_10 = _stride_j_1 * ctr_1 + _data_j_20_38;
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_2m1_36_10 = _stride_j_1 * ctr_1 + _data_j_2m1_36;
      float *RESTRICT _data_j_2m1_310_11 = _stride_j_1 * ctr_1 + _stride_j_1 + _data_j_2m1_310;
      float *RESTRICT _data_j_2m1_38_11 = _stride_j_1 * ctr_1 + _stride_j_1 + _data_j_2m1_38;
      float *RESTRICT _data_j_2m1_312_1m1 = _stride_j_1 * ctr_1 - _stride_j_1 + _data_j_2m1_312;
      float *RESTRICT _data_j_20_33_11 = _stride_j_1 * ctr_1 + _stride_j_1 + _data_j_20_33;
      float *RESTRICT _data_j_20_31_11 = _stride_j_1 * ctr_1 + _stride_j_1 + _data_j_20_31;
      float *RESTRICT _data_j_20_34_1m1 = _stride_j_1 * ctr_1 - _stride_j_1 + _data_j_20_34;
      float *RESTRICT _data_j_21_35_10 = _stride_j_1 * ctr_1 + _data_j_21_35;
      float *RESTRICT _data_j_21_39_11 = _stride_j_1 * ctr_1 + _stride_j_1 + _data_j_21_39;
      float *RESTRICT _data_j_21_37_11 = _stride_j_1 * ctr_1 + _stride_j_1 + _data_j_21_37;
      float *RESTRICT _data_j_21_311_1m1 = _stride_j_1 * ctr_1 - _stride_j_1 + _data_j_21_311;
      float *RESTRICT _data_j_21_32_10 = _stride_j_1 * ctr_1 + _data_j_21_32;
      for (int64_t ctr_0 = 1; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
        _data_rho_20_10[_stride_rho_0 * ctr_0] = -1.0f * _data_j_20_30_10[_stride_j_0 * ctr_0] - 1.0f * _data_j_20_310_10[_stride_j_0 * ctr_0] - 1.0f * _data_j_20_311_10[_stride_j_0 * ctr_0] - 1.0f * _data_j_20_312_10[_stride_j_0 * ctr_0] - 1.0f * _data_j_20_31_10[_stride_j_0 * ctr_0] - 1.0f * _data_j_20_32_10[_stride_j_0 * ctr_0] - 1.0f * _data_j_20_33_10[_stride_j_0 * ctr_0] - 1.0f * _data_j_20_34_10[_stride_j_0 * ctr_0] - 1.0f * _data_j_20_35_10[_stride_j_0 * ctr_0] - 1.0f * _data_j_20_36_10[_stride_j_0 * ctr_0] - 1.0f * _data_j_20_37_10[_stride_j_0 * ctr_0] - 1.0f * _data_j_20_38_10[_stride_j_0 * ctr_0] - 1.0f * _data_j_20_39_10[_stride_j_0 * ctr_0] + _data_j_20_30_10[_stride_j_0 * ctr_0 + _stride_j_0] + _data_j_20_31_11[_stride_j_0 * ctr_0] + _data_j_20_33_11[_stride_j_0 * ctr_0 + _stride_j_0] + _data_j_20_34_1m1[_stride_j_0 * ctr_0 + _stride_j_0] + _data_j_21_311_1m1[_stride_j_0 * ctr_0 + _stride_j_0] + _data_j_21_32_10[_stride_j_0 * ctr_0] + _data_j_21_35_10[_stride_j_0 * ctr_0 + _stride_j_0] + _data_j_21_37_11[_stride_j_0 * ctr_0] + _data_j_21_39_11[_stride_j_0 * ctr_0 + _stride_j_0] + _data_j_2m1_310_11[_stride_j_0 * ctr_0 + _stride_j_0] + _data_j_2m1_312_1m1[_stride_j_0 * ctr_0 + _stride_j_0] + _data_j_2m1_36_10[_stride_j_0 * ctr_0 + _stride_j_0] + _data_j_2m1_38_11[_stride_j_0 * ctr_0] + _data_rho_20_10[_stride_rho_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_990034b4e4dd57d2802b4bcb5f716e46

void ContinuityKernel_single_precision::run(IBlock *block) {
  auto rho = block->getData<field::GhostLayerField<float, 1>>(rhoID);
  auto j = block->getData<field::GhostLayerField<float, 13>>(jID);

  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(j->nrOfGhostLayers()));
  float *RESTRICT const _data_j = j->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho->nrOfGhostLayers()));
  float *RESTRICT _data_rho = rho->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(), int64_t(cell_idx_c(j->xSize()) + 2));
  const int64_t _size_j_0 = int64_t(cell_idx_c(j->xSize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(), int64_t(cell_idx_c(j->ySize()) + 2));
  const int64_t _size_j_1 = int64_t(cell_idx_c(j->ySize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(), int64_t(cell_idx_c(j->zSize()) + 2));
  const int64_t _size_j_2 = int64_t(cell_idx_c(j->zSize()) + 2);
  const int64_t _stride_j_0 = int64_t(j->xStride());
  const int64_t _stride_j_1 = int64_t(j->yStride());
  const int64_t _stride_j_2 = int64_t(j->zStride());
  const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
  const int64_t _stride_rho_0 = int64_t(rho->xStride());
  const int64_t _stride_rho_1 = int64_t(rho->yStride());
  const int64_t _stride_rho_2 = int64_t(rho->zStride());
  internal_990034b4e4dd57d2802b4bcb5f716e46::continuitykernel_single_precision_continuitykernel_single_precision(_data_j, _data_rho, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_rho_0, _stride_rho_1, _stride_rho_2);
}

void ContinuityKernel_single_precision::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto rho = block->getData<field::GhostLayerField<float, 1>>(rhoID);
  auto j = block->getData<field::GhostLayerField<float, 13>>(jID);

  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(j->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(j->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(j->nrOfGhostLayers()));
  float *RESTRICT const _data_j = j->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho->nrOfGhostLayers()));
  float *RESTRICT _data_rho = rho->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 2));
  const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 2));
  const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 2));
  const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 2);
  const int64_t _stride_j_0 = int64_t(j->xStride());
  const int64_t _stride_j_1 = int64_t(j->yStride());
  const int64_t _stride_j_2 = int64_t(j->zStride());
  const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
  const int64_t _stride_rho_0 = int64_t(rho->xStride());
  const int64_t _stride_rho_1 = int64_t(rho->yStride());
  const int64_t _stride_rho_2 = int64_t(rho->zStride());
  internal_990034b4e4dd57d2802b4bcb5f716e46::continuitykernel_single_precision_continuitykernel_single_precision(_data_j, _data_rho, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_rho_0, _stride_rho_1, _stride_rho_2);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif