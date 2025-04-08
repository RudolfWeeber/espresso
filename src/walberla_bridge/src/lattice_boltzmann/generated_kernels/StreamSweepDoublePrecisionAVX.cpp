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
//! \\file StreamSweepDoublePrecisionAVX.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7, lbmpy v1.3.7, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit f36fa0a68bae59f0b516f6587ea8fa7c24a41141

#include <cmath>

#include "StreamSweepDoublePrecisionAVX.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#include <immintrin.h>

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

namespace internal_91e2c9bdb4c4fa8a405803890749bf98 {
static FUNC_PREFIX void streamsweepdoubleprecisionavx_streamsweepdoubleprecisionavx(double *RESTRICT const _data_pdfs, double *RESTRICT _data_pdfs_tmp, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3) {
  for (int64_t ctr_2 = 1; ctr_2 < _size_pdfs_2 - 1; ctr_2 += 1) {
    for (int64_t ctr_1 = 1; ctr_1 < _size_pdfs_1 - 1; ctr_1 += 1) {
      {
        for (int64_t ctr_0 = 1; ctr_0 < (int64_t)((_size_pdfs_0 - 2) / (4)) * (4) + 1; ctr_0 += 4) {
          const __m256d streamed_0 = _mm256_load_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0]);
          const __m256d streamed_1 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0]);
          const __m256d streamed_2 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]);
          const __m256d streamed_3 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1]);
          const __m256d streamed_4 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1]);
          const __m256d streamed_5 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0]);
          const __m256d streamed_6 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]);
          const __m256d streamed_7 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1]);
          const __m256d streamed_8 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1]);
          const __m256d streamed_9 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1]);
          const __m256d streamed_10 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1]);
          const __m256d streamed_11 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0]);
          const __m256d streamed_12 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]);
          const __m256d streamed_13 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]);
          const __m256d streamed_14 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1]);
          const __m256d streamed_15 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]);
          const __m256d streamed_16 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]);
          const __m256d streamed_17 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]);
          const __m256d streamed_18 = _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1]);
          _mm256_store_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + ctr_0], streamed_0);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + _stride_pdfs_tmp_3 + ctr_0], streamed_1);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 2 * _stride_pdfs_tmp_3 + ctr_0], streamed_2);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 3 * _stride_pdfs_tmp_3 + ctr_0], streamed_3);
          _mm256_store_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 4 * _stride_pdfs_tmp_3 + ctr_0], streamed_4);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 5 * _stride_pdfs_tmp_3 + ctr_0], streamed_5);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 6 * _stride_pdfs_tmp_3 + ctr_0], streamed_6);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 7 * _stride_pdfs_tmp_3 + ctr_0], streamed_7);
          _mm256_store_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 8 * _stride_pdfs_tmp_3 + ctr_0], streamed_8);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 9 * _stride_pdfs_tmp_3 + ctr_0], streamed_9);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 10 * _stride_pdfs_tmp_3 + ctr_0], streamed_10);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 11 * _stride_pdfs_tmp_3 + ctr_0], streamed_11);
          _mm256_store_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 12 * _stride_pdfs_tmp_3 + ctr_0], streamed_12);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 13 * _stride_pdfs_tmp_3 + ctr_0], streamed_13);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 14 * _stride_pdfs_tmp_3 + ctr_0], streamed_14);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 15 * _stride_pdfs_tmp_3 + ctr_0], streamed_15);
          _mm256_store_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 16 * _stride_pdfs_tmp_3 + ctr_0], streamed_16);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 17 * _stride_pdfs_tmp_3 + ctr_0], streamed_17);
          _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 18 * _stride_pdfs_tmp_3 + ctr_0], streamed_18);
        }
        for (int64_t ctr_0 = (int64_t)((_size_pdfs_0 - 2) / (4)) * (4) + 1; ctr_0 < _size_pdfs_0 - 1; ctr_0 += 1) {
          const double streamed_0 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0];
          const double streamed_1 = _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0];
          const double streamed_2 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0];
          const double streamed_3 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1];
          const double streamed_4 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1];
          const double streamed_5 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0];
          const double streamed_6 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0];
          const double streamed_7 = _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1];
          const double streamed_8 = _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1];
          const double streamed_9 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1];
          const double streamed_10 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1];
          const double streamed_11 = _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0];
          const double streamed_12 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0];
          const double streamed_13 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1];
          const double streamed_14 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1];
          const double streamed_15 = _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0];
          const double streamed_16 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0];
          const double streamed_17 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1];
          const double streamed_18 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1];
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + ctr_0] = streamed_0;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + _stride_pdfs_tmp_3 + ctr_0] = streamed_1;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 2 * _stride_pdfs_tmp_3 + ctr_0] = streamed_2;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 3 * _stride_pdfs_tmp_3 + ctr_0] = streamed_3;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 4 * _stride_pdfs_tmp_3 + ctr_0] = streamed_4;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 5 * _stride_pdfs_tmp_3 + ctr_0] = streamed_5;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 6 * _stride_pdfs_tmp_3 + ctr_0] = streamed_6;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 7 * _stride_pdfs_tmp_3 + ctr_0] = streamed_7;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 8 * _stride_pdfs_tmp_3 + ctr_0] = streamed_8;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 9 * _stride_pdfs_tmp_3 + ctr_0] = streamed_9;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 10 * _stride_pdfs_tmp_3 + ctr_0] = streamed_10;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 11 * _stride_pdfs_tmp_3 + ctr_0] = streamed_11;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 12 * _stride_pdfs_tmp_3 + ctr_0] = streamed_12;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 13 * _stride_pdfs_tmp_3 + ctr_0] = streamed_13;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 14 * _stride_pdfs_tmp_3 + ctr_0] = streamed_14;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 15 * _stride_pdfs_tmp_3 + ctr_0] = streamed_15;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 16 * _stride_pdfs_tmp_3 + ctr_0] = streamed_16;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 17 * _stride_pdfs_tmp_3 + ctr_0] = streamed_17;
          _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 18 * _stride_pdfs_tmp_3 + ctr_0] = streamed_18;
        }
      }
    }
  }
}
} // namespace internal_91e2c9bdb4c4fa8a405803890749bf98

void StreamSweepDoublePrecisionAVX::run(IBlock *block) {

  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);
  field::GhostLayerField<double, 19> *pdfs_tmp;
  {
    if (cache_pdfs_.find(block) == cache_pdfs_.end()) {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_[block] = pdfs_tmp;
    } else {
      pdfs_tmp = cache_pdfs_[block];
    }
  }

  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT const _data_pdfs = pdfs->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs_tmp->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs_tmp->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(pdfs->xSize()) + 2))
  const int64_t _size_pdfs_0 = int64_t(int64_c(pdfs->xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(pdfs->ySize()) + 2))
  const int64_t _size_pdfs_1 = int64_t(int64_c(pdfs->ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(pdfs->zSize()) + 2))
  const int64_t _size_pdfs_2 = int64_t(int64_c(pdfs->zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
  const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
  const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
  internal_91e2c9bdb4c4fa8a405803890749bf98::streamsweepdoubleprecisionavx_streamsweepdoubleprecisionavx(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3);
  pdfs->swapDataPointers(pdfs_tmp);
}

void StreamSweepDoublePrecisionAVX::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {

  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);
  field::GhostLayerField<double, 19> *pdfs_tmp;
  {
    if (cache_pdfs_.find(block) == cache_pdfs_.end()) {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_[block] = pdfs_tmp;
    } else {
      pdfs_tmp = cache_pdfs_[block];
    }
  }

  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs_tmp->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs_tmp->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
  const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
  const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
  const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
  const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
  const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
  internal_91e2c9bdb4c4fa8a405803890749bf98::streamsweepdoubleprecisionavx_streamsweepdoubleprecisionavx(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3);
  pdfs->swapDataPointers(pdfs_tmp);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif
