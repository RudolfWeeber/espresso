// kernel generated with pystencils v0.3.4+4.g4fecf0c, lbmpy v0.3.4+6.g2faceda,
// lbmpy_walberla/pystencils_walberla from commit
// b17ca5caf00db7d19f86c5f85c6f67fec6c16aff

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
//! \\file StreamSweepDoublePrecisionAVX.h
//! \\author pystencils
//======================================================================================================================

#pragma once
#include "core/DataTypes.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"
#include <set>

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wreorder"
#endif

namespace walberla {
namespace pystencils {

class StreamSweepDoublePrecisionAVX {
public:
  StreamSweepDoublePrecisionAVX(BlockDataID forceID_, BlockDataID pdfsID_,
                                BlockDataID velocityID_)
      : forceID(forceID_), pdfsID(pdfsID_), velocityID(velocityID_){};

  ~StreamSweepDoublePrecisionAVX() {
    for (auto p : cache_pdfs_) {
      delete p;
    }
  }

  void operator()(IBlock *block);
  void runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers, IBlock *block);

  static std::function<void(IBlock *)>
  getSweep(const shared_ptr<StreamSweepDoublePrecisionAVX> &kernel) {
    return [kernel](IBlock *b) { (*kernel)(b); };
  }

  static std::function<void(IBlock *)> getSweepOnCellInterval(
      const shared_ptr<StreamSweepDoublePrecisionAVX> &kernel,
      const shared_ptr<StructuredBlockStorage> &blocks,
      const CellInterval &globalCellInterval, cell_idx_t ghostLayers = 1) {
    return [kernel, blocks, globalCellInterval, ghostLayers](IBlock *b) {
      kernel->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b);
    };
  }

  BlockDataID forceID;
  BlockDataID pdfsID;
  BlockDataID velocityID;

private:
  std::set<field::GhostLayerField<double, 19> *,
           field::SwapableCompare<field::GhostLayerField<double, 19> *>>
      cache_pdfs_;
};

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif