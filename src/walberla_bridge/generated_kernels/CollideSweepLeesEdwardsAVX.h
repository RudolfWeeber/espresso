// kernel generated with pystencils v0.3.3+44.gf2d4181, lbmpy
// v0.3.1+53.g2faceda, lbmpy_walberla/pystencils_walberla from commit
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
//! \\file CollideSweepLeesEdwardsAVX.h
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

class CollideSweepLeesEdwardsAVX {
public:
  CollideSweepLeesEdwardsAVX(BlockDataID densityID_, BlockDataID forceID_,
                             BlockDataID pdfsID_, BlockDataID velocityID_,
                             double omega_bulk, double omega_even,
                             double omega_odd, double omega_shear,
                             bool points_down, bool points_up)
      : densityID(densityID_), forceID(forceID_), pdfsID(pdfsID_),
        velocityID(velocityID_), omega_bulk_(omega_bulk),
        omega_even_(omega_even), omega_odd_(omega_odd),
        omega_shear_(omega_shear), points_down_(points_down),
        points_up_(points_up){};

  void operator()(IBlock *block);
  void runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers, IBlock *block);

  static std::function<void(IBlock *)>
  getSweep(const shared_ptr<CollideSweepLeesEdwardsAVX> &kernel) {
    return [kernel](IBlock *b) { (*kernel)(b); };
  }

  static std::function<void(IBlock *)>
  getSweepOnCellInterval(const shared_ptr<CollideSweepLeesEdwardsAVX> &kernel,
                         const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers = 1) {
    return [kernel, blocks, globalCellInterval, ghostLayers](IBlock *b) {
      kernel->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b);
    };
  }

  BlockDataID densityID;
  BlockDataID forceID;
  BlockDataID pdfsID;
  BlockDataID velocityID;
  double omega_bulk_;
  double omega_even_;
  double omega_odd_;
  double omega_shear_;
  bool points_down_;
  bool points_up_;
};

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif