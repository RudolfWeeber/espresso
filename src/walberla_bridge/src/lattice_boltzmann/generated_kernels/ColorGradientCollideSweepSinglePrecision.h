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
//! \\file ColorGradientCollideSweepSinglePrecision.h
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 17fc54c872bd8ceabf271a7e9e636c7c583f55af


#pragma once
#include "core/DataTypes.h"
#include "core/logging/Logging.h"

#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include <functional>
#include <unordered_map>



#ifdef __GNUC__
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wunused-parameter"
#   pragma GCC diagnostic ignored "-Wreorder"
#endif

namespace walberla {
namespace pystencils {


class ColorGradientCollideSweepSinglePrecision
{
public:
   ColorGradientCollideSweepSinglePrecision( BlockDataID force_aID_, BlockDataID force_bID_, BlockDataID pdfs_aID_, BlockDataID pdfs_bID_, BlockDataID phasefieldID_, BlockDataID rho_aID_, BlockDataID rho_bID_, BlockDataID velocityID_, float beta, float omega_even_a, float omega_even_b, float omega_odd_a, float omega_odd_b, float omega_shear_a, float omega_shear_b, float sigma )
     : force_aID(force_aID_), force_bID(force_bID_), pdfs_aID(pdfs_aID_), pdfs_bID(pdfs_bID_), phasefieldID(phasefieldID_), rho_aID(rho_aID_), rho_bID(rho_bID_), velocityID(velocityID_), beta_(beta), omega_even_a_(omega_even_a), omega_even_b_(omega_even_b), omega_odd_a_(omega_odd_a), omega_odd_b_(omega_odd_b), omega_shear_a_(omega_shear_a), omega_shear_b_(omega_shear_b), sigma_(sigma)
   {}

   

   void run(IBlock * block);

   void runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block);

   
   void operator() (IBlock * block)
   {
     run(block);
   }
   

   static std::function<void (IBlock *)> getSweep(const shared_ptr<ColorGradientCollideSweepSinglePrecision> & kernel)
   {
     return [kernel]
            (IBlock * b)
            { kernel->run(b); };
   }

   static std::function<void (IBlock*)> getSweepOnCellInterval(const shared_ptr<ColorGradientCollideSweepSinglePrecision> & kernel, const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers=1)
   {
     return [kernel, blocks, globalCellInterval, ghostLayers]
            (IBlock * b)
            { kernel->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b); };
   }

   std::function<void (IBlock *)> getSweep()
   {
     return [this]
            (IBlock * b)
            { this->run(b); };
   }

   std::function<void (IBlock *)> getSweepOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers=1)
   {
     return [this, blocks, globalCellInterval, ghostLayers]
            (IBlock * b)
            { this->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b); };
   }

   
   void configure( const shared_ptr<StructuredBlockStorage> & /*blocks*/, IBlock * /*block*/ ){}
   

   

   inline float getBeta() const { return beta_; }
   inline float getOmega_even_a() const { return omega_even_a_; }
   inline float getOmega_even_b() const { return omega_even_b_; }
   inline float getOmega_odd_a() const { return omega_odd_a_; }
   inline float getOmega_odd_b() const { return omega_odd_b_; }
   inline float getOmega_shear_a() const { return omega_shear_a_; }
   inline float getOmega_shear_b() const { return omega_shear_b_; }
   inline float getSigma() const { return sigma_; }
   inline void setBeta(const float value) { beta_ = value; }
   inline void setOmega_even_a(const float value) { omega_even_a_ = value; }
   inline void setOmega_even_b(const float value) { omega_even_b_ = value; }
   inline void setOmega_odd_a(const float value) { omega_odd_a_ = value; }
   inline void setOmega_odd_b(const float value) { omega_odd_b_ = value; }
   inline void setOmega_shear_a(const float value) { omega_shear_a_ = value; }
   inline void setOmega_shear_b(const float value) { omega_shear_b_ = value; }
   inline void setSigma(const float value) { sigma_ = value; }

private:
   
   BlockDataID force_aID;
   BlockDataID force_bID;
   BlockDataID pdfs_aID;
   BlockDataID pdfs_bID;
   BlockDataID phasefieldID;
   BlockDataID rho_aID;
   BlockDataID rho_bID;
   BlockDataID velocityID;
   float beta_;
   float omega_even_a_;
   float omega_even_b_;
   float omega_odd_a_;
   float omega_odd_b_;
   float omega_shear_a_;
   float omega_shear_b_;
   float sigma_;

   

   
};


} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif
