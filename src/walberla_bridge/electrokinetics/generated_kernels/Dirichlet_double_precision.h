// kernel generated with pystencils v0.4.4, lbmpy v0.4.4, lbmpy_walberla/pystencils_walberla from commit 08f04ef64f95609b47838db85862033a1600afa1

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
//! \\file Dirichlet_double_precision.h
//! \\author pystencils
//======================================================================================================================

#pragma once
#include "core/DataTypes.h"

#include "field/GhostLayerField.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "blockforest/StructuredBlockForest.h"
#include "field/FlagField.h"
#include "core/debug/Debug.h"

#include <set>
#include <vector>



#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

namespace walberla {
namespace pystencils {


class Dirichlet_double_precision
{
public:
    struct IndexInfo { 
        int32_t x;
        int32_t y;
        int32_t z;
        int32_t dir;
        double value;
        IndexInfo(int32_t x_, int32_t y_, int32_t z_, int32_t dir_) : x(x_), y(y_), z(z_), dir(dir_), value() {}
        bool operator==(const IndexInfo & o) const {
            return x == o.x && y == o.y && z == o.z && dir == o.dir && floatIsEqual(value, o.value);
        }
    };



    class IndexVectors
    {
    public:
        using CpuIndexVector = std::vector<IndexInfo>;

        enum Type {
            ALL = 0,
            INNER = 1,
            OUTER = 2,
            NUM_TYPES = 3
        };

        IndexVectors() = default;
        bool operator==(IndexVectors & other) { return other.cpuVectors_ == cpuVectors_; }

        CpuIndexVector & indexVector(Type t) { return cpuVectors_[t]; }
        IndexInfo * pointerCpu(Type t)  { return &(cpuVectors_[t][0]); }

        void syncGPU()
        {
            
        }

    private:
        std::vector<CpuIndexVector> cpuVectors_{NUM_TYPES};

        
    };

    Dirichlet_double_precision( const shared_ptr<StructuredBlockForest> & blocks,
                   BlockDataID fieldID_, std::function<double(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)>& dirichletCallback )
        :elementInitaliser(dirichletCallback), fieldID(fieldID_)
    {
        auto createIdxVector = []( IBlock * const , StructuredBlockStorage * const ) { return new IndexVectors(); };
        indexVectorID = blocks->addStructuredBlockData< IndexVectors >( createIdxVector, "IndexField_Dirichlet_double_precision");
    };

    void run (IBlock * block);

    void operator() (IBlock * block)
    {
        run(block);
    }

    void inner (IBlock * block);

    void outer (IBlock * block);

    std::function<void (IBlock *)> getSweep()
    {
        return [this]
               (IBlock * b)
               { this->run(b); };
    }

    std::function<void (IBlock *)> getInnerSweep()
    {
        return [this]
               (IBlock * b)
               { this->inner(b); };
    }

    std::function<void (IBlock *)> getOuterSweep()
    {
        return [this]
               (IBlock * b)
               { this->outer(b); };
    }

    template<typename FlagField_T>
    void fillFromFlagField( const shared_ptr<StructuredBlockForest> & blocks, ConstBlockDataID flagFieldID,
                            FlagUID boundaryFlagUID, FlagUID domainFlagUID)
    {
        for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
            fillFromFlagField<FlagField_T>(blocks, &*blockIt, flagFieldID, boundaryFlagUID, domainFlagUID );
    }


    template<typename FlagField_T>
    void fillFromFlagField( const shared_ptr<StructuredBlockForest> &blocks, IBlock * block, ConstBlockDataID flagFieldID,
                            FlagUID boundaryFlagUID, FlagUID domainFlagUID )
    {
        auto * indexVectors = block->getData< IndexVectors > ( indexVectorID );
        auto & indexVectorAll = indexVectors->indexVector(IndexVectors::ALL);
        auto & indexVectorInner = indexVectors->indexVector(IndexVectors::INNER);
        auto & indexVectorOuter = indexVectors->indexVector(IndexVectors::OUTER);

        auto * flagField = block->getData< FlagField_T > ( flagFieldID );
        

        if( !(flagField->flagExists(boundaryFlagUID) && flagField->flagExists(domainFlagUID) ))
            return;

        auto boundaryFlag = flagField->getFlag(boundaryFlagUID);
        auto domainFlag = flagField->getFlag(domainFlagUID);

        auto inner = flagField->xyzSize();
        inner.expand( cell_idx_t(-1) );

        indexVectorAll.clear();
        indexVectorInner.clear();
        indexVectorOuter.clear();

        
        auto flagWithGLayers = flagField->xyzSizeWithGhostLayer();
        double dot = 0.0; double maxn = 0.0;
        cell_idx_t calculated_idx = 0;
        cell_idx_t dx = 0; cell_idx_t dy = 0;   cell_idx_t dz = 0; 
        cell_idx_t sum_x = 0; cell_idx_t sum_y = 0;  cell_idx_t sum_z = 0; 
        for( auto it = flagField->beginWithGhostLayerXYZ(); it != flagField->end(); ++it )
        {
            
            if( ! isFlagSet(it, boundaryFlag) )
                continue;
            if ( flagWithGLayers.contains(it.x() + cell_idx_c(0), it.y() + cell_idx_c(0), it.z() + cell_idx_c(0)) && isFlagSet( it.neighbor(0, 0, 0 , 0 ), domainFlag ) )
            {
                
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  0 );
                double InitialisatonAdditionalData = elementInitaliser(Cell(it.x(), it.y(), it.z()), blocks, *block);
                element.value = InitialisatonAdditionalData;
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
                
            }
            

        }
        

        indexVectors->syncGPU();
    }

private:
    void run_impl(IBlock * block, IndexVectors::Type type);

    BlockDataID indexVectorID;
    std::function<double(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)> elementInitaliser; 
public:
    BlockDataID fieldID;
};



} // namespace pystencils
} // namespace walberla