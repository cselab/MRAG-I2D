//
//  I2D_CoreFMM_PlanBuilderWim.cpp
//  I2D_ROCKS
//
//  Created by Wim van Rees on 28/10/13.
//
//

#include <tbb/parallel_sort.h>

#include "I2D_CoreFMM_PlanBuilderWim.h"

bool PlanBuilderWim::isCloseBox(Real srcBoxCtr[2],Real srcBoxRad, const Real trgBoxCtr[2],const Real halfTrgBoxWidth) const
{
    
    // minimum distance given that target box is a cube
    const Real x0 = trgBoxCtr[0]-halfTrgBoxWidth;
    const Real x1 = trgBoxCtr[0]+halfTrgBoxWidth;
    const Real y0 = trgBoxCtr[1]-halfTrgBoxWidth;
    const Real y1 = trgBoxCtr[1]+halfTrgBoxWidth;
    const Real x = srcBoxCtr[0];
    const Real y = srcBoxCtr[1];
    
    const Real dx = x < x0 ? (x0-x) : (x > x1 ? (x-x1) : 0); // can replace this with min/max
    const Real dy = y < y0 ? (y0-y) : (y > y1 ? (y-y1) : 0);
    const Real min_p2bDist = std::sqrt(dx*dx+dy*dy);
    
    const Real denom = std::max(min_p2bDist, std::numeric_limits<Real>::epsilon());
    return (srcBoxRad>=_THETA*denom);
}

bool PlanBuilderWim::isIntersectingBox(Real srcBoxCtr[2],Real srcBoxRad, const Real trgBoxCtr[2],const Real halfTrgBoxWidth) const
{
    const Real intersection[3] = {
        std::min(trgBoxCtr[0]+halfTrgBoxWidth, srcBoxCtr[0]+srcBoxRad) - std::max(trgBoxCtr[0]-halfTrgBoxWidth,srcBoxCtr[0]-srcBoxRad),
        std::min(trgBoxCtr[1]+halfTrgBoxWidth, srcBoxCtr[1]+srcBoxRad) - std::max(trgBoxCtr[1]-halfTrgBoxWidth,srcBoxCtr[1]-srcBoxRad),
    };
    
    return intersection[0]>=0 && intersection[1]>=0;
}

void PlanBuilderWim::run() const {
	tbb::parallel_for (blocked_range<int> (0,m_num_target_blocks), *this, auto_partitioner ());
}

void PlanBuilderWim::operator () (const blocked_range<int>& range) const
{
	typedef vector< pair<int,int>, scalable_allocator< pair<int, int> > > Intervals;
	typedef vector< tBox*, scalable_allocator< tBox* > > Expansions;
	
	for (int iblock=range.begin (); iblock != range.end(); ++iblock)
	{
		const BlockInfo info = m_target_blocks [iblock];
		
		HCFMM::BoxIterator<tBox,tbb::scalable_allocator> srcBox(m_root_node);
		
		Real block_org[2];
		info.pos(block_org,0,0);        
        const Real halfTargetWidth = 0.5*(_BLOCKSIZE_-1)*info.h[0];//our destination blocks have size _BLOCKSIZE_
        const Real targetCenter[2] = {block_org[0] + halfTargetWidth,block_org[1] + halfTargetWidth};
        
		while(srcBox!=NULL && (srcBox->nParticles>0))
		{
			bool canRemove=true;
			
            const bool isClose = isCloseBox(srcBox->expansions.Center, srcBox->expansions.Radius,targetCenter,halfTargetWidth);
            const bool isIntersecting = isIntersectingBox(srcBox->expansions.Center, srcBox->expansions.Radius,targetCenter,halfTargetWidth);
            
            const tBox& current = *srcBox;
            
            if(isClose || isIntersecting)
            {
                if(not srcBox->isleaf)
                    canRemove=false;
                else
                    m_plan_ptr->addDirectInteraction (&current, iblock, isIntersecting);
            }
            else
				m_plan_ptr->addIndirectInteraction (&current, iblock);

			if(canRemove)
				srcBox.advanceRemove();
			else
				srcBox++;
		}
		
		m_plan_ptr->merge_direct_intervals(iblock);
	}
}
