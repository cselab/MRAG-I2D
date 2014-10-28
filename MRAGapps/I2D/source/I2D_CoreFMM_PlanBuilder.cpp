/*
 *  I2D_CoreFMM_PlanBuilder.cpp
 *  I2D_ROCKS
 *
 *  Created by Roman Schaerer on 12/26/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include <tbb/parallel_sort.h>

#include "I2D_CoreFMM_PlanBuilder.h"

bool PlanBuilder::_isclose_box(const tBox& box, const BlockInfo& info, tBox::Btype _theta) const
{
	const Real h_block = pow(0.5, info.level);
	
	const Real blockCenter[2] = {
		(0.5 + info.index[0])*h_block,
		(0.5 + info.index[1])*h_block
	};
	
	const Real b2b_dist = sqrt(pow(box.expansions.Center[0]-blockCenter[0], 2) +
							   pow(box.expansions.Center[1]-blockCenter[1], 2));
	
	const Real radius_source = sqrt(pow(box.h[0],2) +
									pow(box.h[1],2));
	
	const Real radius_target = h_block*sqrt(2.0);
	const Real denom = max(b2b_dist-radius_target, std::numeric_limits<Real>::epsilon());
	
	return (radius_source>=_theta*denom);
}

bool PlanBuilder::_is_intersecting(const tBox& box, const BlockInfo& info) const
{
	Real min_pos[2], max_pos[2];
	info.pos(min_pos, 0,0);
	info.pos(max_pos, B::sizeX-1, B::sizeY-1);
	
	const Real intersection[2] = {
		min(max_pos[0], (Real)(box.center[0] + box.h[0]*0.5)) - max(min_pos[0], (Real)(box.center[0] - box.h[0]*0.5)),
		min(max_pos[1], (Real)(box.center[1] + box.h[1]*0.5)) - max(min_pos[1], (Real)(box.center[1] - box.h[1]*0.5))
	};
	
	return intersection[0]>=0 && intersection[1]>=0;
}


bool PlanBuilder::_ws_box_barnes_hut_and (const tBox& source_box, const BlockInfo target_block, tBox::Btype theta) const 
{
	bool bResult = false;

	if (!_is_intersecting (source_box, target_block)) 
	{
		Real minp[2], maxp[2];
		
		target_block.pos(minp, 0,0);
		target_block.pos(maxp, B::sizeX-1,B::sizeX-1);
		
		Real xclosest[2] = {
			max(minp[0], min(maxp[0], source_box.expansions.Center[0])),
			max(minp[1], min(maxp[1], source_box.expansions.Center[1]))
		};
		
		bResult = ws_barnes_hut(&source_box, xclosest, theta);
	}

	assert(bResult==false || bResult == _ws_box_barnes_hut_and_DUMB(source_box, target_block, theta));

	return bResult;
}

bool PlanBuilder::_ws_box_barnes_hut_and_DUMB (const tBox& source_box, const BlockInfo& target_block, tBox::Btype theta) const 
{
	for(int iy=0; iy<B::sizeY; iy++) {
		for(int ix=0; ix<B::sizeX; ix++) {
			Real target_pos[2] = {0,0};
			target_block.pos(target_pos, ix, iy);
			
			if(!ws_barnes_hut(&source_box, target_pos, theta)) {
				return false;
			}
		}
	}
	
	if (_is_intersecting (source_box, target_block))
	{
		printf("Akamon now.\n");
		abort();
		return false;
	}
	
	return true;
}

void PlanBuilder::run() const {
	tbb::parallel_for (blocked_range<int> (0,m_num_target_blocks), *this, auto_partitioner ());
}

void PlanBuilder::operator () (const blocked_range<int>& range) const {	
	typedef vector< pair<int,int>, scalable_allocator< pair<int, int> > > Intervals;
	typedef vector< tBox*, scalable_allocator< tBox* > > Expansions;
	
	for (int iblock=range.begin (); iblock != range.end(); ++iblock)
	{
		const BlockInfo info = m_target_blocks [iblock];
		
		HCFMM::BoxIterator<tBox,tbb::scalable_allocator> it1(m_root_node);
		
		bool canRemove;
		
		while(it1!=NULL && (it1->nParticles>0))
		{
			canRemove=true;
			const tBox& current = *it1;
			const bool is_close = _isclose_box (current, info, _THETA);
			
			if (!is_close && current.parent!=NULL) 
			{
				const bool parent_is_close = _isclose_box (*current.parent, info, _THETA);
				
				if (parent_is_close)
					if (!(_ws_box_barnes_hut_and (*current.parent, info, _THETA))) 
						m_plan_ptr->addIndirectInteraction (&current, iblock);
			}
			else 
			{
				if (current.parent==NULL || !(_ws_box_barnes_hut_and (*current.parent, info, _THETA))) 
				{
					if (_ws_box_barnes_hut_and (current, info, _THETA))
						m_plan_ptr->addIndirectInteraction (&current, iblock);
					else 
					{
						if (!it1->isleaf) 
							canRemove=false;
						else
						{
							if (_is_intersecting(current, info))
								m_plan_ptr->addDirectInteraction (&current, iblock, true);
							else 
								m_plan_ptr->addDirectInteraction (&current, iblock, false);
						}
					}
				}
			}
			
			if(canRemove)
				it1.advanceRemove();
			else
				it1++;
		}
		
		m_plan_ptr->merge_direct_intervals(iblock);
	}
}
