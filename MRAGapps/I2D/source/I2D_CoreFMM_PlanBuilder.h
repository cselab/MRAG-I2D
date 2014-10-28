/*
 *  I2D_CoreFMM_PlanBuilder.h
 *  I2D_ROCKS
 *
 *  Created by Roman Schaerer on 12/26/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once

extern  double _THETA;
#define _FMMSILENT

#include "I2D_CoreFMM_Plan.h"

class PlanBuilder {
public:
	typedef HCFMM::Box<_VortexExpansions<VelocitySourceParticle, _ORDER_>,_FMM_MAX_LEVEL_> tBox;

	PlanBuilder (Plan* _plan, tBox* const _root_node,
			BlockInfo* _target_blocks, int _num_target_blocks) :
				m_plan_ptr (_plan), m_root_node (_root_node),
				m_target_blocks (_target_blocks), m_num_target_blocks (_num_target_blocks),
				m_root_source_particle (&m_root_node->vparticles[0]){
		assert (m_plan_ptr != NULL && m_root_node != NULL && m_target_blocks != NULL);
		assert (m_num_target_blocks >= 0);
	}

	void run () const;
	void operator () (const blocked_range<int>&) const;

	Plan* const m_plan_ptr;
	BlockInfo* const m_target_blocks;
	tBox * const m_root_node;
	VelocitySourceParticle* const m_root_source_particle;
	const int m_num_target_blocks;

protected:

	bool _isclose_box(const tBox& sourceBox, const BlockInfo& info, tBox::Btype _theta) const;
	bool _is_intersecting(const tBox& box, const BlockInfo& info) const;
	bool _ws_box_barnes_hut_and (const tBox& source_box, const BlockInfo target_block, tBox::Btype theta) const;
	bool _ws_box_barnes_hut_and_DUMB (const tBox& source_box, const BlockInfo& target_block, tBox::Btype theta) const;

};
