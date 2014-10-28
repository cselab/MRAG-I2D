//
//  I2D_CoreFMM_PlanBuilderWim.h
//  I2D_ROCKS
//
//  Created by Wim van Rees on 28/10/13.
//
//

#ifndef __I2D_ROCKS__I2D_CoreFMM_PlanBuilderWim__
#define __I2D_ROCKS__I2D_CoreFMM_PlanBuilderWim__

extern  double _THETA;
#define _FMMSILENT

#include "I2D_CoreFMM_Plan.h"

class PlanBuilderWim {
public:
	typedef HCFMM::Box<_VortexExpansions<VelocitySourceParticle, _ORDER_>,_FMM_MAX_LEVEL_> tBox;
    
	PlanBuilderWim (Plan* _plan, tBox* const _root_node,
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
    
    bool isCloseBox(Real srcBoxCtr[2],Real srcBoxRad, const Real trgBoxCtr[2],const Real halfTrgBoxWidth) const;
    bool isIntersectingBox(Real srcBoxCtr[2],Real srcBoxRad, const Real trgBoxCtr[2],const Real halfTrgBoxWidth) const;
};
#endif /* defined(__I2D_ROCKS__I2D_CoreFMM_PlanBuilderWim__) */
