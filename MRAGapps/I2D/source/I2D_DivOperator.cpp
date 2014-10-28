/*
 *  I2D_DivOperator.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include "I2D_DivOperator.h"

#include "I2D_VectorBlockLab.h"
#include "I2D_GradOfVector.h"

struct DivOp: I2D_GradOfVector_4thOrder
{
	Real t;
	int stencil_start[3], stencil_end[3];
	
	DivOp(): t(0)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	DivOp(const DivOp & copy): t(0)
	{
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -2;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	struct TmpAdd { static inline void stream(FluidElement2D& out, Real in) { out.tmp += in; } };
	
	template<typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const 
	{
		//clear all tmp component
		{
			FluidElement2D * const e = &out(0,0);
			
			static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
			for(int i=0; i<n; i++) 
				e[i].tmp = 0;
		}
		
		_dfdx<TmpAdd, 0 >(lab, info, out);
		_dfdy<TmpAdd, 1 >(lab, info, out);
	}
};

void I2D_DivOperator::perform()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	BoundaryInfo& binfo=grid.getBoundaryInfo();  
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	DivOp div;
	block_processing.process< I2D_VectorBlockLab< Streamer_Velocity, 2 >::Lab >(vInfo, coll, binfo, div);
}

