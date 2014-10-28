/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include "I2D_CurlVelocityOperator.h"
#include "I2D_VectorBlockLab.h"
#include "I2D_GradOfVector.h"

struct CurlVel_4thOrder: I2D_GradOfVector_4thOrder
{
	Real t;
	int stencil_start[3], stencil_end[3];
	
	CurlVel_4thOrder(): t(0)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	CurlVel_4thOrder(const CurlVel_4thOrder & copy): t(0)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}


	struct VortAdd
	{ static inline void stream(FluidElement2D& out, Real in) { out.omega += in; } };
	
	struct VortSub
	{ static inline void stream(FluidElement2D& out, Real in) { out.omega -= in; } };
	

	template<typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const 
	{
		//clear all vorticity component
		{
			FluidElement2D * const e = &out(0,0);
			
			static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
			for(int i=0; i<n; i++) 
				e[i].omega = 0.0;
		}
		
		_dfdx<VortAdd, 1 >(lab, info, out);
		_dfdy<VortSub, 0 >(lab, info, out);
	}
};


void I2D_CurlVelocityOperator_4thOrder::perform()
{	
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	BoundaryInfo& binfo=grid.getBoundaryInfo();  
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	CurlVel_4thOrder curl_vel;
	block_processing.process< I2D_VectorBlockLab< Streamer_Velocity, 2 >::Lab >(vInfo, coll, binfo, curl_vel);
}

