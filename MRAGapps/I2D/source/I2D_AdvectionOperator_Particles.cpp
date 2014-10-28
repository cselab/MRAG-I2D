/*
 *  I2D_AdvectionOperator_Particles.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include "I2D_AdvectionOperator_Particles.h"
#include "I2D_VectorBlockLab.h"
#include "I2D_ParticleBlockLab.h"
#include "I2D_GradOfVector.h"
#include "I2D_ParticleBlockLab.h"

struct GetGradUMax: I2D_GradOfVector_4thOrder
{
	map< int, Real>& local_max_gradu;
	Real t;
	int stencil_start[3], stencil_end[3];
	
	GetGradUMax(map< int, Real>& local_max_gradu):local_max_gradu(local_max_gradu), t(0)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	GetGradUMax(const GetGradUMax& c): local_max_gradu(c.local_max_gradu), t(0)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	struct TmpMax
	{ static inline void stream(FluidElement2D& out, Real in) { out.tmp = max((Real)fabs(in), out.tmp); } };
	
	template<typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const 
	{
		//clear all tmps
		{
			FluidElement2D * const e = &out(0,0);
			
			static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
			for(int i=0; i<n; i++) 
				e[i].tmp = 0;
		}
		
		//compute gradu and reduce to one component
		{
			_dfdx<TmpMax, 0 >(lab, info, out);
			_dfdy<TmpMax, 0 >(lab, info, out);
			
			_dfdx<TmpMax, 1 >(lab, info, out);
			_dfdy<TmpMax, 1 >(lab, info, out);
		}

		//reduce to one number
		{
			FluidElement2D * const e = &out(0,0);
			
			Real maxVal = 0;
			
			static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
			for(int i=0; i<n; i++) 
				maxVal = max((Real)fabs(e[i].tmp), maxVal);
			
			map< int, Real>::iterator it = local_max_gradu.find(info.blockID);
			assert(it != local_max_gradu.end());
			it->second = maxVal;
		}
	}
};

Real I2D_AdvectionOperator_Particles::_GTS_CFL()
{
	const Real maxvel = compute_maxvel();
	const Real min_dx = (1./B::sizeX)*pow(0.5,grid.getCurrentMaxLevel());
	
	return min_dx/maxvel * CFL;
}

Real I2D_AdvectionOperator_Particles::estimate_largest_dt()
{
	state = Ready;
	
	map<int, Real> max_grad_u;
	
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	for(vector<BlockInfo>::iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		max_grad_u[it->blockID] = 0;
	
	BoundaryInfo& binfo=grid.getBoundaryInfo();  
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	GetGradUMax compute_max_val(max_grad_u);
	block_processing.process< I2D_VectorBlockLab< Streamer_Velocity, 2 >::Lab >(vInfo, coll, binfo, compute_max_val);
		
	Real maxval = 0;
	for(map< int, Real>::iterator it=max_grad_u.begin(); it!=max_grad_u.end(); it++)
		maxval = max(maxval, it->second);
	
	Real dtLCFL = LCFL/maxval;
	Real dtCFL = _GTS_CFL();
	
	printf("dtCFL=%e, dtLCFL=%e\n", dtCFL, dtLCFL);
	return min(dtCFL, dtLCFL);
}

struct PushRemesh
{
	Real dt,t;
	Real Uinf[2];
	
	int stencil_start[3], stencil_end[3];
	
	PushRemesh(Real dt, Real Uinf[2]): t(0), dt(dt)
	{
		stencil_start[0] = stencil_start[1] = -3;
		stencil_end[0] = stencil_end[1] = +4;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
		
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
	}
	
	PushRemesh(const PushRemesh& c): t(0), dt(c.dt)
	{
		stencil_start[0] = stencil_start[1] = -3;
		stencil_end[0] = stencil_end[1] = +4;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
		
		Uinf[0] = c.Uinf[0];
		Uinf[1] = c.Uinf[1];
	}
	
	template<typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const 
	{
		lab.pcore.push(info, lab, dt, Uinf);
		lab.pcore.remesh(lab);
		
		for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					out.external_data[iy][ix] = lab.pcore.omega_new[iy][ix];
	}
};

struct UpdateOmega
{
	inline void operator() (const BlockInfo& info, FluidBlock2D& b) const
	{	
		FluidElement2D * const dest = &b(0,0);
		
		const Real * const src = &b.external_data[0][0];
		
		const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
		
		for(int i=0; i<n; i++)
			dest[i].omega = src[i];
	}
};

void I2D_AdvectionOperator_Particles::perform_timestep(double dt)
{
	assert(state == Ready);
	state = Done;
	
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	BoundaryInfo& binfo=grid.getBoundaryInfo();  
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	PushRemesh push_remesh(dt, Uinf);
	block_processing.process< I2D_ParticleBlockLab >(vInfo, coll, binfo, push_remesh);
	
	UpdateOmega update;
	block_processing.process(vInfo, coll, update);
	
	rhscounter = vInfo.size();
}