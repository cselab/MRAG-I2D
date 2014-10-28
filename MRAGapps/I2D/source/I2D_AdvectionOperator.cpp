/*
 *  I2D_AdvectionStretchingOperator.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include "I2D_AdvectionOperator.h"
#include "I2D_ScalarBlockLab.h"

struct RHSUpwind3rdOrder
{
	Real t;
	Real Uinf[2];
	
	int stencil_start[3], stencil_end[3];
	
	RHSUpwind3rdOrder(Real t, Real Uinf[2]): t(t)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
		
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
	}
	
	RHSUpwind3rdOrder(const RHSUpwind3rdOrder & copy): t(copy.t)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
		
		Uinf[0] = copy.Uinf[0];
		Uinf[1] = copy.Uinf[1];
	}
	
	template<typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const 
	{
		const Real factor = -1./(6*info.h[0]);
		
		for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
			for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
			{
				const Real dwdx[2] = {
					2*lab(ix+1, iy) + 3*lab(ix, iy) - 6*lab(ix-1, iy) + lab(ix-2, iy),
					-lab(ix+2, iy) + 6*lab(ix+1, iy) - 3*lab(ix, iy) -2*lab(ix-1, iy)};
				
				const Real dwdy[2] = {
					2*lab(ix, iy+1) + 3*lab(ix, iy) - 6*lab(ix, iy-1) + lab(ix, iy-2),
					-lab(ix, iy+2) + 6*lab(ix, iy+1) - 3*lab(ix, iy) -2*lab(ix, iy-1)};
				
				const Real u[2] = {
					Uinf[0] + out(ix, iy).u[0],
					Uinf[1] + out(ix, iy).u[1]};
				
				out.external_data[iy][ix] = factor*(max(u[0], (Real)0)*dwdx[0] + min(u[0], (Real)0)*dwdx[1] + 
													max(u[1], (Real)0)*dwdy[0] + min(u[1], (Real)0)*dwdy[1]);
			}		
		
		//correct for the x-inlet
		if(info.index[0] == 0)
		{
			//overwrite the rhs with the correct values for ix=0 and ix=1
			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<2; ix++)
				{
					//BEGIN CORRECTION
					const Real dwdx[2] = {
						2*lab(ix+1, iy) + 3*lab(ix, iy) - 6*0 + 0,
						-lab(ix+2, iy) + 6*lab(ix+1, iy) - 3*lab(ix, iy) -2*0};
					//END CORRECTION
					
					const Real dwdy[2] = {
						2*lab(ix, iy+1) + 3*lab(ix, iy) - 6*lab(ix, iy-1) + lab(ix, iy-2),
						-lab(ix, iy+2) + 6*lab(ix, iy+1) - 3*lab(ix, iy) -2*lab(ix, iy-1)};
					
					const Real u[2] = {
						Uinf[0] + out(ix, iy).u[0],
						Uinf[1] + out(ix, iy).u[1],
					};
					
					out.external_data[iy][ix] = factor*(
														max(u[0], (Real)0)*dwdx[0] + min(u[0], (Real)0)*dwdx[1] + 
														max(u[1], (Real)0)*dwdy[0] + min(u[1], (Real)0)*dwdy[1] );
				}
		}
		
		if(info.index[0] == pow(2.0, info.level)-1)
		{
			//overwrite the rhs with the correct values for ix=FluidBlock2D::sizeY-1
			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
			{
				const int ix = FluidBlock2D::sizeY-1;
				
				//BEGIN CORRECTION
				const Real dwdx[2] = {
					2*0 + 3*lab(ix, iy) - 6*lab(ix-1, iy) + lab(ix-2, iy),
					-0 + 6*0 - 3*lab(ix, iy) -2*lab(ix-1, iy)};
				//END CORRECTION
				
				const Real dwdy[2] = {
					2*lab(ix, iy+1) + 3*lab(ix, iy) - 6*lab(ix, iy-1) + lab(ix, iy-2),
					-lab(ix, iy+2) + 6*lab(ix, iy+1) - 3*lab(ix, iy) -2*lab(ix, iy-1)};
				
				const Real u[2] = {
					Uinf[0] + out(ix, iy).u[0],
					Uinf[1] + out(ix, iy).u[1],
				};
				
				out.external_data[iy][ix] = factor*(
													max(u[0], (Real)0)*dwdx[0] + min(u[0], (Real)0)*dwdx[1] + 
													max(u[1], (Real)0)*dwdy[0] + min(u[1], (Real)0)*dwdy[1] );
			}
		}
	}
};

struct GetUMax
{
	Real Uinf[2];
	map< int, Real>& local_max_velocities;
	
	GetUMax(map< int, Real>& local_max_velocities, Real Uinf[2]):local_max_velocities(local_max_velocities) 
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
	}
	
	GetUMax(const GetUMax& c): local_max_velocities(c.local_max_velocities)
	{
		Uinf[0] = c.Uinf[0];
		Uinf[1] = c.Uinf[1];
	}
	
	inline void operator() (const BlockInfo& info, FluidBlock2D& b) const
	{
		Real maxVel[2] = {0.0,0.0};
		
		FluidElement2D * e = &b(0,0);
		
		static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
		for(int i=0; i<n; i++) 
		{
			maxVel[0] = max((Real)fabs(Uinf[0] + e[i].u[0]), maxVel[0]);
			maxVel[1] = max((Real)fabs(Uinf[1] + e[i].u[1]), maxVel[1]);
		}
		
		map< int, Real>::iterator it = local_max_velocities.find(info.blockID);
		assert(it != local_max_velocities.end());
		
		it->second = max(maxVel[0],maxVel[1]);
	}
};

Real I2D_AdvectionOperator::compute_maxvel()
{
	map<int, Real> velocities;
	
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	for(vector<BlockInfo>::iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		velocities[it->blockID] = 0;
	
	GetUMax get_velocities(velocities, Uinf);
	block_processing.process(vInfo, grid.getBlockCollection(), get_velocities);
	
	Real maxvel = 0;
	for(map< int, Real>::iterator it=velocities.begin(); it!=velocities.end(); it++)
		maxvel = max(maxvel, it->second);
	
	tmp_maxvel = maxvel;
	
	return maxvel;
}

Real I2D_AdvectionOperator::estimate_largest_dt()
{
	const Real maxvel = compute_maxvel();
	const Real max_dx = (1./B::sizeX)*pow(0.5,grid.getCurrentMinLevel());
	const Real min_dx = (1./B::sizeX)*pow(0.5,grid.getCurrentMaxLevel());
	
	largest_dt = max_dx/maxvel * CFL;
	smallest_dt = min_dx/maxvel * CFL;
	assert(largest_dt >= smallest_dt);
	state = Ready;
	
	return largest_dt;
}

struct AdvectionSimple
{
	Grid<W,B>& grid;
	BlockProcessing& block_processing;
	Real Uinf[2];
	int& rhscounter;
	
	AdvectionSimple(Grid<W,B>& grid, BlockProcessing& block_processing, Real Uinf[2], int& rhscounter_): 
	grid(grid), block_processing(block_processing), rhscounter(rhscounter_)
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
	}
	
	template< template<typename BT> class CorrectLab>
	void integrate(double dt)
	{
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		BoundaryInfo& binfo=grid.getBoundaryInfo();  
		const BlockCollection<B>& coll = grid.getBlockCollection();
		
		RHSUpwind3rdOrder advection_rhs(0, Uinf);
		block_processing.process< CorrectLab >(vInfo, coll, binfo, advection_rhs);
		UpdateScalarRK2<1> stepA(0, dt);
		block_processing.process(vInfo, coll, stepA);
		
		advection_rhs.t = dt;
		block_processing.process< CorrectLab >(vInfo, coll, binfo, advection_rhs);
		UpdateScalarRK2<2> stepB(dt, dt);
		block_processing.process(vInfo, coll, stepB);	
		
		rhscounter += vInfo.size()*2;
	}
};

struct AdvectionLTS
{
	Grid<W,B>& grid;
	BlockProcessing& block_processing;
	SpaceTimeSorter stsorter;
	Real Uinf[2];
	int& rhscounter;
	Real CFL, maxvel;
	const bool bSmartTrick;
	
	AdvectionLTS(Grid<W,B>& grid, BlockProcessing& block_processing, Real Uinf[2], int& rhscounter, Real CFL, Real maxvel, const bool bSmartTrick): 
	grid(grid), block_processing(block_processing), rhscounter(rhscounter), CFL(CFL), maxvel(maxvel), bSmartTrick(bSmartTrick)
	{
		stsorter.connect(grid);
		
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
	}
	
	~AdvectionLTS()
	{
		stsorter.disconnect();
	}
	
	template< template<typename BT> class CorrectLab >
	void integrate(double largest_dt)
	{
		int startlevel = grid.getCurrentMinLevel();
		
		vector<BlockInfo> vEasyBlocks;
		
		if (bSmartTrick)
		{
			const int LMIN = grid.getCurrentMinLevel();
			const int LMAX = grid.getCurrentMaxLevel();
			
			for(int lev=LMIN; lev<=LMAX; lev++)
			{
				const double my_dx = (1./B::sizeX)*pow(0.5,lev);
				const double my_dt = my_dx/maxvel * CFL;
				
				if (my_dt >= largest_dt)
				{
					vector<BlockInfo> v = grid.getBlocksInfoAtLevel(lev);
					vEasyBlocks.insert(vEasyBlocks.end(), v.begin(), v.end());
					
					startlevel = lev;
				}
				else 
					break;
			}
		}
		
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		BoundaryInfo& binfo=grid.getBoundaryInfo();  
		const BlockCollection<B>& coll = grid.getBlockCollection();
		
		stsorter.startSession(largest_dt, 2, 0, startlevel);
		
		while(true)
		{
			SpaceTimeSorter::ETimeInterval type;
			int level;
			double currTime, currDeltaT;
			const bool bContinue = stsorter.getBlocks(level, currDeltaT, currTime, vInfo, type);	
			
			if (bSmartTrick && level == startlevel)
				vInfo = vEasyBlocks;
			
			if (type == SpaceTimeSorter::ETimeInterval_Start)
			{
				RHSUpwind3rdOrder advection_rhs(currTime, Uinf);
				block_processing.process<CorrectLab>(vInfo, coll, binfo, advection_rhs);
				
				UpdateScalarRK2<1> stepA(currTime, currDeltaT);
				block_processing.process(vInfo, coll, stepA);
			}
			else if (type == SpaceTimeSorter::ETimeInterval_End)
			{
				RHSUpwind3rdOrder advection_rhs(currTime, Uinf);
				block_processing.process<CorrectLab>(vInfo, coll, binfo, advection_rhs);
				
				UpdateScalarRK2<2> stepB(currTime, currDeltaT);
				block_processing.process(vInfo, coll, stepB);
			}
			else
				abort();
			
			rhscounter += vInfo.size();
			
			if (!bContinue) break;
		}
		
		stsorter.endSession();
	}
};


void I2D_AdvectionOperator::perform_timestep(double dt)
{
	rhscounter = 0;
	
	if(state != Ready)
	{
		const double val =  estimate_largest_dt();
		if(dt > val)
		{
			printf("oops.  I2D_AdvectionOperator::perform_timestep: dt=%e > largest dt allowed %e\n", dt, val);
		}
		
		assert(dt<= val);
		state = Ready;
	}
	
	state = Done;
	
	const bool bUseLTS = dt > smallest_dt;
	
	if(bUseLTS)
	{
		printf("LTS!\n");		
		AdvectionLTS lts(grid, block_processing, Uinf, rhscounter, CFL, tmp_maxvel, true);		
		lts.integrate<I2D_ScalarBlockLab< Streamer_OmegaLTS >::Lab> (dt);		
		return;
	}
	
	AdvectionSimple advection(grid, block_processing, Uinf, rhscounter);	
	advection.integrate<I2D_ScalarBlockLab< Streamer_OmegaLTS >::Lab>(dt);
}

