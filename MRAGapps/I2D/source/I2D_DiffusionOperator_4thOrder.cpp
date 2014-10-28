/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include "I2D_DiffusionOperator.h"
#include "I2D_ScalarBlockLab.h"
#include "I2D_DivGradOfScalar.h"

struct DiffusionRHS_4thOrder: I2D_DivGradOfScalar_4thOrder
{
	Real t;
	Real viscosity;
	int stencil_start[3], stencil_end[3];
	
	DiffusionRHS_4thOrder(Real viscosity, Real t): viscosity(viscosity), t(t)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	DiffusionRHS_4thOrder(const DiffusionRHS_4thOrder & copy): viscosity(copy.viscosity), t(copy.t)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	struct AddToTmp { static inline void stream(Real& out, Real in) { out += in; } };
	
	template<typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const 
	{
		//clear all tmps
		Real * const e = &out.external_data[0][0];
		
		static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
		
		for(int i=0; i<n; i++) 
			e[i] = 0;
		
		_ddpsiddx_ptr< AddToTmp >(lab, info, e);
		_ddpsiddy_ptr< AddToTmp >(lab, info, e);
		
		for(int i=0; i<n; i++)
			e[i] *= viscosity;
	}
};

double I2D_DiffusionOperator_4thOrder::estimate_largest_dt() const
{
	double max_dx = (1./B::sizeX)*pow(0.5,grid.getCurrentMinLevel());
	
	return FC * pow(max_dx,2)/(viscosity * 2.0 * DIM);
}

double I2D_DiffusionOperator_4thOrder::estimate_smallest_dt() const
{
	double max_dx = (1./B::sizeX)*pow(0.5,grid.getCurrentMaxLevel());
	
	return FC * pow(max_dx,2)/(viscosity * 2.0 * DIM);
}

struct DiffusionSimple_4thOrder
{
	Grid<W,B>& grid;
	BlockProcessing& block_processing;
	int& rhscounter;

	DiffusionSimple_4thOrder(Grid<W,B>& grid, BlockProcessing& block_processing, int& rhscounter): 
	grid(grid), block_processing(block_processing), rhscounter(rhscounter){}
	
	template< template<typename BT> class CorrectLab >
	void integrate(double viscosity, double dt)
	{
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		BoundaryInfo& binfo=grid.getBoundaryInfo();  
		const BlockCollection<B>& coll = grid.getBlockCollection();
		
		DiffusionRHS_4thOrder diffusion_rhs(viscosity, 0);
		block_processing.process< CorrectLab >(vInfo, coll, binfo, diffusion_rhs);
		UpdateScalarRK2<1> stepA(0, dt);
		block_processing.process(vInfo, coll, stepA);
		
		diffusion_rhs.t = dt;
		block_processing.process< CorrectLab >(vInfo, coll, binfo, diffusion_rhs);
		UpdateScalarRK2<2> stepB(dt, dt);
		block_processing.process(vInfo, coll, stepB);	
		
		rhscounter += vInfo.size()*2;
	}
};

struct DiffusionLTS_4thOrder
{
	const double FC, DIM;
	Grid<W,B>& grid;
	BlockProcessing& block_processing;
	SpaceTimeSorter stsorter;
	const bool bSmartTrick;
	int& rhscounter;
	
	DiffusionLTS_4thOrder(Grid<W,B>& grid, BlockProcessing& block_processing, bool bSmartTrick, int& rhscounter, double FC):
	grid(grid), block_processing(block_processing), bSmartTrick(bSmartTrick), rhscounter(rhscounter), FC(FC), DIM(2.0)
	{
		stsorter.connect(grid);
	}
	
	~DiffusionLTS_4thOrder()
	{
		stsorter.disconnect();
	}
	
	template< template<typename BT> class CorrectLab >
	void integrate(double viscosity, double largest_dt)
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
				const double my_dt = FC * pow(my_dx,2)/(viscosity * 2.0 * DIM);
				
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
		
		stsorter.startSession(largest_dt, 4, 0, startlevel);
		
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
				DiffusionRHS_4thOrder diffusion_rhs(viscosity, currTime);
				block_processing.process<CorrectLab>(vInfo, coll, binfo, diffusion_rhs);
					
				UpdateScalarRK2<1> stepA(currTime, currDeltaT);
				block_processing.process(vInfo, coll, stepA);
			}
			else if (type == SpaceTimeSorter::ETimeInterval_End)
			{
				DiffusionRHS_4thOrder diffusion_rhs(viscosity, currTime);
				block_processing.process<CorrectLab>(vInfo, coll, binfo, diffusion_rhs);
				
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

void I2D_DiffusionOperator_4thOrder::perform_timestep(double dt)
{
	rhscounter = 0;
	
	const bool bUseLTS = dt > estimate_smallest_dt();
	
	if (bUseLTS)
	{
		printf("diffusion with LTS!\n");		
		DiffusionLTS_4thOrder lts(grid, block_processing, true, rhscounter,FC);
		lts.integrate<I2D_ScalarBlockLab< Streamer_OmegaLTS >::Lab> (viscosity, dt);		
		return;
	}
	
	DiffusionSimple_4thOrder step(grid, block_processing, rhscounter);
	step.integrate<I2D_ScalarBlockLab<Streamer_OmegaLTS >::Lab >(viscosity, dt);
}


