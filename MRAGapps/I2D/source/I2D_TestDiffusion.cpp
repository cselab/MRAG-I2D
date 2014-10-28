/*
 *  I2D_TestDiffusion.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include "I2D_TestDiffusion.h"

I2D_TestDiffusion::I2D_TestDiffusion(const int argc, const char ** argv):
parser(argc, argv), step_id(0)
{
	printf("////////////////////////////////////////////////////////////\n");
	printf("///////////////     DIFFUSION TEST       ///////////////////\n");
	printf("////////////////////////////////////////////////////////////\n");
	
	const int bpd = parser("-bpd").asInt();
	assert(bpd > 1);
	
	grid = new Grid<W,B>(bpd,bpd);
	assert(grid != NULL);
	
	const int res_jump = max(1, parser("-jump").asInt());
	const int lmax = max(1, parser("-lmax").asInt());
	
	refiner = new Refiner_BlackList(res_jump, lmax);
	
	grid->setRefiner(refiner);
	grid->setCompressor(&compressor);
	
	_ic(*grid);	
	
	if (!parser("-uniform").asBool())
	{
		while(true)
		{
			set<int> boundary_blocks =  _getBoundaryBlockIDs();
			((Refiner_BlackList*)refiner)->set_blacklist(&boundary_blocks);
			const int refinements = Science::AutomaticRefinement<0,0>(*grid, fwt_omega, parser("-rtol").asDouble(), parser("-lmax").asInt(), 1, NULL, _ic, &boundary_blocks);
			
			if (refinements == 0) break;
		}
	}
	
	if (parser("-dumpfreq").asInt() > 0 ) _dump("ic");
	
	diffusion = new I2D_DiffusionOperator_4thOrder(*grid, parser("-nu").asDouble(), parser("-fc").asDouble());
}

void I2D_TestDiffusion::_refine()
{
	if (!parser("-uniform").asBool())
	{
		while(true)
		{
			set<int> boundary_blocks = _getBoundaryBlockIDs();
			((Refiner_BlackList*)refiner)->set_blacklist(&boundary_blocks);
			const int refinements = Science::AutomaticRefinement<0,0>(*grid, fwt_omega, parser("-rtol").asDouble(), parser("-lmax").asInt(), 1, NULL, (void (*)(Grid<W,B>&)) NULL, &boundary_blocks);
			
			if (refinements == 0) break;
		}
	}
}

void I2D_TestDiffusion::run()
{	
	//I2D_ScalarLaplaceOperator lapl(*grid);
	//lapl.perform();
	//_dump("lapl");
	
	//exit(0);

	if ( step_id % max(1, parser("-dumpfreq").asInt()) == 0)
	{
		char buf[500];
		sprintf(buf, "diffusion_%04d", step_id);

		if (parser("-dumpfreq").asInt() > 0 ) _dump(buf);
	}
	
//	_refine();
	
	double dt = diffusion->estimate_largest_dt();
	printf("TIMESTEP IS %e\n", dt);
	diffusion->perform_timestep(dt);
	
	Science::AutomaticCompression<0,0>(*grid, fwt_omega, parser("-ctol").asDouble(), 1);
	
	printf("STEP ID: %d\n", step_id++);
}

void I2D_TestDiffusion::paint()
{
}

void I2D_TestDiffusion::_ic(Grid<W,B>& grid)
{
	printf("IC\n");
	
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		B& b = grid.getBlockCollection()[vInfo[i].blockID];
		
		for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
			{
				Real p[2];
				
				info.pos(p, ix, iy);
				const Real x = p[0];
				const Real y = p[1];
				
				const Real radius = sqrt(pow(p[0]-0.35, 2) + pow(p[1]-0.35, 2) );
				const bool bInside =  radius < 0.1; 
				
				b(ix, iy).omega = bInside;//z*z + y + x*x + z + y*y + x;
				b(ix, iy).tmp = 6;
			}
	}
}

set<int> I2D_TestDiffusion::_getBoundaryBlockIDs()
{
	set<int> result;
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	for(vector<BlockInfo>::iterator it = vInfo.begin(); it != vInfo.end(); it++)
	{
		const bool bX = it->index[0] == 0 || it->index[0] == pow(2, it->level)-1;
		const bool bY = it->index[1] == 0 || it->index[1] == pow(2, it->level)-1;
		
		if (bX || bY) result.insert(it->blockID);
	}
	
	return result;
}

void I2D_TestDiffusion::_dump(string filename)
{
	IO_VTKNative<W,B, 4,0> vtkdumper;
	vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
}