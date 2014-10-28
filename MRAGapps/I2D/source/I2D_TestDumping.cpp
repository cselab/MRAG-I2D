/*
 *  I2D_TestDumping.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */


#include "I2D_TestDumping.h"
#include "I2D_CircularObstacleOperator.h"

ArgumentParser * gparser;

I2D_TestDumping::I2D_TestDumping(const int argc, const char ** argv):
	parser(argc, argv), refiner(2), compressor(2)
{
	gparser = new ArgumentParser(argc, argv);
	
	printf("////////////////////////////////////////////////////////////\n");
	printf("/////////////////     DUMPING TEST       ///////////////////\n");
	printf("////////////////////////////////////////////////////////////\n");
		   
	const int bpd = parser("-bpd").asInt(); assert(bpd > 1);
	
	grid = new Grid<W,B>(bpd,bpd);	assert(grid != NULL);
	grid->setRefiner(&refiner);
	grid->setCompressor(&compressor);
	
	_ic(*grid);	
	
	for(int i=0; i<3;i++)
	{
		set<int> boundary_blocks = _getBoundaryBlockIDs();
		refiner.set_blacklist(&boundary_blocks);
		Science::AutomaticRefinement<0,0>(*grid, fwt_omega, parser("-rtol").asDouble(), parser("-lmax").asInt(), 1, NULL, _ic, &boundary_blocks);
	}
	
	Science::AutomaticCompression<0,0>(*grid, fwt_omega, parser("-ctol").asDouble(), 1, NULL, _ic);
}

void I2D_TestDumping::run()
{	
	IO_VTKNative<W,B, 4,0> vtkdumper;
	vtkdumper.Write(*grid, grid->getBoundaryInfo(), "ic");
	
	exit(0);
}

void I2D_TestDumping::paint()
{
}

void I2D_TestDumping::_ic(Grid<W,B>& grid)
{
	const Real p[3] = {0.5,0.5,0.5};
	const Real radius = 0.25/5;
	const Real smoothing_length = 0.01;
	
	static I2D_CircularObstacleOperator * obstacle = new I2D_CircularObstacleOperator(grid, radius, p, smoothing_length);

	obstacle->characteristic_function();
	
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
			
			const Real my_radius = sqrt(pow(p[0]-0.5, 2) + pow(p[1]-0.5, 2));
			const bool bInside =  my_radius < radius; 
			
			b(ix, iy).omega = b(ix, iy).tmp;
		}
	}
}

set<int> I2D_TestDumping::_getBoundaryBlockIDs()
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