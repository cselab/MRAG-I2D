/*
 *  I2D_TestPenalizationAndOther.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include "I2D_TestPenalizationAndOther.h"

I2D_TestPenalizationAndOther::I2D_TestPenalizationAndOther(const int argc, const char ** argv): parser(argc, argv), step_id(0)
{
	printf("////////////////////////////////////////////////////////////\n");
	printf("/////////////////     PENALIZATION TEST     ////////////////\n");
	printf("////////////////////////////////////////////////////////////\n");
	
	const int bpd = parser("-bpd").asInt();
	assert(bpd > 1);
	
	grid = new Grid<W,B>(bpd,bpd,1);
	assert(grid != NULL);
	
	const int res_jump = max(1, parser("-jump").asInt());
	
	refiner = new Refiner_BlackList(res_jump);
	compressor = new Compressor(res_jump);
	grid->setRefiner(refiner);
	grid->setCompressor(compressor);
	
	_ic_omega(*grid);	
	
	_refine(true);
	_compress(true);
	
	_ic_omega(*grid);
	
	
	const Real lambda = 1e4;
	const Real Uinf[2] = {0.1,0};
	penalization = new I2D_PenalizationOperator(*grid, lambda, Uinf);
	
	div_omega = new I2D_DivOperator(*grid);
	
	Real pos[2] = {0.45,0.5};
	obstacle = new I2D_CircularObstacleOperator(*grid, 1/16., pos, sqrt(2.)*2./FluidBlock2D::sizeX*pow(0.5, parser("-lmax").asInt()));
	obstacle->characteristic_function();
	//
	
	
	
	//_dump("fadeout");
	
	
	if (parser("-dumpfreq").asInt() > 0 ) _dump("ic_penalization");
	
	//	exit(0);
	
}


void I2D_TestPenalizationAndOther::run()
{		
	_refine(false);
	
	double dt = 1e-1;
	printf("Maximum dt is :%e **********************\n", dt);
	
	obstacle->characteristic_function();
	penalization->perform_timestep(dt);
	obstacle->characteristic_function();
	double cor[2] = {0.5,0.5};
	//penalization->compute_dragandstuff(0, obstacle->getD(), cor,"prova.txt");
	
	div_omega->perform();
	//penalization->perform_timestep(dt);
	
	
	
	_compress(false);
	
	printf("STEP ID: %d\n", step_id++);
	
	if ( step_id % max(1, parser("-dumpfreq").asInt()) == 0)
	{
		char buf[500];
		sprintf(buf, "penalization_%04d", step_id);
		
		obstacle->characteristic_function();
		if (parser("-dumpfreq").asInt() > 0 ) _dump(buf);
	}
	
	
}

void I2D_TestPenalizationAndOther::paint()
{
}

set<int> I2D_TestPenalizationAndOther::_getBoundaryBlockIDs()
{
	set<int> result;
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	for(vector<BlockInfo>::iterator it = vInfo.begin(); it != vInfo.end(); it++)
	{
		const bool bX = it->index[0] == 0;// || it->index[0] == pow(2, it->level)-1;
		const bool bY = false;//it->index[1] == 0 || it->index[1] == pow(2, it->level)-1;
		
		if (bX || bY) result.insert(it->blockID);
	}
	
	return result;
}

void I2D_TestPenalizationAndOther::_refine(bool bUseIC)
{
	if (!parser("-uniform").asBool())
	{
		while(true)
		{
			set<int> boundary_blocks = _getBoundaryBlockIDs();
			((Refiner_BlackList*)refiner)->set_blacklist(&boundary_blocks);
			void (*initial_condition)(Grid<W,B>&) = bUseIC ? _ic_omega : NULL;
			const int refinements = Science::AutomaticRefinement<0,0>(*grid, fwt_obstacle, parser("-rtol").asDouble(), parser("-lmax").asInt(), 1, NULL, initial_condition, &boundary_blocks);
			printf("refinements = %d\n", refinements);
			if (refinements == 0) break;
		}
	}
	
	if(bUseIC)
		_ic_velocity(*grid);
}

void I2D_TestPenalizationAndOther::_compress(bool bUseIC)
{
	void (*initial_condition)(Grid<W,B>&) = bUseIC ? _ic_omega : NULL;
	
	if (!parser("-uniform").asBool())
		Science::AutomaticCompression<0,0>(*grid, fwt_obstacle, parser("-ctol").asDouble(), 1, NULL, initial_condition);
	
	if(bUseIC)
		_ic_velocity(*grid);
}

namespace ICpenalization
{
	void _map_space(Real& x, Real& y)
	{
		x = 4*x - 1.5;
		y = 4*y - 1.5;
	}
	
	void _compute_velocity(Real x, Real y, Real& u, Real& v)
	{
		u = 0;//2*pow(sin(M_PI*x),2)*sin(2*M_PI*y)*sin(2*M_PI*z);
		v = 0;//-sin(2*M_PI*x)*pow(sin(M_PI*y),2)*sin(2*M_PI*z);
	}
}
void I2D_TestPenalizationAndOther::_ic_omega(Grid<W,B>& grid)
{
	printf("IC omega\n");
	
	Real pos[2] = {0.45,0.5};
	I2D_CircularObstacleOperator obstacle(grid, 1/16., pos, sqrt(2.)*2./FluidBlock2D::sizeX*pow(0.5, 5));
	obstacle.characteristic_function();
	//exit(0);
	
	/*
	 vector<BlockInfo> vInfo = grid.getBlocksInfo();
	 
	 for(int i=0; i<vInfo.size(); i++)
	 {
	 BlockInfo info = vInfo[i];
	 B& b = grid.getBlockCollection()[vInfo[i].blockID];
	 
	 for(int iz=0; iz<B::sizeZ; iz++)
	 for(int iy=0; iy<B::sizeY; iy++)
	 for(int ix=0; ix<B::sizeX; ix++)
	 {
	 Real p[3];
	 
	 info.pos(p, ix, iy, iz);
	 ICpenalization::_map_space(p[0], p[1], p[2]);
	 
	 const Real radius = sqrt(pow(p[0]-0.35, 2) + pow(p[1]-0.35, 2) + pow(p[2]-0.35, 2));
	 const bool bInside =  radius < 0.15; 
	 
	 b(ix, iy, iz).omega[0] = bInside;
	 b(ix, iy, iz).omega[1] = 1;//sqrt(pow(p[0]-0.25, 2) + pow(p[1]-0.25, 2) + pow(p[2]-0.25, 2))<0.1;
	 b(ix, iy, iz).omega[2] = 0;//p[2];
	 }
	 }*/
}


void I2D_TestPenalizationAndOther::_ic_velocity(Grid<W,B>& grid)
{
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
				ICpenalization::_map_space(p[0], p[1]);
				
				Real u,v,w;
				ICpenalization::_compute_velocity(p[0], p[1], u,v);
				
				b(ix, iy).u[0] = u;
				b(ix, iy).u[1] = v;
			}
	}	
}

void I2D_TestPenalizationAndOther::_dump(string filename)
{
	IO_VTKNative<W,B, 4,0> vtkdumper;
	vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
}
