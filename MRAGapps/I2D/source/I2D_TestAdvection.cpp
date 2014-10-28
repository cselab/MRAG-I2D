/*
 *  I2D_TestAdvectionStretching2D.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include "I2D_TestAdvection.h"
#include "I2D_AdvectionOperator_Particles.h"

static const int maxParticleStencil[2][3] = {
		-3, -3, 0,
		+4, +4, 1
};

I2D_TestAdvection::I2D_TestAdvection(const int argc, const char ** argv): parser(argc, argv), step_id(0)
{
	printf("////////////////////////////////////////////////////////////\n");
	printf("//////////////////       ADVECTION TEST     ////////////////\n");
	printf("////////////////////////////////////////////////////////////\n");

	parser.set_strict_mode();
	const int bpd = parser("-bpd").asInt();
	assert(bpd > 1);

	if (parser("-particles").asBool())
		grid = new Grid<W,B>(bpd,bpd,1, maxParticleStencil);
	else
		grid = new Grid<W,B>(bpd,bpd,1);

	assert(grid != NULL);

	const int res_jump = max(1, parser("-jump").asInt());
	const int lmax = parser("-lmax").asInt();

	refiner = new Refiner(res_jump, lmax);
	compressor = new Compressor(res_jump);
	grid->setRefiner(refiner);
	grid->setCompressor(compressor);

	_ic_omega(*grid);	

	_refine(true);
	_compress(true);

	//if (parser("-dumpfreq").asInt() > 0 ) _dump("ic_advection");

	if (parser("-particles").asBool())
		advection = new I2D_AdvectionOperator_Particles(*grid);
	else
		advection = new I2D_AdvectionOperator(*grid);

	Real Uinf[2] = {0,0};
	advection->set_Uinfinity(Uinf);
}


void I2D_TestAdvection::run()
{	
	char buf[500];
	sprintf(buf, "advection_%04d", step_id);
	//if (step_id%parser("-dumpfreq").asInt()==0 ) _dump(buf);

	_refine(false);

	double dt = advection->estimate_largest_dt();
	printf("Maximum dt is :%e **********************\n", dt);
	advection->perform_timestep(dt);

	_compress(false);

	printf("STEP ID: %d\n", step_id++);
}

void I2D_TestAdvection::paint()
{
}

set<int> I2D_TestAdvection::_getBoundaryBlockIDs()
{
	set<int> result;

	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	for(vector<BlockInfo>::iterator it = vInfo.begin(); it != vInfo.end(); it++)
	{
		const bool bX = it->index[0] == 0 || it->index[0] == pow(2, it->level)-1;
		const bool bY = it->index[1] == 0 || it->index[1] == pow(2, it->level)-1;

		if (bX || bY ) result.insert(it->blockID);
	}

	return result;
}

void I2D_TestAdvection::_refine(bool bUseIC)
{
	if (!parser("-uniform").asBool())
	{
		while(true)
		{
			set<int> boundary_blocks = _getBoundaryBlockIDs();

			void (*initial_condition)(Grid<W,B>&) = bUseIC ? _ic_omega : NULL;
			const int refinements = Science::AutomaticRefinement<0,0>(*grid, fwt_omega, parser("-rtol").asDouble(), parser("-lmax").asInt(), 1, NULL, initial_condition);//, &boundary_blocks);

			if (refinements == 0) break;
		}
	}

	_ic_velocity(*grid);
}

void I2D_TestAdvection::_compress(bool bUseIC)
{
	void (*initial_condition)(Grid<W,B>&) = bUseIC ? _ic_omega : NULL;

	if (!parser("-uniform").asBool())
		Science::AutomaticCompression<0,0>(*grid, fwt_omega, parser("-ctol").asDouble(), 1, NULL, initial_condition);

	_ic_velocity(*grid);
}

namespace AdvectionIC {

void _map_space(Real& x, Real& y)
{
}

void _compute_velocity(Real x, Real y, Real& u, Real& v)
{
	u = 2*pow(sin(M_PI*x),2)*sin(2*M_PI*y);
	v = -sin(2*M_PI*x)*pow(sin(M_PI*y),2);

	/*u = 0;
		 v = 0;
		 w = 0;*/
}

}

void I2D_TestAdvection::_ic_omega(Grid<W,B>& grid)
{
	printf("IC omega\n");

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
				AdvectionIC::_map_space(p[0], p[1]);

				const Real radius = sqrt(pow(p[0]-0.35, 2) + pow(p[1]-0.35, 2));
				const bool bInside =  radius < 0.025;

				b(ix, iy).omega = bInside;
			}
	}
}

void I2D_TestAdvection::_ic_velocity(Grid<W,B>& grid)
{
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
					AdvectionIC::_map_space(p[0], p[1]);

					Real u,v,w;
					AdvectionIC::_compute_velocity(p[0], p[1], u,v);

					b(ix, iy, iz).u[0] = u;
					b(ix, iy, iz).u[1] = v;
				}
	}	
}

void I2D_TestAdvection::_dump(string filename)
{
	IO_VTKNative<W,B, 4,0> vtkdumper;
	vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
}
