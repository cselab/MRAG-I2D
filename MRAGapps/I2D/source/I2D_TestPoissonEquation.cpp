/*
 *  I2D_TestPoissonEquation.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include "I2D_TestPoissonEquation.h"

#include "I2D_VelocitySolver_Mani.h"

#ifdef _I2D_MPI_
#include "I2D_VelocitySolverMPI_Mani.h"
#endif

I2D_TestPoissonEquation::I2D_TestPoissonEquation(const int argc, const char ** argv): parser(argc, argv)
{
	
	printf("////////////////////////////////////////////////////////////\n");
	printf("////////////       POISSON EQUATION TEST     ///////////////\n");
	printf("////////////////////////////////////////////////////////////\n");
	
	const int bpd = parser("-bpd").asInt();
	assert(bpd > 1);
	
	const string fmm_name = parser("-fmm").asString();
#ifdef _I2D_MPI_
	if(fmm_name == "mpi-velocity")
		poisson_solver = new I2D_VelocitySolverMPI_Mani(argc, argv);
#endif
	
	grid = new Grid<W,B>(bpd,bpd,1);
	assert(grid != NULL);
	
	const int res_jump = max(1, parser("-jump").asInt());
	const int lmax = max(1, parser("-lmax").asInt());
	
	refiner = new Refiner_BlackList(res_jump, lmax);
	compressor = new Compressor(res_jump);
	grid->setRefiner(refiner);
	grid->setCompressor(compressor);
	
	_ic_omega(*grid);	
	
	_refine_omega(true);
	_compress(true);
	_dump("poissoneq_ic");
	
	if(fmm_name == "velocity")
		poisson_solver = new I2D_VelocitySolver_Mani(*grid, parser);
	else 
	{
		if(poisson_solver == NULL)
			abort();
#ifdef _I2D_MPI_
		else 	if(fmm_name == "mpi-velocity")
			((I2D_VelocitySolverMPI_Mani *) poisson_solver)->set_grid(*grid);
#endif
	}
	
	curl_velocity = new I2D_CurlVelocityOperator_4thOrder(*grid);
	//divergence = new I2D_DivOperator(*grid);
	//laplace_psi = new I2D_ScalarLaplaceOperator(*grid);
}


void I2D_TestPoissonEquation::run()
{	
	_dump("ic");
	poisson_solver->compute_velocity();
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		B& b = grid->getBlockCollection()[vInfo[i].blockID];
		
		for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
				{
					b(ix,iy).omega = 0.0;
				}
	}
	
	curl_velocity->perform();
	_dump("after");
	
	//_dump("divergence-free1");	
	//poisson_solver->compute_velocity();	
	//_dump("divergence-free2");

	printf("TEST CONCLUDED!\n");
	exit(0);
}

set<int> I2D_TestPoissonEquation::_getBoundaryBlockIDs()
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

void I2D_TestPoissonEquation::_refine_omega(bool bUseIC)
{
	if (!parser("-uniform").asBool())
	{
		while(true)
		{
			set<int> boundary_blocks = _getBoundaryBlockIDs();
			((Refiner_BlackList*)refiner)->set_blacklist(&boundary_blocks);
			void (*initial_condition)(Grid<W,B>&) = bUseIC ? _ic_omega : NULL;
			const int refinements = Science::AutomaticRefinement<0,0>(*grid, fwt_omega, parser("-rtol").asDouble(), parser("-lmax").asInt(), 1, NULL, initial_condition, &boundary_blocks);
			
			if (refinements == 0) break;
		}
	}
}

void I2D_TestPoissonEquation::_refine_vel(bool bUseIC)
{
	if (!parser("-uniform").asBool())
	{
		while(true)
		{
			set<int> boundary_blocks = _getBoundaryBlockIDs();
			((Refiner_BlackList*)refiner)->set_blacklist(&boundary_blocks);
			void (*initial_condition)(Grid<W,B>&) = bUseIC ? _ic_omega : NULL;
			const int refinements = Science::AutomaticRefinement<0,1>(*grid, fwt_velocity, parser("-rtol").asDouble(), parser("-lmax").asInt(), 1, NULL, initial_condition, &boundary_blocks);
			
			if (refinements == 0) break;
		}
	}
}

void I2D_TestPoissonEquation::_compress(bool bUseIC)
{
	void (*initial_condition)(Grid<W,B>&) = bUseIC ? _ic_omega : NULL;
	
	if (!parser("-uniform").asBool())
		Science::AutomaticCompression<0,0>(*grid, fwt_omega, parser("-ctol").asDouble(), 1, NULL, initial_condition);
}

namespace PoissonEquationIC {
	Real _fill(Real x, Real y);
	
	class BS4
	{
	public:
		static inline Real eval(Real x) 
		{
			const Real t = fabs(x);
			
			if (t>2) return 0;
			
			if (t>1) return pow(2-t,3)/6;
			
			return (1 + 3*(1-t)*(1 + (1-t)*(1 - (1-t))))/6;
		}
	};
	
	Real _fill(Real x, Real y)
	{
		const Real r = sqrt(pow(x-0.35, 2) + pow(y-0.35, 2) );
		return BS4::eval(r/0.15);
	}
}

void I2D_TestPoissonEquation::paint()
{
}

void I2D_TestPoissonEquation::_ic_omega(Grid<W,B>& grid)
{
	printf("IC omega\n");
	
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		B& b = grid.getBlockCollection()[vInfo[i].blockID];
		
		for(int iz=0; iz<B::sizeZ; iz++)
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
				{
					Real p[2];
					
					info.pos(p, ix, iy);
					const Real x = p[0];
					const Real y = p[1];
					//PoissonEquationIC::_map_space(p[0], p[1]);
					
					b(ix, iy, iz).omega = PoissonEquationIC::_fill(p[0], p[1]);;//z;//_fill(p[0], p[1], p[2]);//p[2];
					b(ix, iy, iz).u[0] = 0.0;//x*x*x;//PoissonEquationIC::_fill(p[0], p[1], p[2]);
					b(ix, iy, iz).u[1] = 0.0;//y*y;//PoissonEquationIC::_fill(p[0], p[1], p[2]);//sqrt(pow(p[0]-0.25, 2) + pow(p[1]-0.25, 2) + pow(p[2]-0.25, 2))<0.1;
					b(ix, iy, iz).tmp = -313;
				}
	}
}


void I2D_TestPoissonEquation::_dump(string filename)
{
	IO_VTKNative<W,B, 4,0> vtkdumper;
	vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
}

void I2D_TestPoissonEquation::_load_ic(string filename)
{	
	//read status
	{
		FILE * f = fopen((filename + ".status").c_str(), "r");
		assert(f != NULL);
		float val = -1;
		fscanf(f, "time: %e\n", &val);
		assert(val>=0);
		//t=val;
		int step_id = -1;
		fscanf(f, "stepid: %d\n", &step_id);
		assert(step_id >= 0);
		fclose(f);
	}
	
	printf("DESERIALIZATION: time is %f and step id is %d\n", 0, 0);
	
	//read grid
	IO_Binary<W,B> serializer;
	serializer.Read(*grid, filename);
}
