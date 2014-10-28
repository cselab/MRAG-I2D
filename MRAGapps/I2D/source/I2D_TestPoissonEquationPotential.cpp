/*
 *  I2D_TestPoissonEquationPotential.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Mattia Gazzola on 11/05/11.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include "I2D_TestPoissonEquationPotential.h"

#include "I2D_PotentialSolver_Mattia.h"

#ifdef _I2D_MPI_
#include "I2D_VelocitySolverMPI_Mani.h"
#endif

namespace PoissonEquationPotIC {

Real eval(Real x)
{
	const Real t = fabs(x);
	if (t>2) return 0;
	if (t>1) return pow(2-t,3)/6;
	return (1 + 3*(1-t)*(1 + (1-t)*(1 - (1-t))))/6;
}

Real _fill(Real x, Real y)
{
	const Real r = sqrt(pow(x-0.5, 2) + pow(y-0.5, 2) );
	return 0.1*eval(r/0.05);
}

}


I2D_TestPoissonEquationPotential::I2D_TestPoissonEquationPotential(const int argc, const char ** argv): parser(argc, argv)
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

	divergence = new I2D_DivOperator(*grid);
	curl_velocity = new I2D_CurlVelocityOperator_4thOrder(*grid);

	refiner = new Refiner_BlackList(res_jump, lmax);
	compressor = new Compressor(res_jump);
	grid->setRefiner(refiner);
	grid->setCompressor(compressor);

	_ic(*grid);
	_refine(true);

	if(fmm_name == "velocity")
		poisson_solver = new I2D_PotentialSolver_Mattia(*grid, parser);
	else 
	{
		if(poisson_solver == NULL)
			abort();
#ifdef _I2D_MPI_
		else 	if(fmm_name == "mpi-velocity")
			((I2D_VelocitySolverMPI_Mani *) poisson_solver)->set_grid(*grid);
#endif
	}


}

set<int> I2D_TestPoissonEquationPotential::_getBoundaryBlockIDs()
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

void I2D_TestPoissonEquationPotential::_refine(bool bUseIC)
{
	if (parser("-uniform").asBool()) return;

	if (bUseIC)
	{
		while(true)
		{
			const int refinements = Science::AutomaticRefinement<0,3>(*grid, fwt_omegaANDvelocityANDchi, parser("-rtol").asDouble(), parser("-lmax").asInt(), 1, NULL, (void (*)(Grid<W,B>&))NULL, NULL);
			_ic(*grid);
			if (refinements == 0) break;
		}
	}
	else
	{
		while(true)
		{
			const int refinements = Science::AutomaticRefinement<0,3>(*grid, fwt_omegaANDvelocityANDchi, parser("-rtol").asDouble(), parser("-lmax").asInt(), 1, NULL, (void (*)(Grid<W,B>&))NULL, NULL);
			if (refinements == 0) break;
		}
	}
}

void I2D_TestPoissonEquationPotential::_compress(bool bUseIC)
{
	if (parser("-uniform").asBool()) return;
	Science::AutomaticCompression<0,3>(*grid, fwt_omegaANDvelocityANDchi, parser("-ctol").asDouble(), 1, NULL, (void (*)(Grid<W,B>&))NULL);
}

void I2D_TestPoissonEquationPotential::paint()
{
}

void I2D_TestPoissonEquationPotential::_ic(Grid<W,B>& grid)
{
	printf("IC velocity\n");

	vector<BlockInfo> vInfo = grid.getBlocksInfo();

	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		B& b = grid.getBlockCollection()[vInfo[i].blockID];

		for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
			{
				Real p[2];

				info.pos(p,ix,iy);
				const Real x = p[0];
				const Real y = p[1];

				b(ix,iy).omega = 0.0;
				b(ix,iy).u[0] = 1.0;
				b(ix,iy).u[1] = 1.0;
				b(ix,iy).tmp = PoissonEquationPotIC::_fill(x,y);
			}
	}
}


void I2D_TestPoissonEquationPotential::_dump(string filename)
{
	IO_VTKNative<W,B, 4,0> vtkdumper;
	vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
}

Real I2D_TestPoissonEquationPotential::_mollified_heaviside(const double dist, const double eps) const
{
	//Positive outside/negative inside
	const double alpha = M_PI*min(1., max(0., (dist+0.5*eps)/eps));
	return 0.5+0.5*cos(alpha);
}

void I2D_TestPoissonEquationPotential::run()
{
	const int LMAX = parser("-lmax").asInt();
	const int BDP = parser("-bpd").asInt();

	_dump("before");

	poisson_solver->compute_velocity();

	_refine(false);
	_ic(*grid);

	poisson_solver->compute_velocity();

	divergence->perform();

	// Compute Linf
	Real Linf = 0.0;
	Real xerr = 0.0;
	Real yerr = 0.0;
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		B& b = grid->getBlockCollection()[vInfo[i].blockID];

		for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
			{
				Real p[2];
				info.pos(p,ix,iy);
				const Real x = p[0];
				const Real y = p[1];

				const Real error = fabs(PoissonEquationPotIC::_fill(x,y) - b(ix,iy).tmp);
				if(error>Linf)
				{
					Linf = error;
					xerr = p[0];
					yerr = p[1];
				}
				b(ix,iy).tmp = error;
			}
	}

	_dump("after");

	int resolution = 0;

	if( parser("-uniform").asBool() )
		resolution = BDP*32;
	else
		resolution = pow(2,LMAX)*32;

	printf("%d %e\n", resolution, Linf);

	printf("TEST CONCLUDED!\n");

	exit(0);
}


