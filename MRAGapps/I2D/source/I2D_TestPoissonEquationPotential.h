/*
 *  I2D_TestPoissonEquationPotential.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Mattia Gazzola on 11/05/11.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_VelocityOperator.h"
#include "I2D_CurlVelocityOperator.h"
#include "I2D_DivOperator.h"

class I2D_TestPoissonEquationPotential: public I2D_Test
{
	ArgumentParser parser;
	
	Grid<W,B> * grid;
	Refiner * refiner;
	Compressor * compressor;
	
	BlockFWT<W, B, vorticity_projector, false, 1> fwt_omega;
	BlockFWT<W, B, velocity_projector, false, 2> fwt_velocity;
	BlockFWT<W, B, vorticityANDvelocity_projector, false, 3> fwt_omegaANDvelocity;
	BlockFWT<W, B, vorticityANDvelocityANDchi_projector, false, 4> fwt_omegaANDvelocityANDchi;

	
	set<int> _getBoundaryBlockIDs();
	void _dump(string filename);
	void _ic(Grid<W,B>& grid);
	Real _mollified_heaviside(const double dist, const double eps) const;
	
	void _refine(bool bUseIC=false);
	void _compress(bool bUseIC);
	
	I2D_VelocityOperator * poisson_solver;
	I2D_CurlVelocityOperator * curl_velocity;
	I2D_DivOperator * divergence;
	
public:
	
	I2D_TestPoissonEquationPotential(const int argc, const char ** argv);
	
	void run();
	void paint();
};

