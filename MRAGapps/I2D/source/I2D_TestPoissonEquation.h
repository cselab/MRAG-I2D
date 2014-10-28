/*
 *  I2D_TestPoissonEquation.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_VelocityOperator.h"
#include "I2D_CurlVelocityOperator.h"
#include "I2D_DivOperator.h"
//#include "I2D_ScalarLaplaceOperator.h"

class I2D_TestPoissonEquation: public I2D_Test
{
	ArgumentParser parser;
	
	Grid<W,B> * grid;
	Refiner * refiner;
	Compressor * compressor;
	
	BlockFWT<W, B, vorticity_projector, false, 1> fwt_omega;
	BlockFWT<W, B, velocity_projector, false, 2> fwt_velocity;
	
	set<int> _getBoundaryBlockIDs();
	void _dump(string filename);
	static void _ic_omega(Grid<W,B>& grid);
	void _load_ic(string filename);
	void _refine_omega(bool bUseIC);
	void _refine_vel(bool bUseIC=false);
	void _compress(bool bUseIC);
	
	I2D_VelocityOperator * poisson_solver;
	I2D_CurlVelocityOperator * curl_velocity;
	I2D_DivOperator * divergence;
//	I2D_ScalarLaplaceOperator * laplace_psi;
	
public:
	
	I2D_TestPoissonEquation(const int argc, const char ** argv);
	
	void run();
	void paint();
};

