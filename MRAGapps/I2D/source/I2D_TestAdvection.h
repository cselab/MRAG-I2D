/*
 *  I2D_TestAdvection2D.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_AdvectionOperator.h"


class I2D_TestAdvection: public I2D_Test
{
	ArgumentParser parser;
	
	Grid<W,B> * grid;
	Refiner * refiner;
	Compressor * compressor;
	
	BlockFWT<W, B, vorticity_projector, false, 1> fwt_omega;
	
	I2D_AdvectionOperator * advection;
	
	set<int> _getBoundaryBlockIDs();
	void _dump(string filename);
	static void _ic_omega(Grid<W,B>& grid);
	static void _ic_velocity(Grid<W,B>& grid);
	
	int step_id;
	
	void _refine(bool bUseIC);
	void _compress(bool bUseIC);
	
public:
	
	I2D_TestAdvection(const int argc, const char ** argv);
	
	void run();
	void paint();
};

