/*
 *  I2D_TestDumping.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"


class I2D_TestDumping: public I2D_Test
{
	ArgumentParser parser;
	Grid<W,B> * grid;
	Refiner_BlackList refiner;
	Compressor compressor;
	BlockFWT<W, B, vorticity_projector, false, 1> fwt_omega;
	BlockFWT<W, B, velocity_projector, false, 1> fwt_psi;
	
	set<int> _getBoundaryBlockIDs();
	static void _ic(Grid<W,B>& grid);
		
public:
	
	I2D_TestDumping(const int argc, const char ** argv);
	
	void run();
	void paint();
};