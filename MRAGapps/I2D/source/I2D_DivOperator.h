/*
 *  I2D_DivOperator.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 * computes the divergence of the velocity or vorticity field into tmp
 * with a fourth order FD scheme. Tested!
 *
 * velocity():
 * IN: u[0-2]
 * OUT: tmp
 *
 * vorticity():
 * IN: omega[0-2]
 * OUT: tmp
 * 
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"

class I2D_DivOperator
{
	Grid<W,B>& grid;
	BlockProcessing block_processing;
	
public:
	
	I2D_DivOperator(Grid<W,B>& grid): grid(grid){}
	
	void perform();
};
