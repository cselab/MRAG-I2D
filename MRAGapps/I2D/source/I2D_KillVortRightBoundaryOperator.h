/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"

class I2D_KillVortRightBoundaryOperator
{
	int killed_width;
	
	Grid<W,B>& grid;
	BlockProcessing block_processing;
	
	
public:
	
	I2D_KillVortRightBoundaryOperator(Grid<W,B>& grid, Real killed_width): 
	grid(grid), killed_width(killed_width)
	{
	}
		
	void killVorticity();
};



