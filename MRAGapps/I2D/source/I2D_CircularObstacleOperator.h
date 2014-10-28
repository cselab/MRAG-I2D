/*
 *  I2D_SphereObstacleOperator.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 *	compute the characteristic function of an obstacle and
 *	store it into tmp. Tested!
 *
 *	characteristic_function():
 *	IN: 
 *	OUT: tmp
 *
 */
#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"

class I2D_CircularObstacleOperator: public I2D_ObstacleOperator
{
	Real radius, position[2];
	Real smoothing_length, smooth_radius;
	bool SHARP;

	Grid<W,B>& grid;
	BlockProcessing block_processing;
	
	Real _get_smooth_radius(const Real radius, const Real smoothing_length);

public:
	
 I2D_CircularObstacleOperator(Grid<W,B>& grid, const Real radius, const Real position[2], const Real smoothing_length, const bool isSharp=false): 
	grid(grid), radius(radius), SHARP(isSharp)
	{
		this->position[0] = position[0];
		this->position[1] = position[1];
		this->smoothing_length = smoothing_length;

		assert(radius-smoothing_length > 0);
		smooth_radius = _get_smooth_radius(radius, smoothing_length);
	}
	
	void characteristic_function();
	
	static Real mollified_heaviside(const Real x, const Real eps);
	
	Real getD() const {return 2*radius;}
};
