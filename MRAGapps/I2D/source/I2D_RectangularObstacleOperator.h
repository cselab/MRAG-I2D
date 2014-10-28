/*
 *  I2D_RectangularObstacleOperator.h
 *  I2D_ROCKS
 *
 *  Created by Diego Rossinelli on 2/2/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"

class I2D_RectangularObstacleOperator: public I2D_ObstacleOperator
{
	Real aspect_ratioX, position[2];
	Real smoothing_length, eps0, eps1, D;
	
	Grid<W,B>& grid;
	BlockProcessing block_processing;
	
	Real _get_eps1();
	
public:
	
	I2D_RectangularObstacleOperator(Grid<W,B>& grid, const Real aspect_ratioX, const Real position[2], const Real D, const Real smoothing_length): 
	grid(grid), aspect_ratioX(aspect_ratioX), D(D), smoothing_length(smoothing_length)
	{
		this->position[0] = position[0];
		this->position[1] = position[1];
		
		eps0 = smoothing_length;
		assert(D*aspect_ratioX-2*eps0 >= 0);
		assert(D - 2*eps0 >= 0);
		
		eps1 = _get_eps1();
		assert(eps1>=eps0);
	}
	
	void characteristic_function();
	
	static Real mollified_heaviside(const Real x, const Real eps);
	
	Real getD() const {return D;}
};

