/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"

class I2D_EllipticalObstacleOperator: public I2D_ObstacleOperator
{
	Real semiMajorAxis, aspectRatio, angle, epsilon, position[2];
	
	Grid<W,B>& grid;
	BlockProcessing block_processing;

public:
	
	I2D_EllipticalObstacleOperator(Grid<W,B>& grid, const Real semiMajorAxis, const Real aspectRatio, const Real angle, const Real position[2], const Real epsilon): 
	grid(grid), semiMajorAxis(semiMajorAxis), aspectRatio(aspectRatio), angle(angle), epsilon(epsilon)
	{
		assert( (aspectRatio>0.0) && (aspectRatio<=1.0) );
		assert( semiMajorAxis>0.0 );
		
		this->angle = this->angle/180.0*M_PI;
		
		this->position[0] = position[0];
		this->position[1] = position[1];
	}
	
	void characteristic_function();
	
	static Real mollified_heaviside(const Real x, const Real eps);
	
	Real getD() const {return 2*semiMajorAxis;}
};



