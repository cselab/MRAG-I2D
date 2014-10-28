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

class I2D_LinkedBodies: public I2D_ObstacleOperator
{
	Real length, width, epsilon, position[2];
	vector<Real> angles;		
	
	Grid<W,B>& grid;
	BlockProcessing block_processing;
	
public:
	
	I2D_LinkedBodies(Grid<W,B>& grid, const Real length, const Real width, const Real position[2], vector<Real> anglesIn, const Real epsilon): 
	grid(grid), length(length), width(width), epsilon(epsilon)
	{
		this->position[0] = position[0];
		this->position[1] = position[1];
		
		for(unsigned int i=0; i<anglesIn.size(); i++)
			angles.push_back(anglesIn[i]);
	}
	
	void characteristic_function();
	
	static Real mollified_heaviside(const Real x, const Real eps);
	
	Real getD() const {return length;}
};

