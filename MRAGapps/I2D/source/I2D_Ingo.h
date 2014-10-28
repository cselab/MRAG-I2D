/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 *  This class is the same as the I2D_LinkedBodies class. It just adds (if possible) a last plate for which the vertical position 
 *  of the second extremity is constrained to be alligned with the extremity of the leading edge (cf Ingo experiments).
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"

class I2D_Ingo: public I2D_ObstacleOperator
{
	Real length, width, epsilon, position[2];
	vector<Real> angles;		
	
	Grid<W,B>& grid;
	BlockProcessing block_processing;
	
public:
	
	I2D_Ingo(Grid<W,B>& grid, const Real length, const Real width, const Real position[2], vector<Real> anglesIn, const Real epsilon): 
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

