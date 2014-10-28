/*
 *  I2D_RotatingWheelObstacle.h
 *  I2D_ROCKS
 *
 *  Created by Diego Rossinelli on 1/11/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_FloatingObstacleOperator.h"
#include "I2D_PenalizationOperator.h"

class I2D_RotatingWheelObstacle: public I2D_FloatingObstacleOperator
{
	Real smoothing_length, radius;
	Real Uinf[2];
	
	Grid<W,B>& grid;
	BlockProcessing block_processing;
	
	map<int , const VelocityBlock *> desired_velocity;
	I2D_PenalizationOperator& penalization;
		
public:
	
	I2D_RotatingWheelObstacle(ArgumentParser & parser, Grid<W,B>& grid, const Real smoothing_length, const Real radius, const Real Uinf[2], I2D_PenalizationOperator& penalization);
	
	void characteristic_function();		
	Real getD() const {return 2*radius;}
	
	void update(const double dt, const double t, string filename =std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	void computeDesiredVelocity(const double t);
	void create(const double t){}
	
	void save(const double t, string filename = std::string());
	void restart(const double t, string filename = std::string());
	
	void refresh(const double t, string filename = std::string()) {} //i dont need this
};
