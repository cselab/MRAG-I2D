/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#pragma once

#include "I2D_FloatingObstacleOperator.h"
#include "I2D_PenalizationOperator.h"
#include "RL_TabularPolicy.h"

class I2D_ObjectFactory
{	
	const int LMAX;
	const Real eps, charLength, XCM, YCM;
	Real Uinf[2];
	Real modulusMaxV;
	
	Grid<W,B>& grid;
	
	I2D_PenalizationOperator& penalization;
	
	int _lines(string filename);
	
public:
	I2D_ObjectFactory(Grid<W,B>& grid, const Real charLength, const Real XCM, const Real YCM, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int _LMAX);
	~I2D_ObjectFactory();	
	
	void create(ArgumentParser & parser, map< string, vector<I2D_FloatingObstacleOperator *> > & shapesMap, bool smart = false, RL::RL_TabularPolicy ** policy = NULL);
	
	Real getModulusMaxVel(){ return this->modulusMaxV; }
};
