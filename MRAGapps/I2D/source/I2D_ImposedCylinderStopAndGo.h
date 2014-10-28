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
#include "I2D_ImposedCylinder.h"

class I2D_ImposedCylinderStopAndGo: public I2D_ImposedCylinder
{	
	const Real tstop;

public:

	I2D_ImposedCylinderStopAndGo(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real _vx, const Real _vy, const Real _tstop, const Real _D, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization);
	~I2D_ImposedCylinderStopAndGo();
	
	Real getModulusMaxVel();
	void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	
protected:
	void _setMotionPattern(const Real t);
};
