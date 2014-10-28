/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

// TODO: CHECK THIS FILE!

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"
#include "I2D_ImposedEllipse.h"

class I2D_ImposedEllipseRotation: public I2D_ImposedEllipse
{	
	const Real T, XCROT, YCROT, aspectRatio, D;
	
public:

	I2D_ImposedEllipseRotation(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real _xcrot, const Real _ycrot, const Real _T, const Real _D, const Real _aspectRatio, const Real eps,
									const Real Uinf[2], I2D_PenalizationOperator& penalization);
	~I2D_ImposedEllipseRotation();
	
	Real getModulusMaxVel();
	void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);

protected:
	void _setMotionPattern(const Real t);
};
