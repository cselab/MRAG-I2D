/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

// TODO: CHECK THIS FILE!

#include "I2D_ImposedEllipseRotation.h"

I2D_ImposedEllipseRotation::I2D_ImposedEllipseRotation(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real _xcrot, const Real _ycrot, const Real _T, const Real _D, const Real _aspectRatio, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization):
I2D_ImposedEllipse(parser, grid, _xm, _ym, _xcrot, _ycrot, 0.0, 0.0, _D, _aspectRatio, eps, Uinf, penalization), XCROT(_xcrot), YCROT(_ycrot), T(_T), aspectRatio(_aspectRatio), D(_D)
{
}

I2D_ImposedEllipseRotation::~I2D_ImposedEllipseRotation()
{
}

Real I2D_ImposedEllipseRotation::getModulusMaxVel()
{
	assert(T>0.0);
	return fabs( 2.0*M_PI/T*D);
}

void I2D_ImposedEllipseRotation::_setMotionPattern(const Real t)
{
	assert(T>0.0);
	
	shape->vx = 0.0;
	shape->vy = 0.0;
	
	shape->angular_velocity = 2.0*M_PI/T;
}

void I2D_ImposedEllipseRotation::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	assert(T>0.0);

	shape->angle = 2*M_PI*t/T;

	shape->xm = XCROT + cos(shape->angle)*(D/2.0);
	shape->ym = YCROT + sin(shape->angle)*(D/2.0);

	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("update_I2D_ImposedEllipseRotation.txt", t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);		
	}
	else
	{
		// If string is set I open the corresponding file
		ppFile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);
	}

	// Write update data
	fprintf(ppFile, "%e %e %e %e %e %e %e %e %e %e\n", t, shape->xm, shape->ym, shape->vx, shape->vy, shape->angle, shape->angular_velocity, shape->J, shape->m, shape->rho);

	// Cloase file
	fclose(ppFile);
}








