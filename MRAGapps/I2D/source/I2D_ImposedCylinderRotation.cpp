/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#include "I2D_ImposedCylinderRotation.h"

I2D_ImposedCylinderRotation::I2D_ImposedCylinderRotation(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real _radiusRot, const Real _T, const Real _D, const Real _align, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization):
I2D_ImposedCylinder(parser, grid, _xm, _ym, 0.0, 0.0, _D, eps, Uinf, penalization), XCM(_xm), YCM(_ym), radiusRot(_radiusRot), T(_T), align(_align)
{
}

I2D_ImposedCylinderRotation::~I2D_ImposedCylinderRotation()
{
}

Real I2D_ImposedCylinderRotation::getModulusMaxVel()
{
	assert(T>0.0);
	return fabs( 2.0*M_PI*radiusRot/T );
}

void I2D_ImposedCylinderRotation::_setMotionPattern(const Real t)
{
	assert(T>0.0);	
	assert(align==1.0 || align==0.0);

	if( align == 1.0 )
	{
		shape->vx = -2*M_PI/T*sin(2*M_PI/T*t)*(radiusRot);
		shape->vy = 2*M_PI/T*cos(2*M_PI/T*t)*(radiusRot);
	}
	else if( align == 0.0 )
	{
		shape->vx = -2*M_PI/T*cos(2*M_PI/T*t)*(radiusRot);
		shape->vy = -2*M_PI/T*sin(2*M_PI/T*t)*(radiusRot);
	}
	else
	{
		printf("_setMotionPattern\n");
		printf("Set cylinder alignment!\n");
		abort();
	}
	
	shape->angular_velocity = 2.0*M_PI/T;
}

void I2D_ImposedCylinderRotation::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	assert(T>0.0);
	assert(align==1.0 || align==0.0);

	if( align == 1.0 )
	{
		shape->xm = XCM + cos(2*M_PI/T*t)*(radiusRot);
		shape->ym = YCM + sin(2*M_PI/T*t)*(radiusRot);
	}
	else if( align == 0.0 )
	{
		shape->xm = XCM - sin(2*M_PI/T*t)*(radiusRot);
		shape->ym = YCM + cos(2*M_PI/T*t)*(radiusRot);
	}
	else
	{
		printf("update\n");
		printf("Set cylinder alignment!\n");
		abort();
	}
	
	shape->angle = 2*M_PI*t;
	
	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("update_I2D_ImposedCylinderRotation.txt", t == 0.0 ? "w" : "a");
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








