/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#include "I2D_ImposedCylinderOscillation.h"

I2D_ImposedCylinderOscillation::I2D_ImposedCylinderOscillation(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real _A, const Real _T, const string _direction, const Real _D, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization):
I2D_ImposedCylinder(parser, grid, _xm, _ym, 0.0, 0.0, _D, eps, Uinf, penalization), XCM(_xm), YCM(_ym), A(_A), T(_T), direction(_direction)
{
}

I2D_ImposedCylinderOscillation::~I2D_ImposedCylinderOscillation()
{
}

Real I2D_ImposedCylinderOscillation::getModulusMaxVel()
{
	assert(T>0.0);
	return fabs( 2.0*A*M_PI/T );
}

void I2D_ImposedCylinderOscillation::_setMotionPattern(const Real t)
{
	assert(T>0.0);	
	assert(direction=="x" || direction=="y");
	
	if( direction == "x" )
	{
		shape->vx = A*cos(2.0*M_PI/T*t)*(2.0*M_PI/T);
		shape->vy = 0.0;
	}
	else if( direction == "y" )
	{
		shape->vx = 0.0;
		shape->vy = A*cos(2.0*M_PI/T*t)*(2.0*M_PI/T);
	}
	else
	{
		printf("_setMotionPattern\n");
		printf("Set cylinder direction!\n");
		abort();
	}
	
	shape->angular_velocity = 0.0;
}

void I2D_ImposedCylinderOscillation::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	assert(T>0.0);
	assert(direction=="x" || direction=="y");	
	
	if( direction == "x" )
	{
		shape->xm = XCM + A*sin(2.0*M_PI/T*t);
		shape->ym = YCM;
	}
	else if( direction == "y" )
	{
		shape->xm = XCM;
		shape->ym = YCM + A*sin(2.0*M_PI/T*t);
	}
	else
	{
		printf("update\n");
		printf("Set cylinder direction!\n");
		abort();
	}
	
	shape->angle = 0;
	
	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("update_I2D_ImposedCylinderOscillation.txt", t == 0.0 ? "w" : "a");
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








