/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#include "I2D_ImposedCylinderStopAndGo.h"

I2D_ImposedCylinderStopAndGo::I2D_ImposedCylinderStopAndGo(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real _vx, const Real _vy, const Real _tstop, const Real _D, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization):
I2D_ImposedCylinder(parser, grid, _xm, _ym, _vx, _vy, _D, eps, Uinf, penalization), tstop(_tstop)
{
}

I2D_ImposedCylinderStopAndGo::~I2D_ImposedCylinderStopAndGo()
{
}

Real I2D_ImposedCylinderStopAndGo::getModulusMaxVel()
{
	const Real modv = sqrt(vx_imposed*vx_imposed+vy_imposed*vy_imposed);
	assert(modv>0.0);
	return modv;
}

void I2D_ImposedCylinderStopAndGo::_setMotionPattern(const Real t)
{
	const Real modv = sqrt(vx_imposed*vx_imposed+vy_imposed*vy_imposed);
	const Real tfinal = tstop-0.5*eps/modv;

	if(t<tfinal)
	{
		shape->vx = vx_imposed;
		shape->vy = vy_imposed;
	}
	else
	{
		FILE * ppFile = fopen("collisionTime.txt","a");
		assert(ppFile!=NULL);
		fprintf(ppFile,"tstop=%e, tfinal=%e, diff=%e, dt=%e\n",t,tfinal,tfinal-tstop,eps/fabs(vx_imposed));
		fclose(ppFile);

		shape->vx = 0.0;
		shape->vy = 0.0;
	}

	shape->angular_velocity = 0.0;
}

void I2D_ImposedCylinderStopAndGo::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	shape->update_all(dt,t);

	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("update_I2D_ImposedCylinderStopAndGo.txt", t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		ppFile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);
	}

	// Write update data
	fprintf(ppFile, "%e %e %e %e %e %e %e %e %e %e %e %e\n", this->dimT, t, shape->xm, shape->ym, shape->vx, shape->vy, shape->angle, shape->angular_velocity, shape->J, shape->m, shape->rho, this->Cd);

	// Cloase file
	fclose(ppFile);
}








