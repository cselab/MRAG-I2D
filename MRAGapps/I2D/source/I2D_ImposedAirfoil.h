/*
 *  I2D_ImposedAirfoil.h
 *  I2D_ROCKS
 *
 *  Created by Silvio Tödtli on 7/3/12.
 *  Copyright 2012 ETHZ. All rights reserved.
 *
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"
#include "I2D_FloatingObstacleOperator.h"
#include <map>


namespace AirfoilStuff {
	class DiscretizedWing;
}


class I2D_ImposedAirfoil: public I2D_FloatingObstacleOperator
{
public:
	I2D_ImposedAirfoil(Grid<W,B>& grid, ArgumentParser& parser, const Real _xm, const Real _ym, const Real _D, const Real _angle,
			const Real _vx, const Real _vy, const int d1, const int d2, const int d3d4, const Real _eps, const Real Uinf[2],
			I2D_PenalizationOperator& penalization);

	~I2D_ImposedAirfoil();


	void characteristic_function();

	Real getD() const {return 2./16;} // adopted from WingObstacleOperator
	
	void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	void computeDesiredVelocity(const double t);
	void create(const double t){} // remains like this

	void save(const double t, string filename = std::string());
	void restart(const double t, string filename = std::string());
	void refresh(const double t, string filename = std::string());

	AirfoilStuff::DiscretizedWing * wing;
	
protected:
	Real vx_imposed, vy_imposed;
	virtual void _setMotionPattern(const Real t);
};
