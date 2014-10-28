/*
 *  I2D_PitchingAirfoil.h
 *  I2D_ROCKS
 *
 *  Created by Silvio Tödtli on 7/12/12.
 *  Copyright 2012 ETHZ. All rights reserved.
 *
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"
#include "I2D_FloatingObstacleOperator.h"
#include <map>


namespace PitchingStuff {
	class DiscretizedWing;
}


class I2D_PitchingAirfoil: public I2D_FloatingObstacleOperator
{
private:
	const Real plungePeriod, plungeAmplitude, plungePhase;
	const Real pitchPeriod, pitchAmplitude, pitchPhase;
	const Real D_char, l; // D_char: characteristic length, l: x-offset rotation center from leading edge
	const Real xcrot0, ycrot0; // time invariant part of rotation center

public:
	I2D_PitchingAirfoil(Grid<W,B>& grid, ArgumentParser& parser,
		const Real _xm, const Real _ym, const Real _D, const Real _l, const int d1,
		const int d2, const int d3d4, const Real _plungeAmplitude, const Real _plungePeriod,
		const Real _plungePhase, const Real _pitchAmplitude, const Real _pitchPeriod,
		const Real _pitchPhase, const Real _eps, const Real Uinf[2], I2D_PenalizationOperator& penalization);

	~I2D_PitchingAirfoil();

	void characteristic_function();

	Real getD() const {return 2./16;} // adopted from WingObstacleOperator
	Real getModulusMaxVel();
	
	void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	void computeDesiredVelocity(const double t);
	void create(const double t){} // remains like this

	void save(const double t, string filename = std::string());
	void restart(const double t, string filename = std::string());
	void refresh(const double t, string filename = std::string());

	PitchingStuff::DiscretizedWing * wing;
	
protected:
	virtual void _setMotionPattern(const Real t);
};
