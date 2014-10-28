/**
 */
#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"
#include "I2D_ImposedCylinder.h"
#include "I2D_ImposedEllipse.h"

class I2D_ImposedEllipseOscillation : public I2D_ImposedEllipse
{
private:
	const Real plungePeriod;
	const Real plungeAmplitude;
	const Real plungePhase;
	const Real pitchPeriod;
	const Real pitchAmplitude;
	const Real pitchPhase;
	const Real x0; /// initial position
	const Real y0; /// initial position
	const Real D; /// characteristic length


public:
	I2D_ImposedEllipseOscillation(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm,
			const Real _ym, const Real _xcrot, const Real _ycrot,
			const Real _aspectRatio, const Real _plungeAmplitude,
			const Real _plungePeriod, const Real _pitchAmplitude,
			const Real _pitchPeriod, const Real _plungePhase,
			const Real _pitchPhase, const Real _D, const Real eps,
			const Real Uinf[2], I2D_PenalizationOperator& penalization);
	~I2D_ImposedEllipseOscillation();

	Real getModulusMaxVel();
	void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);

protected:
	void _setMotionPattern(const Real t);
};
