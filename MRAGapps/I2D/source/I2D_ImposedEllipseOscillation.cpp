/**
 * Allows you to specify the ellipse plunge and pitch oscillations.
 *
 * @file I2D_ImposedEllipseOscillation.h
 * @author Andrew A. Tchieu
 * @date 21 February 2012
 */

#include "I2D_ImposedEllipseOscillation.h"

/// @todo Need to clean this object factory up, this is getting out of hand...
I2D_ImposedEllipseOscillation::I2D_ImposedEllipseOscillation(ArgumentParser & parser, Grid<W,B>& grid,
const Real _xm, const Real _ym, const Real _xcrot, const Real _ycrot,
		const Real _aspectRatio, const Real _plungeAmplitude,
		const Real _plungePeriod, const Real _pitchAmplitude,
		const Real _pitchPeriod, const Real _plungePhase,
		const Real _pitchPhase, const Real _D, const Real eps,
		const Real Uinf[2], I2D_PenalizationOperator& penalization) :
			I2D_ImposedEllipse(parser, grid, _xm, _ym, _xcrot, _ycrot, 0.0, 0.0, _D, _aspectRatio, eps, Uinf, penalization), // Inherited constructor
			x0(_xm), // define constants
			y0(_ym),
			plungeAmplitude(_plungeAmplitude),
			plungePeriod(_plungePeriod),
			pitchAmplitude(_pitchAmplitude),
			pitchPeriod(_pitchPeriod),
			plungePhase(_plungePhase),
			pitchPhase(_pitchPhase),
			D(_D)
{
	// Initialize
	shape->xm = x0;
	shape->ym = y0 + plungeAmplitude*sin(plungePhase);
	shape->xcrot = shape->xm;
	shape->ycrot = shape->ym;
	shape->angle = pitchAmplitude*sin(pitchPhase);
}

/**
 * Destructor.
 */
I2D_ImposedEllipseOscillation::~I2D_ImposedEllipseOscillation()
{
}

/**
 * Return the maximum velocity from given parameters.
 *
 * @return
 */
Real I2D_ImposedEllipseOscillation::getModulusMaxVel()
{
	Real translation;
	Real rotation;

	translation = fabs(2.0*plungeAmplitude*M_PI/plungePeriod);
	rotation = fabs(pitchAmplitude*M_PI/pitchPeriod*D);
	return (translation > rotation ? translation : rotation);
}

/**
 * Set the function of motion versus time.
 *
 * @param t
 */
void I2D_ImposedEllipseOscillation::_setMotionPattern(const Real t)
{
	Real factor1 = 2.0*M_PI/plungePeriod;
	Real factor2 = 2.0*M_PI/pitchPeriod;
	shape->vx = 0.0;
	shape->vy = plungeAmplitude*factor1*cos(factor1*t+plungePhase);
	shape->angular_velocity = pitchAmplitude*factor2*cos(factor2*t+pitchPhase);
}

/**
 * Update the location of the body explicitly.
 *
 * @param dt timestep
 * @param t time
 * @param filename
 */
void I2D_ImposedEllipseOscillation::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	shape->xm = x0;
	shape->ym = y0 + plungeAmplitude*sin(2.0*M_PI/plungePeriod*t+plungePhase);
	shape->xcrot = shape->xm;
	shape->ycrot = shape->ym;
	shape->angle = pitchAmplitude*sin(2.0*M_PI/pitchPeriod*t+pitchPhase);

	/// Write shape data to file
	FILE * ppFile = NULL;
	if(filename == std::string()) /// String not set, create my file
	{
		ppFile = fopen("update_I2D_ImposedEllipseOscillation.txt",
				t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);
	}
	else /// String set, so write info there
	{
		ppFile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);
	}
	fprintf(ppFile, "%e %e %e %e %e %e %e %e %e %e\n", t, shape->xm, shape->ym,
			shape->vx, shape->vy, shape->angle, shape->angular_velocity, shape->J,
			shape->m, shape->rho);
	fclose(ppFile);
}






