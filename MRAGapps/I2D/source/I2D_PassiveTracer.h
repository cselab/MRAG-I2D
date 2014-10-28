/**
 * @file I2D_PassiveTracer.cpp
 * @author Mattia Gazzola
 * @date Mar 1, 2012
 */

#ifndef I2D_PASSIVETRACER_H_
#define I2D_PASSIVETRACER_H_

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_FloatingObstacleOperator.h"
#include <vector>
#include <map>

class I2D_PassiveTracer: public I2D_FloatingObstacleOperator
{
protected:
	map<long int, vector<Real> > particles;
	map<I3, vector<long int> > block2particles;

public:
	I2D_PassiveTracer(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization);
	~I2D_PassiveTracer();

	void characteristic_function(){};
	Real getD() const {return 0.0;}

	void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	void computeDesiredVelocity(const double t){};
	void create(const double t){}

	void save(const double t, string filename = std::string());
	void restart(const double t, string filename = std::string());
	void refresh(const double t, string filename = std::string());
};

#endif /* I2D_PASSIVETRACER_H_ */
