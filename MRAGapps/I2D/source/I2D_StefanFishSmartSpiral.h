/*
 * I2D_StefanFishSmartSpiral.h
 *
 *  Created on: Mar 6, 2012
 *      Author: mgazzola
 */

#ifndef I2D_STEFANFISHSMARTSPIRAL_H_
#define I2D_STEFANFISHSMARTSPIRAL_H_

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_StefanFishSmart.h"

class I2D_StefanFishSmartSpiral: public I2D_StefanFishSmart
{
protected:
	double spiralT, spiralAngle, omega, A;
	double target[2];
public:
	I2D_StefanFishSmartSpiral(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real T, const Real tau, const Real angle, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID = 0, RL::RL_TabularPolicy ** _policy = NULL, const int seed = 0);
	virtual ~I2D_StefanFishSmartSpiral();

	virtual void mapAction(int action);
	virtual bool mapState(vector<int> & state, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL );
	virtual void reward(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	virtual void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
};

#endif /* I2D_STEFANFISHSMARTSPIRAL_H_ */
