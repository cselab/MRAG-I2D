/*
 * I2D_StefanFishSmartLattice.h
 *
 *  Created on: Feb 29, 2012
 *      Author: mgazzola
 */

#ifndef I2D_STEFANFISHSMARTLATTICE_H_
#define I2D_STEFANFISHSMARTLATTICE_H_

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_StefanFishSmart.h"

class I2D_StefanFishSmartLattice : public I2D_StefanFishSmart
{
protected:
	double scale;
	double cutoff;
	const double T;
	double dir[2];
	double target[2];
	double followedPoint[2];

public:
	I2D_StefanFishSmartLattice(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real _T, const Real tau, const Real angle, const Real _dir[2], const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID = 0, RL::RL_TabularPolicy ** _policy = NULL, const int seed = 0);
	virtual ~I2D_StefanFishSmartLattice();

	virtual void mapAction(int action);
	virtual bool mapState(vector<int> & state, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL );
	virtual void reward(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	virtual void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
};

#endif /* I2D_STEFANFISHSMARTLATTICE_H_ */
