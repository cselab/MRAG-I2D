/*
 * I2D_StefanFishSmart.h
 *
 *  Created on: Feb 8, 2012
 *      Author: mgazzola
 */

#ifndef I2D_STEFANFISHSMART_H_
#define I2D_STEFANFISHSMART_H_

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_StefanFish.h"

class I2D_StefanFishSmart : public I2D_StefanFish
{
protected:
	vector<int> signature;
	double maxDomainRadius;

	bool _mapStateDistance(vector<int> & state, map< string, vector<I2D_FloatingObstacleOperator *> > * _data);

public:
	class StefanFishSmart: public StefanFish
	{
	public:
		StefanFishSmart(double xm, double ym, double _D, double T, double tau, double angle_rad, double eps, const int LMAX);
		virtual ~StefanFishSmart(){};
	};

	I2D_StefanFishSmart(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real T, const Real tau, const Real angle, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID = 0, RL::RL_TabularPolicy ** _policy = NULL, const int seed = 0);
	virtual ~I2D_StefanFishSmart();

	virtual void mapAction(int action);
	virtual bool mapState(vector<int> & state, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL );
	virtual void reward(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
};

#endif /* I2D_STEFANFISHSMART_H_ */
