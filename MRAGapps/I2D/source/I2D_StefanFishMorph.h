/*
 * I2D_CStartLarva.h
 *
 *  Created on: Dec 8, 2011
 *      Author: mgazzola
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"
#include "I2D_FloatingObstacleOperator.h"
#include "I2D_StefanFish.h"

class I2D_StefanFishMorph: public I2D_StefanFish
{
public:

	class StefanFishMorph: public StefanFish
	{
	protected:
		vector<double> WIDTH;

		double _getWidth(const double & ss) const;
		double _coxDeBoorRecursion(const int j, const int d, const double * t, const double u) const;

	public:
		StefanFishMorph(double xm, double ym, double _D, double _T, double phase, double tau, double angle_rad, vector<double> WIDTH, vector<double> BASELINE, vector<double> CURVATURE, double angleInSpace_rad, double eps, const int LMAX, const bool isSharp=false);
		virtual ~StefanFishMorph(){};
	};

	I2D_StefanFishMorph(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real T, const Real phase, const Real tau, const Real angle, vector<double> WIDTH, vector<double> BASELINE, vector<double> CURVATURE, const Real angleInSpace, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID = 0);
	virtual ~I2D_StefanFishMorph(){};
};
