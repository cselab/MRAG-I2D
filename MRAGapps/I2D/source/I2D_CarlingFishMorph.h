/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"
#include "I2D_FloatingObstacleOperator.h"
#include "I2D_CarlingFish.h"

class I2D_CarlingFishMorph: public I2D_CarlingFish
{
public:

	class CarlingFishMorph: public Fish
	{
	protected:
		vector<double> WIDTH;

		double _getWidth(const double & ss) const;
		double _coxDeBoorRecursion(const int j, const int d, const double * t, const double u) const;

	public:
		CarlingFishMorph(double xm, double ym, double _D, double _T, double phase, double angle_rad, vector<double> WIDTH, double angleInSpace_rad, double eps, const int LMAX);
		virtual ~CarlingFishMorph(){};
	};

	I2D_CarlingFishMorph(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real T, const Real phase, const Real angle, vector<double> WIDTH, const Real angleInSpace, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID = 0);
	virtual ~I2D_CarlingFishMorph(){};
};
