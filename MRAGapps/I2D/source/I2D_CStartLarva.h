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

class I2D_CStartLarva: public I2D_StefanFish
{
public:

	class CStartLarva: public StefanFish
	{
	protected:
		double _getWidth(const double & ss) const;
	public:
		double Tprep, Tprop;

		CStartLarva(double xm, double ym, double _D, double Tprep, double Tprop, double phase, double tau, double angle_rad, vector<double> BASELINE, vector<double> CURVATURE, double angleInSpace_rad, double eps, const int LMAX, const bool isSharp=false);
		virtual ~CStartLarva(){};
		void updateInSpace(const double t);
	};

	I2D_CStartLarva(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real Tprep, const Real Tprop, const Real phase, const Real tau, const Real angle, vector<double> BASELINE, vector<double> CURVATURE,
			const Real angleInSpace, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID = 0, RL::RL_TabularPolicy ** policy = NULL);
	virtual ~I2D_CStartLarva(){};
};
