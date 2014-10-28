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

class I2D_StefanFish: public I2D_CarlingFish
{
public:
	class StefanFish: public Fish
	{
	protected:
		double mControlPoints[6];

		double * TANX;
		double * TANY;
		double * VTANX;
		double * VTANY;

		double * DKDT;
		double * KT;
		double * K;
		double * mCurvatureDerivative;
		double mCurvatureStart[6];
		double mCurvatureEnd[6];
		double mCurrentCurvature[6];
		double mCurvatureTimeStart;
		double mCurvatureTimeEnd;

		double * mBaselineDerivative;
		double * mBaselineInterpolated;
		double mBaselineStart[6];
		double mBaselineEnd[6];
		double mCurrentBaseline[6];
		double mBaselineTimeStart;
		double mBaselineTimeEnd;

		double mTauStart;
		double mTauEnd;
		double mTauCurrent;
		double mTauTimeStart;
		double mTauTimeEnd;
		double mTauDerivative;

		double mTStart;
		double mTEnd;
		double mTCurrent;
		double mTTimeStart;
		double mTTimeEnd;
		double mTDerivative;

		double mTimeStartWave;
		bool mBoolWaveStarted;

		void _interpolateCurvature(const double * x, const double * y, const unsigned int n, const double * xx, double * yy, const unsigned int nn);
		void _applyTravellingWave(const double t, const double tau, const double tauDerivative, const double T, const double TDerivative, const double * s, const double * currentCurvature, const double * currentBaseline,
				const double * derivativeCurvature, const double * derivativeBaseline, double * xt, double * xt_dt, const unsigned int nn, bool boolWaveStarted);
		void _solveFrenetInHeadRefFrame(const double ds, const double * k, double * rX, double * rY, double * tanX, double * tanY, double * norX, double * norY, const unsigned int n);
		void _solveVelProfileInHeadRefFrame(const double ds, const double * kt, const double * kt_dt, const double * tanX, const double * tanY, const double * norX, const double * norY,
				double * vX, double * vY, double * vTanX, double * vTanY, double * vNorX, double * vNorY, const unsigned int n);
		void _schedulerCurvature(double t, double timeTransition, double * curvature, int kdim, double * s, int N, double * current_K_profile, double * current_dKdt_profile);
		void _schedulerBaseline(double t, double timeTransition, double * baseline, int kdim, double * s, int N, double * current_B_profile, double * current_dBdt_profile);
		void _schedulerTau(double t, double timeTransition, double tau, double & current_tau, double & current_tau_derivative);
		void _schedulerT(double t, double timeTransition, double _T, double & current_T, double & current_T_derivative);
		bool _checkCurvatureIsTheSame(double const * newControlK, double const * controlK, unsigned int n);
		void _startWave(double time);

	public:
		double S1, S2, S3, S4, tau;
		vector<double> BASELINE;
		vector<double> CURVATURE;
		double mTransitionTime;

		StefanFish(double xm, double ym, double _D, double _T, double phase, double tau, double angle_rad, vector<double> BASELINE, vector<double> CURVATURE, double angleInSpace_rad, double eps, const int LMAX, const bool isSharp=false);
		virtual ~StefanFish();
		void updateInSpace(const double t);
		virtual void save(FILE * f) const;
		virtual void restart(FILE * f);
		void clear();
	};

	I2D_StefanFish(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real T, const Real phase, const Real tau, const Real angle, vector<double> BASELINE, vector<double> CURVATURE,
			const Real angleInSpace, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID = 0, RL::RL_TabularPolicy ** policy = NULL, const int seed = 0);
	virtual ~I2D_StefanFish(){};
};
