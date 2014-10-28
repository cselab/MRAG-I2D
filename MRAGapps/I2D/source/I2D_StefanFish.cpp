/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include "I2D_StefanFish.h"
#include "I2D_Interpolation1D.h"
#include "I2D_Frenet2D.h"
#include "I2D_VectorBlockLab.h"
#include "I2D_GradOfVector.h"
#include <limits>
#include <omp.h>

void I2D_StefanFish::StefanFish::save(FILE * f) const
{
	fprintf(f, "xm: %20.20e\n", xm);
	fprintf(f, "ym: %20.20e\n", ym);
	fprintf(f, "angle: %20.20e\n", angle);
	fprintf(f, "angleInSpace: %20.20e\n", angleInSpace);

	// Curvature transition
	for(unsigned int i=0; i<6; i++)
		fprintf(f, "mCurvatureStart[x]: %20.20e\n", mCurvatureStart[i]);

	for(unsigned int i=0; i<6; i++)
		fprintf(f, "mCurvatureEnd[x]: %20.20e\n", mCurvatureEnd[i]);

	for(unsigned int i=0; i<6; i++)
		fprintf(f, "mCurrentCurvature[x]: %20.20e\n", mCurrentCurvature[i]);

	fprintf(f, "mCurvatureTimeStart: %20.20e\n", mCurvatureTimeStart);
	fprintf(f, "mCurvatureTimeEnd: %20.20e\n", mCurvatureTimeEnd);

	// Baseline transition
	for(unsigned int i=0; i<6; i++)
		fprintf(f, "mBaselineStart[x]: %20.20e\n", mBaselineStart[i]);

	for(unsigned int i=0; i<6; i++)
		fprintf(f, "mBaselineEnd[x]: %20.20e\n", mBaselineEnd[i]);

	for(unsigned int i=0; i<6; i++)
		fprintf(f, "mCurrentBaseline[x]: %20.20e\n", mCurrentBaseline[i]);

	fprintf(f, "mBaselineTimeStart: %20.20e\n", mBaselineTimeStart);
	fprintf(f, "mBaselineTimeEnd: %20.20e\n", mBaselineTimeEnd);

	// Tau transition
	fprintf(f, "mTauTimeStart: %20.20e\n", mTauTimeStart);
	fprintf(f, "mTauTimeEnd: %20.20e\n", mTauTimeEnd);
	fprintf(f, "mTauStart: %20.20e\n", mTauStart);
	fprintf(f, "mTauEnd: %20.20e\n", mTauEnd);
	fprintf(f, "mTauCurrent: %20.20e\n", mTauCurrent);

	// T transition
	fprintf(f, "mTTimeStart: %20.20e\n", mTTimeStart);
	fprintf(f, "mTTimeEnd: %20.20e\n", mTTimeEnd);
	fprintf(f, "mTStart: %20.20e\n", mTStart);
	fprintf(f, "mTEnd: %20.20e\n", mTEnd);
	fprintf(f, "mTCurrent: %20.20e\n", mTCurrent);

	// Wave transition
	fprintf(f, "mTimeStartWave: %20.20e\n", mTimeStartWave);
	fprintf(f, "mBoolWaveStarted: %d\n", mBoolWaveStarted);
}

void I2D_StefanFish::StefanFish::restart(FILE * f)
{
	float val = 0.0;
	int valInt = 0;

	fscanf(f, "xm: %e\n", &val);
	xm = val;
	printf("CarlingFish::restart(): xm is %e\n", xm);

	fscanf(f, "ym: %e\n", &val);
	ym = val;
	printf("CarlingFish::restart(): ym is %e\n", ym);

	fscanf(f, "angle: %e\n", &val);
	angle = val;
	printf("CarlingFish::restart(): angle is %e\n", angle);

	fscanf(f, "angleInSpace: %e\n", &val);
	angleInSpace = val;
	printf("CarlingFish::restart(): angleInSpace is %e\n", angleInSpace);

	// REMEMBER TO DO THE MAPPING OTHERWISE YOU ARE FUCKED!
	_set_ow_mapping();

	// Curvature transition
	for(unsigned int i=0; i<6; i++)
	{
		fscanf(f, "mCurvatureStart[x]: %e\n", &val);
		mCurvatureStart[i] = val;
	}

	for(unsigned int i=0; i<6; i++)
	{
		fscanf(f, "mCurvatureEnd[x]: %e\n", &val);
		mCurvatureEnd[i] = val;
	}

	for(unsigned int i=0; i<6; i++)
	{
		fscanf(f, "mCurrentCurvature[x]: %e\n", &val);
		mCurrentCurvature[i] = val;
	}

	fscanf(f, "mCurvatureTimeStart: %e\n", &val);
	mCurvatureTimeStart = val;
	fscanf(f, "mCurvatureTimeEnd: %e\n", &val);
	mCurvatureTimeEnd = val;

	// Baseline transition
	for(unsigned int i=0; i<6; i++)
	{
		fscanf(f, "mBaselineStart[x]: %e\n", &val);
		mBaselineStart[i] = val;
	}

	for(unsigned int i=0; i<6; i++)
	{
		fscanf(f, "mBaselineEnd[x]: %e\n", &val);
		mBaselineEnd[i] = val;
	}

	for(unsigned int i=0; i<6; i++)
	{
		fscanf(f, "mCurrentBaseline[x]: %e\n", &val);
		mCurrentBaseline[i] = val;
	}

	fscanf(f, "mBaselineTimeStart: %e\n", &val);
	mBaselineTimeStart = val;
	fscanf(f, "mBaselineTimeEnd: %e\n", &val);
	mBaselineTimeEnd = val;

	// Tau transition
	fscanf(f, "mTauTimeStart: %e\n", &val);
	mTauTimeStart = val;
	fscanf(f, "mTauTimeEnd: %e\n", &val);
	mTauTimeEnd = val;
	fscanf(f, "mTauStart: %e\n", &val);
	mTauStart = val;
	fscanf(f, "mTauEnd: %e\n", &val);
	mTauEnd = val;
	fscanf(f, "mTauCurrent: %e\n", &val);
	mTauCurrent = val;

	// T transition
	fscanf(f, "mTTimeStart: %e\n", &val);
	mTTimeStart = val;
	fscanf(f, "mTTimeEnd: %e\n", &val);
	mTTimeEnd = val;
	fscanf(f, "mTStart: %e\n", &val);
	mTStart = val;
	fscanf(f, "mTEnd: %e\n", &val);
	mTEnd = val;
	fscanf(f, "mTCurrent: %e\n", &val);
	mTCurrent = val;

	mTransitionTime = mTEnd;

	// Wave transition
	fscanf(f, "mTimeStartWave: %e\n", &val);
	mTimeStartWave = val;
	fscanf(f, "mBoolWaveStarted: %d\n", &valInt);
	mBoolWaveStarted = valInt;
}

bool I2D_StefanFish::StefanFish::_checkCurvatureIsTheSame(double const * newControlK, double const * controlK, unsigned int n)
{
	double sum = 0.0;
	for(unsigned int i=0;i<n;i++)
	{
		sum += (controlK[i]-newControlK[i]);
	}
	if(sum==0.0){ return true; }
	else{ return false; }
}

void I2D_StefanFish::StefanFish::_applyTravellingWave(const double t, const double tau, const double tauDerivative, const double T, const double TDerivative, const double * s, const double * currentCurvature, const double * currentBaseline,
		const double * derivativeCurvature, const double * derivativeBaseline, double * xt, double * xt_dt, const unsigned int nn, bool boolWaveStarted)
{
	const double offset = 0.0;//-(EXTENSION/2.0)*EPS;
	const double ll = 1.0;
	if( boolWaveStarted )
	{
		for(unsigned int i=0; i<nn; i++)
		{
			const double argument = 2.0*M_PI*( t/T - (s[i]+offset)/ll*tau ) + phase;
			const double argumentDerived = 2.0*M_PI*( (T-t*TDerivative)/(T*T) - (s[i]+offset)/ll*mTauDerivative );
			//const double argumentDerived = 2.0*M_PI*( 1.0/T - (s[i]+offset)/ll*mTauDerivative );
			xt[i] = currentBaseline[i] + currentCurvature[i]*sin(argument);
			xt_dt[i] = derivativeBaseline[i] + derivativeCurvature[i]*sin(argument) + currentCurvature[i]*cos(argument)*argumentDerived;
		}
	}
	else
	{
		for(unsigned int i=0; i<nn; i++)
		{
			xt[i] = currentBaseline[i] + currentCurvature[i];
			xt_dt[i] = derivativeBaseline[i] + derivativeCurvature[i];
		}
	}
}

void I2D_StefanFish::StefanFish::_interpolateCurvature(const double * x, const double * y, const unsigned int n, const double * xx, double * yy, const unsigned int nn)
{
	I2D_Interpolation1D iterpolator;
	const double offset = 0.0;//-(EXTENSION/2.0)*EPS;
	iterpolator.naturalCubicSpline( x, y, n, xx, yy, nn, offset );
}

void I2D_StefanFish::StefanFish::_schedulerCurvature(double t, double timeTransition, double * curvature, int kdim, double * s, int N, double * current_K_profile, double * current_dKdt_profile)
{
	// Curvature checks
	bool curvaturesStartEndEqual = _checkCurvatureIsTheSame(mCurvatureStart, mCurvatureEnd, kdim);
	bool curvaturesEndCurrentEqual = _checkCurvatureIsTheSame(mCurvatureEnd,curvature,kdim);

	// KEEP GOING
	if( curvaturesStartEndEqual && curvaturesEndCurrentEqual )
	{
		// Set time start/end
		mCurvatureTimeStart = 0.0;
		mCurvatureTimeEnd = 0.0;

		// Compute current curvature
		for(int i=0;i<kdim;i++)
		{
			mCurvatureStart[i] = mCurvatureEnd[i];
			mCurrentCurvature[i] =  mCurvatureEnd[i];
		}

		// Compute curvature & curvature derivative profile
		_interpolateCurvature(mControlPoints, mCurrentCurvature, kdim, s, current_K_profile, N);

		for(int i=0;i<N;i++)
			current_dKdt_profile[i] = 0.0;

		std::cout << "KEEP GOING (CURVATURE)" << std::endl;

		return;
	}


	// START TRANSITION FROM STATIONARY
	if( curvaturesStartEndEqual && !curvaturesEndCurrentEqual )
	{
		// Set time start/end
		mCurvatureTimeStart = t;
		mCurvatureTimeEnd = mCurvatureTimeStart + timeTransition;

		// Compute current curvature
		for(int i=0;i<kdim;i++)
		{
			mCurvatureEnd[i] = curvature[i];
			double dummy = 0.0;
			_cubicInterpolation(mCurvatureTimeStart, mCurvatureTimeEnd, t, mCurvatureStart[i], mCurvatureEnd[i], mCurrentCurvature[i], dummy);
		}

		// Compute curvature & curvature derivative profile
		double curvStartLong[N];
		double curvEndLong[N];
		_interpolateCurvature(mControlPoints, mCurvatureStart, kdim, s, curvStartLong, N);
		_interpolateCurvature(mControlPoints, mCurvatureEnd, kdim, s, curvEndLong, N);
		for(int i=0;i<N;i++)
			_cubicInterpolation(mCurvatureTimeStart, mCurvatureTimeEnd, t, curvStartLong[i], curvEndLong[i], current_K_profile[i], current_dKdt_profile[i]);

		std::cout << "START TRANSITION FROM STATIONARY (CURVATURE)" << std::endl;

		return;
	}


	// DURING TRANSITION & END TRANSITION
	if( !curvaturesStartEndEqual && curvaturesEndCurrentEqual )
	{
		if(t<=mCurvatureTimeEnd)
		{
			// Compute current curvature
			for(int i=0;i<kdim;i++)
			{
				double dummy = 0.0;
				_cubicInterpolation(mCurvatureTimeStart, mCurvatureTimeEnd, t, mCurvatureStart[i], mCurvatureEnd[i], mCurrentCurvature[i], dummy);
			}

			// Compute curvature & curvature derivative profile
			double curvStartLong[N];
			double curvEndLong[N];
			_interpolateCurvature(mControlPoints, mCurvatureStart, kdim, s, curvStartLong, N);
			_interpolateCurvature(mControlPoints, mCurvatureEnd, kdim, s, curvEndLong, N);
			for(int i=0;i<N;i++)
				_cubicInterpolation(mCurvatureTimeStart, mCurvatureTimeEnd, t, curvStartLong[i], curvEndLong[i], current_K_profile[i], current_dKdt_profile[i]);

			std::cout << "DURING TRANSITION (CURVATURE)" << std::endl;
		}
		else
		{
			// Set time start/end
			mCurvatureTimeStart = 0.0;
			mCurvatureTimeEnd = 0.0;

			// Compute current curvature
			for(int i=0;i<kdim;i++)
			{
				mCurvatureStart[i] = mCurvatureEnd[i];
				mCurrentCurvature[i] =  mCurvatureEnd[i];
			}

			// Compute curvature & curvature derivative profile
			_interpolateCurvature(mControlPoints, mCurrentCurvature, kdim, s, current_K_profile, N);

			for(int i=0;i<N;i++)
				current_dKdt_profile[i] = 0.0;

			std::cout << "END TRANSITION (CURVATURE)" << std::endl;

			return;
		}
	}


	// TRANSITION WITHIN TRANSITION
	if( !curvaturesStartEndEqual && !curvaturesEndCurrentEqual )
	{
		mCurvatureTimeStart = t;
		mCurvatureTimeEnd = mCurvatureTimeStart + timeTransition;

		for(int i=0;i<kdim;i++)
		{
			mCurvatureStart[i] = mCurrentCurvature[i];
			mCurvatureEnd[i] = curvature[i];
			double dummy = 0.0;
			_cubicInterpolation(mCurvatureTimeStart, mCurvatureTimeEnd, t, mCurvatureStart[i], mCurvatureEnd[i], mCurrentCurvature[i], dummy); // HERE CUBIC INTERPOLATION WITH SLOPE (to be implemented)!!!
		}

		// Compute curvature & curvature derivative profile
		double curvStartLong[N];
		double curvEndLong[N];
		_interpolateCurvature(mControlPoints, mCurvatureStart, kdim, s, curvStartLong, N);
		_interpolateCurvature(mControlPoints, mCurvatureEnd, kdim, s, curvEndLong, N);
		for(int i=0;i<N;i++)
			_cubicInterpolation(mCurvatureTimeStart, mCurvatureTimeEnd, t, curvStartLong[i], curvEndLong[i], current_K_profile[i], current_dKdt_profile[i]); // HERE CUBIC INTERPOLATION WITH SLOPE (to be implemented)!!!

		std::cout << "START TRANSITION WITHIN TRANSITION (CURVATURE)" << std::endl;

		return;
	}
}

void I2D_StefanFish::StefanFish::_schedulerBaseline(double t, double timeTransition, double * baseline, int kdim, double * s, int N, double * current_B_profile, double * current_dBdt_profile)
{
	// Curvature checks
	bool baselinesStartEndEqual = _checkCurvatureIsTheSame(mBaselineStart,mBaselineEnd,kdim);
	bool baselinesEndCurrentEqual = _checkCurvatureIsTheSame(mBaselineEnd,baseline,kdim);

	// KEEP GOING
	if( baselinesStartEndEqual && baselinesEndCurrentEqual )
	{
		// Set time start/end
		mBaselineTimeStart = 0.0;
		mBaselineTimeEnd = 0.0;

		// Compute current curvature
		for(int i=0;i<kdim;i++)
		{
			mBaselineStart[i] = mBaselineEnd[i];
			mCurrentBaseline[i] =  mBaselineEnd[i];
		}

		// Compute curvature & curvature derivative profile
		_interpolateCurvature(mControlPoints, mCurrentBaseline, kdim, s, current_B_profile, N);

		for(int i=0;i<N;i++)
			current_dBdt_profile[i] = 0.0;

		std::cout << "KEEP GOING (BASELINE)" << std::endl;

		return;
	}


	// START TRANSITION FROM STATIONARY
	if( baselinesStartEndEqual && !baselinesEndCurrentEqual )
	{
		// Set time start/end
		mBaselineTimeStart = t;
		mBaselineTimeEnd = mBaselineTimeStart + timeTransition;

		// Compute current curvature
		for(int i=0;i<kdim;i++)
		{
			mBaselineEnd[i] = baseline[i];
			double dummy = 0.0;
			_cubicInterpolation(mBaselineTimeStart, mBaselineTimeEnd, t, mBaselineStart[i], mBaselineEnd[i], mCurrentBaseline[i], dummy);
		}

		// Compute curvature & curvature derivative profile
		double baselineStartLong[N];
		double baselineEndLong[N];
		_interpolateCurvature(mControlPoints, mBaselineStart, kdim, s, baselineStartLong, N);
		_interpolateCurvature(mControlPoints, mBaselineEnd, kdim, s, baselineEndLong, N);
		for(int i=0;i<N;i++)
			_cubicInterpolation(mBaselineTimeStart, mBaselineTimeEnd, t, baselineStartLong[i], baselineEndLong[i], current_B_profile[i], current_dBdt_profile[i]);

		std::cout << "START TRANSITION FROM STATIONARY (BASELINE)" << std::endl;

		return;
	}

	// DURING TRANSITION & END TRANSITION
	if( !baselinesStartEndEqual && baselinesEndCurrentEqual )
	{
		if(t<=mBaselineTimeEnd)
		{
			// Compute current curvature
			for(int i=0;i<kdim;i++)
			{
				double dummy = 0.0;
				_cubicInterpolation(mBaselineTimeStart, mBaselineTimeEnd, t, mBaselineStart[i], mBaselineEnd[i], mCurrentBaseline[i], dummy);
			}

			// Compute curvature & curvature derivative profile
			double baselineStartLong[N];
			double baselineEndLong[N];
			_interpolateCurvature(mControlPoints, mBaselineStart, kdim, s, baselineStartLong, N);
			_interpolateCurvature(mControlPoints, mBaselineEnd, kdim, s, baselineEndLong, N);
			for(int i=0;i<N;i++)
				_cubicInterpolation(mBaselineTimeStart, mBaselineTimeEnd, t, baselineStartLong[i], baselineEndLong[i], current_B_profile[i], current_dBdt_profile[i]);

			std::cout << "DURING TRANSITION (BASELINE)" << std::endl;
		}
		else
		{
			// Set time start/end
			mBaselineTimeStart = 0.0;
			mBaselineTimeEnd = 0.0;

			// Compute current curvature
			for(int i=0;i<kdim;i++)
			{
				mBaselineStart[i] = mBaselineEnd[i];
				mCurrentBaseline[i] =  mBaselineEnd[i];
			}

			// Compute curvature & curvature derivative profile
			_interpolateCurvature(mControlPoints, mCurrentBaseline, kdim, s, current_B_profile, N);

			for(int i=0;i<N;i++)
				current_dBdt_profile[i] = 0.0;

			std::cout << "END TRANSITION (BASELINE)" << std::endl;

			return;
		}
	}

	// TRANSITION WITHIN TRANSITION
	if( !baselinesStartEndEqual && !baselinesEndCurrentEqual )
	{
		mBaselineTimeStart = t;
		mBaselineTimeEnd = mBaselineTimeStart + timeTransition;

		for(int i=0;i<kdim;i++)
		{
			mBaselineStart[i] = mCurrentBaseline[i];
			mBaselineEnd[i] = baseline[i];
			double dummy = 0.0;
			_cubicInterpolation(mBaselineTimeStart, mBaselineTimeEnd, t, mBaselineStart[i], mBaselineEnd[i], mCurrentBaseline[i], dummy);
		}

		// Compute curvature & curvature derivative profile
		double baselineStartLong[N];
		double baselineEndLong[N];
		_interpolateCurvature(mControlPoints, mBaselineStart, kdim, s, baselineStartLong, N);
		_interpolateCurvature(mControlPoints, mBaselineEnd, kdim, s, baselineEndLong, N);
		for(int i=0;i<N;i++)
			_cubicInterpolation(mBaselineTimeStart, mBaselineTimeEnd, t, baselineStartLong[i], baselineEndLong[i], current_B_profile[i], current_dBdt_profile[i]);

		std::cout << "START TRANSITION WITHIN TRANSITION (BASELINE)" << std::endl;

		return;
	}
}

void I2D_StefanFish::StefanFish::_schedulerT(double t, double timeTransition, double _T, double & current_T, double & current_T_derivative)
{
	// Tau checks
	bool TStartEndEqual = (mTStart == mTEnd);
	bool TEndCurrentEqual = (mTEnd == _T);

	// KEEP GOING
	if( TStartEndEqual && TEndCurrentEqual )
	{
		// Set time start/end
		mTTimeStart = 0.0;
		mTTimeEnd = 0.0;

		mTStart = mTEnd;
		mTCurrent = mTEnd;
		mTDerivative = 0.0;

		current_T = mTCurrent;
		current_T_derivative = mTDerivative;

		std::cout << "KEEP GOING (T)" << std::endl;

		return;
	}

	// START TRANSITION FROM STATIONARY
	if( TStartEndEqual && !TEndCurrentEqual )
	{
		// Set time start/end
		mTTimeStart = t;
		mTTimeEnd = mTTimeStart + timeTransition;

		// Compute current curvature
		mTEnd = _T;
		double dummy = 0.0;
		_cubicInterpolation(mTTimeStart, mTTimeEnd, t, mTStart, mTEnd, mTCurrent, mTDerivative);

		current_T = mTCurrent;
		current_T_derivative = mTDerivative;

		std::cout << "START TRANSITION FROM STATIONARY (T)" << std::endl;

		return;
	}


	// DURING TRANSITION & END TRANSITION
	if( !TStartEndEqual && TEndCurrentEqual )
	{
		if(t<=mTTimeEnd)
		{
			// Compute current curvature
			double dummy = 0.0;
			_cubicInterpolation(mTTimeStart, mTTimeEnd, t, mTStart, mTEnd, mTCurrent, mTDerivative);

			current_T = mTCurrent;
			current_T_derivative = mTDerivative;

			std::cout << "DURING TRANSITION (T)" << std::endl;
		}
		else
		{
			// Set time start/end
			mTTimeStart = 0.0;
			mTTimeEnd = 0.0;

			mTStart = mTEnd;
			mTCurrent = mTEnd;

			mTDerivative = 0.0;

			current_T = mTCurrent;
			current_T_derivative = mTDerivative;

			std::cout << "END TRANSITION (T)" << std::endl;

			return;
		}
	}


	// TRANSITION WITHIN TRANSITION
	if( !TStartEndEqual && !TEndCurrentEqual )
	{
		mTTimeStart = t;
		mTTimeEnd = mTTimeStart + timeTransition;

		mTStart = mTCurrent;
		mTEnd = _T;

		double dummy = 0.0;
		_cubicInterpolation(mTTimeStart, mTTimeEnd, t, mTStart, mTEnd, mTCurrent, mTDerivative);

		current_T = mTCurrent;
		current_T_derivative = mTDerivative;

		std::cout << "START TRANSITION WITHIN TRANSITION (T)" << std::endl;

		return;
	}
}

void I2D_StefanFish::StefanFish::_schedulerTau(double t, double timeTransition, double tau, double & current_tau, double & current_tau_derivative)
{
	// Tau checks
	bool tauStartEndEqual = (mTauStart == mTauEnd);
	bool tauEndCurrentEqual = (mTauEnd == tau);

	// KEEP GOING
	if( tauStartEndEqual && tauEndCurrentEqual )
	{
		// Set time start/end
		mTauTimeStart = 0.0;
		mTauTimeEnd = 0.0;

		mTauStart = mTauEnd;
		mTauCurrent = mTauEnd;
		mTauDerivative = 0.0;

		current_tau = mTauCurrent;
		current_tau_derivative = mTauDerivative;

		std::cout << "KEEP GOING (TAU)" << std::endl;

		return;
	}

	// START TRANSITION FROM STATIONARY
	if( tauStartEndEqual && !tauEndCurrentEqual )
	{
		// Set time start/end
		mTauTimeStart = t;
		mTauTimeEnd = mTauTimeStart + timeTransition;

		// Compute current curvature
		mTauEnd = tau;
		double dummy = 0.0;
		_cubicInterpolation(mTauTimeStart, mTauTimeEnd, t, mTauStart, mTauEnd, mTauCurrent, mTauDerivative);

		current_tau = mTauCurrent;
		current_tau_derivative = mTauDerivative;

		std::cout << "START TRANSITION FROM STATIONARY (TAU)" << std::endl;

		return;
	}


	// DURING TRANSITION & END TRANSITION
	if( !tauStartEndEqual && tauEndCurrentEqual )
	{
		if(t<=mTauTimeEnd)
		{
			// Compute current curvature
			double dummy = 0.0;
			_cubicInterpolation(mTauTimeStart, mTauTimeEnd, t, mTauStart, mTauEnd, mTauCurrent, mTauDerivative);

			current_tau = mTauCurrent;
			current_tau_derivative = mTauDerivative;

			std::cout << "DURING TRANSITION (TAU)" << std::endl;
		}
		else
		{
			// Set time start/end
			mTauTimeStart = 0.0;
			mTauTimeEnd = 0.0;

			mTauStart = mTauEnd;
			mTauCurrent = mTauEnd;

			mTauDerivative = 0.0;

			current_tau = mTauCurrent;
			current_tau_derivative = mTauDerivative;

			std::cout << "END TRANSITION (TAU)" << std::endl;

			return;
		}
	}


	// TRANSITION WITHIN TRANSITION
	if( !tauStartEndEqual && !tauEndCurrentEqual )
	{
		mTauTimeStart = t;
		mTauTimeEnd = mTauTimeStart + timeTransition;

		mTauStart = mTauCurrent;
		mTauEnd = tau;

		double dummy = 0.0;
		_cubicInterpolation(mTauTimeStart, mTauTimeEnd, t, mTauStart, mTauEnd, mTauCurrent, mTauDerivative);

		current_tau = mTauCurrent;
		current_tau_derivative = mTauDerivative;

		std::cout << "START TRANSITION WITHIN TRANSITION (TAU)" << std::endl;

		return;
	}
}

void I2D_StefanFish::StefanFish::_solveFrenetInHeadRefFrame(const double ds, const double * k, double * rX, double * rY, double * tanX, double * tanY, double * norX, double * norY, const unsigned int n)
{
	I2D_Frenet2D frenetSolver;
	frenetSolver.solve(ds, k, rX, rY, tanX, tanY, norX, norY, n);
}

void I2D_StefanFish::StefanFish::_solveVelProfileInHeadRefFrame(const double ds, const double * kt, const double * kt_dt, const double * tanX, const double * tanY, const double * norX, const double * norY,
		double * vX, double * vY, double * vTanX, double * vTanY, double * vNorX, double * vNorY, const unsigned int n)
{
	I2D_Frenet2D frenetSolver;
	frenetSolver.solveVelProfile(ds, kt, kt_dt, tanX, tanY, norX, norY, vX, vY, vTanX, vTanY, vNorX, vNorY, n);
}

void I2D_StefanFish::StefanFish::_startWave(double time)
{
	if(mBoolWaveStarted==false)
	{
		mTimeStartWave = time;
		mBoolWaveStarted = true;
	}
}

void I2D_StefanFish::StefanFish::clear()
{
	// Curvature transition
	mCurvatureStart[0] = 0.0;
	mCurvatureStart[1] = 0.0;
	mCurvatureStart[2] = 0.0;
	mCurvatureStart[3] = 0.0;
	mCurvatureStart[4] = 0.0;
	mCurvatureStart[5] = 0.0;

	mCurvatureEnd[0] = 0.0;
	mCurvatureEnd[1] = 0.0;
	mCurvatureEnd[2] = 0.0;
	mCurvatureEnd[3] = 0.0;
	mCurvatureEnd[4] = 0.0;
	mCurvatureEnd[5] = 0.0;

	mCurrentCurvature[0] = 0.0;
	mCurrentCurvature[1] = 0.0;
	mCurrentCurvature[2] = 0.0;
	mCurrentCurvature[3] = 0.0;
	mCurrentCurvature[4] = 0.0;
	mCurrentCurvature[5] = 0.0;

	mCurvatureTimeStart = 0.0;
	mCurvatureTimeEnd = 0.0;

	memset(TANX,0,N*sizeof(double));
	memset(TANY,0,N*sizeof(double));
	memset(VTANX,0,N*sizeof(double));
	memset(VTANY,0,N*sizeof(double));
	memset(K,0,N*sizeof(double));
	memset(KT,0,N*sizeof(double));
	memset(DKDT,0,N*sizeof(double));
	memset(mCurvatureDerivative,0,N*sizeof(double));

	// Baseline transition
	mBaselineStart[0] = 0.0;
	mBaselineStart[1] = 0.0;
	mBaselineStart[2] = 0.0;
	mBaselineStart[3] = 0.0;
	mBaselineStart[4] = 0.0;
	mBaselineStart[5] = 0.0;

	mBaselineEnd[0] = 0.0;
	mBaselineEnd[1] = 0.0;
	mBaselineEnd[2] = 0.0;
	mBaselineEnd[3] = 0.0;
	mBaselineEnd[4] = 0.0;
	mBaselineEnd[5] = 0.0;

	mCurrentBaseline[0] = 0.0;
	mCurrentBaseline[1] = 0.0;
	mCurrentBaseline[2] = 0.0;
	mCurrentBaseline[3] = 0.0;
	mCurrentBaseline[4] = 0.0;
	mCurrentBaseline[5] = 0.0;

	mBaselineTimeStart = 0.0;
	mBaselineTimeEnd = 0.0;

	memset(mBaselineInterpolated,0,N*sizeof(double));
	memset(mBaselineDerivative,0,N*sizeof(double));

	// Tau transition
	mTauTimeStart = 0.0;
	mTauTimeEnd = 0.0;
	mTauStart = 0.0;
	mTauEnd = 0.0;
	mTauCurrent = 0.0;
	mTauDerivative = 0.0;

	// Wave transition
	mTimeStartWave = 0.0;
	mBoolWaveStarted = false;
}

void I2D_StefanFish::StefanFish::updateInSpace(const double t)
{
	clearUniformGrids();

	// Head reference frame quantities
	double hCMX = 0.0;
	double hCMY = 0.0;
	double hVCMX = 0.0;
	double hVCMY = 0.0;

	// Center of mass reference frame quantities
	double cmL = 0.0;
	double cmII = 0.0;
	double omega = 0.0;

	// Frame of reference shift
	const double centerx = 0.5;
	const double centery = 0.5;
	double yyy = 0.0;
	double xxx = 0.0;

	// TEMPORARY QUANTITIES
	double curvature[this->CURVATURE.size()];
	double baseline[this->BASELINE.size()];

	for(unsigned int i=0; i<this->CURVATURE.size(); i++ )
	{
		curvature[i] = this->CURVATURE[i];
		baseline[i] = this->BASELINE[i];
	}

	// Smooth change in configuration
	double current_tau = 0.0;
	double current_tau_derivative = 0.0;
	double current_T = 0.0;
	double current_T_derivative = 0.0;
	assert(mTransitionTime==T);
	_schedulerBaseline(t,mTransitionTime,baseline,6,S,N,mBaselineInterpolated,mBaselineDerivative);
	_schedulerCurvature(t,mTransitionTime,curvature,6,S,N,K,mCurvatureDerivative);
	_schedulerTau(t,mTransitionTime,tau,current_tau,current_tau_derivative);
	_schedulerT(t,mTransitionTime,T,current_T,current_T_derivative);

	// Calculations in the head reference frame
	double wavetime = t;
	_startWave(t);
	if(mBoolWaveStarted==true){ wavetime = t - mTimeStartWave; }
	_applyTravellingWave(wavetime, current_tau, current_tau_derivative, current_T, current_T_derivative, S, K, mBaselineInterpolated, mCurvatureDerivative, mBaselineDerivative, KT, DKDT, N, mBoolWaveStarted);
	_solveFrenetInHeadRefFrame(DS, KT, X, Y, TANX, TANY, NORX, NORY, N);

	// Print midline
	//FILE * ppFile = fopen("MRAGmidline", "a");
	//fprintf(ppFile,"%f ",t);
	//for(unsigned int i=0; i<N; i++) { fprintf(ppFile,"%f ",X[i]); }
	//fprintf(ppFile,"\n");
	//fprintf(ppFile,"%f ",t);
	//for(unsigned int i=0; i<N; i++) { fprintf(ppFile,"%f ",Y[i]); }
	//fprintf(ppFile,"\n");
	//fclose(ppFile);

	_solveVelProfileInHeadRefFrame(DS, KT, DKDT, TANX, TANY, NORX, NORY, VX, VY, VTANX, VTANY, VNORX, VNORY, N);
	_createShapeBoundary(X, Y, W, NORX, NORY, SHAPEX, SHAPEY, N);
	_fillDefGrid();

	// Go into center of mass frame of reference
	_getCenterOfMassFull(hCMX, hCMY); xxx = hCMX; yyy = hCMY;
	if(xxx<0.0){ printf("cazzo xxx negative!!\n"); abort(); }

	_getVelCenterOfMassFull(hVCMX, hVCMY);
	_centerlineCenterOfMassFrameTransform(hCMX, hCMY, hVCMX, hVCMY);
	_defGridCenterOfMassFrameTransform(hCMX, hCMY, hVCMX, hVCMY);

	// Correction for the rotational impulse in the center of mass reference frame
	_getAngularMomentumFull(cmL);
	_getScalarMomentOfInertiaFull(cmII);
	omega = -cmL/cmII;
	_correctCenterlineForRotationalImpulse(omega);
	_correctDefGridForRotationalImpulse(omega);

#ifndef NDEBUG
	// Check, the angular momentum here should be almost zero
	double hL = 0.0;
	_getCenterOfMassFull(hCMX, hCMY);
	_getVelCenterOfMassFull(hVCMX, hVCMY);
	_getAngularMomentumFull(hL);
	if( (fabs(hL) > 1e-10) || (fabs(hCMX) > 1e-10) || (fabs(hCMY) > 1e-10) || (fabs(hVCMX) > 1e-10) || (fabs(hVCMY) > 1e-10) )
	{
		printf("(Everything should be zero!) cmx=%e, cmy=%e\n", hCMX, hCMY);
		printf("(Everything should be zero!) vxcm=%e, vycm=%e\n", hVCMX, hVCMY);
		printf("(Everything should be zero!) L=%e\n", hL);
		abort();
	}
#endif

	// Rotate for the internal angle
	_rotateAboutTheOrigin(angleInSpace, X, Y, N);
	_rotateDefGridAboutTheOrigin(angleInSpace);

	// OPTIONAL
	_rotateAboutTheOrigin(angleInSpace, VX, VY, N);
	_rotateAboutTheOrigin(angleInSpace, NORX, NORY, N);
	_createShapeBoundary(X, Y, W, NORX, NORY, SHAPEX, SHAPEY, N);

	// Compute bounding fish box
	double height = 0.0;
	double width = 0.0;
	_bboxShapeBoundary(SHAPEX, SHAPEY, N, width, height);
	const double spaceX = (1.0-width)/2.0;
	const double spaceY = (1.0-height)/2.0;
	xxx = (xxx<=0.0)?(1.0-spaceX+xxx):xxx+spaceX;
	yyy = (yyy<=0.0)?(1.0-spaceY+yyy):yyy+spaceY;

	// Fish should fit in domain [0,1]
	_rigidTranslation(xxx,yyy);
	traslX = xxx-centerx;
	traslY = yyy-centery;
	angular_velocityInSpace = omega;

	// Compute bounding fish box
	_bbox();

	// Reconstruct shape on uniform fine mesh
	bilinearInterpolation();
}

I2D_StefanFish::StefanFish::StefanFish(double xm, double ym, double _D, double _T, double phase, double tau, double angle_rad, vector<double> BASELINE, vector<double> CURVATURE, double angleInSpace_rad,
					double eps, const int LMAX, const bool isSharp): tau(tau), mTransitionTime(_T), Fish(xm,ym,_D,_T,phase,angle_rad,angleInSpace_rad,eps,LMAX, isSharp)
{
	this->BASELINE = BASELINE;
	this->CURVATURE = CURVATURE;

	printf("--------------SUMMARY STEFANFISH--------------\n");
	printf("xm=%f\n",xm);
	printf("ym=%f\n",ym);
	printf("D=%f\n",_D);
	printf("T=%f\n",_T);
	printf("phase=%f\n",phase);
	printf("tau=%f\n",tau);
	printf("angle_rad=%f\n",angle_rad);
	printf("angleInSpace_rad=%f\n",angleInSpace_rad);
	printf("eps=%f\n",eps);
	printf("LMAX=%d\n",LMAX);
	printf("BASELINE\n");
	for(unsigned int i=0; i<this->BASELINE.size();i++){ printf("%f, ", this->BASELINE[i]); }
	printf("\nCURVATURE\n");
	for(unsigned int i=0; i<this->CURVATURE.size();i++){ printf("%f, ", this->CURVATURE[i]); }
	printf("\n--------------SUMMARY STEFANFISH--------------\n");

	if( this->BASELINE.size()!=this->CURVATURE.size() || this->BASELINE.size()!=6){ printf("Something wrong with BASELINE or CURVATURE!\n"); abort(); }

	S1 = 0.05*1.0;
	S2 = 0.33*1.0;
	S3 = 0.67*1.0;
	S4 = 0.95*1.0;

	mControlPoints[0] = 0.0;
	mControlPoints[1] = S1;
	mControlPoints[2] = S2;
	mControlPoints[3] = S3;
	mControlPoints[4] = S4;
	mControlPoints[5] = 1.0;

	// Curvature transition
	mCurvatureStart[0] = 0.0;
	mCurvatureStart[1] = 0.0;
	mCurvatureStart[2] = 0.0;
	mCurvatureStart[3] = 0.0;
	mCurvatureStart[4] = 0.0;
	mCurvatureStart[5] = 0.0;

	mCurvatureEnd[0] = 0.0;
	mCurvatureEnd[1] = 0.0;
	mCurvatureEnd[2] = 0.0;
	mCurvatureEnd[3] = 0.0;
	mCurvatureEnd[4] = 0.0;
	mCurvatureEnd[5] = 0.0;

	mCurrentCurvature[0] = 0.0;
	mCurrentCurvature[1] = 0.0;
	mCurrentCurvature[2] = 0.0;
	mCurrentCurvature[3] = 0.0;
	mCurrentCurvature[4] = 0.0;
	mCurrentCurvature[5] = 0.0;

	mCurvatureTimeStart = 0.0;
	mCurvatureTimeEnd = 0.0;

	TANX = new double[N];
	TANY = new double[N];
	VTANX = new double[N];
	VTANY = new double[N];
	K = new double[N];
	KT = new double[N];
	DKDT = new double[N];
	mCurvatureDerivative = new double[N];
	memset(TANX,0,N*sizeof(double));
	memset(TANY,0,N*sizeof(double));
	memset(VTANX,0,N*sizeof(double));
	memset(VTANY,0,N*sizeof(double));
	memset(K,0,N*sizeof(double));
	memset(KT,0,N*sizeof(double));
	memset(DKDT,0,N*sizeof(double));
	memset(mCurvatureDerivative,0,N*sizeof(double));

	// Baseline transition
	mBaselineStart[0] = 0.0;
	mBaselineStart[1] = 0.0;
	mBaselineStart[2] = 0.0;
	mBaselineStart[3] = 0.0;
	mBaselineStart[4] = 0.0;
	mBaselineStart[5] = 0.0;

	mBaselineEnd[0] = 0.0;
	mBaselineEnd[1] = 0.0;
	mBaselineEnd[2] = 0.0;
	mBaselineEnd[3] = 0.0;
	mBaselineEnd[4] = 0.0;
	mBaselineEnd[5] = 0.0;

	mCurrentBaseline[0] = 0.0;
	mCurrentBaseline[1] = 0.0;
	mCurrentBaseline[2] = 0.0;
	mCurrentBaseline[3] = 0.0;
	mCurrentBaseline[4] = 0.0;
	mCurrentBaseline[5] = 0.0;

	mBaselineTimeStart = 0.0;
	mBaselineTimeEnd = 0.0;

	mBaselineInterpolated = new double[N];
	mBaselineDerivative = new double[N];
	memset(mBaselineInterpolated,0,N*sizeof(double));
	memset(mBaselineDerivative,0,N*sizeof(double));

	// Tau transition
	mTauTimeStart = 0.0;
	mTauTimeEnd = 0.0;
	mTauStart = tau;
	mTauEnd = tau;
	mTauCurrent = tau;
	mTauDerivative = 0.0;

	// T transition
	mTTimeStart = 0.0;
	mTTimeEnd = 0.0;
	mTStart = _T;
	mTEnd = _T;
	mTCurrent = _T;
	mTDerivative = 0.0;

	// Wave transition
	mTimeStartWave = 0.0;
	mBoolWaveStarted = false;
}

I2D_StefanFish::StefanFish::~StefanFish()
{
	delete [] K;
	delete [] KT;
	delete [] DKDT;
	delete [] TANX;
	delete [] TANY;
	delete [] VTANX;
	delete [] VTANY;
	delete [] mCurvatureDerivative;
	delete [] mBaselineInterpolated;
	delete [] mBaselineDerivative;
}

I2D_StefanFish::I2D_StefanFish(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real T, const Real phase, const Real tau, const Real angle, vector<double> BASELINE, vector<double> CURVATURE,
		const Real angleInSpace, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID, RL::RL_TabularPolicy ** policy, const int seed):
		I2D_CarlingFish(parser, grid, _xm, _ym, D, T, phase, angle, angleInSpace, eps, Uinf, penalization, LMAX, ID, policy, seed)
{
	assert(shape!=NULL);
	if(shape!=NULL)
	{
		delete shape;
		shape = NULL;
	}

	const bool isSharp = parser("-sharp").asBool();

	shape = new StefanFish(_xm, _ym, D, T, phase, tau, angle, BASELINE, CURVATURE, angleInSpace, eps, LMAX, isSharp);
	shape->updateInSpace(0);
	//shape->clear();
}

