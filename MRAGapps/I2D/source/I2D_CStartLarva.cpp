/*
 * I2D_CStartLarva.cpp
 *
 *  Created on: Dec 8, 2011
 *      Author: mgazzola
 */

#include "I2D_CStartLarva.h"
#include "I2D_Interpolation1D.h"
#include "I2D_Frenet2D.h"
#include "I2D_VectorBlockLab.h"
#include "I2D_GradOfVector.h"
#include <limits>
#include <omp.h>

I2D_CStartLarva::CStartLarva::CStartLarva(double xm, double ym, double _D, double Tprep, double Tprop, double phase, double tau, double angle_rad, vector<double> BASELINE, vector<double> CURVATURE, double angleInSpace_rad,
					   double eps, const int LMAX, const bool isSharp):  Tprep(Tprep), Tprop(Tprop), StefanFish(xm, ym, _D, Tprop, phase, tau, angle_rad, BASELINE, CURVATURE, angleInSpace_rad, eps, LMAX, isSharp)
{
	this->BASELINE = BASELINE;
	this->CURVATURE = CURVATURE;

	printf("--------------SUMMARY CSTARTLARVA--------------\n");
	printf("xm=%f\n",xm);
	printf("ym=%f\n",ym);
	printf("D=%f\n",_D);
	printf("T=%f\n",T);
	printf("Tprep=%f\n",Tprep);
	printf("Tprop=%f\n",Tprop);
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
	printf("\n--------------SUMMARY CSTARTLARVA--------------\n");

	if( this->BASELINE.size()!=this->CURVATURE.size() || this->BASELINE.size()!=6){ printf("Something wrong with BASELINE or CURVATURE!\n"); abort(); }

	this->WH = (2.0/31.5)*(1.0-EXTENSION*EPS);
	this->SB = (2.5/29.0)*(1.0-EXTENSION*EPS);
	this->WT = (0.8/31.5)*(1.0-EXTENSION*EPS);
	this->ST = (10.0/29.0)*(1.0-EXTENSION*EPS);
	this->S1 = 0.2*1.0;
	this->S2 = 0.5*1.0;
	this->S3 = 0.75*1.0;
	this->S4 = 0.95*1.0;
	this->mControlPoints[0] = 0.0;
	this->mControlPoints[1] = S1;
	this->mControlPoints[2] = S2;
	this->mControlPoints[3] = S3;
	this->mControlPoints[4] = S4;
	this->mControlPoints[5] = 1.0;

	for(int i = 0; i < N; i++ )
	{
		S[i] = DS*(double)i;
		W[i] = _getWidth(S[i]);
	}
}

double I2D_CStartLarva::CStartLarva::_getWidth(const double & ss) const
{
	const double end = 1.0-EXTENSION*EPS;
	const double s = ss - (EXTENSION/2.0)*EPS;
	const double a = (-2.0*(WT-WH) - WT*(ST-SB)) / ((ST-SB)*(ST-SB)*(ST-SB));
	const double b = (3.0*(WT-WH) + WT*(ST-SB)) / ((ST-SB)*(ST-SB));
	const double d = WH;
	double width = 0.0;

	if( (s>=0.0) && (s<SB) )
		width = WH*sqrt(1.0 - ((SB - s)/SB)*((SB - s)/SB));

	// CUBIC
	if( (s>=SB) && (s<ST) )
		width = a*(s-SB)*(s-SB)*(s-SB)+b*(s-SB)*(s-SB)+d;

	// QUADRATIC
	if( (s>=ST) && (s<=end) )
		width = WT - (WT)*((s-ST)/(end-ST))*((s-ST)/(end-ST));

	return width;
}

void I2D_CStartLarva::CStartLarva::updateInSpace(const double t)
{
	const double time = t;

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
	double timeTransition = this->T;
	double curvature[this->CURVATURE.size()];
	double baseline[this->BASELINE.size()];

	for(unsigned int i=0; i<this->CURVATURE.size(); i++ )
		curvature[i] = this->CURVATURE[i];

	if(time<=Tprep)
	{
		timeTransition = Tprep;
		for(unsigned int i=0; i<this->BASELINE.size(); i++ )
			baseline[i] = this->BASELINE[i];
	}
	else
	{
		timeTransition = Tprop;
		for(unsigned int i=0; i<this->BASELINE.size(); i++ )
			baseline[i] = 0.0;
	}

	// Smooth change in configuration
	double current_tau = 0.0;
	double current_tau_derivative = 0.0;
	double current_T = this->T;
	double current_T_derivative = 0.0;
	_startWave(time);
	_schedulerBaseline(time,timeTransition,baseline,6,S,N,mBaselineInterpolated,mBaselineDerivative);
	_schedulerCurvature(time,timeTransition,curvature,6,S,N,K,mCurvatureDerivative);
	_schedulerTau(time,timeTransition,tau, current_tau, current_tau_derivative);

	// Calculations in the head reference frame
	double wavetime = time;
	if(mBoolWaveStarted==true){ wavetime = time - mTimeStartWave; }
	_applyTravellingWave(wavetime, current_tau, current_tau_derivative, current_T, current_T_derivative, S, K, mBaselineInterpolated, mCurvatureDerivative, mBaselineDerivative, KT, DKDT, N, mBoolWaveStarted);
	_solveFrenetInHeadRefFrame(DS, KT, X, Y, TANX, TANY, NORX, NORY, N);

	// Print midline
	//FILE * ppFile = fopen("MRAGmidline", "a");
	//fprintf(ppFile,"%f ",time);
	//for(unsigned int i=0; i<N; i++) { fprintf(ppFile,"%f ",X[i]); }
	//fprintf(ppFile,"\n");
	//fprintf(ppFile,"%f ",time);
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
	printf("traslY=%e, yyy=%e\n", traslY, yyy);
	angular_velocityInSpace = omega;

	// Compute bounding fish box
	_bbox();

	// Reconstruct shape on uniform fine mesh
	bilinearInterpolation();
}

I2D_CStartLarva::I2D_CStartLarva(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real Tprep, const Real Tprop, const Real phase, const Real tau, const Real angle, vector<double> BASELINE, vector<double> CURVATURE,
		const Real angleInSpace, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID, RL::RL_TabularPolicy ** policy):
		I2D_StefanFish(parser, grid, _xm, _ym, D, Tprop, phase, tau, angle, BASELINE, CURVATURE, angleInSpace, eps, Uinf, penalization, LMAX, ID, policy)
{
	assert(shape!=NULL);
	if(shape!=NULL)
	{
		delete shape;
		shape = NULL;
	}

	const bool isSharp = parser("-sharp").asBool();

	shape = new CStartLarva(_xm, _ym, D, Tprep, Tprop, phase, tau, angle, BASELINE, CURVATURE, angleInSpace, eps, LMAX, isSharp);
	shape->updateInSpace(0.0);
	//shape->clear();
}

