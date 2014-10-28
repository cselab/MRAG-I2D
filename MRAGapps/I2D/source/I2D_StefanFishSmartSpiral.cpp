/*
 * I2D_StefanFishSmartSpiral.cpp
 *
 *  Created on: Mar 6, 2012
 *      Author: mgazzola
 */

#include "I2D_StefanFishSmartSpiral.h"

I2D_StefanFishSmartSpiral::I2D_StefanFishSmartSpiral(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real T, const Real tau, const Real angle, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID, RL::RL_TabularPolicy ** _policy, const int seed) :
I2D_StefanFishSmart(parser, grid, _xm, _ym, D, T, tau, angle, eps, Uinf, penalization, LMAX, ID, _policy, seed)
{
	maxDomainRadius = 0.4;

	// Spiral paramater
	spiralT = 100.0;
	spiralAngle = 3.0*M_PI;
	omega = spiralAngle/spiralT;
	A = 0.1;

	target[0] = (double)_xm;
	target[1] = (double)_ym;

	assert(shape!=NULL);
	if(shape!=NULL)
	{
		delete shape;
		shape = NULL;
	}

	// Set sensors
	signature.clear();
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			signature.push_back(5); // dist from lattice point
			signature.push_back(16); // angle between myVel and lattice point
			signature.push_back(3); // actions

			if(!(*policy)->samedim(signature))
			{
				printf("Policy dimension reset!\n");
				(*policy)->setdim(signature);
			}

			learningInterval = T;
		}

	shape = new StefanFishSmart(_xm, _ym, D, T, tau, angle, eps, LMAX);
	shape->updateInSpace(0.0);
}

I2D_StefanFishSmartSpiral::~I2D_StefanFishSmartSpiral()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

bool I2D_StefanFishSmartSpiral::mapState(vector<int> & state, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	assert(_data!=NULL);
	assert((*_data).size()!=0);

	state.clear();

	// Gather my info (always double check with getInfo!!)
	const double center[2] = {0.5,0.5};
	const double myX[2] = {shape->xm,shape->ym};
	const double myV[2] = {shape->vx,shape->vy};

	// Levels
	const int levelX = signature[0];
	const int levelAngle = signature[1];

	// Constants
	bool valid = true;
	const double cutoff = shape->D/2.0;
	const double h = cutoff / (double)levelX;
	const double dangle = 360.0 / (double)levelAngle;

	double d[2] = {0.0,0.0};
	_dist(target,myX,d);

	const double angle = _angleVectors(myV, d);

	const double IdI = sqrt(d[0]*d[0]+d[1]*d[1]);

	const int idxX = std::max(0,std::min(levelX-1,(int)floor(IdI/h)));
	const int idxangle = std::max(0,std::min(levelAngle-1,(int)floor(angle/dangle)));

	// Prepare state vector
	state.push_back(idxX);
	state.push_back(idxangle);

	{
		const double dist = sqrt( (myX[0]-center[0])*(myX[0]-center[0]) + (myX[1]-center[1])*(myX[1]-center[1]) );
		if( dist > maxDomainRadius || IdI>cutoff )
			valid = false;
	}

	return valid;
}

void I2D_StefanFishSmartSpiral::reward(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	assert(status==Waiting);

	if( t>(learningTimer+learningInterval) && status==Waiting)
	{
		double reward = 0.0;

		// Gather my info (always double check with getInfo!!)
		const double myX[2] = {shape->xm,shape->ym};

		// Constants
		const double cutoff = shape->D/2.0;
		const double spotR = cutoff/2.0;
		const double rewardOnSpot = 5.0;

		double d[2] = {0.0,0.0};
		_dist(target,myX,d);

		const double IdI = sqrt(d[0]*d[0]+d[1]*d[1]);

		if(IdI <= spotR)
			reward += rewardOnSpot*(1-IdI/spotR);
		else
			reward += rewardOnSpot*(1-IdI/spotR);

		integralReward += reward;
	}
}

void I2D_StefanFishSmartSpiral::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	// If string is set I open the corresponding file
	assert(filename!=std::string());
	FILE * ppFile = NULL;
	ppFile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
	assert(ppFile!=NULL);
	fprintf(ppFile, "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", t, shape->xm, shape->ym, shape->vx, shape->vy, shape->angle, shape->angular_velocity, target[0], target[1]);
	fflush(ppFile);
	fclose(ppFile);

	// Cumulative diagnostics
	string cumulative = filename + "_cum";
	ppFile = fopen(cumulative.c_str(),"a");
	assert(ppFile!=NULL);
	fprintf(ppFile, "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", t, shape->xm, shape->ym, shape->vx, shape->vy, shape->angle, shape->angular_velocity, target[0], target[1]);
	fflush(ppFile);
	fclose(ppFile);

	shape->update_all(dt);

	// Update lattice node
	const double center[2] = {0.5,0.5};
	double r = 0.0;
	double theta = 0.0;

	if(t<=spiralT)
	{
		theta = spiralAngle - omega*t;
		r = A*sqrt(theta);
	}
	else
	{
		theta = omega*(t-spiralT);
		r = -A*sqrt(theta);
	}

	const double xt = r*cos(theta);
	const double yt = r*sin(theta);

	target[0] = xt + center[0];
	target[1] = yt + center[1];
}

void I2D_StefanFishSmartSpiral::mapAction(int action)
{
	assert(shape!=NULL);
	StefanFish * myshape = static_cast<StefanFish*>(shape);

	assert(myshape->CURVATURE.size()==6);
	assert(myshape->BASELINE.size()==6);
	assert(learningInterval>=0.0);
	myshape->mTransitionTime = learningInterval;
	myshape->tau = 1.44;
	myshape->CURVATURE[0] = 0.0;
	myshape->CURVATURE[1] = 1.51/3.0;
	myshape->CURVATURE[2] = 0.48/3.0;
	myshape->CURVATURE[3] = 5.74/3.0;
	myshape->CURVATURE[4] = 2.73/3.0;
	myshape->CURVATURE[5] = 0.0;

	switch(action)
	{
	case 0 :
		printf("Go straight!!\n");
		myshape->BASELINE[0] = 0.0;
		myshape->BASELINE[1] = 0.0;
		myshape->BASELINE[2] = 0.0;
		myshape->BASELINE[3] = 0.0;
		myshape->BASELINE[4] = 0.0;
		myshape->BASELINE[5] = 0.0;
		break;
	case 1 :
		printf("Turn left!!\n");
		myshape->BASELINE[0] = -1.0;
		myshape->BASELINE[1] = -1.0;
		myshape->BASELINE[2] = -1.0;
		myshape->BASELINE[3] = -1.0;
		myshape->BASELINE[4] = -1.0;
		myshape->BASELINE[5] = -1.0;
		break;
	case 2 :
		printf("Turn right!!\n");
		myshape->BASELINE[0] = 1.0;
		myshape->BASELINE[1] = 1.0;
		myshape->BASELINE[2] = 1.0;
		myshape->BASELINE[3] = 1.0;
		myshape->BASELINE[4] = 1.0;
		myshape->BASELINE[5] = 1.0;
		break;
	default :
		printf("What da fuck did u map, dude!?\n");
		break;
	}
}

