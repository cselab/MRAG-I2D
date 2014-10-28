/*
 * I2D_StefanFishSmartInline.cpp
 *
 *  Created on: Feb 21, 2012
 *      Author: mgazzola
 */

#include "I2D_StefanFishSmartInline.h"

I2D_StefanFishSmartInline::I2D_StefanFishSmartInline(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real T, const Real tau, const Real angle, const Real _dir[2], const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID, RL::RL_TabularPolicy ** _policy, const int seed) :
I2D_StefanFishSmart(parser, grid, _xm, _ym, D, T, tau, angle, eps, Uinf, penalization, LMAX, ID, _policy, seed)
{
	const double IdirI = sqrt( _dir[0]*_dir[0] + _dir[1]*_dir[1] );
	assert(IdirI!=0.0);

	dir[0] = _dir[0]/IdirI;
	dir[1] = _dir[1]/IdirI;

	target[0] = (double)_xm + D*dir[0];
	target[1] = (double)_ym + D*dir[1];

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
			signature.push_back(20);	// ---> dist from lattice point with cutoff
			signature.push_back(16);	// ---> angle between myVel and lattice point
			signature.push_back(3);		// ---> actions

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

I2D_StefanFishSmartInline::~I2D_StefanFishSmartInline()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

bool I2D_StefanFishSmartInline::mapState(vector<int> & state, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
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
	const double cutoff = shape->D;
	const double h = cutoff / (double)levelX;
	const double dangle = 360.0 / (double)levelAngle;

	double d[2] = {0.0,0.0};
	_dist(target,myX,d);

	const double angle = _angleVectors(myV, d);

	const double IdI = sqrt(d[0]*d[0]+d[1]*d[1]) - shape->D;

	const int idxX = std::max(0,std::min(levelX-1,(int)floor(IdI/h)));
	const int idxangle = std::max(0,std::min(levelAngle-1,(int)floor(angle/dangle)));

	// Prepare state vector
	state.push_back(idxX);
	state.push_back(idxangle);

	bool valid = true;
	{
		const double dist = sqrt( (myX[0]-center[0])*(myX[0]-center[0]) + (myX[1]-center[1])*(myX[1]-center[1]) );
		if( dist > maxDomainRadius || IdI>(5*cutoff) )
			valid = false;
	}

	return valid;
}

void I2D_StefanFishSmartInline::reward(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	assert(shape!=NULL);
	StefanFish * myshape = static_cast<StefanFish*>(shape);

	assert(myshape->mTransitionTime==shape->T);
	assert(myshape->mTransitionTime==learningInterval);

	assert(status==Waiting);

	if( t>(learningTimer+learningInterval) && status==Waiting)
	{
		double reward = 0.0;

		// Gather my info (always double check with getInfo!!)
		const double myX[2] = {shape->xm,shape->ym};

		// Constants
		const double cutoff = shape->D;
		const double rewardOnSpot = 5.0;

		double d[2] = {0.0,0.0};
		_dist(target,myX,d);

		const double IdI2 = d[0]*d[0]+d[1]*d[1] - shape->D*shape->D;
		const double a = (rewardOnSpot/(cutoff*cutoff));
		reward += -a*IdI2 + rewardOnSpot;

		integralReward += reward;
	}
}

void I2D_StefanFishSmartInline::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
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
	const double modv = shape->vx*dir[0] + shape->vy*dir[1];
	target[0] += modv*dir[0]*dt;
	target[1] += modv*dir[1]*dt;
}

void I2D_StefanFishSmartInline::mapAction(int action)
{
	assert(shape!=NULL);
	StefanFish * myshape = static_cast<StefanFish*>(shape);

	assert(myshape->CURVATURE.size()==6);
	assert(myshape->BASELINE.size()==6);
	assert(learningInterval>=0.0);

	myshape->CURVATURE[0] = 0.0;
	myshape->CURVATURE[1] = 1.51/3.0;
	myshape->CURVATURE[2] = 0.48/3.0;
	myshape->CURVATURE[3] = 5.74/3.0;
	myshape->CURVATURE[4] = 2.73/3.0;
	myshape->CURVATURE[5] = 0.0;

	const Real bending = 0.5;

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
		myshape->BASELINE[0] = -bending;
		myshape->BASELINE[1] = -bending;
		myshape->BASELINE[2] = -bending;
		myshape->BASELINE[3] = -bending;
		myshape->BASELINE[4] = -bending;
		myshape->BASELINE[5] = -bending;
		break;
	case 2 :
		printf("Turn right!!\n");
		myshape->BASELINE[0] = bending;
		myshape->BASELINE[1] = bending;
		myshape->BASELINE[2] = bending;
		myshape->BASELINE[3] = bending;
		myshape->BASELINE[4] = bending;
		myshape->BASELINE[5] = bending;
		break;
	case 3 :
		printf("Slow down!!\n");
		myshape->T = myshape->TOrig*2.0;
		break;
	case 4 :
		printf("Accellerate!!\n");
		myshape->T = myshape->TOrig/2.0;
		break;
	default :
		printf("What da fuck did u map, dude!?\n");
		break;
	}

	// Check consistency between timescales
	learningInterval = myshape->T;
	myshape->mTransitionTime = learningInterval;

	assert(myshape->mTransitionTime==myshape->T);
	assert(myshape->mTransitionTime==learningInterval);
}

void I2D_StefanFishSmartInline::restart(const double t, string filename)
{
	// Restart shape
	FILE * ppFile = NULL;
	assert(filename!=std::string());
	ppFile = fopen(filename.c_str(), "r");
	assert(ppFile!=NULL);
	shape->restart(ppFile);

	// Read target information
	float val = 0.0;
	fscanf(ppFile, "target[0]: %e\n", &val);
	target[0] = val;
	fscanf(ppFile, "target[1]: %e\n", &val);
	target[1] = val;

	if(policy!=NULL)
		if((*policy)!=NULL)
			(*policy)->restart();

	// Cloase file
	fclose(ppFile);
}

void I2D_StefanFishSmartInline::save(const double t, string filename)
{
	FILE * ppFile = NULL;
	assert(filename!=std::string());

	// If string is set I open the corresponding file
	ppFile = fopen(filename.c_str(), "w");
	assert(ppFile!=NULL);

	// Actual restart
	shape->save(ppFile);

	// Add target information
	fprintf(ppFile, "target[0]: %20.20e\n", target[0]);
	fprintf(ppFile, "target[1]: %20.20e\n", target[1]);

	// Cloase file
	fclose(ppFile);
}

