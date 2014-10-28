/*
 * I2D_StefanFishSmart.cpp
 *
 *  Created on: Feb 8, 2012
 *      Author: mgazzola
 */

#include "I2D_StefanFishSmart.h"
#include <stdio.h>
#include <math.h>

I2D_StefanFishSmart::StefanFishSmart::StefanFishSmart(double xm, double ym, double _D, double T, double tau, double angle_rad, double eps, const int LMAX) :
StefanFish(xm, ym, _D, T, 0.0, tau, angle_rad, vector<double>(6), vector<double>(6), 0.0, eps, LMAX)
{
}

I2D_StefanFishSmart::I2D_StefanFishSmart(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real T, const Real tau, const Real angle, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID, RL::RL_TabularPolicy ** _policy, const int seed) :
		I2D_StefanFish(parser, grid, _xm, _ym, D, T, 0.0, tau, angle, vector<double>(6), vector<double>(6), 0.0, eps, Uinf, penalization, LMAX, ID, _policy, seed), maxDomainRadius(0.45)
{
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
			signature.push_back(5);
			signature.push_back(2);
			signature.push_back(8);
			signature.push_back(3);

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

I2D_StefanFishSmart::~I2D_StefanFishSmart()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

void I2D_StefanFishSmart::mapAction(int action)
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

bool I2D_StefanFishSmart::_mapStateDistance(vector<int> & state, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	assert(_data!=NULL);
	assert((*_data).size()!=0);

	state.clear();

	// Gather my info (always double check with getInfo!!)
	const double center[2] = {0.5,0.5};
	const double myX[2] = {shape->xm-center[0],shape->ym-center[1]};
	const double myV[2] = {shape->vx,shape->vy};

	// Constants
	const int levelX = signature[0];
	const int levelY = signature[1];
	const int levelAngle = signature[2];
	const double innerR = 0.1;
	const double outerR = 0.25;
	bool valid = true;
	const double cutoff = (outerR-innerR)/2.0;
	const double h = cutoff / (double)levelX;
	const double dangle = 360.0 / (double)levelAngle;

	double dist = sqrt( myX[0]*myX[0] + myX[1]*myX[1] ); dist = (dist==0.0)?1.0:dist;
	double radialDir[2] = {myX[0]/dist,myX[1]/dist};
	radialDir[0] = ((radialDir[0]==0.0) && (radialDir[1]==0.0))?1.0:radialDir[0];
	radialDir[1] = ((radialDir[0]==0.0) && (radialDir[1]==0.0))?0.0:radialDir[1];
	double innerPoint[2] = {innerR*radialDir[0],innerR*radialDir[1]};
	double outerPoint[2] = {outerR*radialDir[0],outerR*radialDir[1]};

	double d[2] = {0.0,0.0};
	if( fabs(dist-innerR)<=fabs(dist-outerR) )
		_dist(innerPoint,myX,d);
	else
		_dist(outerPoint,myX,d);

	const double IdI = sqrt(d[0]*d[0]+d[1]*d[1]);
	const double anglev = atan2(myV[1],myV[0]) /M_PI*180.0;
	const double angled = atan2(d[1],d[0]) /M_PI*180.0;
	double angle = anglev-angled;
	angle = (angle<0.0)?angle+360.0:angle;

	const bool inside = (dist>innerR && dist<outerR)?true:false;
	const int idxX = std::max(0,std::min(levelX-1,(int)floor(IdI/h)));
	const int idxY = inside;
	const int idxangle = std::max(0,std::min(levelAngle-1,(int)floor(angle/dangle)));

	// Prepare state vector
	state.push_back(idxX);
	state.push_back(idxY);
	state.push_back(idxangle);

	{
		const double dist = sqrt( myX[0]*myX[0] + myX[1]*myX[1] );
		if( dist > maxDomainRadius )
			valid = false;
	}

	return valid;
}

bool I2D_StefanFishSmart::mapState(vector<int> & state, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	const bool valid = _mapStateDistance(state,_data);

	return valid;
}

void I2D_StefanFishSmart::reward(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	assert(status==Waiting);

	if( t>(learningTimer+learningInterval) && status==Waiting)
	{
		double reward = 0.0;
		const double rewardRing = 5.0;
		const double innerR = 0.1;
		const double outerR = 0.25;
		const double midR = (outerR+innerR)/2.0;
		const double x[2] = {shape->xm,shape->ym};
		const double dist = sqrt( (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5)  );
		const bool inside = (dist>innerR) && (dist<outerR);

		if(inside)
			reward += rewardRing - (rewardRing/fabs(innerR-midR))*fabs(dist-midR);
		else
		{
			if(dist<=innerR)
				reward += -rewardRing + (rewardRing/innerR)*dist;
			else
				reward += -(rewardRing/fabs(maxDomainRadius-outerR))*fabs(dist-outerR);
		}

		integralReward += reward;
	}
}

