/*
 * I2D_StefanFishSmartLattice.cpp
 *
 *  Created on: Feb 29, 2012
 *      Author: mgazzola
 */

#include "I2D_StefanFishSmartLattice.h"

I2D_StefanFishSmartLattice::I2D_StefanFishSmartLattice(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real _T, const Real tau, const Real angle, const Real _dir[2], const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID, RL::RL_TabularPolicy ** _policy, const int seed) :
T(_T), I2D_StefanFishSmart(parser, grid, _xm, _ym, D, _T, tau, angle, eps, Uinf, penalization, LMAX, ID, _policy, seed)
{
	scale = shape->D;
	cutoff = scale;

	const double IdirI = sqrt( _dir[0]*_dir[0] + _dir[1]*_dir[1] );
	assert(IdirI!=0.0);

	dir[0] = _dir[0]/IdirI;
	dir[1] = _dir[1]/IdirI;

	followedPoint[0] = (double)_xm + scale*dir[0];
	followedPoint[1] = (double)_ym + scale*dir[1];

	target[0] = followedPoint[0];
	target[1] = followedPoint[1];

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
			signature.push_back(10);
			signature.push_back(5);
			signature.push_back(5);

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

I2D_StefanFishSmartLattice::~I2D_StefanFishSmartLattice()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

bool I2D_StefanFishSmartLattice::mapState(vector<int> & state, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
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
	const int levelDiffV = signature[2];

	// Get distance
	double d[2] = {0.0,0.0};
	_dist(followedPoint,myX,d);

	// Constants
	const double h = cutoff / (double)levelX;
	const double dangle = 360.0 / (double)levelAngle;

	// Distance
	const double IdI = fabs(sqrt(d[0]*d[0]+d[1]*d[1]) - scale);
	const int idxX = std::max(0,std::min(levelX-1,(int)floor(IdI/h)));

	// Angle
	const double angle = _angleVectors(myV, d);
	const int idxangle = std::max(0,std::min(levelAngle-1,(int)floor(angle/dangle)));

	// Velocity difference
	double totVx = 0.0;
	double totVy = 0.0;
	string object("I2D_StefanFishSmartLattice");
	vector<I2D_FloatingObstacleOperator *> &agents = (*_data)[object];
	for(vector<I2D_FloatingObstacleOperator *>::iterator it=agents.begin(); it!=agents.end(); ++it)
	{
		I2D_StefanFishSmartLattice * b = static_cast<I2D_StefanFishSmartLattice*>(*it);
		totVx += b->shape->vx;
		totVy += b->shape->vy;
	}
	totVx /= (double)agents.size();
	totVy /= (double)agents.size();
	const double modvLattice = totVx*dir[0] + totVy*dir[1];
	const double IvI = sqrt(shape->vx*shape->vx + shape->vy*shape->vy);
	const double diffVelPercent = fabs(IvI-modvLattice)/modvLattice;
	const int idxDiffVel = _discretizeRange(diffVelPercent, 0, 0.1, levelDiffV);

	// Prepare state vector
	state.push_back(idxX);
	state.push_back(idxangle);
	state.push_back(idxDiffVel);

	bool valid = true;
	{
		const double dist = sqrt( (myX[0]-center[0])*(myX[0]-center[0]) + (myX[1]-center[1])*(myX[1]-center[1]) );
		if( dist > maxDomainRadius || IdI>2.0*cutoff )
			valid = false;
	}

	return valid;
}

void I2D_StefanFishSmartLattice::reward(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	assert(status==Waiting);

	if( t>(learningTimer+learningInterval) && status==Waiting)
	{
		double rewardDist = 0;
		double rewardVel = 0;

		{
			const double myX[2] = {shape->xm,shape->ym};
			double d[2] = {0.0,0.0};
			_dist(followedPoint,myX,d);

			const double rewardOnSpot = 1.0;
			const double IdI2 = fabs(d[0]*d[0]+d[1]*d[1] - scale*scale);
			const double a = (rewardOnSpot/(scale*scale));
			rewardDist += -a*IdI2 + rewardOnSpot;
		}

		{
			// Velocity difference
			double totVx = 0.0;
			double totVy = 0.0;
			string object("I2D_StefanFishSmartLattice");
			vector<I2D_FloatingObstacleOperator *> &agents = (*_data)[object];
			for(vector<I2D_FloatingObstacleOperator *>::iterator it=agents.begin(); it!=agents.end(); ++it)
			{
				I2D_StefanFishSmartLattice * b = static_cast<I2D_StefanFishSmartLattice*>(*it);
				totVx += b->shape->vx;
				totVy += b->shape->vy;
			}
			totVx /= (double)agents.size();
			totVy /= (double)agents.size();
			const double modvLattice = totVx*dir[0] + totVy*dir[1];
			const double IvI = sqrt(shape->vx*shape->vx + shape->vy*shape->vy);
			const double diffVelPercent = fabs(IvI-modvLattice)/modvLattice;

			const double rewardGoodVel = 0.1;
			const double a = (rewardGoodVel/(0.1*0.1));
			rewardVel = -a*diffVelPercent + rewardGoodVel;
		}

		integralReward += rewardDist + rewardVel;
	}
}

void I2D_StefanFishSmartLattice::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
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

	// Update follow point
	const Real modv = shape->vx*dir[0] + shape->vy*dir[1];
	followedPoint[0] += modv*dir[0]*dt;
	followedPoint[1] += modv*dir[1]*dt;

	// Update target node
	double totVx = 0.0;
	double totVy = 0.0;
	string object("I2D_StefanFishSmartLattice");
	vector<I2D_FloatingObstacleOperator *> &agents = (*_data)[object];
	for(vector<I2D_FloatingObstacleOperator *>::iterator it=agents.begin(); it!=agents.end(); ++it)
	{
		I2D_StefanFishSmartLattice * b = static_cast<I2D_StefanFishSmartLattice*>(*it);
		totVx += b->shape->vx;
		totVy += b->shape->vy;
	}
	totVx /= (double)agents.size();
	totVy /= (double)agents.size();

	const Real modvLattice = totVx*dir[0] + totVy*dir[1];
	target[0] += modvLattice*dir[0]*dt;
	target[1] += modvLattice*dir[1]*dt;
}

void I2D_StefanFishSmartLattice::mapAction(int action)
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

	const Real base1 = 0.5;
	const Real deltaSpeed = 0.1;

	switch(action)
	{
	case 0 :
		printf("Go straight!!\n");
		myshape->T = T;
		myshape->BASELINE[0] = 0.0;
		myshape->BASELINE[1] = 0.0;
		myshape->BASELINE[2] = 0.0;
		myshape->BASELINE[3] = 0.0;
		myshape->BASELINE[4] = 0.0;
		myshape->BASELINE[5] = 0.0;
		break;
	case 1 :
		printf("Turn right!!\n");
		myshape->T = T;
		myshape->BASELINE[0] = base1;
		myshape->BASELINE[1] = base1;
		myshape->BASELINE[2] = base1;
		myshape->BASELINE[3] = base1;
		myshape->BASELINE[4] = base1;
		myshape->BASELINE[5] = base1;
		break;
	case 2 :
		printf("Turn left!!\n");
		myshape->T = T;
		myshape->BASELINE[0] = -base1;
		myshape->BASELINE[1] = -base1;
		myshape->BASELINE[2] = -base1;
		myshape->BASELINE[3] = -base1;
		myshape->BASELINE[4] = -base1;
		myshape->BASELINE[5] = -base1;
		break;
	case 3 :
		printf("Accelerate!!\n");
		myshape->T = (1-deltaSpeed)*T;
		myshape->BASELINE[0] = 0.0;
		myshape->BASELINE[1] = 0.0;
		myshape->BASELINE[2] = 0.0;
		myshape->BASELINE[3] = 0.0;
		myshape->BASELINE[4] = 0.0;
		myshape->BASELINE[5] = 0.0;
		break;
	case 4 :
		printf("Slow down!!\n");
		myshape->T = (1+deltaSpeed)*T;
		myshape->BASELINE[0] = 0.0;
		myshape->BASELINE[1] = 0.0;
		myshape->BASELINE[2] = 0.0;
		myshape->BASELINE[3] = 0.0;
		myshape->BASELINE[4] = 0.0;
		myshape->BASELINE[5] = 0.0;
		break;
	default :
		printf("What da fuck did u map, dude!?\n");
		break;
	}
}


