/*
 * RL_SmartyInline.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: mgazzola
 */

#include <vector>
#include <math.h>
#include <assert.h>
#include <limits>
#include "RL_SmartyInline.h"

namespace RL
{

RL_SmartyInline::RL_SmartyInline(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const Real _dir[2], const int _ID, RL_TabularPolicy ** _policy, const int seed) :
										x(_x), y(_y), D(_D),T(1), maxDomainRadius(0.45), RL_Agent(parser, "RL_SmartyInline", _ID, _policy, seed)
{
	modv = D/5.0;

	const Real IdirI = sqrt( _dir[0]*_dir[0] + _dir[1]*_dir[1] );
	assert(IdirI!=0.0);

	dir[0] = _dir[0]/IdirI;
	dir[1] = _dir[1]/IdirI;

	target[0] = x + D*dir[0];
	target[1] = y + D*dir[1];

	vx = modv*dir[0];
	vy = modv*dir[1];

	// Set sensors
	signature.clear();
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			signature.push_back(40);	// ---> dist from lattice point with cutoff
			signature.push_back(20);	// ---> angle between myVel and lattice point
			signature.push_back(3);		// ---> actions

			if(!(*policy)->samedim(signature))
			{
				printf("Policy dimension reset!\n");
				(*policy)->setdim(signature);
			}

			learningInterval = T;
		}
}

RL_SmartyInline::~RL_SmartyInline()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

void RL_SmartyInline::mapAction(int action)
{
	const Real vold[2] = {vx,vy};
	Real v[2] = {vx,vy};
	const Real steeringAngle = M_PI/18.0;
	Real alpha = 0.0;

	switch(action)
	{
	case 0 : // Keep going
		break;
	case 1 : // Rotate 5deg left
		alpha = steeringAngle;
		_rotate(alpha,vold,v);
		break;
	case 2 : // Rotate 5deg right
		alpha = -steeringAngle;
		_rotate(alpha,vold,v);
		break;
	default :
		printf("What da fuck did u map, dude!?\n");
		break;
	}

	//Rescale v to be sure modulus is conserved
	const Real IvI = sqrt(v[0]*v[0]+v[1]*v[1]);
	vx = (v[0]/IvI)*modv;
	vy = (v[1]/IvI)*modv;
}

bool RL_SmartyInline::mapState(vector<int> & state, map< string, vector<RL_Agent *> > * _data)
{
	assert(_data!=NULL);
	assert((*_data).size()!=0);

	state.clear();

	// Gather my info (always double check with getInfo!!)
	const Real center[2] = {0.5,0.5};
	const Real myX[2] = {x,y};
	const Real myV[2] = {vx,vy};

	// Levels
	const int levelX = signature[0];
	const int levelAngle = signature[1];

	// Constants
	const Real cutoff = 2*D;
	const Real h = cutoff / (Real)levelX;
	const Real dangle = 360.0 / (Real)levelAngle;

	Real d[2] = {0.0,0.0};
	_dist(target,myX,d);

	const Real angle = _angleVectors(myV, d);

	const Real IdI = sqrt(d[0]*d[0]+d[1]*d[1]) - D;

	const int idxX = std::max(0,std::min(levelX-1,(int)floor(IdI/h)));
	const int idxangle = std::max(0,std::min(levelAngle-1,(int)floor(angle/dangle)));

	// Prepare state vector
	state.push_back(idxX);
	state.push_back(idxangle);

	bool valid = true;
	{
		const Real dist = sqrt( (myX[0]-center[0])*(myX[0]-center[0]) + (myX[1]-center[1])*(myX[1]-center[1]) );
		if( dist > maxDomainRadius || IdI>(2*cutoff) )
			valid = false;
	}

	return valid;
}

void RL_SmartyInline::reward(const Real t, map< string, vector<RL_Agent *> > * _data)
{
	assert(status==Waiting);

	if( t>(learningTimer+learningInterval) && status==Waiting)
	{
		Real reward = 0.0;

		// Gather my info (always double check with getInfo!!)
		const Real myX[2] = {x,y};

		// Constants
		const Real cutoff = D;
		const Real rewardOnSpot = 5.0;

		Real d[2] = {0.0,0.0};
		_dist(target,myX,d);

		const Real IdI2 = d[0]*d[0]+d[1]*d[1] - D*D;
		const Real a = (rewardOnSpot/(cutoff*cutoff));
		reward += -a*IdI2 + rewardOnSpot;

		integralReward += reward;
	}
}

void RL_SmartyInline::update(const Real dt, const Real t, map< string, vector<RL_Agent *> > * _data, string filename)
{
	x += vx*dt;
	y += vy*dt;

	// Update lattice node
	target[0] += modv*dir[0]*dt;
	target[1] += modv*dir[1]*dt;

	assert(x>=0 && x<=1);
	assert(y>=0 && y<=1);
}

#ifdef _RL_VIZ
void RL_SmartyInline::paint()
{
	assert(x>=0 && x<=1);
	assert(y>=0 && y<=1);
	_paintSphere(x,y,D/2.0,0.0,0.0,1.0);
	_paintSphere(target[0],target[1],D/5.0);
}
#endif

} /* namespace RL */
