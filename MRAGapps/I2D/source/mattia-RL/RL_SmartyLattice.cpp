/*
 * RL_SmartyLattice.cpp
 *
 *  Created on: Apr 5, 2012
 *      Author: mgazzola
 */
#include <vector>
#include <math.h>
#include <assert.h>
#include <limits>
#include <fstream>
#include "RL_SmartyLattice.h"

namespace RL
{

RL_SmartyLattice::RL_SmartyLattice(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const Real _dir[2], const int _ID, RL_TabularPolicy ** _policy, const int seed) :
																x(_x), y(_y), D(_D),T(1.0), maxDomainRadius(0.45), rng(seed), RL_Agent(parser, "RL_SmartyLattice", _ID, _policy, seed)
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
			signature.push_back(20);	// ---> dist from lattice point with cutoff
			signature.push_back(20);	// ---> angle between myVel and lattice point
			signature.push_back(5);		// ---> actions

			if(!(*policy)->samedim(signature))
			{
				printf("Policy dimension reset!\n");
				(*policy)->setdim(signature);
			}

			learningInterval = T;
		}
}

RL_SmartyLattice::~RL_SmartyLattice()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

void RL_SmartyLattice::mapAction(int action)
{
	const Real vold[2] = {vx,vy};
	Real v[2] = {vx,vy};
	const Real steeringAngle = M_PI/18.0;
	Real alpha = 0.0;
	Real referenceV = modv;

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
	case 3 : // Accelerate
		referenceV *= 2.0;
		break;
	case 4 : // Decelerate
		referenceV *= 0.5;
		break;
	default :
		printf("What da fuck did u map, dude!?\n");
		break;
	}

	//Rescale v to be sure modulus is conserved
	const Real IvI = sqrt(v[0]*v[0]+v[1]*v[1]);
	vx = (v[0]/IvI)*referenceV;
	vy = (v[1]/IvI)*referenceV;
	const Real fluctuation = 0.025*modv;
	vx += rng.uniform(-fluctuation,fluctuation);
	vy += rng.uniform(-fluctuation,fluctuation);
}

bool RL_SmartyLattice::mapState(vector<int> & state, map< string, vector<RL_Agent *> > * _data)
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
	const Real cutoff = D;
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
		if( dist > maxDomainRadius )
			valid = false;
	}

	return valid;
}

void RL_SmartyLattice::reward(const Real t, map< string, vector<RL_Agent *> > * _data)
{
	assert(status==Waiting);

	if( t>(learningTimer+learningInterval) && status==Waiting)
	{
		Real reward = 0.0;

		// Gather my info (always double check with getInfo!!)
		const Real myX[2] = {x,y};

		// Constants
		const Real cutoff = D;
		//const Real rewardOnSpot = 0.0;

		Real d[2] = {0.0,0.0};
		_dist(target,myX,d);

		const Real IdI2 = fabs(d[0]*d[0]+d[1]*d[1] - D*D);
		reward += -(1/cutoff*cutoff)*IdI2;

		integralReward += reward;
	}
}

void RL_SmartyLattice::update(const Real dt, const Real t, map< string, vector<RL_Agent *> > * _data, string filename)
{
	x += vx*dt;
	y += vy*dt;

	assert(x>=0 && x<=1);
	assert(y>=0 && y<=1);

	// Update lattice node
	string object("RL_SmartyLattice");
	vector<RL_Agent *> &agents = (*_data)[object];

	// Compute average velocity among all lattice agents
	double totVx = 0.0;
	double totVy = 0.0;
	for(vector<RL_Agent *>::iterator it=agents.begin(); it!=agents.end(); ++it)
	{
		RL_SmartyLattice * b = static_cast<RL_SmartyLattice*>(*it);
		totVx += b->vx;
		totVy += b->vy;
	}
	totVx /= (double)agents.size();
	totVy /= (double)agents.size();

	const double modvLattice = totVx*dir[0] + totVy*dir[1];
	target[0] += modvLattice*dir[0]*dt;
	target[1] += modvLattice*dir[1]*dt;

	// Dump
	//char buf[100];
	//sprintf(buf, "stats_%04d", ID);
	//string namefile(buf);
	//ofstream out(namefile.c_str(), ios_base::app);
	//out << t << "  " << x << "  " << y << "  " << vx << "  " << vy << "  " << target[0] << "  " << target[1] << "  " << modvLattice << endl;
	//out.flush();
	//out.close();
}

#ifdef _RL_VIZ
void RL_SmartyLattice::paint()
{
	assert(x>=0 && x<=1);
	assert(y>=0 && y<=1);
	_paintSphere(x,y,D/2.0,1.0,0.0,0.0);
	_paintSphere(target[0],target[1],D/5.0);
}
#endif

} /* namespace RL */
