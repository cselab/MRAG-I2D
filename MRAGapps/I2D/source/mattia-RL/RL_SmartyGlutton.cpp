/*
 * RL_SmartyGlutton.cpp
 *
 *  Created on: Mar 27, 2012
 *      Author: dalmassg
 */

#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <limits>
#include <iostream>
#include <fstream>
#include "RL_SmartyGlutton.h"
#include "RL_CircularWall.h"
#include "RL_Food.h"
#include "rng.h"

namespace RL {

RL_SmartyGlutton::RL_SmartyGlutton(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const Real _T, const int _nactions, const Real _cutoffGluttons, const Real _cutoffFood, const Real _steer, const Real _noise, const Real _selfAvoid, const int NN, const int _ID, RL_TabularPolicy ** _policy, const int seed):
																				x(_x), y(_y), D(_D),T(_T), nactions(_nactions), cutoffGluttons(_cutoffGluttons), cutoffFood(_cutoffFood), steer(_steer), noise(_noise), selfAvoid(_selfAvoid), rng(seed), myfoodcounter(0), NN(NN), STARVING(true), RL_Agent(parser, "RL_SmartyGlutton", _ID, _policy, seed)
{

	//printf("seed = %f\n", rng.uniform());
	//exit(0);

	if(nactions!=3 && nactions!=4)
	{
		printf("Wrong number of actions(%d)!\n",nactions);
		abort();
	}

	IDcollision = this->ID;

	const Real angleStart = rng.uniform(0.0,2*M_PI);

	//modv = D/5.0;
	modv = D;
	vx = modv*cos(angleStart);
	vy = modv*sin(angleStart);
	Vm[0] = 0.0;
	Vm[1] = 0.0;

	// Set sensors
	signature.clear();
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			signature.push_back(3);		// distance food
			signature.push_back(5);		// angle food
			signature.push_back(6);		// distance gluttons
			signature.push_back(6);		// angle gluttons
			signature.push_back(6);		// angle relative motion gluttons
			signature.push_back(nactions);	// actions

			if(!(*policy)->samedim(signature))
			{
				printf("Policy dimension reset!\n");
				(*policy)->setdim(signature);
			}

			learningInterval = T;
		}
}

RL_SmartyGlutton::~RL_SmartyGlutton()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

void RL_SmartyGlutton::_rotateSimple(const Real theta, Real x[2]) const
{
	const Real a00 = cos(theta);
	const Real a01 = -sin(theta);
	const Real a10 = sin(theta);
	const Real a11 = cos(theta);

	const Real xx = x[0];
	const Real yy = x[1];

	x[0] = a00*xx + a01*yy;
	x[1] = a10*xx + a11*yy;
}

void RL_SmartyGlutton::_mapForFoods(const Real cutoff, const int idxDistLevel, const int idxAngleLevel, int & idxDist, int & idxAngle, map< string, vector<RL_Agent *> > * _data)
{
	// Gather my info (always double check with getInfo!!)
	const Real myX[2] = {x,y};
	const Real myV[2] = {vx,vy};

	// Get foods
	string object("RL_Food");
	vector<RL_Agent *> &foods = (*_data)[object];

	// Find the closest
	RL_Food * closestAgent = NULL;
	Real minDist = numeric_limits<Real>::max();
	for(vector<RL_Agent *>::iterator it=foods.begin(); it!=foods.end(); ++it)
	{
		RL_Food * b = static_cast<RL_Food*>(*it);
		const Real target[2] = {b->x,b->y};
		Real d[2];
		_distPeriodic(target, myX, d);
		const Real IdI = _modv(d);
		closestAgent = (IdI<minDist)?b:closestAgent;
		minDist = (IdI<minDist)?IdI:minDist;
	}

	if(closestAgent!=NULL)
	{
		// Discretize distance from target
		const Real target[2] = {closestAgent->x,closestAgent->y};
		Real d[2];
		_distPeriodic(target, myX, d);
		const Real IdI = _modv(d) - (closestAgent->D/2.0) -D/2;
		const bool deadEnd = IdI>=cutoff;
		idxDist = _discretize1SS(IdI, -D, cutoff, idxDistLevel, deadEnd);

		// Discretize distance to target
		const Real angle = _angleVectors(d, myV);
		idxAngle = _discretize1SS(angle, 0, 360, idxAngleLevel, deadEnd);
	}
}

void RL_SmartyGlutton::_mapForGluttons(const Real cutoff, const int idxDistLevel, const int idxAngleLevel, const int idxAngleNeighLevel, RL_SmartyGlutton * neigh, Real vneigh[2], int & idxDist, int & idxAngle, int & idxAngleNeigh, map< string, vector<RL_Agent *> > * _data)
{
	// Gather my info (always double check with getInfo!!)
	const Real myX[2] = {x,y};
	const Real myV[2] = {vx,vy};
	Real IdI = 0.0;
	Real d[2] = {0.0,0.0};

	RL_SmartyGlutton * closestAgent = neigh;
	const bool deadEnd = (closestAgent==NULL)?true:false;
	if(!deadEnd)
	{
		const Real target[2] = {closestAgent->x,closestAgent->y};
		_distPeriodic(target, myX, d);
		IdI = _modv(d);
	}

	// Discretize distance to closest neighbor
	idxDist = _discretize1SS(IdI, 0, cutoff, idxDistLevel, deadEnd);

	// Discretize angle between distance and my velocity
	const Real angle = _angleVectors(d, myV);
	idxAngle = _discretize1SS(angle, 0, 360, idxAngleLevel, deadEnd);

	// Discretize angle between closest neighbor's velocity and my velocity
	const Real angleNeigh = _angleVectors(vneigh, myV);
	idxAngleNeigh = _discretize1SS(angleNeigh, 0, 360, idxAngleNeighLevel, deadEnd);

}

bool RL_SmartyGlutton::mapState(vector<int> & state, map< string, vector<RL_Agent *> > * _data)
{
	assert(_data!=NULL);
	assert((*_data).size()!=0);

	state.clear();

	const Real cutFood = cutoffFood*D;
	const Real cutGluttons = cutoffGluttons*D;

	RL_SmartyGlutton * neigh = NULL;
	Real vneigh[2] = {0.0,0.0};

	//	Neighborhood velocity and velocity of the closest neighbor (if any) within cutoff
	{
		const Real myX[2] = {x,y};
		Real myVm[2] = {0.0,0.0};
		unsigned int counter1 = 0;

		// Get gluttons
		string object("RL_SmartyGlutton");
		vector<RL_Agent *> &gluttons = (*_data)[object];
		RL_SmartyGlutton * closestAgent = NULL;
		assert((*_data)[object].size() != 0);

		// Collect distances form other gluttons within cutoff
		map<Real,RL_SmartyGlutton *> distToV;
		vector<Real> distances;
		for(vector<RL_Agent *>::iterator it=gluttons.begin(); it!=gluttons.end(); ++it)
		{
			RL_SmartyGlutton * b = static_cast<RL_SmartyGlutton*>(*it);
			const Real target[2] = {b->x,b->y};
			Real d[2];
			_distPeriodic(target, myX, d);
			const Real IdI = _modv(d);

			if(IdI<=cutGluttons && b->ID!=this->ID)
			{
				distances.push_back(IdI);
				distToV[IdI]=b;
				counter1++;
			}
		}

		//Compute vneigh
		sort(distances.begin(), distances.end());
		if(counter1!=0)
		{
			neigh = distToV[ distances[0] ];
			vneigh[0] += neigh->vx;
			vneigh[1] += neigh->vy;
		}

		// Compute vel of neighborhood next to me
		if(nactions==4)
		{
			myVm[0] = vx;
			myVm[1] = vy;
			Real totalWeight = 1.0;
			const int nnInteraction = (NN==-1)?distances.size():NN;
			const int maxNN = (distances.size()<nnInteraction)?distances.size():nnInteraction;
			for(int nn=0; nn<maxNN; nn++)
			{
				RL_SmartyGlutton * b = distToV[ distances[nn] ];
				const Real weight = 1.0;
				totalWeight += weight;
				myVm[0] += weight*b->vx;
				myVm[1] += weight*b->vy;
			}
			myVm[0] /= totalWeight;
			myVm[1] /= totalWeight;
			const Real ImyVmI1 = sqrt(myVm[0]*myVm[0]+myVm[1]*myVm[1]);
			Vm[0] = (myVm[0]/ImyVmI1)*modv;
			Vm[1] = (myVm[1]/ImyVmI1)*modv;

			//Vm[0] = ( counter1 != 0 && ImyVmI1 != 0.0 && maxNN!=0)?(myVm[0]/ImyVmI1)*modv:vx;
			//Vm[1] = ( counter1 != 0 && ImyVmI1 != 0.0 && maxNN!=0)?(myVm[1]/ImyVmI1)*modv:vy;
		}
	}

	// Map foods
	{
		int idxDistMaxFoods = 0;
		int idxAngleMaxFoods = 0;
		_mapForFoods(cutFood, signature[0], signature[1], idxDistMaxFoods, idxAngleMaxFoods, _data);
		state.push_back(idxDistMaxFoods);
		state.push_back(idxAngleMaxFoods);
	}

	// Map gluttons
	{
		int idxDistMaxGluttons = 0;
		int idxAngleMaxGluttons = 0;
		int idxAngleNeigh = 0;
		_mapForGluttons(cutGluttons, signature[2], signature[3], signature[4], neigh, vneigh, idxDistMaxGluttons, idxAngleMaxGluttons, idxAngleNeigh, _data);
		state.push_back(idxDistMaxGluttons);
		state.push_back(idxAngleMaxGluttons);
		state.push_back(idxAngleNeigh);
	}



	// Check validity
	bool valid = true;

	return valid;
}

void RL_SmartyGlutton::mapAction(int action)
{
	myaction = action;

	const Real vold[2] = {vx,vy};
	Real v[2] = {vx,vy};
	const Real steeringAngle = steer/180.0*M_PI;
	Real alpha = 0.0;

	switch(myaction)
	{
	case 0 : 	// Keep going
		break;
	case 1 : 	// Rotate 5deg left
		alpha = steeringAngle;
		_rotate(alpha,vold,v);
		break;
	case 2 : 	// Rotate 5deg right
		alpha = -steeringAngle;
		_rotate(alpha,vold,v);
		break;
	case 3 : 	// Follow mean velocity
		v[0] = Vm[0];
		v[1] = Vm[1];
		break;
	default :
		printf("What da fuck did u map, dude!?\n");
		break;
	}

	_rotateSimple(rng.normal(0.0, noise),v);

	//Rescale v to be sure modulus is conserved
	const Real IvI = sqrt(v[0]*v[0]+v[1]*v[1]);
	vx = (v[0]/IvI)*modv;
	vy = (v[1]/IvI)*modv;
}

void RL_SmartyGlutton::reward(const double t, map< string, vector<RL_Agent *> > * _data)
{
	assert(status==Waiting);

	if( (t>(learningTimer+learningInterval)) && status==Waiting)
	{
		Real reward = 0.0;
		const Real rewardFood = 1.0;
		const Real rewardGluttons = selfAvoid;
		const Real myX[2] = {x,y};

		// Gluttons
		{
			const int IDcollisionBefore = IDcollision;

			string object("RL_SmartyGlutton");
			vector<RL_Agent *> &gluttons = (*_data)[object];
			RL_SmartyGlutton * closestAgent = NULL;

			// Find the closest
			for(vector<RL_Agent *>::iterator it=gluttons.begin(); it!=gluttons.end(); ++it)
			{
				RL_SmartyGlutton * b = static_cast<RL_SmartyGlutton*>(*it);
				const Real target[2] = {b->x,b->y};
				const int otherID = b->ID;
				Real d[2];
				_distPeriodic(target, myX, d);
				const Real IdI = _modv(d) - (b->D/2.0) - D/2;
				const Real localReward = (IdI<=0.0 && otherID!=this->ID)?rewardGluttons:0;
				reward += localReward;

				if(localReward==rewardGluttons)
					break;
			}
		}

		// Food
		{
			string object("RL_Food");
			map< string, vector<RL_Agent *> >::iterator itmap = _data->find(object);

			const bool STARVINGBEFORE = STARVING;

			if( itmap!=(*_data).end() )
			{
				vector<RL_Agent *> &foods = (*_data)[object];
				RL_Food * closestAgent = NULL;

				// Find the closest
				for(vector<RL_Agent *>::iterator it=foods.begin(); it!=foods.end(); ++it)
				{
					RL_Food * b = static_cast<RL_Food*>(*it);
					const Real target[2] = {b->x,b->y};
					Real d[2];
					_distPeriodic(target, myX, d);
					const Real IdI = _modv(d) - (b->D/2.0);
					const Real localReward = (IdI<=0.0)?rewardFood:0.0;
					reward += localReward;
					if(localReward==rewardFood)
						break;
				}
			}
		}

		integralReward += reward;
	}
}

void RL_SmartyGlutton::update(const double dt, const double t, map< string, vector<RL_Agent *> > * _data, string filename)
{
	// Update positions
	x += vx*dt;
	y += vy*dt;

	// Correct for periodic domain [0,1]x[0,1]
	const Real domain_size1D = 1.0;
	x = (x - domain_size1D*floor(x/domain_size1D));
	y = (y - domain_size1D*floor(y/domain_size1D));

	//const Real xold = x;
	//const Real yold = y;
	//x += vx*dt;
	//y += vy*dt;
	//const Real dist = sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) );
	//const bool inside = (dist<0.5);
	//x = inside?x:xold;
	//vx = inside?vx:-vx;
	//y = inside?y:yold;
	//vy = inside?vy:-vy;
}

void RL_SmartyGlutton::dump(const double t, ofstream & out)
{
	out.precision(6);
	out << scientific << t << "\t" << this->x << "\t" << this->y << "\t" << this->vx << "\t" << this->vy << "\t";

	assert(myStateAction.size()!=0);

	for(unsigned int i=0; i<myStateAction.size(); i++)
		out << myStateAction[i] << "\t";

	out << "\n";
	out.flush();
}

#ifdef _RL_VIZ
void RL_SmartyGlutton::paint()
{
	//assert(x>=0 && x<=1);
	//assert(y>=0 && y<=1);
	_paintSphere(x,y,(3.0/2.0)*(D/2.0),0.0,0.741,0.965);
	//_paintSphere(x,y,D/2.0,1.0,0.0,0.0);
}
#endif

} /* namespace RL */

