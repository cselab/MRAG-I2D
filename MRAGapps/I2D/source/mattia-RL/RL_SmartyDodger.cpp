/*
 * RL_SmartyDodger.cpp
 *
 *  Created on: Mar 22, 2012
 *      Author: mgazzola
 */
#include <vector>
#include <math.h>
#include <assert.h>
#include <limits>
#include "RL_SmartyDodger.h"
#include "RL_CircularWall.h"
#include "RL_Column.h"
#include "RL_DynamicColumn.h"

namespace RL {

RL_SmartyDodger::RL_SmartyDodger(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const Real _T, const int _ID, RL_TabularPolicy ** _policy, const int seed) : x(_x), y(_y), D(_D),T(_T), RL_Agent(parser, "RL_SmartyDodger", _ID, _policy, seed)
{
	modv = D/5.0;
	vx = modv;
	vy = 0;

	// Set sensors
	signature.clear();
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			signature.push_back(10);
			signature.push_back(10);
			signature.push_back(10);
			signature.push_back(10);
			signature.push_back(10);
			signature.push_back(10);
			signature.push_back(3);

			if(!(*policy)->samedim(signature))
			{
				printf("Policy dimension reset!\n");
				(*policy)->setdim(signature);
			}

			learningInterval = T;
		}
}

RL_SmartyDodger::~RL_SmartyDodger()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

void RL_SmartyDodger::_mapForMaxRadius(const Real cutoff, const int idxDistLevel, const int idxAngleLevel, int & idxDist, int & idxAngle, map< string, vector<RL_Agent *> > * _data)
{
	// Gather my info (always double check with getInfo!!)
	const Real center[2] = {0.5,0.5};
	const Real myX[2] = {x,y};
	const Real myV[2] = {vx,vy};
	const Real radius = sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) );

	// Get maxDomainRadius
	string object("RL_CircularWall");
	vector<RL_Agent *> &wall = (*_data)[object];
	assert(wall.size()==1);
	RL_CircularWall * b = static_cast<RL_CircularWall*>(*wall.begin());
	assert(b!=NULL);
	const Real maxDomainRadius = b->D/2;

	// Get target point position
	Real dist[2];
	_dist(myX, center, dist);
	Real dir[2];
	const Real IdistI = _normalize(dist,dir);
	const Real target[2] = { maxDomainRadius*dir[0] + center[0], maxDomainRadius*dir[1] + center[0] };

	// Discretize distance from target
	Real d[2];
	_dist(target, myX, d);
	const Real IdI = _modv(d);
	const bool deadEnd = IdI>=cutoff;
	idxDist = _discretize1SS(IdI, 0, cutoff, idxDistLevel, deadEnd);

	// Discretize distance to target
	const Real angle = _angleVectors(d, myV);
	idxAngle = _discretize1SS(angle, 0, 360, idxAngleLevel, deadEnd);
}

void RL_SmartyDodger::_mapForColumns(const Real cutoff, const int idxDistLevel, const int idxAngleLevel, int & idxDist, int & idxAngle, map< string, vector<RL_Agent *> > * _data)
{
	// Gather my info (always double check with getInfo!!)
	const Real myX[2] = {x,y};
	const Real myV[2] = {vx,vy};

	// Get columns
	string object("RL_Column");
	vector<RL_Agent *> &columns = (*_data)[object];

	// Find the closest
	RL_Column * closestAgent = NULL;
	Real minDist = numeric_limits<Real>::max();
	for(vector<RL_Agent *>::iterator it=columns.begin(); it!=columns.end(); ++it)
	{
		RL_Column * b = static_cast<RL_Column*>(*it);
		const Real target[2] = {b->x,b->y};
		Real d[2];
		_dist(target, myX, d);
		const Real IdI = max(0.0,_modv(d) - (b->D/2.0));
		closestAgent = (IdI<minDist)?b:closestAgent;
		minDist = (IdI<minDist)?IdI:minDist;
	}

	if(closestAgent!=NULL)
	{
		// Discretize distance from target
		const Real target[2] = {closestAgent->x,closestAgent->y};
		Real d[2];
		_dist(target, myX, d);
		const Real signIdI = _modv(d) - (closestAgent->D/2.0);
		const Real IdI = fabs(signIdI);
		const bool deadEnd = IdI>=cutoff;
		const bool inside = signIdI<=0;
		idxDist = _discretize2SS(IdI, 0, cutoff, idxDistLevel, deadEnd,inside);

		// Discretize distance to target
		const Real angle = _angleVectors(d, myV);
		idxAngle = _discretize2SS(angle, 0, 360, idxAngleLevel, deadEnd,inside);
	}
}

void RL_SmartyDodger::_mapForDynamicColumns(const Real cutoff, const int idxDistLevel, const int idxAngleLevel, int & idxDist, int & idxAngle, map< string, vector<RL_Agent *> > * _data)
{
	// Gather my info (always double check with getInfo!!)
	const Real myX[2] = {x,y};
	const Real myV[2] = {vx,vy};

	// Get dynamic columns
	string object("RL_DynamicColumn");
	vector<RL_Agent *> &Dcolumns = (*_data)[object];

	// Find the closest
	RL_DynamicColumn * closestAgent = NULL;
	Real minDist = numeric_limits<Real>::max();
	for(vector<RL_Agent *>::iterator it=Dcolumns.begin(); it!=Dcolumns.end(); ++it)
	{
		RL_DynamicColumn * b = static_cast<RL_DynamicColumn*>(*it);
		const Real target[2] = {b->x,b->y};
		Real d[2];
		_dist(target, myX, d);
		const Real IdI = max(0.0,_modv(d) - (b->D/2.0));
		closestAgent = (IdI<minDist)?b:closestAgent;
		minDist = (IdI<minDist)?IdI:minDist;
	}

	if(closestAgent!=NULL)
	{
		// Discretize distance from target
		const Real target[2] = {closestAgent->x,closestAgent->y};
		Real d[2];
		_dist(target, myX, d);
		const Real signIdI = _modv(d) - (closestAgent->D/2.0);
		const Real IdI = fabs(signIdI);
		const bool deadEnd = IdI>=cutoff;
		const bool inside = signIdI<=0;
		idxDist = _discretize2SS(IdI, 0, cutoff, idxDistLevel, deadEnd,inside);

		// Discretize distance to target
		const Real angle = _angleVectors(d, myV);
		idxAngle = _discretize2SS(angle, 0, 360, idxAngleLevel, deadEnd,inside);
	}
}

bool RL_SmartyDodger::mapState(vector<int> & state, map< string, vector<RL_Agent *> > * _data)
{
	assert(_data!=NULL);
	assert((*_data).size()!=0);

	state.clear();

	const Real cutoffWall = 5*D;
	const Real cutoffColumns = 5*D;
	const Real cutoffDynamicColumns = 5*D;

	// Map wall
	{
		int idxDistMaxRadius = 0;
		int idxAngleMaxRadius = 0;
		_mapForMaxRadius(cutoffWall, signature[0], signature[1], idxDistMaxRadius, idxAngleMaxRadius, _data);
		state.push_back(idxDistMaxRadius);
		state.push_back(idxAngleMaxRadius);
	}

	// Map columns
	{
		int idxDistMaxColumns = 0;
		int idxAngleMaxColumns = 0;
		_mapForColumns(cutoffColumns, signature[2], signature[3], idxDistMaxColumns, idxAngleMaxColumns, _data);
		state.push_back(idxDistMaxColumns);
		state.push_back(idxAngleMaxColumns);
	}

	// Map dynamic columns
	{
		int idxDistMaxDynamicColumns = 0;
		int idxAngleMaxDynamicColumns = 0;
		_mapForDynamicColumns(cutoffDynamicColumns, signature[4], signature[5], idxDistMaxDynamicColumns, idxAngleMaxDynamicColumns, _data);
		state.push_back(idxDistMaxDynamicColumns);
		state.push_back(idxAngleMaxDynamicColumns);
	}

	// Check validity
	bool valid = true;
	//if( sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ) > 0.45 )
	//	valid = false;

	return valid;
}

void RL_SmartyDodger::mapAction(int action)
{
	const Real vold[2] = {vx,vy};
	Real v[2] = {vx,vy};
	const Real steeringAngle = 5.0/180.0*M_PI;
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

void RL_SmartyDodger::reward(const double t, map< string, vector<RL_Agent *> > * _data)
{
	assert(status==Waiting);

	if( (t>(learningTimer+learningInterval)) && status==Waiting)
	{
		Real reward = 0.0;
		const Real myX[2] = {x,y};

		// Wall
		{
			string object("RL_CircularWall");
			vector<RL_Agent *> &wall = (*_data)[object];
			assert(wall.size()==1);
			RL_CircularWall * b = static_cast<RL_CircularWall*>(*wall.begin());
			assert(b!=NULL);
			const Real dist = max(0.0,sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)) - b->innerD/2);
			const Real localReward = -1.0/(b->D/2.0 - b->innerD/2)*dist;
			reward += localReward;
		}

		// Columns
		{
			string object("RL_Column");
			vector<RL_Agent *> &columns = (*_data)[object];

			// Find the closest
			for(vector<RL_Agent *>::iterator it=columns.begin(); it!=columns.end(); ++it)
			{
				RL_Column * b = static_cast<RL_Column*>(*it);
				const Real target[2] = {b->x,b->y};
				Real d[2];
				_dist(target, myX, d);
				const Real IdI = _modv(d) - (b->D/2.0);
				const Real localReward = (IdI<=0.0)?-1:0;
				reward += localReward;
				if(localReward==-1)
					break;
			}
		}

		// DynamicColumns
		{
			string object("RL_DynamicColumn");
			vector<RL_Agent *> &Dcolumns = (*_data)[object];

			// Find the closest
			for(vector<RL_Agent *>::iterator it=Dcolumns.begin(); it!=Dcolumns.end(); ++it)
			{
				RL_DynamicColumn * b = static_cast<RL_DynamicColumn*>(*it);
				const Real target[2] = {b->x,b->y};
				Real d[2];
				_dist(target, myX, d);
				const Real IdI = _modv(d) - (b->D/2.0);
				const Real localReward = (IdI<=0.0)?-1:0;
				reward += localReward;
				if(localReward==-1)
					break;
			}
		}

		//integralReward += (reward==0)?1:reward;
		integralReward += reward;
	}
}

void RL_SmartyDodger::update(const double dt, const double t, map< string, vector<RL_Agent *> > * _data, string filename)
{
	Real const xold = x;
	Real const yold = y;
	x += vx*dt;
	y += vy*dt;
	const Real dist = sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) );
	const bool inside = (dist<0.5);

	x = inside?x:xold;
	vx = inside?vx:-vx;
	y = inside?y:yold;
	vy = inside?vy:-vy;

	//const Real domain_size1D = 1.0;
	//x = (x - domain_size1D*floor(x/domain_size1D));
	//y = (y - domain_size1D*floor(y/domain_size1D));
	//assert(x>=0 && x<=1);
	//assert(y>=0 && y<=1);
}

#ifdef _RL_VIZ
void RL_SmartyDodger::paint()
{
	assert(x>=0 && x<=1);
	assert(y>=0 && y<=1);
	//_drawCircle(D/2.0,x,y);
	_paintSphere(x,y,D/2.0,0.0,1.0,0.0);
}
#endif

} /* namespace RL */
