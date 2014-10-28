/*
 * RL_SmartyCircle.cpp
 *
 *  Created on: May 27, 2011
 *      Author: mgazzola
 */

#include <vector>
#include <math.h>
#include <assert.h>
#include <limits>
#include "RL_SmartyCircle.h"

namespace RL
{

RL_SmartyCircle::RL_SmartyCircle(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const int _ID, RL_TabularPolicy ** _policy, const int seed) :
		x(_x), y(_y), D(_D),T(1), maxDomainRadius(0.45), RL_Agent(parser, "RL_SmartyCircle", _ID, _policy, seed)
{
	modv = D/5.0;
	vx = modv;
	vy = 0;

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
}

RL_SmartyCircle::~RL_SmartyCircle()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

void RL_SmartyCircle::mapAction(int action)
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

bool RL_SmartyCircle::mapState(vector<int> & state, map< string, vector<RL_Agent *> > * _data)
{
	assert(_data!=NULL);
	assert((*_data).size()!=0);

	state.clear();

	// Gather my info (always double check with getInfo!!)
	const Real center[2] = {0.5,0.5};
	const Real myX[2] = {x-center[0],y-center[1]};
	const Real myV[2] = {vx,vy};

	// Constants
	const int levelX = signature[0];
	const int levelY = signature[1];
	const int levelAngle = signature[2];
	const Real innerR = 0.1;
	const Real outerR = 0.25;
	bool valid = true;
	const Real cutoff = (outerR-innerR)/2.0;
	const Real h = cutoff / (Real)levelX;
	const Real dangle = 360.0 / (Real)levelAngle;

	Real dist = sqrt( myX[0]*myX[0] + myX[1]*myX[1] ); dist = (dist==0.0)?1.0:dist;
	Real radialDir[2] = {myX[0]/dist,myX[1]/dist};
	radialDir[0] = ((radialDir[0]==0.0) && (radialDir[1]==0.0))?1.0:radialDir[0];
	radialDir[1] = ((radialDir[0]==0.0) && (radialDir[1]==0.0))?0.0:radialDir[1];
	Real innerPoint[2] = {innerR*radialDir[0],innerR*radialDir[1]};
	Real outerPoint[2] = {outerR*radialDir[0],outerR*radialDir[1]};

	Real d[2] = {0.0,0.0};
	if( fabs(dist-innerR)<=fabs(dist-outerR) )
		_dist(innerPoint,myX,d);
	else
		_dist(outerPoint,myX,d);

	const Real IdI = sqrt(d[0]*d[0]+d[1]*d[1]);
	const Real anglev = atan2(myV[1],myV[0]) /M_PI*180.0;
	const Real angled = atan2(d[1],d[0]) /M_PI*180.0;
	Real angle = anglev-angled;
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
		const Real dist = sqrt( myX[0]*myX[0] + myX[1]*myX[1] );
		if( dist > maxDomainRadius )
			valid = false;
	}

	return valid;
}

void RL_SmartyCircle::reward(const Real t, map< string, vector<RL_Agent *> > * _data)
{
	assert(status==Waiting);

	if( (t>(learningTimer+learningInterval)) && status==Waiting)
	{
		Real reward = 0.0;
		const Real rewardRing = 5.0;
		const Real innerR = 0.1;
		const Real outerR = 0.25;
		const Real midR = (outerR+innerR)/2.0;
		const Real xx[2] = {x,y};
		const Real dist = sqrt( (xx[0]-0.5)*(xx[0]-0.5) + (xx[1]-0.5)*(xx[1]-0.5)  );

		const Real lenght = outerR - innerR;
		const Real IdI2 = fabs(dist-midR);
		const Real a = (rewardRing/(lenght*lenght));
		reward += -a*IdI2 + rewardRing;

		integralReward += reward;

		if(integralReward>rewardRing)
		{
			printf("reward=%e\n",reward);
			printf("integralReward=%e\n",integralReward);
		}

		assert(integralReward<=rewardRing);
	}
}

void RL_SmartyCircle::update(const Real dt, const Real t, map< string, vector<RL_Agent *> > * _data, string filename)
{
	x += vx*dt;
	y += vy*dt;

	const Real domain_size1D = 1.0;
	x = (x - domain_size1D*floor(x/domain_size1D));
	y = (y - domain_size1D*floor(y/domain_size1D));

	assert(x>=0 && x<=1);
	assert(y>=0 && y<=1);
}

#ifdef _RL_VIZ
void RL_SmartyCircle::paint()
{
	assert(x>=0 && x<=1);
	assert(y>=0 && y<=1);
	_paintSphere(x,y,D/2.0,0.0,0.0,1.0);
	_drawCircle(maxDomainRadius,0.5,0.5,1,0,0);
	const Real innerR = 0.1;
	const Real outerR = 0.25;
	const Real midR = (outerR+innerR)/2.0;
	_drawCircle(midR,0.5,0.5);
}
#endif

}

