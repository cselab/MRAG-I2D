/*
 * RL_Food.cpp
 *
 *  Created on: Mar 27, 2012
 *      Author: dalmassg
 */

#include "rng.h"
#include <assert.h>
#include "RL_Food.h"

namespace RL
{

RL_Food::RL_Food(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const int _ID, const int seed, RL_TabularPolicy ** _policy) : x(_x), y(_y), D(_D), count(1), rng(seed), RL_Agent(parser, "RL_Food", _ID, NULL, seed)
{
	//printf("seed = %f\n", rng.uniform());
	//exit(0);
}

RL_Food::~RL_Food()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

void RL_Food::update(const double dt, const double t, map< string, vector<RL_Agent *> > * _data, string filename)
{
	//const Real stepLength = D*dt;
	//const Real angle = rng.uniform( 0.0,2.0*M_PI );
	//x += stepLength*cos(angle);
	//y += stepLength*sin(angle);

	if(count%500000 == 0)
	{
		//const Real angleJump = rng.uniform(0.0, 2.0*M_PI);
		//const Real radiusJump = rng.uniform(0.0, 0.5);
		//x = 0.5 + radiusJump*cos(angleJump);
		//y = 0.5 + radiusJump*sin(angleJump);
		x += rng.uniform( -0.5, 0.5 );
		y += rng.uniform( -0.5, 0.5 );
	}

	// Periodic domain reprojection
	const Real domain_size1D = 1.0;
	x = (x - domain_size1D*floor(x/domain_size1D));
	y = (y - domain_size1D*floor(y/domain_size1D));

	// Arena reprojection
	//const Real rad = sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) );
	//x = (rad>0.5)?0.5:x;
	//y = (rad>0.5)?0.5:y;

	assert(x>=0 && x<=1);
	assert(y>=0 && y<=1);
	count++;
}

#ifdef _RL_VIZ
void RL_Food::paint()
{
	assert(x>=0 && x<=1);
	assert(y>=0 && y<=1);
	_paintSphere(x,y,D/2.0,0.8,0.5,0.0);
}
#endif

} /* namespace RL */
