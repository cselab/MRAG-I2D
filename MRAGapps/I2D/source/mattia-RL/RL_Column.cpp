/*
 * RL_Column.cpp
 *
 *  Created on: Mar 21, 2012
 *      Author: mgazzola
 */
#include <assert.h>
#include "RL_Column.h"

namespace RL
{

RL_Column::RL_Column(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const int _ID, RL_TabularPolicy ** _policy) : x(_x), y(_y), D(_D), RL_Agent(parser, "RL_Column")
{
}

RL_Column::~RL_Column()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

#ifdef _RL_VIZ
void RL_Column::paint()
{
	assert(x>=0 && x<=1);
	assert(y>=0 && y<=1);
	_paintSphere(x,y,D/2.0,0.752941,0.752941,0.752941);
}
#endif

} /* namespace RL */
