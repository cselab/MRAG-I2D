/*
 * RL_CircularWall.cpp
 *
 *  Created on: Mar 23, 2012
 *      Author: mgazzola
 */
#include <assert.h>
#include "RL_CircularWall.h"

namespace RL
{

RL_CircularWall::RL_CircularWall(MRAG::ArgumentParser & parser) : x(0.5), y(0.5), D(1), innerD(0.95), RL_Agent(parser, "RL_CircularWall")
{
}

RL_CircularWall::~RL_CircularWall()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

#ifdef _RL_VIZ
void RL_CircularWall::paint()
{
	assert(x>=0 && x<=1);
	assert(y>=0 && y<=1);
	_drawFullCircle(D/2.0,0.5,0.5,1,1,1);
}
#endif

} /* namespace RL */
