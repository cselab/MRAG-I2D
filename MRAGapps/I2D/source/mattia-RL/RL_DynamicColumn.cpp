/*
 * RL_DynamicColumn.cpp
 *
 *  Created on: Apr 5, 2012
 *      Author: dalmassg
 */

#include "RL_DynamicColumn.h"

namespace RL
{

RL_DynamicColumn::RL_DynamicColumn(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const int _ID, RL_TabularPolicy ** _policy) : x(_x), y(_y), D(_D), RL_Agent(parser, "RL_Column")
{
}

RL_DynamicColumn::~RL_DynamicColumn()
{
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

void RL_DynamicColumn::update(const Real dt, const Real t, map< string, vector<RL_Agent *> > * _data, string filename)
{
	const double angleOsc = 2.0/180.0*M_PI;
	const Real xC = x-0.5-sin(angleOsc);
	const Real yC = y-0.5-sin(angleOsc);

//	if(count >= 50000)
//	{
		const double angleRot = 90.0/180.0*M_PI;
		x += ( xC*cos(angleRot)-yC*sin(angleRot) )*0.005*dt;
		y += ( xC*sin(angleRot)+yC*cos(angleRot) )*0.005*dt;
//	}
//	count++;
}

#ifdef _RL_VIZ
void RL_DynamicColumn::paint()
{
	_paintSphere(x,y,D/2.0,0.502941,0.502941,0.502941);
}
#endif

} /* namespace RL */
