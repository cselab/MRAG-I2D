/*
 * RL_Food.h
 *
 *  Created on: Mar 27, 2012
 *      Author: dalmassg
 */

#ifndef RL_FOOD_H_
#define RL_FOOD_H_

#include <stdio.h>
#include <stdlib.h>
#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "RL_Agent.h"

namespace RL
{

class RL_Food : public RL_Agent
{
protected:
	RNG rng;
	long unsigned int count;

public:
	const Real D;
	Real x,y;
	RL_Food(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const int _ID = 0, const int seed = 0, RL_TabularPolicy ** _policy = NULL);
	virtual ~RL_Food();

	virtual void update(const double dt, const double t, map< string, vector<RL_Agent *> > * _data = NULL, string filename = string());

#ifdef _RL_VIZ
	virtual void paint();
#endif
};

} /* namespace RL */

#endif /* RL_FOOD_H_ */
