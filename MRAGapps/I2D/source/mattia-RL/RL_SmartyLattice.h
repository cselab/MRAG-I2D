/*
 * RL_SmartyLattice.h
 *
 *  Created on: Apr 5, 2012
 *      Author: mgazzola
 */

#ifndef RL_SMARTYLATTICE_H_
#define RL_SMARTYLATTICE_H_

#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "RL_Environment.h"
#include "RL_Agent.h"
#include "rng.h"

namespace RL
{

class RL_SmartyLattice : public RL_Agent
{
protected:
	Real dir[2];
	Real target[2];
	vector<int> signature;
	const Real D,T,maxDomainRadius;
	Real x,y,vx,vy,modv;
	RNG rng;

public:

	RL_SmartyLattice(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const Real _dir[2], const int _ID = 0, RL_TabularPolicy ** _policy = NULL, const int seed = 0);
	virtual ~RL_SmartyLattice();

	virtual void update(const Real dt, const Real t, map< string, vector<RL_Agent *> > * _data = NULL, string filename = string());
	virtual void mapAction(int action);
	virtual bool mapState(vector<int> & state, map< string, vector<RL_Agent *> > * _data = NULL );
	virtual void reward(const Real t, map< string, vector<RL_Agent *> > * _data = NULL);

#ifdef _RL_VIZ
	virtual void paint();
#endif
};

} /* namespace RL */
#endif /* RL_SMARTYLATTICE_H_ */
