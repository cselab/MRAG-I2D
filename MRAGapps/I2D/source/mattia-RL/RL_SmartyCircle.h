/*
 * RL_SmartyCircle.h
 *
 *  Created on: May 27, 2011
 *      Author: mgazzola
 */

#ifndef RL_SMARTPOINT_H_
#define RL_SMARTPOINT_H_

#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "RL_Environment.h"
#include "RL_Agent.h"

namespace RL
{

class RL_SmartyCircle : public RL_Agent
{
protected:
	vector<int> signature;
	const Real D,T,maxDomainRadius;
	Real x,y,vx,vy,modv;

public:

	RL_SmartyCircle(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const int _ID = 0, RL_TabularPolicy ** _policy = NULL, const int seed = 0);
	virtual ~RL_SmartyCircle();

	virtual void update(const Real dt, const Real t, map< string, vector<RL_Agent *> > * _data = NULL, string filename = string());
	virtual void mapAction(int action);
	virtual bool mapState(vector<int> & state, map< string, vector<RL_Agent *> > * _data = NULL );
	virtual void reward(const Real t, map< string, vector<RL_Agent *> > * _data = NULL);

#ifdef _RL_VIZ
	virtual void paint();
#endif
};

}

#endif /* RL_SMARTPOINT_H_ */
