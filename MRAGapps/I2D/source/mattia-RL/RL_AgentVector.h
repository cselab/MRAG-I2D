/*
 * RL_AgentVector.h
 *
 *  Created on: Mar 21, 2012
 *      Author: mgazzola
 */

#ifndef RL_AGENTVECTOR_H_
#define RL_AGENTVECTOR_H_

#include <vector>
#include "RL_Environment.h"
#include "RL_Agent.h"

namespace RL {

class RL_AgentVector : public RL_Agent
{
	long unsigned int counterAgents;
	bool shared;
	vector<RL_Agent *> agents;
	vector<bool> valids;
	vector<string> name_learning;
	vector<string> name_saveQ;
	MRAG::Profiler profiler;

public:
	map< string, vector<RL_Agent *> > data;

	RL_AgentVector(MRAG::ArgumentParser & parser, map< string, vector<RL_Agent *> > _data );
	~RL_AgentVector();

	// Online learning
	void savePolicy(string name=string());
	void restartPolicy(string name=string());

	void update(const double dt, const double t, map< string, vector<RL_Agent *> > * _data = NULL, string filename = string());
	bool choose(const double t, map< string, vector<RL_Agent *> > * _data = NULL );
	void learn(const double t, map< string, vector<RL_Agent *> > * _data = NULL, string name = string());
	void reward(const double t, map< string, vector<RL_Agent *> > * _data = NULL);

#ifdef _RL_VIZ
	virtual void paint();
#endif
};

} /* namespace RL */
#endif /* RL_AGENTVECTOR_H_ */
