/*
 * RL_AgentVector.cpp
 *
 *  Created on: Mar 21, 2012
 *      Author: mgazzola
 */
#include <stdio.h>
#include <stdlib.h>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "RL_AgentVector.h"

namespace RL {

struct TBB_start_choose
{
	double t;
	vector<RL_Agent *> & agents;
	map< string, vector<RL_Agent *> > * data;
	vector<bool> & valids;

	TBB_start_choose(double t, vector<RL_Agent *> & agents, map< string, vector<RL_Agent *> > * data, vector<bool> & valids): t(t),agents(agents),data(data),valids(valids)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			valids[i] = agents[i]->startChoose(t,data);
	}
};

struct TBB_choose_no_shared
{
	double t;
	vector<RL_Agent *> & agents;
	map< string, vector<RL_Agent *> > * data;
	vector<bool> & valids;

	TBB_choose_no_shared(double t, vector<RL_Agent *> & agents, map< string, vector<RL_Agent *> > * data, vector<bool> & valids): t(t),agents(agents),data(data),valids(valids)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
		{
			valids[i] = agents[i]->startChoose(t,data);
			agents[i]->choose(t,data);
		}
	}
};

struct TBB_choose
{
	double t;
	vector<RL_Agent *> & agents;
	map< string, vector<RL_Agent *> > * data;

	TBB_choose(double t, vector<RL_Agent *> & agents, map< string, vector<RL_Agent *> > * data): t(t),agents(agents),data(data)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->choose(t,data);
	}
};

struct TBB_learn
{
	double t;
	vector<RL_Agent *> & agents;
	vector<string> & names;
	map< string, vector<RL_Agent *> > * data;

	TBB_learn(double t, vector<RL_Agent *> & agents, map< string, vector<RL_Agent *> > * data,vector<string> & names): t(t),agents(agents),data(data),names(names)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->learn(t,data,names[i]);
	}
};

struct TBB_save
{
	vector<RL_Agent *> & agents;
	vector<string> & names;

	TBB_save(vector<RL_Agent *> & agents, vector<string> & names): agents(agents), names(names)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->savePolicy(names[i]);
	}
};

struct TBB_restart
{
	vector<RL_Agent *> & agents;
	vector<string> & names;

	TBB_restart(vector<RL_Agent *> & agents, vector<string> & names): agents(agents), names(names)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->restartPolicy(names[i]);
	}
};

struct TBB_stop_test
{
	double t;
	vector<RL_Agent *> & agents;
	map< string, vector<RL_Agent *> > * data;

	TBB_stop_test(double t, vector<RL_Agent *> & agents, map< string, vector<RL_Agent *> > * data): t(t),agents(agents),data(data)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->stopTest(t,data);
	}
};

struct TBB_set_test
{
	double t;
	vector<RL_Agent *> & agents;

	TBB_set_test(double t, vector<RL_Agent *> & agents): t(t),agents(agents)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->setReward(t);
	}
};

struct TBB_learn_no_shared_test
{
	double t;
	vector<RL_Agent *> & agents;
	map< string, vector<RL_Agent *> > * data;
	vector<string> & names;

	TBB_learn_no_shared_test(double t, vector<RL_Agent *> & agents, map< string, vector<RL_Agent *> > * data, vector<string> & names): t(t), agents(agents), data(data), names(names)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
		{
			RL_Agent * a = agents[i];
			a->setReward(t);
			a->mapTest(t,data);
			a->stopTest(t,data);
			a->learn(t,data,names[i]);
		}
	}
};

struct TBB_map_test
{
	double t;
	vector<RL_Agent *> & agents;
	map< string, vector<RL_Agent *> > * data;

	TBB_map_test(double t, vector<RL_Agent *> & agents, map< string, vector<RL_Agent *> > * data): t(t),agents(agents),data(data)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->mapTest(t,data);
	}
};

struct TBB_reward
{
	double t;
	vector<RL_Agent *> & agents;
	map< string, vector<RL_Agent *> > * data;

	TBB_reward(double t, vector<RL_Agent *> & agents, map< string, vector<RL_Agent *> > * data): t(t),agents(agents),data(data)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->reward(t,data);
	}
};

struct TBB_update
{
	double t,dt;
	vector<RL_Agent *> & agents;
	map< string, vector<RL_Agent *> > * data;

	TBB_update(double dt, double t, vector<RL_Agent *> & agents, map< string, vector<RL_Agent *> > * data): dt(dt),t(t),agents(agents),data(data)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->update(dt,t,data);
	}
};

RL_AgentVector::RL_AgentVector(MRAG::ArgumentParser & parser, map< string, vector<RL_Agent *> > _data): RL_Agent(parser), data(_data), counterAgents(0)
{
	shared = parser("-shared").asBool();

	// Prepare stuff for tbb
	for( map< string, vector<RL_Agent *> >::iterator it = data.begin(); it!=data.end(); ++it)
		for( vector<RL_Agent *>::iterator it2 = it->second.begin(); it2!=it->second.end(); ++it2)
			agents.push_back(*it2);

	valids.resize(agents.size(),true);

	for( map< string, vector<RL_Agent *> >::iterator it = data.begin(); it!=data.end(); ++it)
	{
		string type = it->first;
		unsigned int counter = 0;
		for( vector<RL_Agent *>::iterator it2 = it->second.begin(); it2!=it->second.end(); ++it2)
		{
			char ending[100];

			//if(!shared)
			sprintf(ending, "_%04d", counter);
			//else
			//	sprintf(ending, "_%04d", 0);

			string nameL("learning_" + type + ending);
			string nameS("saveQ_" + type + ending);
			name_learning.push_back(nameL);
			name_saveQ.push_back(nameS);
			counter++;
		}
	}
}

RL_AgentVector::~RL_AgentVector()
{
	for( vector<RL_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		if( (*it)!=NULL )
		{
			delete (*it);
			(*it) = NULL;
		}

	agents.clear();
	data.clear();

	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

void RL_AgentVector::update(const double dt, const double t, map< string, vector<RL_Agent *> > * _data, string filename )
{
	TBB_update tbb_update(dt,t,agents,&data);
	tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_update, tbb::auto_partitioner());

	// Common dump serial
	if(counterAgents%10000==0)
	{
		ofstream out("commondump", ios_base::app);
		for( vector<RL_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
			(*it)->dump(t,out);

		out.close();
	}

	counterAgents++;
}

bool RL_AgentVector::choose(const double t, map< string, vector<RL_Agent *> > * _data )
{
	bool valid = true;

	if(!shared)
	{
		TBB_choose_no_shared choose_no_shared(t,agents,&data,valids);
		tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), choose_no_shared, tbb::auto_partitioner());
	}
	else
	{
		TBB_start_choose tbb_start_choose(t,agents,&data,valids);
		tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_start_choose, tbb::auto_partitioner());

		for( vector<RL_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
			(*it)->choose(t,&data);
	}

	for( vector<bool>::iterator it=valids.begin(); it!=valids.end(); ++it)
		valid *= (*it);

	return valid;
}

void RL_AgentVector::learn(const double t, map< string, vector<RL_Agent *> > * _data, string name )
{
	if(!shared)
	{
		// Everything parallel PARALLEL
		TBB_learn_no_shared_test tbb_learn_no_shared_test(t,agents,&data,name_learning);
		tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_learn_no_shared_test, tbb::auto_partitioner());
	}
	else
	{
		// Set reward SERIAL
		for( vector<RL_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
			(*it)->setReward(t);

		// Map end state PARALLEL
		TBB_map_test tbb_map_test(t,agents,&data);
		tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_map_test, tbb::auto_partitioner());

		// Stop test SERIAL
		for( vector<RL_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
			(*it)->stopTest(t,&data);

		// Learn SERIAL
		unsigned int counter = 0;
		for( vector<RL_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		{
			(*it)->learn(t,&data,name_learning[counter]);
			counter++;
		}
	}
}

void RL_AgentVector::reward(const double t, map< string, vector<RL_Agent *> > * _data )
{
	TBB_reward tbb_reward(t,agents,&data);
	tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_reward, tbb::auto_partitioner());
}

void RL_AgentVector::savePolicy(string name)
{
	if(!shared)
	{
		TBB_save tbb_save(agents,name_saveQ);
		tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_save, tbb::auto_partitioner());
	}
	else
	{
		unsigned int totalCounter = 0;
		for( map< string, vector<RL_Agent *> >::iterator it = data.begin(); it!=data.end(); ++it)
		{
			string type = it->first;
			unsigned int counter = 0;
			for( vector<RL_Agent *>::iterator it2 = it->second.begin(); it2!=it->second.end(); ++it2)
			{
				if(counter==0)
					(*it2)->savePolicy(name_saveQ[totalCounter]);

				counter++;
				totalCounter++;
			}
		}
	}
}

void RL_AgentVector::restartPolicy(string name)
{
	if(!shared)
	{
		// Everything parallel PARALLEL
		TBB_restart tbb_restart(agents,name_saveQ);
		tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_restart, tbb::auto_partitioner());
	}
	else
	{
		unsigned int totalCounter = 0;
		for( map< string, vector<RL_Agent *> >::iterator it = data.begin(); it!=data.end(); ++it)
		{
			string type = it->first;
			unsigned int counter = 0;
			for( vector<RL_Agent *>::iterator it2 = it->second.begin(); it2!=it->second.end(); ++it2)
			{
				if(counter==0)
					(*it2)->restartPolicy(name_saveQ[totalCounter]);

				counter++;
				totalCounter++;
			}
		}
	}
}

#ifdef _RL_VIZ
void RL_AgentVector::paint()
{
	for( vector<RL_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->paint();
}
#endif

} /* namespace RL */
