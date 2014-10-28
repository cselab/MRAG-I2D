/*
 *  PF_AgentVector.cpp
 *  DipoleCode
 *
 *  Created by Alexia on 3/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <complex>
#include <math.h>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "PF_AgentVector.h"

using namespace PF;

struct TBB_map_start
{
	Real t;
	vector<PF_Agent *> & agents;
	map< string, vector<PF_Agent *> > * agentCollection;
	vector<bool> & valids;

	TBB_map_start(Real t, vector<PF_Agent *> & agents, map< string, vector<PF_Agent *> > * agentCollection, vector<bool> & valids): t(t),agents(agents),agentCollection(agentCollection),valids(valids)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			valids[i] = agents[i]->mapStart(t,agentCollection);
	}
};

struct TBB_choose_no_shared
{
	Real t;
	vector<PF_Agent *> & agents;
	map< string, vector<PF_Agent *> > * agentCollection;
	vector<bool> & valids;

	TBB_choose_no_shared(Real t, vector<PF_Agent *> & agents, map< string, vector<PF_Agent *> > * agentCollection, vector<bool> & valids): t(t),agents(agents),agentCollection(agentCollection),valids(valids)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
		{
			valids[i] = agents[i]->mapStart(t,agentCollection);
			agents[i]->choose(t,agentCollection);
		}
	}
};

struct TBB_save
{
	vector<PF_Agent *> & agents;
	vector<string> & names;

	TBB_save(vector<PF_Agent *> & agents, vector<string> & names): agents(agents), names(names)
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
	vector<PF_Agent *> & agents;
	vector<string> & names;

	TBB_restart(vector<PF_Agent *> & agents, vector<string> & names): agents(agents), names(names)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->restartPolicy(names[i]);
	}
};

struct TBB_learn_no_shared_test
{
	Real t;
	vector<PF_Agent *> & agents;
	map< string, vector<PF_Agent *> > * agentCollection;
	vector<string> & names;

	TBB_learn_no_shared_test(Real t, vector<PF_Agent *> & agents, map< string, vector<PF_Agent *> > * agentCollection, vector<string> & names): t(t), agents(agents), agentCollection(agentCollection), names(names)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
		{
			PF_Agent * a = agents[i];
			a->setReward(t);
			a->mapTest(t,agentCollection);
			a->stopTest(t,agentCollection);
			a->learn(t,agentCollection,names[i]);
		}
	}
};

struct TBB_map_test
{
	Real t;
	vector<PF_Agent *> & agents;
	map< string, vector<PF_Agent *> > * agentCollection;

	TBB_map_test(Real t, vector<PF_Agent *> & agents, map< string, vector<PF_Agent *> > * agentCollection): t(t),agents(agents),agentCollection(agentCollection)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->mapTest(t,agentCollection);
	}
};

struct TBB_reward
{
	Real t;
	vector<PF_Agent *> & agents;
	map< string, vector<PF_Agent *> > *agentCollection;

	TBB_reward(Real t, vector<PF_Agent *> & agents, map< string, vector<PF_Agent *> > * agentCollection): t(t),agents(agents),agentCollection(agentCollection)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->reward(t,agentCollection);
	}
};

struct TBB_fitness
{
	vector<PF_Agent *> & agents;

	TBB_fitness(vector<PF_Agent *> & agents) :
			agents(agents)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->fitness();
	}
};

struct TBB_error
{
	Real t;
	vector<PF_Agent *> & agents;

	TBB_error(vector<PF_Agent *> & agents) :
			agents(agents)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->error();
	}
};


struct TBB_update
{
	Real t,dt;
	vector<PF_Agent *> & agents;
	map < string, vector <PF_Agent*> > *agentCollection;

	TBB_update(Real dt, Real t, vector<PF_Agent *> & agents, map < string, vector <PF_Agent*> > *agentCollection): dt(dt),t(t),agents(agents),agentCollection(agentCollection)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->update(dt,t,agentCollection);
	}
};

struct TBB_updatePostVelocityUpdate
{
	Real t,dt;
	vector<PF_Agent *> & agents;
	map < string, vector <PF_Agent*> > *agentCollection;

	TBB_updatePostVelocityUpdate(Real dt, Real t, vector<PF_Agent *> & agents, map < string, vector <PF_Agent*> > *agentCollection): dt(dt),t(t),agents(agents),agentCollection(agentCollection)
	{
	}

	void operator() ( const tbb::blocked_range<size_t> &r ) const
	{
		for (size_t i=r.begin(); i!=r.end();++i)
			agents[i]->updatePostVelocityUpdate(dt,t,agentCollection);
	}
};

PF_AgentVector::PF_AgentVector(MRAG::ArgumentParser &parser, const Real lr, const Real greedyEps, bool isSmooth, bool isControlled, map< string, vector<PF_Agent*> > collection): PF_Agent(parser, lr, greedyEps, isSmooth, isControlled), agentCollection(collection)
{
	shared = parser("-shared").asBool();

	// Prepare stuff for tbb
	for( map< string, vector<PF_Agent *> >::iterator it = agentCollection.begin(); it!=agentCollection.end(); ++it)
		for( vector<PF_Agent *>::iterator it2 = it->second.begin(); it2!=it->second.end(); ++it2)
			agents.push_back(*it2);

	valids.resize(agents.size(),true);

	for( map< string, vector<PF_Agent *> >::iterator it = agentCollection.begin(); it!=agentCollection.end(); ++it)
	{
		string type = it->first;
		unsigned int counter = 0;
		for( vector<PF_Agent *>::iterator it2 = it->second.begin(); it2!=it->second.end(); ++it2)
		{
			char ending[100];
			if(!shared)
				sprintf(ending, "_%04d", counter);
			else
				sprintf(ending, "_%04d", 0);

			string nameL("learning_" + type + ending);
			string nameS("saveQ_" + type + ending);
			name_learning.push_back(nameL);
			name_saveQ.push_back(nameS);
			counter++;
		}
	}
}

PF_AgentVector::~PF_AgentVector()
{	
	for(vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		if((*it) != NULL)
		{
			delete (*it);
			(*it) = NULL;
		}

	agentCollection.clear();
	agents.clear();

	if(policy != NULL)
		if((*policy) != NULL)
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

int PF_AgentVector::numberOfAgents()
{
	return (int)agents.size();
}

void PF_AgentVector::update(const Real dt, const Real time, map < string, vector <PF_Agent*> > *collection)
{	
	TBB_update tbb_update(dt,time,agents,&agentCollection);
	tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_update, tbb::auto_partitioner());
}

void PF_AgentVector::updatePostVelocityUpdate(const Real dt, const Real time, map<string, vector<PF_Agent*> >* collection)
{
	TBB_updatePostVelocityUpdate tbb_updatePostVelocityUpdate(dt, time, agents, &agentCollection);
	tbb::parallel_for(tbb::blocked_range<size_t>(0, agents.size()), tbb_updatePostVelocityUpdate, tbb::auto_partitioner());
}

void PF_AgentVector::saveData(Real time)
{
	for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->saveData(time);
}

bool PF_AgentVector::choose(const Real time, map < string, vector <PF_Agent*> > *collection)
{
	bool valid = true;

	if(!shared)
	{
		TBB_choose_no_shared choose_no_shared(time,agents,&agentCollection,valids);
		tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), choose_no_shared, tbb::auto_partitioner());
	}
	else
	{
		TBB_map_start tbb_map_start(time,agents,&agentCollection,valids);
		tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_map_start, tbb::auto_partitioner());

		for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
			(*it)->choose(time,&agentCollection);
	}

	for( vector<bool>::iterator it=valids.begin(); it!=valids.end(); ++it)
		valid *= (*it);

	return valid;
}

void PF_AgentVector::learn(const Real time, map < string, vector <PF_Agent*> > *collection, string filename)
{
	if(!shared)
	{
		// Everything parallel PARALLEL
		TBB_learn_no_shared_test tbb_learn_no_shared_test(time,agents,&agentCollection,name_learning);
		tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_learn_no_shared_test, tbb::auto_partitioner());
	}
	else
	{
		// Set reward SERIAL
		for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
			(*it)->setReward(time);

		// Map end state PARALLEL
		TBB_map_test tbb_map_test(time,agents,&agentCollection);
		tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_map_test, tbb::auto_partitioner());

		// Stop test SERIAL
		for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
			(*it)->stopTest(time,&agentCollection);

		// Learn SERIAL
		unsigned int counter = 0;
		for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		{
			(*it)->learn(time,&agentCollection,name_learning[counter]);
			counter++;
		}
	}
}

void PF_AgentVector::reward(const Real time, map < string, vector <PF_Agent*> > *collection)
{
	TBB_reward tbb_reward(time,agents,&agentCollection);
	tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_reward, tbb::auto_partitioner());
}

void PF_AgentVector::fitness()
{
	TBB_fitness tbb_fitness(agents);
	tbb::parallel_for(tbb::blocked_range<size_t>(0, agents.size()), tbb_fitness, tbb::auto_partitioner());
}

void PF_AgentVector::error()
{
	TBB_error tbb_error(agents);
	tbb::parallel_for(tbb::blocked_range<size_t>(0, agents.size()), tbb_error, tbb::auto_partitioner());
}

void PF_AgentVector::savePolicy(string name)
{
	if(!shared)
	{
		TBB_save tbb_save(agents,name_saveQ);
		tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_save, tbb::auto_partitioner());
	}
	else
	{
		unsigned int counter = 0;

		vector<PF_Agent *>::iterator it = agents.begin();
		(*it)->savePolicy(name_saveQ[counter]);

//		for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
//		{
//			(*it)->savePolicy(name_saveQ[counter]);
//			counter++;
//		}
	}
}

void PF_AgentVector::restartPolicy(string name)
{
	if(!shared)
	{
		// Everything parallel PARALLEL
		TBB_restart tbb_restart(agents,name_saveQ);
		tbb::parallel_for(tbb::blocked_range<size_t>(0,agents.size()), tbb_restart, tbb::auto_partitioner());
	}
	else
	{
		unsigned int counter = 0;
		vector<PF_Agent *>::iterator it = agents.begin();
		(*it)->restartPolicy(name_saveQ[counter]);

//		for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
//		{
//			(*it)->restartPolicy(name_saveQ[counter]);
//			counter++;
//		}
	}
}

void PF_AgentVector::resetToInitialCondition()
{
	for (vector<PF_Agent *>::iterator it = agents.begin(); it != agents.end(); ++it)
		(*it)->resetToInitialCondition();
}

void PF_AgentVector::storePosition(vector <pair <Real,Real> > &position, map < pair <Real,Real>, PF_Agent*> *positionAgent)
{
	position.clear();
	(*positionAgent).clear();
	int size = 0;

	for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
	{
		if (!(*it)->isDead)
		{
			(*it)->storePosition(position);
			int size_next = position.size();
			int diff = size_next-size;
			for (int i = diff; i > 0; i--)
			{
				pair <Real,Real> coord = position[size_next-i];
				(*positionAgent)[coord] = (*it);
			}
			size = size_next;
		}
	}
}

void PF_AgentVector::getCoordinates(vector<pair<Real, Real> >& coordinates)
{
	coordinates.clear();
	for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
			(*it)->getCoordinates(coordinates);
}

void PF_AgentVector::getVortices(vector <Vortex> &vortices)
{
	vortices.clear();

	for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->getVortices(vortices);

}

void PF_AgentVector::getNominalGammas(vector <Real> &nominalGammas)
{
	nominalGammas.clear();
	for(vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->getNominalGammas(nominalGammas);
}

void PF_AgentVector::getTargets(vector <pair <Real,Real> > &targets, map < pair <Real,Real>, PF_Agent*> *targetsAgent)
{
	targets.clear();
	(*targetsAgent).clear();
	int size = 0;

	for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
	{
		(*it)->getTargets(targets);
		int size_next = targets.size();
		int diff = size_next-size;
		for (int i = diff; i > 0; i--)
		{
			pair <Real,Real> coord = targets[size_next-i];
			(*targetsAgent)[coord] = (*it);
		}
		size = size_next;
	}
}

void PF_AgentVector::getForwardVelocities(vector <Real> &forwardVelocities)
{
	forwardVelocities.clear();
	for(vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->getForwardVelocities(forwardVelocities);
}


void PF_AgentVector::getFitnessValues(vector<Real> &fitnesses)
{
	fitnesses.clear();
	for (vector<PF_Agent *>::iterator it = agents.begin(); it != agents.end(); ++it)
		(*it)->getFitnessValues(fitnesses);
}

void PF_AgentVector::getErrorValues(vector<Real> &errors)
{
	errors.clear();
	for (vector<PF_Agent *>::iterator it = agents.begin(); it != agents.end(); ++it)
		(*it)->getErrorValues(errors);
}

void PF_AgentVector::getTargetPoints(vector<pair<Real, Real> > &targetPoints)
{
	for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->getTargetPoints(targetPoints);
}

#ifdef _RL_VIZ
void PF_AgentVector::paint()
{
	for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->paint();
}

void PF_AgentVector::paint(Real scale, const Real center[2])
{
	for( vector<PF_Agent *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->paint(scale, center);
}
#endif
