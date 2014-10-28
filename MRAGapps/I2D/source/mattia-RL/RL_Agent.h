/*
 * RL_Agent.h
 *
 *  Created on: May 27, 2011
 *      Author: mgazzola
 */

#ifndef RL_AGENT_H_
#define RL_AGENT_H_

#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "RL_Environment.h"
#include "RL_TabularPolicy.h"
#include "RL_QLearning.h"

#ifdef _RL_VIZ
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include "GLUT/glut.h"
#else
#include <GL/gl.h>
#endif
#endif

using namespace std;

namespace RL
{

class RL_Agent
{
protected:
	// Online learning
	enum QLearningStatus{ Waiting, Ready };
	QLearningStatus status;
	RL::RL_TabularPolicy ** policy;
	RL::RL_TabularPolicy * policystore;
	Real integralReward;
	Real learningInterval;
	Real learningTimer;
	vector<int> myState;
	vector<int> myStateAction;
	int myAction;
	bool shared;

	// Gives back the vector between x2 and x1   x2 * <------- * x1
	void _dist(const Real x2[2], const Real x1[2], Real d[2]) const;
	void _distPeriodic(const Real x2[2], const Real x1[2], Real d[2]) const;

	// Gives back the angle between v1 and v2
	// angle(v1,v2) = -angle(v2,v1)
	//											v1
	//                                        /
	// 										 /____ v2
	Real _angleVectors(const Real v1[2], const Real v2[2]) const;

	// Tells you if x2 is to your right given a certain direction
	//  x1 *------> dir
	//			|
	//			|
	//			*x2
	bool _isRight(const Real x1[2], const Real dir[2], const Real x2[2]) const;

	// Rotate a vector
	void _rotate(const Real alpha, const Real v1[2], Real v2[2]) const;

	Real _modv(const Real v[2]) const;

	Real _normalize(const Real v[2], Real n[2]) const;

	int _discretize(const Real signal, const Real minsignal, const Real maxsignal, const int nlevels) const;

	int _discretize1SS(const Real signal, const Real minsignal, const Real maxsignal, const int nlevels, const bool deadEnd) const;

	// Discretize and account for two special states when the signal is beyond the cutoff and when it s below a certain value,
	// i.e. distance < 0 and therefore agent inside another object. The priority is given to the inside obstacle state.
	int _discretize2SS(const Real signal, const Real minsignal, const Real maxsignal, const int nlevels, const bool deadEnd = false, const bool inside = false) const;

public:

	int ID;
	string name;

	RL_Agent(MRAG::ArgumentParser & parser, string _name = string(), const int _ID = 0, RL::RL_TabularPolicy ** _policy = NULL, const int seed = 0) :
		name(_name), ID(_ID), policy(_policy), policystore(NULL), shared(true), integralReward(0.0), status(Ready), learningInterval(0.0), learningTimer(0.0)
	{
		if(policy!=NULL)
			if((*policy)==NULL)
			{
				const Real LR = parser("-lr").asDouble();
				const Real GAMMA = parser("-gamma").asDouble();
				const Real GREEDYEPS = parser("-greedyEps").asDouble();
				const int learnDump = parser("-learnDump").asInt();
				shared = parser("-shared").asBool();

				assert( LR>=0.0 && LR<=1.0);
				assert( GAMMA>=0.0 && GAMMA<=1.0 );
				assert( GREEDYEPS>=0.0 && GREEDYEPS<=1.0 );
				assert( (LR+GAMMA)<=1.0 );

				if(LR>0.0 && LR<=1.0 && GAMMA>=0.0 && GAMMA<=1.0 && GREEDYEPS>=0.0 && GREEDYEPS<=1.0 && (LR+GAMMA)>1.0)
				{
					printf("RL parameters are fucked up!\n");
					abort();
				}

				policystore = new RL_QLearning(LR,GAMMA,GREEDYEPS,seed,learnDump);
				policy = &policystore;
				assert(policy!=NULL);
			}
	};

	virtual ~RL_Agent()
	{
		if(policystore!=NULL)
		{
			delete policystore;
			policystore = NULL;
		}

		if(policy!=NULL)
			if( (*policy)!=NULL )
			{
				delete (*policy);
				(*policy) = NULL;
			};
	}


	// Methods implemented here once for all
	bool startChoose(const double t, map< string, vector<RL_Agent *> > * _data = NULL);
	void setReward(const double t);
	void mapTest(const double t, map< string, vector<RL_Agent *> > * _data);
	void stopTest(const double t, map< string, vector<RL_Agent *> > * _data = NULL);

	// Methods re-implemented ONLY by agent vector class
	virtual void savePolicy(string name=string());
	virtual void restartPolicy(string name=string());
	virtual bool choose(const double t, map< string, vector<RL_Agent *> > * _data = NULL);
	virtual void learn(const double t, map< string, vector<RL_Agent *> > * _data = NULL, string name = string());

	// Methods re-implementable by all agents
	virtual void update(const double dt, const double t, map< string, vector<RL_Agent *> > * _data = NULL, string filename = string()){};
	virtual void reward(const double t, map< string, vector<RL_Agent *> > * _data = NULL){};
	virtual void mapAction(int action){};
	virtual bool mapState(vector<int> & state, map< string, vector<RL_Agent *> > * _data = NULL ){ return true; };
	virtual void dump(const double t, ofstream & out){};

#ifdef _RL_VIZ
	void _paintSphere(Real x, Real y, Real radius, Real r = 0.0, Real g = 0.0, Real b = 0.0) const;
	void _drawCircle(Real radius, Real xc = 0.5, Real yc = 0.5, Real r = 0.0, Real g = 0.0, Real b = 0.0) const;
	void _drawFullCircle(Real radius, Real xc = 0.5, Real yc = 0.5, Real r = 0.0, Real g = 0.0, Real b = 0.0) const;
	virtual void paint() = 0;
#endif
};

}

#endif /* RL_AGENT_H_ */
