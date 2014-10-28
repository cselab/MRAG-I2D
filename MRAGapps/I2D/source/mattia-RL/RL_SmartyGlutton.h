/*
 * RL_SmartyGlutton.h
 *
 *  Created on: Mar 27, 2012
 *      Author: dalmassg
 */

#ifndef RL_SMARTYGLUTTON_H_
#define RL_SMARTYGLUTTON_H_
#include <iostream>
#include <fstream>
#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "RL_Environment.h"
#include "RL_Agent.h"

namespace RL
{

class RL_SmartyGlutton : public RL_Agent
{
protected:
	bool STARVING;
	int IDcollision;
	vector<int> signature;
	const Real D,T;
	Real x,y,vx,vy,modv,Vm[2];
	const int nactions;
	const int NN;
	int myaction;
	int myfoodcounter;
	const Real cutoffGluttons;
	const Real cutoffFood;
	const Real steer;
	const Real noise;
	const Real selfAvoid;
	RNG rng;

	void _rotateSimple(const Real theta, Real x[2]) const;
	void _mapForMaxRadius(const Real cutoff, const int idxDistLevel, const int idxAngleLevel, int & idxDist, int & idxAngle, map< string, vector<RL_Agent *> > * _data);
	void _mapForFoods(const Real cutoff, const int idxDistLevel, const int idxAngleLevel, int & idxDist, int & idxAngle, map< string, vector<RL_Agent *> > * _data);
	void _mapForGluttons(const Real cutoff, const int idxDistLevel, const int idxAngleLevel, const int idxAngleNeighLevel, RL_SmartyGlutton * neigh, Real vneigh[2], int & idxDist, int & idxAngle, int & idxAngleNeigh, map< string, vector<RL_Agent *> > * _data);

public:

	RL_SmartyGlutton(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const Real _T, const int _nactions, const Real _cutoffGluttons, const Real _cutoffFood, const Real _steer, const Real _noise, const Real _selfAvoid, const int NN, const int _ID, RL_TabularPolicy ** _policy = NULL, const int seed = 0);
	virtual ~RL_SmartyGlutton();

	virtual void update(const double dt, const double t, map< string, vector<RL_Agent *> > * _data = NULL, string filename = string());
	virtual void mapAction(int action);
	virtual bool mapState(vector<int> & state, map< string, vector<RL_Agent *> > * _data = NULL );
	virtual void reward(const double t, map< string, vector<RL_Agent *> > * _data = NULL);
	virtual void dump(const double t, ofstream & out);

#ifdef _RL_VIZ
	virtual void paint();
#endif
};

} /* namespace RL */
#endif /* RL_SMARTYGLUTTON_H_ */
