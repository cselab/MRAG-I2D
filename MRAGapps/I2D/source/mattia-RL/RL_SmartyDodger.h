/*
 * RL_SmartyDodger.h
 *
 *  Created on: Mar 22, 2012
 *      Author: mgazzola
 */

#ifndef RL_SMARTYDODGER_H_
#define RL_SMARTYDODGER_H_
#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "RL_Environment.h"
#include "RL_Agent.h"

namespace RL
{

class RL_SmartyDodger : public RL_Agent
{
protected:
	vector<int> signature;
	const Real D,T;
	Real x,y,vx,vy,modv;

	void _mapForMaxRadius(const Real cutoff, const int idxDistLevel, const int idxAngleLevel, int & idxDist, int & idxAngle, map< string, vector<RL_Agent *> > * _data);
	void _mapForColumns(const Real cutoff, const int idxDistLevel, const int idxAngleLevel, int & idxDist, int & idxAngle, map< string, vector<RL_Agent *> > * _data);
	void _mapForDynamicColumns(const Real cutoff, const int idxDistLevel, const int idxAngleLevel, int & idxDist, int & idxAngle, map< string, vector<RL_Agent *> > * _data);

public:

	RL_SmartyDodger(MRAG::ArgumentParser & parser, const Real _x, const Real _y, const Real _D, const Real _T, const int _ID = 0, RL_TabularPolicy ** _policy = NULL, const int seed = 0);
	virtual ~RL_SmartyDodger();

	virtual void update(const double dt, const double t, map< string, vector<RL_Agent *> > * _data = NULL, string filename = string());
	virtual void mapAction(int action);
	virtual bool mapState(vector<int> & state, map< string, vector<RL_Agent *> > * _data = NULL );
	virtual void reward(const double t, map< string, vector<RL_Agent *> > * _data = NULL);

#ifdef _RL_VIZ
	virtual void paint();
#endif
};

} /* namespace RL */
#endif /* RL_SMARTYDODGER_H_ */
