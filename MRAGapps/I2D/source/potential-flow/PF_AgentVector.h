/*
 *  PF_AgentVector.h
 *  DipoleCode
 *
 *  Created by Alexia on 3/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#include <vector>
#include "PF_Agent.h"

namespace PF
{

class PF_AgentVector : public PF_Agent
{
public:
	//constructor & destructor
	PF_AgentVector(MRAG::ArgumentParser &parser, const Real lr, const Real greedyEps, bool isSmooth, bool isControlled, map< string, vector<PF_Agent*> > collection);
	~PF_AgentVector();
	
	//method
	int numberOfAgents();

	void update(const Real dt, const Real time, map < string, vector <PF_Agent*> > *collection = NULL);
	void updatePostVelocityUpdate(const Real dt, const Real time, map < string, vector <PF_Agent*> > *collection = NULL);
	void storePosition(vector <pair <Real,Real> > &position, map < pair <Real,Real>, PF_Agent*> *positionAgent = NULL);

	void getCoordinates(vector<pair<Real, Real> >& coordinates);
	void getVortices(vector <Vortex> &vortices);
	void getNominalGammas(vector <Real> & nominalGammas);
	void getTargets(vector <pair <Real,Real> > &targets, map < pair <Real,Real>, PF_Agent*> *targetsAgent = NULL);
	void getForwardVelocities(vector <Real> &forwardVelocities);
	void getTargetPoints(vector<pair<Real, Real> > &targetPoints);

	void saveData(Real time);

	void getFitnessValues(vector<Real> &fitnesses);
	void getErrorValues(vector<Real> &errors);

	bool choose(const Real time, map < string, vector <PF_Agent*> > *collection = NULL);
	void reward(const Real time, map < string, vector <PF_Agent*> > *collection = NULL);
	void fitness();
	void error();
	void learn(const Real time, map < string, vector <PF_Agent*> > *collection = NULL, string filename = string());
	void savePolicy(string name = string());
	void restartPolicy(string name = string());
	void resetToInitialCondition();

#ifdef _RL_VIZ
	void paint();
	void paint(Real scale, const Real center[2]);
#endif

private:
	//attributs
	bool shared;
	vector<string> name_learning;
	vector<string> name_saveQ;
	vector<PF_Agent*> agents;
	vector<bool> valids;
	map < string, vector<PF_Agent*> > agentCollection;
};

} /*namespace PF*/
	
