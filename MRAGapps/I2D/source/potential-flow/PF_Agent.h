/*
 *  PF_Agent.h
 *  DipoleCode
 *
 *  Created by Alexia on 3/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#pragma once

#include <vector>
#include <complex>
#include <map>
#include <string>
#include <numeric>
#include "RL_Environment.h"
#include "RL_TabularPolicy.h"

#ifdef _RL_VIZ
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include "GLUT/glut.h"
#else
#include <GL/gl.h>
#endif
#endif

const complex <Real> I(0,1);

namespace PF
{

struct Vortex
{
	Real x;
	Real y;
	Real gamma;
};

class PF_Agent
{
public:
	//constructor & destructor
	PF_Agent(MRAG::ArgumentParser & parser, const Real lr, const Real greedyEps, bool isSmooth, bool isControlled, string _name = string(), const int _ID = 0, RL::RL_TabularPolicy **_policy = NULL, const int seed = 0);
	virtual ~PF_Agent();
	
	//methods	
	virtual int numberOfAgents(){return 0;};

	virtual void update(const Real dt, const Real time, map < string, vector <PF_Agent*> > *collection = NULL) = 0;
	virtual void updatePostVelocityUpdate(const Real dt, const Real time, map < string, vector <PF_Agent*> > *collection = NULL){} ;
	virtual void storePosition(vector <pair <Real,Real> > &position, map < pair <Real,Real>, PF_Agent*> *positionAgent = NULL){};
	virtual void setVelocity(complex <Real> velocityAgent){};
	virtual void reachedDtMin(){};
	virtual void saveData(Real time) = 0;

	virtual void getCoordinates(vector<pair<Real, Real> >& coordinates){};
	virtual void getVortices(vector <Vortex> &vortices){};
	virtual void getNominalGammas(vector <Real> & nominalGammas){};
	virtual void getTargets(vector <pair <Real,Real> > &targets, map < pair <Real,Real>, PF_Agent*> *targetsAgent = NULL){};
	virtual void getForwardVelocities(vector <Real> &forwardVelocities){};
	virtual void getTargetPoints(vector<pair<Real, Real> > &targetPoints){};

	virtual void getFitnessValues(vector<Real> &fitnesses){};
	virtual void getErrorValues(vector<Real> &errors){};
	
	bool mapStart(const Real time, map< string, vector<PF_Agent*> > *collection = NULL);
	void setReward(const Real time);
	void mapTest(const Real time, map < string, vector <PF_Agent*> > *collection);
	void stopTest(const Real time, map < string, vector <PF_Agent*> > *collection = NULL);

	virtual bool choose(const Real time, map < string, vector <PF_Agent*> > *collection = NULL);
	virtual void learn(const Real time, map < string, vector <PF_Agent*> > *collection = NULL, string filename = string());
	virtual void savePolicy(string name = string());
	virtual void restartPolicy(string name = string());
	virtual void resetToInitialCondition(){};

	virtual void reward(const Real time, map < string, vector <PF_Agent*> > *collection = NULL){};
	virtual void fitness(){};
	virtual void error(){};
	virtual bool mapState(vector <int> &state, map < string, vector <PF_Agent*> > *collection = NULL){ return true; };
	virtual void mapAction(int action, Real time = 0){};
	virtual void perform(const Real time = 0){};

#ifdef _RL_VIZ
	void _paintSphere(Real x, Real y, Real radius, Real r = 0.0, Real g = 0.0, Real b = 0.0) const;
	void _paintX(Real x, Real y, Real radius, Real r, Real g, Real b) const;
	void _paintTriangle(Real xc, Real yc, Real radius, Real direction, Real r, Real g, Real b) const;
	void _paintCone(Real xc, Real yc, Real radius, Real direction, Real r, Real g, Real b) const;
	void _drawCircle(Real radius, Real xc = 0.5, Real yc = 0.5, Real r = 0.0, Real g = 0.0, Real b = 0.0) const;
	void _drawFullCircle(Real radius, Real xc = 0.5, Real yc = 0.5, Real r = 0.0, Real g = 0.0, Real b = 0.0) const;
	void _rotate(Real point[2], Real angle) const;
	virtual void paint() = 0;
	virtual void paint(Real scale, const Real center[2]) = 0;
#endif

	//public attributs
	int ID;
	string name;
	bool isDead;
	bool isAveraged;

protected:
	//attributs	
	RL::RL_TabularPolicy * policystore;
	RL::RL_TabularPolicy **policy;
	enum QLearningStatus{ Waiting, Ready };
	QLearningStatus status;
	Real integralReward;
	Real learningInterval;
	Real learningTimer;
	Real fitnessValue;
	Real errorValue;
	vector<int> myState;
	int myAction;
	int myPreviousAction;
	
	bool SMOOTH;
	bool ISCONTROLLED;

	void _dist(const Real x2[2], const Real x1[2], Real d[2]) const;
	Real _angleVectors(const Real v1[2], const Real v2[2]) const;
	Real _modv(const Real v[2]) const;
	void _normalize(const Real v[2], Real n[2]) const;
	int _discretize(const Real signal, const Real minsignal, const Real maxsignal, const int nlevels, const bool deadEnd = false) const;
	int _discretizeRange(const Real value, const Real minvalue, const Real maxvalue, const int levels);
	int _discretizeAngle(Real angle, const Real sightangle, const int levels);
};

} /*namespace PF*/
