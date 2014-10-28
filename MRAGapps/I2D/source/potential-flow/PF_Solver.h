/*
 *  PF_Solver.h
 *  DipoleCode
 *
 *  Created by Alexia on 3/19/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <complex>
#include <deque>
#include "IF2D_Test.h"
#include "RL_Environment.h"
#include "RL_TabularPolicy.h"
#include "PF_Agent.h"
//#include "PF_Tracer.h"

namespace PF
{

class PF_Solver : public IF2D_Test
{
protected:

	enum
	{
		learning, averaging
	} solverState;

	Real CHARLENGTH, XPOS, YPOS;
	bool RESTART;
	int SAVEFREQ;
	int FITNESSSAVEFREQ; /// number of iterations before dumping a fitness file
	bool ISCONTROLLED; /// if any of the agents are controlled by learning
	bool SMOOTH;
	Real LEARNINGTIME; /// amount of time given to all agents to initially learn a 'good' policy
	Real FITNESSTIME; /// amount of time given to average fitness function over
	Real FITNESSBUFFER; /// take (FINTESSBUFFER - 1)% longer to get rid of transient after a refresh
	int NAVG; /// number of final FITNESSTIME runs to average over
	int NLEARNINGLEVELS;
	bool INDIVIDUALFITNESS;
	Real LR;
	Real GREEDYEPS;
	int FITSELECT;
	bool ISLABFRAME;
	bool ISUSINGTRACERS;

	MRAG::ArgumentParser parser;
	PF_Agent * collection;
	PF_Agent * tracers;
	RL::RL_TabularPolicy *policy;
	MRAG::Profiler profiler;
	
	bool needsRefreshing;
	bool isUsingSoftReset;
	int nTimesLearned;
	int nTimesAveraged;
	int nCompletedAveraging;
	int nTeribleRefreshes;
	int nAveragingScrewUps;
	int crashPenalty; // a measure of robustness of the formation + policy
	double totalLearningTime;
	Real unitVelocity;
	Real averagedSchoolFitness;
	Real timeAveragedSchoolFitness;
	Real currentFitness;
	Real schoolError;
	Real penaltyForBeingTooClose;
	deque<vector<Real> > individualFitnessesInWindow;
	deque<Real> schoolFitnessInWindow;
	deque<double> timeInWindow;
	deque<double> timeStepsInWindow;

	vector<Real> currentFitnesses;
	vector<vector<double> > currentIndividualFitnesses;

	Real schoolCenter[2];
	Real schoolEdgeLength;
	vector<pair<Real, Real> > schoolPoints;
	vector<double> avgFitness;
	vector<double> currentIndividualFitness;

	void _dispose();
	void _refresh();
	void _softRefresh();
	void _prepareAgents();
	void _prepareTracers();
	void _computeVelocity();
	void _computeVelocityTracers();
	double _setTimeStep();
	void _save(double time);
	void _saveResetTime(double time);

	void _computeFitness(double time, double dt);
	void _computeTimeAveragedFitness(double time);
	void _computeError(double time);
	void _checkStatus(double time, double & timeSinceValidRefresh);
	void _checkStatusPostAction(double time, double dt, bool valid, double & timeAtLastValidRefresh, double & timeAtLastNonValidRefresh);
	void _checkSettings();
	Real _computeTooClosePenalty();

	void _getTargetSchoolInfo(Real x[2], Real& edgeLength);

#ifdef _RL_VIZ
	void _paint(const double time);
	void _paintSquareCanvas(double edgeLength, double bottomLeftCorner[2], double bgColor[3], bool isOutlinedBlack);
	void _paintBox(Real edgeLength, Real bottomLeftCorner[2], Real bgColor[3]);
	void _paintText(float x, float y, const char* text, float r, float g, float b, float a);
	void _paintSphere(Real x, Real y, Real radius, Real r, Real g, Real b) const;
#endif

public:
	PF_Solver(int argc, const char ** argv);
	virtual ~PF_Solver();

	void run();
	void paint(){};
};

} /*namespace PF*/
