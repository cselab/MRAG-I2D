/*
 * RL_QLearning.h
 *
 *  Created on: May 26, 2011
 *      Author: mgazzola
 */

#ifndef RL_QLEARNING_H_
#define RL_QLEARNING_H_

#include <vector>
#include <string>

#include "rng.h"

#include "MRAGcore/MRAGProfiler.h"
#include "RL_TabularPolicy.h"

namespace RL
{

class RL_QLearning : public RL_TabularPolicy
{
protected:
	double LR, gamma, epsilon;
	int seed;
	unsigned long long plays;
	int freq, freqUsage;
	bool isGreedy;
	double totalIntegralReward;
	vector<double> rewards;
	vector< double > singleIntegralRewards;
	vector< vector<int> > stateActionStarts;
	vector< vector<int> > stateEnds;
	vector<unsigned long long> histPlays;
	vector<double> histRewards;
	vector<double> histUsage;

	RNG rng;
	MRAG::Profiler profiler;

	int _epsilonGreedySelection(const vector<int> & state);
	double _getValue(const vector<int> & idx);
	double _getMaxValue(const vector<int> & state);
	bool _jumpedToDifferentState(const vector<int> & stateActionStart, const vector<int> & stateEnd);

public:
	RL_QLearning(double LR, double gamma, double epsilon, int seed=0, int freq = 1);
	virtual ~RL_QLearning();

	int selectAction(const vector<int> & state);
	void setReward(double reward);
	void setStateActionStart( const vector<int> & idx );
	void setStateEnd( const vector<int> & idx );
	void update(string name = "learning");
};

}

#endif /* RL_QLEARNING_H_ */
