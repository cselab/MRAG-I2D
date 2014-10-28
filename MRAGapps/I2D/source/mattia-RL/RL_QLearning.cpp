/*
 * RL_QLearning.cpp
 *
 *  Created on: May 26, 2011
 *      Author: mgazzola
 */

#include <iostream>
#include <assert.h>
#include <fstream>
#include "RL_QLearning.h"

namespace RL
{

RL_QLearning::RL_QLearning(double LR, double gamma, double epsilon, int seed, int freq): rng(seed), LR(LR), gamma(gamma), epsilon(epsilon), plays(0), freq(freq), totalIntegralReward(0.0)
{
	assert( LR>=0.0 && LR<=1.0);
	assert( gamma>=0.0 && gamma<=1.0 );
	assert( epsilon>=0.0 && epsilon<=1.0 );
	assert( (LR+gamma)<=1.0 );

	if(LR>0.0 && LR<=1.0 && gamma>=0.0 && gamma<=1.0 && epsilon>=0.0 && epsilon<=1.0 && (LR+gamma)>1.0)
	{
		printf("RL parameters are fucked up!\n");
		abort();
	}
}

RL_QLearning::~RL_QLearning()
{
}

int RL_QLearning::_epsilonGreedySelection(const vector<int> & state)
{
	const int DIM = dim.size();
	assert( (int)state.size() == (DIM-1) );
	vector<int> idx(state);
	idx.push_back(0.0);

	double maxQ = Q.read(idx);
	int maxQIndex = 0;

	if( rng.uniform(0.0,1.0) < epsilon )
		maxQIndex = floor( rng.uniform( 0.0, (double)dim[DIM-1] ) );
	else
	{
		//vector< double > qs;
		//vector< int > qsIdx;

		for(int k=0; k < dim[DIM-1]; k++)
		{
			idx[DIM-1] = k;
			const double qread = Q.read(idx);
			//qs.push_back(qread);
			maxQIndex = (qread>maxQ)?k:maxQIndex;
			maxQ = max(maxQ, qread);
		}

		//for(int k=0; k < dim[DIM-1]; k++)
		//{
		//	idx[DIM-1] = k;
		//	if(qs[k]==maxQ)
		//		qsIdx.push_back(k);
		//}


		//printf("shit=%d\n",(int)floor(rng.uniform(0.0,(double)qsIdx.size())));
		//return qsIdx[(int)floor(rng.uniform(0.0,(double)qsIdx.size()))];
	}

	return maxQIndex;
}

int RL_QLearning::selectAction(const vector<int> & state)
{
	assert( (int)state.size() == ((int)dim.size()-1) );
	return _epsilonGreedySelection(state);
}

double RL_QLearning::_getValue(const vector<int> & idx)
{
	return Q.read(idx);
}

double RL_QLearning::_getMaxValue(const vector<int> & state)
{
	const int DIM = dim.size();
	assert( (int)state.size() == (DIM-1) );
	vector<int> idx(state);
	idx.push_back(0.0);

	double maxQ = Q.read(idx);

	for(int k=0; k < dim[DIM-1]; k++)
	{
		idx[DIM-1] = k;
		maxQ = max(maxQ, Q.read(idx));
	}

	return maxQ;
}

void RL_QLearning::setReward(double reward)
{
	rewards.push_back(reward);
}

void RL_QLearning::setStateActionStart( const vector<int> & idx )
{
	assert( idx.size() == dim.size() );
	Qusage(idx) += 1;
	stateActionStarts.push_back(idx);
}

void RL_QLearning::setStateEnd( const vector<int> & idx )
{
	assert( idx.size() == (dim.size()-1) );
	stateEnds.push_back(idx);
}

bool RL_QLearning::_jumpedToDifferentState(const vector<int> & stateActionStart, const vector<int> & stateEnd)
{
	const int sizeSA = stateActionStart.size();
	const int sizeSE = stateEnd.size();
	assert(sizeSA==(sizeSE+1));

	bool same = true;
	for( unsigned int i=0; i<stateEnd.size(); i++)
		same *= (stateActionStart[i]==stateEnd[i]);

	return same;
}

void RL_QLearning::update(string name)
{
	const unsigned int rsize = rewards.size();
	assert( rsize == stateActionStarts.size() );
	assert( stateEnds.size() == stateActionStarts.size() );

	if( rsize==0 )
		return;

	// Calculate delta state action value for update
	vector<double> deltaStateActionValues(rsize,0);
	for(unsigned int i=0; i<rsize; ++i)
	{
		const double reward = rewards[i];
		assert(rewards[i]==rewards[i]);
		const double Qsa = _getValue( stateActionStarts[i] );
		const double QsaMax = _getMaxValue( stateEnds[i] );
		deltaStateActionValues[i] = reward + gamma*QsaMax - Qsa;
		assert(deltaStateActionValues[i]==deltaStateActionValues[i]);

		totalIntegralReward += reward;
	}

	// Update state action values
	for(unsigned int i=0; i<rsize; i++)
		Q(stateActionStarts[i]) += LR * deltaStateActionValues[i];

	// Dump statistics
	histPlays.push_back(plays);
	histRewards.push_back(totalIntegralReward);
	histUsage.push_back(Q.usage());

	if(histPlays.size()%freq==0)
	{
		ofstream out(name.c_str(), ios_base::app);

		for(unsigned int i=0; i<histRewards.size(); ++i)
			out << histPlays[i] << "  " << histRewards[i] << "  " << histUsage[i] << endl;

		out.flush();
		out.close();

		histPlays.clear();
		histRewards.clear();
		histUsage.clear();
	}

	// Clear vectors
	rewards.clear();
	stateActionStarts.clear();
	stateEnds.clear();

	// Clean up integral quantities
	assert(totalIntegralReward==totalIntegralReward);
	totalIntegralReward = 0.0;

	// Update counter
	plays += rsize;
}

}
