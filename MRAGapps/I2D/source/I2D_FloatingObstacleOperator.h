/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"
#include "I2D_PenalizationOperator.h"
#include "RL_TabularPolicy.h"
#include "RL_QLearning.h"

using namespace std;

class I2D_FloatingObstacleOperator : public I2D_ObstacleOperator
{		
protected:
	// Hydrodynamics
	Real eps, D, Cd, dimT;
	Real Uinf[2];

	Grid<W,B>& grid;
	BlockProcessing block_processing;

	map<int, const VelocityBlock *> desired_velocity;
	I2D_PenalizationOperator& penalization;

	// Online learning
	enum QLearningStatus{ Waiting, Ready };
	QLearningStatus status;
	RL::RL_TabularPolicy ** policy;
	RL::RL_TabularPolicy * policystore;
	double integralReward;
	double learningInterval;
	double learningTimer;

	// Gives back the vector between x2 and x1   x2 * <------- * x1
	void _dist(const double x2[2], const double x1[2], double d[2]) const;

	// Gives back the angle between v1 and v2
	// angle(v1,v2) = -angle(v2,v1)
	//											v1
	//                                        /
	// 										 /____ v2
	double _angleVectors(const double v1[2], const double v2[2]) const;

	//  x1 *------> dir
	//			|
	//			|
	//			*x2
	bool _isRight(const double x1[2], const double dir[2], const double x2[2]) const;

	int _discretizeRange(const double value, const double minvalue, const double maxvalue, const int levels);

public:
	I2D_FloatingObstacleOperator(ArgumentParser & parser, Grid<W,B>& grid, const Real D, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int _ID = 0, RL::RL_TabularPolicy ** _policy = NULL, const int seed = 0) :
		grid(grid), eps(eps), D(D), Cd(0.0), dimT(0.0), penalization(penalization), ID(_ID), policy(_policy), policystore(NULL), integralReward(0.0), status(Ready), learningInterval(0.0), learningTimer(0.0)
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];

		if(policy!=NULL)
			if((*policy)==NULL)
			{
				const Real LR = parser("-lr").asDouble();
				const Real GAMMA = parser("-gamma").asDouble();
				const Real GREEDYEPS = parser("-greedyEps").asDouble();

				assert( LR>=0.0 && LR<=1.0);
				assert( GAMMA>=0.0 && GAMMA<=1.0 );
				assert( GREEDYEPS>=0.0 && GREEDYEPS<=1.0 );
				assert( (LR+GAMMA)<=1.0 );

				if(LR>0.0 && LR<=1.0 && GAMMA>=0.0 && GAMMA<=1.0 && GREEDYEPS>=0.0 && GREEDYEPS<=1.0 && (LR+GAMMA)>1.0)
				{
					printf("RL parameters are fucked up!\n");
					abort();
				}

				srand ( time(NULL) );

				policystore = new RL::RL_QLearning(LR,GAMMA,GREEDYEPS,seed);
				policy = &policystore;
				assert(policy!=NULL);
			}
	}

	virtual ~I2D_FloatingObstacleOperator()
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
			}
	}

	// Hydrodynamics
	virtual void computeDragAndStuff(const Real time, const Real charLength = 1.0, const Real charVel = 1.0);
	virtual void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL) = 0;
	virtual void computeDesiredVelocity(const double t) = 0;
	virtual void create(const double t) = 0;
	virtual void save(const double t, string filename = std::string()) = 0;
	virtual void restart(const double t, string filename = std::string()) = 0;
	virtual void refresh(const double t, string filename = std::string()) = 0;
	map<int , const VelocityBlock *> getDesiredVelocity(){ return desired_velocity; }
	virtual vector<Real> getMass();

	// Online learning
	enum Label { FISH };
	int ID;
	virtual void savePolicy(string name=string());
	virtual void restartPolicy(string name=string());
	virtual bool choose(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	virtual void stopTest(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	virtual void learn(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL,string name=string());
	virtual void reward(const double t, map< string, vector<I2D_FloatingObstacleOperator* > > * _data = NULL){};
	virtual void mapAction(int action){};
	virtual bool mapState(vector<int> & state, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL ){ return true; };

	// accessor
	virtual Real getDrag() const { return this->Cd; }
};
