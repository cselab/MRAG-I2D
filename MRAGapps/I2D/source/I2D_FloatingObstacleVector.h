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
#include <string>
#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_FloatingObstacleOperator.h"

using namespace std;

class I2D_FloatingObstacleVector : public I2D_FloatingObstacleOperator
{
	const Real charLength, charVel;
	vector<I2D_FloatingObstacleOperator *> agents;
	VelocityBlock * arrayVelBlock;
public:
	map< string, vector<I2D_FloatingObstacleOperator *> > data;

	I2D_FloatingObstacleVector(ArgumentParser & parser, Grid<W,B>& grid, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, map< string, vector<I2D_FloatingObstacleOperator *> > _data, Real charLength, Real charVel);
	~I2D_FloatingObstacleVector(); // _data will be deallocated by this destructor 
	
	void characteristic_function();
	Real getD() const;
	Real getDrag() const;
	
	// Hydrodynamics
	void computeDragAndStuff(const Real time, const Real charLength, const Real charVel);
	void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	void computeDesiredVelocity(const double t);
	void create(const double t);
	void save(const double t, string filename = std::string());
	void restart(const double t, string filename = std::string());
	void refresh(const double t, string filename = std::string());
	vector<Real> getMass();


	// Online learning
	void savePolicy(string name=string());
	void restartPolicy(string name=string());
	bool choose(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL );
	void learn(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL, string name=string());
	void reward(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
};
