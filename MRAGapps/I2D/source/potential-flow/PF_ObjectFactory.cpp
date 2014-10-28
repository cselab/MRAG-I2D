/*
 *  PF_ObjectFactory.cpp
 *  DipoleCode
 *
 *  Created by Alexia on 4/4/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdlib.h>
#include <fstream>
#include <assert.h>
#include <map>
#include <vector>
#include <complex>
#include <math.h>
#include <algorithm>
#include <iostream>
#include "rng.h"
#include "RL_QLearning.h"
#include "PF_ObjectFactory.h"
#include "PF_Dipole.h"
#include "PF_DipoleInline.h"
#include "PF_SwarmCylinders.h"
#include "PF_DipoleInlineBravais.h"
#include "PF_DipoleEight.h"

using namespace PF;

PF_ObjectFactory::PF_ObjectFactory(const Real charLength, const Real XCM, const Real YCM)
		: charLength(charLength), XCM(XCM), YCM(YCM), unitVelocity(0.0)
{
}

int PF_ObjectFactory::_lines(string filename)
{
	// Open file stream
	ifstream filestream(filename.c_str());

	// Check that the filestream is correctly open
	if (!filestream.good())
	{
		printf("** Ooops: file not found. Exiting now.\n");
		abort();
	}

	// Count number of lines contained in the file
	int c = 0;
	string line;
	while (getline(filestream, line))
		c++;
	filestream.close();

	// Return number of lines
	return c;
}

void PF_ObjectFactory::maximumVelocity(Real v)
{
	const Real V = v / charLength;
	unitVelocity = (unitVelocity < V) ? V : unitVelocity;
}

Real PF_ObjectFactory::getUnitVelocity()
{
	if (unitVelocity <= 0.0)
	{
		printf("error : unitVelocity == %f\n", unitVelocity);
		abort();
	}
	return unitVelocity;
}

void PF_ObjectFactory::create(MRAG::ArgumentParser & parser, map<string, vector<PF_Agent *> >& shapesMap, const Real lr, const Real greedyEps, RL::RL_TabularPolicy ** policy)
{
	shapesMap.clear();

	// Read parser information
	Real LR = lr;
	Real GREEDYEPS = greedyEps;

	const Real GAMMA = parser("-gamma").asDouble();
	const bool SHARED = parser("-shared").asBool();
	const string factoryFile = parser("-factory").asString();
	const bool SMOOTH = parser("-smooth").asBool();
	const bool ISCONTROLLED = parser("-isControlled").asBool();
	const int learnDump = parser("-learnDump").asInt();
	const int FITSELECT = parser("-fitselection").asInt();

	assert(factoryFile != "");
	assert((LR>=0.0 && LR<=1.0));
	assert(GAMMA>=0.0 && GAMMA<=1.0);
	assert((GREEDYEPS>=0.0 && GREEDYEPS<=1.0));

	srand(time(NULL));
//	srand ( 100 );

	RNG rng(rand());

	// Count number of shapes
	const int N = _lines(factoryFile);

	if (SHARED) if (LR >= 0.0 && LR <= 1.0 && GAMMA >= 0.0 && GAMMA <= 1.0 && GREEDYEPS >= 0.0 && GREEDYEPS <= 1.0) (*policy) = new RL::RL_QLearning((double) LR, (double) GAMMA, (double) GREEDYEPS, rand(), learnDump);

	// Open factory file and retrieve information
	FILE * ppFile;
	ppFile = fopen(factoryFile.c_str(), "r");
	if (ppFile == NULL)
	{
		printf("could not open ctrl file %s!\n", factoryFile.c_str());
		abort();
	}

	printf("\n--------------PF OBJECT FACTORY: START--------------\n");
	printf("** Number of objects = %d\n", N);

	unsigned int counterID = 1;
	for (int i = 0; i < N; i++)
	{
		char nameObject[1000];
		fscanf(ppFile, "%s", nameObject);
		string name(nameObject);

		if (name == "PF_Dipole")
		{
			float variable = 0.0;
			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength;
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength;
			fscanf(ppFile, " xm=%f", &variable);
			const Real xm = XCM + variable * charLength;
			fscanf(ppFile, " ym=%f", &variable);
			const Real ym = YCM + variable * charLength;
			const complex<Real> location(xm, ym);
			fscanf(ppFile, " angle=%f", &variable);
			const Real angle = variable * M_PI;
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable;

			PF_Dipole* object = new PF_Dipole(parser, d, v, location, angle, lr, greedyEps, SMOOTH, ISCONTROLLED, T, counterID, policy);

			map<string, vector<PF_Agent*> >::iterator it = shapesMap.find(name);

			if (it == shapesMap.end()) shapesMap[name] = vector<PF_Agent*>();

			shapesMap[name].push_back(object);

			maximumVelocity(v);

			counterID++;
		}

		else if (name == "PF_Dipoles")
		{
			name = "PF_Dipole";
			float variable = 0.0;
			int variableInt = 0;
			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength;
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength;
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable;
			//assert(T>=0.01);
			fscanf(ppFile, " n=%d", &variableInt);
			const int number = variableInt;
			assert(number>=0);

			//			vector <complex <Real> > locations;
			//			locations.clear();
			const Real theta = 2 * M_PI / number;
			for (int j = 0; j < number; j++)
			{
				// set dipoles around a circle
				const Real x = XCM + 0.25 * cos(theta * j);
				const Real y = YCM + 0.25 * sin(theta * j);
				complex<Real> location(x, y);
				const Real angle = 0.5 * M_PI + theta * j;

				PF_Dipole* object = new PF_Dipole(parser, d, v, location, angle, lr, greedyEps, SMOOTH, ISCONTROLLED, T, counterID, policy);

				map<string, vector<PF_Agent*> >::iterator it = shapesMap.find(name);

				if (it == shapesMap.end()) shapesMap[name] = vector<PF_Agent*>();

				shapesMap[name].push_back(object);

				maximumVelocity(v);
				counterID++;
			}
		}
		else if (name == "PF_DipoleInline")
		{
			float variable = 0.0;
			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength;
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength;
			fscanf(ppFile, " radius=%f", &variable);
			const Real radius = variable;
			fscanf(ppFile, " acc=%f", &variable);
			const Real acc = variable;
			fscanf(ppFile, " xm=%f", &variable);
			const Real xm = XCM + variable * charLength;
			fscanf(ppFile, " ym=%f", &variable);
			const Real ym = YCM + variable * charLength;
			const complex<Real> location(xm, ym);
			Real dir[2];
			fscanf(ppFile, " dirx=%f", &variable);
			dir[0] = variable;
			fscanf(ppFile, " diry=%f", &variable);
			dir[1] = variable;
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable;

			PF_DipoleInline * object = new PF_DipoleInline(parser, d, v, radius, acc, location, dir, lr, greedyEps, SMOOTH, ISCONTROLLED, FITSELECT, T, counterID, policy, rand());

			map<string, vector<PF_Agent*> >::iterator it = shapesMap.find(name);

			if (it == shapesMap.end())
			{
				shapesMap[name] = vector<PF_Agent*>();
				maximumVelocity(v);
			}

			shapesMap[name].push_back(object);

			counterID++;
		}
		/**
		 * A single dipole that takes assigned actions at specified times
		 */
		else if (name == "PF_DipoleRandomActions")
		{
			float variable = 0.0;
			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength;
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength;
			fscanf(ppFile, " radius=%f", &variable);
			const Real radius = variable;
			fscanf(ppFile, " acc=%f", &variable);
			const Real acc = variable;
			fscanf(ppFile, " xm=%f", &variable);
			const Real xm = XCM + variable * charLength;
			fscanf(ppFile, " ym=%f", &variable);
			const Real ym = YCM + variable * charLength;
			const complex<Real> location(xm, ym);
			Real dir[2];
			fscanf(ppFile, " dirx=%f", &variable);
			dir[0] = variable;
			fscanf(ppFile, " diry=%f", &variable);
			dir[1] = variable;
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable;

			PF_DipoleInline * object = new PF_DipoleInline(parser, d, v, radius, acc, location, dir, lr, greedyEps, SMOOTH, ISCONTROLLED, FITSELECT, T, counterID, policy, rand());

			map<string, vector<PF_Agent*> >::iterator it = shapesMap.find(name);

			if (it == shapesMap.end())
			{
				shapesMap[name] = vector<PF_Agent*>();
				maximumVelocity(v);
			}

			shapesMap[name].push_back(object);

			counterID++;
		}
		/**
		 * A dipole that follows a specified path figure eight path.
		 */
		else if (name == "PF_DipoleEight")
		{
			name = "PF_DipoleEight";
			float variable = 0.0;
			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength;
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength;
			fscanf(ppFile, " radius=%f", &variable);
			const Real radius = variable;
			fscanf(ppFile, " acc=%f", &variable);
			const Real acc = variable;
			fscanf(ppFile, " xm=%f", &variable);
			const Real xm = XCM + variable * charLength;
			fscanf(ppFile, " ym=%f", &variable);
			const Real ym = YCM + variable * charLength;
			const complex<Real> location(xm, ym);
			Real dir[2];
			fscanf(ppFile, " dirx=%f", &variable);
			dir[0] = variable;
			fscanf(ppFile, " diry=%f", &variable);
			dir[1] = variable;
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable;

			PF_DipoleEight * object = new PF_DipoleEight(parser, d, v, radius, acc, location, dir, lr, greedyEps, SMOOTH, ISCONTROLLED, T, counterID, policy, rand());

			map<string, vector<PF_Agent*> >::iterator it = shapesMap.find(name);

			if (it == shapesMap.end())
			{
				shapesMap[name] = vector<PF_Agent*>();
				maximumVelocity(v);
			}

			shapesMap[name].push_back(object);

			counterID++;
		}
		/**
		 * A general 2D Bravais lattice that is extended to a circle of specified
		 * radius.
		 */
		else if (name == "PF_DipoleInlineBravais")
		{
			name = "PF_DipoleInline";
			Real variable = 0.0;
			int variableInt = 0;

			/// Properties of an individual swimmer
			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength; /// scaled swimmer 'size'
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength; /// scaled velocity
			fscanf(ppFile, " radius=%f", &variable);
			const Real radius = variable; /// turn radius over swimmer 'size'
			fscanf(ppFile, " acc=%f", &variable);
			const Real acc = variable; /// ratio defining dvel/(nominal vel)
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable; /// time between decisions

			assert(d > 0);
			assert(v > 0);
			assert(radius > 0);
			assert(acc > 0);
			assert(T > 0);

			/// Properties of the group/swarm
			fscanf(ppFile, " a1=%f", &variable);
			const Real a1 = variable * charLength;
			fscanf(ppFile, " a2=%f", &variable);
			const Real a2 = variable * charLength;
			fscanf(ppFile, " phi=%f", &variable);
			const Real phi = variable * M_PI / 180.0;
			fscanf(ppFile, " theta=%f", &variable);
			const Real theta = variable * M_PI / 180.0;
			fscanf(ppFile, " nDesired=%d", &variableInt);
			const Real nDesired = variableInt;

			assert(a1 > 0);
			assert(a2 > 0);
			assert(nDesired > 0);

			/// Create lattice
			DipoleInlineBravais swarm(a1, a2, phi, theta, nDesired);

			FILE* swarmFid;
			swarmFid = fopen("swarm.txt", "r");
			if (swarmFid == NULL)
			{
				swarm.computeRadius();
				printf("** Writing coordinates to swarm.txt!\n");
				swarm.printPositions("swarm.txt");
			}
			else /// read in previously created coordinates
			{
				fclose(swarmFid);
				printf("** Reading in coordinates (swarm.txt exists)!\n");
				swarm.readPositions("swarm.txt");
			}

			vector<double> xs = swarm.getX();
			vector<double> ys = swarm.getY();
			vector<bool> isInterior = swarm.getIsInterior();

			/// MOVE DIPOLES
			const Real scaleRatio = 1.0;
			Real const yMin = *(std::min_element(ys.begin(), ys.end()));

			for (int i = 0; i < (int) xs.size(); i++)
			{
				xs[i] = xs[i] * scaleRatio;
				ys[i] = (ys[i] - yMin) * scaleRatio;
			}

			/// ALLOCATE PF_DipoleInline SWIMMERS GIVEN XS AND YS
			for (int i = 0; i < swarm.getn(); i++)
			{
				const Real xm = XCM + xs[i];
				const Real ym = YCM + ys[i];
				assert(xm > 0 && xm < 1);
				assert(ym > 0 && ym < 1);

				const complex<Real> location(xm, ym);
				const Real dir[2] = { 0, 1 };

				PF_DipoleInline * object = new PF_DipoleInline(parser, d, v, radius, acc, location, dir, lr, greedyEps, SMOOTH, ISCONTROLLED, FITSELECT, T, counterID, policy, rand());

				if (isInterior[i] == false)
					object->isAveraged = false;

				map<string, vector<PF_Agent *> >::iterator it = shapesMap.find(name);
				if (it == shapesMap.end())
				{
					shapesMap[name] = vector<PF_Agent *>();
					maximumVelocity(v);
				}

				shapesMap[name].push_back(object);

				counterID++;
			}
			printf("** School successfully created.\n");
		}
		/**
		 * TODO: Bravais diamond school.
		 */
		else if (name == "PF_DipoleInlineBravaisDiamond")
		{
			name = "PF_DipoleInline";
			Real variable = 0.0;
			int variableInt = 0;

			/// Properties of an individual swimmer
			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength; /// scaled swimmer 'size'
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength; /// scaled velocity
			fscanf(ppFile, " radius=%f", &variable);
			const Real radius = variable; /// turn radius over swimmer 'size'
			fscanf(ppFile, " acc=%f", &variable);
			const Real acc = variable; /// ratio defining dvel/(nominal vel)
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable; /// time between decisions

			assert(d > 0);
			assert(v > 0);
			assert(radius > 0);
			assert(acc > 0);
			assert(T > 0);

			/// Properties of the group/swarm
			fscanf(ppFile, " b=%f", &variable);
			const Real b = variable * charLength;
			fscanf(ppFile, " h=%f", &variable);
			const Real h = variable * charLength;
			fscanf(ppFile, " theta=%f", &variable);
			const Real theta = variable * M_PI / 180.0;
			fscanf(ppFile, " nDesired=%d", &variableInt);
			const Real nDesired = variableInt;

			const Real phi = atan2(2 * h, b);
			const Real a1 = b;
			const Real a2 = h / sin(phi);

			assert(a1 > 0);
			assert(a2 > 0);
			assert(nDesired > 0);

			/// Create lattice
			DipoleInlineBravais swarm(a1, a2, phi, theta, nDesired);

//			FILE* swarmFid;
//			swarmFid = fopen("swarm.txt", "r");
//			if (swarmFid == NULL)
//			{
				swarm.computeRadius();
				printf("** Writing coordinates to swarm.txt!\n");
				swarm.printPositions("swarm.txt");
//			}
//			else /// read in previously created coordinates
//			{
//				fclose(swarmFid);
//				printf("** Reading in coordinates (swarm.txt exists)!\n");
//				swarm.readPositions("swarm.txt");
//			}

			vector<double> xs = swarm.getX();
			vector<double> ys = swarm.getY();
			vector<bool> isInterior = swarm.getIsInterior();

			/// MOVE DIPOLES
			const Real scaleRatio = 1.0;
			Real const yMin = *(std::min_element(ys.begin(), ys.end()));

			for (int i = 0; i < (int) xs.size(); i++)
			{
				xs[i] = xs[i] * scaleRatio;
				ys[i] = (ys[i] - yMin) * scaleRatio;
			}

			/// ALLOCATE PF_DipoleInline SWIMMERS GIVEN XS AND YS
			for (int i = 0; i < swarm.getn(); i++)
			{
				const Real xm = XCM + xs[i];
				const Real ym = YCM + ys[i];

//				assert(xm > 0 && xm < 1);
//				assert(ym > 0 && ym < 1);

				const complex<Real> location(xm, ym);
				const Real dir[2] = { 0, 1 };

				PF_DipoleInline * object = new PF_DipoleInline(parser, d, v, radius, acc, location, dir, lr, greedyEps, SMOOTH, ISCONTROLLED, FITSELECT, T, counterID, policy, rand());

				if (isInterior[i] == false) object->isAveraged = false;

				map<string, vector<PF_Agent *> >::iterator it = shapesMap.find(name);
				if (it == shapesMap.end())
				{
					shapesMap[name] = vector<PF_Agent *>();
					maximumVelocity(v);
				}

				shapesMap[name].push_back(object);

				counterID++;
			}
			printf("** School successfully created.\n");
		}
		/**
		 * TODO: Bravais square school.
		 */
		else if (name == "PF_DipoleInlineBravaisSquare")
		{
			name = "PF_DipoleInline";
			Real variable = 0.0;
			int variableInt = 0;

			/// Properties of an individual swimmer
			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength; /// scaled swimmer 'size'
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength; /// scaled velocity
			fscanf(ppFile, " radius=%f", &variable);
			const Real radius = variable; /// turn radius over swimmer 'size'
			fscanf(ppFile, " acc=%f", &variable);
			const Real acc = variable; /// ratio defining dvel/(nominal vel)
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable; /// time between decisions

			assert(d > 0);
			assert(v > 0);
			assert(radius > 0);
			assert(acc > 0);
			assert(T > 0);

			/// Properties of the group/swarm
			fscanf(ppFile, " b=%f", &variable);
			const Real b = variable * charLength;
			fscanf(ppFile, " h=%f", &variable);
			const Real h = variable * charLength;
			fscanf(ppFile, " theta=%f", &variable);
			const Real theta = variable * M_PI / 180.0;
			fscanf(ppFile, " nDesired=%d", &variableInt);
			const Real nDesired = variableInt;

			const Real phi = 90 * M_PI / 180.0;
			const Real a1 = b;
			const Real a2 = h;

			assert(a1 > 0);
			assert(a2 > 0);
			assert(nDesired > 0);

			/// Create lattice
			DipoleInlineBravais swarm(a1, a2, phi, theta, nDesired);

//			FILE* swarmFid;
//			swarmFid = fopen("swarm.txt", "r");
//			if (swarmFid == NULL)
//			{
				swarm.computeRadius();
				printf("** Writing coordinates to swarm.txt!\n");
				swarm.printPositions("swarm.txt");
//			}
//			else /// read in previously created coordinates
//			{
//				fclose(swarmFid);
//				printf("** Reading in coordinates (swarm.txt exists)!\n");
//				swarm.readPositions("swarm.txt");
//			}

			vector<double> xs = swarm.getX();
			vector<double> ys = swarm.getY();
			vector<bool> isInterior = swarm.getIsInterior();

			/// MOVE DIPOLES
			const Real scaleRatio = 1.0;
			Real const yMin = *(std::min_element(ys.begin(), ys.end()));

			for (int i = 0; i < (int) xs.size(); i++)
			{
				xs[i] = xs[i] * scaleRatio;
				ys[i] = (ys[i] - yMin) * scaleRatio;
				std::cout << ys[i];
			}

			/// ALLOCATE PF_DipoleInline SWIMMERS GIVEN XS AND YS
			for (int i = 0; i < swarm.getn(); i++)
			{
				const Real xm = XCM + xs[i];
				const Real ym = YCM + ys[i];

//				assert(xm >= 0 && xm <= 1);
//				assert(ym >= 0 && ym <= 1);

				const complex<Real> location(xm, ym);
				const Real dir[2] = { 0, 1 };

				PF_DipoleInline * object = new PF_DipoleInline(parser, d, v, radius, acc, location, dir, lr, greedyEps, SMOOTH, ISCONTROLLED, FITSELECT, T, counterID, policy, rand());

				if (isInterior[i] == false) object->isAveraged = false;

				map<string, vector<PF_Agent *> >::iterator it = shapesMap.find(name);
				if (it == shapesMap.end())
				{
					shapesMap[name] = vector<PF_Agent *>();
					maximumVelocity(v);
				}

				shapesMap[name].push_back(object);

				counterID++;
			}
			printf("** School successfully created.\n");
		}
		/**
		 * TODO: Bravais equilateral lattice.
		 */
		else if (name == "PF_DipoleInlineBravaisHexagon")
		{
			name = "PF_DipoleInline";
			Real variable = 0.0;
			int variableInt = 0;

			/// Properties of an individual swimmer
			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength; /// scaled swimmer 'size'
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength; /// scaled velocity
			fscanf(ppFile, " radius=%f", &variable);
			const Real radius = variable; /// turn radius over swimmer 'size'
			fscanf(ppFile, " acc=%f", &variable);
			const Real acc = variable; /// ratio defining dvel/(nominal vel)
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable; /// time between decisions

			assert(d > 0);
			assert(v > 0);
			assert(radius > 0);
			assert(acc > 0);
			assert(T > 0);

			/// Properties of the group/swarm
			fscanf(ppFile, " b=%f", &variable);
			const Real b = variable * charLength;
			fscanf(ppFile, " theta=%f", &variable);
			const Real theta = variable * M_PI / 180.0;
			fscanf(ppFile, " nDesired=%d", &variableInt);
			const Real nDesired = variableInt;

			const Real phi = 60 * M_PI / 180.0;
			const Real a1 = b;
			const Real a2 = b;

			assert(a1 > 0);
			assert(a2 > 0);
			assert(nDesired > 0);

			/// Create lattice
			DipoleInlineBravais swarm(a1, a2, phi, theta, nDesired);

//			FILE* swarmFid;
//			swarmFid = fopen("swarm.txt", "r");
//			if (swarmFid == NULL)
//			{
				swarm.computeRadius();
				printf("** Writing coordinates to swarm.txt!\n");
				swarm.printPositions("swarm.txt");
//			}
//			else /// read in previously created coordinates
//			{
//				fclose(swarmFid);
//				printf("** Reading in coordinates (swarm.txt exists)!\n");
//				swarm.readPositions("swarm.txt");
//			}

			vector<double> xs = swarm.getX();
			vector<double> ys = swarm.getY();
			vector<bool> isInterior = swarm.getIsInterior();

			/// MOVE DIPOLES
			const Real scaleRatio = 1.0;
			Real const yMin = *(std::min_element(ys.begin(), ys.end()));

			for (int i = 0; i < (int) xs.size(); i++)
			{
				xs[i] = xs[i] * scaleRatio;
				ys[i] = (ys[i] - yMin) * scaleRatio;
			}

			/// ALLOCATE PF_DipoleInline SWIMMERS GIVEN XS AND YS
			for (int i = 0; i < swarm.getn(); i++)
			{
				const Real xm = XCM + xs[i];
				const Real ym = YCM + ys[i];
				assert(xm >= 0 && xm <= 1);
				assert(ym >= 0 && ym <= 1);

				const complex<Real> location(xm, ym);
				const Real dir[2] = { 0, 1 };

				PF_DipoleInline * object = new PF_DipoleInline(parser, d, v, radius, acc, location, dir, lr, greedyEps, SMOOTH, ISCONTROLLED, FITSELECT, T, counterID, policy, rand());

				if (isInterior[i] == false) object->isAveraged = false;

				map<string, vector<PF_Agent *> >::iterator it = shapesMap.find(name);
				if (it == shapesMap.end())
				{
					shapesMap[name] = vector<PF_Agent *>();
					maximumVelocity(v);
				}

				shapesMap[name].push_back(object);

				counterID++;
			}
			printf("** School successfully created.\n");
		}
		/**
		 * This is a model swarm based on repulsive far field interactions between
		 * a SPLINE shaped boundary and N objects distributed inside the boundary.
		 * See Diego Rossinelli's thesis for more details.
		 *
		 * Please specify the details for an individual dipole swimmer and the
		 * parameters associated with the spline shape (area, k1, k2, and alpha as
		 * presented in Rossinelli's thesis.
		 */
		else if (name == "PF_SwarmCylinders")
		{
			name = "PF_DipoleInline";
			Real variable = 0.0;
			int variableInt = 0;

			/// Properties of an individual swimmer

			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength; /// scaled swimmer 'size'
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength; /// scaled velocity
			fscanf(ppFile, " radius=%f", &variable);
			const Real radius = variable; /// turn radius over swimmer 'size'
			fscanf(ppFile, " acc=%f", &variable);
			const Real acc = variable; /// ratio defining dvel/(nominal vel)
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable; /// time between decisions

			assert(d > 0);
			assert(v > 0);
			assert(radius > 0);
			assert(acc > 0);
			assert(T > 0);

			/// Properties of the group/swarm

			fscanf(ppFile, " n=%d", &variableInt);
			const int nd = variableInt; /// number of individuals
			fscanf(ppFile, " avgDist=%f", &variable);
			const Real avgDistToDipole = variable * charLength; /// area enclosed by spline
			fscanf(ppFile, " k1=%f", &variable);
			const Real k1 = variable; /// spline parameter
			fscanf(ppFile, " k2=%f", &variable);
			const Real k2 = variable; /// spline parameter
			fscanf(ppFile, " alpha=%f", &variable);
			const Real alpha = variable; /// spline paremter

			assert(nd > 0);
			assert(k1 > 0);
			assert(k2 > 0);
			assert(alpha > 0 && alpha < M_PI);

			double rmin = 0.05; /// only used in initial placement
			FILE* fid;
			fid = fopen("parameters.txt", "w");
			fprintf(fid, "%e", radius);
			fprintf(fid, "%e", rmin);
			fclose(fid);

			SwarmCylinders swarm(nd, 0.5 * d, rmin, 1.0, k1, k2, alpha);
			FILE* swarmFid;
			swarmFid = fopen("swarm.txt", "r");
			if (swarmFid == NULL)
			{
				printf("** Calculating new coordinates!\n");
				swarm.placeCylinders();
				swarm.printSpline();
				swarm.findEquilibrium();
				swarm.rotateClockwiseQuarterTurn();

				printf("** Writing coordinates to swarm.txt!\n");
				swarm.printPositions("swarm.txt");
			}
			else /// read in previously created coordinates
			{
				fclose(swarmFid);
				printf("** Reading in coordinates (swarm.txt exists)!\n");
				swarm.readPositions("swarm.txt");
			}

			/// CALCULATE EQUILIBRIUM POSITIONS
			vector<double> xs;
			vector<double> ys;

			xs = swarm.getX();
			ys = swarm.getY();

			/// SCALE DIPOLES AND MOVE DIPOLES
			Real const areaPerDipole = avgDistToDipole * avgDistToDipole;
			Real const totalSwarmArea = areaPerDipole * Real(nd);
			Real const scaleRatio = sqrt(totalSwarmArea);

			Real const yMin = *(std::min_element(ys.begin(), ys.end()));

			for (int i = 0; i < (int) xs.size(); i++)
			{
				xs[i] = xs[i] * scaleRatio;
				ys[i] = (ys[i] - yMin) * scaleRatio;
			}

			/// ALLOCATE PF_DipoleInline SWIMMERS GIVEN XS AND YS
			for (int i = 0; i < nd; i++)
			{
				const Real xm = XCM + xs[i];
				const Real ym = YCM + ys[i];
				assert(xm > 0 && xm < 1);
				assert(ym > 0 && ym < 1);

				const complex<Real> location(xm, ym);
				const Real dir[2] = { 0, 1 };

				PF_DipoleInline * object = new PF_DipoleInline(parser, d, v, radius, acc, location, dir, lr, greedyEps, SMOOTH, ISCONTROLLED, FITSELECT, T, counterID, policy, rand());

				map<string, vector<PF_Agent *> >::iterator it = shapesMap.find(name);
				if (it == shapesMap.end())
				{
					shapesMap[name] = vector<PF_Agent *>();
					maximumVelocity(v);
				}

				shapesMap[name].push_back(object);

				counterID++;
			}
			printf("** School successfully created.\n");

		}
		else if (name == "PF_DipoleInlineSquaredSchool")
		{
			name = "PF_DipoleInline";
			float variable = 0.0;
			int variableInt = 0;
			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength;
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength;
			fscanf(ppFile, " radius=%f", &variable);
			const Real radius = variable;
			fscanf(ppFile, " acc=%f", &variable);
			const Real acc = variable;
			fscanf(ppFile, " deltax=%f", &variable);
			const Real deltax = variable * charLength;
			fscanf(ppFile, " deltay=%f", &variable);
			const Real deltay = variable * charLength;
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable;
			fscanf(ppFile, " n=%d", &variableInt);
			const int number = variableInt;
			assert(number>=0);

			int side = (int) sqrt(number);
			const Real xcorrect = 0.5 * (side - 1) * deltax; //to center the school according to the axe x=0.5
			for (int ix = 0; ix < side; ix++)
				for (int iy = 0; iy < side; iy++)
				{
					const Real xm = XCM - xcorrect + ix * deltax;
					const Real ym = YCM + iy * deltay;
					const complex<Real> location(xm, ym);
					const Real dir[2] = { 0, 1 };

					PF_DipoleInline * object = new PF_DipoleInline(parser, d, v, radius, acc, location, dir, lr, greedyEps, SMOOTH, ISCONTROLLED, FITSELECT, T, counterID, policy, rand());

					map<string, vector<PF_Agent *> >::iterator it = shapesMap.find(name);
					if (it == shapesMap.end())
					{
						shapesMap[name] = vector<PF_Agent *>();
						maximumVelocity(v);
					}

					shapesMap[name].push_back(object);

					counterID++;
				}
		}
		else if (name == "PF_DipoleInlineDiamondSchool")
		{
			name = "PF_DipoleInline";
			float variable = 0.0;
			int variableInt = 0;
			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength;
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength;
			fscanf(ppFile, " radius=%f", &variable);
			const Real radius = variable;
			fscanf(ppFile, " acc=%f", &variable);
			const Real acc = variable;
			fscanf(ppFile, " deltax=%f", &variable);
			const Real deltax = variable * charLength;
			fscanf(ppFile, " deltay=%f", &variable);
			const Real deltay = variable * charLength;
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable;
			fscanf(ppFile, " n=%d", &variableInt);
			const int number = variableInt;
			assert(number>=0);

			const Real fx = deltax * cos(0.25 * M_PI);
			const Real fy = deltay * cos(0.25 * M_PI);
			int side = (int) sqrt(number);
			for (int ix = 0; ix < side; ix++)
				for (int iy = 0; iy < side; iy++)
				{
					const Real xm = XCM + ix * fx - iy * fy;
					const Real ym = YCM + iy * fy + ix * fx;
					const complex<Real> location(xm, ym);
					const Real dir[2] = { 0, 1 };

					PF_DipoleInline * object = new PF_DipoleInline(parser, d, v, radius, acc, location, dir, lr, greedyEps, SMOOTH, ISCONTROLLED, FITSELECT, T, counterID, policy, rand());

					map<string, vector<PF_Agent *> >::iterator it = shapesMap.find(name);

					if (it == shapesMap.end())
					{
						shapesMap[name] = vector<PF_Agent *>();
						maximumVelocity(v);
					}

					shapesMap[name].push_back(object);

					counterID++;
				}
		}
		else if (name == "PF_DipoleInlineHexagonalSchool")
		{
			name = "PF_DipoleInline";
			float variable = 0.0;
			int variableInt = 0;
			fscanf(ppFile, " d=%f", &variable);
			const Real d = variable * charLength;
			fscanf(ppFile, " v=%f", &variable);
			const Real v = variable * charLength;
			fscanf(ppFile, " radius=%f", &variable);
			const Real radius = variable;
			fscanf(ppFile, " acc=%f", &variable);
			const Real acc = variable;
			fscanf(ppFile, " delta=%f", &variable);
			const Real delta = variable * charLength;
			Real dir[2];
			fscanf(ppFile, " dirx=%f", &variable);
			dir[0] = variable;
			fscanf(ppFile, " diry=%f", &variable);
			dir[1] = variable;
			fscanf(ppFile, " T=%f", &variable);
			const Real T = variable;
			fscanf(ppFile, " n=%d", &variableInt);
			const int number = variableInt;
			assert(number>=0);

			//			int side = (int)sqrt(number);
			vector<complex<Real> > locations;
			locations.clear();
			const Real theta = M_PI / 3.;
			Real distCenter = 2 * delta * cos(M_PI / 6.);
			Real center[2];
			pair<Real, Real> prCenter(XCM, YCM);
			vector<pair<Real, Real> > centers;
			centers.clear();
			centers.push_back(prCenter);
			vector<pair<Real, Real> > allCenters;
			allCenters.clear();
			allCenters.push_back(prCenter);
			int ccl = 0;
			int nbDipole = 0;
			const Real error = 1e-6;
			while (nbDipole < number)
			{
				Real length = (int) centers.size();

				for (int tr = 0; tr < length; tr++)
				{
					center[0] = centers[tr].first;
					center[1] = centers[tr].second;

					for (int hexa = 0; hexa < 6; hexa++)
					{
						const Real xm = center[0] + delta * sin(theta * hexa);
						const Real ym = center[1] + delta * cos(theta * hexa);
						const complex<Real> location(xm, ym);

						bool same = false;
						Real realLoc = real(location);
						Real imagLoc = imag(location);
						for (int k = 0; k < (int) locations.size(); k++)
						{
							Real realLocK = real(locations[k]);
							Real imagLocK = imag(locations[k]);
							if (realLocK > realLoc - error && realLocK < realLoc + error && imagLocK > imagLoc - error && imagLocK < imagLoc + error)
							{
								same = true;
								break;
							}
						}

						if (!same)
						{
							locations.push_back(location);
							nbDipole++;
						}

						if (nbDipole == number) break;
					}
					if (nbDipole == number) break;
				}

				ccl++;
				vector<pair<Real, Real> > interCenters;
				interCenters.clear();
				interCenters = centers;
				centers.clear();
				for (int gs = 0; gs < length; gs++)
					for (int h = 0; h < ccl * 6; h++)
					{
						Real cx = interCenters[gs].first;
						Real cy = interCenters[gs].second;
						center[0] = cx + distCenter * cos(theta * h);
						center[1] = cy + distCenter * sin(theta * h);

						bool same = false;
						for (int k = 0; k < (int) allCenters.size(); k++)
							if (allCenters[k].first > center[0] - error && allCenters[k].first < center[0] + error && allCenters[k].second > center[1] - error && allCenters[k].second < center[1] + error)
							{
								same = true;
								break;
							}

						if (!same)
						{
							prCenter.first = center[0];
							prCenter.second = center[1];
							centers.push_back(prCenter);
							allCenters.push_back(prCenter);
						}

					}
			}

			for (int j = 0; j < (int) locations.size(); j++)
			{
				const complex<Real> location = locations[j];

				PF_DipoleInline * object = new PF_DipoleInline(parser, d, v, radius, acc, location, dir, lr, greedyEps, SMOOTH, ISCONTROLLED, FITSELECT, T, counterID, policy, rand());

				map<string, vector<PF_Agent *> >::iterator it = shapesMap.find(name);

				if (it == shapesMap.end())
				{
					shapesMap[name] = vector<PF_Agent *>();
					maximumVelocity(v);
				}

				shapesMap[name].push_back(object);

				counterID++;
			}
		}
		else
		{
			printf("Case not defined!! One of your shapes' name in the factory file is wrong!\n");
			abort();
		}

	}

	fclose(ppFile);

	printf("--------------PF OBJECT FACTORY: END----------------\n");
}
