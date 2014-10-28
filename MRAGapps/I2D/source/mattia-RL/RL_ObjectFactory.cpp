/*
 * RL_ObjectFactory.cpp
 *
 *  Created on: Mar 21, 2012
 *      Author: mgazzola
 */

#include <stdlib.h>
#include <fstream>
#include <assert.h>
#include "rng.h"
#include "RL_ObjectFactory.h"
#include "RL_QLearning.h"
//#include "RL_SmartyCircle.h"
//#include "RL_Column.h"
//#include "RL_DynamicColumn.h"
//#include "RL_SmartyDodger.h"
//#include "RL_CircularWall.h"
#include "RL_Food.h"
#include "RL_SmartyGlutton.h"
//#include "RL_SmartyInline.h"
//#include "RL_SmartyLattice.h"

namespace RL
{

RL_ObjectFactory::RL_ObjectFactory(const Real charLength, const Real XCM, const Real YCM) : charLength(charLength), XCM(XCM), YCM(YCM)
{
}

RL_ObjectFactory::~RL_ObjectFactory()
{
}

int RL_ObjectFactory::_lines(string filename)
{
	// Open file stream
	ifstream filestream(filename.c_str());

	// Check that the filestream is correctly open
	if(!filestream.good())
	{
		printf("ooops: file not found. Exiting now.\n");
		abort();
	}

	// Count number of lines contained in the file
	int c = 0;
	string line;
	while( getline(filestream, line) ) c++;
	filestream.close();

	// Return number of lines
	return c;
}

void RL_ObjectFactory::create(MRAG::ArgumentParser & parser, map< string, vector<RL_Agent *> >& shapesMap, RL_TabularPolicy ** policy)
{
	shapesMap.clear();

	// Read parser information
	const Real LR = parser("-lr").asDouble();
	const Real GAMMA = parser("-gamma").asDouble();
	const Real GREEDYEPS = parser("-greedyEps").asDouble();
	const bool SHARED = parser("-shared").asBool();
	const string factoryFile = parser("-factory").asString();
	const int learnDump = parser("-learnDump").asInt();

	assert(factoryFile != "");
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
	//printf("seed %d!\n", rand());
	//exit(0);

	RNG rng( rand() );

	// Count number of shapes
	const int N = _lines(factoryFile);

	if(SHARED)
		(*policy) = new RL_QLearning(LR,GAMMA,GREEDYEPS,rand(),learnDump);

	// Open factory file and retrieve information
	FILE * ppFile;
	ppFile = fopen(factoryFile.c_str(),"r");
	if(ppFile==NULL){ printf("could not open ctrl file %s!\n", factoryFile.c_str()); abort(); }

	printf("\n--------------RL OBJECT FACTORY: START--------------\n");
	printf("number of objects = %d\n",N);

	unsigned int counterID = 1;
	for(int i=0; i<N; i++)
	{
		char nameObject[1000];
		fscanf(ppFile,"%s",nameObject);
		string name(nameObject);

		/*
		if( name == "RL_CircularWall" )
		{
			RL_CircularWall * object = new RL_CircularWall(parser);

			map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<RL_Agent *>();

			shapesMap[name].push_back(object);

			counterID++;
		}
		/*
		else if( name == "RL_SmartyCircle" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;

			RL_SmartyCircle * object = new RL_SmartyCircle(parser, xm, ym, d, counterID, policy, rand());

			map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<RL_Agent *>();

			shapesMap[name].push_back(object);

			counterID++;
		}
		*/
		/*
		else if( name == "RL_SmartyInline" )
		{
			Real dir[2];
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," dirx=%f",&variable);
			dir[0] = variable;
			fscanf(ppFile," diry=%f",&variable);
			dir[1] = variable;

			RL_SmartyInline * object = new RL_SmartyInline(parser, xm, ym, d, dir, counterID, policy, rand());

			map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<RL_Agent *>();

			shapesMap[name].push_back(object);

			counterID++;
		}
		*/
		/*
		else if( name == "RL_SmartyLattice" )
		{
			Real dir[2];
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," dirx=%f",&variable);
			dir[0] = variable;
			fscanf(ppFile," diry=%f",&variable);
			dir[1] = variable;

			RL_SmartyLattice * object = new RL_SmartyLattice(parser, xm, ym, d, dir, counterID, policy, rand());

			map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<RL_Agent *>();

			shapesMap[name].push_back(object);

			counterID++;
		}
		*/
		/*
		else if( name == "RL_SmartyLattices" )
		{
			name = "RL_SmartyLattice";
			Real dir[2];
			int variableInt = 0.0;
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," dirx=%f",&variable);
			dir[0] = variable;
			fscanf(ppFile," diry=%f",&variable);
			dir[1] = variable;
			fscanf(ppFile," n=%d",&variableInt);
			const int number = variableInt;
			assert(number>=0);

			int side = (int)sqrt(number);
			for(int ix=0; ix<side; ix++)
				for(int iy=0; iy<side; iy++)
				{
					const Real xm = XCM + ix*charLength;
					const Real ym = YCM + iy*charLength;

					RL_SmartyLattice * object = new RL_SmartyLattice(parser, xm, ym, d, dir, counterID, policy, rand());

					map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

					if(it==shapesMap.end())
						shapesMap[name] = vector<RL_Agent *>();

					shapesMap[name].push_back(object);

					counterID++;
				}
		}
		*/
		/*
		else if( name == "RL_SmartyDodger" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;
			assert(T>=0.01);

			RL_SmartyDodger * object = new RL_SmartyDodger(parser, xm, ym, d, T, counterID, policy, rand());

			map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<RL_Agent *>();

			shapesMap[name].push_back(object);

			counterID++;
		}*/
		/*
		else if( name == "RL_SmartyDodgers" )
		{
			name = "RL_SmartyDodger";
			float variable = 0.0;
			int variableInt = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;
			assert(T>=0.01);
			fscanf(ppFile," n=%d",&variableInt);
			const int number = variableInt;
			assert(number>=0);

			for(int j=0; j<number; j++)
			{
				const Real radius = rng.uniform(0.0,0.3);
				const Real angle = rng.uniform(0.0,2*M_PI);
				const Real xx = radius*cos(angle);
				const Real yy = radius*sin(angle);

				RL_SmartyDodger * object = new RL_SmartyDodger(parser, xx+0.5, yy+0.5, d, T, counterID, policy, rand());

				map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

				if(it==shapesMap.end())
					shapesMap[name] = vector<RL_Agent *>();

				shapesMap[name].push_back(object);

				counterID++;
			}
		}*/
		/*
		else if( name == "RL_Column" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;

			RL_Column * object = new RL_Column(parser, xm, ym, d, counterID);

			map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<RL_Agent *>();

			shapesMap[name].push_back(object);

			counterID++;
		}*/
		/*
		else if( name == "RL_Columns" )
		{
			name = "RL_Column";
			float variable = 0.0;
			int variableInt = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," n=%d",&variableInt);
			const int number = variableInt;
			assert(number>=0);

			for(int j=0; j<number; j++)
			{
				const Real radius = rng.uniform(0.0,0.45);
				const Real angle = rng.uniform(0.0,2*M_PI);
				const Real xx = radius*cos(angle);
				const Real yy = radius*sin(angle);
				//RL_Column * object = new RL_Column(parser, xx+0.5, yy+0.5, d, counterID);
				RL_Column * object = new RL_Column(parser, rng.uniform(0,1), rng.uniform(0,1), d, counterID);

				map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

				if(it==shapesMap.end())
					shapesMap[name] = vector<RL_Agent *>();

				shapesMap[name].push_back(object);

				counterID++;
			}
		}
		*/
		/*
		else if( name == "RL_DynamicColumn" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;

			RL_DynamicColumn * object = new RL_DynamicColumn(parser, xm, ym, d, counterID);

			map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<RL_Agent *>();

			shapesMap[name].push_back(object);

			counterID++;
		}
		else if( name == "RL_DynamicColumns" )
		{
			name = "RL_DynamicColumn";
			float variable = 0.0;
			int variableInt = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," n=%d",&variableInt);
			const int number = variableInt;
			assert(number>=0);

			for(int j=0; j<number; j++)
			{
				const Real radius = rng.uniform(0.0,0.45);
				const Real angle = rng.uniform(0.0,2*M_PI);
				const Real xx = radius*cos(angle);
				const Real yy = radius*sin(angle);
				RL_DynamicColumn * object = new RL_DynamicColumn(parser, rng.uniform(0,1), rng.uniform(0,1), d, counterID);

				map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

				if(it==shapesMap.end())
					shapesMap[name] = vector<RL_Agent *>();

				shapesMap[name].push_back(object);

				counterID++;
			}
		}*/


		///else if( name == "RL_Food" )
		if( name == "RL_Food" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			const Real radius = rng.uniform(0.0,0.45);
			const Real angle = rng.uniform(0.0,2*M_PI);
			const Real xx = radius*cos(angle);
			const Real yy = radius*sin(angle);


			RL_Food * object = new RL_Food(parser, 0.5, 0.5, d, counterID);

			map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<RL_Agent *>();

			shapesMap[name].push_back(object);

			counterID++;
		}
		else if( name == "RL_Foods" )
		{
			name = "RL_Food";
			float variable = 0.0;
			int variableInt = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," n=%d",&variableInt);
			const int number = variableInt;
			assert(number>=0);

			for(int j=0; j<number; j++)
			{
				const Real radius = rng.uniform(0.0,0.3);
				const Real angle = rng.uniform(0.0,2*M_PI);
				const Real xx = radius*cos(angle);
				const Real yy = radius*sin(angle);

				RL_Food * object = new RL_Food(parser, xx+0.5, yy+0.5, d, counterID, rand());

				map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

				if(it==shapesMap.end())
					shapesMap[name] = vector<RL_Agent *>();

				shapesMap[name].push_back(object);

				counterID++;
			}
		}
		else if( name == "RL_SmartyGlutton" )
		{
			float variable = 0.0;
			int variableInt = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;
			assert(T>=0.01);
			fscanf(ppFile," nactions=%d",&variableInt);
			const int nactions = variableInt;
			fscanf(ppFile," cutoffGluttons=%f",&variable);
			const Real cutoffGluttons = variable;
			fscanf(ppFile," cutoffFood=%f",&variable);
			const Real cutoffFood = variable;
			fscanf(ppFile," steer=%f",&variable);
			const Real steer = variable;
			fscanf(ppFile," noise=%f",&variable);
			const Real noise = variable;
			fscanf(ppFile," selfavoid=%f",&variable);
			const Real selfavoid = variable;
			fscanf(ppFile," NN=%d",&variableInt);
			const int NN = variableInt;

			RL_SmartyGlutton * object = new RL_SmartyGlutton(parser, xm, ym, d, T, nactions, cutoffGluttons, cutoffFood, steer, noise, selfavoid, NN, counterID, policy, rand());

			map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<RL_Agent *>();

			shapesMap[name].push_back(object);

			counterID++;
		}
		else if( name == "RL_SmartyGluttons" )
		{
			name = "RL_SmartyGlutton";
			float variable = 0.0;
			int variableInt = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			assert(d>0.0);
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;
			assert(T>=0.01);
			fscanf(ppFile," nactions=%d",&variableInt);
			const int nactions = variableInt;
			fscanf(ppFile," cutoffGluttons=%f",&variable);
			const Real cutoffGluttons = variable;
			fscanf(ppFile," cutoffFood=%f",&variable);
			const Real cutoffFood = variable;
			fscanf(ppFile," steer=%f",&variable);
			const Real steer = variable;
			fscanf(ppFile," noise=%f",&variable);
			const Real noise = variable;
			fscanf(ppFile," selfavoid=%f",&variable);
			const Real selfavoid = variable;
			fscanf(ppFile," NN=%d",&variableInt);
			const int NN = variableInt;
			fscanf(ppFile," n=%d",&variableInt);
			const int number = variableInt;
			assert(number>=0);

			const int sidemore = (int)ceil(sqrt((Real)number));
			const Real spacing = 1.0/(Real)sidemore;

			int countside = 0;
			for(int ix=0; ix<sidemore; ix++)
				for(int iy=0; iy<sidemore; iy++)
				{
					if(countside<number)
					{
						//const Real radius = rng.uniform(0.0,0.25);
						//const Real angle = rng.uniform(0.0,2*M_PI);
						//const Real xx = 0.5+radius*cos(angle);
						//const Real yy = 0.5+radius*sin(angle);

						const Real xx = ix*spacing + spacing/2.0 + rng.uniform(0.0,spacing/2.0);
						const Real yy = iy*spacing + spacing/2.0 + rng.uniform(0.0,spacing/2.0);

						RL_SmartyGlutton * object = new RL_SmartyGlutton(parser, xx, yy, d, T, nactions, cutoffGluttons, cutoffFood, steer, noise, selfavoid, NN, counterID, policy, rand());

						map< string, vector<RL_Agent *> >::iterator it = shapesMap.find(name);

						if(it==shapesMap.end())
							shapesMap[name] = vector<RL_Agent *>();

						shapesMap[name].push_back(object);

						counterID++;
						countside++;
					}
				}
		}
		else
		{
			printf("Case not defined!! One of your shapes' name in the factory file is wrong!\n");
			abort();
		}

	}

	fclose(ppFile);

	printf("--------------RL OBJECT FACTORY: END--------------\n");
}

} /* namespace RL */
