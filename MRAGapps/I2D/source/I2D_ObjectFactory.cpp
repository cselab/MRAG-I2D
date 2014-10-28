/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#include <iostream>
#include <fstream>
#include "RL_QLearning.h"
#include "I2D_FloatingObstacleOperator.h"
#include "I2D_ObjectFactory.h"
#include "I2D_FloatingCylinder.h"
#include "I2D_FloatingEllipse.h"
#include "I2D_ImposedCylinder.h"
#include "I2D_ImposedEllipse.h"
#include "I2D_ImposedEllipseRotation.h"
#include "I2D_ImposedEllipseOscillation.h"
#include "I2D_ImposedCylinderRotation.h"
#include "I2D_ImposedCylinderOscillation.h"
#include "I2D_ImposedCylinderStopAndGo.h"
#include "I2D_CarlingFish.h"
#include "I2D_StefanFish.h"
#include "I2D_CStartLarva.h"
#include "I2D_CarlingFishMorph.h"
#include "I2D_StefanFishMorph.h"
#include "I2D_StefanFishSmart.h"
#include "I2D_CarlingFishAirfoil.h"
#include "I2D_StefanFishSmartInline.h"
#include "I2D_StefanFishSmartLattice.h"
#include "I2D_PassiveTracer.h"
#include "I2D_StefanFishSmartSpiral.h"
#include "I2D_ThreeLinkFish.h"
#include "I2D_ImposedAirfoil.h"
#include "I2D_PitchingAirfoil.h"
#include "I2D_FloatingRotatingCylinderPair.h"

I2D_ObjectFactory::I2D_ObjectFactory(Grid<W,B>& grid, const Real charLength, const Real XCM, const Real YCM, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int _LMAX):
grid(grid), eps(eps), charLength(charLength), XCM(XCM), YCM(YCM), penalization(penalization), LMAX(_LMAX)
{
	this->Uinf[0] = Uinf[0];
	this->Uinf[1] = Uinf[1];
}


I2D_ObjectFactory::~I2D_ObjectFactory()
{
}

int I2D_ObjectFactory::_lines(string filename)
{	
	// Open file stream 
	ifstream filestream(filename.c_str());

	// Check that the filestream is correctly open
	if(!filestream.good())
	{
		cout << "ooops: file not found. Exiting now." << endl;
		exit(-1);
	}

	// Count number of lines contained in the file
	int c = 0;
	string line;
	while( getline(filestream, line) ) c++;
	filestream.close();	

	// Return number of lines
	return c;
}

void I2D_ObjectFactory::create(ArgumentParser & parser, map< string, vector<I2D_FloatingObstacleOperator *> > & shapesMap, bool smart, RL::RL_TabularPolicy ** policy)
{
	shapesMap.clear();
	modulusMaxV = sqrt(Uinf[0]*Uinf[0]+Uinf[1]*Uinf[1]);

	// Read parser information
	const Real LR = parser("-lr").asDouble();
	const Real GAMMA = parser("-gamma").asDouble();
	const Real GREEDYEPS = parser("-greedyEps").asDouble();
	const bool SHARED = parser("-shared").asBool();
	const string factoryFile = parser("-factory").asString();

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

	if(smart)
		if(SHARED)
			(*policy) = new RL::RL_QLearning(LR,GAMMA,GREEDYEPS,rand());

	// Count number of shapes
	const int N = _lines(factoryFile);

	// Open factory file and retrieve information
	FILE * ppFile;
	ppFile = fopen(factoryFile.c_str(),"r");
	if(ppFile==NULL){ printf("could not open ctrl file %s!\n", factoryFile.c_str()); abort(); }

	printf("\n--------------OBJECT FACTORY: START--------------\n");	
	printf("number of objects = %d\n",N);

	for(int i=0; i<N; i++)
	{
		char nameObject[1000];
		fscanf(ppFile,"%s",nameObject);
		string name(nameObject);

		if( name == "I2D_FloatingCylinder" )
		{
			float variable = 0.0;			
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;

			I2D_FloatingCylinder * object = new I2D_FloatingCylinder(parser, grid, xm, ym, d, eps, Uinf, penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);
		}
		else if( name == "I2D_PassiveTracer" )
		{
			float variable = 0.0;
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;

			I2D_PassiveTracer * object = new I2D_PassiveTracer(parser, grid, xm, ym, 0.0, eps, Uinf, penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);
		}
		else if( name == "I2D_FloatingEllipse" )
		{
			float variable = 0.0;			
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," aspectRatio=%f",&variable);
			const Real aspectRatio = variable;

			I2D_FloatingEllipse * object = new I2D_FloatingEllipse(parser, grid, xm, ym, d, aspectRatio, eps, Uinf, penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);
		}
		else if( name == "I2D_ImposedCylinder" )
		{
			float variable = 0.0;			
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," vx=%f",&variable);
			const Real vx = variable*charLength;
			fscanf(ppFile," vy=%f",&variable);
			const Real vy = variable*charLength;

			I2D_ImposedCylinder * object = new I2D_ImposedCylinder(parser, grid, xm, ym, vx, vy, d, eps, Uinf, penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,(Real)sqrt(vx*vx+vy*vy));
		}
		else if(name == "I2D_ImposedCylinder_DiamArray") // diamond array
		{
			// read values form factory file
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable); // physical x-coord of array center
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable); // physical y-coord of array center
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," N=%f",&variable); // number of objects to create, must be square
			const Real Na = variable;
			if (sqrt(Na) - int(sqrt(Na)) != 0){ printf("Wrong number of objects\n"); abort(); }
			fscanf(ppFile," b=%f",&variable); // horizontal spacing between objects
			const Real b = variable*charLength;
			fscanf(ppFile," h=%f",&variable); // vertical spacing between objects
			const Real h = variable*charLength;

			// ensure correct array dimensions
			int s = int (sqrt(Na)-1); // spacing number array center - outmost cylinder, characteristic number
			if (xm-s*b<d/2.0 || xm+s*b>1.0-d/2.0){ printf("Wrong array dimension\n"); abort(); }
			if (ym-s*h<d/2.0 || ym+s*b>1.0-d/2.0){ printf("Wrong array dimension\n"); abort(); }

			// create array
			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);
			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();
			for(int cx=-s; cx<=s; cx++)
				for(int cy=-(s-abs(cx)); cy<=s-abs(cx); cy+=2)
				{
					I2D_ImposedCylinder * object = new I2D_ImposedCylinder(parser, grid, xm+cx*b, ym+cy*h, 0, 0, d, eps, Uinf, penalization);
					shapesMap[name].push_back(object);
				}
		}
		else if(name == "I2D_ImposedCylinder_RectArray") // rectangle array
		{
			// read values from factory file
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable); // physical x-coord of array center
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable); // physical y-coord of array center
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," N=%f",&variable); // number of objects to create, must be square
			const Real Na = variable;
			if (sqrt(Na) - int(sqrt(Na)) != 0){ printf("Wrong number of objects\n"); abort(); }
			fscanf(ppFile," b=%f",&variable); // horizontal spacing between objects
			const Real b = variable*charLength;
			fscanf(ppFile," h=%f",&variable); // vertical spacing between objects
			const Real h = variable*charLength;

			// ensure correct array dimensions
			float s = (sqrt(Na)-1.0)/2.0; // spacing number array center - outmost cylinder, characteristic number
			if (xm-s*b<d/2.0 || xm+s*b>1.0-d/2.0){ printf("Wrong array dimension\n"); abort(); }
			if (ym-s*h<d/2.0 || ym+s*b>1.0-d/2.0){ printf("Wrong array dimension\n"); abort(); }

			// create array
			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);
			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();
			for(float cx=-s; cx<=s; cx++)
				for(float cy=-s; cy<=s; cy++)
				{
					I2D_ImposedCylinder * object = new I2D_ImposedCylinder(parser, grid, xm+cx*b, ym+cy*h, 0, 0, d, eps, Uinf, penalization);
					shapesMap[name].push_back(object);
				}
		}
		else if( name == "I2D_ImposedAirfoil" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," d1=%f",&variable); // maximum camber as percentage of the chord
			const int d1 = int (variable);
			if(d1<0 || d1>9){ printf("Wrong maximum camber\n"); abort(); }
			fscanf(ppFile," d2=%f",&variable); // distance of maximum camber from leading edge in tens of percent's from the chord
			const int d2 = int (variable);
			if(d2<0 || d2>9){ printf("Wrong distance of maximum camber\n"); abort(); }
			fscanf(ppFile," d3d4=%f",&variable); // maximum thickness of the airfoil as percentage of the chord
			const int d3d4 = int (variable);
			if(d3d4<0 || d3d4>99){ printf("Wrong maximum thickness\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable); // physical x-coordinate leading edge
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable); // physical y-coordinate leading edge
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," angle=%f",&variable); // orientation of NACA profile, [ °]
			const Real angle = variable;
			fscanf(ppFile," vx=%f",&variable);
			const Real vx = variable*charLength;
			fscanf(ppFile," vy=%f",&variable);
			const Real vy = variable*charLength;

			I2D_ImposedAirfoil * object = new I2D_ImposedAirfoil(grid, parser, xm, ym, d, angle, vx, vy, d1, d2, d3d4, eps, Uinf, penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,(Real)sqrt(vx*vx+vy*vy));
		}
		else if( name == "I2D_PitchingAirfoil" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," d1=%f",&variable); // maximum camber as percentage of the chord
			const int d1 = int (variable);
			if(d1<0 || d1>9){ printf("Wrong maximum camber\n"); abort(); }
			fscanf(ppFile," d2=%f",&variable); // distance of maximum camber from leading edge in tens of percent's from the chord
			const int d2 = int (variable);
			if(d2<0 || d2>9){ printf("Wrong distance of maximum camber\n"); abort(); }
			fscanf(ppFile," d3d4=%f",&variable); // maximum thickness of the airfoil as percentage of the chord
			const int d3d4 = int (variable);
			if(d3d4<0 || d3d4>99){ printf("Wrong maximum thickness\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable); // physical x-coordinate leading edge
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable); // physical y-coordinate leading edge
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," l=%f",&variable); // distance in x between leading edge airfoil and center of rotation
			const Real l = variable*charLength;
			if(l<0.0 || l>d){ printf("Wrong center of rotation\n"); abort(); }
			fscanf(ppFile," plAmp=%f",&variable); // plunge Amplitude
			const Real plAmp = variable*charLength;
			fscanf(ppFile," plPer=%f",&variable); // plunge Period
			const Real plPer = variable;
			fscanf(ppFile," plPhase=%f",&variable); // plunge Phase
			const Real plPhase = variable;
			fscanf(ppFile," piAmp=%f",&variable); // pitch Amplitude, [rad], rotation: [-piAmp, piAmp]
			const Real piAmp = variable;
			fscanf(ppFile," piPer=%f",&variable); // pitch Period
			const Real piPer = variable;
			fscanf(ppFile," piPhase=%f",&variable); // pitch Phase
			const Real piPhase = variable;

			I2D_PitchingAirfoil * object = new I2D_PitchingAirfoil(grid, parser, xm, ym, d,
					l, d1, d2, d3d4, plAmp, plPer, plPhase, piAmp, piPer, piPhase, eps, Uinf,
					penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			//assert(object->getModulusMaxVel() > 0.0); // remove comment later
		}
		else if(name == "I2D_ImposedAirfoil_DiamArray") // diamond array
		{
			// read values from factory file
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," d1=%f",&variable); // maximum camber as percentage of the chord
			const int d1 = int (variable);
			if(d1<0 || d1>9){ printf("Wrong maximum camber\n"); abort(); }
			fscanf(ppFile," d2=%f",&variable); // distance of maximum camber from leading edge in tens of percent's from the chord
			const int d2 = int (variable);
			if(d2<0 || d2>9){ printf("Wrong distance of maximum camber\n"); abort(); }
			fscanf(ppFile," d3d4=%f",&variable); // maximum thickness of the airfoil as percentage of the chord
			const int d3d4 = int (variable);
			if(d3d4<0 || d3d4>99){ printf("Wrong maximum thickness\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable); // physical x-coord of array center
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable); // physical y-coord of array center
			const Real ym = YCM + variable*charLength;
			fscanf (ppFile," angle=%f",&variable); // orientation of NACA profiles
			const Real angle = variable;
			fscanf(ppFile," N=%f",&variable); // number of objects to create, must be square
			const Real Na = variable;
			if (sqrt(Na) - int(sqrt(Na)) != 0){ printf("Wrong number of objects\n"); abort(); }
			fscanf(ppFile," b=%f",&variable); // horizontal spacing between objects
			const Real b = variable*charLength;
			if (b < d/2.0){ printf("Wrong horizontal spacing\n"); abort(); }
			fscanf(ppFile," h=%f",&variable); // vertical spacing between objects
			const Real h = variable*charLength;
			if (h < d3d4/100.0*d){ printf("Wrong vertical spacing\n"); abort(); }

			// ensure correct array dimensions
			int s = int (sqrt(Na)-1); // spacing number array center - outmost airfoil, characteristic number
			if (xm-s*b<0.0 || xm+s*b>1.0-d){ printf("Wrong array dimension\n"); abort(); }
			if (ym-s*h<d3d4/100.0*d || ym+s*b>1.0-d3d4/100.0*d){ printf("Wrong array dimension\n"); abort(); }

			// create array
			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);
			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();
			for(int cx=-s; cx<=s; cx++)
				for(int cy=-(s-abs(cx)); cy<=s-abs(cx); cy+=2)
				{
					I2D_ImposedAirfoil * object = new I2D_ImposedAirfoil(grid, parser, xm+cx*b, ym+cy*h, d, angle, 0.0, 0.0, d1, d2, d3d4, eps, Uinf, penalization);
					shapesMap[name].push_back(object);
				}
		}
		else if( name == "I2D_ImposedEllipse" )
		{
			float variable = 0.0;			
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," vx=%f",&variable);
			const Real vx = variable*charLength;
			fscanf(ppFile," vy=%f",&variable);
			const Real vy = variable*charLength;
			fscanf(ppFile," aspectRatio=%f",&variable);
			const Real aspectRatio = variable;

			I2D_ImposedEllipse * object = new I2D_ImposedEllipse(parser, grid, xm, ym, xm, ym, vx, vy, d, aspectRatio, eps, Uinf, penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,(Real)sqrt(vx*vx+vy*vy));
		}
		else if( name == "I2D_ImposedEllipseRotation" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xcrot = XCM + variable*charLength;
			const Real xm = xcrot + d/2.0;
			fscanf(ppFile," ym=%f",&variable);
			const Real ycrot = YCM + variable*charLength;
			const Real ym = ycrot;
			fscanf(ppFile," aspectRatio=%f",&variable);
			const Real aspectRatio = variable;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;

			I2D_ImposedEllipseRotation * object = new I2D_ImposedEllipseRotation(parser, grid, xm, ym, xcrot, ycrot, T, d, aspectRatio, eps, Uinf, penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_ImposedEllipseOscillation" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," aspectRatio=%f",&variable);
			const Real aspectRatio = variable;
			fscanf(ppFile," plAmp=%f",&variable);
			const Real plAmp = variable*charLength;
			fscanf(ppFile," plPer=%f",&variable);
			const Real plPer = variable;
			fscanf(ppFile," plPhase=%f",&variable);
			const Real plPhase = variable;
			fscanf(ppFile," piAmp=%f",&variable);
			const Real piAmp = variable;
			fscanf(ppFile," piPer=%f",&variable);
			const Real piPer = variable;
			fscanf(ppFile," piPhase=%f",&variable);
			const Real piPhase = variable;
			I2D_ImposedEllipseOscillation * object = new I2D_ImposedEllipseOscillation(parser, grid, xm, ym, xm, ym, aspectRatio, plAmp, plPer, piAmp, piPer, plPhase, piPhase, d, eps, Uinf, penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV, object->getModulusMaxVel());
			assert(object->getModulusMaxVel() > 0.0);
		}
		else if( name == "I2D_ImposedCylinderRotation" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," radiusRot=%f",&variable);
			const Real radiusRot = variable*charLength;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;
			fscanf(ppFile," alignment=%f",&variable);
			const Real align = variable;
			I2D_ImposedCylinderRotation * object = new I2D_ImposedCylinderRotation(parser, grid, xm, ym, radiusRot, T, d, align, eps, Uinf, penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_ImposedCylinderOscillation" )
		{
			float variable = 0.0;			
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," A=%f",&variable);
			const Real A = variable*charLength;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;
			assert(T>0.0);
			char dir[1000];
			fscanf(ppFile," dir=%s",dir);
			string direction(dir);
			I2D_ImposedCylinderOscillation * object = new I2D_ImposedCylinderOscillation(parser, grid, xm, ym, A, T, direction, d, eps, Uinf, penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_ImposedCylinderStopAndGo" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," vx=%f",&variable);
			const Real vx = variable*charLength;
			fscanf(ppFile," vy=%f",&variable);
			const Real vy = variable*charLength;
			fscanf(ppFile," tstop=%f",&variable);
			const Real tstop = variable;
			assert(tstop>0.0);
			I2D_ImposedCylinderStopAndGo * object = new I2D_ImposedCylinderStopAndGo(parser, grid, xm, ym, vx, vy, tstop, d, eps, Uinf, penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_ThreeLinkFish" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," aspectRatio=%f",&variable);
			const Real aspectRatio = variable;

			I2D_ThreeLinkFish * object = new I2D_ThreeLinkFish(parser, grid, xm, ym, d, aspectRatio, eps, Uinf, penalization);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_CarlingFish" )
		{
			const Real angleInSpace = 0.0;

			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," phase=%f",&variable);
			const Real phase = variable;
			fscanf(ppFile," angle=%f",&variable);
			const Real angle = variable/180.0*M_PI;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;

			I2D_CarlingFish * object = new I2D_CarlingFish(parser, grid, xm, ym, d, T, phase, angle, angleInSpace, eps, Uinf, penalization, LMAX, i);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
			printf("Factory just made a CarlingFish\n");
		}
		else if( name == "I2D_CarlingFishAirfoil" )
		{
			const Real angleInSpace = 0.0;

			char buf[4];
			fscanf(ppFile," naca=%s",buf);
			string naca(buf);
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," phase=%f",&variable);
			const Real phase = variable;
			fscanf(ppFile," angle=%f",&variable);
			const Real angle = variable/180.0*M_PI;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;

			I2D_CarlingFishAirfoil * object = new I2D_CarlingFishAirfoil(parser, grid, naca, xm, ym, d, T, phase, angle, angleInSpace, eps, Uinf, penalization, LMAX, i);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_StefanFish" )
		{
			const Real angleInSpace = 0.0;
			vector<double> BASELINE;
			vector<double> CURVATURE;

			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," phase=%f",&variable);
			const Real phase = variable;
			fscanf(ppFile," tau=%f",&variable);
			const Real tau = variable;
			fscanf(ppFile," angle=%f",&variable);
			const Real angle = variable/180.0*M_PI;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;
			fscanf(ppFile," b1=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b2=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b3=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b4=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b5=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b6=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," k1=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k2=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k3=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k4=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k5=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k6=%f",&variable);
			CURVATURE.push_back(variable);

			I2D_StefanFish * object = new I2D_StefanFish(parser, grid, xm, ym, d, T, phase, tau, angle, BASELINE, CURVATURE, angleInSpace, eps, Uinf, penalization, LMAX, i);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_CStartLarva" )
		{
			const Real angleInSpace = 0.0;
			vector<double> BASELINE;
			vector<double> CURVATURE;

			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," phase=%f",&variable);
			const Real phase = variable;
			fscanf(ppFile," tau=%f",&variable);
			const Real tau = variable;
			fscanf(ppFile," angle=%f",&variable);
			const Real angle = variable/180.0*M_PI;
			fscanf(ppFile," Tprep=%f",&variable);
			const Real Tprep = variable;
			fscanf(ppFile," Tprop=%f",&variable);
			const Real Tprop = variable;
			fscanf(ppFile," b1=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b2=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b3=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b4=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b5=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b6=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," k1=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k2=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k3=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k4=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k5=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k6=%f",&variable);
			CURVATURE.push_back(variable);

			I2D_CStartLarva * object = new I2D_CStartLarva(parser, grid, xm, ym, d, Tprep, Tprop, phase, tau, angle, BASELINE, CURVATURE, angleInSpace, eps, Uinf, penalization, LMAX, i);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_CarlingFishMorph" )
		{
			const Real angleInSpace = 0.0;
			vector<double> WIDTH;

			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," phase=%f",&variable);
			const Real phase = variable;
			fscanf(ppFile," angle=%f",&variable);
			const Real angle = variable/180.0*M_PI;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;
			fscanf(ppFile," w1=%f",&variable);
			WIDTH.push_back(variable);
			fscanf(ppFile," w2=%f",&variable);
			WIDTH.push_back(variable);
			fscanf(ppFile," w3=%f",&variable);
			WIDTH.push_back(variable);
			fscanf(ppFile," w4=%f",&variable);
			WIDTH.push_back(variable);

			I2D_CarlingFishMorph * object = new I2D_CarlingFishMorph(parser, grid, xm, ym, d, T, phase, angle, WIDTH, angleInSpace, eps, Uinf, penalization, LMAX,i);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_StefanFishMorph" )
		{
			const Real angleInSpace = 0.0;
			vector<double> WIDTH;
			vector<double> BASELINE;
			vector<double> CURVATURE;

			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," phase=%f",&variable);
			const Real phase = variable;
			fscanf(ppFile," tau=%f",&variable);
			const Real tau = variable;
			fscanf(ppFile," angle=%f",&variable);
			const Real angle = variable/180.0*M_PI;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;
			fscanf(ppFile," w1=%f",&variable);
			WIDTH.push_back(variable);
			fscanf(ppFile," w2=%f",&variable);
			WIDTH.push_back(variable);
			fscanf(ppFile," w3=%f",&variable);
			WIDTH.push_back(variable);
			fscanf(ppFile," w4=%f",&variable);
			WIDTH.push_back(variable);
			fscanf(ppFile," b1=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b2=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b3=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b4=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b5=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," b6=%f",&variable);
			BASELINE.push_back(variable);
			fscanf(ppFile," k1=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k2=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k3=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k4=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k5=%f",&variable);
			CURVATURE.push_back(variable);
			fscanf(ppFile," k6=%f",&variable);
			CURVATURE.push_back(variable);

			I2D_StefanFishMorph * object = new I2D_StefanFishMorph(parser, grid, xm, ym, d, T, phase, tau, angle, WIDTH, BASELINE, CURVATURE, angleInSpace, eps, Uinf, penalization, LMAX, i);

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_StefanFishSmart" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," tau=%f",&variable);
			const Real tau = variable;
			fscanf(ppFile," angle=%f",&variable);
			const Real angle = variable/180.0*M_PI;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;

			I2D_StefanFishSmart * object = new I2D_StefanFishSmart(parser, grid, xm, ym, d, T, tau, angle, eps, Uinf, penalization, LMAX, i, policy, rand());

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_StefanFishSmartInline" )
		{
			Real dir[2] = {0.0,0.0};
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," tau=%f",&variable);
			const Real tau = variable;
			fscanf(ppFile," angle=%f",&variable);
			const Real angle = variable/180.0*M_PI;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;
			fscanf(ppFile," dirx=%f",&variable);
			dir[0] = variable;
			fscanf(ppFile," diry=%f",&variable);
			dir[1] = variable;
			I2D_StefanFishSmartInline * object = new I2D_StefanFishSmartInline(parser, grid, xm, ym, d, T, tau, angle, dir, eps, Uinf, penalization, LMAX, i, policy, rand());

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_StefanFishSmartLattice" )
		{
			Real dir[2] = {0.0,0.0};
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," tau=%f",&variable);
			const Real tau = variable;
			fscanf(ppFile," angle=%f",&variable);
			const Real angle = variable/180.0*M_PI;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;
			fscanf(ppFile," dirx=%f",&variable);
			dir[0] = variable;
			fscanf(ppFile," diry=%f",&variable);
			dir[1] = variable;
			I2D_StefanFishSmartLattice * object = new I2D_StefanFishSmartLattice(parser, grid, xm, ym, d, T, tau, angle, dir, eps, Uinf, penalization, LMAX, i, policy, rand());

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_StefanFishSmartSpiral" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," tau=%f",&variable);
			const Real tau = variable;
			fscanf(ppFile," angle=%f",&variable);
			const Real angle = variable/180.0*M_PI;
			fscanf(ppFile," T=%f",&variable);
			const Real T = variable;
			I2D_StefanFishSmartSpiral * object = new I2D_StefanFishSmartSpiral(parser, grid, xm, ym, d, T, tau, angle, eps, Uinf, penalization, LMAX, i, policy, rand());

			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);

			if(it==shapesMap.end())
				shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();

			shapesMap[name].push_back(object);

			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}
		else if( name == "I2D_FloatingRotatingCylinderPair" )
		{
			float variable = 0.0;
			fscanf(ppFile," d=%f",&variable);
			const Real d = variable*charLength;
			if(d<=0.0 || d>charLength){ printf("Wrong object dimension\n"); abort(); }
			fscanf(ppFile," xm=%f",&variable);
			const Real xm = XCM + variable*charLength;
			fscanf(ppFile," ym=%f",&variable);
			const Real ym = YCM + variable*charLength;
			fscanf(ppFile," angle=%f",&variable);
			const Real angle = variable;
			fscanf(ppFile," width=%f",&variable);
			const Real width = variable;
			fscanf(ppFile," g1=%f",&variable);
			const Real g1 = variable;
			fscanf(ppFile," g2=%f",&variable);
			const Real g2 = variable;
            
			I2D_FloatingRotatingCylinderPair * object = new I2D_FloatingRotatingCylinderPair(parser, grid, xm, ym, d, angle, width, g1, g2, eps, Uinf, penalization, LMAX, i, policy, rand());
            
			map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = shapesMap.find(name);
            
			if(it==shapesMap.end())
            shapesMap[name] = vector<I2D_FloatingObstacleOperator *>();
            
			shapesMap[name].push_back(object);
            
			modulusMaxV = max(modulusMaxV,object->getModulusMaxVel());
			assert(object->getModulusMaxVel()>0.0);
		}        
		else 
		{
			printf("Case not defined!! One of your shapes' name in the factory file is wrong!\n");
			abort();
		}

	}

	fclose(ppFile);

	printf("--------------OBJECT FACTORY: END--------------\n");

	assert(modulusMaxV != 0.0);
}


