/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include "I2D_MRAGOptimisation.h"
#include "I2D_FlowPastFixedObstacle.h"
#include "I2D_FluidMediatedInteractions.h"
#include "I2D_FlowPastFloatingObstacle.h"
#include "I2D_RotatingCylinderPair.h"
#include "I2D_FTLE.h"
#include "I2D_SmartInteractions.h"

using namespace MRAG;
using namespace std;

I2D_Test * test = NULL;
#undef _MRAG_GLUT_VIZ

#ifdef _MRAG_GLUT_VIZ 
struct VisualSupport
{	
	static void display()
	{
	}

	static void idle(void)
	{
		glClear(GL_COLOR_BUFFER_BIT);
		test->run();
		glutSwapBuffers();
	}

	static void run(int argc, const char ** argv)
	{
		static bool bSetup = false;

		if (!bSetup)
		{
			setup(argc, argv);
			bSetup = true;
		}

		glutDisplayFunc(display);
		glutIdleFunc(idle);

		glutMainLoop();
	}

	static void setup(int argc,  const char ** argv)
	{
		glutInit(&argc, const_cast<char **>(argv));
		glutInitWindowSize(800,800);
		glutInitWindowPosition(0, 0);
		glutInitDisplayMode(GLUT_DEPTH| GLUT_STENCIL |GLUT_RGBA | GLUT_DOUBLE );

		glutCreateWindow("Fish");

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		glOrtho(0.0, 1.0, 0.0, 1.0, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glEnable(GL_TEXTURE_2D);

		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	}
};
#endif

int main (int argc, const char ** argv) 
{
	printf("///////////////////////////////////////////////////////////////\n");
	printf("////////////        THIS IS CHLOE'S MAIN        ///////////////\n");
	printf("///////////////////////////////////////////////////////////////\n");

	ArgumentParser parser(argc, argv);

	Environment::setup(max(1, parser("-nthreads").asInt()));

	if( parser("-study").asString() == "FLOW_PAST_FLOATING_OBSTACLE" )
		test = new I2D_FlowPastFloatingObstacle(argc, argv);
	else if( parser("-study").asString() == "OPTIMIZATION_FLOW_PAST_FIXED_OBSTACLE" )
		test = new I2D_MRAGOptimisation(argc, argv);
	else if( parser("-study").asString() == "FLOW_PAST_FIXED_OBSTACLE" )
		test = new I2D_FlowPastFixedObstacle(argc, argv);
	else if( parser("-study").asString() == "FLUID_MEDIATED_INTERACTIONS" )
		test = new I2D_FluidMediatedInteractions(argc, argv);
	else if( parser("-study").asString() == "SMART_INTERACTIONS" )
		test = new I2D_SmartInteractions(argc, argv);
	else if( parser("-study").asString() == "FTLE" )
		test = new I2D_FTLE(argc, argv);
	else if( parser("-study").asString() == "ROTATING_CYLINDER_PAIR" )
        test = new I2D_RotatingCylinderPair(argc, argv);
	else
	{
		printf("Study case not defined!\n"); 
		abort();
	}

	tbb::tick_count t1,t0;

	{
		t0=tbb::tick_count::now();

#ifdef _MRAG_GLUT_VIZ 
		VisualSupport::run(argc, argv);
#else
		test->run();
#endif

		t1=tbb::tick_count::now();
	}

	printf("we spent: %2.2f \n",(t1-t0).seconds());

	return 0;
}
