/*
 *  I2D_main.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */


#include "I2D_TestDumping.h"
#include "I2D_TestDiffusion.h"
#include "I2D_TestAdvection.h"
#include "I2D_TestPenalizationAndOther.h"
#include "I2D_TestPoissonEquation.h"
#include "I2D_TestPoissonEquationPotential.h"
#include "I2D_TestMultipole.h"

using namespace MRAG;
using namespace std;

I2D_Test * test = NULL;

#ifdef _MRAG_GLUT_VIZ 
struct VisualSupport
{	
	static void display()
	{
		glClear(GL_COLOR_BUFFER_BIT);
		
		test->run();
		test->paint();
		
		glutSwapBuffers();
	}
	
	static void idle(void) { glutPostRedisplay(); }
	
	static void run(int argc,  char ** argv)
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
	
	static void setup(int argc,  char ** argv)
	{
		glutInit(&argc, argv);
		glutInitWindowSize(800,800);
		glutInitWindowPosition(0, 0);
		glutInitDisplayMode(GLUT_DEPTH| GLUT_STENCIL |GLUT_RGBA | GLUT_DOUBLE );
		
		glutCreateWindow("HYPRE Poisson Solvers");
		
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		
		glOrtho(-0.2, 1.2, -0.2, 1.2, -1, 1);
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

int main (int argc,  const char ** argv) 
{
	ArgumentParser parser(argc, argv);
	
	Environment::setup(max(1, parser("-nthreads").asInt()));
	
	printf("INPUT IS %s\n", parser("-study").asString().data());
	if(parser("-study").asString() == "diffusion")
		test = new I2D_TestDiffusion(argc, (const char **)argv);
	else if(parser("-study").asString() == "advection")
		test = new I2D_TestAdvection(argc, (const char **)argv);
	else if(parser("-study").asString() == "penalization")
		test = new I2D_TestPenalizationAndOther(argc, (const char **)argv);
	else if(parser("-study").asString() == "poissoneq")
		test = new I2D_TestPoissonEquation(argc, (const char **)argv);
	else if(parser("-study").asString() == "poissoneqpot")
			test = new I2D_TestPoissonEquationPotential(argc, (const char **)argv);
	else if(parser("-study").asString() == "multipole")
        test = new I2D_TestMultipole(argc, (const char **)argv);
	else if(parser("-study").asString() == "dumping")
		test = new I2D_TestDumping(argc, (const char **)argv);
	else
	{
		printf("Study case is not set!\n");
		abort();
	}
	
	tbb::tick_count t1,t0;
	
	{
		t0=tbb::tick_count::now();
		
#ifdef _MRAG_GLUT_VIZ 
		VisualSupport::run(argc, argv);
#else
		const int nsteps = max(50, parser("-nsteps").asInt());
		for(int i=0; i<nsteps; i++) 
			test->run();
#endif
		
		t1=tbb::tick_count::now();
	}
	
	printf("we spent: %2.2f \n",(t1-t0).seconds());
	
    return 0;
}
