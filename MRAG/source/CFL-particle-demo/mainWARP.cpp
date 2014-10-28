/*
 *  mainWARP.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/27/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include "mainWARP.h"

#include <stdlib.h>
#include <iostream>
#include "../MRAGcore/MRAGEnvironment.h"

#ifdef	_MRAG_GLUT_VIZ
#ifndef __apple__
#include "GL/glut.h"
#else
#include "GLUT/glut.h"
#endif
//#include "screenshot.h"
extern "C"	
GLint gltWriteTGA(const char *szFileName);
#endif

const int nMaxSteps = 100000;
int iCurrStep = 0;
DemoWARP * demo = NULL;
#ifdef	_MRAG_GLUT_VIZ

static void display(void)
{
	//glClearColor(1,0,1,0);
	glClear(GL_COLOR_BUFFER_BIT);
	
	demo->simulation_render(true);
	
	
	
	const bool bStop = (iCurrStep++>=nMaxSteps);
	if (bStop) exit(0);
	glFinish();
	
	
	
	static int iFrameCounter = 0;
	char buf[300];
	sprintf(buf, "asdasd%05d.tga", iFrameCounter++);
	printf("Writing to %s\n", buf);
	gltWriteTGA(buf);
	glutSwapBuffers();
	demo->simulation_run(1);
	
}

static void idle(void)
{
	glutPostRedisplay();
}
#endif

int main(int argc, char ** argv)
{
	srand(3290);
	
#ifdef _MRAG_GLUT_VIZ
	{
		glutInit(&argc, argv);
		glutInitWindowSize(800,800);
		glutInitWindowPosition(0, 0);
		glutInitDisplayMode(GLUT_DEPTH| GLUT_STENCIL |GLUT_RGBA | GLUT_DOUBLE );
		
		glutCreateWindow("WARP");
		
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		
		glOrtho(-0.2, 1.2, -0.2, 1.2, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glEnable(GL_TEXTURE_2D);
		
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
		
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		
		glutDisplayFunc(display);
		glutIdleFunc(idle);
	}
#endif	
	Environment::setup();
	demo = new DemoWARP;
	demo->simulation_init();
	
#ifdef _MRAG_GLUT_VIZ			
	glutMainLoop();
#else
	const int nsteps = 4;
	for(int i=0; i<nsteps; i++)
		demo->simulation_run(1);
#endif
}

