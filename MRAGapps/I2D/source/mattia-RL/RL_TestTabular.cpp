#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <iostream>
#include <vector>

#include "rng.h"
#include "RL_TestTabular.h"
#include "RL_QLearning.h"
#include "RL_AgentVector.h"
#include "RL_ObjectFactory.h"

#ifdef _RL_VIZ
#include "FotoCamera.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include "GLUT/glut.h"
#else
#include <GL/gl.h>
#endif
#endif

namespace RL
{

RL_TestTabular::RL_TestTabular(const int argc, const char ** argv):	parser(argc, argv), policy(NULL)
{
	printf("//////////////////////////////////////////////////////////////////////\n");
	printf("////////////            REINFORCEMENT LEARNING         ///////////////\n");
	printf("//////////////////////////////////////////////////////////////////////\n");

	parser.set_strict_mode();
	SAVEFREQ = parser("-savefreq").asInt();
	XPOS = parser("-xpos").asDouble();
	YPOS = parser("-ypos").asDouble();
	CHARLENGTH = parser("-D").asDouble();
	parser.save_options();

	assert( SAVEFREQ >= 0.0 );
	assert( XPOS>=0.0 && XPOS<=1.0 );
	assert( YPOS>=0.0 && YPOS<=1.0 );
	assert( CHARLENGTH>=0.0 && CHARLENGTH<=1.0 );

	_prepareAgents();
}

void RL_TestTabular::_refresh()
{
	_dispose();
	_prepareAgents();
}

void RL_TestTabular::_prepareAgents()
{
	assert(agent==NULL);

	map< string, vector<RL_Agent *> > shapesMap;
	RL_ObjectFactory factory(CHARLENGTH, XPOS, YPOS);
	factory.create(parser,shapesMap,&policy);
	agent = new RL_AgentVector(parser,shapesMap);

	assert(agent!=NULL);

	agent->restartPolicy();
}

void RL_TestTabular::_dispose()
{
	if(agent!=NULL)
	{
		delete agent;
		agent = NULL;
	}

	assert(agent==NULL);

	if(policy!=NULL)
	{
		delete policy;
		policy = NULL;
	}

	assert(policy==NULL);
}

void RL_TestTabular::_save()
{
	agent->savePolicy();
	printf("Policies saved successfully!\n");
}

RL_TestTabular::~RL_TestTabular()
{
	_dispose();
}

#ifdef _RL_VIZ
void RL_TestTabular::_paint()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushAttrib(GL_ENABLE_BIT);
	agent->paint();
	glPopAttrib();

	glutSwapBuffers();
}
#endif

void RL_TestTabular::run()
{
/*
	const Real myX[2] = {0.5,0.5};
	const Real target[2] = {0.55,0.55};

	const Real myV[2] = {0.0,1.0};

	Real d[2] = {0.0,.0};

	//_distPeriodic(target, myX, d);
	const Real x2[2] = {target[0],target[1]};
	const Real x1[2] = {myX[0],myX[1]};
	const Real domain_size = 1.0;
	const Real inv_domain_size = 1.0/domain_size;
	const Real _rx = x2[0] - x1[0];
	const Real _ry = x2[1] - x1[1];
	d[0] = _rx - domain_size*floor(0.5+_rx*inv_domain_size);
	d[1] = _ry - domain_size*floor(0.5+_ry*inv_domain_size);
	cout << "d[0]=" << d[0] << "  d[1]=" << d[1] << std::endl;


	//_angleVectors(d, myV);
	//const Real v1[2] = {d[0],d[1]};
	//const Real v2[2] = {myV[0],myV[1]};
	//const Real anglev = atan2(v1[1],v1[0]) /M_PI*180.0;
	//const Real angled = atan2(v2[1],v2[0]) /M_PI*180.0;
	//const Real angle = anglev-angled;
	//const Real value = (angle<0.0)?angle+360.0:angle;
	//cout << "angle=" << value << std::endl;

	//right
	const Real v1[2] = {-1,0};
	Real v2[2] = {0.0,0.0};
	const Real alpha = 5/180.0*M_PI;
	//void RL_Agent::_rotate(const Real alpha, const Real v1[2], Real v2[2]) const
	v2[0] = v1[0]*cos(alpha) - v1[1]*sin(alpha);
	v2[1] = v1[0]*sin(alpha) + v1[1]*cos(alpha);
	cout << "v2[0]=" << v2[0] << "  v2[1]=" << v2[1] << std::endl;

	exit(0);
	*/

#ifdef _RL_VIZ
	const unsigned int framePerCycle = 2;
	const Real cycle = 1.0;
	const Real fotoDT = cycle/((Real)framePerCycle-1.0);
	Real fotoTimer = 0.0;
	FotoCamera foto;
#endif

	// Retrieve number of gluttons
	RL_AgentVector * b = static_cast<RL_AgentVector*>(agent);
	map< string, vector<RL_Agent *> >::iterator itmap = b->data.find("RL_SmartyGlutton");
	const long unsigned int ngluttons = (itmap!=b->data.end())?itmap->second.size():0;
	const long unsigned int totNumPolicyEval = 1e15;
	const long unsigned int nplays = (ngluttons==0)?totNumPolicyEval:((int)ceil((Real)totNumPolicyEval/(Real)ngluttons));

	// Start simulation
	double t = 0.0;
	const double dt = 0.005;
	long unsigned int step_id = 0;
	while(true)
	{
		// Choose an action
		profiler.push_start("CHOOSE");
		const bool valid = agent->choose(t);
		profiler.pop_stop();
		if(!valid){ printf("NON VALID: REFRESH!\n"); agent->savePolicy(); _refresh(); continue; }

		// Update agents
		profiler.push_start("UPDATE");
		agent->update(dt,t);
		profiler.pop_stop();

		// Compute rewards
		profiler.push_start("REWARD");
		agent->reward(t);
		profiler.pop_stop();

		profiler.push_start("LEARN");
		agent->learn(t);
		profiler.pop_stop();

		if(step_id % SAVEFREQ == 0)
		{
			profiler.printSummary();
			_save();
		}

		t += dt;
		step_id++;

		//if(t>300)
		//	exit(0);

		if(step_id==nplays)
		{
			printf("Simulation succesfully finished!\n");
			exit(0);
		}

#ifdef _RL_VIZ
		profiler.push_start("PAINT");
		fotoTimer += dt;
		if(fotoTimer > fotoDT || step_id==0)
		{
			fotoTimer = 0.0;
			_paint();
			//foto.shoot();
		}
		profiler.pop_stop();
#endif
	}
}

}
