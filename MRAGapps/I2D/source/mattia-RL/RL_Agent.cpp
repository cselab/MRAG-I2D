/*
 * RL_Agent.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: mgazzola
 */

#include <math.h>
#include <assert.h>
#include "RL_Agent.h"


using namespace std;

namespace RL
{

void RL_Agent::_rotate(const Real alpha, const Real v1[2], Real v2[2]) const
{
	v2[0] = v1[0]*cos(alpha) - v1[1]*sin(alpha);
	v2[1] = v1[0]*sin(alpha) + v1[1]*cos(alpha);
}

void RL_Agent::_dist(const Real x2[2], const Real x1[2], Real d[2]) const
{
	d[0] = x2[0] - x1[0];
	d[1] = x2[1] - x1[1];
}

void RL_Agent::_distPeriodic(const Real x2[2], const Real x1[2], Real d[2]) const
{
	const Real domain_size = 1.0;
	const Real inv_domain_size = 1.0/domain_size;

	const Real _rx = x2[0] - x1[0];
	const Real _ry = x2[1] - x1[1];

	d[0] = _rx - domain_size*floor(0.5+_rx*inv_domain_size);
	d[1] = _ry - domain_size*floor(0.5+_ry*inv_domain_size);
}

Real RL_Agent::_modv(const Real v[2]) const
{
	return sqrt( v[0]*v[0] + v[1]*v[1] );
}

Real RL_Agent::_normalize(const Real v[2], Real n[2]) const
{
	const Real IvI = _modv(v);
	n[0] = v[0] / IvI;
	n[1] = v[1] / IvI;
	return IvI;
}

Real RL_Agent::_angleVectors(const Real v1[2], const Real v2[2]) const
{
	const Real anglev = atan2(v1[1],v1[0]) /M_PI*180.0;
	const Real angled = atan2(v2[1],v2[0]) /M_PI*180.0;
	const Real angle = anglev-angled;
	return (angle<0.0)?angle+360.0:angle;
}

int RL_Agent::_discretize(const Real signal, const Real minsignal, const Real maxsignal, const int nlevels) const
{
	assert(maxsignal>=minsignal);
	const Real interval = maxsignal - minsignal;
	const Real delta = interval / (Real)nlevels;
	return max(0,min(nlevels-1,(int)floor((signal - minsignal)/delta)));
}

int RL_Agent::_discretize2SS(const Real signal, const Real minsignal, const Real maxsignal, const int nlevels, const bool deadEnd, const bool inside) const
{
	assert(maxsignal>=minsignal);

	const int numberOfLevels = nlevels-2;

	if (inside)
		return nlevels-1;

	if (deadEnd)
		return nlevels-2;

	const Real interval = maxsignal - minsignal;
	const Real delta = interval / (Real)numberOfLevels;
	return max(0,min(numberOfLevels-1,(int)floor((signal - minsignal)/delta)));
}

int RL_Agent::_discretize1SS(const Real signal, const Real minsignal, const Real maxsignal, const int nlevels, const bool deadEnd) const
{
	assert(maxsignal>=minsignal);

	const int numberOfLevels = nlevels-1;

	if (deadEnd)
		return nlevels-1;

	const Real interval = maxsignal - minsignal;
	const Real delta = interval / (Real)numberOfLevels;
	return max(0,min(numberOfLevels-1,(int)floor((signal - minsignal)/delta)));
}

bool RL_Agent::_isRight(const Real x1[2], const Real dir[2], const Real x2[2]) const
{
	const Real alpha = -M_PI/2.0;
	const Real n[2] = { dir[0]*cos(alpha) - dir[1]*sin(alpha), dir[0]*sin(alpha) + dir[1]*cos(alpha) };
	Real d[2];
	_dist(x2,x1,d);
	return ((d[0]*n[0]+d[1]*n[1])>=0.0);
}

void RL_Agent::savePolicy(string name)
{
	if(policy!=NULL)
		if((*policy)!=NULL)
			(*policy)->save(name);
}

void RL_Agent::restartPolicy(string name)
{
	if(policy!=NULL)
		if((*policy)!=NULL)
			(*policy)->restart(name);
}

bool RL_Agent::startChoose(const double t, map< string, vector<RL_Agent *> > * _data)
{
	bool valid = true;

	if(policy!=NULL && learningTimer==0.0 && status==Ready)
	{
		assert((*policy)!=NULL);
		assert(status == Ready);
		valid = mapState(myState,_data);
	}

	return valid;
}

bool RL_Agent::choose(const double t, map< string, vector<RL_Agent *> > * _data)
{
	bool valid = true;

	if(policy!=NULL && learningTimer==0.0 && status==Ready)
	{
		assert((*policy)!=NULL);
		assert(status == Ready);
		int action = (*policy)->selectAction(myState);
		mapAction(action);
		myState.push_back(action);
		(*policy)->setStateActionStart(myState);
		myStateAction = myState;
		integralReward = 0.0;
		learningTimer = t;
		status = Waiting;
	}

	return valid;
}

void RL_Agent::setReward(const double t)
{
	if(policy!=NULL && (t>(learningTimer+learningInterval)) && status==Waiting )
	{
		assert(status==Waiting);
		assert((*policy)!=NULL);
		(*policy)->setReward(integralReward);
	}
}

void RL_Agent::mapTest(const double t, map< string, vector<RL_Agent *> > * _data)
{
	if(policy!=NULL && (t>(learningTimer+learningInterval)) && status==Waiting )
	{
		assert(status==Waiting);
		assert((*policy)!=NULL);
		mapState(myState,_data);
	}
}

void RL_Agent::stopTest(const double t, map< string, vector<RL_Agent *> > * _data)
{
	if(policy!=NULL && (t>(learningTimer+learningInterval)) && status==Waiting )
	{
		assert(status==Waiting);
		assert((*policy)!=NULL);

		(*policy)->setStateEnd(myState);
		status = Ready;
	}
}

void RL_Agent::learn(const double t, map< string, vector<RL_Agent *> > * _data, string name)
{
	if(policy!=NULL && (t>(learningTimer+learningInterval)) && status==Ready)
	{
		assert((*policy)!=NULL);
		assert(status == Ready);
		(*policy)->update(name);
		learningTimer = 0.0;
	}
}

#ifdef _RL_VIZ
void RL_Agent::_paintSphere(Real x, Real y, Real radius, Real r, Real g, Real b) const
{
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	GLfloat lightpos[] = {0, 0, -1.0, 0};
	glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
	GLfloat lightColor[] = {r,g,b,1};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);

	glPushMatrix();
	glColor3f(r,g,b);
	glTranslated(x,y,0);
	glutSolidSphere(radius, 20,20);
	glPopMatrix();

	/*
	const Real delta = M_PI/180;

	glPushMatrix();

	glColor3f(r,g,b);
	glBegin(GL_POLYGON);
	for (int i=0; i<360; i++)
	{
		const Real degInRad = i*delta;
		glVertex2f(x+cos(degInRad)*radius,y+sin(degInRad)*radius);
	}
	glEnd();

	glColor3f(0,0,0);
	glBegin(GL_LINE_LOOP);
	for (int i=0; i<360; i++)
	{
		const Real degInRad = i*delta;
		glVertex2f(x+cos(degInRad)*radius,y+sin(degInRad)*radius);
	}
	glEnd();

	glPopMatrix();
	*/
}



void RL_Agent::_drawCircle(Real radius, Real xc, Real yc, Real r, Real g, Real b) const
{
	const Real deg2rad = M_PI/90;

	GLfloat lightColorExit[] = {r*1.3,g*1.3,b*1.3,1};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColorExit);

	glPushMatrix();

	glColor3f(r,g,b);
	glBegin(GL_LINE_LOOP);
	for (int i=0; i<360; i++)
	{
		Real degInRad = i*deg2rad;
		glVertex2f(xc+cos(degInRad)*radius,yc+sin(degInRad)*radius);
	}
	glEnd();

	glPopMatrix();
}

void RL_Agent::_drawFullCircle(Real radius, Real xc, Real yc, Real r, Real g, Real b) const
{
	const Real deg2rad = M_PI/90;

	glPushMatrix();

	glColor3f(r,g,b);
	glBegin(GL_POLYGON);
	for (int i=0; i<360; i++)
	{
		Real degInRad = i*deg2rad;
		glVertex2f(xc+cos(degInRad)*radius,yc+sin(degInRad)*radius);
	}
	glEnd();

	glPopMatrix();
}

#endif

}



