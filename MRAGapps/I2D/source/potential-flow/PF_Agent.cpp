/*
 *  PF_Agent.cpp
 *  DipoleCode
 *
 *  Created by Alexia on 3/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#include "assert.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include "RL_QLearning.h"
#include "PF_Agent.h"

using namespace PF;

PF_Agent::PF_Agent(MRAG::ArgumentParser & parser, const Real lr, const Real greedyEps, bool isSmooth, bool isControlled, string _name, const int _ID, RL::RL_TabularPolicy **_policy, const int seed)
		: name(_name), ID(_ID), policy(_policy), policystore(NULL), integralReward(0.0), status(Ready), learningInterval(0.0), learningTimer(0.0), fitnessValue(0.0), errorValue(0.0), SMOOTH(isSmooth), ISCONTROLLED(isControlled), myAction(0), myPreviousAction(0), isAveraged(true)
{
	isDead = false;

	if (policy != NULL)
	{
		if ((*policy) == NULL)
		{
			Real LR = lr;
			Real GREEDYEPS = greedyEps;

			const Real GAMMA = parser("-gamma").asDouble();

			assert((LR >= 0.0 && LR <= 1.0));
			assert(GAMMA >= 0.0 && GAMMA <= 1.0);
			assert((GREEDYEPS >= 0.0 && GREEDYEPS <= 1.0));
			assert((LR+GAMMA) <= 1.0);

			srand(time(NULL));

			if (LR >= 0.0 && LR <= 1.0 && GAMMA >= 0.0 && GAMMA <= 1.0 && GREEDYEPS >= 0.0 && GREEDYEPS <= 1.0)
			{
				policystore = new RL::RL_QLearning(LR, GAMMA, GREEDYEPS, seed, 10);
				policy = &policystore;
				assert(policy!=NULL);
			}
		}
		assert((*policy)!=NULL);
	}
}

PF_Agent::~PF_Agent()
{
	if (policy != NULL)
	{
		if ((*policy) != NULL)
		{
			delete (*policy);
			(*policy) = NULL;
		}
	}
}

bool PF_Agent::mapStart(const Real time, map<string, vector<PF_Agent*> > *collection)
{
	bool valid = true;

	if (policy != NULL && learningTimer == 0.0 && status == Ready)
	{
		assert((*policy)!=NULL);
		assert(status == Ready);
		valid = mapState(myState, collection);
	}

	return valid;
}

bool PF_Agent::choose(const Real time, map<string, vector<PF_Agent*> > *collection)
{
	bool valid = true;

	if (policy != NULL)
	{
		assert((*policy)!=NULL);
		if (learningTimer == 0.0 && status == Ready)
		{
			int action = 0;
			if (ISCONTROLLED) action = (*policy)->selectAction(myState);
			mapAction(action, time);
			myState.push_back(action);
			myPreviousAction = myAction;
			myAction = action;
			(*policy)->setStateActionStart(myState);
			integralReward = 0.0;
			learningTimer = time;
			status = Waiting;
		}

		if (SMOOTH) perform(time);
	}

	return valid;
}

void PF_Agent::setReward(const Real time)
{
	if (policy != NULL && (time > (learningTimer + learningInterval)) && status == Waiting)
	{
		assert((*policy)!=NULL);
		assert(integralReward==integralReward);
		(*policy)->setReward(integralReward);
	}
}

void PF_Agent::mapTest(const Real time, map<string, vector<PF_Agent*> > *collection)
{
	if (policy != NULL && (time > (learningTimer + learningInterval)) && status == Waiting)
	{
		assert((*policy)!=NULL);
		mapState(myState, collection);
	}
}

void PF_Agent::stopTest(const Real time, map<string, vector<PF_Agent*> > *collection)
{
	if (policy != NULL && (time > (learningTimer + learningInterval)) && status == Waiting)
	{
		assert((*policy)!=NULL);
		(*policy)->setStateEnd(myState);
		status = Ready;
	}
}

void PF_Agent::learn(const Real time, map<string, vector<PF_Agent*> > *collection, string filename)
{
	if (policy != NULL && (time > (learningTimer + learningInterval)) && status == Ready)
	{
		assert((*policy)!=NULL);
		(*policy)->update(filename);
		learningTimer = 0.0;
	}
}

void PF_Agent::savePolicy(string name)
{
	if (policy != NULL) if ((*policy) != NULL) (*policy)->save(name);
}

void PF_Agent::restartPolicy(string name)
{
	if (policy != NULL) if ((*policy) != NULL) (*policy)->restart(name);
}

void PF_Agent::_dist(const Real x2[2], const Real x1[2], Real d[2]) const
{
	d[0] = x2[0] - x1[0];
	d[1] = x2[1] - x1[1];
}

Real PF_Agent::_modv(const Real v[2]) const
{
	return sqrt(v[0] * v[0] + v[1] * v[1]);
}

Real PF_Agent::_angleVectors(const Real v1[2], const Real v2[2]) const
{
	const Real anglev = atan2(v1[1], v1[0]) / M_PI * 180.0;
	const Real angled = atan2(v2[1], v2[0]) / M_PI * 180.0;
	const Real angle = anglev - angled;
	return (angle < 0.0) ? angle + 360.0 : angle;
}

void PF_Agent::_normalize(const Real v[2], Real n[2]) const
{
	const Real IvI = _modv(v);
	n[0] = v[0] / IvI;
	n[1] = v[1] / IvI;
}

int PF_Agent::_discretize(const Real signal, const Real minsignal, const Real maxsignal, const int nlevels, const bool deadEnd) const
{
	assert(maxsignal >= minsignal);

	const int numberOfLevels = (deadEnd) ? nlevels - 1 : nlevels;

	if (deadEnd) return numberOfLevels;

	const Real interval = maxsignal - minsignal;
	const Real delta = interval / (Real) numberOfLevels;
	return max(0, min(numberOfLevels - 1, (int) floor((signal - minsignal) / delta)));
}

int PF_Agent::_discretizeRange(const Real value, const Real minvalue, const Real maxvalue, const int levels)
{
	assert(maxvalue>=minvalue);
	const Real h = (maxvalue - minvalue) / levels;
	return max(0, min(levels - 1, (int) floor((value - minvalue) / h)));
}

int PF_Agent::_discretizeAngle(Real angle, const Real sightangle, const int levels)
{
	assert(sightangle>0.0);
	assert(levels % 2 == 0);

	const Real dangle = sightangle / (levels - 1);
	angle += dangle / 2.0;
	angle = (angle > 360.0) ? angle - 360.0 : angle;
	if (angle < 180.0)
	{
		return max(0, min(levels / 2, (int) floor(angle / dangle)));
	}
	else
	{
		return max(levels / 2, min(levels - 1, levels - 1 - (int) floor((360.0 - angle) / dangle)));
	}
}

#ifdef _RL_VIZ
void PF_Agent::_paintSphere(Real x, Real y, Real radius, Real r, Real g, Real b) const
{
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	GLfloat lightpos[] = { 0, 0, -1.0, 0 };
	glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
	GLfloat lightColor[] = { r, g, b, 1 };
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);

	glPushMatrix();
	glColor3f(r, g, b);
	glTranslated(x, y, 0);
	glutSolidSphere(radius, 50, 50);
	glPopMatrix();
}

/**
 * Paint the fuck out of an X on the screen.
 */
void PF_Agent::_paintX(Real x, Real y, Real radius, Real r, Real g, Real b) const
{
	const int nx = 13;
	const Real t = 0.4;
	const Real xx[nx] = { t, 1, 1 - t, 0, t - 1, -1, -t, -1, t - 1, 0, 1 - t, 1, t };
	const Real yx[nx] = { 0, 1 - t, 1, t, 1, 1 - t, 0, t - 1, -1, -t, -1, t - 1, 0 };

	glPushMatrix();

	/// Fill
	glColor3f(r, g, b);
	glBegin(GL_POLYGON);
	for (int i = 0; i < nx; i++)
		glVertex2f(x + xx[i] * radius, y + yx[i] * radius);
	glEnd();

	/// Outline
	glColor3f(0, 0, 0);
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < nx; i++)
		glVertex2f(x + xx[i] * radius, y + yx[i] * radius);
	glEnd();

	glPopMatrix();
}

void PF_Agent::_drawCircle(Real radius, Real xc, Real yc, Real r, Real g, Real b) const
{
	const Real deg2rad = M_PI / 180;
	GLfloat red = 1.3 * r;
	GLfloat green = 1.3 * g;
	GLfloat blue = 1.3 * b;
	GLfloat lightColorExit[] = { red, green, blue, 1 };
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColorExit);

	glPushMatrix();

	glColor3f(r, g, b);
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < 360; i++)
	{
		Real degInRad = i * deg2rad;
		glVertex2f(xc + cos(degInRad) * radius, yc + sin(degInRad) * radius);
	}
	glEnd();

	glPopMatrix();
}

void PF_Agent::_paintTriangle(Real xc, Real yc, Real radius, Real direction, Real r, Real g, Real b) const
{
	GLfloat red = 1.3 * r;
	GLfloat green = 1.3 * g;
	GLfloat blue = 1.3 * b;
	GLfloat lightColorExit[] = { red, green, blue, 1 };

	Real halfLength = 1.0 * radius;
	Real vertex1[2] = { (Real) 0.5 * radius, -halfLength };
	Real vertex2[2] = { -(Real) 0.5 * radius, -halfLength };
	Real vertex3[2] = { (Real) 0.0, halfLength };

	direction -= 0.5 * M_PI;
	_rotate(vertex1, direction);
	_rotate(vertex2, direction);
	_rotate(vertex3, direction);

	glPushMatrix();
	glColor3f(r, g, b);

	glBegin(GL_TRIANGLES);
	glVertex2f(vertex3[0] + xc, vertex3[1] + yc);
	glVertex2f(vertex2[0] + xc, vertex2[1] + yc);
	glVertex2f(vertex1[0] + xc, vertex1[1] + yc);
	glEnd();

	glColor3f(0, 0, 0);
	glBegin(GL_LINE_LOOP);
	glVertex2f(vertex1[0] + xc, vertex1[1] + yc);
	glVertex2f(vertex2[0] + xc, vertex2[1] + yc);
	glVertex2f(vertex3[0] + xc, vertex3[1] + yc);
	glEnd();

	glPopMatrix();
}

void PF_Agent::_paintCone(Real xc, Real yc, Real radius, Real direction, Real r, Real g, Real b) const
{
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);

	glEnable(GL_LIGHT0);
	GLfloat lightpos[] = { 0, 0, 1, 0 };
	GLfloat lightColor[] = { r, g, b, 1 };
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);
	glLightfv(GL_LIGHT0, GL_POSITION, lightpos);

//	glEnable(GL_LIGHT1); // not sure why i have to turn this light on!
//	lightpos[2] = 1;
//	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor);
//	glLightfv(GL_LIGHT1, GL_POSITION, lightpos);

	GLfloat directionDegrees = 180.0 * direction / M_PI;

	glPushMatrix();
	glColor3f(r, g, b);
	glTranslated(xc, yc, 0);
	glTranslated(-radius * cos(direction), -radius * sin(direction), 0); // move back a little so that the COM coincides with (xc, yc)
	glRotated(90, 0.0, 1.0, 0.0); // rotate about
	glRotated(-directionDegrees, 1.0, 0.0, 0.0);
	glutSolidCone(0.5 * radius, 2 * radius, 50, 50);
	glPopMatrix();
}

void PF_Agent::_drawFullCircle(Real radius, Real xc, Real yc, Real r, Real g, Real b) const
{
	const Real deg2rad = M_PI / 180;

	glPushMatrix();

	glColor3f(r, g, b);
	glBegin(GL_POLYGON);
	for (int i = 0; i < 360; i++)
	{
		Real degInRad = i * deg2rad;
		glVertex2f(xc + cos(degInRad) * radius, yc + sin(degInRad) * radius);
	}
	glEnd();

	glPopMatrix();
}

void PF_Agent::_rotate(Real point[2], Real angle) const
{
	Real xnew = point[0] * cos(angle) - point[1] * sin(angle);
	Real ynew = point[0] * sin(angle) + point[1] * cos(angle);
	point[0] = xnew;
	point[1] = ynew;
}

#endif
