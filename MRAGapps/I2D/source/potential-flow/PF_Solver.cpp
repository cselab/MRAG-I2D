/*
 *  PF_Solver.cpp
 *  DipoleCode
 *
 *  Created by Alexia on 3/19/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#include <assert.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <numeric>
#include <algorithm>
#include "RL_QLearning.h"
#include "PF_AgentVector.h"
#include "PF_ObjectFactory.h"
#include "PF_Solver.h"
//#include "PF_Tracer.h"

#ifdef _RL_VIZ
#include "FotoCamera.h"
#ifdef __APPLE__
#include "GLUT/glut.h"
#else
#include <GL/gl.h>
#endif
#endif

using namespace PF;

PF_Solver::PF_Solver(int argc, const char ** argv) :
		parser(argc, argv), policy(NULL), collection(NULL), unitVelocity(0.0), ISLABFRAME(true), FITSELECT(1)
{
	// Parse input
	parser.set_strict_mode();

	SAVEFREQ = parser("-savefreq").asInt();
	RESTART = parser("-restart").asBool();
	XPOS = parser("-xpos").asDouble();
	YPOS = parser("-ypos").asDouble();
	CHARLENGTH = parser("-D").asDouble();
	ISCONTROLLED = parser("-isControlled").asBool();
	SMOOTH = parser("-smooth").asBool();
	LEARNINGTIME = parser("-learningTime").asDouble();
	FITNESSTIME = parser("-fitnessTime").asDouble();
	FITNESSSAVEFREQ = parser("-fitnessSaveFreq").asInt();
	FITNESSBUFFER = parser("-fitnessBuffer").asDouble();
	NAVG = parser("-navg").asInt();
	NLEARNINGLEVELS = parser("-nLearningLevels").asInt();
	INDIVIDUALFITNESS = parser("-individualFitness").asBool();
	LR = parser("-lr").asDouble();
	GREEDYEPS = parser("-greedyEps").asDouble();

	parser.unset_strict_mode();
	FITSELECT = parser("-fitselection").asInt();
	if (FITSELECT == 0) FITSELECT = 1;
	ISLABFRAME = parser("-islabframe").asBool();
	ISUSINGTRACERS = 0;
	ISUSINGTRACERS = parser("-isUsingTracers").asBool();

	parser.save_options();

	// Write fitness (overwrite, one number)
	FILE* fid;
	fid = fopen("fitness", "w");
	fprintf(fid, "%e \n", 1e10); // PARAMETERS OKAY -> BUILD INITIAL SCHOOL
	fclose(fid);

	// Decrease learning time at each level by half
	totalLearningTime = 0.0;
	for (int i = 0; i < NLEARNINGLEVELS; i++)
		totalLearningTime += LEARNINGTIME/(pow(2.0,(double) i));

	// Initialize counters
	nTimesLearned = 0;
	nTimesAveraged = 0;
	nTeribleRefreshes = 0;
	nCompletedAveraging = 0;
	nAveragingScrewUps = 0;

	// Initialize values
	currentFitness = 0;
	crashPenalty = 0;
	penaltyForBeingTooClose = 0.0;
	needsRefreshing = false;
	isUsingSoftReset = true;

	_prepareAgents();
	//if (ISUSINGTRACERS) _prepareTracers();
	_checkSettings();

	avgFitness = vector<double>(collection->numberOfAgents(), 0.0);
	currentIndividualFitness = vector<double>(collection->numberOfAgents(), 0.0);
}

PF_Solver::~PF_Solver()
{
	_dispose();
}

/**
 * Prepare agents from the factory file.
 */
void PF_Solver::_prepareAgents()
{
	assert(collection==NULL);

	map<string, vector<PF_Agent*> > shapesMap;
	PF_ObjectFactory factory(CHARLENGTH, XPOS, YPOS);
	factory.create(parser, shapesMap, LR, GREEDYEPS, &policy);

	collection = new PF_AgentVector(parser, LR, GREEDYEPS, SMOOTH, ISCONTROLLED, shapesMap);

	assert(collection!=NULL);
	collection->restartPolicy();

	unitVelocity = factory.getUnitVelocity();
}

/*
void PF_Solver::_prepareTracers()
{
	assert(tracers == NULL);
	map<string, vector<PF_Agent*> > shapesMap;

	// Create a vector of agents by the name of tracer
	shapesMap.clear();

	string name("PF_Tracer");
	int counterID = 1;

	{
		const Real dx = 0.025;
		const int nGrid = 1.0 / dx;

		for (int i = 0; i < nGrid; i++)
		{
			const Real x = (i + 0.5) * dx;
			for (int j = 0; j < nGrid; j++)
			{
				const Real y = (j + 0.5) * dx;
				const complex<Real> location(x, y);
				PF_Tracer * object = new PF_Tracer(parser, location, counterID);
				map<string, vector<PF_Agent*> >::iterator it = shapesMap.find(name);
				if (it == shapesMap.end())
					shapesMap[name] = vector<PF_Agent*>();
				shapesMap[name].push_back(object);
				counterID++;
			}
		}
	}

	// Assign tracers
	tracers = new PF_AgentVector(parser, 0, 0, false, false, shapesMap);
	assert(tracers != NULL);
}
*/

/**
 * Refresh when you reach an invalid situation.
 */
void PF_Solver::_refresh()
{
	individualFitnessesInWindow.clear();
	schoolFitnessInWindow.clear();
	timeInWindow.clear();
	timeStepsInWindow.clear();

	_dispose();
	_prepareAgents();
}

void PF_Solver::_softRefresh()
{
	individualFitnessesInWindow.clear();
	schoolFitnessInWindow.clear();
	timeInWindow.clear();
	timeStepsInWindow.clear();

	collection->resetToInitialCondition();

	//printf("REFRESH\n");
	//exit(0);
}

void PF_Solver::_dispose()
{
	if (collection != NULL)
	{
		delete collection;
		collection = NULL;
	}
	assert(collection==NULL);

	if (policy != NULL)
	{
		delete policy;
		policy = NULL;
	}

	assert(policy==NULL);
}

/**
 * Compute the velocity of all agents.
 */
void PF_Solver::_computeVelocity()
{
	vector<complex<Real> > velocities;
	vector<pair<Real, Real> > targets;
	map<pair<Real, Real>, PF_Agent*> targetsAgent;
	vector<Vortex> vortices;

	collection->getVortices(vortices);
	collection->getTargets(targets, &targetsAgent);

	for (int n = 0; n < (int) targets.size(); n++)
	{
		complex<Real> velocity(0, 0);
		Real xn = targets[n].first;
		Real yn = targets[n].second;
		for (int j = 0; j < (int) vortices.size(); j++)
		{
			Real xj = vortices[j].x;
			Real yj = vortices[j].y;
			Real gammaj = vortices[j].gamma;
			Real X = xn - xj;
			Real Y = yn - yj;
			velocity += (j == n) ? (Real) (0) : -gammaj * (Y + I * X) / (Real) (2 * M_PI * (X * X + Y * Y));
		}
		velocities.push_back(velocity); /// beware, these are conjugate velocities
	}

	vector<complex<Real> >::iterator itVel = velocities.begin();
	for (vector<pair<Real, Real> >::iterator itTarg = targets.begin(); itTarg != targets.end(); ++itTarg, ++itVel)
	{
		PF_Agent *agent = targetsAgent[*itTarg];
		agent->setVelocity(*itVel);
	}
}

/**
 * Compute the velocity of all tracers.
 */
void PF_Solver::_computeVelocityTracers()
{
	const Real exagerateFactor = 10.0;
	vector<complex<Real> > velocities;
	vector<pair<Real, Real> > tracerPoints;
	map<pair<Real, Real>, PF_Agent*> targetsTracer;
	vector<Vortex> vortices;

	collection->getVortices(vortices);
	tracers->getTargets(tracerPoints, &targetsTracer); /// @todo change this

	for (int n = 0; n < (int) tracerPoints.size(); n++)
	{
		complex<Real> velocity(0, 0);
		Real xn = tracerPoints[n].first;
		Real yn = tracerPoints[n].second;
		for (int j = 0; j < (int) vortices.size(); j++)
		{
			Real xj = vortices[j].x;
			Real yj = vortices[j].y;
			Real gammaj = vortices[j].gamma;
			Real X = xn - xj;
			Real Y = yn - yj;
			velocity += - exagerateFactor * gammaj * (Y + I * X) / (Real) (2 * M_PI * (X * X + Y * Y));
		}
		velocities.push_back(velocity); /// beware, these are conjugate velocities
	}

	vector<complex<Real> >::iterator itVel = velocities.begin();
	for (vector<pair<Real, Real> >::iterator itTarg = tracerPoints.begin(); itTarg != tracerPoints.end(); ++itTarg, ++itVel)
	{
		PF_Agent *agent = targetsTracer[*itTarg];
		agent->setVelocity(*itVel);
	}
}

/**
 * Set the time step based on minimum separation.
 *
 * @return
 */
double PF_Solver::_setTimeStep()
{
	// Create positions-agents map
	vector<pair<Real, Real> > position;
	map<pair<Real, Real>, PF_Agent*> positionAgent;
	collection->storePosition(position, &positionAgent);

	//	const Real dtMax(1e-3);
	const Real dtMax = 0.005 / unitVelocity;
	const Real dtMin = 5e-4 / unitVelocity ;
	const Real DMIN = CHARLENGTH / 2.0 / (2 * sqrt(2 * M_PI));
	const Real DMAX = 2.0 * CHARLENGTH / (2 * sqrt(2 * M_PI));

	vector<pair<Real, Real> > storeIndex;
	pair<Real, Real> coord;

	int size = position.size();
	Real minDist = numeric_limits<Real>::max();
	for (int n = 0; n < size - 1; n++)
	{
		Real x1 = position[n].first;
		Real y1 = position[n].second;
		for (int s = n + 1; s < size; s++)
		{
			Real x2 = position[s].first;
			Real y2 = position[s].second;
			Real dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
			if (dist < minDist)
			{
				storeIndex.clear();
				minDist = dist;
				coord.first = x1;
				coord.second = y1;
				storeIndex.push_back(coord);
				coord.first = x2;
				coord.second = y2;
				storeIndex.push_back(coord);
			}
		}
	}

	double dt = max(dtMin, min(dtMax, max((Real) 0, minDist - DMIN) * (dtMax - dtMin) / (DMAX - DMIN)));

	if ((Real) dt == dtMin)
		for (vector<pair<Real, Real> >::iterator it = storeIndex.begin(); it != storeIndex.end(); ++it)
		{
			PF_Agent *agent = positionAgent[*it];
			agent->reachedDtMin();
		}

	dt = dtMax;
	return (dt);
}

/**
 * Save the policies and the last time where we are able restart time.
 *
 * @param time
 */
void PF_Solver::_save(double time)
{
	collection->saveData(time);

	// Save time as well
	ofstream out("RestartTime.txt", ios::out | ios::trunc);
	out.precision(10);
	out << time << " ";
	out << nTimesLearned;
	out.close();
}

/**
 * Save all resets (different from restarts).
 * @param time
 */
void PF_Solver::_saveResetTime(double time)
{
	FILE* fid;
	fid = fopen("resets.txt", "a");
	fprintf(fid, "%e\n", time);
	fclose(fid);
}

#ifdef _RL_VIZ
void PF_Solver::_paintSquareCanvas(double edgeLength, double bottomLeftCorner[2], double bgColor[3], bool isOutlinedBlack)
{
	const double square[5][2] = {
			{bottomLeftCorner[0], bottomLeftCorner[1]},
			{bottomLeftCorner[0] + edgeLength, bottomLeftCorner[1]},
			{bottomLeftCorner[0] + edgeLength, bottomLeftCorner[1] + edgeLength},
			{bottomLeftCorner[0], bottomLeftCorner[1] + edgeLength},
			{bottomLeftCorner[0], bottomLeftCorner[1]}
	};

//	GLfloat lightColorExit[] = {(GLfloat)bgColor[0], (GLfloat)bgColor[1], (GLfloat)bgColor[2], 1};
//	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColorExit);
	glPushMatrix();

	glPushMatrix();
	glColor3f(bgColor[0], bgColor[1], bgColor[2]);
	glBegin(GL_POLYGON);
	for (int i = 0; i < 5; i++)
	{
		glVertex2f(square[i][0], square[i][1]);
	}
	glEnd();
	glPopMatrix();

	if (isOutlinedBlack)
	{
		glPushMatrix();
		glColor3f(0, 0, 0);
		glBegin(GL_LINE_LOOP);
		for (int i = 0; i < 5; i++)
		{
			glVertex2f(square[i][0], square[i][1]);
		}
		glEnd();
		glPopMatrix();
	}
}

void PF_Solver::_paintBox(Real edgeLength, Real bottomLeftCorner[2], Real bgColor[3])
{
	const double square[5][2] = {
			{bottomLeftCorner[0], bottomLeftCorner[1]},
			{bottomLeftCorner[0] + edgeLength, bottomLeftCorner[1]},
			{bottomLeftCorner[0] + edgeLength, bottomLeftCorner[1] + edgeLength},
			{bottomLeftCorner[0], bottomLeftCorner[1] + edgeLength},
			{bottomLeftCorner[0], bottomLeftCorner[1]}
	};

//	glPushMatrix();
	glLineWidth(1.0f);
	glColor3f(bgColor[0], bgColor[1], bgColor[2]);
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < 5; i++)
	{
		glVertex2f(square[i][0], square[i][1]);
	}
	glEnd();
//	glPopMatrix();
}

void PF_Solver::_paintSphere(Real x, Real y, Real radius, Real r, Real g, Real b) const
{
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
}

/**
 * Draw the agents on screen.
 *
 * @param dt
 * @param timeHere
 * @param initialTime
 */
void PF_Solver::_paint(const double time)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushAttrib(GL_ENABLE_BIT);

	Real edgeLength;
	Real centerOfSchool[2];
	_getTargetSchoolInfo(centerOfSchool, edgeLength); // todo remove this from here to prevent extra computation when painting
	Real bottomLeftCorner[2] = { centerOfSchool[0] - 0.5 * edgeLength, centerOfSchool[1] - 0.5 * edgeLength };
	Real lineColor[3] = {0.5, 0.5, 0.5};

	if (ISLABFRAME)
	{
//		_paintBox(edgeLength, bottomLeftCorner, lineColor);
		collection->paint();
		if (ISUSINGTRACERS) tracers->paint();
	}
	else
	{
		collection->paint(1.0/edgeLength, centerOfSchool);
		if (ISUSINGTRACERS) tracers->paint(1.0/edgeLength, centerOfSchool);
	}

	char * timeStr = new char[20];
	sprintf(timeStr, "t = %6.2f", time);
//	_paintText(0.05f, 0.05f, timeStr, 0.7f, 0.7f, 0.7f, 1.0f);

	glPopAttrib();
	glutSwapBuffers();
}

/**
 * Draw text on the screen.
 *
 * @param dt
 * @param timeHere
 * @param initialTime
 */
void PF::PF_Solver::_paintText(float x, float y, const char* text, float r, float g, float b, float a)
{
//	glShadeModel(GL_SMOOTH);
//	glEnable(GL_LIGHTING);
//	glEnable(GL_LIGHT0);
//	GLfloat lightpos[] = { 0, 0, 0.5, 0 };
//	glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
//	GLfloat lightColor[] = { r, g, b, 1 };
//	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);

	glPushMatrix();
	glColor3f(r, g, b);
	glTranslated(x, y, 0);

	if (!text || !strlen(text))
		return;
	bool blending = false;
	if (glIsEnabled(GL_BLEND))
		blending = true;
	glEnable(GL_BLEND);
	glColor4f(r, g, b, a);
	glRasterPos2f(x, y);
	while (*text)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *text);
		text++;
	}
	if (!blending)
		glDisable(GL_BLEND);

	glPopMatrix();
}
#endif


void PF::PF_Solver::_checkStatus(double time, double & timeAtLastValidRefresh)
{
	// LEARNING STATE
	if (nTimesLearned < NLEARNINGLEVELS)
	{
		if (time - timeAtLastValidRefresh > LEARNINGTIME)
		{
			nTimesLearned++;
			printf("** Saving policy!\n");
			collection->savePolicy();
			_saveResetTime(time);
			timeAtLastValidRefresh = time;

			if ((nTimesLearned == NLEARNINGLEVELS))
			{
				LR = 0.0;
				GREEDYEPS = 0.0;
			}
			else
			{
				LR /= 2.0;
				GREEDYEPS /= 2.0;
			}

//			LEARNINGTIME /= 2.0;

			_refresh(); // hard refresh

			printf("** Finished learning stage %d/%d. New values for LR = %e, GREEDYEPS = %e, LEARNINGTIME = %0.2f.\n", nTimesLearned, NLEARNINGLEVELS, LR, GREEDYEPS, LEARNINGTIME);

			needsRefreshing = true;
		}
		solverState = learning;
	}

	// AVERAGING STATE
	else if (nTimesAveraged < NAVG && ISCONTROLLED)
	{
		LR = 0.0;
		GREEDYEPS = 0.0;

		if (time - timeAtLastValidRefresh > FITNESSTIME * FITNESSBUFFER) // finished with one averaging run, save fitness
		{
			_computeTimeAveragedFitness(time);

			currentFitness = (nCompletedAveraging * currentFitness + timeAveragedSchoolFitness) / (nCompletedAveraging + 1); // running average
			crashPenalty = nAveragingScrewUps;

			// Write fitness (overwrite, one number)
			FILE* fid;
			fid = fopen("fitness", "w");
			fprintf(fid, "%e \n", currentFitness + (double)crashPenalty + penaltyForBeingTooClose);
			fclose(fid);

			// Calculate fitness using another way and print history to a file
			Real sum = 0.0;
			currentFitnesses.push_back(timeAveragedSchoolFitness);
			const int nAveragesSoFar = currentFitnesses.size();
			for (int i = 0; i < nAveragesSoFar; ++i)
			{
				sum += currentFitnesses[i];
			}
			const Real thisAverage = sum / nAveragesSoFar;

			Real sigma = 0;
			Real sigmasum = 0;
			for (int i = 0; i < nAveragesSoFar; ++i)
			{
				sigmasum += (currentFitnesses[i] - thisAverage) * (currentFitnesses[i] - thisAverage);
			}
			sigma = sqrt(sigmasum / nAveragesSoFar);

			// Write fitness history
			fid = fopen("fitnessStatsHistory", "a");
			fprintf(fid, "%d %e %e\n", nAveragesSoFar, thisAverage, sigma);
			fclose(fid);

			// Calculate individual fitness and write to a single file
			const int n = avgFitness.size();
			vector<pair<Real, Real> > coordinates;
			collection->getCoordinates(coordinates);

			fid = fopen("individualFitnesses", "w");
			for (int i = 0; i < n; ++i)
			{
				currentIndividualFitness[i] = (nCompletedAveraging * currentIndividualFitness[i] + avgFitness[i]) / (nCompletedAveraging + 1); // running average
				fprintf(fid, "%e %e %e\n", coordinates[i].first, coordinates[i].second, currentIndividualFitness[i]);
			}
			fclose(fid);

			// Calculate individual fitness statistics
			fid = fopen("individualFitnessStatsHistory", "a");
			currentIndividualFitnesses.push_back(avgFitness);
			fprintf(fid, "** Average %d\n", nAveragesSoFar);
			for (int j = 0; j < n; j++)
			{
				Real sum = 0.0;
				for (int i = 0; i < nAveragesSoFar; ++i)
				{
					const vector<double> &fitnessForThisSet = currentIndividualFitnesses[i];
					sum += fitnessForThisSet[j];
				}

				Real thisAverage = sum / nAveragesSoFar;
				Real sigma = 0;
				Real sigmasum = 0;
				for (int i = 0; i < nAveragesSoFar; ++i)
				{
					sigmasum += (currentFitnesses[i] - thisAverage) * (currentFitnesses[i] - thisAverage);
				}
				sigma = sqrt(sigmasum / nAveragesSoFar);
				fprintf(fid, "%e %e %e %e\n", coordinates[j].first, coordinates[j].second, thisAverage, sigma);
			}
			fclose(fid);

			_saveResetTime(time);
			timeAtLastValidRefresh = time;
			_refresh(); // hard resets here

			needsRefreshing = true;
			nTimesAveraged++;
			nCompletedAveraging++;

			printf("** Finished averaging completely %d times (out of %d, NAVG = %d).\n", nCompletedAveraging, nTimesAveraged, NAVG);
			printf("** Screwed up %d times so far (out of %d, NAVG = %d) \n", nAveragingScrewUps, nTimesAveraged, NAVG);
			printf("** LR = %e. FITNESSTIME = %0.2f, FITNESSBUFFER = %0.2f.\n", LR, FITNESSTIME, FITNESSBUFFER);
			printf("** currentFitness = %e, crashPenalty = %e, penaltyForBeingTooClose = %e.\n", currentFitness, (double)crashPenalty, penaltyForBeingTooClose);
			printf("** currentFitnessCheck = %e.\n", thisAverage);
			printf("** USING FITNESS %d!!! \n", FITSELECT);
		}

		solverState = averaging;
	}

	// EXIT STATES
	else if (!ISCONTROLLED) abort();

	else if (nTimesAveraged == NAVG)
	{
		printf("** EXITING: reached expected sim time %e\n", time);
		// Write fitness (overwrite, one number)
		crashPenalty = NAVG - nCompletedAveraging;
		FILE* fid;
		fid = fopen("fitness", "w");
		fprintf(fid, "%e \n", currentFitness + (double)crashPenalty + penaltyForBeingTooClose);
		fclose(fid);
		abort();
	}
}

void PF::PF_Solver::_checkStatusPostAction(double time, double dt, bool valid, double & timeAtLastValidRefresh, double & timeAtLastNonValidRefresh)
{
	if (ISCONTROLLED)
	{
		if (!valid || dt <= 1e-4) // || time step is too small -> probably collision
		{
			if (solverState == averaging)
			{
				nTimesAveraged++;
				nAveragingScrewUps++;
				_saveResetTime(time);
				timeAtLastValidRefresh = time;
				printf("** WARNING: Screwed up %d times so far (out of %d, NAVG = %d) \n", nAveragingScrewUps, nTimesAveraged, NAVG);
			}

			printf("** NON VALID: Refresh at t = %e using %s.\n", time, isUsingSoftReset ? "soft resets" : "hard resets");

			// CHECK IF CRASHING TOO MUCH
			if (solverState == learning)
			{
				const int nMaxTerribleRefresh = 500;
				const double minTimeBetweenRefresh = 10.0/unitVelocity;

				if (time - timeAtLastNonValidRefresh < minTimeBetweenRefresh)
				{
					nTeribleRefreshes++;
					printf("** That is %d terrible refreshes with minTimeBetweenRefresh = %0.4f\n", nTeribleRefreshes, minTimeBetweenRefresh);
				}

				if (nTeribleRefreshes == nMaxTerribleRefresh)
				{
					printf("** Reached nMaxTerribleRefresh = %d. Agents are probably close given the decision making time.\n", nMaxTerribleRefresh);
					FILE* fid;
					fid = fopen("fitness", "w");
					fprintf(fid, "%e \n", 1e8 - time); // CAN'T LEARN AT ALL -> ABORT (subtract time to rank better, if the simulation goes longer = better)
					fclose(fid);
					abort();
				}
			}

			if (isUsingSoftReset) // to reduce dumping to files
			{
				_softRefresh();
			}
			else
			{
				printf("** Saving policy!\n");
				collection->savePolicy();
				_saveResetTime(time);
				_refresh();
			}

			timeAtLastNonValidRefresh = time;
		}
	}
}

void PF::PF_Solver::_checkSettings()
{
	assert(SAVEFREQ >= 0.0);
	assert(XPOS >= 0.0 && XPOS <= 1.0);
	assert(YPOS >= 0.0 && YPOS <= 1.0);
	assert(CHARLENGTH >= 0.0 && CHARLENGTH <= 1.0);

	const double timeToReachOtherSide = 1.0 / (unitVelocity * CHARLENGTH);
	if (FITNESSTIME * FITNESSBUFFER > timeToReachOtherSide)
	{
		printf("** ABORTING: FITNESSTIME * FITNESSBUFFER > timeToReachOtherSide! Check your settings!\n");
		printf("** Lower FITNESSTIME,FITNESSBUFFER, or characteristic length.\n");
		printf("** unitVelocity = %e, CHARLENGTH = %e, FITNESSTIME = %e.\n", unitVelocity * CHARLENGTH, CHARLENGTH, FITNESSTIME);
		abort();
	}
}

Real PF_Solver::_computeTooClosePenalty()
{
	const Real distanceFactor = 1.0;
	const Real minDist = distanceFactor * CHARLENGTH;
	const Real penalizationFactor = 0.0;
	Real penalty = 0.0;

	vector<pair<Real, Real> > points;
	collection->getCoordinates(points);

	const int n = points.size();
	for (int i = 0; i < n; i++)
	{
		const Real xi = points[i].first;
		const Real yi = points[i].second;

		for (int j = i + 1; j < n; j++) // only count the upper corner
		{
			const Real xj = points[j].first;
			const Real yj = points[j].second;
			const Real dist = (xi - xj) * (xi - xj) + (yi - yj) * (yi - yj);

			if (dist < minDist * minDist)
				penalty += penalizationFactor * CHARLENGTH * CHARLENGTH / dist;
		}
	}
	return (penalty);
}

/**
 *
 * @param time
 * @param dt
 */
void PF::PF_Solver::_computeFitness(double time, double dt)
{
	collection->fitness();

	if (solverState == averaging) // Start saving fitnesses in deque if averaging
	{
		// Get all fitnesses
		vector<Real> allFitnesses;
		collection->getFitnessValues(allFitnesses);

		// Add to deques
		individualFitnessesInWindow.push_back(allFitnesses);
		timeInWindow.push_back(time);
		timeStepsInWindow.push_back(dt);

		// Check if front of deque is outside window, if yes, pop off
		while (time - timeInWindow.front() > FITNESSTIME)
		{
			individualFitnessesInWindow.pop_front();
			timeInWindow.pop_front();
			timeStepsInWindow.pop_front();
		}
	}
}


void PF::PF_Solver::_computeTimeAveragedFitness(double time)
{
	// Count the number of guys that you want to average
	const double actualWindow = timeInWindow.back() - timeInWindow.front();
	const int n = collection->numberOfAgents();

//	vector<double> avgFitness(n, 0.0);

	for (int i = 0; i < (int)individualFitnessesInWindow.size(); i++)
	{
		const vector<Real> fitnessHere = individualFitnessesInWindow[i];

		for (int j = 0; j < n; j++)
		{
			avgFitness[j] += fitnessHere[j]*timeStepsInWindow[i];
		}
	}

	int nInAverage = n;
	Real avgSchoolFitness = 0.0;
	for (int i = 0; i < n; i++)
	{
		avgFitness[i] /= actualWindow;

		if (collection[i].isAveraged == true)
		{
			avgSchoolFitness += avgFitness[i];
		}
		else
		{
			nInAverage--;
		}
	}
	avgSchoolFitness /= nInAverage;
	timeAveragedSchoolFitness = avgSchoolFitness;

//	// Write fitness history
//	FILE* fid;
//	fid = fopen("schoolFitnesshistory", "a");
//	fprintf(fid, "%e \n", avgSchoolFitness);
//	fclose(fid);
//
//	vector<pair<Real, Real> > coordinates;
//	collection->getCoordinates(coordinates);
//
//	// Write individual fitnesses
//	fid = fopen("individualFitnessesHistory", "a");
//	for (int i = 0; i < n; i++)
//		fprintf(fid, "%e %e %e\n", coordinates[i].first, coordinates[i].second, avgFitness[i]);
//	fclose(fid);

	/** COMMENT ME OUT! **/
//	fid = fopen("blah.txt", "a");
//	fprintf(fid, "%e ", time);
//	for (int i = 0; i < n; i++)
//		fprintf(fid, "%e ", avgFitness[i]);
//	fprintf(fid, "\n");
//	fclose(fid);
//
//	fid = fopen("actions.txt", "a");
//	fprintf(fid, "%e ", time);
//	const vector<Real> printMe = individualFitnessesInWindow.back();
//	for (int i = 0; i < n; i++)
//		fprintf(fid, "%e ", printMe[i]);
//	fprintf(fid, "\n");
//	fclose(fid);
//
//	fid = fopen("avgFitness.txt", "a");
//	fprintf(fid, "%e ", time);
//	fprintf(fid, "%e \n", avgSchoolFitness);
//	fclose(fid);
	/** COMMENT ME OUT! **/
}


void PF::PF_Solver::_computeError(double time)
{
	vector<Real> errors;
	collection->error();
	collection->getErrorValues(errors);

	// Compute average
	schoolError = accumulate(errors.begin(), errors.end(), 0.0);
	schoolError /= collection->numberOfAgents();

	FILE *fid;
	fid = fopen("error.txt", "a");
	fprintf(fid, "%e %e \n", time, schoolError);
	fclose(fid);
}


void PF::PF_Solver::_getTargetSchoolInfo(Real x[2], Real & edgeLength)
{
	const Real windowScale = 1.5; // todo

	vector<pair<Real, Real> > targetPoints;

	collection->getTargetPoints(targetPoints);
	vector<Real> xT;
	vector<Real> yT;

	int nT = targetPoints.size();

	for (int i = 0; i < nT; i++)
	{
		const pair<Real, Real> point = targetPoints[i];
		xT.push_back(point.first);
		yT.push_back(point.second);
	}

	const Real xmin = *(min_element(xT.begin(), xT.end()));
	const Real xmax = *(max_element(xT.begin(), xT.end()));
	const Real ymin = *(min_element(yT.begin(), yT.end()));
	const Real ymax = *(max_element(yT.begin(), yT.end()));
	const Real xrange = xmax - xmin;
	const Real yrange = ymax - ymin;

	edgeLength = xrange > yrange ? (windowScale * xrange) : (windowScale * yrange);

	x[0] = 0.5 * (xmax + xmin);
	x[1] = 0.5 * (ymax + ymin);

	schoolEdgeLength = edgeLength;
	schoolCenter[0] = x[0];
	schoolCenter[1] = x[1];
	schoolPoints = targetPoints;
}

/**
 * Main guts of the solver...
 */
void PF_Solver::run()
{
#ifdef _RL_VIZ
	const unsigned int framePerCycle = 3;
	const Real cycle = 1.0;
	const Real fotoDT = cycle/((Real)framePerCycle-1.0);
	Real fotoTimer = 0.0;
	FotoCamera foto;
#endif

	double time = 0.0;
	double timeAtLastValidRefresh = 0.0;
	double timeAtLastNonValidRefresh = 0.0;
	long unsigned int step_id = 0;

	penaltyForBeingTooClose = _computeTooClosePenalty();
	printf("** Penalty for being too close = %e.\n", penaltyForBeingTooClose);

	_save(time);

	while (true)
	{
		_checkStatus(time, timeAtLastValidRefresh);

		if (needsRefreshing) {
			needsRefreshing = false;
			continue;
		}

		profiler.push_start("DT");
		const double dt = _setTimeStep();
		profiler.pop_stop();

		profiler.push_start("CHOOSE");
		const bool valid = collection->choose(time);
		profiler.pop_stop();

		profiler.push_start("CHECK");
		_checkStatusPostAction(time, dt, valid, timeAtLastValidRefresh, timeAtLastNonValidRefresh);
		profiler.pop_stop();

		profiler.push_start("VEL");
		_computeVelocity();
#ifdef _RL_VIZ
		if (ISUSINGTRACERS) _computeVelocityTracers();
#endif
		profiler.pop_stop();

		profiler.push_start("UPDATE");
		collection->update(dt, time);
		collection->updatePostVelocityUpdate(dt, time);
		profiler.pop_stop();

		if (ISCONTROLLED)
		{
			profiler.push_start("REWARD");
			collection->reward(time);
			profiler.pop_stop();

			profiler.push_start("LEARN");
			collection->learn(time);
			profiler.pop_stop();
		}

		if (step_id % SAVEFREQ == 0)
		{
			printf("** Saving policy!\n");
			collection->savePolicy();
			_save(time);
			profiler.printSummary();
		}

		time += dt;
		step_id++;

		profiler.push_start("FITNESS");
		_computeFitness(time, dt);
		profiler.pop_stop();

		// Average only when needed
		if (step_id % FITNESSSAVEFREQ == 0 && solverState == averaging)
		{
			_save(time);
			_computeError(time);
			_computeTimeAveragedFitness(time);
		}

#ifdef _RL_VIZ
		if (ISUSINGTRACERS) tracers->update(dt, time);
		profiler.push_start("PAINT");
		fotoTimer += dt;
		if(fotoTimer > fotoDT || step_id==0)
		{
			fotoTimer = 0.0;
			_paint(time);
//			foto.shoot();
		}
		profiler.pop_stop();
#endif
	}
}

