#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <complex>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>

#include "PF_DipoleEight.h"

using namespace PF;

PF_DipoleEight::PF_DipoleEight(MRAG::ArgumentParser & parser, const Real _D, const Real _V, const Real _radius, const Real _acc, const complex<Real> _locationCenter, const Real _dir[2], const Real lr, const Real greedyEps, bool isSmooth, bool isControlled, const Real T, const int _ID, RL::RL_TabularPolicy **_policy, const int seed)
		: PF_Agent(parser, lr, greedyEps, isSmooth, isControlled, "PF_DipoleEight", _ID, _policy, seed), D(_D / (2 * sqrt(2 * M_PI))), V(_V), locationCenter(_locationCenter), locationCenter0(_locationCenter), maxDomainRadius(0.5)
{
	cutoff = 10 * D; /// note that D is actually \ell_n in the paper, see above, in the factory, D = \ell 2*sqrt(2*\pi)
	gamma = 2 * M_PI * V * D;
	deltaStir = 1. / (_radius * 4 * sqrt(2 * M_PI));
	deltaSpeed = _acc;
	followFactor = 2.0;

	Real IdirI = sqrt(_dir[0] * _dir[0] + _dir[1] * _dir[1]);
	assert(IdirI != 0.);
	dir[0] = _dir[0] / IdirI;
	dir[1] = _dir[1] / IdirI;

	setupInitialCondition();

	wForBeingOnTarget = 0.9;
	wForGoodAction = 0.1;

	// Assign the number levels to discretize the state
	signature.clear();
	if (policy != NULL) if ((*policy) != NULL)
	{
		signature.push_back(30); // n levels to discretize dipole distance to target point
		signature.push_back(30); // n levels to discretize dipole orientation with respect to target point
		signature.push_back(5); // number of actions

		if (!(*policy)->samedim(signature))
		{
			printf("** Policy dimension reset!\n");
			(*policy)->setdim(signature);
		}

		learningInterval = T;
		myAction = 0;
	}
}

PF_DipoleEight::~PF_DipoleEight()
{
	if (policy != NULL) if ((*policy) != NULL)
	{
		delete (*policy);
		(*policy) = NULL;
	}
}

void PF_DipoleEight::setupInitialCondition()
{
	isDead = false;

	locationCenter = locationCenter0;
	gammaRight = -gamma;
	gammaLeft = gamma;

	vx = V * dir[0];
	vy = V * dir[1];

	alpha = atan2(dir[1], dir[0]);

	_positionVortices();

	followedPoint[0] = real(locationCenter);
	followedPoint[1] = imag(locationCenter);

	circleRadius = 0.15;
	vphi = V / circleRadius;
	circleCenter[0] = 0.5;
	circleCenter[1] = 0.5;

	targetPoint[0] = real(locationCenter) + followFactor * dir[0] * D;
	targetPoint[1] = imag(locationCenter) + followFactor * dir[1] * D;

	velocity.clear();

	pathAction = straight;
}

/**
 * Place vortices after updating the center location.
 */
void PF_DipoleEight::_positionVortices()
{
	locationRightVortex = locationCenter - I * D * exp(I * alpha) / (Real) 2.;
	locationLeftVortex = locationCenter + I * D * exp(I * alpha) / (Real) 2.;
	assert(real(locationCenter) == real(locationCenter));
	assert(imag(locationCenter) == imag(locationCenter));
	assert(real(locationRightVortex) == real(locationRightVortex));
	assert(imag(locationRightVortex) == imag(locationRightVortex));
	assert(real(locationLeftVortex) == real(locationLeftVortex));
	assert(imag(locationLeftVortex) == imag(locationLeftVortex));
}

/**
 * Update velocities and positions of vortices.
 *
 * @param dt
 * @param time
 * @param collection
 */
void PF_DipoleEight::update(const Real dt, const Real time, map<string, vector<PF_Agent*> > *collection)
{
	complex<Real> newLocation;
	Real newAlpha;

	complex<Real> velocityRight = velocity.front();
	complex<Real> velocityLeft = velocity.back();

	complex<Real> conjVelocity = (Real) 0.5 * (velocityRight + velocityLeft);
	Real alphaDot = real((velocityRight - velocityLeft) * exp(I * alpha)) / D;

	complex<Real> vel = conj(conjVelocity);
	vx = real(vel);
	vy = imag(vel);

	if (isDead)
	{
		newLocation = locationCenter;
		newAlpha = alpha;
	}
	else
	{
		complex<Real> conjLocation = conj(locationCenter) + dt * conjVelocity;
		newLocation = conj(conjLocation);
		newAlpha = alpha + dt * alphaDot;
	}

	locationCenter = newLocation;

	assert(real(locationCenter) == real(locationCenter));
	assert(imag(locationCenter) == imag(locationCenter));

	alpha = newAlpha;
	while (alpha > M_PI)
		alpha -= 2 * M_PI;
	while (alpha < -M_PI)
		alpha += 2 * M_PI;

	_positionVortices();
}

/**
 * Update post velocity update (because you need velocities from all dipoles)
 *
 * @param dt
 * @param time
 * @param collection
 */
void PF_DipoleEight::updatePostVelocityUpdate(const Real dt, const Real time, map<string, vector<PF_Agent*> >* collection)
{
	// Previous update determined the average velocity of the school thus
	// allowing the school as a whole to move faster than the speed of the agent

	// isInBox
	boxLength = 0.35;
	{
		const Real boxCenter = 0.5;
		if (abs(targetPoint[0] - boxCenter) > 0.5 * boxLength || abs(targetPoint[1] - boxCenter) > 0.5 * boxLength) isInBox = false;
		else isInBox = true;
	}

	//	selectAction();
	if (isInBox)
	{
		pathAction = straight;
		circleCenter[0] = 0;
		circleCenter[1] = 0;
	}
	else if (!isInBox && pathAction == straight)
	{
		int randomness = rand();
		if (randomness % 2 == 0) pathAction = turnLeft;
		else pathAction = turnRight;

		// initializeTurn()
		{
			Real normal[2];
			if (pathAction == turnRight)
			{
				normal[0] = dir[1];
				normal[1] = -dir[0];
			}
			if (pathAction == turnLeft)
			{
				normal[0] = -dir[1];
				normal[1] = dir[0];
			}

			circleCenter[0] = targetPoint[0] + circleRadius * normal[0];
			circleCenter[1] = targetPoint[1] + circleRadius * normal[1];
		}

		const Real dx = targetPoint[0] - circleCenter[0];
		const Real dy = targetPoint[1] - circleCenter[1];

		phi = atan2(dy, dx);
	}

	//	peformAction();
	if (pathAction == straight)
	{
		targetPoint[0] += V * dir[0] * dt;
		targetPoint[1] += V * dir[1] * dt;

		followedPoint[0] = targetPoint[0] - followFactor * dir[0] * D; // follow one \ell behind at all times
		followedPoint[1] = targetPoint[1] - followFactor * dir[1] * D;
	}
	else if (pathAction == turnRight)
	{
		const Real lagDistance = followFactor * D;
		const Real deltaPhi = lagDistance / circleRadius;

		targetPoint[0] = circleCenter[0] + circleRadius * cos(phi);
		targetPoint[1] = circleCenter[1] + circleRadius * sin(phi);

		dir[0] = sin(phi);
		dir[1] = -cos(phi);

		followedPoint[0] = circleCenter[0] + circleRadius * cos(phi + deltaPhi);
		followedPoint[1] = circleCenter[1] + circleRadius * sin(phi + deltaPhi);

		phi -= vphi * dt;
	}
	else if (pathAction == turnLeft)
	{
		const Real lagDistance = followFactor * D;
		const Real deltaPhi = lagDistance / circleRadius;

		targetPoint[0] = circleCenter[0] + circleRadius * cos(phi);
		targetPoint[1] = circleCenter[1] + circleRadius * sin(phi);

		dir[0] = -sin(phi);
		dir[1] = cos(phi);

		followedPoint[0] = circleCenter[0] + circleRadius * cos(phi - deltaPhi);
		followedPoint[1] = circleCenter[1] + circleRadius * sin(phi - deltaPhi);

		phi += vphi * dt;
	}
	else
	{
		printf("No appropriate action chose for the path! Aborting...\n");
		abort();
	}
}

/**
 * Store position of the center.
 *
 * @param position
 * @param positionAgent
 */
void PF_DipoleEight::storePosition(vector<pair<Real, Real> > &position, map<pair<Real, Real>, PF_Agent*> *positionAgent)
{
	pair<Real, Real> coordinates;
	coordinates.first = real(locationCenter);
	coordinates.second = imag(locationCenter);
	position.push_back(coordinates);
}

/**
 *
 * @param coordinates
 */
void PF_DipoleEight::getCoordinates(vector<pair<Real, Real> >& coordinates)
{
	pair<Real, Real> point;
	point.first = real(locationCenter);
	point.second = imag(locationCenter);
	coordinates.push_back(point);
}

/**
 * Get the individual vortices comprising the dipole.
 *
 * @param vortices
 */
void PF_DipoleEight::getVortices(vector<Vortex> &vortices)
{
	Vortex pointR;
	pointR.x = real(locationRightVortex);
	pointR.y = imag(locationRightVortex);
	pointR.gamma = gammaRight;
	vortices.push_back(pointR);

	Vortex pointL;
	pointL.x = real(locationLeftVortex);
	pointL.y = imag(locationLeftVortex);
	pointL.gamma = gammaLeft;
	vortices.push_back(pointL);
}

/**
 * Get the nominal/initial strength of the vortex as defined by the factory.
 *
 * @param nominalGammas
 */
void PF_DipoleEight::getNominalGammas(vector<Real> & nominalGammas)
{
	nominalGammas.push_back(gamma);
}

/**
 * Get the locations of individual vortices of the dipole.
 *
 * @param targets
 * @param targetsAgent
 */
void PF_DipoleEight::getTargets(vector<pair<Real, Real> > &targets, map<pair<Real, Real>, PF_Agent*> *targetsAgent)
{
	pair<Real, Real> point;
	point.first = real(locationRightVortex);
	point.second = imag(locationRightVortex);
	targets.push_back(point);

	point.first = real(locationLeftVortex);
	point.second = imag(locationLeftVortex);
	targets.push_back(point);
}

/**
 * Get the velocities in the direction of desired travel.
 *
 * @param forwardVelocities
 */
void PF_DipoleEight::getForwardVelocities(vector<Real> &forwardVelocities)
{
	forwardVelocities.push_back(vx * dir[0] + vy * dir[1]);
}

/**
 *
 * @param velocityAgent
 */
void PF_DipoleEight::setVelocity(complex<Real> velocityAgent)
{
	if (velocity.size() == 2) velocity.clear();

	velocity.push_back(velocityAgent);
}

/**
 * Set the dipole dead and turn off its strengths.
 */
void PF_DipoleEight::reachedDtMin()
{
	isDead = true; // put dipole to dead and avoid integrating
	gammaRight = 0.0;
	gammaLeft = 0.0;
}

bool PF_DipoleEight::mapState(vector<int> &state, map<string, vector<PF_Agent*> > *collection)
{
	assert(collection!=NULL);
	assert((*collection).size()!=0);

	state.clear();

	// Gather my info (always double check with getInfo!!)
	const Real myX[2] = { real(locationCenter), imag(locationCenter) };
	const Real myV[2] = { vx, vy };

	// Levels
	const int levelX = signature[0];
	const int levelAngle = signature[1];

	// Get distance
	Real d[2] = { 0.0, 0.0 };
	_dist(targetPoint, myX, d); // remem

	// Constants
	const Real h = cutoff / (Real) levelX;
	const Real dangle = 360.0 / (Real) levelAngle;

	// Distance
	const Real IdI = fabs(sqrt(d[0] * d[0] + d[1] * d[1]));
	const int idxX = std::max(0, std::min(levelX - 1, (int) floor(IdI / h)));

	// Angle
	Real angle = _angleVectors(myV, d) + dangle / 2.0;
	angle = (angle > 360.0) ? angle - 360.0 : angle;
	const int idxangle = std::max(0, std::min(levelAngle - 1, (int) floor(angle / dangle)));

	// Prepare state vector
	state.push_back(idxX);
	state.push_back(idxangle);

	bool valid = true;
	if (ISCONTROLLED)
		if (_isOutsideDomain(myX) || IdI > (20 * cutoff)) /// || far far away from the followed point
			valid = false;

	return valid;
}

/**
 * Map the action given the case.
 * @param action
 * @param time
 */
void PF_DipoleEight::mapAction(int action, Real time)
{
	if (!isDead) switch (action)
	{
	case 0: // go straight
		  gammaRight = -gamma;
		  gammaLeft = gamma;
		  isTakingBadAction = false;
		  break;
	  case 1: // turn right
		  gammaRight = -(1.0 + deltaStir) * gamma;
		  gammaLeft = (1.0 - deltaStir) * gamma;
		  isTakingBadAction = true;
		  break;
	  case 2: // turn left
		  gammaRight = -(1.0 - deltaStir) * gamma;
		  gammaLeft = (1.0 + deltaStir) * gamma;
		  isTakingBadAction = true;
		  break;
	  case 3: // travel faster
		  gammaRight = -(1 + deltaSpeed) * gamma;
		  gammaLeft = (1 + deltaSpeed) * gamma;
		  isTakingBadAction = true;
		  break;
	  case 4: // travel slower
		  gammaRight = -(1 - deltaSpeed) * gamma;
		  gammaLeft = (1 - deltaSpeed) * gamma;
		  isTakingBadAction = false;
		  break;
	  default:
		  printf("Wrong mapAction\n");
		  abort();
		  break;
	  }
}

void PF_DipoleEight::perform(const Real time)
{
	Real diffTime = min(1.0, 0.1 + (time - learningTimer) / (0.5 * learningInterval));

	Real pgR = 0.;
	Real pgL = 0.;
	Real gR = 0.;
	Real gL = 0.;

	switch (myPreviousAction)
	{
	case 0: // go straight
		pgR = -gamma;
		pgL = gamma;
		break;
	case 1: // turn right
		pgR = -(1.0 + deltaStir) * gamma;
		pgL = (1.0 - deltaStir) * gamma;
		break;
	case 2: // turn left
		pgR = -(1.0 - deltaStir) * gamma;
		pgL = (1.0 + deltaStir) * gamma;
		break;
	case 3: // travel faster
		pgR = -(1.0 + deltaSpeed) * gamma;
		pgL = (1.0 + deltaSpeed) * gamma;
		break;
	case 4: // travel slower
		pgR = -(1.0 - deltaSpeed) * gamma;
		pgL = (1.0 - deltaSpeed) * gamma;
		break;
	default:
		printf("** Wrong previous mapAction!!!\n");
		abort();
		break;
	}

	switch (myAction)
	{
	case 0: // go straight
		gR = -gamma;
		gL = gamma;
		break;
	case 1: // turn right
		gR = -(1.0 + deltaStir) * gamma;
		gL = (1.0 - deltaStir) * gamma;
		break;
	case 2: // turn left
		gR = -(1.0 - deltaStir) * gamma;
		gL = (1.0 + deltaStir) * gamma;
		break;
	case 3: // travel faster
		gR = -(1.0 + deltaSpeed) * gamma;
		gL = (1.0 + deltaSpeed) * gamma;
		break;
	case 4: // travel slower
		gR = -(1.0 - deltaSpeed) * gamma;
		gL = (1.0 - deltaSpeed) * gamma;
		break;
	default:
		printf("** Wrong mapAction!!!\n");
		abort();
		break;
	}

	const Real diffGammaR = pgR - gR;
	const Real diffGammaL = pgL - gL;

	gammaRight = pgR - diffGammaR * diffTime;
	gammaLeft = pgL - diffGammaL * diffTime;
}

/**
 * Assign rewards.
 * @param time
 * @param collection
 */
void PF_DipoleEight::reward(const Real time, map<string, vector<PF_Agent*> > *collection)
{
	assert(status == Waiting);

	if (time > (learningTimer + learningInterval) && status == Waiting)
	{
		Real rewardDist = 0;
		{ /// Reward based on agent == follow point
			const Real myX[2] = { real(locationCenter), imag(locationCenter) };
			Real d[2] = { 0.0, 0.0 };
			_dist(followedPoint, myX, d);
			const Real rewardOnSpot = 1.0;
			const Real IdI2 = fabs(d[0] * d[0] + d[1] * d[1]);
			const Real a = (rewardOnSpot / (D * D));
			rewardDist = -a * IdI2 + rewardOnSpot;
		}

		Real rewardForGoodAction = 0;
		{ /// Reward based on taking good actions
			switch (myAction)
			{
			case 0: // go straight
				rewardForGoodAction = 0;
				break;
			case 1: // turn right
				rewardForGoodAction = -1;
				break;
			case 2: // turn left
				rewardForGoodAction = -1;
				break;
			case 3: // travel faster
				rewardForGoodAction = -1;
				break;
			case 4: // travel slower
				rewardForGoodAction = 1;
				break;
			default:
				printf("** Wrong mapAction!!!\n");
				abort();
				break;
			}
		}

		integralReward += wForBeingOnTarget * rewardDist + wForGoodAction * rewardForGoodAction;
	}
}

#ifdef _RL_VIZ
/**
 * Draw the dipoles onto the screen
 * @param time
 */
void PF_DipoleEight::paint()
{
	const Real x = real(locationCenter);
	const Real y = imag(locationCenter);

	glPushMatrix();
	glColor3f(0.9, 0.9, 0.9);
	glBegin(GL_LINE_LOOP);
	glVertex2f(0.5 + 0.5 * boxLength, 0.5 + 0.5 * boxLength);
	glVertex2f(0.5 - 0.5 * boxLength, 0.5 + 0.5 * boxLength);
	glVertex2f(0.5 - 0.5 * boxLength, 0.5 - 0.5 * boxLength);
	glVertex2f(0.5 + 0.5 * boxLength, 0.5 - 0.5 * boxLength);
	glVertex2f(0.5 + 0.5 * boxLength, 0.5 + 0.5 * boxLength);
	glEnd();
	glPopMatrix();

//	_paintSphere(circleCenter[0], circleCenter[1], 2.0 * D, 0.7, 0.2, 0.2);
//	_paintSphere(targetPoint[0], targetPoint[1], 2.0 * D, 1.0, 1.0, 1.0);
	_drawFullCircle(2.0 * D, targetPoint[0], targetPoint[1], 1.0, 0.0, 0.0);
	_drawCircle(2.0 * D, targetPoint[0], targetPoint[1], 0.0, 0.0, 0.0);

//	_paintCone(x, y, 5.0 * D, alpha, 0.0, 0.741, 0.965); // swimmer size is actually approx 5 times the separating distince D, see constructor
	_paintTriangle(x, y, 5.0 * D, alpha, 0.0, 1.0, 0.0);
}

void PF_DipoleEight::paint(Real scale, const Real center[2])
{
	const Real windowCenter[] =
	{	0.5, 0.5};
	const Real x = scale * (real(locationCenter) - center[0]) + windowCenter[0];
	const Real y = scale * (imag(locationCenter) - center[1]) + windowCenter[1];

	_paintTriangle(x, y, 5.0 * scale * D, alpha, 0.0, 1.0, 0.0);
//	_paintCone(x, y, 5.0 * scale * D, alpha, 0.0, 0.741, 0.965); // swimmer size is actually approx 5 times the separating distince D, see constructor
}

#endif

/**
 * Save data to data.txt.
 * @param time
 */
void PF_DipoleEight::saveData(Real time)
{
	//save the position of the centers of the dipoles
	ofstream out("data.txt", ios::out | ios::app);

	if (out)
	{
		Real velDipole = sqrt(vx * vx + vy * vy);
		out << time << " " << real(locationCenter) << " " << imag(locationCenter) << " " << real(locationRightVortex) << " " << imag(locationRightVortex) << " " << real(locationLeftVortex) << " "
			<< imag(locationLeftVortex) << " " << alpha << " " << gammaRight << " " << gammaLeft << " " << D << " " << velDipole << " " << targetPoint[0] << " " << targetPoint[1] << endl;
		out.close();
	}
	else cerr << "** Problem opening the file !" << endl;
}

/**
 * Compute instantaneous fitness for this fucker.
 */
void PF::PF_DipoleEight::fitness()
{
	fitnessValue = (abs(gammaLeft + gammaRight) + abs(gammaLeft) + abs(gammaRight) - 2 * gamma) / gamma;
}

void PF::PF_DipoleEight::getFitnessValues(vector<Real>& fitnesses)
{
	fitnesses.push_back(fitnessValue);
}

/**
 * Compute the instantaneous error.
 */
void PF::PF_DipoleEight::error()
{
	const Real myX[2] = { real(locationCenter), imag(locationCenter) };
	Real d[2] = { 0.0, 0.0 };
	_dist(followedPoint, myX, d);
	errorValue = (d[0] * d[0] + d[1] * d[1]) / (D * D);
}

void PF::PF_DipoleEight::getErrorValues(vector<Real>& errors)
{
	errors.push_back(errorValue);
}

/**
 * Get target points and store in vector.
 * @param targetPoints
 */
void PF::PF_DipoleEight::getTargetPoints(vector<pair<Real, Real> > &targetPoints)
{
	pair<Real, Real> point;
	point.first = targetPoint[0];
	point.second = targetPoint[1];
	targetPoints.push_back(point);
}

/**
 * Used for soft resets to avoid dumping a policy.
 */
void PF_DipoleEight::resetToInitialCondition()
{
	setupInitialCondition();
}

/**
 * Check if outside the specified domain.
 *
 * @param myX
 * @return
 */
bool PF_DipoleEight::_isOutsideDomain(const Real myX[2])
{
	if ((myX[0] > 1 || myX[0] < 0) || (myX[1] > 1 || myX[1] < 0)) return (true);
	else return (false);
}

