#include <assert.h>
#include <math.h>

#include <complex>
#include <iostream>
#include <fstream>

#include "PF_DipoleInline.h"

using namespace PF;

PF_DipoleInline::PF_DipoleInline(MRAG::ArgumentParser & parser, const Real _D, const Real _V, const Real _radius, const Real _acc, const complex<Real> _locationCenter, const Real _dir[2], const Real lr, const Real greedyEps, bool isSmooth, bool isControlled, int fitselect, const Real T, const int _ID, RL::RL_TabularPolicy **_policy, const int seed)
		: PF_Agent(parser, lr, greedyEps, isSmooth, isControlled, "PF_DipoleInline", _ID, _policy, seed), D(_D / (2 * sqrt(2 * M_PI))), V(_V), locationCenter(_locationCenter), locationCenter0(_locationCenter), maxDomainRadius(0.5), fitselect(fitselect)
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
//			signature.push_back(5); // n levels to discretize speed between dipole and the group
//			signature.push_back(5); // n levels to discretize distance from follow to target point
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

PF_DipoleInline::~PF_DipoleInline()
{
	if (policy != NULL) if ((*policy) != NULL)
	{
		delete (*policy);
		(*policy) = NULL;
	}
}

void PF_DipoleInline::setupInitialCondition()
{
	isDead = false;

	locationCenter = locationCenter0;
	gammaRight = -gamma;
	gammaLeft = gamma;

	vx = V * dir[0];
	vy = V * dir[1];

	alpha = atan2(dir[1], dir[0]);
	;

	_positionVortices();

	followedPoint[0] = real(locationCenter);
	followedPoint[1] = imag(locationCenter);

	targetPoint[0] = real(locationCenter) + followFactor * dir[0] * D;
	targetPoint[1] = imag(locationCenter) + followFactor * dir[1] * D;

	velocity.clear();
}

/**
 * Place vortices after updating the center location.
 */
void PF_DipoleInline::_positionVortices()
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
void PF_DipoleInline::update(const Real dt, const Real time, map<string, vector<PF_Agent*> > *collection)
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
void PF_DipoleInline::updatePostVelocityUpdate(const Real dt, const Real time, map<string, vector<PF_Agent*> >* collection)
{
	// Previous update determined the average velocity of the school thus
	// allowing the school as a whole to move faster than the speed of the agent

//	Real totVx = 0.0;
//	Real totVy = 0.0;
//	vector<PF_Agent *> &agents = (*collection)["PF_DipoleInline"];
//
//	for (vector<PF_Agent *>::iterator it = agents.begin(); it != agents.end(); ++it)
//	{
//		PF_DipoleInline * b = static_cast<PF_DipoleInline*>(*it);
//		totVx += b->vx;
//		totVy += b->vy;
//	}
//
//	totVx /= (Real) agents.size();
//	totVy /= (Real) agents.size();
//
//	const Real modvLattice = totVx * dir[0] + totVy * dir[1];
//	targetPoint[0] += modvLattice * dir[0] * dt; // not used here
//	targetPoint[1] += modvLattice * dir[1] * dt; // not used here

	targetPoint[0] += V * dir[0] * dt;
	targetPoint[1] += V * dir[1] * dt;

	followedPoint[0] = targetPoint[0] - followFactor * dir[0] * D; // follow one \ell behind at all times
	followedPoint[1] = targetPoint[1] - followFactor * dir[1] * D;
}

/**
 * Store position of the center.
 *
 * @param position
 * @param positionAgent
 */
void PF_DipoleInline::storePosition(vector<pair<Real, Real> > &position, map<pair<Real, Real>, PF_Agent*> *positionAgent)
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
void PF_DipoleInline::getCoordinates(vector<pair<Real, Real> >& coordinates)
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
void PF_DipoleInline::getVortices(vector<Vortex> &vortices)
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
void PF_DipoleInline::getNominalGammas(vector<Real> & nominalGammas)
{
	nominalGammas.push_back(gamma);
}

/**
 * Get the locations of individual vortices of the dipole.
 *
 * @param targets
 * @param targetsAgent
 */
void PF_DipoleInline::getTargets(vector<pair<Real, Real> > &targets, map<pair<Real, Real>, PF_Agent*> *targetsAgent)
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
void PF_DipoleInline::getForwardVelocities(vector<Real> &forwardVelocities)
{
	forwardVelocities.push_back(vx * dir[0] + vy * dir[1]);
}

/**
 *
 * @param velocityAgent
 */
void PF_DipoleInline::setVelocity(complex<Real> velocityAgent)
{
	if (velocity.size() == 2) velocity.clear();

	velocity.push_back(velocityAgent);
}

/**
 * Set the dipole dead and turn off its strengths.
 */
void PF_DipoleInline::reachedDtMin()
{
	isDead = true; // put dipole to dead and avoid integrating
	gammaRight = 0.0;
	gammaLeft = 0.0;
}

bool PF_DipoleInline::mapState(vector<int> &state, map<string, vector<PF_Agent*> > *collection)
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
		if (_isOutsideDomain(myX) || IdI > (8 * cutoff)) /// || far far away from the followed point
			valid = false;

	return (valid);
}

/**
 * Map the action given the case.
 * @param action
 * @param time
 */
void PF_DipoleInline::mapAction(int action, Real time)
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

void PF_DipoleInline::perform(const Real time)
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
void PF_DipoleInline::reward(const Real time, map<string, vector<PF_Agent*> > *collection)
{
	assert(status==Waiting);

	if (time > (learningTimer + learningInterval) && status == Waiting)
	{
		Real rewardDist = 0;
		{ // Reward based on agent == follow point
			const Real myX[2] = { real(locationCenter), imag(locationCenter) };
			Real d[2] = { 0.0, 0.0 };
			_dist(followedPoint, myX, d);
			const Real rewardOnSpot = 1.0;
			const Real IdI2 = fabs(d[0] * d[0] + d[1] * d[1]);
			const Real a = (rewardOnSpot / (D * D));
			rewardDist = -a * IdI2 + rewardOnSpot;
		}

//		Real rewardVel = 0;
//		{ // Reward based on having the same speed as your neighbors (speeds are discontinuous!!)
//			Real totVx = 0.0;
//			Real totVy = 0.0;
//
//			// Compute average velocity among all lattice agents
//			vector<PF_Agent *> &agents = (*collection)["PF_DipoleInline"];
//
//			for (vector<PF_Agent *>::iterator it = agents.begin(); it != agents.end(); ++it)
//			{
//				PF_DipoleInline * b = static_cast<PF_DipoleInline*>(*it);
//				totVx += b->vx;
//				totVy += b->vy;
//			}
//
//			totVx /= (Real) agents.size();
//			totVy /= (Real) agents.size();
//
//			const Real modvLattice = totVx * dir[0] + totVy * dir[1];
//			const Real IvI = vx * dir[0] + vy * dir[1];
//			const Real diffVelPercent = (fabs(IvI - modvLattice) / modvLattice) * (fabs(IvI - modvLattice) / modvLattice);
//
//			const Real rewardGoodVel = 1.0;
//			const Real a = (rewardGoodVel / (0.1 * 0.1));
//			rewardVel = -a * diffVelPercent + rewardGoodVel;
//		}

		Real rewardForGoodAction = 0;
		{ /// Reward based on taking good actions
//			rewardForGoodAction = isTakingBadAction ? 0.0 : 1.0; // before testing

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

//		Real rewardInRightDirection = 0;
//		{ /// Reward based on orientation with respect to where it wants to go
//			const Real angleDifference = alpha - alpha0; // or use actual bearing atan2(vy,vx);
//			const Real angleFactor = 10.0;
//			rewardInRightDirection = -angleFactor * abs(angleDifference) + 1.0;
//		}

		integralReward += wForBeingOnTarget * rewardDist + wForGoodAction * rewardForGoodAction;

//		// Sum up integral reward as a combination of the ones calculated above
//		const Real w[5] = { 1.0, 0, 0, 0.1, 0 }; // weights are hard coded right now
//		integralReward += w[0] * rewardDist + w[1] * rewardVel + w[2] * rewardTarget + w[3] * rewardForBadAction + w[4] * rewardInRightDirection;
//		printf("Rewards (dist, vel, target, action): (%7.4f, %7.4f, %7.4f, %7.4f, %7.4f) \n", rewardDist, rewardVel, rewardTarget, rewardNoAction, rewardInRightDirection);
	}
}

#ifdef _RL_VIZ
/**
 * Draw the dipoles onto the screen
 * @param time
 */
void PF_DipoleInline::paint()
{
	const Real x = real(locationCenter);
	const Real y = imag(locationCenter);

	if (isDead)
	{
//		_paintSphere(x, y, 2.5 * D, 1.0, 0.2, 0.1); // swimmer size is actually approx 5 times the separating distince D, see constructor
		_paintX(x, y, 5.0 * D, 1.0, 0.2, 0.1);
	}
	else
	{
		if (isAveraged == false)
		{
			_paintCone(x, y, 5.0 * D, alpha, 0.7, 0.7, 0.7); // swimmer size is actually approx 5 times the separating distince D, see constructor
//			_paintTriangle(x, y, 5.0 * D, alpha, 0.7, 0.7, 0.7); // swimmer size is actually approx 5 times the separating distince D, see constructor
		}
		else if (isAveraged == true)
		{
//			_paintCone(x, y, 5. * D, alpha, 0.0, 1.0, 0.0); // swimmer size is actually approx 5 times the separating distince D, see constructor
			_drawFullCircle(1.5 * D, targetPoint[0], targetPoint[1], 1.0, 0.0, 0.0);
			_drawCircle(1.5 * D, targetPoint[0], targetPoint[1], 0.0, 0.0, 0.0); // swimmer size is actually approx 5 times the separating distince D, see constructor
			_paintTriangle(x, y, 5.0 * D, alpha, 0.0, 1.0, 0.0); // swimmer size is actually approx 5 times the separating distince D, see constructor
		}
	}
}

void PF_DipoleInline::paint(Real scale, const Real center[2])
{
	const Real windowCenter[] = { 0.5, 0.5 };
	const Real x = scale * (real(locationCenter) - center[0]) + windowCenter[0];
	const Real y = scale * (imag(locationCenter) - center[1]) + windowCenter[1];

	if (isDead)
	{
		_paintSphere(x, y, 2.5 * scale * D, 1.0, 0.2, 0.1); // swimmer size is actually approx 5 times the separating distince D, see constructor
//		_paintX(x, y, 5.0 * scale * D, 1.0, 0.2, 0.1);
	}
	else if (isAveraged == false)
	{
		_paintCone(x, y, 5.0 * scale * D, alpha, 0.7, 0.7, 0.7); // swimmer size is actually approx 5 times the separating distince D, see constructor
//			_paintTriangle(x, y, 5.0 * D, alpha, 0.7, 0.7, 0.7); // swimmer size is actually approx 5 times the separating distince D, see constructor
	}
	else if (isAveraged == true)
	{
		_paintCone(x, y, 5. * scale * D, alpha, 0.0, 0.741, 0.965); // swimmer size is actually approx 5 times the separating distince D, see constructor
//					_paintTriangle(x, y, 5.0 * D, alpha, 0.0, 0.741, 0.965); // swimmer size is actually approx 5 times the separating distince D, see constructor
	}
}

#endif

/**
 * Save data to data.txt.
 * @param time
 */
void PF_DipoleInline::saveData(Real time)
{
	//save the position of the centers of the dipoles
	ofstream out("data.txt", ios::out | ios::app);

	if (out)
	{
		Real velDipole = sqrt(vx * vx + vy * vy);
		out << time << " " << real(locationCenter) << " " << imag(locationCenter) << " " << real(locationRightVortex) << " " << imag(locationRightVortex) << " " << real(locationLeftVortex) << " " << imag(locationLeftVortex) << " " << alpha << " " << gammaRight << " " << gammaLeft << " " << D << " " << velDipole << endl;
		out.close();
	}
	else cerr << "** Problem opening the file !" << endl;
}

/**
 * Compute instantaneous fitness for this fucker.
 */
void PF::PF_DipoleInline::fitness()
{
	Real fitnessSpeed = (abs(gammaLeft) + abs(gammaRight) - 2 * gamma) / gamma;
	Real fitnessTurn = abs(gammaLeft + gammaRight) / gamma;
	if (fitselect == 1) fitnessValue = fitnessSpeed + fitnessTurn;
	else if (fitselect == 2) fitnessValue = fitnessSpeed;
	else if (fitselect == 3) fitnessValue = fitnessTurn;
}

void PF::PF_DipoleInline::getFitnessValues(vector<Real>& fitnesses)
{
	fitnesses.push_back(fitnessValue);
}

/**
 * Compute the instantaneous error.
 */
void PF::PF_DipoleInline::error()
{
	const Real myX[2] = { real(locationCenter), imag(locationCenter) };
	Real d[2] = { 0.0, 0.0 };
	_dist(followedPoint, myX, d);
	errorValue = (d[0] * d[0] + d[1] * d[1]) / (D * D);
}

void PF::PF_DipoleInline::getErrorValues(vector<Real>& errors)
{
	errors.push_back(errorValue);
}

/**
 * Get target points and store in vector.
 * @param targetPoints
 */
void PF::PF_DipoleInline::getTargetPoints(vector<pair<Real, Real> > &targetPoints)
{
	pair<Real, Real> point;
	point.first = targetPoint[0];
	point.second = targetPoint[1];
	targetPoints.push_back(point);
}

/**
 * Used for soft resets to avoid dumping a policy.
 */
void PF_DipoleInline::resetToInitialCondition()
{
	setupInitialCondition();
}

/**
 * Check if outside the specified domain.
 *
 * @param myX
 * @return
 */
bool PF_DipoleInline::_isOutsideDomain(const Real myX[2])
{
	if ((myX[0] > 1 || myX[0] < 0) || (myX[1] > 1 || myX[1] < 0)) return (true);
	else return (false);
}

