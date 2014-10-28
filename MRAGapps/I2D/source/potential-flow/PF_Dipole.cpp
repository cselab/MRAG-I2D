/*
 *  PF_Dipole.cpp
 *  DipoleCode
 *
 *  Created by Alexia on 3/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#include "assert.h"
#include <iostream>
#include <fstream>
#include <complex>
#include <math.h>
#include "PF_Dipole.h"

using namespace PF;

PF_Dipole::PF_Dipole(MRAG::ArgumentParser & parser, const Real _D, const Real _V, const complex<Real> _locationCenter, const Real _alpha, const Real lr, const Real greedyEps, bool isSmooth, bool isControlled, const Real T, const int _ID, RL::RL_TabularPolicy **_policy)
		: PF_Agent(parser, lr, greedyEps, isSmooth, isControlled, "PF_Dipole", _ID, _policy), D(_D / (2 * sqrt(2 * M_PI))), V(_V), locationCenter(_locationCenter), alpha(_alpha), gammaRight(0), gammaLeft(0), deadDipole(false), maxDomainRadius(0.5)
{			
	cutoff = 10*D;
	gamma = 2*M_PI*V*D;

	vx = V * cos(alpha);
	vy = V * sin(alpha);

	_positionVortices();

	velocity.clear();
	gammaRight = -gamma;
	gammaLeft = gamma;

	// For the learning
	signature.clear();
	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			signature.push_back(5);
			signature.push_back(2);
			signature.push_back(8);
			//			signature.push_back(10);
			signature.push_back(3);

			if(!(*policy)->samedim(signature))
			{
				printf("** Policy dimension reset!\n");
				(*policy)->setdim(signature);
			}

			learningInterval = T;
			myAction = 0;
		}
}

PF_Dipole::~PF_Dipole()
{
	if(policy != NULL)
		if( (*policy) != NULL)
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

void PF_Dipole::_positionVortices()
{
	locationRightVortex = locationCenter-I*D*exp(I*alpha)/(Real)(2);
	locationLeftVortex = locationCenter+I*D*exp(I*alpha)/(Real)(2);
}

void PF_Dipole::update(const Real dt, const Real time, map < string, vector <PF_Agent*> > *collection)
{
	complex <Real> newLocation;
	Real newAlpha;

	complex <Real> velocityRight = velocity.front();
	complex <Real> velocityLeft = velocity.back();

	complex <Real> conjVelocity = (Real)(0.5)*(velocityRight+velocityLeft);
	Real alphaDot = real((velocityRight-velocityLeft)*exp(I*alpha))/D;

	complex <Real> vel = conj(conjVelocity);
	vx = real(vel);
	vy = imag(vel);

	if (deadDipole) 
	{
		newLocation = locationCenter;
		newAlpha = alpha;
	}
	else 
	{
		complex <Real> conjLocation = conj(locationCenter) + dt*conjVelocity;
		newLocation = conj(conjLocation);
		newAlpha = alpha + dt*alphaDot;
	}

	//update
	locationCenter = newLocation;
	alpha = newAlpha;
	while (alpha > M_PI)
		alpha -= 2*M_PI;
	while (alpha < -M_PI)
		alpha += 2*M_PI;

	//	Real alpha0 = 0.000001*M_PI;
	//	alpha += alpha0;
	_positionVortices();	

	velocity.clear();
}

void PF_Dipole::storePosition(vector <pair <Real,Real> > &position, map < pair <Real,Real>, PF_Agent*> *positionAgent)
{	
	pair <Real,Real> coordinates;

	if (deadDipole)
	{
		coordinates.first = INFINITY;
		coordinates.second = INFINITY;
	}
	else
	{
		coordinates.first = real(locationCenter);
		coordinates.second = imag(locationCenter);	
	}

	position.push_back(coordinates);
}

void PF_Dipole::getVortices(vector <Vortex> &vortices)
{
	Vortex point;
	point.x = real(locationRightVortex);
	point.y = imag(locationRightVortex);
	point.gamma = gammaRight;
	vortices.push_back(point);

	point.x = real(locationLeftVortex);
	point.y = imag(locationLeftVortex);
	point.gamma = gammaLeft;
	vortices.push_back(point);
}

void PF_Dipole::getTargets(vector <pair <Real,Real> > &targets, map < pair <Real,Real>, PF_Agent*> *targetsAgent)
{
	pair <Real,Real> point;
	point.first = real(locationRightVortex);
	point.second = imag(locationRightVortex);
	targets.push_back(point);

	point.first = real(locationLeftVortex);
	point.second = imag(locationLeftVortex);
	targets.push_back(point);
}

void PF_Dipole::setVelocity(complex <Real> velocityAgent)
{
	velocity.push_back(velocityAgent);
}

void PF_Dipole::reachedDtMin()
{
	assert(!deadDipole);

	//put dipole to dead
	deadDipole = true;
	gammaRight = 0.;
	gammaLeft = 0.;
}

bool PF_Dipole::mapState(vector <int> &state, map < string, vector <PF_Agent*> > *collection)
{
	assert(collection!=NULL);
	assert((*collection).size()!=0);

	state.clear();

	// Gather my info
	const Real center[2] = {0.5,0.5};
	Real x = real(locationCenter);
	Real y = imag(locationCenter);
	const Real myX[2] = {x-center[0],y-center[1]};

	// Constants
	const int levelX = signature[0];
	const int levelAngle = signature[2];
	const Real innerR = 0.1;
	const Real outerR = 0.3;
	bool valid = true;
	const Real newcutoff = (outerR-innerR)/2.0;
	const Real h = newcutoff / (Real)levelX;
	const Real dangle = 360.0 / (Real)levelAngle;

	Real dist = sqrt( myX[0]*myX[0] + myX[1]*myX[1] ); dist = (dist==0.0)?1.0:dist;
	Real radialDir[2] = {myX[0]/dist,myX[1]/dist};
	radialDir[0] = ((radialDir[0]==0.0) && (radialDir[1]==0.0))?1.0:radialDir[0];
	radialDir[1] = ((radialDir[0]==0.0) && (radialDir[1]==0.0))?0.0:radialDir[1];
	Real innerPoint[2] = {innerR*radialDir[0],innerR*radialDir[1]};
	Real outerPoint[2] = {outerR*radialDir[0],outerR*radialDir[1]};

	Real d[2] = {0.0,0.0};
	if( fabs(dist-innerR)<=fabs(dist-outerR) )
		_dist(innerPoint,myX,d);
	else
		_dist(outerPoint,myX,d);

	const Real IdI = sqrt(d[0]*d[0]+d[1]*d[1]);
	const Real anglev = alpha /M_PI*180.0;
	const Real angled = atan2(d[1],d[0]) /M_PI*180.0;
	Real angle = anglev-angled;
	angle = (angle<0.0)?angle+360.0:angle;

	const bool inside = (dist>innerR && dist<outerR)?true:false;
	const int idxX = std::max(0,std::min(levelX-1,(int)floor(IdI/h)));
	const int idxY = inside;
	const int idxangle = std::max(0,std::min(levelAngle-1,(int)floor(angle/dangle)));

	// Prepare state vector
	state.push_back(idxX);
	state.push_back(idxY);
	state.push_back(idxangle);

	//	// Get other dipoles
	//	vector<PF_Agent*> &agents = (*collection)["PF_Dipole"];
	//	vector <pair <Real,Real> > coord;
	//
	//	// Find the closest
	//	PF_Dipole * closestAgent = NULL;
	//	Real minDist = numeric_limits<Real>::max();
	//	for(vector<PF_Agent *>::iterator it=agents.begin(); it!=agents.end(); ++it)
	//	{
	//		PF_Dipole *b = static_cast<PF_Dipole*>(*it);
	//		coord.clear();
	//		b->storePosition(coord);
	//		const Real target[2] = {coord[0].first-center[0],coord[0].second-center[1]};
	//		Real d[2];
	//		_dist(target, myX, d);
	//		const Real IdI = _modv(d);
	//		if (IdI == 0.)
	//			continue;
	//		closestAgent = (IdI<minDist)?b:closestAgent;
	//		minDist = (IdI<minDist)?IdI:minDist;
	//	}
	//
	//	Real idxDist = 0;
	//
	//	if(closestAgent!=NULL)
	//	{
	//		// Discretize distance from target
	//		coord.clear();
	//		closestAgent->storePosition(coord);
	//		const Real target[2] = {coord[0].first,coord[0].second};
	//		Real d[2];
	//		_dist(target, myX, d);
	//		const Real signIdI = _modv(d);
	//		const Real IdI = fabs(signIdI);
	//		const bool deadEnd = IdI >= cutoff;
	//		idxDist = _discretize(IdI, 0, cutoff, signature[2], deadEnd);
	//	}
	//
	//	state.push_back(idxDist);

	{
		const Real dist = sqrt( myX[0]*myX[0] + myX[1]*myX[1] );
		if( dist > maxDomainRadius )
			valid = false;

		//		if(deadDipole)
		//			valid = false;
	}

	return valid;
}

void PF_Dipole::_mapForMaxRadius(int &idxDist, int &idxAngle, map< string, vector<PF_Agent*> > *collection)
{
	// Gather my info
	const Real center[2] = {0.5,0.5};
	Real x = real(locationCenter);
	Real y = imag(locationCenter);
	const Real myX[2] = {x,y};

	// Get target point position
	Real dist[2];
	_dist(myX, center, dist);
	Real dir[2];
	_normalize(dist,dir);
	const Real target[2] = { maxDomainRadius*dir[0] + center[0], maxDomainRadius*dir[1] + center[1] };

	// Discretize distance from target
	Real d[2];
	_dist(target, myX, d);
	const Real IdI = _modv(d);
	const bool deadEnd = IdI >= cutoff;
	idxDist = _discretize(IdI, 0, cutoff, signature[0], deadEnd);

	// Discretize distance to target
	const Real angle = _angleVectors(d);
	idxAngle = _discretize(angle, 0, 360, signature[1], deadEnd);
}

Real PF_Dipole::_angleVectors(const Real v1[2]) const
{
	Real anglev = atan2(v1[1],v1[0]) /M_PI*180.0;
	Real angled = alpha/M_PI*180.;
	Real angle = anglev-angled;

	return (angle < 0.0)?angle+360.0:angle;
}

void PF_Dipole::mapAction(int action, Real time)
{	
	//	int aa = 3;
	//	Real velDipole = sqrt(vx*vx+vy*vy);
	//	Real ampF = gamma*0.01;
	//	Real freqF = 0.047747;
	//	Real phaseF = 3*M_PI/4.;
	//	Real gammaF = ampF*sin(2*M_PI*freqF*time+phaseF);
	//	//gammaF = ampF*sin(phaseF);
	//
	//	Real st = freqF*(gamma+gammaF)/(velDipole*velDipole);
	//	printf("gamma %f\n", gamma);
	//	printf("gammaF %f\n", gammaF);
	//	printf("inside sinus %f\n", 2*M_PI*freqF*time+phaseF);
	//	printf("time %f\n", time);
	//	//assert(2*M_PI*freqF*time+phaseF < 3*M_PI/4.);
	//	printf("strouhal number %f\n", st);
	//	printf("circulation + %f\n", gamma+gammaF);
	//	printf("circulation - %f\n", -(gamma-gammaF));
	//	printf("dipole velocity %f\n", velDipole);
	//	printf("optimal frequency %f\n", 0.3*velDipole*velDipole/(gamma));
	//	printf("alpha %f\n", alpha);
	//	printf("*****************************\n");
	//	//assert(st >= 0.2 && st <= 0.4);

	const Real deltaStir = 0.5;

	if (!deadDipole)
		switch(action)
		{
		case 0 : //go straigth
			gammaRight = -gamma;
			gammaLeft = gamma;
			break;
		case 1 : //turn right
			gammaRight = -(1.0+deltaStir)*gamma;
			gammaLeft = (1.0-deltaStir)*gamma;
			break;
		case 2 : //turn left
			gammaRight = -(1.0-deltaStir)*gamma;
			gammaLeft = (1.0+deltaStir)*gamma;
			break;
			//		case 3 : //test flapping
			//			gammaRight = -(gamma-gammaF);
			//			gammaLeft = gamma+gammaF;
			//			break;
		default :
			printf("Wrong mapAction\n");
			abort();
			break;
		}

//	ofstream out("action.txt", ios::out | ios ::app);
//
//	if(out)
//	{
//		out << time << " " << action << " " << gammaRight << " " << gammaLeft << " " << myAction << " " << myPreviousAction << endl;
//		out.close();
//	}
//	else
//		cerr << "Problem opening the file !" << endl;
}

void PF_Dipole::perform(const Real time)
{
	Real diffTime = min(1.0,0.1+(time-learningTimer)/(0.5*learningInterval));

	Real pgR = 0.;
	Real pgL = 0.;
	Real gR = 0.;
	Real gL = 0.;

	const Real deltaStir = 0.5;

	switch(myPreviousAction)
	{
	case 0 : //go straight
		pgR = -gamma;
		pgL = gamma;
		break;
	case 1 : //turn right
		pgR = -(1.0+deltaStir)*gamma;
		pgL = (1.0-deltaStir)*gamma;
//		pgR = -(1.0+deltaSpeed)*gamma;
//		pgL = (1.0+deltaSpeed)*gamma;
		break;
	case 2 : //turn left
		pgR = -(1.0-deltaStir)*gamma;
		pgL = (1.0+deltaStir)*gamma;
//		pgR = -(1.0-deltaSpeed)*gamma;
//		pgL = (1.0-deltaSpeed)*gamma;
		break;
	default :
		printf("Wrong previous mapAction\n");
		abort();
		break;
	}

	switch(myAction)
	{
	case 0 : //go straight
		gR = -gamma;
		gL = gamma;
		break;
	case 1 : //turn right
		gR = -(1.0+deltaStir)*gamma;
		gL = (1.0-deltaStir)*gamma;
//		gR = -(1.0+deltaSpeed)*gamma;
//		gL = (1.0+deltaSpeed)*gamma;
		break;
	case 2 : //turn left
		gR = -(1.0-deltaStir)*gamma;
		gL = (1.0+deltaStir)*gamma;
//		gR = -(1.0-deltaSpeed)*gamma;
//		gL = (1.0-deltaSpeed)*gamma;
		break;
	default :
		printf("Wrong mapAction\n");
		abort();
		break;
	}

	const Real diffGammaR = pgR-gR;
	const Real diffGammaL = pgL-gL;

	gammaRight = pgR-diffGammaR*diffTime;
	gammaLeft = pgL-diffGammaL*diffTime;

//	ofstream out("action.txt", ios::out | ios ::app);
//
//	if(out)
//	{
//		out << time << " " << gammaRight << " " << gammaLeft << " " << diffGammaR << " " << myAction << " " << myPreviousAction << endl;
//		out.close();
//	}
//	else
//		cerr << "Problem opening the file !" << endl;

}

void PF_Dipole::reward(const Real time, map < string, vector <PF_Agent*> > *collection)
{
	assert(status==Waiting);

	if( time > (learningTimer+learningInterval) && status==Waiting)
	{
		Real reward = 0.0;
		const Real rewardRing = 1.0;
		const Real innerR = 0.1;
		const Real outerR = 0.3;
		const Real midR = (outerR+innerR)/2.0;
		Real x = real(locationCenter);
		Real y = imag(locationCenter);
		const Real xx[2] = {x,y};
		const Real dist = sqrt( (xx[0]-0.5)*(xx[0]-0.5) + (xx[1]-0.5)*(xx[1]-0.5)  );

		const Real lenght = outerR - innerR;
		const Real IdI2 = fabs(dist-midR);
		const Real a = (rewardRing/(lenght*lenght));
		reward += -a*IdI2 + rewardRing;

		//		const Real velDipole = sqrt(vx*vx+vy*vy);
		//		reward += (velDipole >= 0.6*V)?0.5:0;

		//		printf("velocity dipole %f\n", velDipole);
		//		printf("velocity min %f\n", 0.4*V);

		// Get other dipoles
		//		vector<PF_Agent*> &agents = (*collection)["PF_Dipole"];
		//		vector <pair <Real,Real> > coord;
		//		Real turn = 0.;
		//
		//		for(vector<PF_Agent *>::iterator it=agents.begin(); it!=agents.end(); ++it)
		//		{
		//			PF_Dipole *b = static_cast<PF_Dipole*>(*it);
		//			coord.clear();
		//			b->storePosition(coord);
		//			Real angle = b->alpha;
		//			const Real turnSens = ((coord[0].first > 0.5 && angle > 0.) || (coord[0].first < 0.5 && angle < 0.))?-1:1;
		//			turn += turnSens;
		//		}
		//
		//		const Real sensRight = (turn <= 0)?false:true;
		//		const Real sensDipole = ((xx[0] > 0.5 && alpha > 0.) || (xx[0] < 0.5 && alpha < 0.))?false:true;
		//		reward += (sensRight == sensDipole)?0.5:0;
		//		const Real rewardSens = (sensRight == sensDipole)?0.5:0;
		//
		//		printf("sens onthe cercle %f\n", turn);
		//		if (rewardSens == 0.)
		//		{
		//			printf("reward sens on the cercle %f\n", 0.);
		//			printf("sensDipole %b\n", sensDipole);
		//		}

		//		// Find the closest
		//		PF_Dipole * closestAgent = NULL;
		//		Real minDist = numeric_limits<Real>::max();
		//		for(vector<PF_Agent *>::iterator it=agents.begin(); it!=agents.end(); ++it)
		//		{
		//			PF_Dipole *b = static_cast<PF_Dipole*>(*it);
		//			coord.clear();
		//			b->storePosition(coord);
		//			const Real target[2] = {coord[0].first,coord[0].second};
		//			Real d[2];
		//			_dist(target, xx, d);
		//			const Real IdI = _modv(d);
		//			if (IdI == 0.)
		//				continue;
		//			closestAgent = (IdI<minDist)?b:closestAgent;
		//			minDist = (IdI<minDist)?IdI:minDist;
		//		}
		//
		//		if(closestAgent!=NULL)
		//		{
		//			// Discretize distance from target
		//			const Real rewardAvoid = 4.;
		//			coord.clear();
		//			closestAgent->storePosition(coord);
		//			const Real target[2] = {coord[0].first,coord[0].second};
		//			Real d[2];
		//			_dist(target, xx, d);
		//			const Real signIdI = _modv(d);
		//			reward += (signIdI<=D)?(rewardAvoid)/(D/3.)*signIdI+2*rewardAvoid/3.:rewardAvoid;
		//		}

		integralReward += reward;

		//		if(integralReward>rewardRing)
		//		{
		//			printf("reward=%e\n",reward);
		//			printf("integralReward=%e\n",integralReward);
		//		}
		//
		//		assert(integralReward<=rewardRing);
	}
}

#ifdef _RL_VIZ
void PF_Dipole::paint()
{
		const Real innerR = 0.1;
		const Real outerR = 0.3;
		const Real midR = (outerR+innerR)/2.0;

	_drawCircle(maxDomainRadius,0.5,0.5,0.0,0.0,0.0);
	_drawCircle(midR,0.5,0.5,1.0,0.0,0.0);

	const Real xr = real(locationRightVortex);
	const Real yr = imag(locationRightVortex);
	const Real xl = real(locationLeftVortex);
	const Real yl = imag(locationLeftVortex);

	_paintSphere(xl,yl,D/4.0,0.0,1.0,0.0);
	_paintSphere(xr,yr,D/4.0,0.0,1.0,0.0);

	//	_paintSphere(xr,yr,D/8.0,0.0,1.0,0.0);
	//	_paintSphere(xl,yl,D/8.0,0.0,1.0,0.0);
}
#endif

void PF_Dipole::saveData(Real time)
{
	//save the position of the centers of the dipoles
	ofstream out("data.txt", ios::out | ios ::app);

	if(out)
	{
		Real velDipole = sqrt(vx*vx+vy*vy);
		//		out << time << " " << real(locationCenter) << " " << imag(locationCenter) << " " << real(locationRightVortex) << " " << imag(locationRightVortex) << " " << real(locationLeftVortex) << " " << imag(locationLeftVortex) << " " << alpha << " " << D << " " << velDipole << endl;
		out << time << " " << real(locationCenter) << " " << imag(locationCenter) << " " << real(locationRightVortex) << " " << imag(locationRightVortex) << " " << real(locationLeftVortex) << " " << imag(locationLeftVortex) << " " << alpha << " " << gammaRight << " " << gammaLeft << " " << D << " " << velDipole << endl;
		//		out << real(locationCenter) << " " << imag(locationCenter) << endl;
		out.close();
		//printf("Data successfully saved!\n");
	}
	else
		cerr << "Problem opening the file !" << endl;
}
