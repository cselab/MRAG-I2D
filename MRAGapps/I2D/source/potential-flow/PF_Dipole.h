/*
 *  PF_Dipole.h
 *  DipoleCode
 *
 *  Created by Alexia on 3/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#include <complex>
#include <vector>
#include "PF_Agent.h"

namespace PF
{

class PF_Dipole : public PF_Agent
{
public:
	//constructor & destructor
	PF_Dipole(MRAG::ArgumentParser & parser, const Real _D, const Real _V, const complex <Real> _locationCenter, const Real _alpha, const Real lr, const Real greedyEps, bool isSmooth, bool isControlled, const Real T = 0, const int _ID = 0, RL::RL_TabularPolicy **_policy = NULL);
	virtual ~PF_Dipole();
	
	//methods
	virtual void update(const Real dt, const Real time, map < string, vector <PF_Agent*> > *collection = NULL);
	virtual void storePosition(vector <pair <Real,Real> > &position, map < pair <Real,Real>, PF_Agent*> *positionAgent = NULL);
	virtual void getVortices(vector <Vortex> &vortices);
	virtual void getTargets(vector <pair <Real,Real> > &targets, map < pair <Real,Real>, PF_Agent*> *targetsAgent = NULL);
	virtual void setVelocity(complex <Real> velocityAgent);
	virtual void reachedDtMin();
	virtual void saveData(Real time);

	virtual void reward(const Real time, map < string, vector <PF_Agent*> > *collection = NULL);
	virtual bool mapState(vector <int> &state, map < string, vector <PF_Agent*> > *collection = NULL);
	virtual void mapAction(int action, Real time = 0);
	virtual void perform(const Real time = 0);
	
#ifdef _RL_VIZ
	virtual void paint();
	virtual void paint(Real scale, const Real center[2]) {};
#endif

protected:
	//attributs
	const Real D;
	const Real V;
	Real vx, vy;
	Real gamma;
	complex <Real> locationCenter;
	complex <Real> locationRightVortex;
	complex <Real> locationLeftVortex;
	Real alpha;
	Real gammaRight;
	Real gammaLeft;
	vector <complex <Real> > velocity;
	bool deadDipole;
	
	vector<int> signature;
	const Real maxDomainRadius;
	Real cutoff;
	
	void _positionVortices();
	void _mapForMaxRadius(int &idxDist, int &idxAngle, map< string, vector<PF_Agent*> > *collection);
	Real _angleVectors(const Real v1[2]) const;
};

}/*namespace PF*/
