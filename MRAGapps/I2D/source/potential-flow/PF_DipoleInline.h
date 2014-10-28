#include <complex>
#include <vector>
#include <map>

#include "PF_Agent.h"

namespace PF
{

class PF_DipoleInline: public PF_Agent
{
public:
	PF_DipoleInline(MRAG::ArgumentParser & parser, const Real _D, const Real _V, const Real _acc, const Real _radius, const complex<Real> _locationCenter, const Real _dir[2], const Real lr, const Real greedyEps, bool isSmooth, bool isControlled, int fitselect, const Real T = 0, const int _ID = 0, RL::RL_TabularPolicy **_policy = NULL, const int seed = 0);
	virtual ~PF_DipoleInline();

	virtual void update(const Real dt, const Real time, map<string, vector<PF_Agent*> > *collection = NULL);
	virtual void updatePostVelocityUpdate(const Real dt, const Real time, map<string, vector<PF_Agent*> > *collection = NULL);
	virtual void storePosition(vector<pair<Real, Real> > & position, map<pair<Real, Real>, PF_Agent*> *positionAgent = NULL);

	virtual void getCoordinates(vector<pair<Real, Real> >& coordinates);
	virtual void getVortices(vector<Vortex> & vortices);
	virtual void getNominalGammas(vector<Real> & nominalGammas);
	virtual void getTargets(vector<pair<Real, Real> > & targets, map<pair<Real, Real>, PF_Agent*> *targetsAgent = NULL);
	virtual void getForwardVelocities(vector<Real> & forwardVelocities);
	virtual void getFitnessValues(vector<Real> &fitnesses);
	virtual void getErrorValues(vector<Real> &errors);
	virtual void getTargetPoints(vector<pair<Real, Real> > &targetPoints);

	virtual void setupInitialCondition();
	virtual void resetToInitialCondition();

	virtual void setVelocity(complex<Real> velocityAgent);
	virtual void reachedDtMin();
	virtual void saveData(Real time);

	virtual void reward(const Real time, map<string, vector<PF_Agent*> > *collection = NULL);
	virtual bool mapState(vector<int> & state, map<string, vector<PF_Agent*> > *collection = NULL);
	virtual void mapAction(int action, Real time = 0);
	virtual void perform(const Real time);

	virtual void fitness();
	virtual void error();

#ifdef _RL_VIZ
	virtual void paint();
	virtual void paint(Real scale, const Real center[2]);
#endif

protected:
	Real D;
	const Real V;
	Real vx, vy;
	Real gamma;
	complex<Real> locationCenter;
	const complex<Real> locationCenter0;
	complex<Real> locationRightVortex;
	complex<Real> locationLeftVortex;
	Real alpha;
	Real gammaRight;
	Real gammaLeft;
	vector<complex<Real> > velocity;
	Real followFactor; // distance behind the target lattice point to follow
	vector<int> signature;
	Real cutoff;
	const Real maxDomainRadius;
	Real followedPoint[2]; // point swimmer creates for itself to travel so that he does not do something out of its capabilities
	Real targetPoint[2]; // point corresponding to RIGID advection of the lattice in the specified direction of travel
	Real dir[2];
	Real deltaStir;
	Real deltaSpeed;
	bool isTakingBadAction; // tells me whether agent took action right now or not
	int fitselect; // chooses a certain fitness

	Real wA, wD, wR, wL; // weights for cost function (accel, decel, right, left)
	Real wForBeingOnTarget, wForGoodAction; // weights for reward

	void _positionVortices();
	bool _isOutsideDomain(const Real myX[2]);
};

}
