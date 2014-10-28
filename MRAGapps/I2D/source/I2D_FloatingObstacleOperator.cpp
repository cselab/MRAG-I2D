/*
 * I2D_FloatingObstacleOperator.cpp
 *
 *  Created on: Sep 7, 2011
 *      Author: mgazzola
 */

#include "I2D_FloatingObstacleOperator.h"
#include "I2D_Clear.h"

namespace FloatingObstacleOperatorStuff
{

struct Mass
{
	Real density;
	map< int, Real>& local_mass;

	Mass(Real density, map< int, Real>& local_mass): density(density), local_mass(local_mass)
	{
	}

	Mass(const Mass& c): density(c.density), local_mass(c.local_mass)
	{
	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& out) const
	{
		const int n = FluidBlock2D::sizeX*FluidBlock2D::sizeY;

		FluidElement2D * const b = &out(0,0);
		Real total = 0;
		for(int i=0; i<n; i++)
			total += b[i].tmp;

		total *= density*info.h[0]*info.h[0];

		map< int, Real>::iterator it = local_mass.find(info.blockID);
		assert(it != local_mass.end());
		it->second = total;
	}
};

struct Diagnostic
{
	Real force[2];
	Real area;
	Real gridpoints;

	Diagnostic():area(0), gridpoints(0)
	{
		force[0] = force[1] = 0;
	}

	const Diagnostic& operator += (const Diagnostic& c)
	{
		force[0] += c.force[0];
		force[1] += c.force[1];
		area += c.area;
		gridpoints += c.gridpoints;

		return *this;
	}
};

struct ComputeDiagnostics
{
	Real lambda, Uinf[2];
	map<int, Diagnostic>& diagnostics;
	map<int, const VelocityBlock *>& desiredVels;

	ComputeDiagnostics(map< int, Diagnostic>& diagnostics, map<int, const VelocityBlock *>& desiredVels, Real lambda, const Real _Uinf[2]):
		diagnostics(diagnostics), desiredVels(desiredVels), lambda(lambda)
	{
		this->Uinf[0] = _Uinf[0];
		this->Uinf[1] = _Uinf[1];
	}

	ComputeDiagnostics(const ComputeDiagnostics& c): diagnostics(c.diagnostics), desiredVels(c.desiredVels), lambda(c.lambda)
	{
		this->Uinf[0] = c.Uinf[0];
		this->Uinf[1] = c.Uinf[1];
	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{
		map<int, const VelocityBlock *>::iterator it = desiredVels.find(info.blockID);

		if( it!=desiredVels.end() )
		{
			const VelocityBlock * targetVelBlock = it->second;

			map< int, Diagnostic>::iterator it = diagnostics.find(info.blockID);
			assert(it != diagnostics.end());
			Diagnostic& diag = it->second;

			diag.force[0] = diag.force[1]  = 0.0;
			diag.area = 0.0;
			diag.gridpoints = 0.0;

			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
				{
					Real p[2];
					info.pos(p,ix,iy);

					const Real Xs = b(ix,iy).tmp;
					const Real ux = b(ix,iy).u[0];
					const Real uy = b(ix,iy).u[1];
					const Real ux_target = targetVelBlock->u[0][iy][ix];
					const Real uy_target = targetVelBlock->u[1][iy][ix];

					diag.force[0] += (ux-ux_target)*Xs;
					diag.force[1] += (uy-uy_target)*Xs;

					diag.area += Xs;
					diag.gridpoints += Xs;
				}

			const Real dA = pow(info.h[0], 2);

			diag.force[0] *= dA*lambda;
			diag.force[1] *= dA*lambda;

			diag.area *= dA;

			diag.force[0] += diag.area*lambda*Uinf[0];
			diag.force[1] += diag.area*lambda*Uinf[1];
		}
	}
};

}


void I2D_FloatingObstacleOperator::_dist(const double x2[2], const double x1[2], double d[2]) const
{
	d[0] = x2[0] - x1[0];
	d[1] = x2[1] - x1[1];
}

double I2D_FloatingObstacleOperator::_angleVectors(const double v1[2], const double v2[2]) const
{
	const double anglev = atan2(v1[1],v1[0]) /M_PI*180.0;
	const double angled = atan2(v2[1],v2[0]) /M_PI*180.0;
	const double angle = anglev-angled;
	return (angle<0.0)?angle+360.0:angle;
}

bool I2D_FloatingObstacleOperator::_isRight(const double x1[2], const double dir[2], const double x2[2]) const
{
	const double alpha = -M_PI/2.0;
	const double n[2] = { dir[0]*cos(alpha) - dir[1]*sin(alpha), dir[0]*sin(alpha) + dir[1]*cos(alpha) };
	double d[2];
	_dist(x2,x1,d);
	return ((d[0]*n[0]+d[1]*n[1])>=0.0);
}

int I2D_FloatingObstacleOperator::_discretizeRange(const double value, const double minvalue, const double maxvalue, const int levels)
{
	assert(maxvalue>=minvalue);
	const double h = (maxvalue-minvalue)/levels;
	return  max(0,min(levels-1,(int)floor((value+minvalue)/h)));
}

void I2D_FloatingObstacleOperator::savePolicy(string name)
{
	if(policy!=NULL)
		if((*policy)!=NULL)
			(*policy)->save(name);
}

void I2D_FloatingObstacleOperator::restartPolicy(string name)
{
	if(policy!=NULL)
		if((*policy)!=NULL)
			(*policy)->restart(name);
}

void I2D_FloatingObstacleOperator::computeDragAndStuff(const Real time, const Real charLength, const Real charVel)
{
	if(charVel<=0.0){ printf("Something wrong with characteristic velocity, charVel=%e!\n", charVel); abort(); }
	if(charLength<=0.0){ printf("Something wrong with characteristic lenght, charLength=%e!\n", charLength); abort(); }

	typedef FloatingObstacleOperatorStuff::Diagnostic DiagType; // stores force, area and grid points of each object

	map<int, DiagType> diagnostics; // stores all objects

	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	for(map<int, const VelocityBlock *>::iterator it=desired_velocity.begin(); it!=desired_velocity.end(); it++)
		diagnostics[it->first] = DiagType(); // create information container (force, area, grid points) for each object

	FloatingObstacleOperatorStuff::ComputeDiagnostics getDiag(diagnostics, desired_velocity, penalization.getLambda(), Uinf);
	block_processing.process(vInfo, grid.getBlockCollection(), getDiag); // compute force for each object and assign it to object

	DiagType global;
	for(map< int, DiagType>::iterator it=diagnostics.begin(); it!=diagnostics.end(); it++)
		global += it->second;

	this->Cd = 2*global.force[0]/(pow(charVel, 2)*charLength);
	this->dimT = 2*charVel*time/charLength;
}

bool I2D_FloatingObstacleOperator::choose(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	bool valid = true;

	if(policy!=NULL)
	{
		assert((*policy)!=NULL);
		if(learningTimer==0.0 && status==Ready)
		{
			assert(status == Ready);
			vector<int> state;
			valid = mapState(state,_data);
			int action = (*policy)->selectAction(state);
			mapAction(action);
			state.push_back(action);
			(*policy)->setStateActionStart(state);
			integralReward = 0.0;
			learningTimer = t;
			status = Waiting;
		}
	}

	return valid;
}

void I2D_FloatingObstacleOperator::stopTest(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	if(policy!=NULL)
	{
		assert(status==Waiting);

		assert((*policy)!=NULL);
		if( t>(learningTimer+learningInterval) && status==Waiting)
		{
			assert(status == Waiting);
			(*policy)->setReward(integralReward);
			vector<int> state;
			mapState(state,_data);
			(*policy)->setStateEnd(state);
			status = Ready;
		}
	}
}

void I2D_FloatingObstacleOperator::learn(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data, string name)
{
	if(policy!=NULL)
	{
		assert((*policy)!=NULL);
		if( t>(learningTimer+learningInterval) && status==Ready)
		{
			assert(status == Ready);
			(*policy)->update(name);
			learningTimer = 0.0;
		}
	}
}

std::vector<Real> I2D_FloatingObstacleOperator::getMass()
{
	I2D_Clear cleaner;
	cleaner.clearTmp(grid);

	characteristic_function();

	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	map<int, Real> local_data;
	for(vector<BlockInfo>::iterator it=vInfo.begin(); it!=vInfo.end(); it++)
			local_data[it->blockID] = 0;

	// Compute mass
	const Real density = 1;
	FloatingObstacleOperatorStuff::Mass mass(density,local_data);
	block_processing.process(vInfo, coll, mass);

	Real totalmass = 0;
	for(map< int, Real>::iterator it=local_data.begin(); it!=local_data.end(); it++)
		totalmass += it->second;

	std::vector<Real> result;
	result.push_back(totalmass);

	return result;
}

