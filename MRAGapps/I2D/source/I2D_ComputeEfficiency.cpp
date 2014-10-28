/*
 * I2D_ComputeEfficiency.cpp
 *
 *  Created on: Mar 8, 2012
 *      Author: mgazzola
 */

#include "I2D_ComputeEfficiency.h"
#include "I2D_VectorBlockLab.h"
#include "I2D_GradOfVector.h"
#include "I2D_FloatingObstacleVector.h"
#include "I2D_CarlingFishMorph.h"

namespace ComputeEfficiency
{

struct DissipationPower: I2D_GradOfVector_4thOrder
{
	Real viscosity, density, t;
	int stencil_start[3], stencil_end[3];
	map< int, Real>& local_dissipation;

	DissipationPower(Real viscosity, Real density, map< int, Real>& local_dissipation): viscosity(viscosity), density(density), t(0), local_dissipation(local_dissipation)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}

	DissipationPower(const DissipationPower & copy): viscosity(copy.viscosity), density(copy.density), t(0), local_dissipation(copy.local_dissipation)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}

	struct Assign { static inline void stream(Real& out, Real in) { out = in; } };
	struct Add2Square { static inline void stream(Real& out, Real in) { out += 2.0*in*in; } };


	template<typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const
	{
		const int n = FluidBlock2D::sizeX*FluidBlock2D::sizeY;

		// Make my little copy :)
		Real e[n];
		Real * const ext = &out.external_data[0][0];

		_dfdx_ptr<Assign, 1 >(lab, info, e); //grad_vx
		_dfdy_ptr<Assign, 0 >(lab, info, ext); //grad_uy

		for(int i=0; i<n; i++)
		{
			const Real grad_vx = e[i];
			const Real grad_uy = ext[i];
			e[i] = grad_vx*grad_vx + grad_uy*grad_uy + 2.0*grad_vx*grad_uy;
		}

		_dfdx_ptr<Add2Square, 0 >(lab, info, e); //grad_ux
		_dfdy_ptr<Add2Square, 1 >(lab, info, e); //grad_vy

		{
			FluidElement2D * const b = &out(0,0);
			Real total = 0;
			for(int i=0; i<n; i++)
				total += (1.0-b[i].tmp)*e[i];

			total *= (viscosity*density)*info.h[0]*info.h[0];

			map< int, Real>::iterator it = local_dissipation.find(info.blockID);
			assert(it != local_dissipation.end());
			it->second = total;
		}
	}
};


struct FluidKineticEnergy
{
	Real density;
	map< int, Real>& local_fluidKinetic;

	FluidKineticEnergy(Real density, map< int, Real>& local_fluidKinetic): density(density), local_fluidKinetic(local_fluidKinetic)
	{
	}

	FluidKineticEnergy(const FluidKineticEnergy& c): density(c.density), local_fluidKinetic(c.local_fluidKinetic)
	{
	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& out) const
	{
		const int n = FluidBlock2D::sizeX*FluidBlock2D::sizeY;

		FluidElement2D * const b = &out(0,0);
		Real total = 0;
		for(int i=0; i<n; i++)
		{
			const Real vx = b[i].u[0];
			const Real vy = b[i].u[1];
			total += (1.0-b[i].tmp)*(vx*vx+vy*vy);
		}

		total *= 0.5*density*info.h[0]*info.h[0];

		map< int, Real>::iterator it = local_fluidKinetic.find(info.blockID);
		assert(it != local_fluidKinetic.end());
		it->second = total;
	}
};


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

}


void I2D_ComputeEfficiency::compute(I2D_FloatingObstacleOperator * floatingObstacle, double t, double dt, double startTime, double endTime)
{
	assert(endTime >= startTime);

	if(t>=startTime && t<=endTime)
	{
		I2D_FloatingObstacleVector * vec = static_cast<I2D_FloatingObstacleVector*>(floatingObstacle);

		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		BoundaryInfo& binfo=grid.getBoundaryInfo();
		const BlockCollection<B>& coll = grid.getBlockCollection();
		map<int, Real> local_data;
		for(vector<BlockInfo>::iterator it=vInfo.begin(); it!=vInfo.end(); it++)
			local_data[it->blockID] = 0;

		// Get Masses
		vector<Real> masses = floatingObstacle->getMass();

		// Get positions
		string object("I2D_CarlingFishMorph");
		vector<I2D_FloatingObstacleOperator *> &agents = (vec->data)[object];
		assert(agents.size()==1);
		vector<I2D_FloatingObstacleOperator *>::iterator it=agents.begin();
		I2D_CarlingFishMorph * b = static_cast<I2D_CarlingFishMorph*>(*it);

		vector< vector<double> > infoShapes;
		vector<double> positions;
		positions.push_back(b->shape->xm);
		positions.push_back(b->shape->ym);
		infoShapes.push_back(positions);

		assert( masses.size() == (infoShapes.size()+1) );

		// Compute dissipated power
		for(map<int, Real>::iterator it=local_data.begin(); it!=local_data.end(); it++)
			it->second = 0;

		ComputeEfficiency::DissipationPower dissipation(viscosity,density,local_data);
		block_processing.process< I2D_VectorBlockLab< Streamer_Velocity, 2 >::Lab >(vInfo, coll, binfo, dissipation);

		Real dissipatedPower = 0;
		for(map< int, Real>::iterator it=local_data.begin(); it!=local_data.end(); it++)
			dissipatedPower += it->second;

		// Compute fluid kinetic energy
		for(map<int, Real>::iterator it=local_data.begin(); it!=local_data.end(); it++)
			it->second = 0;

		ComputeEfficiency::FluidKineticEnergy kinetic(density,local_data);
		block_processing.process(vInfo, coll, kinetic);

		Real fluidKineticEnergy = 0;
		for(map< int, Real>::iterator it=local_data.begin(); it!=local_data.end(); it++)
			fluidKineticEnergy += it->second;

		FILE * ppFile = fopen("efficiencies.txt", "a");
		assert(ppFile!=NULL);
		fprintf(ppFile,"%e %e %e ",t,dissipatedPower,fluidKineticEnergy);
		assert(masses.size()>=2);
		for(int i=0; i<masses.size()-1;++i)
			fprintf(ppFile,"%e %e %e ",masses[i],infoShapes[i][0],infoShapes[i][1]);

		fprintf(ppFile,"%e\n",masses[masses.size()-1]);

		// Cloase file
		fclose(ppFile);
	}
}



