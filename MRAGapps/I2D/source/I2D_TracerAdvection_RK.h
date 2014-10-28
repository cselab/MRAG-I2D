/**
 * @file I2D_TracerAdvection_RK.cpp
 * @author Mattia Gazzola
 * @date Mar 1, 2012
 */
#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ParticleBlockLab.h"

struct I2D_TracerAdvection_RK
{
	map<long int, vector<Real> > & particles;
	map<I3, vector<long int> > & block2particles;

	Real dt, t;
	Real Uinf[2];

	int stencil_start[3], stencil_end[3];

	// Constructor
	// @param x
	// @param y
	// @param dt
	// @param Uinf
	I2D_TracerAdvection_RK(map<long int, vector<Real> > & particles, map<I3, vector<long int> > & block2particles, Real dt, Real Uinf[2]): t(0), particles(particles), block2particles(block2particles), dt(dt)
	{
		stencil_start[0] = stencil_start[1] = -4;
		stencil_end[0] = stencil_end[1] = +5;
		stencil_start[2] = 0;
		stencil_end[2] = 1;

		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
	}

	// Copy constructor
	// @param c
	I2D_TracerAdvection_RK(const I2D_TracerAdvection_RK& c): t(0), particles(c.particles), block2particles(c.block2particles), dt(c.dt)
	{
		stencil_start[0] = stencil_start[1] = -4;
		stencil_end[0] = stencil_end[1] = +5;
		stencil_start[2] = 0;
		stencil_end[2] = 1;

		Uinf[0] = c.Uinf[0];
		Uinf[1] = c.Uinf[1];
	}

	// Template to overload () operator for use in block_processing.process, see
	// I2D_CoreParticles.cpp > push for more information.
	// @param lab
	// @param info
	// @param out
	template <typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const
	{
		/// Interpolation stencil info
		const int VSX = -3;
		const int VSY = -3;
		const int VEX = B::sizeX + 3;
		const int VEY = B::sizeY + 3;
		const int KS = -1;
		const int KE = +3;

		/// Runge-Kutta factors
		const int nSubSteps = 4;
		const Real kFactor[nSubSteps] = {1./6., 1./3., 1./3., 1./6.};
		const Real hFactor[nSubSteps] = {0.5, 0.5, 1.0, 0.0};
		//const int nSubSteps = 2;
		//const Real kFactor[nSubSteps] = {0.0, 1.0};
		//const Real hFactor[nSubSteps] = {0.5, 0.0};

		/// Block info
		const Real dx = 1.0/pow(2.0, info.level); /// dx of the block
		const Real h = info.h[0];
		const int index[2] = {info.index[0], info.index[1]};
		const Real startBlock[2] = {index[0]*dx, index[1]*dx}; /// actual location of the start of the block
		const Real endBlock[2] = {startBlock[0]+dx, startBlock[1]+dx}; /// actual location of the end limits of the block
		const Real ori[2] = {info.origin[0],info.origin[1]};

		I3 currentBlock(info.index[0],info.index[1],info.level);
		const vector<long int> & currentParticles = block2particles[currentBlock];

		for(unsigned int j=0; j<currentParticles.size(); j++)
		{
			Real & x = particles[currentParticles[j]][0];
			Real & y = particles[currentParticles[j]][1];

			/// Check if this passive tracer is in the block, exit if not in block
			assert( x>=startBlock[0] && x<endBlock[0] && y>=startBlock[1] && y<endBlock[1] );

			Real xNext = x;
			Real yNext = y;
			for(int i = 0; i < nSubSteps; i++)
			{
				/// Get grid index where x lies relative to block (0,0) particle
				const Real xp[2] = { (xNext - ori[0])/h, (yNext - ori[1])/h }; /// x/h of the anchor point relative to the block
				const Real ap[2] = { floor(xp[0]), floor(xp[1]) };
				const int iap[2] = { (int)ap[0], (int)ap[1] };
				assert(iap[0] >= -2 && iap[0] < 33 && iap[1] >= -2 && iap[1] < 33);

				/// Compute weights and interpolate to get velocity
				Real weights[2][4];
				lab.pcore._computeWeights(xp, ap, weights);
				const int start[2] = { max(KS, VSX - iap[0]), max(KS, VSY - iap[1]) };
				const int end[2] = { min(KE, VEX - iap[0]), min(KE, VEY - iap[1]) };
				const Real vx = Uinf[0] + lab.pcore.template _sample<1>(weights, start, end, iap, lab);
				const Real vy = Uinf[1] + lab.pcore.template _sample<2>(weights, start, end, iap, lab);

				const Real kvx = dt*vx;
				const Real kvy = dt*vy;

				/// Update position this substep
				xNext = x + kvx*hFactor[i];
				yNext = y + kvy*hFactor[i];

				x += kvx*kFactor[i];
				y += kvy*kFactor[i];
			}
		}
	}
};
