/*
 *  I2D_CoreParticles.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ParticleKernels.h"

static const Real a3[4] = {-1./2.,3./2.,3./2.,-1./2.};
static const Real a2[4] = {5./2.,-5./2.,-5./2.,5./2.};
static const Real a1[4] = {-4.,0.,-0.,-4.};
static const Real a0[4] = {2.,1.,1.,2.};

class I2D_CoreParticles
{
	static const int KS = -1;
	static const int KE = +3;
	
	static const int VSX = -3;
	static const int VSY = -3;
	
	static const int VEX = B::sizeX + 3;
	static const int VEY = B::sizeY + 3;
	
	static const int PSX = -2;
	static const int PSY = -2;
	
	static const int PEX = B::sizeX + 2;
	static const int PEY = B::sizeY + 2;
	
	static const int NPX = PEX-PSX;
	static const int NPY = PEY-PSY;
	
	inline bool _is_close(Real a, Real b) const
	{
		const Real tol= sizeof(Real)<8 ? 5e-4 : 1e-4;
		
		if (!(fabs(a-b)<tol))
			printf("not ! %e %e\n", a-b,tol);
		//a/ssert(fabs(a-b)<tol);
		return true;//fabs(a-b)<tol;
	}
	
public:
	
	Real xparticles[NPY][NPX][2];
	Real omega_new[B::sizeY][B::sizeX];
#ifndef NDEBUG	
	Real pou[B::sizeY][B::sizeX];
#endif
	void inline _computeWeights(const Real xp[2], const Real ap[2], Real (weights[2])[4]) 
	{			
		for(int c=0; c<2; c++)
		{
			Real t[4] = {
				fabs(xp[c] - (ap[c] + -1)),
				fabs(xp[c] - (ap[c] + +0)),
				fabs(xp[c] - (ap[c] + +1)),
				fabs(xp[c] - (ap[c] + +2))
			};
			
			for(int i=0; i<4; i++)
				weights[c][i] = a0[i] + t[i]*(a1[i] + t[i]*(a2[i] + t[i]*a3[i]));
			
#ifndef NDEBUG
			assert(_is_close(weights[c][0] + weights[c][1] + weights[c][2] + weights[c][3], 1));
			assert(_is_close(weights[c][0], MP4::eval(t[0])));
			assert(_is_close(weights[c][1], MP4::eval(t[1])));
			assert(_is_close(weights[c][2], MP4::eval(t[2])));
			assert(_is_close(weights[c][3], MP4::eval(t[3])));
#endif
		}
		
#ifndef NDEBUG
		{
			Real unit_partition = 0;
			
			for(int sy=KS; sy<KE; sy++)
				for(int sx=KS; sx<KE; sx++)
					unit_partition += weights[0][sx-KS]*weights[1][sy-KS];
			
			assert (_is_close(unit_partition, 1));
		}
#endif
	}
	
	template< int component, typename LabVel>
	Real inline _sample(const Real (w[2])[4], const int start[2], const int end[2], const int iap[2], LabVel& lab)
	{
		static const int KS = -1;
		static const int KE = +3;
		
		Real sample = 0;
		
		for(int sy=start[1]; sy<end[1]; sy++)
		{
			const Real wy = w[1][sy-KS];
			
			for(int sx=start[0]; sx<end[0]; sx++)
				sample += wy*w[0][sx-KS] * lab.template get<component>(iap[0] + sx, iap[1] + sy);
		}
		
#ifndef NDEBUG
		{
			Real unit_partition = 0;
			
			for(int sy=KS; sy<KE; sy++)
				for(int sx=KS; sx<KE; sx++)
					unit_partition += w[0][sx-KS]*w[1][sy-KS];
			
			assert (_is_close(unit_partition, 1));
		}
#endif
		return sample;
	}
	
	void inline _scatter(const Real (w[2])[4], const int start[2], const int end[2], const int iap[2], Real value)
	{
		static const int KS = -1;
		static const int KE = +3;
		
		for(int sy=start[1]; sy<end[1]; sy++)
		{
			const Real wy = w[1][sy-KS];
			
			for(int sx=start[0]; sx<end[0]; sx++)
			{
				const Real wxwy =  wy*w[0][sx-KS];
				
				omega_new[iap[1]+sy][iap[0]+sx] += wxwy*value;
#ifndef NDEBUG			
				assert(iap[0]+sx>=0 && iap[0]+sx<B::sizeX);
				assert(iap[1]+sy>=0 && iap[1]+sy<B::sizeY);
				
				pou[iap[1]+sy][iap[0]+sx] += wxwy;
#endif
			}
		}
		
#ifndef NDEBUG
		{
			double unit_partition = 0;
			
			for(int sy=KS; sy<KE; sy++)
				for(int sx=KS; sx<KE; sx++)
					unit_partition += w[0][sx-KS]*w[1][sy-KS];
			
			assert (_is_close(unit_partition, 1));
		}
#endif
	}
	
	template<typename LabVel>
	void push(const BlockInfo& info, LabVel& lab, double dt, const Real Uinf[3])
	{
		static const int KS = -1;
		static const int KE = +3;
		
		const Real factor1 = dt*0.5/info.h[0];
		const Real factor2 = dt/info.h[0];
		
		for(int iy=PSY; iy<PEY; iy++)
			for(int ix=PSX; ix<PEX; ix++)
			{
				Real xp[2] = {
					ix + factor1*(Uinf[0] + lab.template get<1>(ix, iy)),
					iy + factor1*(Uinf[1] + lab.template get<2>(ix, iy))};
				
				Real ap[2] = {
					floor(xp[0]),
					floor(xp[1])};
				
				int iap[2] = {
					(int)ap[0],
					(int)ap[1] };
				
				Real weights[2][4];
				_computeWeights(xp, ap, weights);
				
				const int start[2] = {
					max(KS, VSX - iap[0]), 
					max(KS, VSY - iap[1])
				};
				
				const int end[2] = {
					min(KE, VEX - iap[0]),  
					min(KE, VEY - iap[1])
				};
				
				Real * const final_xp = xparticles[iy-PSY][ix-PSX];
				
				final_xp[0] = ix + factor2*(Uinf[0] + _sample<1>(weights, start, end, iap, lab));
				final_xp[1] = iy + factor2*(Uinf[1] + _sample<2>(weights, start, end, iap, lab));
			}
	}
	
	template<typename LabVel>
	void remesh(LabVel& lab)
	{
		static const int KS = -1;
		static const int KE = +3;
		
		for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
				omega_new[iy][ix] = 0;
#ifndef NDEBUG
		for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
				pou[iy][ix] = 0;
#endif
		
		for(int iy=PSY; iy<PEY; iy++)
			for(int ix=PSX; ix<PEX; ix++)
			{
				const Real * const xp = xparticles[iy-PSY][ix-PSX];
				
				Real ap[2] = {
					floor(xp[0]),
					floor(xp[1]) };
				
				int iap[2] = {
					(int)ap[0],
					(int)ap[1] };
				
				Real weights[2][4];
				_computeWeights(xp, ap, weights);
				
				const int start[2] = {
					max(KS, 0 - iap[0]), 
					max(KS, 0 - iap[1])
				};
				
				const int end[2] = {
					min(KE, B::sizeX - iap[0]),  
					min(KE, B::sizeY - iap[1])
				};
				
				const Real omega = lab.template get<0>(ix, iy);
				
				
				_scatter(weights, start, end, iap, omega);
			}
		
#ifndef NDEBUG
		for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
				assert(_is_close(pou[iy][ix], 1.));
#endif
	}
};