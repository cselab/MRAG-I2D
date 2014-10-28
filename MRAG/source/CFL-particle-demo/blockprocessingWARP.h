/*
 *  blockprocessingWARP.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/28/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include <limits>
using namespace std;

class MP4
{
public:
	
	inline static const Real _f0(Real t)
	{
#ifndef NDEBUG
		assert(t>=0);
#endif
		return 1.+t*t*(-5./2.+3./2.*t);
	}
	
	inline static const Real _f1(Real t)
	{
#ifndef NDEBUG
		assert(t>=0);
#endif
		return 2.+t*(-4. + t*(5./2.-1./2.*t));
	}
	
	static const Real eval(Real x)
	{
		Real t = fabs(x);
		switch(min((int)2, (int)t))
		{
			case 2:
				return 0;
			case 1:
				return _f1(t);
			case 0:
				return _f0(t);
		}
		
		abort();
		return 0;
	}
};

#pragma mark -
#pragma mark CFLParticle classes
template<int blockSize, int maximumDistance=1>
struct CLFParticle
{
	static const int startParticle = -2;
	static const int endParticle = blockSize+2;
	static const int nParticles = endParticle - startParticle;

	static const int startVelocity = -3;
	static const int endVelocity = blockSize+3;
	static const int nVelocityPoints = endVelocity-startVelocity;

#pragma mark FindMaximumDT
	struct FindMaximumDT
	{
		FindMaximumDT(){}
		FindMaximumDT(const FindMaximumDT& c){}

		template<typename B>
		inline void operator()(const BlockInfo& info, B& b) const
		{
			typedef typename B::ElementType E;
			
			const int n = B::sizeZ*B::sizeY*B::sizeX;
			
			Real dt = HUGE_VALF;
			const Real h = info.h[0]; 
			E* ptrE = &(b[0]);
			
			for(int iE=0; iE<n; iE++, ptrE++)
			{
				const Real maxVel = max(fabs(ptrE->u[0]), fabs(ptrE->u[1]));
				dt = min(dt, h*maximumDistance/maxVel);
			}

			b.maximumDT = dt;
		}
	};
	
#pragma mark UpdateOmega
	struct UpdateOmega
	{
		UpdateOmega(){}
		UpdateOmega(const UpdateOmega& c){}
		
		template<typename B>
		inline void operator()(const BlockInfo& info, B& b) const
		{
			typedef typename B::ElementType E;
			
			const int n = B::sizeZ*B::sizeY*B::sizeX;
			
			E* ptrE = &(b[0]);
			
			for(int iE=0; iE<n; iE++, ptrE++)
				ptrE->w = ptrE->tmp;
		}
	};


#pragma mark ParticleSampler
	struct ParticleSampler
	{
		Real wx[4], wy[4];
		int iap[2];
		
		void _computeWeights(const Real xp[2]) 
		{			
			const Real ap[2] = { floor(xp[0]), floor(xp[1]) };
			
			iap[0] = (int)(xp[0]>=ap[0]+0.5) + (int)ap[0];
			iap[1] = (int)(xp[1]>=ap[1]+0.5) + (int)ap[1];
			
			Real tx[4] = {
				fabs(xp[0] - iap[0] - (0.5+-2)),
				fabs(xp[0] - iap[0] - (0.5+-1)),
				fabs(xp[0] - iap[0] - (0.5+ 0)),
				fabs(xp[0] - iap[0] - (0.5+ 1))
			};
			
			wx[0] = MP4::_f1(tx[0]);
			wx[1] = MP4::_f0(tx[1]);
			wx[2] = MP4::_f0(tx[2]);
			wx[3] = MP4::_f1(tx[3]);
			
			Real ty[4] = {
				fabs(xp[1] - iap[1] - (0.5+-2)),
				fabs(xp[1] - iap[1] - (0.5+-1)),
				fabs(xp[1] - iap[1] - (0.5+ 0)),
				fabs(xp[1] - iap[1] - (0.5+ 1))
			};
			
			wy[0] = MP4::_f1(ty[0]);
			wy[1] = MP4::_f0(ty[1]);
			wy[2] = MP4::_f0(ty[2]);
			wy[3] = MP4::_f1(ty[3]);
		}
		
		template<int startLayer, int endLayer, int components>
		void m2p(Real x[2], Layer<startLayer, endLayer, components>& workingplace, Real sample[components])
		{
			static const int kernel_start = -2;
			static const int kernel_end = 2;
			
			_computeWeights(x);
			
			for(int ic=0; ic<components; ic++)
				sample[ic] = 0;
			
			const int start[2] = {max(kernel_start, startLayer-iap[0]), max(kernel_start, startLayer-iap[1]) };
			const int end[2] = {min(kernel_end, endLayer-iap[0]),  min(kernel_end, endLayer-iap[1])};

			Real unit_partition = 0;
			
			for(int ic=0; ic<components; ic++)
			for(int sy=start[1]; sy<end[1]; sy++)
			for(int sx=start[0]; sx<end[0]; sx++)
			{
				unit_partition += wx[sx-kernel_start]*wy[sy-kernel_start];
				sample[ic] += wx[sx-kernel_start]*wy[sy-kernel_start] * workingplace(iap[0] + sx, iap[1] + sy, ic);
			}
			
			if (!(unit_partition > 1. - 1e-6))
				printf("not ! %f\n", 1.-unit_partition);
			assert(unit_partition > 1. - 1e-6);
		}
		
		template<int startLayer, int endLayer, int components>
		void p2m(Real value[components], Real x[2], Layer<startLayer, endLayer, components>& workingplace)
		{
			static const int kernel_start = -2;
			static const int kernel_end = 2;
		
			_computeWeights(x);
			
			const int start[2] = {max(kernel_start, startLayer-iap[0]), max(kernel_start, startLayer-iap[1]) };
			const int end[2] = {min(kernel_end, endLayer-iap[0]),  min(kernel_end, endLayer-iap[1])};
						
			for(int ic=0; ic<components; ic++)
			for(int sy=start[1]; sy<end[1]; sy++)
			for(int sx=start[0]; sx<end[0]; sx++)
				workingplace(iap[0] + sx, iap[1] + sy, ic) += value[ic]*wx[sx-kernel_start]*wy[sy-kernel_start];
		}
	};
	

#pragma mark WorkingPlaceWARP

	struct WorkingPlaceWARP
	{
		Layer<startVelocity, endVelocity, 2> u;
		Layer<startParticle, endParticle, 2> x; 
		Layer<startParticle, endParticle> w_carried;
		Layer<startParticle, endParticle, 2> tmp;
		Layer<0, blockSize, 2> w;
		
		ParticleSampler sampler;
	};

#pragma mark BlockLabWARP
	template<typename BlockType>
	class BlockLabWARP: public BlockLab<BlockType>
	{		
		WorkingPlaceWARP* m_workingplace;

	public:
		
		BlockLabWARP(): BlockLab<BlockType>()
		{
			m_workingplace = allocator<WorkingPlaceWARP >().allocate(1);
		}
		
		~BlockLabWARP()
		{
			allocator<WorkingPlaceWARP >().deallocate(m_workingplace, 1);
			m_workingplace = NULL;
		}
		
		WorkingPlaceWARP& workingPlace(){return *m_workingplace;}
		
		const int * getStencilStart(){return this->m_stencilStart;}
	};

#pragma mark AdvectParticles
	struct AdvectParticles
	{
		Real dt;
		int stencil_start[3], stencil_end[3];
		
		AdvectParticles(float dt_): dt(dt_) 
		{
			stencil_start[0] = stencil_start[1] = startVelocity;
			stencil_start[2] = 0;
			
			stencil_end[0] = stencil_end[1] = endVelocity - blockSize + 1;
			stencil_end[2] = +1; 
		}
		
		AdvectParticles(const AdvectParticles & copy) : dt(copy.dt)
		{
			stencil_start[0] = stencil_start[1] = startVelocity;
			stencil_start[2] = 0;
			
			stencil_end[0] = stencil_end[1] = endVelocity - blockSize + 1;
			stencil_end[2] = +1;    
		}
		
		template<typename Lab, typename WP>
		void prepare(const BlockInfo& info, Lab& lab, WP& wp) const
		{
			const Real factor = 1/info.h[0];
			
			for(int iy=startVelocity; iy<endVelocity; iy++)
			for(int ix=startVelocity; ix<endVelocity; ix++)
			{
				wp.u(ix, iy, 0) = lab(ix, iy).u[0]*factor;
				wp.u(ix, iy, 1) = lab(ix, iy).u[1]*factor;
			}
			
			for(int iy=startParticle; iy<endParticle; iy++)
			for(int ix=startParticle; ix<endParticle; ix++)
			{
				wp.x(ix, iy,0) = ix + 0.5;
				wp.x(ix, iy,1) = iy + 0.5;
				wp.w_carried(ix, iy) = lab(ix,iy).w;
			}
			
			for(int iy=0; iy<blockSize; iy++)
			for(int ix=0; ix<blockSize; ix++)
			{
				wp.w(ix,iy, 0) = 0;
				wp.w(ix,iy, 1) = 0;
			}
		}
		
		template<typename WP>
		void advectGridPoints(WP& wp) const
		{
			const Real factor = dt*0.5;
			
			for(int iy=startParticle; iy<endParticle; iy++)
			for(int ix=startParticle; ix<endParticle; ix++)
			{
				Real x_mesh[2] = {ix + 0.5, iy + 0.5};
				
				for(int ic=0; ic<2; ic++)
					wp.x(ix, iy, ic) = x_mesh[ic] + wp.u(ix, iy, ic)*factor;
			}
		}
		
		template<typename WP>
		void advectParticles(WP& wp) const
		{
			for(int iy=startParticle; iy<endParticle; iy++)
			for(int ix=startParticle; ix<endParticle; ix++)
			{
				Real x_current[2] = {
					wp.x(ix, iy, 0),
					wp.x(ix, iy, 1)
				};
				
				Real new_velocity[2];
				
				wp.sampler.m2p(x_current, wp.u, new_velocity);
				
				wp.tmp(ix, iy, 0) = new_velocity[0];
				wp.tmp(ix, iy, 1) = new_velocity[1];
				wp.x(ix, iy, 0) = ix + 0.5 + new_velocity[0]*dt;
				wp.x(ix, iy, 1) = iy + 0.5 + new_velocity[1]*dt;
			}
		}
		
		template<typename WP, typename Block>
		void remeshParticles(Block& block, WP& wp) const 
		{
			for(int iy=startParticle; iy<endParticle; iy++)
			for(int ix=startParticle; ix<endParticle; ix++)
			{
				Real x_current[2] = {
					wp.x(ix, iy, 0),
					wp.x(ix, iy, 1)
				};

				Real value[2] = {wp.w_carried(ix, iy), 1} ;
				wp.sampler.p2m(value, x_current, wp.w);
			}
		}
		
		template<typename WP, typename Block>
		void finalize(Block& block, WP& wp) const
		{
			for(int iy=0; iy<Block::sizeY; iy++)
			for(int ix=0; ix<Block::sizeX; ix++)
			{
				block(ix, iy).tmp = wp.w(ix, iy,0);
				
				//if (wp.w(ix,iy,1) < 1 - 1e-5)
				//	printf("partition of unit corrupted %d %d %e\n", ix, iy, wp.w(ix,iy,1));
				//assert(wp.w(ix,iy,1) > 1 - 1e-5);
			}
		}
		
		template<typename BlockType>
		inline void operator()(BlockLabWARP<BlockType>& lab, const BlockInfo& info, BlockType& o) const 
		{
			prepare(info, lab, lab.workingPlace());
			advectGridPoints(lab.workingPlace());
			advectParticles(lab.workingPlace());
			remeshParticles(o,lab.workingPlace());
			finalize(o,lab.workingPlace());
		}
	};

	struct AdvectParticlesRK3: public AdvectParticles
	{
		AdvectParticlesRK3(float dt_): AdvectParticles(dt_) { }
		AdvectParticlesRK3(const AdvectParticlesRK3 & copy) : AdvectParticles(copy.dt){  }
		
		template<typename WP>
		void thirdStep(WP& wp) const
		{
			const Real dt_ = this->dt;
			//1. copy velocity 2 to tmp
			//2. construct position 3: x_mesh - u_mesh + 3*tmp
			//3. sample velocity, store it in u_p
			//4. assemble final position: x_mesh + 1/6(u_mesh + 4*tmp + u_p)
			
			//1. already done in advectParticles
			
			//2. - 3. - 4.
			for(int iy=startParticle; iy<endParticle; iy++)
			for(int ix=startParticle; ix<endParticle; ix++)
			{
				Real x_current[2] = {
					ix + 0.5 + dt_*(2*wp.tmp(ix,iy,0)-wp.u(ix,iy,0)), 
					iy + 0.5 + dt_*(2*wp.tmp(ix,iy,1)-wp.u(ix,iy,1))
				};
				
				Real new_velocity[2];
				
				wp.sampler.m2p(x_current, wp.u, new_velocity);
				
				wp.x(ix, iy, 0) = ix + 0.5 + 1./6*(wp.u(ix,iy,0) + 4*wp.tmp(ix,iy,0) + new_velocity[0])*dt_;
				wp.x(ix, iy, 1) = iy + 0.5 + 1./6*(wp.u(ix,iy,1) + 4*wp.tmp(ix,iy,1) + new_velocity[1])*dt_;
			}
		}

		
		template<typename BlockType>
		inline void operator()(BlockLabWARP<BlockType>& lab, const BlockInfo& info, BlockType& o) const 
		{
			this->prepare(info, lab, lab.workingPlace());
			this->advectGridPoints(lab.workingPlace());
			this->advectParticles(lab.workingPlace());
			thirdStep(lab.workingPlace());
			this->remeshParticles(o,lab.workingPlace());
			this->finalize(o,lab.workingPlace());
		}
	};
};

