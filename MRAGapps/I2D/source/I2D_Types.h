/*
 *  I2D_Types.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "I2D_Test.h"

using namespace std;
using namespace MRAG;

struct VelocityBlock
{
	Real u[2][_BLOCKSIZE_][_BLOCKSIZE_];
	
	void clear()
	{
		memset(u, 0, sizeof(Real)*_BLOCKSIZE_*_BLOCKSIZE_*2);
	}
	
	static VelocityBlock * allocate(const int nitems)
	{
		VelocityBlock * ptr = (VelocityBlock*)_mm_malloc(sizeof(VelocityBlock)*nitems, 16);
		
		return new (ptr) VelocityBlock[nitems];
	}
	
	static void deallocate(const VelocityBlock *& velblocks)
	{
		_mm_free(const_cast<VelocityBlock *>(velblocks));
		
		velblocks = NULL;
	}
	
	static void deallocate(VelocityBlock *& velblocks)
	{
		_mm_free(velblocks);
		
		velblocks = NULL;
	}
};

//===================================================
//		FLUID ELEMENT
//===================================================
struct FluidElement2D
{
	Real omega;
	Real u[2];
	Real tmp;
	
	operator Real() const 
	{
		abort();
		return (Real)omega;
	}
	
	void operator += (FluidElement2D f)
	{
		omega += f.omega;
		u[0] += f.u[0];
		u[1] += f.u[1];
		tmp += f.tmp;
	}
	
	Real giveMe(int i, Real h=0)
	{
		switch(i)
		{
			case 0: return omega;
			case 1: return u[0];
			case 2: return u[1];
			case 3: return tmp;
				
			default: abort(); return 0;	
		}
	}
	
	//scalar reconstruction in time
	Real time_rec(Real t) 
	{ 
		return omega + t*tmp;
	}
	
	/*void partA_ScalarRK2(Real t, Real dt)
	 {
	 psi = tmp;
	 omega -= tmp*t;
	 }
	 
	 void partB_ScalarRK2(Real t, Real dt)
	 {
	 omega = time_rec(t-dt/2.0) + tmp*(dt/2.0);
	 psi = 0.0;
	 }*/
};

inline FluidElement2D operator*(const FluidElement2D& p, Real v)
{
	FluidElement2D res;
	
	res.omega = p.omega*v;
	res.u[0] = p.u[0]*v;
	res.u[1] = p.u[1]*v;
	res.tmp = p.tmp*v;
	
	return res;
}


//===================================================
//		FLUID BLOCK
//===================================================
struct FluidBlock2D: public Block<FluidElement2D, _BLOCKSIZE_, _BLOCKSIZE_, 1>
{
	Real external_data[_BLOCKSIZE_][_BLOCKSIZE_];
	
	typedef FluidElement2D ElementType;
	
	FluidBlock2D(ElementType e = ElementType()): Block<FluidElement2D, _BLOCKSIZE_, _BLOCKSIZE_, 1>(e){}
	
	void serialize(FILE* outputBinary)
	{
		fwrite(reinterpret_cast<char*> (this),sizeof(FluidBlock2D),1,outputBinary);
	}
	
	void deserialize(FILE* inputBinary)
	{
		fread(reinterpret_cast<char*> (this),sizeof(FluidBlock2D),1,inputBinary);
	}
};


//===================================================
//		STREAMERSs
//===================================================
struct Streamer_OmegaLTS
{
	inline static Real operate(FluidElement2D& input, Real time) { return input.time_rec(time);}
};

struct Streamer_Omega
{
	inline static Real operate(FluidElement2D& input, Real time) { return input.omega;}
};

struct Streamer_Velocity
{
	template <int components>
	inline static void operate(FluidElement2D& input, Real time, Real output[components]){ abort(); }
};

template <> inline void Streamer_Velocity::operate<2>(FluidElement2D& input, Real time, Real output[2])
{
	output[0] = input.u[0];
	output[1] = input.u[1];
}

//===================================================
//		RK2
//===================================================
template <int stage>
struct UpdateScalarRK2
{
	Real t, dt;
	
	UpdateScalarRK2(Real t_, Real dt_): t(t_), dt(dt_) {}
	UpdateScalarRK2(const UpdateScalarRK2& task): t(task.t), dt(task.dt){}
	
	inline void operator() (const BlockInfo& info, FluidBlock2D& b) const
	{	
		static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
		FluidElement2D *const e = &b(0,0);
		const Real * const rhs = &b.external_data[0][0];
		
		switch(stage)
		{
			case 1:
				
				for(int i=0; i<n; i++) 
				{
					e[i].tmp = rhs[i];
					e[i].omega -= rhs[i]*t;
				}
				
				break;
			case 2:
				
				for(int i=0; i<n; i++) 
				{
					e[i].omega = e[i].time_rec(t-dt/2) + rhs[i]*(dt/2);
					e[i].tmp = 0;
				}
				
				break;
			default:
				abort();
		}
	}
};		

template <int stage>
struct UpdateScalarRK2_Simple
{
	Real dt;
	
	UpdateScalarRK2_Simple(Real dt_): dt(dt_) {}
	UpdateScalarRK2_Simple(const UpdateScalarRK2_Simple& c): dt(c.dt){}
	
	inline void operator() (const BlockInfo& info, FluidBlock2D& b) const
	{	
		static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
		FluidElement2D *const e = &b(0,0);
		const Real * const rhs = &b.external_data[0][0];
		
		switch(stage)
		{
			case 1:
				
				for(int i=0; i<n; i++) 
				{
					e[i].tmp = rhs[i];
					e[i].omega += rhs[i]*dt;
				}
				
				break;
			case 2:
				
				for(int i=0; i<n; i++) 
				{
					e[i].omega += e[i].tmp * (-dt/2.0) + rhs[i]*(dt/2.0);
					e[i].tmp = 0;
				}
				
				break;
			default:
				abort();
		}
	}
};		


//===================================================
//		PROJECTORs
//===================================================
template <typename T, int i> inline Real velocity_projector_impl(const T&t);
template <> inline Real velocity_projector_impl<FluidElement2D, 0>(const FluidElement2D&t) {return t.u[0];}
template <> inline Real velocity_projector_impl<FluidElement2D, 1>(const FluidElement2D&t) {return t.u[1];}
make_projector(velocity_projector, velocity_projector_impl)

template <typename T, int i> inline Real vorticity_projector_impl(const T&t);
template <> inline Real vorticity_projector_impl<FluidElement2D, 0>(const FluidElement2D&t) {return t.omega;}
make_projector(vorticity_projector, vorticity_projector_impl)

template <typename T, int i> inline Real vorticityANDvelocityANDchi_projector_impl(const T&t);
template <> inline Real vorticityANDvelocityANDchi_projector_impl<FluidElement2D, 0>(const FluidElement2D&t) {return t.omega;}
template <> inline Real vorticityANDvelocityANDchi_projector_impl<FluidElement2D, 1>(const FluidElement2D&t) {return t.u[0];}
template <> inline Real vorticityANDvelocityANDchi_projector_impl<FluidElement2D, 2>(const FluidElement2D&t) {return t.u[1];}
template <> inline Real vorticityANDvelocityANDchi_projector_impl<FluidElement2D, 3>(const FluidElement2D&t) {return t.tmp;}
make_projector(vorticityANDvelocityANDchi_projector, vorticityANDvelocityANDchi_projector_impl)

template <typename T, int i> inline Real vorticityANDvelocity_projector_impl(const T&t);
template <> inline Real vorticityANDvelocity_projector_impl<FluidElement2D, 0>(const FluidElement2D&t) {return t.omega;}
template <> inline Real vorticityANDvelocity_projector_impl<FluidElement2D, 1>(const FluidElement2D&t) {return t.u[0];}
template <> inline Real vorticityANDvelocity_projector_impl<FluidElement2D, 2>(const FluidElement2D&t) {return t.u[1];}
make_projector(vorticityANDvelocity_projector, vorticityANDvelocity_projector_impl)

template <typename T, int i> inline Real obstacle_projector_impl(const T&t);
template <> inline Real obstacle_projector_impl<FluidElement2D, 0>(const FluidElement2D&t) {return t.tmp;}
make_projector(obstacle_projector, obstacle_projector_impl)


//===================================================
//		COMPUTATION MAIN CLASS
//===================================================

//typedef Wavelets_Haar W;
//typedef Wavelets_Interp4thOrder W;
//typedef Wavelets_Interp2ndOrder W;
//typedef Wavelets_AverageInterp3rdOrder W;
typedef Wavelets_AverageInterp5thOrder W;
typedef FluidBlock2D B;
typedef Multithreading::BlockProcessing_TBB<B> BlockProcessing;

