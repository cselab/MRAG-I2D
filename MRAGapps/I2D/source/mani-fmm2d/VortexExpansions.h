/*
 *  VortexExpansions.h
 *  hcfmm
 *
 *  Created by Manfred on 1/29/09.
 *  Copyright 2009 ETHZ. All rights reserved.
 *
 */
#pragma once

#include <complex>
//#include "hcfmm_complex.h"
#include "hcfmm_types.h"
#include "binomial.h"

//size_t global_box;

//------------------------------VortexExpansions--------------------------//
//This is the physics that is required to calculate the Expansions.        //
//require 3 Functions: 
//void calculateExpansions (tParticle* in_particles, tExpansions* out_expansions)
//void gatherExpansions (tExpansions* in_expansions, tExpansions* out_expansions)
//void evaluateExpansions (tExpansions* in_expansions, tRHS* out_RHS)

template <typename Particle, int _order>
struct _VortexExpansions {
	typedef _VortexExpansions<Particle,_order> SelfType;
	typedef Particle ParticleType;
	typedef typename Particle::RHSType tRHS;
	typedef typename Particle::BaseType Btype;
	typedef	std::complex<Btype> ExpansionsValueType;
	//typedef HCFMM::fmmcomplex<Btype> ExpansionsValueType;
	static const int order=_order;
	
	//The Data:
	ExpansionsValueType values[Particle::pdim][_order+1];
	Btype Center[Particle::dim];
	Btype Radius; //New: using the radius instead of h.
	
	void calculateExpansions (ParticleType* in_particles, int nParticles);
	void gatherExpansions (SelfType* in_expansions);
    void evaluateExpansions (Btype *location,tRHS* out_RHS);
	
	

	
	void clear()
	{
		for (int pd=0;pd<Particle::pdim;++pd) //HOLLY CRAP: when you put _order here instead of Particle::pdim, you can look for a long time (bless valgrind)
		{
			for (int n=0;n<_order;++n)
			{
				
				ExpansionsValueType dummy=ExpansionsValueType(0,0);
				values[pd][n]=dummy;
			}
		}

		
	}

	_VortexExpansions():Radius(0)
	{
		clear();
		
		for (int d=0; d<Particle::dim;++d)
		{
			Center[d]=0;
		}
	}
	
	_VortexExpansions(const _VortexExpansions &B)
	{
		
		for (int pd=0;pd<Particle::pdim;++pd)
		{
			for (int n=0;n<_order;++n)
			{
				
					values[pd][n]=B.values[pd][n];
			
			}
		}
		
	
		for (int d=0; d<Particle::dim;++d)
		{
			Center[d]=B.Center[d];
		}
		
	}
	
	void setcenter(Btype* in_center)
	{
		for (int d=0; d<Particle::dim;++d)
		{
			Center[d]=in_center[d];
		}
	}
	
	void print()
	{
		for (int n=0;n<_order;++n)
		{
				for (int d=0;d<Particle::pdim;++d)
				{
					//std::cout << "expansion [" << d << "][" << n << "]="<<values[d][n] <<std::endl;
					std::cout << "expansion [" << d << "][" << n << "]="<<values[d][n].real() << " " << values[d][n].imag()  <<std::endl;
					if(isnan(values[1][1].real()))
					{
						std::cout << "warning: we have nans" << std::endl;
					}
				}
		}
	}
	
	
};


#if 0 // this is an implementation based on Rick Beatson's and Leslie Gerards "Short course on Fast Multipole Methods" 
//needed some adaptation for the usage for the Vortex case instead of gravity case. but should  only be for the evaulation actually.

template <typename Particle, int _order>
void _VortexExpansions<Particle,_order>::calculateExpansions(Particle* in_particles, int nParticles)
{
	//ExpansionsValueType Cnm[Particle::pdim][_order*(_order+2)];
	ExpansionsValueType rp,prod,ZM; //Complex Distance relative to center: (Z_m-Z_M)
	
	Particle* p=in_particles;
	
	clear(); //initialize to zero.
	
	ZM=ExpansionsValueType(this->Center[0],this->Center[1]);
	
	for (int i=0;i<nParticles;++i)
	{
		
		//calculate distance:
		rp=ExpansionsValueType(p[i].x[0],p[i].x[1])-ZM;
		
		for (int pd=0;pd<Particle::pdim;++pd)
		{
		this->values[pd][0]+=p[i].w[pd];
		//prod=ExpansionsValueType(1,0);  //power of rp.
		for (int n=1;n<_order;++n)
		{
			//TODO: can improve this in efficiency
			
				//this->values[pd][n]+=p[i].w[pd]*prod;
				this->values[pd][n]-=p[i].w[pd]/Btype(n)*pow(rp,n);
				
			
			//prod=prod*rp;
		}
		}
	} //end particle-loop.
	
	ExpansionsValueType prefac=ExpansionsValueType(1/(2*M_PI),1);
	for (int pd=0;pd<Particle::pdim;++pd)
		for (int n=1;n<_order;++n)
			this->values[pd][n]*=prefac;
	
	
	
	
	
}


/*
 * Shift the Coefficients of a child to the location of the parent and adds it to the parent.
 **/

template <typename Particle, int _order>
void _VortexExpansions<Particle,_order>::gatherExpansions (_VortexExpansions<Particle,_order>* in_expansions)
{
	std::cout << "SHOULD NOT BE HERE" <<std::endl;
	unsigned long int bcoeff(0);
	ExpansionsValueType rb,prod,csum;
	rb=ExpansionsValueType(in_expansions->Center[0],in_expansions->Center[1])-ExpansionsValueType(this->Center[0],this->Center[1]);
	
	
	for (int j=0; j<_order;++j)//looping to get all the parents expansions filled up.
	{
		//for each value a_k wee need to loop from 0..k, and sum up. 
		for (int pd=0;pd<Particle::pdim;++pd)
		{
			csum=ExpansionsValueType(0,0);
			prod=ExpansionsValueType(1,0);
			for (int k=j;k>=0;--k) //looping reverse because of term rp^(l-k)
			{
				bcoeff=binomial(j,k);
				csum+=in_expansions->values[pd][k]*prod*ExpansionsValueType(bcoeff,0);
				
				prod*=rb; //rb^(j-k)
			}
			this->values[pd][j]+=csum;
		}
		
	}
	
	
	
}

template <typename Particle, int _order>
void _VortexExpansions<Particle,_order>::evaluateExpansions (typename Particle::BaseType *location, typename Particle::RHSType *out_RHS)
{
	
    ExpansionsValueType csum,rp,prefac,prod;
	
	rp=ExpansionsValueType(location[0],location[1])-ExpansionsValueType(this->Center[0],this->Center[1]);
	
	prefac=ExpansionsValueType(1.0/(2.0*M_PI),1); //factor before sum: i/2Pi
	
	//calculation of sum: Sum_{k=0}^{P}(a_k/(z-zm)^k)
	for (int pd=0;pd<Particle::pdim;++pd)
	{
		
		prod=ExpansionsValueType(1,0)/rp;
		csum=ExpansionsValueType(std::real(this->values[pd][0]),0)/rp*prefac;
		for (int n=1;n<_order;++n)
		{
			prod/=rp; //1/(rp^n)
			csum-=prod*this->values[pd][n]*Btype(n)*prod;
			
		}
		//csum*=prefac;
		//TODO: ok this really makes no sense for pd>0;
		out_RHS->x[0]+=std::real(csum);
		out_RHS->x[1]+=-std::imag(csum);
		
	}
	
	
	
}


#else  //This is an implementation based on petros Thesis

template <typename Particle, int _order>
void _VortexExpansions<Particle,_order>::calculateExpansions(Particle* in_particles, int nParticles)
{
	//ExpansionsValueType Cnm[Particle::pdim][_order*(_order+2)];
	ExpansionsValueType rp,prod,ZM; //Complex Distance relative to center: (Z_m-Z_M)
	
	Particle* p=in_particles;
	
	clear(); //initialize to zero.
	
	ZM=ExpansionsValueType(this->Center[0],this->Center[1]);
	for (int i=0;i<nParticles;++i)
	{
		
		//calculate distance:
		rp=ExpansionsValueType(p[i].x[0],p[i].x[1])-ZM;
		prod=ExpansionsValueType(1,0);  //power of rp.
		for (int n=0;n<_order;++n)
		{
			
			for (int pd=0;pd<Particle::pdim;++pd)
			{
				this->values[pd][n]+=p[i].w[pd]*prod;
				//this->values[pd][n]+=p[i].w[pd]*pow(rp,n);
				
			}
			prod=prod*rp;
		}
		
	} //end particle-loop.
	

}


/*
 * Shift the Coefficients of a child to the location of the parent and adds it to the parent.
 **/

template <typename Particle, int _order>
void _VortexExpansions<Particle,_order>::gatherExpansions (_VortexExpansions<Particle,_order>* in_expansions)
{
	//std::cout << "SHOULD NOT BE HERE" <<std::endl;
	unsigned long int bcoeff(0);
	ExpansionsValueType rb,prod,csum;
	rb=ExpansionsValueType(in_expansions->Center[0],in_expansions->Center[1])-ExpansionsValueType(this->Center[0],this->Center[1]);
	
	
	for (int j=0; j<_order;++j)//looping to get all the parents expansions filled up.
	{
		//for each value a_k wee need to loop from 0..k, and sum up. 
		for (int pd=0;pd<Particle::pdim;++pd)
		{
			csum=ExpansionsValueType(0,0);
			prod=ExpansionsValueType(1,0);
			for (int k=j;k>=0;--k) //looping reverse because of term rp^(l-k)
			{
				bcoeff=binomial(j,k);
				csum+=in_expansions->values[pd][k]*prod*ExpansionsValueType(bcoeff,0);
		
				prod*=rb; //rb^(j-k)
			}
			this->values[pd][j]+=csum;
		}
		
	}
	

	
}

template <typename Particle, int _order>
void _VortexExpansions<Particle,_order>::evaluateExpansions (typename Particle::BaseType *location, typename Particle::RHSType *out_RHS)
{
    ExpansionsValueType csum,rp,prefac,prod;
	
	rp=ExpansionsValueType(location[0],location[1])-ExpansionsValueType(this->Center[0],this->Center[1]);
	
	prefac=-ExpansionsValueType(0,1.0/(2.0*M_PI)); 
	prefac=prefac/rp;                          //factor before sum: i/2Pi*1/(z-zm)
	
	//calculation of sum: Sum_{k=0}^{P}(a_k/(z-zm)^k)
	for (int pd=0;pd<Particle::pdim;++pd)
	{
		csum=ExpansionsValueType(0,0);
		prod=ExpansionsValueType(1,0);
#pragma unroll
	for (int n=0;n<_order;++n)
	{
		csum+=prod*this->values[pd][n];
		prod/=rp; //1/(rp^n)
		//csum+=this->values[pd][n]*ExpansionsValueType(1,0)/pow(rp,n);
	}
		csum*=prefac;
		//TODO: ok this really makes no sense for pd>0;
		//out_RHS->x[0]+=std::real(csum);
		//out_RHS->x[1]-=std::imag(csum);
		out_RHS->x[0]+=csum.real();
		out_RHS->x[1]-=csum.imag();
	}
	
	
	//global_box++;


}
			
			
#endif

			
			


