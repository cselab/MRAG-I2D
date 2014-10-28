/*
 *  PotentialExpansions.h
 *  hcfmm
 *
 *  Created by Manfred on 1/29/09.
 *  Copyright 2009 ETHZ. All rights reserved.
 *
 */

#include <fstream>
#include <complex>
#include "hcfmm_types.h"






//maping function for the negative indices
inline int _map(int n, int m)
{
	int res=n*n+(m+n);
	assert(res<_ORDER_*(_ORDER_+2));
	return res;
}


//debug: printing to file:
template <typename Ntype>
void dumpData(Ntype* data, char* fname )
{
	std::ofstream outf(fname);
	outf.setf(std::ios_base::scientific);
	outf.width(15);
	for (int n=0;n<_ORDER_;++n)
	{
		for (int m=-n;m<=n;++m)
		{
			outf  << data[_map(n,m)]<<std::endl;
		}
	}
	
	outf.close();
	
}

//NEW implementation based on Kautzleben 1965 (Aeronomie, BandI) (digged out from Bosch 2000)
template <class Particle, int _order>
void calcPnm(typename Particle::BaseType Pnm_out[_order][_order], const typename Particle::BaseType  &theta_in)
{
	typedef typename Particle::BaseType Btype;
	Btype sine,cosine;
	sine=sin(theta_in);
	cosine=cos(theta_in);
	
	
	//first calculate without (-1)^m
	Pnm_out[0][0]=Btype(1);
	Pnm_out[1][0]=cosine;
	Pnm_out[1][1]=sine;
	
	for (int n=2;n<_order;++n)
	{
		Pnm_out[n][0]=(2*n-1)/Btype(n)*cosine*Pnm_out[n-1][0]-(n-1)/Btype(n)*Pnm_out[n-2][0];
		Pnm_out[n][n-1]=(2*n-1)*sine*Pnm_out[n-1][n-2];
		
		for (int m=1;m<n-1;++m)
		{
			Pnm_out[n][m]=(2*n-1)*sine*Pnm_out[n-1][m-1]+Pnm_out[n-2][m];
		}
		Pnm_out[n][n]=(2*n-1)*sine*Pnm_out[n-1][n-1];
		
		
	}
	
	
	
	
	for (int n=1;n<_order;++n)
	{
		for (int m=1;m<=n;m+=2)
		{
			Pnm_out[n][m]*=(-1);  //add factor (-1)^m (only for odd m)
		}
	}
	
	
}


template<typename NumType, int _order>
void calcYnm(NumType phi, NumType Pnm[_order][_order], std::complex<NumType> Ynm[_order*(_order+2)] )
{
	for (int n=0;n<_order;++n)
	{
		Ynm[_map(n,0)]=Pnm[n][0];  //special case n=0: sqrtfac=1, cos(0)=1, sin(0)=0
		
		for (int m=1;m<=n;++m)
		{
			//Ynm[_map(n,m)]=sqrtfac<NumType>(n,m)*Pnm[n][m]*exp(m*phi*std::complex<NumType>(0,1));
			Ynm[_map(n,m)]=sqrtfac<NumType,_order>(n,m)*Pnm[n][m]*exp(m*phi*std::complex<NumType>(0,1));

			Ynm[_map(n,-m)]=std::conj(Ynm[_map(n,m)]);

		}
	}
}


//------------------------------ParticleExpansions--------------------------//
//This is the physics that is required to calculate the Expansions.        //
//require 3 Functions: 
//void calculateExpansions (tParticle* in_particles, tExpansions* out_expansions)
//void gatherExpansions (tExpansions* in_expansions, tExpansions* out_expansions)
//void evaluateExpansions (tExpansions* in_expansions, tRHS* out_RHS)

template <typename Particle, int _order>
struct _PotentialExpansions {
	typedef _PotentialExpansions<Particle,_order> SelfType;
	typedef Particle ParticleType;
	typedef typename Particle::RHSType tRHS;
	typedef typename Particle::BaseType Btype;
	typedef	std::complex<Btype> ExpansionsValueType;
	static const int order=_order;
	
	//The Data:
	ExpansionsValueType values[Particle::pdim][_order*(_order+2)];
	Btype Center[Particle::dim];
	Btype Radius; //New: using the radius instead of h.
	
	void calculateExpansions (ParticleType* in_particles, int nParticles);
	void gatherExpansions (SelfType* in_expansions);
    void evaluateExpansions (Btype *location,tRHS* out_RHS);
	
	

	
	void clear()
	{
		for (int pd=0;pd<Particle::pdim;++pd)
			for (int n=0;n<_order;++n)
				for (int m=-n;m<=n;++m)
					values[pd][_map(n,m)]=0;
	}

	
	
	_PotentialExpansions():Radius(0)
	{
		clear();
		
		for (int d=0; d<Particle::dim;++d)
			Center[d]=0;
	}
	
	_PotentialExpansions(const _PotentialExpansions &B)
	{
		for (int pd=0;pd<Particle::pdim;++pd)
			for (int n=0;n<_order;++n)
				for (int m=-n;m<=n;++m)
					values[pd][_map(n,m)]=B.values[pd][_map(n,m)];

		for (int d=0; d<Particle::dim;++d)
			Center[d]=B.Center[d];
	}
	
	void setcenter(Btype* in_center)
	{
		for (int d=0; d<Particle::dim;++d)
			Center[d]=in_center[d];
	}
	
	void print()
	{
		for (int n=0;n<_order;++n)
			for (int m=-n;m<=n;++m)
				for (int d=0;d<Particle::pdim;++d)
					std::cout << "expansion [" << d << "][" << n << "][" <<m << "]="<<values[d][_map(n,m)] <<std::endl; 
	}
};

template <typename Particle, int _order>
void _PotentialExpansions<Particle,_order>::calculateExpansions(Particle* in_particles, int nParticles)
{
	ExpansionsValueType Ynm[_order*(_order+2)]; 
	//ExpansionsValueType Cnm[Particle::pdim][_order*(_order+2)];
	Btype Pnm[_order][_order];
	Btype rp[Particle::dim]; //spherical coordinates. relative to center.
	Particle* p=in_particles;
	Btype cosine,sine,val,prod;
	
	clear(); //initialize to zero.
	
	//debug:
//	std::cout << "Center of expansion is: " << std::endl;
//	std::cout << "[ " << this->Center[0] << "," <<this->Center[1] << "," <<this->Center[2] << std::endl;  
	
	for (int i=0;i<nParticles;++i)
	{
		cart2sph(p[i].x,this->Center, rp); //conversion to spherical coords.
		cosine=cos(rp[1]);
		sine=sin(rp[1]);
		val=-sine;
		prod=Btype(1);
		
		//calc: Pnm
		calcPnm<Particle,_order>(Pnm,rp[1]);
		//calc Ynm
		calcYnm(rp[2], Pnm,Ynm);
        //debug:
		//std::cout << "Ynm("<<rp[1]<<","<<rp[2]<<")"<< std::endl;
//		for (int n=0;n<_order;++n)
//			for (int m=-n;m<=n;++m)
//			{
//				std::cout << n << " " << m  << " " << Ynm[_map(n,m)] <<std::endl;
//			}
		//calc Cnm (expansions) //there seems to be a difference in Epton/Dembart to Greengard. the whole term i^-m*Anm/(-1)^n is missing in Eptons formulation for Cnm.
		prod=Btype(1);
		val=rp[0];
		for (int n=0;n<_order;++n)
		{
			for (int m=-n;m<=n;++m) // doing -n,n to follow Matlab Matlab implementation (might undo for performance, later)
			{
				for (int d=0;d<Particle::pdim;++d)
				{
					Btype A=Anm(n,m);
					int pm1=pm1n(n);
					//this->values[d][_map(n,m)]+=(p[i].w[d]*prod*Ynm[_map(n,m)]*A)/Btype(pm1);
					//Epton?:
					//this->values[d][_map(n,m)]+=(p[i].w[d]*prod*Ynm[_map(n,m)]);
					//afterMatlab:
					this->values[d][_map(n,m)]+=(p[i].w[d]*prod*Ynm[_map(n,-m)]);

				}
			}
			prod=prod*val;
			//debug:
			//std::cout << "debugging expansions:  i=" << i << ", n=" << n << ", prod(rp^n)= " <<prod << " w[i]" << p[i].w[0] << std::endl;
		}
		
		
	} //end particle-loop.
	
}

//Another Implementation, after reading deeply into the papers, doing a Matlab-prototype ... Wed. 1.April.2009
template <typename Particle, int _order>
void _PotentialExpansions<Particle,_order>::gatherExpansions (_PotentialExpansions<Particle,_order>* in_expansions)
{
	Btype Pnm[_order][_order];
	ExpansionsValueType Ynm[_order*(_order+2)]; 
	std::complex<Btype> csum;
	Btype rb[Particle::dim]; //spherical coordinates of childbox. relative to center.
	Btype cosine,sine,val,prod,rinv;
	//some temporary variables:
	Btype tmp0(0),tmp1(0);
	int tmpexp=0;
	
	ExpansionsValueType CI=ExpansionsValueType(0,1);
	
	//calculate spherical coordinates of child
	cart2sph(in_expansions->Center, this->Center,rb);
	//TODO: Should we compute the radius of the Box again?
	//compute Polynomial:
	rinv=Btype(1)/rb[0];
	sine=sin(rb[1]);
	cosine=cos(rb[1]);
	val=-sine;
	prod=Btype(1);
	//compute Pnm: 
	calcPnm<Particle,_order>(Pnm, rb[1]);
	//compute Ynm 
	calcYnm(rb[2], Pnm,Ynm);
	
	
	//compute sum of expansions. (using Greengard,Rohklin 1997, note that a naive implementation of their formula would lead to negative factorials!)
	//outer loops for Mjk
	val=rb[0];
	for (int j=0; j<_order;++j)
	{
		for (int k=-j;k<=j;++k)
		{
			for (int pd=0;pd<Particle::pdim;++pd)
			{
			// now the inner loops for n,-m (looping over the childrens part)
			csum=ExpansionsValueType(0,0);
			prod=1;	
			for(int n=0; n<=j;++n)
			{
				for(int m=-n;m<=n;++m)
				{
					int pr=j-n;
					int qr=k-m;
					if(abs(pr)>=abs(qr))
					{
						//tmp0=Anm(pr,qr);
						//tmp1=Anm(n,m)*tmp0/Anm(j,k);
						tmp1=calcFracAnm<Btype,_order>(j,k,n,m);
						tmpexp=abs(k)-abs(m)-abs(k-m);
						//csum+=in_expansions->values[pd][_map(pr,qr)]*pow(CI,tmpexp)*pow(rb[0],n)*Ynm[_map(n,-m)]*tmp1;
						csum+=in_expansions->values[pd][_map(pr,qr)]*pow(CI,tmpexp)*prod*Ynm[_map(n,-m)]*tmp1;
						
					}
				}
				prod*=val;
			}
			this->values[pd][_map(j,k)]+=csum;
			//debug:
			//std::cout << "j= " << j <<",k=" << k << ": csum=" << csum <<std::endl;
					
			}
			
		}
	}
	
}
//corrected r,theta,phi

template <typename Particle, int _order>
void _PotentialExpansions<Particle,_order>::evaluateExpansions (typename Particle::BaseType *location, typename Particle::RHSType *out_RHS)
{
	ExpansionsValueType Ynm[_order*(_order+2)]; 
	Btype Pnm[_order][_order];
	Btype rp[Particle::dim]; //spherical coordinates. relative to center.
	
	//calculate spherical coordinates relative to center of Box.
	cart2sph(location, this->Center, rp);
	
	//compute Legendre-Polynomial.	
	//calc Pnm
    calcPnm<Particle,_order>(Pnm,rp[1]);
	calcYnm(rp[2], Pnm, Ynm);
	
	Btype prod=rp[0];
	
	for (int n=0;n<_order;++n)
	{
		for (int m=-n;m<=n;++m)
			for (int d=0;d<Particle::pdim;++d)
				out_RHS->x[d]+=std::real(this->values[d][_map(n,m)]/prod*Ynm[_map(n,m)]);

		prod=prod*rp[0];
	}
}
			
			


			
			


