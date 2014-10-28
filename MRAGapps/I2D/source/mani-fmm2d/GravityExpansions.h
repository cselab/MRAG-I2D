/*
 *  GravityExpansions.h
 *  hcfmm
 *
 *  Created by Manfred on 4/21/09.
 *  Copyright 2009 ETHZ. All rights reserved.
 *
 */

#include <complex>
#include "hcfmm_types.h"
//#include "boost/math/special_functions/spherical_harmonic.hpp"






//maping function for the negative indices
template <int _order>
inline int _map(const int &n, const int &m)
{
	int res=n*n+(m+n);
	assert(res<_order*(_order+2));
	assert(res>=0);
	return res;
}


//debug: printing to file:
template <typename Ntype,int _order>
void dumpData(Ntype* data, char* fname )
{
	std::ofstream outf(fname);
	outf.setf(std::ios_base::scientific);
	outf.width(15);
	for (int n=0;n<_order;++n)
	{
		for (int m=-n;m<=n;++m)
		{
			outf  << data[_map<_order>(n,m)]<<std::endl;
		}
	}
	
	outf.close();
	
}


/*
//calc Pnm and Ynm should be ok, since checked with boost-implementation of Ynm.
template <class Particle, int _order>
void calcPnm(typename Particle::BaseType Pnm_out[_order][_order], const typename Particle::BaseType  &theta_in)
{
	const int order=_order;
	typedef typename Particle::BaseType Btype;
	Btype val=-sin(theta_in);
	Btype cosine=cos(theta_in);
	Btype prod=Btype(1);
	//calc Pnm
	for (int m=0;m<order;++m)
	{
		Pnm_out[m][m]=fracfac<Btype>(m)*prod;
		prod=prod*val;
	}
	
	for (int m=0;m<order-1;++m)
	{
		Pnm_out[m+1][m]=cosine*(2*m+1)*Pnm_out[m][m];
	}
	
	for (int n=2;n<order;++n)
	{
		val=cosine*(2*n-1);
		for (int m=0;m<=n-2;++m)
		{
			Pnm_out[n][m]=(val*Pnm_out[n-1][m]-(n+m-1)*Pnm_out[n-2][m])/Btype(n-m);
		}
	}
	
	//debug:
	//	std::cout << "Pnm is:" << std::endl;
	//	for (int j=0;j<_order;++j)
	//	{
	//    std::cout << "[ ";
	//		for (int i=0;i<_order;++i)
	//		{
	//			std::cout << Pnm_out[i][j] <<" " ;
	//		}
	//		std::cout << "]" << std::endl;
	//	}
}
*/
/*
//Boost Version:
template <class Particle, int _order>
void calcPnm(typename Particle::BaseType Pnm_out[_order][_order], const typename Particle::BaseType  &theta_in)
{
	typedef typename Particle::BaseType Btype;
	Btype cosine;
	cosine=cos(theta_in);
	
	for (int n=0;n<_order;++n)
		for (int m=0;m<=n;++m)
		{
			Pnm_out[n][m]=boost::math::legendre_p(n,m,cosine);
		}
	
}
*/


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


//TODO: improve this withoug using a temporary Pnmtmp. (derive formula directly including pm1n)
template<typename Particle, int _order>
void calcDPnm(typename Particle::BaseType DPnm_out[_order][_order], const typename Particle::BaseType Pnm_in[_order][_order])
    {
		typename Particle::BaseType Pnmtmp[_order][_order];
		typename Particle::BaseType prefac;
		
		for (int n=0;n<_order;++n)
		{
			for (int m=0;m<=n;++m)
			{
				Pnmtmp[n][m]=pm1n(m)*Pnm_in[n][m]; //(-1^m-factor) //copying to Pnmtmp which included (-1)^m
			}
		}
		
		DPnm_out[0][0]=0;
		for (int n=1;n<_order;++n)
		{
			DPnm_out[n][0]=-Pnmtmp[n][1];
			for (int m=1; m<n;++m)
			{
				prefac=(n+m)*(n-m+1);
				DPnm_out[n][m]=pm1n(m)*0.5*(prefac*Pnmtmp[n][m-1]-Pnmtmp[n][m+1]);
			}
			DPnm_out[n][n]=pm1n(n)*n*Pnmtmp[n][n-1];
			
		}
		//debug:
//		std::cout << "Pnm is: " << std::endl;
//		for (int i=0;i<_order;++i)
//		{
//		std::cout << "[ ";
//			for (int j=0;j<=i;++j)
//			{
//				std::cout << Pnm_in[i][j] <<" " ;
//			}
//			std::cout << "]" << std::endl;
//		}
//		std::cout << "DPnm is: " << std::endl;
//		for (int i=0;i<_order;++i)
//		{
//			std::cout << "[ ";
//			for (int j=0;j<=i;++j)
//			{
//				std::cout << DPnm_out[i][j] <<" " ;
//			}
//			std::cout << "]" << std::endl;
//		}
		
		
	}
	




template<typename NumType, int _order>
void calcYnm(NumType phi, NumType Pnm[_order][_order], std::complex<NumType> Ynm[_order*(_order+2)] )
{
	
	for (int n=0;n<_order;++n)
	{
		Ynm[_map<_order>(n,0)]=Pnm[n][0];  //special case n=0: sqrtfac=1, cos(0)=1, sin(0)=0
		
		for (int m=1;m<=n;++m)
		{
			Ynm[_map<_order>(n,m)]=sqrtfac<NumType,_order>(n,m)*Pnm[n][m]*exp(m*phi*std::complex<NumType>(0,1));
			Ynm[_map<_order>(n,-m)]=std::conj(Ynm[_map<_order>(n,m)]);

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
struct _GravityExpansions {
	typedef _GravityExpansions<Particle,_order> SelfType;
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
	
	void _calcNablaCnm(ExpansionsValueType NablaCnm[3],Btype Pnm[_order][_order], Btype DPnm[_order][_order], ExpansionsValueType Ynm[_order*(_order+2)],   const Btype rp[3], const int &n, const int& m );

	
	void clear()
	{
		for (int pd=0;pd<Particle::pdim;++pd)
		{
			for (int n=0;n<_order;++n)
			{
				for (int m=-n;m<=n;++m)
				{
					values[pd][_map<_order>(n,m)]=0;
				}
			}
		}

		
	}

	
	
	_GravityExpansions():Radius(0)
	{
		clear();
		
		for (int d=0; d<Particle::dim;++d)
		{
			Center[d]=0;
		}
	}
	
	_GravityExpansions(const SelfType &B)
	{
		
		for (int pd=0;pd<Particle::pdim;++pd)
		{
			for (int n=0;n<_order;++n)
			{
				for (int m=-n;m<=n;++m)
				{
					values[pd][_map<_order>(n,m)]=B.values[pd][_map<_order>(n,m)];
				}
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
			for (int m=-n;m<=n;++m)
			{
				for (int d=0;d<Particle::pdim;++d)
				{
					std::cout << "expansion [" << d << "][" << n << "][" <<m << "]="<<values[d][_map<_order>(n,m)] <<std::endl; 
				}
			}
		}
	}
	
	
};



//a new inline function: calcNablaCnm: 
//it calculates the coefficients to express the derivative d()/dx,d()/dy,d()/dz of Phi=laplace(U) which can then be multiplied with the coefficients Mnm[:]
//perhaps a macro would be better, and storing sin(phi),sin(theta) etc..
template<typename Particle, int _order>
inline void _GravityExpansions<Particle,_order>::_calcNablaCnm(ExpansionsValueType NablaCnm[3],typename Particle::BaseType Pnm[_order][_order], typename Particle::BaseType DPnm[_order][_order], ExpansionsValueType Ynm[_order*(_order+2)],   const typename Particle::BaseType rp[3],const int &n,const int &m  )
{
	//abort();
	const Btype epstheta=1e-5; //epsilon to catch singularity in 1/sin(theta)
	const typename Particle::BaseType& R=rp[0];
	const typename Particle::BaseType& theta=rp[1];
	const typename Particle::BaseType& phi=rp[2];
	typename Particle::BaseType C0;
	ExpansionsValueType dUdr,dUdphi,dUdtheta;
	ExpansionsValueType imag(0,1),tmpZ;
	
	if((fabs(theta) < epstheta) || (fabs(M_PI-theta) < epstheta) )
	{
		/*std::cout << "not yet implemented" <<std::endl;
		exit(1);*/
		
		//Check if singularity is at pi or at 0 (has effect on sign of some of the terms)
		int signum;
		if(fabs(theta) < epstheta)
		{
			signum=1;
		}
		else {
			signum=-1;
		}

		C0=1.0/(pow(R,Btype(n+2)));
		typename Particle::BaseType pfac;
		const typename Particle::BaseType nmterm=0.5*sqrt(n*(1.0+n));
		switch (m) {
			case 0:
				NablaCnm[0]=0;
				NablaCnm[1]=0;
				NablaCnm[2]=-C0*cos(theta)*(n+1)*pow((Real)signum,n);
				break;
			case 1:
				pfac=pow((Real)signum,n+1);
				NablaCnm[0]=-pfac*nmterm*C0;
				NablaCnm[1]=-imag*pfac*nmterm*C0;
				NablaCnm[2]=0;
				break;
			case -1:
				pfac=pow((Real)signum,n+1);
				NablaCnm[0]=-pfac*nmterm*C0;
				NablaCnm[1]=imag*pfac*nmterm*C0;
				NablaCnm[2]=0;
				break;
			default:
				NablaCnm[0]=0;
				NablaCnm[1]=0;
				NablaCnm[2]=0;
				break;
		}

		
	}
	else
	{
		C0=1.0/(pow(R,Btype(n+2)));
		dUdr=Ynm[_map<_order>(n,m)]*Btype(-n-1);
		tmpZ=Btype(m)*imag;
		dUdphi=(tmpZ)*Ynm[_map<_order>(n,m)]/(ExpansionsValueType)sin(theta); 
		dUdtheta=sqrtfac<typename Particle::BaseType,_order>(n,m)*exp(phi*tmpZ)*DPnm[n][abs(m)];
		NablaCnm[0]=(ExpansionsValueType)(sin(theta)*cos(phi))*dUdr+(ExpansionsValueType)(cos(theta)*cos(phi))*dUdtheta-(ExpansionsValueType)(sin(phi))*dUdphi;
		NablaCnm[1]=(ExpansionsValueType)(sin(theta)*sin(phi))*dUdr+(ExpansionsValueType)(cos(theta)*sin(phi))*dUdtheta+(ExpansionsValueType)(cos(phi))*dUdphi;	
		NablaCnm[2]=(ExpansionsValueType)(cos(theta))*dUdr-(ExpansionsValueType)(sin(theta))*dUdtheta;
		NablaCnm[0]*=C0;
		NablaCnm[1]*=C0;
		NablaCnm[2]*=C0;
		//std::cout << "NablaC (r,theta,phi,n,m): (" <<R << "," << theta << "," << phi << ","<< n << "," << m << " ) " <<  " ( " << NablaCnm[0] <<  "," << NablaCnm[1] << "," << NablaCnm[2] << ")" <<std::endl;
	}
	

}



template <typename Particle, int _order>
void _GravityExpansions<Particle,_order>::calculateExpansions(Particle* in_particles, int nParticles)
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
//				std::cout << n << " " << m  << " " << Ynm[_map<_order>(n,m)] <<std::endl;
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
					//Btype A=Anm(n,m);
					//int pm1=pm1n(n);
					//this->values[d][_map<_order>(n,m)]+=(p[i].w[d]*prod*Ynm[_map<_order>(n,m)]*A)/Btype(pm1);
					//Epton?:
					//this->values[d][_map<_order>(n,m)]+=(p[i].w[d]*prod*Ynm[_map<_order>(n,m)]);
					//afterMatlab:
					this->values[d][_map<_order>(n,m)]+=(p[i].w[d]*prod*Ynm[_map<_order>(n,-m)]);

				}
			}
			prod=prod*val;
			//debug:
			//std::cout << "debugging expansions:  i=" << i << ", n=" << n << ", prod(rp^n)= " <<prod << " w[i]" << p[i].w[0] << std::endl;
		}
		
		
	} //end particle-loop.
	
	//calculate (-m) expansions.
	//not-doing for after matlab
//	for (int n=0;n<_order;++n)
//	{
//		for (int m=1; m<=n;++m)
//		{
//			for (int d=0;d<Particle::pdim;++d)
//			{
//				this->values[d][_map<_order>(n,-m)]=std::conj(this->values[d][_map<_order>(n,m)]);
//			}
//		}
//	}

// Epton DOES NOT DO THAT!	
//	for (int m=1;m<_order;++m)
//	{
//		//CI=(0,1) (not -1)
//		std::complex<Btype> tmp=pow(std::complex<Btype>(0,1),-m);
//		for (int n=m;n<_order;++n)
//		{
//			for (int d=0;d<Particle::pdim;++d)
//			{
//				this->values[d][_map<_order>(n,m)]*=tmp;
//				this->values[d][_map<_order>(n,-m)]*=tmp;
//			}
//		}
//	}
}
//reread...done
//reread2--> found bug in cart2sph->changed order of r,phi,theta to r,theta,phi!



//Another Implementation, after reading deeply into the papers, doing a Matlab-prototype ... Wed. 1.April.2009
template <typename Particle, int _order>
void _GravityExpansions<Particle,_order>::gatherExpansions (_GravityExpansions<Particle,_order>* in_expansions)
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
						//csum+=in_expansions->values[pd][_map<_order>(pr,qr)]*pow(CI,tmpexp)*pow(rb[0],n)*Ynm[_map<_order>(n,-m)]*tmp1;
						csum+=in_expansions->values[pd][_map<_order>(pr,qr)]*pow(CI,tmpexp)*prod*Ynm[_map<_order>(n,-m)]*tmp1;
						
					}
				}
				prod*=val;
			}
			this->values[pd][_map<_order>(j,k)]+=csum;
			//debug:
			//std::cout << "j= " << j <<",k=" << k << ": csum=" << csum <<std::endl;
					
			}
			
		}
	}
	
}
//corrected r,theta,phi





//THIS IS DIFFERENT for GravityParticle (Force Evaluation)
template <typename Particle, int _order>
void _GravityExpansions<Particle,_order>::evaluateExpansions (typename Particle::BaseType *location, typename Particle::RHSType *out_RHS)
{
	
#ifndef TEST_0
	ExpansionsValueType Ynm[_order*(_order+2)]; 
	ExpansionsValueType Outer[_order*(_order+2)];
	Btype Pnm[_order][_order];
	Btype DPnm[_order][_order];
	ExpansionsValueType NablaCnm[3]; //the Nablaoperator in 3D-Cartesian coordinates  (simplifies to coefficients in Spherical harmonics)
	Btype rp[Particle::dim]; //spherical coordinates. relative to center.
	Btype cosine,sine,val;
	//Btype contris[_order*(_order+2)];
	char buf[22];
	const double factorPI =  0.25/M_PI;

	
	//calculate spherical coordinates relative to center of Box.
	//debug:
	//std::cout << "calculating spherical coordinates relative to center: [" << this->Center[0] << " " << this->Center[1] <<" "<<this->Center[2] << "]" <<std::endl;
//	std::cout << "evaluating at location: " << location[0] << ", "<< location[1] <<", "<< location[2] <<std::endl;
	cart2sph(location, this->Center, rp);
//	std::cout << "rp: [" << rp[0] << " " << rp[1] <<" "<<rp[2]<< "]" <<std::endl;
	
	//compute Legendre-Polynomial.
	sine=sin(rp[1]);
	cosine=cos(rp[1]);
	val=-sine;
	
	//calc Pnm
    calcPnm<Particle,_order>(Pnm,rp[1]);
	calcDPnm<Particle,_order>(DPnm,Pnm); //derivatives of legendre polynomials
	calcYnm(rp[2], Pnm, Ynm);
	
	
	
	//for (int n=_order-1;n>=0;--n)
	for (int n=0;n<_order;++n)
	{
		//prod=pow(rp[0],n+1);
		for (int m=-n;m<=n;++m)
		{
			
			/*SelfType::*/
			_calcNablaCnm(NablaCnm,Pnm,DPnm,Ynm,rp,n,m);
			//void _calcNablaCnm(Btype NablaCnm[3],Btype Pnm[_order][_order], Btype DPnm[_order][_order], Btype Ynm[_order*(_order+2)],   const Btype rp[3], const int &n, const int& m );

			out_RHS->x[0]+=factorPI*std::real(NablaCnm[0]*this->values[0][_map<_order>(n,m)]); //assume only one field quantitiy.
			out_RHS->x[1]+=factorPI*std::real(NablaCnm[1]*this->values[0][_map<_order>(n,m)]);
			out_RHS->x[2]+=factorPI*std::real(NablaCnm[2]*this->values[0][_map<_order>(n,m)]);
			//debug:
//			std::cout << "NablaCnm (r,theta,phi,n,m): (" << rp[0] << "," << rp[1] << "," << rp[2] << ","<< n << "," << m << " ) " <<  " ( " << NablaCnm[0] <<  "," << NablaCnm[1] << "," << NablaCnm[2] << ")" <<std::endl;
//			std::cout <<  "Contrib. (r,theta,phi,n,m): (" <<rp[0] << "," << rp[1] << "," << rp[2] << ","<< n << "," << m << " ) " <<  " ( " << std::real(NablaCnm[0]*this->values[0][_map<_order>(n,m)]) <<  "," << std::real(NablaCnm[1]*this->values[0][_map<_order>(n,m)]) << "," << std::real(NablaCnm[2]*this->values[0][_map<_order>(n,m)]) << ")" <<std::endl;
//			std::cout << "Psum. (r,theta,phi,n,m): (" <<rp[0] << "," << rp[1] << "," << rp[2] << ","<< n << "," << m << " ) " <<  " ( " << out_RHS->x[0] <<  "," << out_RHS->x[1] << "," << out_RHS->x[2] << ")" <<std::endl;

		}
		
		
	}
	
	
	
	//debug:
	//sprintf(buf,"contributions_o%2.2i.dat",_order);
	//dumpData(contris, buf);
	
	
	//calc outer expansion: //checked against Epton/Dembart it is ok.
//	prod=Btype(1);
//    for (int n=0;n<_order;++n)
//	{
//		prod=prod*rp[0]; //calculate r^(n+1)
//		for (int m=-n;m<=n;++m)
//		{
//			//Outer[n][mp]=pow(-1,n)*pow(ExpansionsValueType(0,1),abs(m))*Ynm[n][mp]/(Anm(n,m)*prod);
//			Btype pw=Btype(pm1n(n)); //calculate (-1)^n
//			Outer[_map<_order>(n,m)]=pw*pow(ExpansionsValueType(0,1),abs(m))*Ynm[_map<_order>(n,m)]/(Anm(n,m)*prod);
//			
//		}
//	}
//	

//	//debug:
//	std::cout << "Outer Expansion Onm(theta,phi): theta: " << rp[1] << " phi:" << rp[2] <<std::endl; 
//	for (int n=0;n<_order;++n)
//	{
//		for (int m=-n;m<=n;++m)
//		{
//			std::cout << n << ", " << m << ": " << Outer[_map<_order>(n,m)] <<std::endl;
//		}
//	}
	
	//Sum up the different moments:
//	for (int n=0;n<_order;++n)
//	{
//		for (int m=-n;m<=n;++m)
//		{
//			for (int d=0;d<Particle::pdim;++d)
//			{
//				//ToDo: why is it Outer(n,-m)-->it is like that in Epton/dembart eq. 2.39.
//				out_RHS->x[d]+=std::real(this->values[d][_map<_order>(n,m)]*Outer[_map<_order>(n,-m)]);
//			}
//		}
//	}
	
#else	
	for (int pd=0;pd<Particle::pdim;++pd)
	{
		out_RHS->x[pd]+=std::real(this->values[pd][0]); //for testing reasons:mass summation
	}
#endif
}
			
			


			
			


