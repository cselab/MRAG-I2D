/*
 *  hcfmm_types.h
 *  hcfmm
 *
 *  Created by Manfred on 12/23/08.
 *  Copyright 2008 ETHZ. All rights reserved.
 *
 */
#pragma once

#include "math.h"
#include "precompute.h"
#include <cassert>


//min-max funcs:

/*template <class T>
T const& min(T const& x, T const& y)
{
	return x < y ? x : y;
}

template <class T>
T const& max(T const& x, T const& y)
{
	return x > y ? x : y;
}*/


//Bounding box type:
template <class BaseType, int sDim>
struct bbox
{
	
	BaseType lower[sDim];
	BaseType upper[sDim];
	
	bbox()
	{
		for (int i=0;i<sDim;++i)
		{
			lower[i]=BaseType(HUGE);
			upper[i]=-BaseType(HUGE);
		}
	}
	
	
	bbox(const bbox<BaseType,sDim> &incoming)
	{
		for (int i=0;i<sDim;++i)
		{
			lower[i]=incoming.lower[i];
			upper[i]=incoming.upper[i];
		}
	}
	
	void compare (const bbox<BaseType,sDim> &incoming)
	{
		for (int i=0;i<sDim;++i)
		{
			lower[i]=min(lower[i],incoming.lower[i]);
			upper[i]=max(upper[i],incoming.upper[i]);
		}
	}
	
	void print()
	{
		printf("Bounding Box: [  ");
		for (int d=0; d<sDim;++d)
		{
			printf("%f, ",BaseType(lower[d]));
		}
		
		printf(" / ");
		
		for (int d=0; d<sDim;++d)
		{
			printf("%f, ",BaseType(upper[d]));
		}
	    
	    printf(" ] \n");
		
	}
	
	bbox & operator=(const bbox& incoming )
	{
		for (int i=0;i<sDim;++i)
		{
			lower[i]=incoming.lower[i];
			upper[i]=incoming.upper[i];
		}
		return this;
	}
	
};

//particle index stuff:
//typedef unsigned int p_key;
//typedef ot_hkey p_key;
typedef unsigned short int tBoxKey;


template <class Particle>
struct p_ind
{
	tBoxKey key;
	Particle* p;
	
	//Bless Valgrind (don't need to change anything here):
	//ok I don't like this, but because I had problems with using vector:push_back(), I am overwriting the constructor /destructors:
	//	
	//	p_ind(): key(0), p(NULL) {
	//	//debug:
	//	//	std::cout << "creating a new index " <<std::endl;
	//	};
	//	
	//	~p_ind() {
	//		//debug:
	//		std::cout << "destroying an index, pointing to: " << p <<std::endl;
	//
	//		key=0; p=NULL;
	//		
	//	}
	//	
	//	//man do I really need this?:
	//	p_ind<Particle>& operator=(const p_ind<Particle> &b)
	//	{
	//		//debug:
	//		std::cout << "using the equal operator" <<std::endl;
	//		key=b.key;
	//		p=b.p;
	//		return *this;
	//	}
	//	
	
	
};


template <class Particle>
bool operator<(p_ind<Particle> P1, p_ind<Particle> P2)
{
	return (P1.key<P2.key);
}

template <class Particle>
bool operator<=(p_ind<Particle> P1, p_ind<Particle> P2)
{
	return (P1.key<=P2.key);
}


template <class Particle>
tBoxKey getboxKey(Particle* p, bbox<typename Particle::BaseType,Particle::dim>& bounds)
{
	tBoxKey res=0;
	for (int d=0; d<Particle::dim;++d)
	{
		(p->x[d]<(bounds.upper[d]-bounds.lower[d])/(typename Particle::BaseType(2)))?res=res:res+=(1<<d);
	}
	
	return res;
}

/*
 template <class Particle>
 p_ind<Particle> createpIndex(Particle* in_p, bbox<typename Particle::BaseType,Particle::dim>& bounds)
 {
 p_ind<Particle> res;
 res.p=in_p;
 res.key=0;
 for (int d=0; d<Particle::dim;++d)
 {
 (in_p->x[d]<(bounds.upper[d]-bounds.lower[d])/(typename Particle::BaseType(2)))?res.key=res.key:res.key+=(1<<d);
 }
 
 return res;
 }
 */
template <class Particle>
p_ind<Particle> createpIndex(Particle* in_p, typename Particle::BaseType *center, typename Particle::BaseType *h)
{
	
	typedef typename Particle::BaseType Btype; 
	p_ind<Particle> res;
	res.p=in_p;
	res.key=0;
	for (int d=0; d<Particle::dim;++d)
	{
		(in_p->x[d]<center[d])?res.key=res.key:res.key+=(1<<d);
	}
	
	
	return res;
}


inline void lsfkey2bits(int key, int dim, int* bits)
{
	for (int d=0; d<dim; ++d)
	{
		bits[d]=(int)((key>>d) & 1);
	}
}


/*
 inline unsigned long int factorial(unsigned long int n)
 {
 assert(n>=0);
 assert(n<13); //only valid for n<13 (overflow?)
 if(n!=0)
 {
 return n*factorial(n-1);
 }
 else
 {
 return 1;
 }
 }
 /*
 
 
 /*
 struct factorial
 {
 static const int minval=10;  //minimum precomputed values
 typedef unsigned long int ntype;
 //private:
 bool initialized;
 int &returnvalue;
 int mxvalue;
 ntype* values;
 
 factorial(int n=0):initialized(false),mxvalue(0),values(0)
 {
 int mxvalue=max(n,minval);
 values=new ntype[mxvalue+1];
 values[0]=1;
 values[1]=1;
 for (int i=2;i<=mxvalue;++i)
 {
 values[i]=values[i-1]*i;
 }
 initialized=true;
 }
 
 ntype operator()(ntype n)
 {
 if (n<mxvalue)
 {
 return values[n];
 }
 else
 {
 ntype* tmp=values;
 //create new array:
 values=new ntype[n+1];
 //copy new values into array.
 for (int i=0;i<=mxvalue;++i)
 {
 values[i]=tmp[i];
 }
 //calculate new values.
 for (int i=mxvalue+1;i<=n;++i)
 {
 values[i]=values[i-1]*i;
 }
 return values[n];
 }
 }
 };
 
 
 
 
 template <typename NumType>
 inline NumType Anm(int n, int m)
 {
 return pow(NumType(-1),n)/sqrt(factorial(n-m)*factorial(n+m));
 }
 */


//IMPROVE: AVOID THIS (only need it in the bad implementation of Pnm)
//template <typename NumType>
//inline unsigned long int fracfac(unsigned long int m)
//{
//	//assert(m<7);
//	unsigned long int tmp1=factorial(2*m);
//	unsigned long int tmp2=factorial(m);
//	unsigned long int tmp3=(1<<m);
//	//std::cout  << "[m,fac(2m),fac(m),(1<<m)]: [" << m <<"," << tmp1 <<","<< tmp2 <<","<< tmp3 <<"]" <<std::endl; 
//	return NumType(tmp1)/(NumType(tmp2*tmp3));
//	//return factorial(2*m)/NumType((1<<m)*factorial(m));
//}

template <typename NumType,int _order>
inline NumType sqrtfac(int n, int m)
{
	
	/*
	 const int max_int=2*_order;
	 const int Nn=1;
	 const int Nd=1;
	 int FinNumerators[Nn]={n-abs(m)};
	 int FinDenominators[Nd]={n+abs(m)};
	 int cNum[max_int+1]={}; //include 0, allow negative numbers
	 NumType Result(0);
	 
	 int cmax=0;
	 
	 for (int i=0;i<Nn;++i){ //loop over elements in numerator
	 for (int j=0;j<=FinNumerators[i];++j){
	 assert(j<=max_int);
	 cNum[j]+=1;
	 cmax=max(cmax,j);
	 }
	 }
	 
	 for (int i=0;i<Nd;++i){ //loop over elements in denominator
	 for (int j=0;j<=FinDenominators[i];++j){
	 assert(j<=max_int);
	 cNum[j]-=1;
	 cmax=max(cmax,j);
	 }		
	 }
	 
	 
	 //Further clean up the integerlist (something like a prime factorization)
	 for (int k=cmax;k>3;k--) //Starting from the biggest numbers
	 {
	 int p=2; //smallest divisor is 2 
	 while((cNum[k]!=0) && (p<k/2)  ) //if we have an non-zero entry and p is between 2 and k/2
	 {
	 if (k%p==0)                  //if it k is dividable by p
	 {
	 cNum[p]+=cNum[k];  //split it up into pieces and add it to the other stuff.
	 cNum[k/p]+=cNum[k];
	 cNum[k]=0;
	 }
	 p++;
	 }
	 }
	 
	 while(cNum[cmax]==0) //find new cmax, if it changed
	 {cmax--;}
	 
	 //This implementation uses the logarithmic.--> have always small number 
	 //a different way would be to find the ggT and eliminate as many terms as possible -->would only use integer as long as possible
	 for (int i=1;i<=cmax;++i){
	 Result+=cNum[i]*log(NumType(i));
	 }
	 //std::cout << "FoF:  before exp: " << Result << std::endl;
	 return sqrt(exp(Result));
	 
	 */
	
	//IMPROVE THIS: don't really need to calculate the high factorials (worst case up to 2*n) but only n
	return sqrt(factorial(n-abs(m))/NumType(factorial(n+abs(m))));
	
}



//warning: this is only for 3d!
template <typename NumType>
inline void cart2sph(NumType* in_x, NumType* x0,NumType* out_x)//r,theta,phi !!!!! it is r,theta,phi NOT r,phi,theta
{
	NumType R(0);
	NumType dx[3]={0,0,0};
	for (int d=0;d<3;++d)
	{
		out_x[d]=0;
		dx[d]=in_x[d]-x0[d];
		R+=(in_x[d]-x0[d])*(in_x[d]-x0[d]);
	}
	if (R>NumType(0))
	{
		R=sqrt(R);
		out_x[0]=R;
		//out_x[1]=atan2(dx[2],sqrt(dx[0]*dx[0]+dx[1]*dx[1]));//acos((in_x[2]-x0[2])/R); //theta, the angle between 0->pi
		out_x[1]=acos((in_x[2]-x0[2])/R);
		out_x[2]=atan2(dx[1],dx[0]); //phi, the angle between 0->2pi
	}
}

template <typename NumType>
inline void sph2cart(NumType* in_sph ,NumType* out_x)//r,theta,phi !!!!! it is r,theta,phi NOT r,phi,theta
{
	out_x[0]=0;
	out_x[1]=0;
	out_x[2]=0;
	NumType R=in_sph[0];
	if (R>NumType(0))
	{
		out_x[0]=R*sin(in_sph[1])*cos(in_sph[2]);
		out_x[1]=R*sin(in_sph[1])*sin(in_sph[2]);
		out_x[2]=R*cos(in_sph[1]);
	}
}



template <typename tint>
inline tint pm1n(const tint &n) //(return -1^n in a more efficient way (hopefully))
{
	return 1-2*(n%2);
}


//calculates the Term involving Fractions of Factorials, which is used in the shifting formula: 
//FracAnm=A_n^m*A_(j-n)^(k-m)/A_j^k
//Have simplified the term as much as possible with Mathematica:
//Sqrt[(j - k)! (j + k)!]/ ( Sqrt[(j + k - m - n)! (j - k + m - n)!] * Sqrt[(-m + n)! (m + n)!]) 
//template <class ReturnType,int _order>
//ReturnType calcFracAnm(int j, int k, int n, int m)
//{
//	
//	const int max_int=2*_order;
//	const int Nn=2;
//	const int Nd=4;
//	int FinNumerators[Nn]={j-k,j+k};
//	int FinDenominators[Nd]={j+k-m-n,j-k+m-n,n-m,n+m};
//	int cNum[max_int+1]={}; //include 0, allow negative numbers
//	ReturnType Result(0);
//	
//	int cmax=0;
//	
//	for (int i=0;i<Nn;++i){ //loop over elements in numerator
//		for (int j=0;j<=FinNumerators[i];++j){
//			assert(j<=max_int);
//			cNum[j]+=1;
//			cmax=max(cmax,j);
//		}
//	}
//	
//	for (int i=0;i<Nd;++i){ //loop over elements in denominator
//		for (int j=0;j<=FinDenominators[i];++j){
//			assert(j<=max_int);
//			cNum[j]-=1;
//			cmax=max(cmax,j);
//		}		
//	}
//	
//	
//	//Further clean up the integerlist (something like a prime factorization)
//	for (int k=cmax;k>3;k--) //Starting from the biggest numbers
//	{
//		int p=2; //smallest divisor is 2 
//		while((cNum[k]!=0) && (p<k/2)  ) //if we have an non-zero entry and p is between 2 and k/2
//		{
//			if (k%p==0)                  //if it k is dividable by p
//			{
//				cNum[p]+=cNum[k];  //split it up into pieces and add it to the other stuff.
//				cNum[k/p]+=cNum[k];
//				cNum[k]=0;
//			}
//			p++;
//		}
//	}
//	
//    while(cNum[cmax]==0) //find new cmax, if it changed
//	{cmax--;}
//	
//	//This implementation uses the logarithmic.--> have always small number 
//	//a different way would be to find the ggT and eliminate as many terms as possible -->would only use integer as long as possible
//	for (int i=1;i<=cmax;++i){
//		Result+=cNum[i]*log(ReturnType(i));
//	}
//	//std::cout << "FoF:  before exp: " << Result << std::endl;
//	return sqrt(exp(Result));
//}


//Trivial implementation (now uses the double-precision factorial.)
template <class ReturnType, int _order>
ReturnType calcFracAnm(int j, int k, int n, int m)
{
	int pr=j-n;
	int qr=k-m;
	assert(abs(pr)>=abs(qr));
	
	ReturnType Result;
	Result=Anm(n,m);
	Result/=Anm(j,k);
	Result*=Anm(j-n,k-m);
	return Result;
}
