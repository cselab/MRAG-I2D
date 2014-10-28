/*
 *  factorial.h
 *  hcfmm
 *
 *  Created by Manfred on 2/4/09.
 *  Copyright 2009 ETHZ. All rights reserved.
 *
 */

#ifndef _PRECOMPUTE_
#define _PRECOMPUTE_
#include <assert.h>
#include <math.h>
#include <iostream>

//#include "boost/math/special_functions/factorials.hpp"



struct tfactorial
{
	//typedef long unsigned  int ntype;
	typedef double ntype;
	//private:
	int mxvalue;
	ntype* values;
	
	tfactorial(int n=0):mxvalue(0),values(0)
	{
		(n>2)?mxvalue=n:mxvalue=2;
		values=new ntype[mxvalue+1];
		values[0]=1;
		values[1]=1;
		for (ntype i=2;i<=mxvalue;++i)
		{
			values[int(i)]=values[int(i-1)]*i;
			std::cout << " i " << i << " "  << values[int(i)] <<std::endl;
			//values[i]=boost::math::factorial<ntype>(i);
		}
	}
	
	ntype operator()(const int &n)
	{
		assert(n>=0);
		//assert(n<24);
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
			for (unsigned int i=0;i<=mxvalue;++i)
			{
				values[i]=tmp[i];
			}
			delete[] tmp;
			//calculate new values.
			for (unsigned int i=mxvalue+1;i<=n;++i)
			{
				values[i]=values[i-1]*i;
				//values[i]=boost::math::factorial<ntype>(i);
			}
			mxvalue=n;
			return values[n];
		}
	}
	
	~tfactorial()
	{
		delete[] values;
	}
};

extern tfactorial factorial;

/*
template <typename NumType>
inline NumType Anm(int n, int m)
{
	return pow(NumType(-1),n)/sqrt(factorial(n-m)*factorial(n+m));
}
*/

class tAnm
{
	//TODO: make this somehow generic (not easy)
	typedef double store_type;
	typedef double ntype;
	
private:
	store_type* values;
	int maxn;
	
	inline int map(int n, int m)
	{
		return n*n+(m+n);
	}
	
public:
	tAnm(int _n):maxn(_n),values(NULL)
	{
		
		values=new store_type[map(_n+1,_n+1)];
		
		for (int n=0;n<=_n;++n)
		{
			for(int m=-n;m<=n;++m)
			{
				values[map(n,m)]=pow(-1.,n)/sqrt(factorial(n-m)*factorial(n+m));
			}
		}
		//printMe();
	}
	
	ntype operator()(int _n, int _m)
	{
		//assert(abs(_m)<=_n);

		if((_n<=maxn))
		{
			return values[map(_n,_m)];
		}
		else
		{
			store_type* tmp= values;
			values=new store_type[map(_n+1,_m+1)];
			//copy old values:
			for (int i=0;i<=map(maxn,maxn);++i)
			{
				values[i]=tmp[i];
			}
			delete[] tmp;
			for (int n=maxn+1;n<=_n;++n)
			{
				for(int m=-n;m<=n;++m)
				{
					values[map(n,m)]=pow(-1.,n)/sqrt(factorial(n-m)*factorial(n+m));
				}
			}
			maxn=_n;
			return ntype(values[map(_n,_m)]);
		}
	}
	
	void printMe()
	{
		for (int n=0;n<=maxn;++n)
		{
			for(int m=-n;m<=n;++m)
			{
				std::cout << "Anm(" <<n<<","<<m<<")="<<values[map(n,m)]<<std::endl;
			}
		}
	}
	
	~tAnm()
	{
		delete[] values;
	}

	
};




extern tAnm Anm;
#endif