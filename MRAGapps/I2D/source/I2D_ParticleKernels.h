/*
 *  I2D_ParticleKernels.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

class MP4
{
public:
	
	inline static Real _f0(Real t)
	{
#ifndef NDEBUG
		assert(t>=0);
#endif
		return 1.+t*t*(-5./2.+3./2.*t);
	}
	
	inline static Real _f1(Real t)
	{
#ifndef NDEBUG
		assert(t>=0);
#endif
		return 2.+t*(-4. + t*(5./2.-1./2.*t));
	}
	
	static Real eval(Real x)
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
