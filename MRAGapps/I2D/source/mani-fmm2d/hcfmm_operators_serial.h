/*
 *  hcfmm_operators_serial.h
 *  hcfmm
 *
 *  Created by Manfred on 12/23/08.
 *  Copyright 2008 ETHZ. All rights reserved.
 *
 */

#include "hcfmm_types.h"
#pragma once
template <class BaseType, class Particle, int sDim>
bbox<BaseType,sDim> getBoundingBox(Particle* vParticles, int n)
{
	bbox<BaseType,sDim> res;
	
	for (int i=0; i<n;++i)
	{
		for (int d=0; d<sDim;++d)
		{
			res.upper[d]=max(res.upper[d],vParticles[i].x[d]);
			res.lower[d]=min(res.lower[d],vParticles[i].x[d]);
		}
	}
	return res;
}
