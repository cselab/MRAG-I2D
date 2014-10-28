/*
 *  I2D_CoreFMM_Check.cpp
 *  I2D_ROCKS
 *
 *  Created by Roman Schaerer on 12/26/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "I2D_CoreFMM_AggressiveVel.h"

class I2D_CoreFMM_Check : public I2D_CoreFMM_AggressiveVel
{
public:

	virtual void solve(const Real theta, const Real inv_scaling, BlockInfo * dest, const int nblocks, VelocitySourceParticle * srcparticles, const int nparticles);

};
