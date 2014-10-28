/*
 *  I2D_CoreFMM_AggressiveVel.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_FMMTypes.h"

class I2D_CoreFMM_AggressiveVel
{
protected:
	int timestamp;
	
public:

	void setVerbose (bool _b_verbose) {
		b_verbose = _b_verbose;
	}

	bool isVerbose () {
		return b_verbose;
	}

	I2D_CoreFMM_AggressiveVel (bool _b_verbose = false) : b_verbose (_b_verbose) {}

	virtual void solve(const Real theta, const Real inv_scaling, BlockInfo * dest, const int nblocks, VelocitySourceParticle * srcparticles, const int nparticles);

	bool b_verbose;
};
