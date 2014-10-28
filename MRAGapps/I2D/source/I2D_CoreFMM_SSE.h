/*
 *  I2D_CoreFMM_SSE.h
 *  I2D_ROCKS
 *
 *  Created by Roman Schaerer on 12/26/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "I2D_CoreFMM_AggressiveVel.h"

class I2D_CoreFMM_SSE : public I2D_CoreFMM_AggressiveVel {
public:

	void setVerbose (bool _b_verbose) {
		b_verbose = _b_verbose;
	}

	bool isVerbose () {
		return b_verbose;
	}

	I2D_CoreFMM_SSE (bool _b_verbose = false) : b_verbose (_b_verbose) {}

	virtual void solve(const Real theta, const Real inv_scaling, BlockInfo * dest, const int nblocks, VelocitySourceParticle * srcparticles, const int nparticles);

protected:

	bool b_verbose;

};
