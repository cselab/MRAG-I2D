/*
 * IF2DSmartInteractions.h
 *
 *  Created on: Feb 7, 2012
 *      Author: mgazzola
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_FluidMediatedInteractions.h"
#include "RL_TabularPolicy.h"
#include "I2D_ComputeEfficiency.h"

class I2D_SmartInteractions : public I2D_FluidMediatedInteractions
{
protected:
	Real TSTARTFTLE, TENDFTLE;
	Real TSTARTEFF, TENDEFF;

	I2D_ComputeEfficiency * efficiency;
	RL::RL_TabularPolicy * policy;

	void _prepareAgents();
	void _dispose();
	Real _initial_dt(int nsteps);

public:
	I2D_SmartInteractions(const int argc, const char ** argv);
	virtual ~I2D_SmartInteractions();
	void refresh();
	void run();
};

