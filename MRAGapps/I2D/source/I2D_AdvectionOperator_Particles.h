/*
 *  I2D_AdvectionOperator_Particles.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 *	this operator advect omega[0-2] with the velocity field u[0-2]
 *	with CFL particles. Tested!
 *
 *	IN: omega[0-2], u[0-2]
 *	OUT: omega[0-2]
 */
#pragma once

#include "I2D_AdvectionOperator.h"

class I2D_AdvectionOperator_Particles: public I2D_AdvectionOperator
{
	Real LCFL;
	
	Real _GTS_CFL(); //estimate the dt-CFL based on a global time stepping
	
public:
	I2D_AdvectionOperator_Particles(Grid<W,B>& grid, double CFL=0.25, double LCFL=0.25):
	I2D_AdvectionOperator(grid, CFL), LCFL(LCFL){}
	
	Real estimate_largest_dt();
	
	void perform_timestep(double dt);
};