/*
 *  I2D_AdvectionOperator.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 *	this operator integrate the advection term using a third-order upwind scheme.
 *	Tested!
 *
 *	IN: omega[0-2], u[0-2]
 *	OUT: omega[0-2]
 *
 */
#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"

class I2D_AdvectionOperator
{
protected:
	enum State {Initialized, Ready, Done};
	
	State state;
	int rhscounter;
	const double CFL;
	double largest_dt, smallest_dt, tmp_maxvel;

	Grid<W,B>& grid;
	BlockProcessing block_processing;
		
	Real Uinf[2];
		
public:
	
	I2D_AdvectionOperator(Grid<W,B>& grid, double CFL=0.25): 
		grid(grid), state(Initialized), CFL(CFL), rhscounter(0)
	{
		Uinf[0] = Uinf[1] = 0.0;
	}
	
	virtual Real estimate_largest_dt();
	
	virtual void perform_timestep(double dt);
	
	virtual void set_Uinfinity(const Real Uinf[2])
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
	}
	
	Real compute_maxvel();
	
	int get_nofrhs() const { return rhscounter; }
};
