/*
 *  I2D_DiffusionOperator.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 *	This operator integrate the diffusion term with a fourth order finite-difference scheme and RK2.
 *	It uses LTS (if the dt is big enough). Tested!
 *
 *	perform_timestep():
 *	IN: omega[0-2]
 *	TMP: tmp
 *	OUT: omega[0-2]
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"

class I2D_DiffusionOperator
{
protected:
	const double FC, DIM;
	double viscosity;
	int rhscounter;
	Grid<W,B>& grid;
	BlockProcessing block_processing;
	
public:
	
	I2D_DiffusionOperator(Grid<W,B>& grid, double viscosity, double FC): grid(grid), viscosity(viscosity), FC(FC), DIM(2.0)
	{
		//don't forget to specify nu!
		assert(viscosity>0);
		assert(FC>0.0 && FC<1.0);
	}
	
	virtual double estimate_largest_dt() const = 0;

	virtual double estimate_smallest_dt() const = 0;

	virtual void perform_timestep(double dt) = 0;
	
	void set_viscosity(double nu){ viscosity = nu; };
	
	int get_nofrhs() const { return rhscounter; }
};

class I2D_DiffusionOperator_4thOrder: public I2D_DiffusionOperator
{
	
public:
	
	I2D_DiffusionOperator_4thOrder(Grid<W,B>& grid, double viscosity, double FC): I2D_DiffusionOperator(grid,viscosity,FC) {}
	
	double estimate_largest_dt()  const;

	double estimate_smallest_dt() const;

	void perform_timestep(double dt);	
};
