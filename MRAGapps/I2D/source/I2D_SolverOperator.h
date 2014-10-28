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
#include "I2D_VelocityOperator.h"

class I2D_SolverOperator
{
protected:
	double viscosity;
	Grid<W,B>& grid;
	BlockProcessing block_processing;
	
	I2D_VelocityOperator * velsolver;
	
public:
	
	I2D_SolverOperator(Grid<W,B>& grid, double viscosity, I2D_VelocityOperator * velsolver_): grid(grid), viscosity(viscosity), velsolver(velsolver_)
	{
		//don't forget to specify nu!
		assert(viscosity>0);
	}
	
	virtual double estimate_largest_dt()  const =0;

	virtual void perform_timestep(double dt)=0;
};

class I2D_SolverOperator_4thOrder: public I2D_SolverOperator
{
	double _estimate_smallest_dt() const;
	
public:
	
	I2D_SolverOperator_4thOrder(Grid<W,B>& grid, double viscosity, I2D_VelocityOperator * velsolver_): I2D_SolverOperator(grid,viscosity,velsolver_) {}
	
	double estimate_largest_dt()  const;
	void perform_timestep(double dt);
};