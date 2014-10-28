/*
 *  I2D_PenalizationOperator.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 *	perform_timestep() integrates the penalization term with an implicit Euler, updates u[0-2] and
 *	omega[0-2]. compute_dragandstuff() computes the forces using tmp and u[0-2]. Tested!
 *
 *	perform_timestep():
 *	IN: u[0-2], omega[0-2], tmp as Xs
 *	TMP: tmp, psi, external
 *	OUT: u[0-2], omega[0-2]
 *
 *	compute_dragandstuff():
 *	IN: u[0-2], tmp as Xs
 *	OUT: (nothing)
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"

class I2D_PenalizationOperator
{
	Real lambda, Cd,Cl;
	Real Uinf[2];
	
	Grid<W,B>& grid;
	BlockProcessing block_processing;
	
	bool bAppendToFile;
	
	const map<int, const VelocityBlock*> * desired_velocities;
	
public:
	
	I2D_PenalizationOperator(Grid<W,B>& grid, Real lambda, const Real Uinf[2], bool bAppendToFile=false): 
		grid(grid), lambda(lambda), bAppendToFile(bAppendToFile), desired_velocities(NULL)
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
		
	}
	
	void set_desired_velocity(const map<int, const VelocityBlock*> * desired_velocities_)
	{
		this->desired_velocities = desired_velocities_;
	}
	
	void perform_timestep(Real dt);
	
	void compute_dragandstuff(Real time, const Real D, const Real cor[2], string filename);
	
	void set_Uinfinity(const Real _Uinf[2])
	{
		this->Uinf[0] = _Uinf[0];
		this->Uinf[1] = _Uinf[1];
	};
	
	Real getDrag(){ return this->Cd; }
	Real getLift(){ return this->Cl; }
	Real getLambda(){ return this->lambda; }
};
