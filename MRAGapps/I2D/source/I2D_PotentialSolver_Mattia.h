/*
 *  I2D_PotentialSolver_Mattia.h
 *  IncompressibleFluids2D
 *
 *  Created by Mattia Gazzola on 11/05/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 */

#pragma once

#include "I2D_VelocitySolver_Mani.h"

class I2D_PotentialSolver_Mattia: public I2D_VelocitySolver_Mani
{
protected:

	void _updateBlocks();
	void _collect_sourceparticles();
	void _count_sourceparticles();

public:

	I2D_PotentialSolver_Mattia(Grid<W,B>& grid, ArgumentParser& parser_): I2D_VelocitySolver_Mani(grid,parser_) {}

	void compute_velocity();
};
