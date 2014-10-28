/*
 *  I2D_VelocitySolverMPI_Mani.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 *	This operator takes omega[0-2] and compute the velocity
 *	in u[0-2]. It allocates and uses external memory. Tested!
 *
 *	IN:		omega[0-2]
 *	OUT:	velocity[0-2]
 *
 */

#pragma once

#ifdef _I2D_MPI_

#include "I2D_VelocitySolver_Mani.h"
#include "I2D_Headers.h"
#include "I2D_Types.h"

class I2D_VelocitySolverMPI_Mani: public I2D_VelocitySolver_Mani
{
	vector<int> workIDstart2node, work2node;
	int comm_rank, comm_size, stepid, ntotaldest_blocks;
	
	vector<BlockInfo> mydestinfo;
	vector< VelocityBlock *> alldestblocks;
	VelocityBlock * mydestblocks;
	
	void _wait_for_work();
	void _bcast_header();
	void _bcast_sourcedata();
	void _compute();
	void _master_update_blocks(int rank);
	void _collect_results();
	
public:
	
	I2D_VelocitySolverMPI_Mani(const int argc, const char ** argv);
	
	void set_grid(Grid<W,B>& grid)
	{
		set_grid_ptr(&grid);
	}
	
	void compute_velocity();
};

#endif

