/*
 *  I2D_CoreFMM_Check.cpp
 *  I2D_ROCKS
 *
 *  Created by Roman Schaerer on 12/26/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#include <stdio.h>
#include <tbb/parallel_sort.h>

#include "I2D_CoreFMM_Check.h"

#include "I2D_CoreFMM_AggressiveVel.h"
#include "I2D_CoreFMM_PlanBuilder.h"

//#include "I2D_CoreFMM_GPU.h"
#include "I2D_CoreFMM_SSE.h"

typedef HCFMM::Box<_VortexExpansions<VelocitySourceParticle, _ORDER_>,_FMM_MAX_LEVEL_> tBox;
typedef HCFMM::boxBuilder_serial<_VortexExpansions<VelocitySourceParticle, _ORDER_>,_FMM_MAX_LEVEL_> tBoxBuilder;


void I2D_CoreFMM_Check::solve(const Real theta, const Real inv_scaling, BlockInfo * dest, const int nblocks, VelocitySourceParticle * srcparticles, const int nparticles) {
	std::cout << "\nrun I2D_CoreFMM_Check::solve (...)\n";

	//I2D_CoreFMM_GPU gpu_solver;
	I2D_CoreFMM_SSE sse_solver;
	I2D_CoreFMM_AggressiveVel cpu_solver;

	/*gpu_solver.solve(theta, inv_scaling, dest, nblocks, srcparticles, nparticles);

	VelocityBlock* gpu_result = VelocityBlock::allocate(nblocks);

	for (int i=0;i<nblocks;++i) {
		assert (dest[i].ptrBlock != NULL);
		gpu_result[i] = *(VelocityBlock*)dest[i].ptrBlock;
	}*/
	VelocityBlock* gpu_result = VelocityBlock::allocate(nblocks);


	sse_solver.solve(theta, inv_scaling, dest, nblocks, srcparticles, nparticles);

	VelocityBlock* sse_result = VelocityBlock::allocate(nblocks);

	for (int i=0;i<nblocks;++i) {
		assert (dest[i].ptrBlock != NULL);
		sse_result[i] = *(VelocityBlock*)dest[i].ptrBlock;
	}

	cpu_solver.solve(theta, inv_scaling, dest, nblocks, srcparticles, nparticles);

	double max_rel_error_gpu, l1_rel_error_gpu;
	double max_rel_error_sse, l1_rel_error_sse;
	double tmp;
	unsigned long counter_gpu, counter_sse;
	counter_gpu = counter_sse = 0;

	max_rel_error_gpu = l1_rel_error_gpu = max_rel_error_sse = l1_rel_error_sse = 0.0;

	for (int i=0;i<nblocks;++i) {
		assert (dest[i].ptrBlock != NULL);
		VelocityBlock& b = *(VelocityBlock*)dest[i].ptrBlock;
		for (int y=0; y<_BLOCKSIZE_; ++y) {
			for (int x=0; x<_BLOCKSIZE_; ++x) {
				if (std::abs(b.u[0][y][x])  > std::numeric_limits <double>::epsilon ()) {
					tmp = std::abs (gpu_result[i].u[0][y][x] - b.u[0][y][x])/std::abs(b.u[0][y][x]);
					max_rel_error_gpu = std::max (max_rel_error_gpu, tmp);
					l1_rel_error_gpu += tmp;
					counter_gpu += 1;
				}
				if (std::abs(b.u[1][y][x])  > std::numeric_limits <double>::epsilon ()) {
					tmp = std::abs (gpu_result[i].u[1][y][x] - b.u[1][y][x])/std::abs(b.u[1][y][x]);
					max_rel_error_gpu = std::max (max_rel_error_gpu, tmp);
					l1_rel_error_gpu += tmp;
					counter_gpu += 1;
				}
				if (std::abs(b.u[0][y][x])  > std::numeric_limits <double>::epsilon ()) {
					tmp = std::abs (sse_result[i].u[0][y][x] - b.u[0][y][x])/std::abs(b.u[0][y][x]);
					max_rel_error_sse = std::max (max_rel_error_sse, tmp);
					l1_rel_error_sse += tmp;
					counter_sse += 1;
				}
				if (std::abs(b.u[1][y][x])  > std::numeric_limits <double>::epsilon ()) {
					tmp = std::abs (sse_result[i].u[1][y][x] - b.u[1][y][x])/std::abs(b.u[1][y][x]);
					max_rel_error_sse = std::max (max_rel_error_sse, tmp);
					l1_rel_error_sse += tmp;
					counter_sse += 1;
				}

			}
		}
	}
	l1_rel_error_sse /= counter_sse;
	l1_rel_error_gpu /= counter_gpu;

	std::cout << "\n\n//////////////// Error Analysis //////////////////\n\n" \
			"rel max error sse: " << setprecision (4) << scientific << max_rel_error_sse << "\trel l1 error sse: " << l1_rel_error_sse << "\n\n" \
			"//////////////////////////////////////////////////\n\n";

	VelocityBlock::deallocate(gpu_result);
	VelocityBlock::deallocate(sse_result);

}
