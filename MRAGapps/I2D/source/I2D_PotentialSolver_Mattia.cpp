/*
 *  I2D_PotentialSolver_Mattia.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Mattia Gazzola on 11/05/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <vector>
#include <xmmintrin.h>

#include "I2D_PotentialSolver_Mattia.h"

using namespace MRAG;
using namespace std;

struct GetTmp { static inline Real stream(FluidElement2D& out) { return out.tmp; } };

void I2D_PotentialSolver_Mattia::_updateBlocks()
{
	const BlockCollection<B>& coll = grid_ptr->getBlockCollection();
	
	vector<B*> destblocks;
	for(unsigned int iblock=0; iblock<vDest.size(); iblock++)
		destblocks.push_back( &coll[vDest[iblock].blockID] );
	
	UpdateBlocksPot update;
	update.blocks = destblocks;
	update.myvelblocks = my_velBlocks;	
	tbb::parallel_for(blocked_range<int>(0, vDest.size()), update, auto_partitioner());
	
	vDest.clear();
}

void I2D_PotentialSolver_Mattia::_count_sourceparticles()
{
	vector<BlockInfo> vInfo = grid_ptr->getBlocksInfo();
	const BlockCollection<B>& coll = grid_ptr->getBlockCollection();

	blockid2info.clear();
	for(vector<BlockInfo>::iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		blockid2info[it->blockID] = SourceParticlesInfo();

	ThresholdParticles<GetTmp,0> countparticles(tolParticle,scaling_factor, blockid2info);
	block_processing.process(vInfo, coll , countparticles);

	int curr = 0;
	for(map<int, SourceParticlesInfo>::iterator it=blockid2info.begin(); it!=blockid2info.end(); ++it)
	{
		it->second.start = curr;
		curr += it->second.nsource_particles;
	}

	nsource_particles = curr;
}

void I2D_PotentialSolver_Mattia::_collect_sourceparticles()
{
	vector<BlockInfo> vInfo = grid_ptr->getBlocksInfo();
	const BlockCollection<B>& coll = grid_ptr->getBlockCollection();

	srcparticles = new VelocitySourceParticle[nsource_particles];
	assert(srcparticles != NULL);

	ThresholdParticles<GetTmp,1> collectparticles(tolParticle,scaling_factor,blockid2info);
	collectparticles.destptr = srcparticles;
	block_processing.process(vInfo, coll, collectparticles);
}

void I2D_PotentialSolver_Mattia::compute_velocity()
{
	_count_sourceparticles();

	if (nsource_particles == 0) return;

	_collect_sourceparticles();
	_compute();
	_updateBlocks();
	_cleanup();
}

