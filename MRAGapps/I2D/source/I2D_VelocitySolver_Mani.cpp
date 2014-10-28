/*
 *  I2D_VelocitySolver_Mani.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <vector>
#include <xmmintrin.h>

#include "I2D_VelocitySolver_Mani.h"
#include "I2D_Clear.h"

using namespace MRAG;
using namespace std;

struct GetOmega { static inline Real stream(FluidElement2D& out) { return out.omega; } };

void I2D_VelocitySolver_Mani::_count_sourceparticles()
{
	vector<BlockInfo> vInfo = grid_ptr->getBlocksInfo();
	const BlockCollection<B>& coll = grid_ptr->getBlockCollection();
	
	blockid2info.clear();
	for(vector<BlockInfo>::iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		blockid2info[it->blockID] = SourceParticlesInfo();
	
	ThresholdParticles<GetOmega,0> countparticles(tolParticle,scaling_factor,blockid2info);
	block_processing.process(vInfo,coll,countparticles);
	
	int curr = 0;
	for(map<int, SourceParticlesInfo>::iterator it=blockid2info.begin(); it!=blockid2info.end(); ++it)
	{
		it->second.start = curr;
		curr += it->second.nsource_particles;
	}
	
	nsource_particles = curr;	
}

void I2D_VelocitySolver_Mani::_collect_sourceparticles()
{
	vector<BlockInfo> vInfo = grid_ptr->getBlocksInfo();
	const BlockCollection<B>& coll = grid_ptr->getBlockCollection();
	
	srcparticles = new VelocitySourceParticle[nsource_particles];
	assert(srcparticles != NULL);
	
	ThresholdParticles<GetOmega,1> collectparticles(tolParticle,scaling_factor,blockid2info);
	collectparticles.destptr = srcparticles;
	block_processing.process(vInfo, coll, collectparticles);
}

void I2D_VelocitySolver_Mani::_compute()
{		
	if (!bSKIPBLOCKS || nsource_particles == 0)
		vDest = grid_ptr->getBlocksInfo();
	else
	{
		vector<BlockInfo> vInfo = grid_ptr->getBlocksInfo();
				
		for(unsigned int i=0; i<vInfo.size(); i++)
		{
			assert(blockid2info.find(vInfo[i].blockID) != blockid2info.end());
			if (blockid2info[vInfo[i].blockID].nsource_particles > 0)
				vDest.push_back(vInfo[i]);
		}
		
		vector<BlockInfo> vNeighbors = grid_ptr->getNeighborsInfo(vDest);
	
		vDest.insert(vDest.end(), vNeighbors.begin(), vNeighbors.end());
		
		printf("FMM-SKIP ACTIVE. : %zu instead of %zu (%f)\n", vDest.size(), vInfo.size(), vDest.size()/(double)vInfo.size());
	}
	
	my_velBlocks = VelocityBlock::allocate(vDest.size());	
	
	for(unsigned int i=0;i<vDest.size();i++)
		vDest[i].ptrBlock = &my_velBlocks[i]; 

	coreFMM->solve(theta, 1.0/scaling_factor, &vDest.front(), vDest.size(), srcparticles, nsource_particles);
	
	for(unsigned int i=0;i<vDest.size();i++)
		vDest[i].ptrBlock = NULL;
}

void I2D_VelocitySolver_Mani::_updateBlocks()
{
	const BlockCollection<B>& coll = grid_ptr->getBlockCollection();
	
	vector<B*> destblocks;
	for(unsigned int iblock=0; iblock<vDest.size(); iblock++)
		destblocks.push_back( &coll[vDest[iblock].blockID] );
	
	UpdateBlocks update;
	update.blocks = destblocks;
	update.myvelblocks = my_velBlocks;	
	tbb::parallel_for(blocked_range<int>(0, vDest.size()), update, auto_partitioner());
	
	vDest.clear();
}

void I2D_VelocitySolver_Mani::_cleanup()
{
	delete [] srcparticles; srcparticles = NULL;
	nsource_particles = 0;
	
	VelocityBlock::deallocate(my_velBlocks);
}

void I2D_VelocitySolver_Mani::compute_velocity()
{	
	I2D_Clear cleaner;
	cleaner.clearVel(*grid_ptr);

	_count_sourceparticles();
	
	if (nsource_particles == 0) return;
	
	_collect_sourceparticles();
	_compute();
	_updateBlocks();
	_cleanup();
}
