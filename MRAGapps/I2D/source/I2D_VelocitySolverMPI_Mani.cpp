/*
 *  I2D_VelocitySolverMPI_Mani.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#ifdef _I2D_MPI_

#include "I2D_VelocitySolverMPI_Mani.h"

#include <limits>
#include <mpi.h>
#include <tbb/parallel_for.h>
#include "I2D_CoreFMM_AggressiveVel.h"

namespace Velocity_MPI {
	vector<MPI_Request> requests;
}

struct ClearVelocity
{
	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{		
		FluidElement2D * const e = &b(0,0);
		
		static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
		for(int i=0; i<n; i++)			
			e[i].u[0] = e[i].u[1] = 0;
	}
};

void print(vector<int>& v, string s)
{
	printf("%s:\t", s.data());
	for(vector<int>::iterator it=v.begin(); it!=v.end(); it++)
		printf("%d ", *it);
	printf("\n");
}

void I2D_VelocitySolverMPI_Mani::_wait_for_work()
{
	while(true)
	{
		_bcast_header();
		_bcast_sourcedata();
		_compute();
		_collect_results();
		
		stepid++;
	}
}

void I2D_VelocitySolverMPI_Mani::_bcast_header()
{
	const int message_bytes = (2 + 2*comm_size)*sizeof(message_bytes);
	char * const message = new char[message_bytes];
	assert((((unsigned long int)message) & 0x7) == 0);
	
	if(comm_rank == 0)
	{
		int * ptr = (int*)message;
		
		*ptr = nsource_particles; ptr ++;
		*ptr = ntotaldest_blocks; ptr ++;
		
		for(int i=0;i<comm_size; i++, ptr++)
			*ptr = work2node[i];
		
		for(int i=0;i<comm_size; i++, ptr++)
			*ptr = workIDstart2node[i];
		
		//printf("sending...\n");
		//printf("nsource_particles=%d\n", nsource_particles);
		//print(work2node, "sending");
		//print(workIDstart2node, "sending");
	}
	
	MPI_Bcast(message, message_bytes, MPI_BYTE, 0, MPI_COMM_WORLD);
	
	if (comm_rank !=0)
	{
		//printf("receiving...\n");
		
		int * ptr = (int *)message;
		
		nsource_particles = *ptr; ptr++;
		ntotaldest_blocks = *ptr; ptr++;
		
		work2node = vector<int>(comm_size);
		for(int i=0;i<comm_size; i++, ptr++)
			work2node[i] = *ptr; 
		
		workIDstart2node = vector<int> (comm_size);
		for(int i=0;i<comm_size; i++, ptr ++)
			workIDstart2node[i] = *ptr;
		
		//printf("receiving nsource_particles=%d\n", nsource_particles);
		//print(work2node, "received");
		//print(workIDstart2node, "received");
	}
	
	delete [] message;
}

I2D_VelocitySolverMPI_Mani::I2D_VelocitySolverMPI_Mani(const int argc, const char ** argv):
I2D_VelocitySolver_Mani(argc, argv), stepid(0)
{
	MPI_Init(const_cast<int *>(&argc), const_cast<char***>(&argv));
	MPI_Comm_rank (MPI_COMM_WORLD, &comm_rank); /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &comm_size);
	//printf( "Hello world from process %d of %d\n", comm_rank, comm_size );
	
	ArgumentParser parser(argc, argv);
	
	parser.set_strict_mode();
	
	const int lmax = parser("-lmax").asInt();
	const double min_dV = pow(pow(0.5, lmax)/B::sizeX, 2);
	
	scaling_factor = max(1., 1e-4/min_dV);
	theta = parser("-fmm-theta").asDouble();
	tolParticle = 2*numeric_limits<Real>::epsilon();
	printf("I2D_VelocitySolverMPI_Mani: setting scaling_factor to %20.20e\n", scaling_factor);
	printf("I2D_VelocitySolverMPI_Mani: setting theta to %20.20e\n", theta);
	printf("I2D_VelocitySolverMPI_Mani: setting tolParticle to %20.20e\n", tolParticle);
	
	if (comm_rank != 0)
		_wait_for_work();
}

void I2D_VelocitySolverMPI_Mani::_compute()
{	
	assert(coreFMM != NULL);
	coreFMM->solve(theta, 1./scaling_factor, &mydestinfo.front(), work2node[comm_rank], srcparticles, nsource_particles);
}

void I2D_VelocitySolverMPI_Mani::_collect_results()
{
	if(comm_rank != 0)
	{
		Velocity_MPI::requests.resize(1);
		
		MPI_Isend(mydestblocks, mydestinfo.size()*sizeof(VelocityBlock), MPI_BYTE, 0, stepid, MPI_COMM_WORLD, &Velocity_MPI::requests.front());
		//printf("slave sent!\n");
		
		{
			MPI_Status status;
			int flag = 0;
			
			while(!flag)
			{
				int res  = MPI_Test(&Velocity_MPI::requests.front(), &flag, &status);
				assert(res == MPI_SUCCESS);
			}
			
			//printf("rank %d: cleanup now!\n", comm_rank);
		}
	}
	else 
	{
		Velocity_MPI::requests.resize(comm_size);
		set<int> pending_slaves;
		
		for(int slave=1; slave<comm_size; slave++)
		{
			////printf("master: from slave %d i expect
			MPI_Irecv(alldestblocks[slave], work2node[slave]*sizeof(VelocityBlock), MPI_BYTE, slave, stepid, MPI_COMM_WORLD, &Velocity_MPI::requests[slave]);
			pending_slaves.insert(slave);
		}
		
		for(int c=1; c<comm_size; c++)
		{
			int slave = -1;
			bool bReceived = false;
			
			while(!bReceived)
				for(set<int>::iterator it = pending_slaves.begin(); it!=pending_slaves.end(); it++)
				{
					MPI_Status status;
					int flag = 0;
					
					int res  = MPI_Test(&Velocity_MPI::requests[*it], &flag, &status);
					assert(res == MPI_SUCCESS);
					
					if (flag)
					{
						bReceived = true;
						slave = *it;
						break;
					}
				}
			
			pending_slaves.erase(slave);
			_master_update_blocks(slave);
		}
		
		assert(pending_slaves.size() == 0);
		
		//printf("rank %d: cleanup now!\n", comm_rank);
	}
	
	//cleanup
	delete [] srcparticles; srcparticles = NULL;
	nsource_particles = 0;
	
	//delete [] mydestblocks; mydestblocks = NULL;
	VelocityBlock::deallocate(mydestblocks);
	ntotaldest_blocks = 0;
	
	if(comm_rank == 0)
		for(int i=1; i<comm_size; i++)
		{
			VelocityBlock::deallocate(alldestblocks[i]);
			//delete [] alldestblocks[i];
			//alldestblocks[i] = NULL;
		}
	
	alldestblocks.clear();
	mydestinfo.clear();
	workIDstart2node.clear();
	work2node.clear();
}

void I2D_VelocitySolverMPI_Mani::_bcast_sourcedata()
{	 
	const long int blockinfo_bytes = 4*sizeof(short int);
	const long int particle_bytes = nsource_particles*sizeof(VelocitySourceParticle);
	const long int message_bytes = particle_bytes + ntotaldest_blocks*blockinfo_bytes;
	
	char * const message = new char[message_bytes];
	assert((((unsigned long int)message) & 0x7) == 0);
	assert(message != NULL);
	
	if (comm_rank == 0)
	{
		char * ptr = message;
		
		memcpy(ptr, srcparticles, particle_bytes);
		
		ptr += particle_bytes;
		
		vector<BlockInfo> vInfo = grid_ptr->getBlocksInfo();
		for(int iblock=0; iblock<ntotaldest_blocks; iblock++)
		{
			short int * info = (short int *)ptr;
			
			info[0] = vInfo[iblock].index[0];
			info[1] = vInfo[iblock].index[1];
			info[2] = 0;
			info[3] = vInfo[iblock].level;
			
			ptr += blockinfo_bytes;
		}
		
		//printf("sending data...\n");
	}
	
	MPI_Bcast(message, message_bytes, MPI_BYTE, 0, MPI_COMM_WORLD);
	
	mydestinfo.clear();
	
	if (comm_rank != 0)
	{
		//printf("receiving data for %d particles..\n", nsource_particles);
		
		const char * ptr = message;
		
		srcparticles = new VelocitySourceParticle[nsource_particles];
		assert(srcparticles != NULL);
		
		memcpy(srcparticles, ptr, particle_bytes );
		ptr += particle_bytes;
		
		const int start = workIDstart2node[comm_rank];
		const int end = start + work2node[comm_rank];
		
		mydestblocks = VelocityBlock::allocate(end-start);//new VelocityBlock[end-start];
		
		ptr += start*blockinfo_bytes;
		for(int iblock=start; iblock<end; iblock++)
		{
			const short int * const infoptr = (short int *)ptr;
			assert(message + message_bytes > ptr);
			
			const int idx[3] = {infoptr[0], infoptr[1], infoptr[2]};
			const int level = infoptr[3];
			
			ptr += blockinfo_bytes;
			
			BlockInfo info(-1, idx, level);
			
			const double dilate = pow(2.0, -info.level);
			
			const Real h[3] = {
				dilate/B::sizeX, 
				dilate/B::sizeY,
				dilate/B::sizeZ
			};
			
			info.h[0] = h[0];
			info.h[1] = h[1];
			info.h[2] = h[2];
			
			const Real p[3] = {
				info.index[0]*dilate + (W::bIsCellCentered ? 0.5*h[0] : 0),
				info.index[1]*dilate + (W::bIsCellCentered ? 0.5*h[1] : 0),
				info.index[2]*dilate + (W::bIsCellCentered ? 0.5*h[2] : 0),
			};
			
			info.origin[0] = p[0];
			info.origin[1] = p[1];
			info.origin[2] = p[2];	
			
			info.ptrBlock = &mydestblocks[iblock-start];
			
			mydestinfo.push_back(info);
		}
		
		//printf("received data\n");
	}
	else 
	{
		vector<BlockInfo> vInfo = grid_ptr->getBlocksInfo();
		
		const int start = workIDstart2node[comm_rank];
		const int end = start + work2node[comm_rank];
		mydestinfo.clear();
		mydestblocks = VelocityBlock::allocate(end-start);//new VelocityBlock[end-start];
		
		for(int iblock=start; iblock<end; iblock++)
		{
			BlockInfo info = vInfo[iblock];
			
			info.ptrBlock = &mydestblocks[iblock - start];
			
			mydestinfo.push_back(info);
		}
		
		alldestblocks = vector< VelocityBlock *>(comm_size);
		alldestblocks[0] = mydestblocks;
		
		for (int i=1; i<comm_size; i++) 
			alldestblocks[i] = VelocityBlock::allocate(work2node[i]);//new VelocityBlock[work2node[i]];
	}
	
	delete [] message;
}

void I2D_VelocitySolverMPI_Mani::_master_update_blocks(int rank)
{
	vector<BlockInfo> vInfo = grid_ptr->getBlocksInfo();
	const BlockCollection<B>& coll = grid_ptr->getBlockCollection();
	
	const int start = workIDstart2node[rank];
	const int end = start + work2node[rank];
	
	vector<B*> destblocks;
	for(int iblock=start; iblock<end; iblock++)
		destblocks.push_back( &coll[vInfo[iblock].blockID] );
	
	UpdateBlocks update;
	update.blocks = destblocks;
	update.myvelblocks = alldestblocks[rank];
	
	tbb::parallel_for(blocked_range<int>(0, work2node[rank]), update, auto_partitioner());
}

void I2D_VelocitySolverMPI_Mani::compute_velocity()
{
	Profiler profiler;
	
	assert(grid_ptr != NULL);
	assert(comm_rank == 0);
	
	int byte_size = -1;
	MPI_Type_size(MPI_BYTE, &byte_size);
	assert(byte_size  == sizeof(char));
	
	assert(srcparticles == NULL);
	assert(nsource_particles == 0);
	
	//estimate the partition of work:
	vector<BlockInfo> vInfo = grid_ptr->getBlocksInfo();
	ntotaldest_blocks = vInfo.size();
	
	const int slaves = comm_size - 1;
	int master_work = -1, slave_homogenous_work = -1, slave_remainder_work = -1;
	
	if(slaves > 0)
	{
		master_work =  (int)(3./7.*vInfo.size()/comm_size);
		slave_homogenous_work = (vInfo.size()-master_work)/slaves;
		slave_remainder_work = vInfo.size() - master_work - slaves*slave_homogenous_work;
	}
	else
	{
		master_work = vInfo.size();
		slave_homogenous_work = slave_remainder_work = 0;
	}
	
	work2node = vector<int>(comm_size);
	work2node[0] = master_work;
	for(int i=1; i<comm_size; i++)
		work2node[i] = slave_homogenous_work;
	
	for(int i=0; i<slave_remainder_work; i++)
		work2node[i+1]++;
	
	print(work2node, "decomposition");	
	workIDstart2node = vector<int> (comm_size);
	workIDstart2node[0]=0;
	for(int i=1; i<comm_size; i++)
		workIDstart2node[i] = workIDstart2node[i-1] + work2node[i-1];
	
	profiler.push_start("MPI-FMM");
	profiler.push_start("count source particles");
	_count_sourceparticles();
	profiler.pop_stop();
	
	if(nsource_particles == 0)
	{
		const BlockCollection<B>& coll = grid_ptr->getBlockCollection();
		
		ClearVelocity clear;
		block_processing.process(vInfo, coll ,clear);
		
		return;
	}
	
	profiler.push_start("broad cast header");
	_bcast_header();
	profiler.pop_stop();
	
	profiler.push_start("collect source particles");
	_collect_sourceparticles();
	profiler.pop_stop();

	profiler.push_start("broadcast source particles");
	_bcast_sourcedata();
	profiler.pop_stop();
	
	profiler.push_start("the master part");
	_compute();
	_master_update_blocks(0);
	profiler.pop_stop();//("compute and collect results");
	
	profiler.push_start("compute and collect results");
	_collect_results();
	profiler.pop_stop();
	
	stepid++;
	profiler.pop_stop();
	profiler.printSummary();
}
#endif

