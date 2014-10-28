/*
 *  I2D_CoreFMM_AggressiveVel.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#include <tbb/parallel_sort.h>

#include "I2D_CoreFMM_AggressiveVel.h"

extern  double _THETA;
#define _FMMSILENT

#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <iomanip>

#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGEnvironment.h"
#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGProfiler.h"

#include "mani-fmm2d/VortexExpansions.h"
#include "mani-fmm2d/hcfmm_box.h"
#include "mani-fmm2d/hcfmm_boxBuilder_serial.h"
#ifndef _MRAG_TBB
#include "mani-fmm2d/hcfmm_evaluator_serial.h"
#else
#include "mani-fmm2d/hcfmm_evaluator_tbb.h"
#endif

struct AggressiveExpansion:  _VortexExpansions<VelocitySourceParticle, _ORDER_>
{
	const static int order = _ORDER_;
	typedef VelocitySourceParticle ParticleType;
	
	void evaluateExpansions (Real *location, VelocityRHS *out_RHS)
	{		
		const ExpansionsValueType rp = ExpansionsValueType(location[0],location[1])-ExpansionsValueType(this->Center[0],this->Center[1]);
				
		std::complex<double> csum = std::complex<double>(0,0);
		std::complex<double> prod = std::complex<double>(1,0);
		
#pragma unroll
		for (int n=0;n<_ORDER_;++n)
		{
			csum+=(prod*(std::complex<double>)this->values[0][n]);
			prod/=(std::complex<double>)rp; 
		}
		
		csum*= -ExpansionsValueType(0,1)/rp;
		
		out_RHS->x[0]+=csum.real();
		out_RHS->x[1]-=csum.imag();
	}
};

typedef HCFMM::Box<AggressiveExpansion,_FMM_MAX_LEVEL_> tBox;
typedef HCFMM::boxBuilder_serial<AggressiveExpansion,_FMM_MAX_LEVEL_> tBoxBuilder;

class Measurements {
public:
	Measurements (int num_target_blocks) :
		num_direct_evals (num_target_blocks), num_indirect_evals (num_target_blocks) {
	}

	std::vector <double> num_direct_evals, num_indirect_evals;
};

struct VelocityEvaluator
{
	BlockInfo * destblocks;
	
	tBox * const rootNode;
	const Real inv_scaling;
	
	Measurements& measurements;

	VelocityEvaluator(tBox * rootNode, const Real inv_scaling, Measurements& _measurements) :
	rootNode(rootNode), inv_scaling(inv_scaling), destblocks(NULL), measurements (_measurements)
	{
	}
	
	VelocityEvaluator(const VelocityEvaluator& c) :
	rootNode(c.rootNode), inv_scaling(c.inv_scaling), destblocks(c.destblocks), measurements (c.measurements)
	{
	}
	
	//we want to be able to distinguish between far(0) close(1) and very close(2)
	bool _isclose_box(tBox* sourceBox, const BlockInfo& info, tBox::Btype _theta) const
	{
		const Real h_block = pow(0.5, info.level);
		
		const Real blockCenter[2] = {
			(0.5 + info.index[0])*h_block, 
			(0.5 + info.index[1])*h_block
		};
		
		const Real b2b_dist = sqrt(pow(sourceBox->expansions.Center[0]-blockCenter[0], 2) +
								   pow(sourceBox->expansions.Center[1]-blockCenter[1], 2));
		
		const Real radius_source = sqrt(pow(sourceBox->h[0],2) + 
										pow(sourceBox->h[1],2));
		
		const Real radius_target = h_block*sqrt(2.0);
		const Real denom = max(b2b_dist-radius_target, std::numeric_limits<Real>::epsilon());
		
		return (radius_source>=_theta*denom);
	}		
	
	bool _is_intersecting(tBox& box, const BlockInfo& info) const
	{			
		Real min_pos[2], max_pos[2];
		info.pos(min_pos, 0,0);
		info.pos(max_pos, B::sizeX-1, B::sizeY-1);
		
		const Real intersection[2] = {
			min(max_pos[0], (Real)(box.center[0] + box.h[0]*0.5)) - max(min_pos[0], (Real)(box.center[0] - box.h[0]*0.5)),
			min(max_pos[1], (Real)(box.center[1] + box.h[1]*0.5)) - max(min_pos[1], (Real)(box.center[1] - box.h[1]*0.5))
		};
		
		return intersection[0]>=0 && intersection[1]>=0;
	}
	
	//LOOPC-STYLE:
	void operator()(blocked_range<int> range) const
	{
		for(int iblock=range.begin(); iblock<range.end(); iblock++)
		{
			
			const BlockInfo info = destblocks[iblock];
			assert(destblocks[iblock].ptrBlock != NULL);
			VelocityBlock& my_b = *(VelocityBlock*)destblocks[iblock].ptrBlock;
			
			HCFMM::BoxIterator<tBox,tbb::scalable_allocator> it1(rootNode);
			
			tBox* current;
			bool canRemove;
			
			//clean velocity field:
			my_b.clear();
			
			while(it1!=NULL && (it1->nParticles>0))
			{
				canRemove=true;
				current=it1;
				
				const bool is_close = _isclose_box(current, info, _THETA);
				
				if(!is_close && current->parent!=NULL)
				{
					const bool parent_is_close = _isclose_box(current->parent, info, _THETA);
					
					if(parent_is_close) //if we were already well separated from parent, we should skipt this.

							for(int iy=0; iy<B::sizeY; iy++)
								for(int ix=0; ix<B::sizeX; ix++)
								{
									Real target_pos[2] = {0,0};
									info.pos(target_pos, ix, iy);
									
									if(current->parent==NULL || !(ws_barnes_hut(current->parent, target_pos, _THETA)))
									{
										VelocityRHS rhs;
										
										it1->expansions.evaluateExpansions(target_pos, &rhs);
										
										my_b.u[0][iy][ix] += rhs.x[0];
										my_b.u[1][iy][ix] += rhs.x[1];

										measurements.num_indirect_evals[iblock] += 1;
									}
								}
				}
				else  
				{
					//treat each point separately
						for(int iy=0; iy<B::sizeY; iy++)
							for(int ix=0; ix<B::sizeX; ix++)
							{
								Real target_pos[2] = {0,0};
								info.pos(target_pos,ix,iy);
								
								if(current->parent==NULL || !(ws_barnes_hut(current->parent, target_pos, _THETA))) 
								{	
									//when we were already well separated from the parent, we should not further interact.
									if (ws_barnes_hut(current, target_pos, _THETA))
									{
										measurements.num_indirect_evals[iblock] += 1;

										VelocityRHS rhs;
										
										it1->expansions.evaluateExpansions(target_pos, &rhs);
										
										my_b.u[0][iy][ix] += rhs.x[0];
										my_b.u[1][iy][ix] += rhs.x[1];
									}
									else
									{
										if(!it1->isleaf) 
											//Case 2: its not a leaf so we further traverse into the tree (summing up the children and bailing out)
											canRemove=false;
										else  //its close and a leaf ->interactDirectly with particles
										{
											Real u[2] = {0,0};
											
											const int nof_sources = current->nParticles;
											const VelocitySourceParticle * const p = current->vparticles;
											
											measurements.num_direct_evals[iblock] += nof_sources;

											if (_is_intersecting(*current, info)) //tBox current and this block are intersecting
												for (int i=0;i<nof_sources;++i)
												{
													Real r[2] = {
														target_pos[0] - p[i].x[0],
														target_pos[1] - p[i].x[1]
													};
													
													const Real distance_2 = r[0]*r[0] + r[1]*r[1];
													
													if (distance_2==0) continue;
													
													const Real factor = 1/distance_2;
													
													u[0] -= factor*p[i].w[0]*r[1];
													u[1] += factor*p[i].w[0]*r[0];
												}												
											else
												for (int i=0;i<nof_sources;++i)
												{
													Real r[2] = {
														target_pos[0] - p[i].x[0],
														target_pos[1] - p[i].x[1]
													};
													
													const Real factor = 1/(r[0]*r[0] + r[1]*r[1]);
													
													u[0] -= factor*p[i].w[0]*r[1];
													u[1] += factor*p[i].w[0]*r[0];
												}													
											
											assert(!isnan(u[0]));
											assert(!isnan(u[1]));
											
											my_b.u[0][iy][ix] += u[0];
											my_b.u[1][iy][ix] += u[1];
										}
									} 
								}
							}
				}
				
				if(canRemove)
					it1.advanceRemove();
				else
					it1++;
			}
			
			//multiply by scaling factor
			{
				Real * const ue = (Real *)my_b.u[0];
				Real * const ve = (Real *)my_b.u[1];
				const Real scale = 1./(2.0*M_PI)*inv_scaling;
				static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
				for(int i=0; i<n; i++) 
				{
					ue[i] *= scale;
					ve[i] *= scale;
				}
			}
#ifndef _FMMSILENT
			printf("Done with %d %d %d l=%d\n", info.index[0], info.index[1], info.index[2], info.level);
#endif
		}
		
	}
};

void I2D_CoreFMM_AggressiveVel::solve(const Real theta, const Real inv_scaling, BlockInfo * dest, const int nblocks, VelocitySourceParticle * srcparticles, const int nparticles)
{
	Profiler profiler;
	
	_THETA = theta;
	
	tBox * rootBox=new tBox;
	tick_count start_before_tree = tick_count::now ();
	profiler.push_start("tree");
	tBoxBuilder::buildBoxes(srcparticles, nparticles, rootBox);
	profiler.pop_stop();
	
	tick_count start_before_expansions = tick_count::now ();
	profiler.push_start("expansions");
	tBoxBuilder::generateExpansions(rootBox);
	profiler.pop_stop();
	
	tick_count start_before_evaluations = tick_count::now ();
	profiler.push_start("evaluations");
	Measurements measurements (nblocks);
	VelocityEvaluator evaluator(rootBox, inv_scaling, measurements);
	evaluator.destblocks = dest;
	tbb::parallel_for(blocked_range<int>(0, nblocks), evaluator, auto_partitioner());
	profiler.pop_stop();
	tick_count end = tick_count::now ();

	if (timestamp++ % 5 == 0) profiler.printSummary();
	
	if (b_verbose) {
	
		const double tree_wallclock_time = (start_before_expansions - start_before_tree).seconds();
		const double evaluations_wallclock_time = (end - start_before_evaluations).seconds();
		const double expansions_wallclock_time = (start_before_evaluations - start_before_expansions).seconds();

		const double num_effective_interactions = (double)nparticles*nblocks*_BLOCKSIZE_*_BLOCKSIZE_;

		double effective_direct_gflops, effective_indirect_gflops, effective_total_gflops;
		double actual_direct_gflops, actual_indirect_gflops, actual_total_gflops;

		double num_direct_evals, num_indirect_evals;
		num_direct_evals = num_indirect_evals = 0;
		for (int i=0; i<nblocks; ++i) {
			num_direct_evals += (double)measurements.num_direct_evals[i];
			num_indirect_evals += (double)measurements.num_indirect_evals[i] ;
		}

		effective_direct_gflops = num_direct_evals * 11./1e9;
		actual_direct_gflops = num_direct_evals * 15./1e9;
		effective_indirect_gflops = num_indirect_evals * (16 + 14*_ORDER_) * 1./1e9;
		actual_indirect_gflops = num_indirect_evals * (20 + 14*_ORDER_) * 1./1e9;

		effective_total_gflops = effective_direct_gflops + effective_indirect_gflops;
		actual_total_gflops = actual_direct_gflops + actual_indirect_gflops;

		//Write data for break even plot
		std::fstream file;


		file.open ("measurements.txt",  std::fstream::app | std::fstream::out);
		file << setprecision (6) << scientific << num_effective_interactions << "\t" << effective_total_gflops << "\t" \
				<< actual_total_gflops << "\t" << evaluations_wallclock_time << "\t" << \
				expansions_wallclock_time << "\t" << tree_wallclock_time << "\n";
		file.close ();


	}
	delete rootBox;
}
