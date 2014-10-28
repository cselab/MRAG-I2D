/*
 *  hcfmm_evaluator_serial.h
 *  hcfmm
 *
 *  Created by Manfred on 1/27/09.
 *  Copyright 2009 ETHZ. All rights reserved.
 *
 */

#ifndef _HCFMM_EVAL_TBB_
#define _HCFMM_EVAL_TBB_
#include <cassert>
#include "hcfmm_evaluator_serial.h"
#include "well_separated.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"



namespace HCFMM{
	
template <typename tExpansions, typename TargetPoint,  int _maxlevel>
class Evaluator_tbb: public Evaluator_serial<tExpansions,TargetPoint,_maxlevel>
	{
	
	public:
		typedef Evaluator_serial<tExpansions,TargetPoint,_maxlevel> MotherClass;
		typedef Evaluator_tbb<tExpansions,TargetPoint,_maxlevel> SelfClass;
		typedef typename tExpansions::ParticleType Particle;
		typedef Box<tExpansions,_maxlevel> tBox;
		typedef typename TargetPoint::RHSType tRHS;
		
		//Constructor
		Evaluator_tbb(tBox* in_box):MotherClass(in_box)
		{
		}
		
		//overloading 	//overloading these now: (with tbb-functionality)

		void evaluate(TargetPoint* inTargets, int nTargets);

	protected:

		struct tbb_interact{
			typedef typename TargetPoint::RHSType tRHS;
			TargetPoint* const inTargets;
			SelfClass* const thisEvaluator;
			tRHS* const local_RHS;
			typename TargetPoint::BaseType theta;

			
#ifndef _OPT_Split
			//operator (needed for tbb)
			void operator() ( const tbb::blocked_range<size_t> &r ) const {
				for (size_t i=r.begin(); i!=r.end();++i)
				{
					//ok: this might be problematic: due to memory access of different threads to the same region.???
					thisEvaluator->_interact(&(inTargets[i]),&(inTargets[i].RHS),theta);
				}
				
			}
#else
			void operator() ( const tbb::blocked_range<size_t> &r ) const {
				
				for (size_t i=r.begin(); i!=r.end();++i)
				{
					//ok: this might be problematic: due to memory access of different threads to the same region.???
					thisEvaluator->_interact(&(inTargets[i]),&(local_RHS[i]),theta);
				}
				
				
			}
#endif
			
			//constructor:
			tbb_interact(TargetPoint* in_Targets, SelfClass* in_thisEvaluator, typename Particle::BaseType in_theta, tRHS* outRHS=NULL):
			inTargets(in_Targets),thisEvaluator(in_thisEvaluator),local_RHS(outRHS),theta(in_theta)
			{}
			

			
		};
		
#ifdef _OPT_Split
		struct tbb_writeback{
			typedef typename TargetPoint::RHSType tRHS;
			tRHS* const local_RHS;
            TargetPoint* local_Targets;
			
			tbb_writeback(TargetPoint* in_Targets, tRHS* in_RHS):
			local_RHS(in_RHS),local_Targets(in_Targets)
			{}
			
			void operator() ( const tbb::blocked_range<size_t> &r ) const {
				for (size_t i=r.begin(); i!=r.end();++i)
				{
					local_Targets[i].RHS=local_RHS[i];
				}

			}

		};
#endif
		
		
		struct tbb_interact_recursive{
			typedef typename TargetPoint::RHSType tRHS;
			TargetPoint* const inTargets;
			SelfClass* const thisEvaluator;
			tRHS* const local_RHS;
			typename TargetPoint::BaseType theta;

			
			//operator (needed for tbb)
			void operator() ( const tbb::blocked_range<size_t> &r ) const {
				for (size_t i=r.begin(); i!=r.end();++i)
				{
					tRHS tmpRHS;
					tmpRHS=thisEvaluator->traverseRecursive(&inTargets[i],thisEvaluator->rootBox,theta);
					inTargets[i].RHS=tmpRHS;
				}
				
			}
			
			//constructor:
			tbb_interact_recursive(TargetPoint* in_Targets, SelfClass* in_thisEvaluator, typename TargetPoint::BaseType in_theta, tRHS* outRHS=NULL):
			inTargets(in_Targets),thisEvaluator(in_thisEvaluator),local_RHS(outRHS),theta(in_theta)
			{}
			
			
			
		};		
		
		
		struct tbb_interact_loopchanged{
			typedef typename TargetPoint::RHSType tRHS;
			TargetPoint* const inTargets;
			SelfClass* const thisEvaluator;
			tRHS* const local_RHS;
			typename TargetPoint::BaseType theta;
			
			//operator (needed for tbb)
			void operator() ( const tbb::blocked_range<size_t> &r ) const {
				
				BoxIterator<tBox> it1(thisEvaluator->rootBox);
				bool canRemove;

//for j=0;j<D;j+=chunk				
				
				while(it1!=NULL)
				{
					
					canRemove=true; //we assume that we can remove the current branch
					tBox* current=it1;

					for (size_t i=r.begin(); i!=r.end();++i)
						//for i=i*h
					{
						if(current->parent==NULL || !(ws_barnes_hut(current->parent, inTargets[i].x, theta))) //when we were already well separated from the parent, we should not further interact.
						{
						//Case One: Well Separated: Use the Multipoles to calculate rhs and continue
						if (ws_barnes_hut(current, inTargets[i].x, theta))
						{
							it1->expansions.evaluateExpansions(&(inTargets[i].x[0]),&(inTargets[i].RHS));
						}
						else// if not Well Separated
						{
							if(!it1->isleaf) //Case 2: its not a leaf so we further traverse into the tree (summing up the children and bailing out)
							{
								canRemove=false; //need to traverse further into tree.
							}
							else  //its close and a leaf ->interactDirectly with particles
							{
								thisEvaluator->_interactDirect(&(inTargets[i]),it1,&(inTargets[i].RHS));
								
							}
						}
						
						}
					} //Particle-loop.
					
					if(canRemove)
					{
						it1.advanceRemove();
					}
					else{
						it1++;
					}
					
				}  //While-Loop over boxes.
			}
			
			//constructor:
			tbb_interact_loopchanged(TargetPoint* in_Targets, SelfClass* in_thisEvaluator, typename Particle::BaseType in_theta=0.75 , tRHS* outRHS=NULL):
			inTargets(in_Targets),thisEvaluator(in_thisEvaluator),local_RHS(outRHS),theta(in_theta)
			{}
			
			
			
		};
		
			
		//tbboperator for TBOX
		struct tbb_evaluate_tbox{
			typedef typename TargetPoint::RHSType tRHS;
			SelfClass* const thisEvaluator;
			typename TargetPoint::BaseType theta;
			recursiveListBuilder<tBox>* const iList; 
			
			
			//operator (needed for tbb)
			void operator() ( const tbb::blocked_range<size_t> &r ) const {
				
				typename recursiveListBuilder<tBox>::tContainer::iterator it2;

				
				for (size_t iL=r.begin(); iL!=r.end();++iL)
				{
					
					for(it2=iList->InteractionList[iL].farBoxList.begin();it2!=iList->InteractionList[iL].farBoxList.end();++it2)
					{
						for(int i=0;i<iList->InteractionList[iL].Leaf->nParticles;++i)
						{
							(*it2)->expansions.evaluateExpansions(iList->InteractionList[iL].Leaf->vparticles[i].x,&(iList->InteractionList[iL].Leaf->vparticles[i].RHS));
						}
					}
					
					for (it2=iList->InteractionList[iL].closeBoxList.begin();it2!=iList->InteractionList[iL].closeBoxList.end();++it2)
					{
						for(int i=0;i<iList->InteractionList[iL].Leaf->nParticles;++i)
						{
							if(ws_barnes_hut(*it2, iList->InteractionList[iL].Leaf->vparticles[i].x, theta))
							{
								(*it2)->expansions.evaluateExpansions(iList->InteractionList[iL].Leaf->vparticles[i].x,&(iList->InteractionList[iL].Leaf->vparticles[i].RHS));
							}
							else
							{
								thisEvaluator->_interactDirect(&(iList->InteractionList[iL].Leaf->vparticles[i]),(*it2),&(iList->InteractionList[iL].Leaf->vparticles[i].RHS));
							}
						}
					}
					
					for (it2=iList->InteractionList[iL].directBoxList.begin();it2!=iList->InteractionList[iL].directBoxList.end();++it2)
					{
						for(int i=0;i<iList->InteractionList[iL].Leaf->nParticles;++i)
						{
							thisEvaluator->_interactDirect(&(iList->InteractionList[iL].Leaf->vparticles[i]),(*it2),&(iList->InteractionList[iL].Leaf->vparticles[i].RHS));
						}
						
					}
					
					
				}
				
			}
			
			//constructor:
			tbb_evaluate_tbox(recursiveListBuilder<tBox>* in_iList, SelfClass* in_thisEvaluator, typename Particle::BaseType in_theta):
			thisEvaluator(in_thisEvaluator),theta(in_theta),iList(in_iList)
			{}
			
			
			
		};
		
		
        //void _interact(TargetPoint* inTargets, tRHS* rhs_out);
		//inline void _interactDirect(TargetPoint* inTargets, tBox* p_Box, tRHS* rhs_out);
		//inline void _interactBox(TargetPoint* inTargets, tBox* p_Box, tRHS* rhs_out);

	};
	
template <typename tExpansions, typename TargetPoint,  int _maxlevel>
void Evaluator_tbb<tExpansions,TargetPoint,_maxlevel>::evaluate(TargetPoint* inTargets, int nTargets)
{

	if(this->rootBox->nParticles<1)
		return;
	
#ifndef _TUSEBOX
#ifndef _CHUNKSIZE
#ifndef _TLOOPC
#ifndef _TRECURSIVE	
#ifndef _OPT_Split
	typename Particle::BaseType theta=_THETA;
#ifndef _FMMSILENT
	std::cout << "YOU ARE CALLING: Evaluator_tbb:evaluate, theta is: " << _THETA <<std::endl;
#endif

	//for (int i=0; i<nTargets;++i) //this is probably the most expensive Loop
//	{
//		int k=nTargets/10;
//		(k==0)?k=1:k=k;
//		if(i%k==0) std::cout << "Calculating... Passed " << i << " of " << nTargets <<std::endl;
//		_interact(&(inTargets[i]),&(inTargets[i].RHS));
//	}
	tbb::parallel_for(tbb::blocked_range<size_t>(0,nTargets), tbb_interact(inTargets,this,theta), tbb::auto_partitioner());
#else
#ifndef _FMMSILENT
	std::cout << "YOU ARE CALLING: Evaluator_tbb:evaluate --> using a copy for the RHSs" <<std::endl;
#endif

	typedef typename TargetPoint::RHSType tRHS;
	tRHS* tmpRHS=new tRHS[nTargets];
	tbb::parallel_for(tbb::blocked_range<size_t>(0,nTargets), tbb_interact(inTargets,this,tmpRHS), tbb::auto_partitioner());
	
	//also parallelize this.
//	for(int i=0;i<nTargets;++i)
//	{
//		inTargets[i].RHS=tmpRHS[i];
//	}

	tbb::parallel_for(tbb::blocked_range<size_t>(0,nTargets), tbb_writeback(inTargets,tmpRHS), tbb::auto_partitioner());

	delete[] tmpRHS;
	
	
#endif
#else
	typename Particle::BaseType theta=_THETA;
#ifndef _FMMSILENT
	std::cout << "YOU ARE CALLING: Evaluator_tbb:evaluate RECURSIVE, theta is: " << theta <<std::endl;
#endif
   	tbb::parallel_for(tbb::blocked_range<size_t>(0,nTargets), tbb_interact_recursive(inTargets,this,theta), tbb::auto_partitioner());

	
#endif
#else
	
	typename Particle::BaseType theta=_THETA;
#ifndef _FMMSILENT
	std::cout << "YOU ARE CALLING: Evaluator_tbb:evaluate LOOP_CHANGED, theta=" << theta <<std::endl;
#endif

#ifndef _TBBGRAINSIZE
   	tbb::parallel_for(tbb::blocked_range<size_t>(0,nTargets), tbb_interact_loopchanged(inTargets,this,theta), tbb::auto_partitioner());
#else
#ifndef _FMMSILENT
	std::cout << "GRAINSIZE Is: " << _TBBGRAINSIZE <<std::endl;
#endif
	tbb::parallel_for(tbb::blocked_range<size_t>(0,nTargets,_TBBGRAINSIZE), tbb_interact_loopchanged(inTargets,this,theta));
#endif
	
#endif
	
	
#else
	typename Particle::BaseType theta=_THETA;
#ifndef _FMMSILENT
	std::cout << "YOU ARE CALLING: Evaluator_tbb:evaluate LOOP_CHANGED theta=" << theta <<std::endl;
	std::cout << "CHUNKSIZE is: " << _CHUNKSIZE << std::endl;
#endif
	tbb::parallel_for(tbb::blocked_range<size_t>(0,nTargets), tbb_interact_loopchanged_chunked<_CHUNKSIZE>(inTargets,this,theta), tbb::auto_partitioner());

	
#endif
#else
	typename Particle::BaseType theta=_THETA;
#ifndef _FMMSILENT
	std::cout << "You are calling: Evaluator_tbb:evaluate USE BOX for TARGETS TOO, theta is: " << theta <<std::endl;
#endif
	recursiveListBuilder<tBox> iList(this->rootBox,this->rootBox,theta);
	tbb::parallel_for(tbb::blocked_range<size_t>(0,iList.InteractionList.size()),tbb_evaluate_tbox(&iList,this,theta),tbb::auto_partitioner());
	
	
#endif
	
}
	



	

}

#endif

