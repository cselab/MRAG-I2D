/*
 *  hcfmm_evaluator_serial.h
 *  hcfmm
 *
 *  Created by Manfred on 1/27/09.
 *  Copyright 2009 ETHZ. All rights reserved.
 *
 */

#ifndef _HCFMM_EVAL_SER_
#define _HCFMM_EVAL_SER_

#include <cassert>
#include "hcfmm_evaluator.h"
#include "well_separated.h"
//for timings:
#include <tbb/tick_count.h>





namespace HCFMM{
	
	
template <typename tExpansions, typename TargetPoint,  int _maxlevel>
class Evaluator_serial: public Evaluator<tExpansions,TargetPoint,_maxlevel>
	{
	
	public:
		typedef Evaluator<tExpansions,TargetPoint,_maxlevel> MotherClass;
		typedef Evaluator_serial<tExpansions,TargetPoint,_maxlevel> SelfClass;
		typedef typename tExpansions::ParticleType Particle;
		typedef Box<tExpansions,_maxlevel> tBox;
		typedef typename TargetPoint::RHSType tRHS;
		
		Evaluator_serial(tBox* in_box):MotherClass(in_box){}
		void evaluate(TargetPoint* inTargets, int nTargets, const typename TargetPoint::BaseType theta=_THETA);
		
		
	//private:
		void _interact(TargetPoint* inTargets, tRHS* rhs_out,const typename Particle::BaseType theta) const; //, size_t &cnt_box, size_t& cnt_direct
		inline void _interactDirect(TargetPoint* inTargets, tBox* p_Box, tRHS* rhs_out) const;
		inline void _interactBox(TargetPoint* inTargets, tBox* p_Box, tRHS* rhs_out) const;
		tRHS traverseRecursive(TargetPoint* inTarget, tBox* node, const typename tExpansions::Btype theta) const;

	};
	
	

template <typename tBox>
struct recursiveListBuilder
	{
		
	typedef tBox* p2Box;
	typedef std::list<p2Box> tContainer;	
	static const int _nc=1<<tBox::Particle::dim;
		
		
	struct DirectBoxComparison
		{
			bool operator()(const p2Box & a, const p2Box& b)
			{
				return (a->vparticles < b->vparticles);
			}
		};


		struct InteractionListElement {
			p2Box Leaf; //pointer to leaf box
			tContainer farBoxList;
			tContainer nearBoxList;
			tContainer directBoxList;
			tContainer closeBoxList;
			
//			InteractionListElement( const InteractionListElement& rhs):
//			Leaf(rhs.Leaf),farBoxList(rhs.farBoxList),nearBoxList(rhs.nearBoxList),directBoxList(rhs.directBoxList)
//			{}
//			
//			InteractionListElement():
//			Leaf(NULL),farBoxList(),nearBoxList(),directBoxList(),closeBoxList()
//			{}
			
			bool operator==(const InteractionListElement &other) const {
				return this->Leaf==other.Leaf;
			}
			
			bool operator<(const InteractionListElement &other) const {
				return (this->Leaf<other.Leaf);
			}
			
		};
		std::vector<InteractionListElement> InteractionList;
		
		struct binpred {
			
			bool operator()(const InteractionListElement& lhs, const InteractionListElement& rhs)
			{
				return(lhs.Leaf==rhs.Leaf);
			}
		};
		
	typedef typename std::vector<InteractionListElement>::iterator iterator;

		
	bool recursiveBuild(tBox *target_rootBox, InteractionListElement iElement, typename tBox::Particle::BaseType theta)
	{
	
		//std::cout << "Size of incoming nearBoxList= " <<iElement.nearBoxList.size() <<std::endl; 
		
		tBox* current;
		tBox* child;
		int isclose;
		//if we have reached the leaf of the target-Tree:
		//go through all near_list and traverse down to source-leafs and put into direct or far, at the end, near must be empty
		if(target_rootBox->isleaf)
		{
			while(!iElement.nearBoxList.empty())
			{
				current=iElement.nearBoxList.front();
				isclose=isclose_box_barnes_hut(current, target_rootBox, theta);
				if(isclose==0)
				{
					iElement.farBoxList.push_back(current);  //it is far now (because we are in a different targetBox)
				}
				else {
					if(current->isleaf)
					{
						if(isclose==2) //it is very close (all have to do direct interactions)
						{
							iElement.directBoxList.push_back(current);  //it is close and source is also leaf -> direct
						}
						else {
							iElement.closeBoxList.push_back(current);
						}
					}
					else {
						for(int kb=0;kb<_nc;++kb)   //it is near but, not a leaf -> put children in near list.
						{
							if(current->children[kb]!=NULL)
							{
							  iElement.nearBoxList.push_back(current->children[kb]);
							}
						}
					}
				}
				iElement.nearBoxList.pop_front();
			}
			iElement.Leaf=target_rootBox;
			
			//reprocess->close or direct
			typename tContainer::iterator it_close=iElement.closeBoxList.begin();
			typename tContainer::iterator tmp_it;
			int cnt_corr(0);
			while(it_close!=iElement.closeBoxList.end())
			{
				isclose=isclose_box_barnes_hut(iElement.Leaf,*it_close, theta);
				
				if(isclose==2)
				{
					iElement.directBoxList.push_back(*it_close);
					tmp_it=it_close;
					iElement.closeBoxList.erase(tmp_it);
					cnt_corr++;
				}
				
				it_close++;
			}
			//std::cout << "reprocessed close-> corrected " << cnt_corr << " size direct: " <<iElement.directBoxList.size() << " size close: " << iElement.closeBoxList.size()  <<std::endl;
			//put it in the Leaf-List
			InteractionList.push_back(iElement);
			return true;
		}
		
		//if we have not reached the leaf of the target-Tree:
		//check things in near-list, traverse one level down in the source-Tree.
		//call this function on all children handing over the current iElement.
		else
		{
            int cSize=iElement.nearBoxList.size();
			int cnt=0;
			typename tContainer::iterator tmpit;
			typename tContainer::iterator it1=iElement.nearBoxList.begin();
			while(cnt<cSize && it1!=iElement.nearBoxList.end()) //loop once through near list.
			{
				isclose=isclose_box_barnes_hut(*it1, target_rootBox, theta);
				if(isclose==0)
				{
					iElement.farBoxList.push_back(*it1);
				}
				else {
					if ((*it1)->isleaf) {
						if(isclose==2)
						{
							iElement.directBoxList.push_back(*it1);
						}
						else
						{
							iElement.closeBoxList.push_back(*it1);
						}
					}
					else
					{
						for (int kb=0;kb<_nc;++kb)
						{
						  if ((*it1)->children[kb]!=NULL)
						  {
						    iElement.nearBoxList.push_back((*it1)->children[kb]);
						  }
					    }
					}
				}
				
				tmpit=it1;
				it1++;
				iElement.nearBoxList.erase(tmpit);
				cnt++;
			}
			
			while(it1!=iElement.nearBoxList.end())
			{
				isclose=isclose_box_barnes_hut(*it1, target_rootBox, theta);
				if(!isclose)
				{
					iElement.farBoxList.push_back(*it1);
					tmpit=it1;
					it1++;
					iElement.nearBoxList.erase(tmpit);
				}
				else {
					if ((*it1)->isleaf) {
						if (isclose==2)
						{
						   iElement.directBoxList.push_back(*it1);
						}
						else
						{
							iElement.closeBoxList.push_back(*it1);
						}
						
						tmpit=it1;
						it1++;
						iElement.nearBoxList.erase(tmpit);
					}
					else{
						it1++;
					}
				}
				
			}
			
			bool tmpbol=true;
			bool tmpbol2=true;
			for (int kb=0; kb<_nc;++kb)
			{
				if(target_rootBox->children[kb]!=NULL)
				{
					tmpbol2=recursiveBuild(target_rootBox->children[kb],iElement,theta);
					tmpbol=(tmpbol && tmpbol2);
				}
			}
			return tmpbol;
			
		}
		
	}
		
		recursiveListBuilder(tBox *target_rootBox, tBox* source_rootBox, typename tBox::Particle::BaseType theta)
		{
			InteractionList.clear();
			InteractionListElement startElement;
			startElement.farBoxList.clear();
			startElement.nearBoxList.clear();
			startElement.directBoxList.clear();
			startElement.Leaf=NULL;
			startElement.nearBoxList.push_back(source_rootBox);
			
			recursiveBuild(target_rootBox,startElement,theta);

			size_t cbox(0),cdir(0);
			iterator it1=InteractionList.begin();
			typename tContainer::iterator it2;
			while (it1!=InteractionList.end())
			{
				cbox+=it1->Leaf->nParticles*it1->farBoxList.size();
				for (it2=it1->directBoxList.begin();it2!=it1->directBoxList.end();++it2)
				{
					cdir+=(it1->Leaf->nParticles)*(*it2)->nParticles;
				}
				for (it2=it1->closeBoxList.begin();it2!=it1->closeBoxList.end();++it2)
				{
					cdir+=(it1->Leaf->nParticles)*(*it2)->nParticles;
				}
				
#ifdef _DOSORT
				//std::sort(it1->directBoxList.begin(),it1->directBoxList.end(),DirectBoxComparison());
				//std::sort(it1->farBoxList.begin(),it1->farBoxList.end());
				it1->farBoxList.sort();
				it1->directBoxList.sort(DirectBoxComparison());
#endif
				
				it1++;
			}
			std::cout << "Planned Interactions: box: " << cbox << "direct: " <<cdir <<std::endl; 


			
		}
		
	};


#ifndef _TUSEBOX2
#ifndef _TUSEBOX
#ifndef _TLOOPC
#ifndef _TRECURSIVE
template <typename tExpansions, typename TargetPoint,  int _maxlevel>
void Evaluator_serial<tExpansions,TargetPoint,_maxlevel>::evaluate(TargetPoint* inTargets, int nTargets,const typename TargetPoint::BaseType theta)
{
	if(this->rootBox->nParticles<1)
		return;
	
	//size_t cnt_box(0);
	//size_t cnt_direct(0);
//	global_box=0;
//	global_direct=0;
//	global_self=0;

	for (int i=0; i<nTargets;++i) //this is probably the most expensive Loop
	{
//		int k=nTargets/10;
//		(k==0)?k=1:k=k;
//		if(i%k==0) std::cout << "Calculating... Passed " << i << " of " << nTargets <<std::endl;
		_interact(&(inTargets[i]),&(inTargets[i].RHS),theta);
	}
//	std::cout << "Standard evaluate: I had " << cnt_box << " box_interactions " <<  cnt_direct << " direct interactions" <<std::endl;
//	std::cout << "global counters: box: " <<global_box << " direct: " << global_direct  << " self: " << global_self<< std::endl;


}
#else
	template <typename tExpansions, typename TargetPoint,  int _maxlevel>
	void Evaluator_serial<tExpansions,TargetPoint,_maxlevel>::evaluate(TargetPoint* inTargets, int nTargets,const typename TargetPoint::BaseType theta)
	{
		if(this->rootBox->nParticles<1)
			return;
		
		tRHS res;
		for (int i=0; i<nTargets;++i) //this is probably the most expensive Loop
		{
//			int k=nTargets/10;
//			(k==0)?k=1:k=k;
//			if(i%k==0) std::cout << "Calculating... Passed " << i << " of " << nTargets <<std::endl;
			res=traverseRecursive(&inTargets[i],this->rootBox,theta);
			inTargets[i].RHS=res;
		}
	}
	
#endif
#else

	template <typename tExpansions, typename TargetPoint,  int _maxlevel>
	void Evaluator_serial<tExpansions,TargetPoint,_maxlevel>::evaluate(TargetPoint* inTargets, int nTargets, const typename TargetPoint::BaseType theta)
	{
		if(this->rootBox->nParticles<1)
			return;
		
		//global_box=0;
		//global_direct=0;
#ifndef _FMMSILENT
		std::cout << "Evaluating SERIAL, LOOP CHANGED -> theta is: " << theta <<std::endl;
#endif
		BoxIterator<tBox> it1(this->rootBox);
		tBox* current;
		bool canRemove;
#ifdef VERBOSE
		int cnt_box(0),cnt_direct(0);
#endif
		while(it1!=NULL)
		{
			canRemove=true; //we assume that we can remove the current branch
			current=it1;
			
			for (int i=0; i<nTargets;++i) //this is probably the most expensive Loop
			{
				if(current->parent==NULL || !(ws_barnes_hut(current->parent, inTargets[i].x, theta))) //when we were already well separated from the parent, we should not further interact.
				{
					if (ws_barnes_hut(current, inTargets[i].x, theta))
					{
						it1->expansions.evaluateExpansions(&(inTargets[i].x[0]),&(inTargets[i].RHS));
#ifdef VERBOSE
						cnt_box++;
#endif
					}
					else// if not Well Separated
					{
						if(!it1->isleaf) //Case 2: its not a leaf so we further traverse into the tree (summing up the children and bailing out)
						{
							canRemove=false;
						}
						else  //its close and a leaf ->interactDirectly with particles
						{
							this->_interactDirect(&(inTargets[i]),it1,&(inTargets[i].RHS));
#ifdef VERBOSE
							cnt_direct+=it1->nParticles;
#endif
						}
					}
				}
			}
			
			if(canRemove)
			{
				it1.advanceRemove();
			}
			else{
				it1++;
			}
			
		}
#ifdef #VERBOSE
		std::cout << "Evaluated with Tloopc-> box interactions: " <<cnt_box << " direct: " << cnt_direct <<std::endl;
#endif
//		std::cout << "global counters: box: " <<global_box << " direct: " << global_direct << std::endl;
	}
	
	
	
#endif
#else
	template <typename tExpansions, typename TargetPoint,  int _maxlevel>
	void Evaluator_serial<tExpansions,TargetPoint,_maxlevel>::evaluate(TargetPoint* inTargets, int nTargets, const typename TargetPoint::BaseType theta)
	{
		if(this->rootBox->nParticles<1)
			return;
//		global_box=0;
//		global_direct=0;
//		global_self=0;
#ifndef _FMMSILENT
		std::cout << "Evaluating SERIAL,USING BOX STRUCTURE FOR TARGETS , theta is: " << theta  <<std::endl;
#endif
#ifdef VERBOSE
		int cnt_direct(0);
		int cnt_box(0);
		int cnt_corrected(0);
#endif
		typedef tBox* p2Box;
		double tdur;
		tbb::tick_count tstart,tend;
		
		tstart=tbb::tick_count::now();
		recursiveListBuilder<tBox> iList(this->rootBox,this->rootBox,theta);
		tend=tbb::tick_count::now();
		tdur=(tend-tstart).seconds();
		std::cout << "recursive interaction list Build took: " << tdur <<std::endl;
		
		typename recursiveListBuilder<tBox>::iterator it1;
		typename recursiveListBuilder<tBox>::tContainer::iterator it2;
		
		for (it1=iList.InteractionList.begin();it1!=iList.InteractionList.end();++it1)
		{
			for (it2=it1->farBoxList.begin();it2!=it1->farBoxList.end();++it2)
			{
				for(int i=0;i<it1->Leaf->nParticles;++i)
				{
					(*it2)->expansions.evaluateExpansions(it1->Leaf->vparticles[i].x,&(it1->Leaf->vparticles[i].RHS));
					cnt_box++;
				}
			}
			
			for (it2=it1->closeBoxList.begin();it2!=it1->closeBoxList.end();++it2)
			{
				for(int i=0;i<it1->Leaf->nParticles;++i)
				{
					if(ws_barnes_hut(*it2, it1->Leaf->vparticles[i].x, theta))
					{
						
						(*it2)->expansions.evaluateExpansions(it1->Leaf->vparticles[i].x,&(it1->Leaf->vparticles[i].RHS));
						cnt_corrected++;
						cnt_box++;
					}
					else
					{
					_interactDirect(&(it1->Leaf->vparticles[i]),(*it2),&(it1->Leaf->vparticles[i].RHS));
						cnt_direct+=(*it2)->nParticles;
					}
				}
			}
			
			for (it2=it1->directBoxList.begin();it2!=it1->directBoxList.end();++it2)
			{
				for(int i=0;i<it1->Leaf->nParticles;++i)
				{
				   _interactDirect(&(it1->Leaf->vparticles[i]),(*it2),&(it1->Leaf->vparticles[i].RHS));
					cnt_direct+=(*it2)->nParticles;
				}
				
			}
			
	
			//cnt_direct+=it1->directBoxList.size()*it1->Leaf->nParticles+it1->closeBoxList.size()*it1->Leaf->nParticles;
			//cnt_box+=it1->farBoxList.size()*it1->Leaf->nParticles;
			
		}
		
		std::cout << "finished with evaluate-serial, Interactions with Boxes: " << cnt_box << " direct: " <<cnt_direct << " corrected: " <<cnt_corrected << std::endl;
//		std::cout << "global counters: box: " <<global_box << " direct: " << global_direct  << " self: " << global_self<< std::endl;

//		static const int nc=1<<Particle::dim;
//		
//		
//		BoxIterator<tBox> it1(this->rootBox);
//		size_t nLeaves=0;
//		p2Box current,currentLeaf;
//		
//		while (it1!=NULL)
//		{
//			if(it1->isleaf)
//			nLeaves++;
//			it1++;
//		}
//		std::cout << "found" << nLeaves << "Leafs" <<std::endl;
//
//		p2Box* leaves=new p2Box[nLeaves];
//		
//		size_t iBox=0;
//		BoxIterator<tBox> it2(this->rootBox);
//
//		while (it2!=NULL)
//		{
//			current=it2;
//			if(current->isleaf)
//			{
//				leaves[iBox]=current;
//				iBox++;
//
//			}
//			it2++;
//		}
//		assert(iBox==nLeaves);
//		
//				
//		for (iBox=0;iBox<nLeaves;++iBox)
//		{
//			BoxIterator<tBox> it3(this->rootBox);
//			currentLeaf=leaves[iBox];
//			while(it3!=NULL)
//			{
//				current=it3;
//			//check if BOXES are well separated:
//			if(ws_box_barnes_hut(current,currentLeaf,theta))
//			{
//				for(int i=0;i<currentLeaf->nParticles;++i)
//				{
//					if(!ws_barnes_hut(current, currentLeaf->vparticles[i].x, theta))
//					{
//						std::cout << "Something went wrong:  " <<std::endl;
//						exit(1);
//					}
//					current->expansions.evaluateExpansions(currentLeaf->vparticles[i].x,&(currentLeaf->vparticles[i].RHS));
//				}
//				cnt_box+=currentLeaf->nParticles;
//				it3.advanceRemove();
//			}
//			else
//			{
//				if(current->isleaf)
//				{
//				//Direct interactions with particles in here and brothers particles:
//				for (int i=0;i<currentLeaf->nParticles;++i)
//				{
//					//this->_interactDirect(&(inTargets[i]),it1,&(inTargets[i].RHS));
//
//					_interactDirect(&(currentLeaf->vparticles[i]),current,&(currentLeaf->vparticles[i].RHS));
//				}
//					cnt_direct+=currentLeaf->nParticles;
//					it3.advanceRemove();
//				}
//				else {
//					it3++;
//				}
//			}
//			
//			}
//		}
//		
//		std::cout << "finished with evaluate-serial, Interactions with Boxes: " << cnt_box << " direct: " <<cnt_direct << std::endl; 
//		
	}
	
	
#endif
#else //TUSEBOX2:
	
	
	template <typename tBox>
	bool recursivelyPutLists(tBox *node, std::list<tBox*> toHandle, const typename tBox::Btype theta)
	{
		//std::list<tBox* > giveDown;
		static const int nc=1<<tBox::Particle::dim;
		typename std::list<tBox* >::iterator it1=toHandle.begin();
		typename std::list<tBox* >::iterator tmpit;
		char isclose;
		bool res=false;
		
		std::cout << "Level: " << node->level <<  " toHandle.size= " << toHandle.size() << std::endl; 
		
		if(toHandle.empty())
		{
			return true; //nothing to do, everything is done on the parents basis already
		}
		
		if(node->isleaf)
		{
			while(it1!=toHandle.end())
			{
				isclose=isclose_box_barnes_hut(node, *it1, theta);
				if(!isclose) //it is far away -> so these guys (in it1) can interact with me as a box.
				{
					node->toDo_far.push_back(*it1);
//					tmpit=it1;
//					it1++;
//					toHandle.erase(tmpit);
					//std::cout << "Check one 601" <<std::endl;
				}
				else  //it is close two cases: 
				{
					if((*it1)->isleaf)
					{
						node->toDo_direct.push_back(*it1);
						//std::cout << "Check two 608" <<std::endl;

					}
					else {  //still has children need to go down there.
						BoxIterator<tBox> it_targets(*it1);
						it_targets++;
						//std::cout << "Check three 614" <<std::endl;
        				while(it_targets!=NULL)
						{
							tBox * current=it_targets;
							isclose=isclose_box_barnes_hut(node, current, theta);
							if(!isclose)
							{
								node->toDo_far.push_back(current);
								it_targets.advanceRemove();
							}
							else {
								if(it_targets->isleaf)
								{
									node->toDo_direct.push_back(current);
									it_targets.advanceRemove();
								}
								else
								{
									it_targets++;
								}
							}
						}

					}
				}
				it1++; //need to go through.
			}
			
			return true;
		}
		
		else {  //ok, this is difficult: we can either go down in the targets tree or in the source tree, criteria should probably be the level together with isclose==1
			
			//since we are now a level down compared to the call from above, the first thing is to check the "tohandle" again:
			while(it1!=toHandle.end())
			{
				isclose=isclose_box_barnes_hut(node, *it1,theta);
				if(!isclose)
				{
					node->toDo_far.push_back(*it1);
					tmpit=it1;
					it1++;
					toHandle.erase(tmpit);
				}
				else {
					if((*it1)->isleaf)
					{
						node->toDo_direct.push_back(*it1);
						tmpit=it1;
						it1++;
						toHandle.erase(tmpit);
					}
					//most difficult: it is close but still has children:
					if(node->level>=(*it1)->level && isclose==2) //no way that i am going to be far for them, so treat this at my childrens level.
					//if(node->level>(*it1)->level)
					{
						it1++;
						//std::cout << "Check four 670" <<std::endl;

					}
					else {  //mmh, might be that i appear far to some of those childrens boxes. if we do that, we take it out from toHandle, put possibly put the children into the list.
					
					  BoxIterator<tBox> it_targets(*it1);
					  it_targets++; //wanna only get the children.
					  while(it_targets!=NULL)
					  {
						  tBox* current=it_targets;
						  isclose=isclose_box_barnes_hut(node, current, theta);
						  if(!isclose)
						  {
							  node->toDo_far.push_back(current);
							  it_targets.advanceRemove();
							  //std::cout << "Check five 685" <<std::endl;

						  }
						  else 
						  {
							  if(current->isleaf) //i appear close to a leaf -> all those particles have to interact directly with mines.
							  {
								  node->toDo_direct.push_back(current);
								  it_targets.advanceRemove();
								  //std::cout << "Check six 694" <<std::endl;

							  }
							  else{
							  if(node->level>=current->level && isclose==2)//no way that i am going to be far for them, so treat this at my childrens level.
							  //if(isclose==2)
							  //if(node->level>current->level)
							  {
								  toHandle.push_back(current);
								  it_targets.advanceRemove();
								  //std::cout << "Check seven 702" <<std::endl;

							  }
							  else {
								  it_targets++; //further subdivide in he targets tree.
								  //std::cout << "Check eight 707	" <<std::endl;

							  }
							  }
						  }
						  
					  }
					  //done with the children, can continue:
					  //since we have split the parent into childrens, we need to delete it from the toHandle list.
						tmpit=it1;
						it1++;
						toHandle.erase(tmpit);
						
					  

					}
					
					
				}
			}
			
			//checked everything here in to handle->so lets give it down to the children.
			res=true;
			bool tmpbol;
			for (int kb=0;kb<nc;++kb)
			{
				tmpbol=recursivelyPutLists(node->children[kb],toHandle,theta);
				res=(tmpbol&&res);
			}
			return res;
			
		}
		
		
		
	}
	
	
	template <typename tExpansions, typename TargetPoint,  int _maxlevel>
	void Evaluator_serial<tExpansions,TargetPoint,_maxlevel>::evaluate(TargetPoint* inTargets, int nTargets, const typename TargetPoint::BaseType theta)
	{
		
		if(this->rootBox->nParticles<1)
			return;
		//		global_box=0;
		//		global_direct=0;
		//		global_self=0;
		
		size_t cnt_box(0);
		size_t cnt_direct(0);
		tbb::tick_count tstart,tend;
		double tdur(0);
		
		
		
#ifndef _FMMSILENT
		std::cout << "Evaluating SERIAL,USING BOX STRUCTURE Building lists at the sources , theta is: " << theta  <<std::endl;
#endif
		
		tstart=tbb::tick_count::now();
		std::list<tBox* > toHandle;
		toHandle.push_back(this->rootBox);
		recursivelyPutLists(this->rootBox,toHandle,theta);
		tend=tbb::tick_count::now();
		tdur=(tend-tstart).seconds();
#ifdef VERBOSE
		std::cout << "recursive interaction listPutting (put them at the treenodes) took: " << tdur <<std::endl;
#endif
		
		
		typename std::list<tBox* >::iterator listIterator;
		
		BoxIterator<tBox> it1(this->rootBox);
		
		while(it1!=NULL)
		{
			
			listIterator=it1->toDo_far.begin();
			while(listIterator!=it1->toDo_far.end())
			{
				Particle* ptarget=(*listIterator)->vparticles;
				for(int i=0;i<(*listIterator)->nParticles;++i)
				{
					it1->expansions.evaluateExpansions(&(ptarget->x[0]),&(ptarget->RHS));
					ptarget++;
				}
				cnt_box+=(*listIterator)->nParticles;
				listIterator++;
			}
			
			listIterator=it1->toDo_direct.begin();
			while(listIterator!=it1->toDo_direct.end())
			{
				
				for(int i=0;i<(*listIterator)->nParticles;++i)
				{
					_interactDirect(&((*listIterator)->vparticles[i]),it1,&((*listIterator)->vparticles[i].RHS));
				}
				cnt_direct+=(*listIterator)->nParticles*it1->nParticles;

				listIterator++;
			}
			
			it1++;
		}
#ifdef VERBOSE
		std::cout << "finished with evaluate-serial, Interactions with Boxes: " << cnt_box << " direct: " <<cnt_direct << std::endl;
#endif
	}
	
#endif

template <typename tExpansions, typename TargetPoint,  int _maxlevel>
	void Evaluator_serial<tExpansions,TargetPoint,_maxlevel>::_interact(TargetPoint* inTarget, typename TargetPoint::RHSType *rhs_out, const typename Particle::BaseType theta) const //,size_t &cnt_box, size_t &cnt_direct
{
	if(this->rootBox->nParticles<1)
		return;
	BoxIterator<tBox> it1(this->rootBox);
	
	
	while(it1!=NULL)
	{
		tBox* current=it1;
		if(ws_barnes_hut(current, inTarget->x,theta)) // it is far away  we can use the multipoles
		{
				it1->expansions.evaluateExpansions(inTarget->x,rhs_out);
				it1.advanceRemove(); //remove this branch from todo-List.
				//cnt_box++;
		}
		else //it is close:
		{
		    if(it1->isleaf) //it is a leaf so we interact directly
			{
				_interactDirect(inTarget,it1,rhs_out);
				//cnt_direct+=it1->nParticles;
				it1.advanceRemove(); //remove this branche from toDo list
			}
			else //it is not a leaf, so we further traverse here
			{
				it1++;
			}
		}
	}
	
	
	
	
	
	
}
	
	
	
////obsolete
//template <typename tExpansions, typename TargetPoint,  int _maxlevel>
//void Evaluator_serial<tExpansions,TargetPoint,_maxlevel>::_interact(TargetPoint* inTarget, typename TargetPoint::RHSType *rhs_out)
//{
//	typedef typename box_walker<Box<tExpansions,_maxlevel> >::t_p2Box_list t_p2Box_list;
//	typedef typename box_walker<Box<tExpansions,_maxlevel> >::t_p2Box_it t_p2Box_it;
//	
//	t_p2Box_it close_it,far_it;
//	
//	box_walker<Box<tExpansions,_maxlevel> > my_walker(MotherClass::rootBox);
//	my_walker._collectSplit(&(inTarget->x[0]));
//	
////	std::cout << "looping through the interaction lists:" << std::endl;
////	std::cout << "Close (direct) : " << my_walker.direct_list.size() << std::endl;
////	std::cout << "Far (use expansions) : " << my_walker.ws_list.size() << std::endl;
//
//	
//	close_it=my_walker.direct_list.begin();
//	while(close_it!=my_walker.direct_list.end()) //loop close particles.
//	{
//		//_interactDirect(inTarget,*close_it,&(inTarget->RHS));
//		_interactDirect(inTarget,*close_it,rhs_out);
//		close_it++;
//	}
//	far_it=my_walker.ws_list.begin();
//	while(far_it!=my_walker.ws_list.end()) //loop far boxes.
//	{
//		//(*far_it)->expansions.evaluateExpansions(inTarget->x,&(inTarget->RHS));
//        (*far_it)->expansions.evaluateExpansions(inTarget->x,rhs_out);
//		far_it++;
//	}
//	
//	
//	
//
//}

template <typename tExpansions, typename TargetPoint,  int _maxlevel>
void Evaluator_serial<tExpansions,TargetPoint,_maxlevel>::_interactDirect(TargetPoint* inTarget, Box<tExpansions,_maxlevel>* p_Box, typename TargetPoint::RHSType *rhs_out)const
{
	if(this->rootBox->nParticles<1)
		return;
	//std::cout << "interacting directly" <<std::endl;
	//assert(p_Box->isleaf==true);
	
	for (int i=0;i<p_Box->nParticles;++i)
	{

		*(rhs_out)+=inTarget->computeRHS(&(p_Box->vparticles[i]));
		
	}
}
	
	
//recursive function:
template <typename tExpansions, typename TargetPoint,  int _maxlevel>
	typename TargetPoint::RHSType Evaluator_serial<tExpansions,TargetPoint,_maxlevel>::traverseRecursive(TargetPoint* inTarget, tBox* node, const typename tExpansions::Btype theta) const
	{
		assert(node!=NULL);
		static const int nc=(1<<Particle::dim);
		tRHS res;
		//Case One: Well Separated: Use the Multipoles to calculate res and bail out
		if (ws_barnes_hut(node, inTarget->x, theta))
		{
			node->expansions.evaluateExpansions(inTarget->x,&res);
			return res;
		}
		else  // if not Well Separated
		{
			if (!node->isleaf)  //Case 2: its not a leaf so we further traverse into the tree (summing up the children and bailing out)
			{
				for( int kb=0;kb<nc;++kb)
				{
				   if (node->children[kb]!=NULL)
				   {
					 tRHS tmp;
					 tmp=traverseRecursive(inTarget,node->children[kb],theta);
					 res+=tmp;
				   }
				}
				return res;
			}
			else // Case 3: we have reached the bottom of the tree carry out direct interactions.
			{
				_interactDirect(inTarget,node,&res);
				return res;
			}
		}
		
	}


	

}

#endif

