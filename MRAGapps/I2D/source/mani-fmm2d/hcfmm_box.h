/*
 *  hcfmm_box.h
 *  hcfmm
 *
 *  Created by Manfred on 12/22/08.
 *  Copyright 2008 ETHZ. All rights reserved.
 *
 */

#ifndef _HCFMM_BOX_
#define _HCFMM_BOX_

#include <cassert>
#include <vector>
#include "hcfmm_types.h"
#include <list>
#include <vector>
#include <algorithm>
#include <complex>
#include <queue>

namespace HCFMM{
	
	template <typename tExpansions, int _maxlevel>
	class Box 
	{
	public:
		//typedefs:
		typedef typename tExpansions::ParticleType Particle;
		typedef typename Particle::BaseType Btype;
		//private variables:
		bool isleaf;
		int level;
		Particle* vparticles;
		int nParticles; //if leaf: number of particles in vparticles, else, number of particles contained in children.
		//Btype expansions[Particle::pdim][_order][_order];
		//std::complex<Btype> expansions[Particle::pdim][_order][_order*2+1];
		tExpansions expansions;
		bool got_expansions;
		bool got_COM;
		Btype h[Particle::dim];
		Btype center[Particle::dim]; //geometrical center
		Btype COM[Particle::dim];    //center of mass
		Btype TotalMass;
		int maxlevelinuse;
		
		unsigned int id; // Added by Roman
		unsigned int getId () const {
			return id;
		}
		Real getExpansionsCenter (int d) const {
			return expansions.Center [d];
		}

		
		//This is a special try: storing the interaction lists at the nodes:
		std::list< Box<tExpansions,_maxlevel>* > toDo_far;
		std::list< Box<tExpansions,_maxlevel>* > toDo_direct;
		
		
		//bbox<Btype,Particle::dim> bounds;
		
		//Pointer Variables (careful here when deriving from this class)
		Box* children[1<<Particle::dim];
		Box* parent;
		
		//real public stuff:
		static const int maxlevel=_maxlevel;
		static const int order=tExpansions::order;
		
		
		
		//class management:
		//constructor:
		Box():
		isleaf(false),level(0), nParticles(0),TotalMass(0), vparticles(NULL),parent(NULL), got_expansions(false), got_COM(false),expansions(),maxlevelinuse(0),
		toDo_direct(),toDo_far(), id (0)
		{
			for (int j=0;j<(1<<Particle::dim);++j)
			{
				children[j]=NULL;
			}
			
			for (int d=0;d<Particle::dim;++d)
			{
				h[d]=Btype(0);
				COM[d]=Btype(0);
				center[d]=Btype(0);
			}
			
			
			
			
		}
		
		//Copy-Constructor:
		Box(const Box<tExpansions,maxlevel> &B):
		level(B.level),isleaf(B.isleaf),vparticles(B.vparticles),parent(B.parent),nParticles(B.nParticles),expansions(B.expansions),
		got_COM(B.got_COM),got_expansions(B.got_expansions)
		{
			for (int d=0;d<Particle::dim;++d)
			{
				h[d]=B.h[d];
				center[d]=B.center[d];
				COM[d]=B.COM[d];
			}
			//DIEGO COMMENTED THIS
			//for (int pd=0;pd<Particle::pdim;++pd)
			//{
			//	TotalMass[pd]=B.TotalMass[pd];
			//}
			//AND WROTE THIS
			TotalMass = B.TotalMass;
			
			for (int i=0; i<(1<<Particle::dim);++i)
			{
				children[i]=B.children[i];
			}
			
		}
		
		//assign operator:
		Box & operator=(const Box &B)
		{
			level(B.level);
			isleaf(B.isleaf);
			vparticles(B.vparticles);
			parent(B.parent);
			nParticles(B.nParticles);
			expansions(B.expansions),
			got_COM(B.got_COM);
			got_expansions(B.got_expansions);
			
			for (int d=0;d<Particle::dim;++d)
			{
				h[d]=B.h[d];
				center[d]=B.center[d];
				COM[d]=B.COM[d];
			}
			//DIEGO COMMENTED THIS
			//for (int pd=0;pd<Particle::pdim;++pd)
			//{
			//	TotalMass[pd]=B.TotalMass[pd];
			//}
			//AND WROTE THIS
			
			TotalMass = B.TotalMass;
			
			for (int i=0; i<(1<<Particle::dim);++i)
			{
				children[i]=B.children[i];
			}
		}
		
		
		//destructor:
		~Box()
		{
			//std::cerr << "Calling the Box-Destructor, level: " << level <<std::endl;

			if (!isleaf)
			{
				for (int i=0; i<(1<<Particle::dim);++i)
				{
					if (children[i]!=NULL)
					{
						delete children[i];
						children[i] = NULL;
					}
					
				}
			}
			else
			{
				vparticles=NULL;
			}
			
		}
		
		
		
	};
	
	
	////////////////box-iterator/////////////////
	template <typename tBox, template <typename T> class tAllocator = std::allocator>
	class BoxIterator
	{
	public:
		typedef typename tBox::Particle Particle;
		typedef unsigned short int ts;
		typedef tBox* p2Box;
		typedef std::deque< p2Box,tAllocator<p2Box> > tContainerQueue;
		typedef std::queue< p2Box, tContainerQueue> tQueue;
		ts nc;
		
		BoxIterator(tBox* in_root) : nc((1<<Particle::dim))
		{
			if (in_root!=NULL)
				toDo.push(in_root);
		}
		
		
		inline tBox& operator*() const {
			assert (!toDo.empty());
			return *(toDo.front());
		}
		
		inline operator tBox*() const {
			return toDo.front();
		}
		
		inline void operator++(int)
		{
			assert(!toDo.empty());
			tBox* current;
			current=toDo.front();
			toDo.pop();
			_expand(current);
		}
		
		//this one removes the top entry without expanding children
		inline void advanceRemove()
		{
			if(!toDo.empty())
				toDo.pop();
		}
		
		inline tBox* operator->() const {
			if (!toDo.empty())
			{
				return toDo.front();
			}
			else{
				return NULL;
			}
		}
		
		inline bool operator==(const tBox *other) const {
			return ( (toDo.empty() && other==NULL ) || (toDo.front()==other) );
		}
		
		inline bool operator!=(const tBox *other) const {
			return !(*this==other);
		}
		
		inline void reset()
		{
			//toDo.clear();
			//TODO: ther must be a nicer way than this:
			while(!toDo.empty())
			{
				toDo.pop();
			}
			toDo.push(root);
		}
		
	private:
		//data-members:
		//int pos[_maxlevel][1<<Particle::dim];
		tQueue toDo;
		p2Box root;
		
		
		//functions:
		inline void _expand(tBox* in_node)
		{
			for (ts i=0;i<nc;++i)
			{
				if(in_node->children[i]!=NULL)
				{
					toDo.push(in_node->children[i]);
				}
			}
		}
		
	};
	
	
	
	
	
	////////////////box-walker//////////////////
	//inefficient, uses lists
	template <class tBox>  //should only depend on a box.
	struct box_ref{
		
		//data members:
		tBox* p_Box;
		bool expanded;
		
		//constructor:
		box_ref():p_Box(NULL),expanded(true){}
		//copy-constructor:
		box_ref(const box_ref &B):p_Box(B.p_Box),expanded(B.expanded){}
		//destructor:
		~box_ref(){p_Box=NULL;}
		//operator=:
		box_ref & operator=(const box_ref &B){p_Box=B.p_Box;expanded=B.expanded; return *this;}
		
		
		inline bool operator<(const box_ref &b)
		{
			return this->level<b.level;
		}
		
		
		
	};
	
	template <class tBox>
	inline bool operator<(const box_ref<tBox> &a, const box_ref<tBox> &b)
	{
		//(a.level<b.level) return true;
		//else return(a.p_Box<b.p_Box);
		return(a.p_Box->level<b.p_Box->level);
	}
	
	template <class tBox>
	inline bool operator>(const box_ref<tBox> &a, const box_ref<tBox> &b)
	{
		//if(a.level>b.level) return true;
		//else return(a.p_Box>b.p_Box);
		return(a.p_Box->level>b.p_Box->level);
		
	}
	
	
	template <class tBox>
	struct box_walker
	{
		typedef typename tBox::Particle Particle;
		
		
		
		typedef std::list<box_ref<tBox> > t_ref_list;
		typedef typename std::list<box_ref<tBox> >::iterator t_ref_it;	
		//	typedef std::list<tBox* > t_p2Box_list;
		//	typedef typename std::list<tBox* >::iterator t_p2Box_it;
		
		//	typedef std::vector<box_ref<tBox> > t_ref_list;
		//	typedef typename std::vector<box_ref<tBox> >::iterator t_ref_it;
		typedef std::vector<tBox* > t_p2Box_list;
		typedef typename std::vector<tBox* >::iterator t_p2Box_it;
		
		
		
		
		box_walker(tBox* rootNode):
		_rootBox(rootNode),direct_list(),ws_list(),work_list()
		{
			size_t maxLength=1<<(Particle::dim*rootNode->maxlevelinuse);
			//		work_list.reserve(maxLength);
			ws_list.reserve(maxLength);
			direct_list.reserve(maxLength);
			
		}
		
		box_walker(const box_walker &B):
		_rootBox(B.rootBox),direct_list(B.direct_list),ws_list(B.direct_list),work_list(B.work_list)
		{
		}
		
		~box_walker()
		{
			_rootBox=NULL; //we don't want to destroy the rootBox.
		}
		
		box_walker& operator=(const box_walker& B)
		{
			work_list=B.work_list;
			ws_list=B.ws_list;
			direct_list=B.direct_list;
			_rootBox=B._rootBox;
		}
		
		void _collect();
		void _collectSplit(typename tBox::Btype *tpx);
		void callbottomup(void (*p_func)(tBox*));
		void calltopbottom(void (*p_func)(tBox*));
		void printStats();
		void checkExpansions();
		t_ref_list work_list;
		//how big are empty lists?
		t_p2Box_list ws_list;
		t_p2Box_list direct_list;
		
		
		//private:
		tBox* _rootBox;
	};
	
	template <class tBox>
	void box_walker<tBox>::_collectSplit(typename tBox::Btype *tpx) //tpoint x.
	{
		const int minParticles=1;
		work_list.clear();
		ws_list.clear();
		direct_list.clear();
		//std::cout << "[box_walker: _collectSplit: ] collecting elements" <<std::endl;
		tBox* current=_rootBox;
		box_ref<tBox> tmp;
		//std::list<box_ref<tBox> > res;
		t_ref_it it1;
		
		tmp.p_Box=current;
		tmp.expanded=false;
		//res.push_back(tmp);
		//it1=res.begin();
		//while(it1!=res.end())
		work_list.push_back(tmp);
		it1=work_list.begin();
		while(it1!=work_list.end())
		{
			//if well separated add to wslist, else pushback and check children.
			if(ws_barnes_hut(it1->p_Box,tpx) && it1->p_Box->nParticles>=minParticles)
			{
				ws_list.push_back(it1->p_Box);
			}
			else
			{
				if(it1->expanded==false)//we have not looked at the children, and it this box is not well separated
				{
					//std::cout << "I am at level " << it1->level << "adding children to toDo-list" <<std::endl;
					current=it1->p_Box;
					for (int kb=0;kb<(1<<Particle::dim);++kb)
					{
						if (current->children[kb]!=NULL)
						{
							if(ws_barnes_hut(current->children[kb],tpx) && it1->p_Box->nParticles>=minParticles) //the kid is well separated
							{
								ws_list.push_back(current->children[kb]); //-->we add it to the list and continue with next kid
							}
							else 
							{
								if(current->children[kb]->isleaf) //it is not well-separated and a leaf
								{
									direct_list.push_back(current->children[kb]); //we want to interact with the particles of that leaf directly
								}
								else  //it is not well-separated and not a leaf, so we need to loop over the children of this kid. -->put it in the toDo list.
								{
									tmp.p_Box=current->children[kb];
									tmp.expanded=false;
									//res.push_back(tmp);
									work_list.push_back(tmp);
								}
							}
						}
						
					}
					it1->expanded=true;
				}
			}
			
			it1++;
		}
		
		//std::cout << "collected " << res.size() << " items." << std::endl;
		//std::cout << "Close: " << direct_list.size() << " Far: " << ws_list.size() << std::endl;
		//work_list=res;
		
		//TODO: should I do this?
		//	size_t s=direct_list.size();
		//	direct_list.resize(s);
		//	s=ws_list.size();
		//	ws_list.resize(s);
		//	s=work_list.size();
		//	work_list.resize(s);
		
		
	}
	
	
	template <class tBox>
	void box_walker<tBox>::_collect()
	{
		
		work_list.clear();
#ifndef _FMMSILENT
		std::cout << "collecting elements" <<std::endl;
#endif
		tBox* current=_rootBox;
		box_ref<tBox> tmp;
		//	std::list<box_ref<tBox> > res;
		t_ref_it it1;
		tmp.p_Box=current;
		tmp.expanded=false;
		//	res.push_back(tmp);
		work_list.push_back(tmp);
		//	it1=res.begin();
		it1=work_list.begin();
		//	while(it1!=res.end())
		while(it1!=work_list.end())
		{
			
			if(it1->expanded==false)
			{
#ifdef VERBOSE
				std::cout << "I am at level " << it1->p_Box->level << "adding children to toDo-list" <<std::endl;
#endif
				current=it1->p_Box;
				for (int kb=0;kb<(1<<Particle::dim);++kb)
				{
					if (current->children[kb]!=NULL)
					{
						tmp.p_Box=current->children[kb];
						tmp.expanded=false;
						//res.push_back(tmp);
						work_list.push_back(tmp);
					}
					
				}
			}
			it1->expanded=true;
			it1++;
		}
#ifdef VERBOSE
		std::cout << "collected " << work_list.size() << " items." << std::endl;
#endif
		//work_list=res;
		
		
		
	}
	
	template <class tBox>
	void box_walker<tBox>::callbottomup(void (*p_func)(tBox*))
	{
		typename std::list<box_ref<tBox> >::iterator it1;
		_collect();
		//note: the list should be automatically sorted in the way we collect the work.
		it1=work_list.end();
		while(it1!=work_list.begin())
		{
			it1--;
			p_func(it1->p_Box);
		}
	}
	
	template <class tBox>
	void box_walker<tBox>::calltopbottom(void (*p_func)(tBox*))
	{
		typename std::list<box_ref<tBox> >::iterator it1;
		_collect();
		//note: the list should be automatically sorted in the way we collect the work.
		it1=work_list.begin();
		while(it1!=work_list.end())
		{
			p_func(it1->p_Box);
			it1++;
		}
	}
	
	
	
	
} //namespace

#endif
