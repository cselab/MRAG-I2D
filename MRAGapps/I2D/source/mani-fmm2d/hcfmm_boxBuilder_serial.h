/*
 *  hcfmm_box_serial.h
 *  hcfmm
 *
 *  Created by Manfred on 1/14/09.
 *  Copyright 2009 ETHZ. All rights reserved.
 *
 */


#ifndef HCFMM_BOXBUILDER_SER
#define HCFMM_BOXBUILDER_SER
#include <cassert>
#include "hcfmm_boxBuilder.h"
#include "hcfmm_operators_serial.h"
#include <algorithm>

namespace HCFMM{
	
	template <typename tExpansions, int _maxlevel>
	struct boxBuilder_serial: public boxBuilder<tExpansions,_maxlevel>
	{
		//typdefs:
		typedef boxBuilder<tExpansions,_maxlevel>MotherClass;
		typedef boxBuilder_serial<tExpansions,_maxlevel> SelfClass;
		typedef typename tExpansions::ParticleType Particle;
		typedef typename Particle::BaseType Btype;
		typedef Box<tExpansions,_maxlevel> tBox;
		//public functions
		static void buildBoxes(Particle * in_p, int nParticles, tBox* rootBox);
		static void generateExpansions(tBox* rootBox);
		boxBuilder_serial():MotherClass(){};
		~boxBuilder_serial(){};
		//private functions:
		static void _split(tBox* rootBox, int levels=_maxlevel);
		static void _calculateCOM(tBox* rootBox);
		static void _computeRadius(tBox* rootBox);
		static int _getMaxLevelinuse(tBox* rootBox);

	};
	
	template <typename tExpansions, int _maxlevel>
    int boxBuilder_serial<tExpansions,_maxlevel>::_getMaxLevelinuse(Box<tExpansions,_maxlevel>* rootBox)
	{
		tBox* current=rootBox;
		int cmax=0;
		if(current!=NULL)
		{
			cmax=current->level;
			for (int kb=0;kb<(1<<Particle::dim);++kb)
			{
				if(current->children[kb]!=NULL)
				{
					cmax=max(cmax,_getMaxLevelinuse(current->children[kb]));
				}
			}
		}
		
		return cmax;
		
	}
	
	
	template <typename tExpansions, int _maxlevel>
    void boxBuilder_serial<tExpansions,_maxlevel>::buildBoxes(Particle * in_p, int in_nParticles,	Box<tExpansions,_maxlevel>* rootBox)
	{

#ifndef _FMMSILENT
		std::cout<<"Called box_builder Serial " <<std::endl;
#endif
		assert(rootBox->level==0 && rootBox->vparticles==NULL);
		if (in_nParticles<1)
		{
#ifndef _FMMSILENT
			std::cout << "[boxBuilder]: no Source particles found. creating empty box" <<std::endl;
#endif
			return;
		}
		bbox<Btype,Particle::dim> bounds=getBoundingBox<Btype,Particle,Particle::dim>(in_p,in_nParticles);
		//rootBox->bounds=getBoundingBox<Btype,Particle,Particle::dim>(in_p,in_nParticles);
		for (int d=0;d<Particle::dim;++d)
		{
			rootBox->h[d]=bounds.upper[d]-bounds.lower[d];
			rootBox->center[d]=bounds.lower[d]+rootBox->h[d]/Btype(2.0);
		}
		rootBox->nParticles=in_nParticles;
		rootBox->vparticles=in_p;
#ifndef _FMMSILENT
		std::cout << "Having " << rootBox->nParticles << " items in initial Box"  <<std::endl;
		std::cout << "bBox is: " <<std::endl;
		bounds.print();
#endif
		SelfClass::_split(rootBox);
		rootBox->maxlevelinuse=_getMaxLevelinuse(rootBox);
#ifndef _FMMSILENT
		std::cout << "Max Level in use is: " << rootBox->maxlevelinuse << "of " << _maxlevel << "possible" <<std::endl;
#endif
		
	}
	
	

	template <typename tExpansions, int _maxlevel>
    void boxBuilder_serial<tExpansions,_maxlevel>::_split(Box<tExpansions,_maxlevel>* rootBox, int levels)
	{
		//implement splitting here:
		//1)create a  list of keys with pointers to the particles
		//2)sort the keys
		//3)rearrange the particles
		//4)create children with pointers to the sorted particles
		//5) split the kid!
		
		
		//hardcoded: //minimum particles per box.
#ifndef _SMAX
		const int smax=20;
#else
		const int smax=_SMAX;
#endif
		p_ind<Particle> *pmap;
		
		assert(rootBox->children[0]==NULL);
		
		if (rootBox->level+1<levels)
		{
		
		//std::vector<p_ind<Particle> > pmap;
		//new:
		//pmap.reserve(rootBox->nParticles);	
		
			//TODO: Improve this, its is not necc. to copy anyway.
		//dont need a vector, I know the size:
		pmap=new p_ind<Particle>[rootBox->nParticles];
			

		int children_count[(1<<Particle::dim)];
		int kb=0;
		for (int k=0;k<(1<<Particle::dim);++k)
		{
			children_count[k]=0;
		}
		
		
		//1
		//we create a full blown copy:
		Particle* tmpArr=new Particle[rootBox->nParticles];
		memcpy(tmpArr, rootBox->vparticles,rootBox->nParticles*sizeof(Particle));
		for (int i=0; i<rootBox->nParticles;++i)
		{
			p_ind<Particle> dummy=createpIndex(&(tmpArr[i]), rootBox->center,rootBox->h);
			//pmap.push_back(dummy);
			pmap[i]=dummy;
			kb=pmap[i].key;
			assert(kb<(1<<Particle::dim));
			children_count[kb]++;
		}
		
		//2
		//std::sort(pmap.begin(),pmap.end());
		//std::sort(pmap,&pmap[rootBox->nParticles]);
		if (rootBox->nParticles >= 20000)
		   tbb::parallel_sort(pmap,&pmap[rootBox->nParticles]);
		else
		   std::sort(pmap,&pmap[rootBox->nParticles]);	
		
#ifdef VERBOSE
		std::cout << "sorted keys" << std::endl;
#endif
		//3
				
		for (int i=0;i<rootBox->nParticles;++i)
		{
			//memcpy(&(rootBox->vparticles[i]),pmap[i].p,sizeof(Particle));
			rootBox->vparticles[i]=*(pmap[i].p);
#ifdef VERBOSE
			std::cout  << "particle " << i << ": pos: ["<< rootBox->vparticles[i].x[0] << " " << rootBox->vparticles[i].x[1] << " " << rootBox->vparticles[i].x[2] << " ] key: " << pmap[i].key <<std::endl;  
#endif
		}
		//delete[] vtmp; //THIS DOES INVALIDATE THE CONTIGUITY of THE BASe-ARRAY!:
		delete[] tmpArr; //This should be ok, because it is an array that we temporarily generated.
		delete[] pmap;
		
		//checking only:
		/*
		std::cout << "checking sorted array:" <<std::endl;
		for (int i=1; i<rootBox->nParticles;++i)
		{
			p_ind<Particle> tmpk1=createpIndex(&(rootBox->vparticles[i-1]),rootBox->center,rootBox->h);
			p_ind<Particle> tmpk2=createpIndex(&(rootBox->vparticles[i]), rootBox->center,rootBox->h);
			//tmpk1.key=ot_hkey_TruncateAtDepth(tmpk1.key,1);
			//tmpk2.key=ot_hkey_TruncateAtDepth(tmpk2.key,1);
			std::cout << i << ": " << rootBox->vparticles[i].x[0] << "," << rootBox->vparticles[i].x[1] << "," << rootBox->vparticles[i].x[2] << " : " << pmap[i].key << "---";
			if ( tmpk1 <= tmpk2 )
				std::cout << "ok. "<<std::endl;
			else
				std::cout << "error!: prev: "<<tmpk1.key << " cur: " <<tmpk2.key <<std::endl;
			
		}
		*/
		
		
		
		//4:create children:
#ifdef VERBOSE
		std::cout << "starting to create children:" <<std::endl;
#endif
		kb=0;
		int cntkids=0;
		int kids_st(0),kids_nd(0);
		int kbits[Particle::dim];
		rootBox->isleaf=false;

		static unsigned int current_id = 1;

		while(kb<(1<<Particle::dim))
		{
			if(children_count[kb]==0)
			{
				rootBox->children[kb]=NULL;
			}
			else
			{
				
				//this->children[kb]->makeChild(children_count[kb], &(this->vparticles[kids_st]), this);

				cntkids+=1;
				kids_nd+=children_count[kb];

				rootBox->children[kb]=new Box<tExpansions,_maxlevel>;

				////////// Added by Roman
				rootBox->children[kb]->id = current_id;
				++current_id;
				//////////

				rootBox->children[kb]->parent=rootBox;
				rootBox->children[kb]->level=rootBox->level+1;
				rootBox->children[kb]->isleaf=true;
				rootBox->children[kb]->vparticles=&(rootBox->vparticles[kids_st]);
				rootBox->children[kb]->nParticles=children_count[kb];
				lsfkey2bits(kb, Particle::dim, kbits);
				for (int d=0; d<Particle::dim;++d)
				{
				  (kbits[d]==0)?kbits[d]=-1:kbits[d]=1;
				  rootBox->children[kb]->h[d]=rootBox->h[d]/Btype(2);
				  rootBox->children[kb]->center[d]=rootBox->center[d]+Btype(kbits[d])*rootBox->h[d]/Btype(4);
				}

				
				//std::cout << "Bounding box of new child:" <<std::endl;
				//rootBox->children[kb]->bounds.print();
#ifdef VERBOSE
				std::cout << "Center of new child: [ "; 
				for (int d=0;d<Particle::dim;++d)
				{
					std::cout << rootBox->children[kb]->center[d]  <<" ";
				}
				std::cout << " ] " <<std::endl;
				std::cout << "number of particles of child" << kb << ": " << rootBox->children[kb]->nParticles << std::endl;
#endif
				kids_st=kids_nd;
			}
			kb++;
		}
		assert(kids_nd==rootBox->nParticles);
#ifdef VERBOSE
		std::cout << "created " << cntkids << "children" <<std::endl; 
#endif	
		//testing purpose:
		int sm=0;
		int chkNKids=0;
		for (kb=0;kb<(1<<Particle::dim);++kb)
		{
			if (rootBox->children[kb]!=NULL)
			{
#ifdef VERBOSE
				std::cout<<rootBox->children[kb]->nParticles<< std::endl;
#endif
				sm+=rootBox->children[kb]->nParticles;
				chkNKids+=1;
			}
		}
#ifdef VERBOSE
		std::cout << "in total counted: " << sm << std::endl;
#endif
		assert(rootBox->nParticles==sm);
		assert(cntkids==chkNKids);
		
			
		
		//5 split children again:
		if (rootBox->level+1<levels)
		{
#ifdef VERBOSE
		std::cout << "splitting again..." <<std::endl;
#endif

			for (int kb=0; kb<(1<<Particle::dim);++kb)
			{
				if (rootBox->children[kb]!=NULL)
				{
					if (rootBox->children[kb]->nParticles>smax)
					{
						_split(rootBox->children[kb],levels);
					}
					else
					{
#ifdef VERBOSE
						std::cout << "contains less than " << smax << "Particles, not splitting." <<std::endl;
#endif
					}
				}
				else
				{
#ifdef VERBOSE

					std::cout << kb << " is a empty box(no particles in this sector):  not splitting." <<std::endl;
#endif
				}
			}
		 }
		else
		{
#ifdef VERBOSE
			std::cout << "reached Maxlevel, stopping to split" <<std::endl;
#endif 
		}
		}
		else
		{
#ifdef VERBOSE
			std::cout << "not splitting, reached maxlevel (from the beginning)" <<std::endl;
#endif
			rootBox->isleaf=true;
#ifdef VERBOSE
			std::cout << "setting rootBox->isleaf to: " << rootBox->isleaf <<std::endl;
#endif
		}
	}
 
	
	template <typename tExpansions, int _maxlevel>
    void boxBuilder_serial<tExpansions,_maxlevel>::_computeRadius(Box<tExpansions,_maxlevel>* rootBox)
	{
		typename Particle::BaseType tmpR=0;
		
		if(rootBox->isleaf) //compute based on particles
		{
			
			for (int i=0;i<rootBox->nParticles;++i)
			{
				tmpR=0;
				for (int d=0;d<Particle::dim;++d)
				{
					tmpR+=(rootBox->COM[d]-rootBox->vparticles[i].x[d])*(rootBox->COM[d]-rootBox->vparticles[i].x[d]);
				}
				tmpR=sqrt(tmpR);
				rootBox->expansions.Radius=max(rootBox->expansions.Radius,tmpR);
			}
		}
		else               //compute based on children.
		{
			for (int kb=0;kb<(1<<Particle::dim);++kb)
			{
				if(rootBox->children[kb]!=NULL)
				{
				tmpR=0;
				for (int d=0;d<Particle::dim;++d)
				{
					//tmpR+=(rootBox->COM[d]-rootBox->children[kb]->COM[d])*(rootBox->COM[d]-rootBox->children[kb]->COM[d]);
                    // wim: replaced with correct expression: need to take radius of children into account when computing parent radius
                    tmpR+=(rootBox->COM[d]-rootBox->children[kb]->COM[d]-rootBox->children[kb]->expansions.Radius)*(rootBox->COM[d]-rootBox->children[kb]->COM[d]-rootBox->children[kb]->expansions.Radius);
                }
				tmpR=sqrt(tmpR);
				rootBox->expansions.Radius=max(rootBox->expansions.Radius,tmpR);
                }
			}
		}
		
	}
	
  	
	template <typename tExpansions, int _maxlevel>
	void boxBuilder_serial<tExpansions,_maxlevel>::_calculateCOM(Box<tExpansions,_maxlevel>* rootBox)
	{
	 
		assert(rootBox->got_COM==false);
		typename Particle::BaseType currentMass = 0; //, totalMass = 0;
		//Nullify:
		for (int d=0;d<Particle::dim;++d)
		{
			rootBox->COM[d]=0;
		}
		
		rootBox->TotalMass=0;
		
	   //check leaf-status
       if(!(rootBox->isleaf)) //gather from children
	   {
		   
		   for (int kb=0;kb<(1<<Particle::dim);++kb)
		   {
			   if(rootBox->children[kb]!=NULL)
			   { 
				   currentMass=rootBox->children[kb]->TotalMass;
				   assert(currentMass>=Btype(0)); //TODO: what happens if mass is zero-> the center will be zero?
				   assert(rootBox->children[kb]->got_COM==true);
				   for(int d=0;d<Particle::dim;++d)
				   {
					   rootBox->COM[d]+=rootBox->children[kb]->TotalMass*rootBox->children[kb]->COM[d];
				   }
				   rootBox->TotalMass+=currentMass;
				   
			   }
			   
		   }
	   }
	   else  //calculate based on particles.
	   {
		   for (int i=0;i<rootBox->nParticles;++i)
		   {
			   currentMass=0;
			   for (int pd=0;pd<Particle::pdim;++pd)
			   {
				   currentMass+=(rootBox->vparticles[i].w[pd])*(rootBox->vparticles[i].w[pd]);
			   }
			   currentMass=sqrt(currentMass);
			   rootBox->TotalMass+=currentMass;
			   for(int d=0;d<Particle::dim;++d)
			   {
				   rootBox->COM[d]+=currentMass*rootBox->vparticles[i].x[d];
			   }
			   
		   }
		}
		
		//ToDo: Check is this correct or should we keep a separate COM for each pd?

		//correct with totalMass.
		for(int d=0;d<Particle::dim;++d)
		{
			rootBox->COM[d]/=rootBox->TotalMass;
		}
		
		rootBox->got_COM=true;
		
		//debug...:
//		if (rootBox->parent==NULL)
//		{
//			std::cout << "[boxBuilder_serial:] got COM: [ " << rootBox->COM[0]<< " " << rootBox->COM[1] << " " << rootBox->COM[2] << " ]" << std::endl;
//			std::cout << "[boxBuilder_serial:] TotalMass: " << rootBox->TotalMass << std::endl;
//		}
//		
//		//debug:...:
//		double COM_serial[3];
//		computeCOM_SingleCPU(rootBox->vparticles, rootBox->nParticles, COM_serial);
//		std::cout << "[boxBuilder_serial:] got COM(based on children): [ " << rootBox->COM[0]<< " " << rootBox->COM[1] << " " << rootBox->COM[2] << " ]" << std::endl;
//        std::cout << "[boxBuilder_serial:] SingleCPU: [ " << COM_serial[0]<< " " << COM_serial[1] << " " << COM_serial[2] << " ]" << std::endl;
//        
		
	}
	
	template <typename tExpansions, int _maxlevel>
	void boxBuilder_serial<tExpansions,_maxlevel>::generateExpansions(Box<tExpansions,_maxlevel>* rootBox)
	{
		
		//1) get an upward iterator
		//2) calculate expansions at lowest level
		//3) shift expansions upward
		//-------------------------//
#ifndef _FMMSILENT
		std::cout << "Generating Expansions" <<std::endl;
#endif
#ifdef _VERBOSE
		std::cout << "RootBox->isleaf? " << rootBox->isleaf << std::endl;
#endif
		
		typename box_walker<Box<tExpansions,_maxlevel> >::t_ref_it up_it;
		
		Btype myc[3]={0,0,0};
		//Btype myc[3]={0.00967776, -0.0131336, 0.0376924};

		
		//1) get an iterator:
		box_walker<Box<tExpansions,_maxlevel> > box_iterator(rootBox);
		box_iterator._collect();
		up_it=box_iterator.work_list.end();
		//ToDo: check whether the work_list assures that the expansions of the children are always calculated
		//(should be like that, since thats the way we build up the work_list. to be sure I added a boolean.)
		
		while (up_it != box_iterator.work_list.begin())
		{
			
			up_it--;
			//2) calculate expansions at lowest level (->isleaf==true)
			if(up_it->p_Box->isleaf==true)
			{
#ifdef VERBOSE
				std::cout << "[generateExpansions]handling leaf" << std::endl;
				std::cout << "leaf-status: " << up_it->p_Box->isleaf <<std::endl;
#endif
				_calculateCOM(up_it->p_Box);
//TODO SET CENTER

				//up_it->p_Box->expansions.setcenter(up_it->p_Box->center);
//				if (_TESTCASE_==1 || _TESTCASE_==4 )
//				{
//					up_it->p_Box->expansions.setcenter(myc);
//				}
//				else
//				{
                   up_it->p_Box->expansions.setcenter(&(up_it->p_Box->COM[0]));
//				}

				_computeRadius(up_it->p_Box);
				up_it->p_Box->expansions.calculateExpansions(up_it->p_Box->vparticles,up_it->p_Box->nParticles);
				up_it->p_Box->got_expansions=true;
				
			}
			//3) shift expansions of children to center of current Box and sum up.
			else
			{
//				if (_TESTCASE_==4)
//				{
//				std::cout << "[generateExpansions]handling branch" << std::endl;
//				std::cout << "I SHOULD NOT BE HERE in TEST 4 (set maxlevel=0)" <<std::endl;
//				std::cout << "leaf-status: " << up_it->p_Box->isleaf <<std::endl;
//				std::cout << "test-true: "  << true << std::endl;
//				}

				_calculateCOM(up_it->p_Box);
				//gatherExpansions(up_it->p_Box);
				assert(!up_it->p_Box->isleaf); //can't gather as leaf.
				assert(!up_it->p_Box->got_expansions); //don't do it twice.
				up_it->p_Box->expansions.clear();
//TODO SET CENTER
				//up_it->p_Box->expansions.setcenter(up_it->p_Box->center);

//				if (_TESTCASE_==1 || _TESTCASE_==4)
					up_it->p_Box->expansions.setcenter(myc);
//				else
					up_it->p_Box->expansions.setcenter(&(up_it->p_Box->COM[0]));
				_computeRadius(up_it->p_Box);
				for (int kb=0;kb<(1<<Particle::dim);++kb)
				{
					if(up_it->p_Box->children[kb]!=NULL)
					{
						assert(up_it->p_Box->children[kb]->got_expansions);
						assert(up_it->p_Box->children[kb]->got_COM);
						up_it->p_Box->expansions.gatherExpansions(&(up_it->p_Box->children[kb]->expansions));
					}
				}
				up_it->p_Box->got_expansions=true;
				
			}

			
		}
		
		up_it=box_iterator.work_list.end();
#ifndef _FMMSILENT
		if (rootBox->nParticles>0)
		{
		std::cout << "I have recursively calculated the expansions." <<std::endl;
		std::cout << "Expansions at root: " << std::endl;
        rootBox->expansions.print();
		std::cout << "COM of Root: [ ";
		for (int d=0;d<Particle::dim;++d)
		{
			std::cout << rootBox->COM[d] << " " ;
		}
		std::cout << "]" <<std::endl;
		
		std::cout << "TotalMass of Root: [ " <<rootBox->TotalMass << "]" <<std::endl;
		}
#endif
 		
	}

	

	
} //namespace

#endif
