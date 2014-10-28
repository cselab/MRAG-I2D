/*
 *  hcfmm_boxBuilder.h
 *  hcfmm
 *
 *  Created by Manfred on 1/21/09.
 *  Copyright 2009 ETHZ. All rights reserved.
 *
 */

#ifndef _HCFMM_BOX_BUILDER
#define _HCFMM_BOX_BUILDER
#include "hcfmm_box.h"

namespace HCFMM{

template <class tExpansions,int _maxlevel>
struct boxBuilder {
	typedef typename tExpansions::ParticleType Particle;
	void buildBoxes(Particle * in_p, int nParticles,	Box<tExpansions,_maxlevel>* rootBox);
    void generateExpansions(Box<tExpansions,_maxlevel>* rootBox);
	static void printStats(Box<tExpansions,_maxlevel>* rootBox);
	
	virtual ~boxBuilder();
	
};
	
template <class tExpansions,int _maxlevel>
	void boxBuilder<tExpansions,_maxlevel>::printStats(Box<tExpansions,_maxlevel>* rootBox)
	{
		typedef  Box<tExpansions,_maxlevel> tBox;
		BoxIterator<tBox> it1(rootBox);
		
		std::cout << "Tree Stats: " << std::endl;
		while(it1!=NULL)
		{
			std::cout << "Level: " << it1->level << " nParticles: " << it1->nParticles << " Leaf: " << it1->isleaf << std::endl;
			std::cout << "geometrical Center: [ ";
			for (int d=0;d<Particle::dim;++d)
			{ 
				std::cout << it1->center[d] << " ";
			}
			std::cout  << "]" << std::endl;
			std::cout << "COM: [ ";
			for (int d=0;d<Particle::dim;++d)
			{ 
				std::cout << it1->COM[d] << " ";
			}
			std::cout  << "]" << std::endl;
			std::cout << "h: [ ";
			for (int d=0;d<Particle::dim;++d)
			{ 
				std::cout << it1->h[d] << " ";
			}
			std::cout  << "]" << std::endl;
			std::cout << "Radius: " << it1->expansions.Radius << " TotalMass: " << it1->TotalMass << std::endl;
			it1++;
		}
	}
	
	
}

#endif