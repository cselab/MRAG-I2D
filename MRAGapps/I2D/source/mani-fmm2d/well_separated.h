/*
 *  well_separated.h
 *  hcfmm
 *
 *  Created by Manfred on 1/27/09.
 *  Copyright 2009 ETHZ. All rights reserved.
 *
 */

#ifndef _HCFMM_WS_
#define _HCFMM_WS_

#include "hcfmm_box.h"
#include <cmath>
#include <limits>

namespace HCFMM{
	
	template<class tBox>
	inline	bool ws_barnes_hut(tBox* checkBox, typename tBox::Btype* tpx, typename tBox::Btype theta)
	{	
		typename tBox::Btype p2b_distance(0),rBox(0);
		
		for (int d=0;d<tBox::Particle::dim;++d)
		{
			p2b_distance += (checkBox->expansions.Center[d]-tpx[d])*(checkBox->expansions.Center[d]-tpx[d]); //center of mass
			rBox += checkBox->h[d]*checkBox->h[d];
		}
		
		return (sqrt(rBox/p2b_distance)<theta);
	}
	
	//we want to be able to distinguish between far(0) close(1) and very close(2)
	template<class tBox>
	inline int isclose_box_barnes_hut(tBox* sourceBox, tBox* targetBox, typename tBox::Btype theta)
	{
		//std::cout << "calling isclose... " <<std::endl;
		int res;
		typename tBox::Btype b2b_dist(0),tmp(0),rBox(0),rTarget(0);
		for (int d=0;d<tBox::Particle::dim;++d)
		{
			//tmp=min(std::numeric_limits<typename tBox::Particle::BaseType>::epsilon(),fabs(sourceBox->center[d]-targetBox->center[d])-(sourceBox->h[d]+targetBox->h[d]));
			tmp=sourceBox->center[d]-targetBox->center[d];
			b2b_dist+=tmp*tmp;
			rBox+=sourceBox->h[d]*sourceBox->h[d];
			rTarget+=targetBox->h[d]*targetBox->h[d];
		}
        tmp=sqrt(rTarget);
		//res=(sqrt(rBox)/sqrt(b2b_dist-tmp)<theta);
		
		tmp=max(sqrt(b2b_dist)-tmp,std::numeric_limits<typename tBox::Btype>::epsilon());
	    //std::cout << (sqrt(rBox)/tmp) << " vs " << theta <<std::endl;
		
	    if((sqrt(rBox)/tmp)<theta)
		{
			res=0;
			return res;
		}
		else 
		{  //close or very close?
			//std::cout << "box is close" <<std::endl;
			res=2; 
			if(sqrt(rBox)/sqrt(b2b_dist)<theta)
				res=1;
			return res;
			
		}
		
	}
	
	
	
	
	
	
	template<class tBox>
	inline bool ws_box_barnes_hut(tBox* sourceBox, tBox* targetBox, typename tBox::Btype theta)
	{
		bool res;
		typename tBox::Btype b2b_dist(0),tmp(0),rBox(0),rTarget(0);
		for (int d=0;d<tBox::Particle::dim;++d)
		{
			//tmp=min(std::numeric_limits<typename tBox::Particle::BaseType>::epsilon(),fabs(sourceBox->center[d]-targetBox->center[d])-(sourceBox->h[d]+targetBox->h[d]));
			tmp=sourceBox->center[d]-targetBox->center[d];
			b2b_dist+=tmp*tmp;
			rBox+=sourceBox->h[d]*sourceBox->h[d];
			rTarget+=targetBox->h[d]*targetBox->h[d];
		}
        tmp=sqrt(rTarget);
		//res=(sqrt(rBox)/sqrt(b2b_dist-tmp)<theta);
		
		tmp=max(sqrt(b2b_dist)-tmp,std::numeric_limits<typename tBox::Btype>::epsilon());
		res=(sqrt(rBox)/tmp<theta);
		
		return res;
	}
	
	
}

#endif
