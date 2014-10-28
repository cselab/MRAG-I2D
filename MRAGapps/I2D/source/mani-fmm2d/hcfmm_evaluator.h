/*
 *  hcfmm_evaluator.h
 *  hcfmm
 *
 *  Created by Manfred on 1/26/09.
 *  Copyright 2009 ETHZ. All rights reserved.
 *
 */
#include "hcfmm_box.h"

#ifndef _HCFMM_EVAL_
#define _HCFMM_EVAL_

namespace HCFMM{
	
template <typename tExpansions, typename TargetPoint, int _maxlevel>
class Evaluator
	{
		
	public:
		typedef Box<tExpansions,_maxlevel> tBox;
		Evaluator(tBox* in_rootBox):rootBox(in_rootBox)
		{
			assert (rootBox!=NULL);
			assert (rootBox->got_expansions==true);
		}
		//copy-constr.
		Evaluator(const Evaluator &E):rootBox(E.rootBox){}
		//destructor:
		virtual ~Evaluator(){rootBox=NULL;}
		//assign-operator:
		Evaluator & operator=(const Evaluator &E)
		{
			rootBox=E.rootBox;
		}
		
	protected:
		tBox* rootBox;
		
	};



}


#endif
