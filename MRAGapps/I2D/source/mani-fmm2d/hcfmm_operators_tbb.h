/*
 *  hcfmm_operators_tbb.h
 *  hcfmm
 *
 *  Created by Manfred on 12/23/08.
 *  Copyright 2008 ETHZ. All rights reserved.
 *
 */

#include "hcfmm_types.h"

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>


#pragma once
template <class BaseType, class Particle, int sDim>
struct getBoundingBox_TBB
{
	Particle* vParticles;
	bbox<BaseType,sDim> curBbox;
	int nParticles;
	
	getBoundingBox_TBB(Particle* _vParticles, bbox<BaseType,sDim> _BBox, int _nParticles):
	vParticles(_vParticles),
	curBbox(_BBox),
	nParticles(_nParticles)
	{}
	
	getBoundingBox_TBB(const getBoundingBox_TBB& body)
	{
		vParticles=body.vParticles;
		nParticles=body.nParticles;
		curBbox=body.curBbox;
	}
	
	void operator()(const tbb::blocked_range<int> & range)
    {
		Particle* work_Particles=vParticles;
		bbox<BaseType,sDim> work_bbox=curBbox;
		int iStart=range.begin();
		int iEnd=range.end();
		for ( int i=iStart;i!=iEnd;++i)
		{
			for (int d=0; d<sDim;++d)
			{
				work_bbox.upper[d]=max(work_bbox.upper[d],work_Particles[i].x[d]);
				work_bbox.lower[d]=min(work_bbox.lower[d],work_Particles[i].x[d]);
			}
		}
		
		//write back:
		curBbox=work_bbox;
		
	}
	
	
	void print() {curBbox.print();}
	
	getBoundingBox_TBB(getBoundingBox_TBB& x, tbb::split) : nParticles(x.nParticles), vParticles(x.vParticles), curBbox(x.curBbox) {}
	
	void join (const getBoundingBox_TBB& y ) {curBbox.compare(y.curBbox);}
	
};