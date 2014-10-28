/*
 *  I2D_FMMTypes.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"

#include <xmmintrin.h>

struct SourceParticlesInfo
{
	int nsource_particles;
	int start;
};

struct VelocityRHS
{
	static const int dim=2;
	Real x[2];
	
	VelocityRHS()
	{
		x[0] = x[1] = 0.0;
	}
};

struct VelocitySourceParticle
{
	typedef VelocityRHS RHSType;
	typedef Real BaseType;
	
	static const int dim = 2;
	static const int pdim = 1;
	
	Real x[2];
	Real w[1];
	
	VelocitySourceParticle()
	{
		x[0] = x[1] = 0.0;
		w[0] = 0.0;
	}
};

struct UpdateBlocks
{
	VelocityBlock * myvelblocks;
	vector<FluidBlock2D *> blocks;
	
	UpdateBlocks(){}
	UpdateBlocks(const UpdateBlocks& c): myvelblocks(c.myvelblocks), blocks(c.blocks){}
	
	void operator()(blocked_range<int> range) const
	{
		for(int iblock=range.begin(); iblock!=range.end(); ++iblock)
		{
			FluidElement2D * const dst = &(*blocks[iblock])(0,0);

			const Real * const srcu = (const Real *)myvelblocks[iblock].u[0];
			const Real * const srcv = (const Real *)myvelblocks[iblock].u[1];
			
			static const int BS = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
			for(int i=0; i<BS; i++)
			{
				dst[i].u[0] = srcu[i];
				dst[i].u[1] = srcv[i];
			}
		}
	}
};

struct UpdateBlocksPot
{
	VelocityBlock * myvelblocks;
	vector<FluidBlock2D *> blocks;

	UpdateBlocksPot(){}
	UpdateBlocksPot(const UpdateBlocksPot& c): myvelblocks(c.myvelblocks), blocks(c.blocks){}

	void operator()(blocked_range<int> range) const
	{
		for(int iblock=range.begin(); iblock!=range.end(); ++iblock)
		{
			FluidElement2D * const dst = &(*blocks[iblock])(0,0);

			const Real * const srcu = (const Real *)myvelblocks[iblock].u[0];
			const Real * const srcv = (const Real *)myvelblocks[iblock].u[1];

			static const int BS = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
			for(int i=0; i<BS; i++)
			{
				dst[i].u[0] += srcv[i];
				dst[i].u[1] -= srcu[i];

				//dst[i].u[0] = srcv[i];
				//dst[i].u[1] = -srcu[i];
			}
		}
	}
};
