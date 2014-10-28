/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include "I2D_Clear.h"
#include "I2D_Headers.h"
#include "I2D_Types.h"

struct I2D_TmpToZero
{	
	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{	
		FluidElement2D * const e = &b(0,0);
		static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
		for(int i=0; i<n; i++) 
			e[i].tmp = 0.0;
	}
};

struct I2D_VelToZero
{
	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{
		FluidElement2D * const e = &b(0,0);
		static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
		for(int i=0; i<n; i++)
		{
			e[i].u[0] = 0.0;
			e[i].u[1] = 0.0;
		}
	}
};

void I2D_Clear::clearTmp(Grid<W,B> & grid)
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	I2D_TmpToZero clean;
	block_processing.process(vInfo, coll, clean);
}

void I2D_Clear::clearVel(Grid<W,B> & grid)
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	I2D_VelToZero clean;
	block_processing.process(vInfo, coll, clean);
}

