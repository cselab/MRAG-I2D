/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include <limits>
#include "I2D_KillVortRightBoundaryOperator.h"

namespace KillStuff
{	
	struct FillBlocks2
	{
		const int killed_width;
		const Real max_dx;
		
		FillBlocks2(int killed_width, Real max_dx): 
		killed_width(killed_width), max_dx(max_dx)
		{
		}
		
		FillBlocks2(const FillBlocks2& c): 
		killed_width(c.killed_width),max_dx(c.max_dx)
		{
		}
		
		bool _is_touching(const BlockInfo& info) const // true: block within killing zone, false: block not within killing zone
		{
			Real max_pos[2];			
			info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1); // max position in physical coordinates
			
			
			return ( max_pos[0]>(1.0-(5+killed_width)*max_dx));
		}
		
		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{	
			if(_is_touching(info))
			{				
				Real factor = 0.0;
				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						Real p[2];
						info.pos(p, ix, iy);						
						factor = cos( (M_PI/(2*killed_width*max_dx))*max(0.0,(p[0]-1.0+(5+killed_width)*max_dx)) );
						factor = max(Real(0.0),factor); // smooth within killing zone (factor <= 1) and kill at very boundaries (factor < 0)
						b(ix, iy).omega = b(ix,iy).omega*factor;
					}				
			}		
		}
	};
}


void I2D_KillVortRightBoundaryOperator::killVorticity()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();
	const Real max_dx= (1./B::sizeX)*pow(0.5,grid.getCurrentMinLevel());;
	
	KillStuff::FillBlocks2 fill(killed_width,max_dx);
	block_processing.process(vInfo, coll, fill);
}