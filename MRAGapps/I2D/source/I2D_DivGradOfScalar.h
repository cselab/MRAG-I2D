/*
 *  I2D_DivGradOfScalar_2ndOrder.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once
struct I2D_DivGradOfScalar_2ndOrder
{	
	template<int direction, int BSsize>
	void _get_regions(const BlockInfo& info, int (region[3])[2]) const
	{
		const int max_index = (1<<info.level)-1;
		
		const bool touching_low = info.index[direction]==0;
		const bool touching_high = info.index[direction]==max_index;
		
		region[0][0] = 0;	
		region[0][1] = touching_low;
		
		region[1][0] = touching_low;	
		region[1][1] = BSsize - touching_high;
		
		region[2][0] = BSsize - touching_high; 
		region[2][1] = BSsize;
	}
	
	template<typename streamer, typename Lab>
	void _ddpsiddx(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const
	{
		const Real factor = 1/pow(info.h[0],2);
		
		int region[3][2];
		_get_regions<0, B::sizeX>(info, region);
		
		if (region[0][1]>region[0][0])
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=region[0][0]; ix<region[0][1]; ix++)
					streamer::stream(out(ix, iy), factor*(lab(ix, iy) - 2*lab(ix+1, iy) + lab(ix+2, iy)));
		
		for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=region[1][0]; ix<region[1][1]; ix++)
				streamer::stream(out(ix, iy), factor*(lab(ix+1, iy) -2*lab(ix, iy) + lab(ix-1, iy)));
		
		if (region[2][1]>region[2][0])
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=region[2][0]; ix<region[2][1]; ix++)
					streamer::stream(out(ix, iy), factor*(lab(ix, iy) - 2*lab(ix-1, iy) + lab(ix-2, iy)));
	}
	
	template<typename streamer, typename Lab>
	void _ddpsiddy(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const
	{
		const Real factor = 1/pow(info.h[0],2);
		
		int region[3][2];
		_get_regions<1, B::sizeY>(info, region);
		
		if (region[0][1]>region[0][0])
			for(int iy=region[0][0]; iy<region[0][1]; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
					streamer::stream(out(ix, iy), factor*(lab(ix, iy) - 2*lab(ix, iy+1) + lab(ix, iy+2)));
		
		for(int iy=region[1][0]; iy<region[1][1]; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
				streamer::stream(out(ix, iy), factor*(lab(ix, iy+1) -2*lab(ix, iy) + lab(ix, iy-1)));
		
		if (region[2][1]>region[2][0])
			for(int iy=region[2][0]; iy<region[2][1]; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
					streamer::stream(out(ix, iy), factor*(lab(ix, iy) - 2*lab(ix, iy-1) + lab(ix, iy-2)));		
	}
};

struct I2D_DivGradOfScalar_4thOrder
{	
	template<int direction, int BSsize>
	void _get_regions(const BlockInfo& info, int (region[5])[2]) const
	{
		const int max_index = (1<<info.level)-1;
		
		const bool touching_low = info.index[direction]==0;
		const bool touching_high = info.index[direction]==max_index;
		
		for(int i=0; i<5; i++)
			region[i][0] = region[i][1] = 0;
		
		region[2][1] = BSsize;
		
		if (!touching_low && !touching_high) return;
		
		if (touching_low) 
		{
			region[0][0] = 0;	region[0][1] = 1;
			region[1][0] = 1;	region[1][1] = 2;
			region[2][0] = 2;
		}
		
		if (touching_high) 
		{
			region[2][1] = BSsize-2;
			region[3][0] = BSsize-2;	region[3][1] = BSsize-1;
			region[4][0] = BSsize-1;	region[4][1] = BSsize;
		}
	}
	
	template<typename streamer, typename Lab>
	void _ddpsiddx(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const
	{
		const Real factor = 1./(pow(info.h[0],2)*12);
		
		int region[5][2];
		_get_regions<0, B::sizeX>(info, region);
		
		if (region[0][1]>region[0][0])
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=region[0][0]; ix<region[0][1]; ix++)
					streamer::stream(out(ix, iy), factor*(35*lab(ix, iy) - 104*lab(ix+1, iy) + 114*lab(ix+2, iy) - 56*lab(ix+3, iy) + 11*lab(ix+4, iy)));
		
		if (region[1][1]>region[1][0])
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=region[1][0]; ix<region[1][1]; ix++)
					streamer::stream(out(ix, iy), factor*(11*lab(ix-1, iy) -20*lab(ix, iy) + 6*lab(ix+1, iy) +4*lab(ix+2, iy) - lab(ix+3, iy)));
		
		for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=region[2][0]; ix<region[2][1]; ix++)
					streamer::stream(out(ix, iy), factor*(-lab(ix-2, iy) + 16*lab(ix-1, iy) - 30*lab(ix, iy) + 16*lab(ix+1, iy) - lab(ix+2, iy)));
		
		if (region[3][1]>region[3][0])
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=region[3][0]; ix<region[3][1]; ix++)
					streamer::stream(out(ix, iy), factor*(-lab(ix-3, iy) + 4*lab(ix-2, iy) + 6*lab(ix-1, iy) - 20*lab(ix, iy) + 11*lab(ix+1, iy)));
		
		if (region[4][1]>region[4][0])
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=region[4][0]; ix<region[4][1]; ix++)
					streamer::stream(out(ix, iy), factor*(11*lab(ix-4, iy) - 56*lab(ix-3, iy) + 114*lab(ix-2, iy) - 104*lab(ix-1, iy) + 35*lab(ix, iy)));
	}
	
	template<typename streamer, typename Lab>
	void _ddpsiddx_ptr(Lab& lab, const BlockInfo& info, Real *  const ptr) const
	{
		const Real factor = 1./(pow(info.h[0],2)*12);
		
		int region[5][2];
		_get_regions<0, B::sizeX>(info, region);
		
		if (region[0][1]>region[0][0])
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=region[0][0]; ix<region[0][1]; ix++)
					streamer::stream(ptr[iy*B::sizeX + ix], factor*(35*lab(ix, iy) - 104*lab(ix+1, iy) + 114*lab(ix+2, iy) - 56*lab(ix+3, iy) + 11*lab(ix+4, iy)));
		
		if (region[1][1]>region[1][0])
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=region[1][0]; ix<region[1][1]; ix++)
					streamer::stream(ptr[iy*B::sizeX + ix], factor*(11*lab(ix-1, iy) -20*lab(ix, iy) + 6*lab(ix+1, iy) +4*lab(ix+2, iy) - lab(ix+3, iy)));
		
		for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=region[2][0]; ix<region[2][1]; ix++)
				streamer::stream(ptr[iy*B::sizeX + ix], factor*(-lab(ix-2, iy) + 16*lab(ix-1, iy) - 30*lab(ix, iy) + 16*lab(ix+1, iy) - lab(ix+2, iy)));
		
		if (region[3][1]>region[3][0])
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=region[3][0]; ix<region[3][1]; ix++)
					streamer::stream(ptr[iy*B::sizeX + ix], factor*(-lab(ix-3, iy) + 4*lab(ix-2, iy) + 6*lab(ix-1, iy) - 20*lab(ix, iy) + 11*lab(ix+1, iy)));
		
		if (region[4][1]>region[4][0])
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=region[4][0]; ix<region[4][1]; ix++)
					streamer::stream(ptr[iy*B::sizeX + ix], factor*(11*lab(ix-4, iy) - 56*lab(ix-3, iy) + 114*lab(ix-2, iy) - 104*lab(ix-1, iy) + 35*lab(ix, iy)));
	}
	
	template<typename streamer, typename Lab>
	void _ddpsiddy(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const
	{
		const Real factor = 1./(pow(info.h[0],2)*12);
		
		int region[5][2];
		_get_regions<1, B::sizeY>(info, region);
		
		if (region[0][1]>region[0][0])
			for(int iy=region[0][0]; iy<region[0][1]; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
					streamer::stream(out(ix, iy), factor*(35*lab(ix, iy) - 104*lab(ix, iy+1) + 114*lab(ix, iy+2) - 56*lab(ix, iy+3) + 11*lab(ix, iy+4)));
		
		if (region[1][1]>region[1][0])
			for(int iy=region[1][0]; iy<region[1][1]; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
					streamer::stream(out(ix, iy), factor*(11*lab(ix, iy-1) -20*lab(ix, iy) + 6*lab(ix, iy+1) +4*lab(ix, iy+2) - lab(ix, iy+3)));
		
		for(int iy=region[2][0]; iy<region[2][1]; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
				streamer::stream(out(ix, iy), factor*(-lab(ix, iy-2) + 16*lab(ix, iy-1) - 30*lab(ix, iy) + 16*lab(ix, iy+1) - lab(ix, iy+2)));
		
		if (region[3][1]>region[3][0])
			for(int iy=region[3][0]; iy<region[3][1]; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
					streamer::stream(out(ix, iy), factor*(-lab(ix, iy-3) + 4*lab(ix, iy-2) + 6*lab(ix, iy-1) - 20*lab(ix, iy) + 11*lab(ix, iy+1)));
		
		if (region[4][1]>region[4][0])
			for(int iy=region[4][0]; iy<region[4][1]; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
					streamer::stream(out(ix, iy), factor*(11*lab(ix, iy-4) - 56*lab(ix, iy-3) + 114*lab(ix, iy-2) - 104*lab(ix, iy-1) + 35*lab(ix, iy)));
	}
	
	template<typename streamer, typename Lab>
	void _ddpsiddy_ptr(Lab& lab, const BlockInfo& info, Real * const ptr) const
	{
		const Real factor = 1./(pow(info.h[0],2)*12);
		
		int region[5][2];
		_get_regions<1, B::sizeY>(info, region);
		
		if (region[0][1]>region[0][0])
			for(int iy=region[0][0]; iy<region[0][1]; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
					streamer::stream(ptr[iy*B::sizeX + ix], factor*(35*lab(ix, iy) - 104*lab(ix, iy+1) + 114*lab(ix, iy+2) - 56*lab(ix, iy+3) + 11*lab(ix, iy+4)));
		
		if (region[1][1]>region[1][0])
			for(int iy=region[1][0]; iy<region[1][1]; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
					streamer::stream(ptr[iy*B::sizeX + ix], factor*(11*lab(ix, iy-1) -20*lab(ix, iy) + 6*lab(ix, iy+1) +4*lab(ix, iy+2) - lab(ix, iy+3)));
		
		for(int iy=region[2][0]; iy<region[2][1]; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
				streamer::stream(ptr[iy*B::sizeX + ix], factor*(-lab(ix, iy-2) + 16*lab(ix, iy-1) - 30*lab(ix, iy) + 16*lab(ix, iy+1) - lab(ix, iy+2)));
		
		if (region[3][1]>region[3][0])
			for(int iy=region[3][0]; iy<region[3][1]; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
					streamer::stream(ptr[iy*B::sizeX + ix], factor*(-lab(ix, iy-3) + 4*lab(ix, iy-2) + 6*lab(ix, iy-1) - 20*lab(ix, iy) + 11*lab(ix, iy+1)));
		
		if (region[4][1]>region[4][0])
			for(int iy=region[4][0]; iy<region[4][1]; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
					streamer::stream(ptr[iy*B::sizeX + ix], factor*(11*lab(ix, iy-4) - 56*lab(ix, iy-3) + 114*lab(ix, iy-2) - 104*lab(ix, iy-1) + 35*lab(ix, iy)));
	}
};