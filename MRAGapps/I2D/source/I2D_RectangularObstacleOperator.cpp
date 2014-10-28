/*
 *  I2D_RectangularObstacleOperator.cpp
 *  I2D_ROCKS
 *
 *  Created by Diego Rossinelli on 2/2/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include "I2D_RectangularObstacleOperator.h"

#include <limits>

Real I2D_RectangularObstacleOperator::mollified_heaviside(const Real x, const Real eps)
{
	const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));
	
	return 0.5+0.5*cos(alpha);
}

inline double area(double LX, double LY, double eps0, double eps1)
{
	const double Asmooth1 = LX*LY + 2*eps0*(LX+LY) + M_PI*eps0*eps0; 
	const double Asmooth2 = (LX+LY)*(eps1-eps0) + 0.5*M_PI*(eps1-eps0)*(eps1+eps0) - 2./M_PI*pow(eps0-eps1,2);
	
	return Asmooth1+Asmooth2;
}

Real I2D_RectangularObstacleOperator::_get_eps1()
{	
	const double LX = aspect_ratioX*D-2*eps0;
	const double LY = D-2*eps0;
	const double desired_area = aspect_ratioX*D*D;
	
	double l = eps0;
	double r = 2*eps0;
	
	assert(area(LX, LY, eps0, l)<desired_area);
	assert(area(LX, LY, eps0, r)>=desired_area);
	
	int counter = 0;
	while(true)
	{
		const double rel_diff = fabs((area(LX, LY, eps0, l)-desired_area)/desired_area);		
	
		printf("iter %d: %e -> eps0 %e eps1 %e \n", counter, rel_diff, eps0, l);
		
		if(rel_diff < 1e-13) return l;
		
		double m = l*0.5 + r*0.5;
		
		if(area(LX, LY, eps0,m)<desired_area)
			l=m;
		else
			r=m;
		
		counter++;
	}
	
	abort();
	return 0;
}


struct FillBlocksPlate
{
	const Real LX, LY, eps0, eps1;
	Real position[2];
	Real bbox[2][2];
	
	void _find_bbox()
	{
		const double WX = LX/2 + eps0 + eps1;
		const double WY = LY/2 + eps0 + eps1;
		
		bbox[0][0] = position[0] - WX;
		bbox[0][1] = position[0] + WX;
		bbox[1][0] = position[1] - WY;
		bbox[1][1] = position[1] + WY;
	}
	
	FillBlocksPlate(Real LX, Real LY, Real eps0, Real eps1, Real position[2]): 
	LX(LX), LY(LY), eps0(eps0), eps1(eps1)
	{
		this->position[0] = position[0];
		this->position[1] = position[1];
		
		_find_bbox();
	}
	
	FillBlocksPlate(const FillBlocksPlate& c): 
	LX(c.LX), LY(c.LY), eps0(c.eps0), eps1(c.eps1)
	{
		position[0] = c.position[0];
		position[1] = c.position[1];
		
		_find_bbox();
	}
	
	bool _is_touching(const BlockInfo& info) const
	{
		Real min_pos[2], max_pos[2];
		
		info.pos(min_pos, 0,0);
		info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1);
		
		Real intersection[2][2] = {
			max(min_pos[0], bbox[0][0]), min(max_pos[0], bbox[0][1]),
			max(min_pos[1], bbox[1][0]), min(max_pos[1], bbox[1][1]),
		};
		
		return 
		intersection[0][1]-intersection[0][0]>0 && 
		intersection[1][1]-intersection[1][0]>0 ;
	}
	
	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{		
		if(_is_touching(info))
		{
			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
				{
					Real p[2];
					info.pos(p, ix, iy);
					
					const Real pclst[2] = {
						max(position[0]-LX/2, min(position[0]+LX/2, p[0])),
						max(position[1]-LY/2, min(position[1]+LY/2, p[1]))
					};
					
					
					const Real r = sqrt(pow(p[0]-pclst[0],2) + 
										pow(p[1]-pclst[1],2) );
					
					b(ix, iy).tmp = I2D_RectangularObstacleOperator::mollified_heaviside(r-(eps0+eps1)*0.5, eps1-eps0);
				}
		}
		else 
		{
			FluidElement2D * const e = &b(0,0);
			
			static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
			for(int i=0; i<n; i++) 
				e[i].tmp = 0;
		}
		
	}
};

void I2D_RectangularObstacleOperator::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	const double LX = D*aspect_ratioX - 2*eps0;
	const double LY = D - 2*eps0;
	FillBlocksPlate fill(LX, LY, eps0, eps1*2, position);
	block_processing.process(vInfo, coll, fill);
}