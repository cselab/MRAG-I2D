/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include <limits>
#include "I2D_EllipticalObstacleOperator.h"

namespace EllyStuff
{
	
	struct FillBlocksElly
	{
		const Real semiMajorAxis, aspectRatio, angle, epsilon;
		Real position[2];
		Real box[2][2];
		
		void _rotate(Real v[2], Real angle)
		{
			Real a00 = cos(angle);
			Real a01 = -sin(angle);
			Real a10 = sin(angle);
			Real a11 = cos(angle);
			
			Real vx = v[0];
			Real vy = v[1];
			
			v[0] = a00*(vx-position[0]) + a01*(vy-position[1]);
			v[1] = a10*(vx-position[0]) + a11*(vy-position[1]);		
			
			v[0] += position[0];
			v[1] += position[1];
		}
		
		void _find_sphere_box()
		{
			Real semiMinorAxis = semiMajorAxis*aspectRatio;
			
			Real minX = position[0] - semiMajorAxis;
			Real maxX = position[0] + semiMajorAxis;		
			Real minY = position[1] - semiMinorAxis;
			Real maxY = position[1] + semiMinorAxis;
			
			Real v1[2] = { minX, minY };
			Real v2[2] = { maxX, maxY };
			Real v3[2] = { minX, maxY };
			Real v4[2] = { maxX, minY };
			
			_rotate(v1,angle);
			_rotate(v2,angle);
			_rotate(v3,angle);
			_rotate(v4,angle);
			
			box[0][0] = min((Real)min((Real)min(v1[0],v2[0]),v3[0]),v4[0]); // min x
			box[0][1] = max((Real)max((Real)max(v1[0],v2[0]),v3[0]),v4[0]); // max x
			box[1][0] = min((Real)min((Real)min(v1[1],v2[1]),v3[1]),v4[1]); // min y
			box[1][1] = max((Real)max((Real)max(v1[1],v2[1]),v3[1]),v4[1]); // max y
		}
		
		FillBlocksElly(Real semiMajorAxis, Real aspectRatio, Real angle, Real epsilon,  Real position[2]): 
		semiMajorAxis(semiMajorAxis), aspectRatio(aspectRatio), angle(angle), epsilon(epsilon)
		{
			this->position[0] = position[0];
			this->position[1] = position[1];
			
			_find_sphere_box();
		}
		
		FillBlocksElly(const FillBlocksElly& c): 
		semiMajorAxis(c.semiMajorAxis), aspectRatio(c.aspectRatio), angle(c.angle), epsilon(c.epsilon)
		{
			position[0] = c.position[0];
			position[1] = c.position[1];
			
			_find_sphere_box();
		}
		
		bool _is_touching(const BlockInfo& info) const
		{
			Real min_pos[2], max_pos[2];
			
			info.pos(min_pos, 0,0);
			info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1);
			
			Real intersection[2][2] = {
				max(min_pos[0], box[0][0]), min(max_pos[0], box[0][1]),
				max(min_pos[1], box[1][0]), min(max_pos[1], box[1][1]),
			};
			
			return 
			intersection[0][1]-intersection[0][0]>0 && 
			intersection[1][1]-intersection[1][0]>0 ;
		}
		
		inline void operator()(const BlockInfo& info, FluidBlock2D& o) const
		{				
			if(_is_touching(info))
			{
				const Real a = semiMajorAxis;
				const Real b = semiMajorAxis*aspectRatio;
				Real dx = 0.0;
				Real dy = 0.0;
				Real alpha = 0.0;
				Real radius = 0.0;
				Real dist = 0.0;
				
				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						Real p[2];
						info.pos(p, ix, iy);
						
						dx = p[0]-position[0];
						dy = p[1]-position[1];
						
						alpha = atan2(dy,dx) - angle;
						radius = a*b/sqrt(b*b*cos(alpha)*cos(alpha) + a*a*sin(alpha)*sin(alpha));	
						dist = sqrt( dx*dx + dy*dy ) - radius;
						
						o(ix,iy).tmp = I2D_EllipticalObstacleOperator::mollified_heaviside(dist, epsilon);
					}
			}
			else 
			{
				FluidElement2D * const e = &o(0,0);
				
				static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
				for(int i=0; i<n; i++) 
					e[i].tmp = 0.0;
			}
		}
	};
	
}


Real I2D_EllipticalObstacleOperator::mollified_heaviside(const Real x, const Real eps)
{
	const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));
	
	return 0.5+0.5*cos(alpha);
}

void I2D_EllipticalObstacleOperator::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	EllyStuff::FillBlocksElly fill(semiMajorAxis, aspectRatio, angle, epsilon, position);
	block_processing.process(vInfo, coll, fill);
}


