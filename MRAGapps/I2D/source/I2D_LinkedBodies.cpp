/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#include <limits>
#include "I2D_LinkedBodies.h"

namespace LinksStuff
{
	struct CleanTmp
	{
		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{		
			FluidElement2D * const e = &b(0,0);
				
				static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
				for(int i=0; i<n; i++) 
					e[i].tmp = 0;			
		}
	};
	
	struct FillBlocksRectangle
	{
		Real m_length;
		Real m_width;
		Real m_epsilon;
		Real m_cm[2], v1[2], v2[2], v3[2], v4[2], nor1[2], nor2[2];
		Real m_angle;
		Real rectangle_box[2][2];
		
		FillBlocksRectangle(Real length, Real width, Real epsilon,  Real cm[2], Real angle)
		{
			m_length = length;
			m_width = width;
			m_epsilon = epsilon;
			m_angle = angle;
			
			m_cm[0] = cm[0];
			m_cm[1] = cm[1];
			
			_compute_vertices(v1,v2,v3,v4,nor1,nor2,m_cm,m_length,m_width,m_angle);
			_find_rectangle_box();
		}
		
		FillBlocksRectangle(const FillBlocksRectangle& c):
		m_length(c.m_length), m_width(c.m_width), m_epsilon(c.m_epsilon), m_angle(c.m_angle)
		{
			m_cm[0] = c.m_cm[0];
			m_cm[1] = c.m_cm[1];
			
			v1[0] = c.v1[0];
			v1[1] = c.v1[1];
			v2[0] = c.v2[0];
			v2[1] = c.v2[1];
			v3[0] = c.v3[0];
			v3[1] = c.v3[1];
			v4[0] = c.v4[0];
			v4[1] = c.v4[1];
			nor1[0] = c.nor1[0];
			nor1[1] = c.nor1[1];
			nor2[0] = c.nor2[0];
			nor2[1] = c.nor2[1];
			
			rectangle_box[0][0] = c.rectangle_box[0][0];
			rectangle_box[0][1] = c.rectangle_box[0][1];
			rectangle_box[1][0] = c.rectangle_box[1][0];
			rectangle_box[1][1] = c.rectangle_box[1][1];
		}
		
		bool _is_touching(const BlockInfo& info) const
		{
			Real min_pos[2], max_pos[2];
			
			info.pos(min_pos, 0,0);
			info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1);
			
			Real intersection[2][2] = {
				max(min_pos[0], rectangle_box[0][0]), min(max_pos[0], rectangle_box[0][1]),
				max(min_pos[1], rectangle_box[1][0]), min(max_pos[1], rectangle_box[1][1]),
			};
			
			return 
			intersection[0][1]-intersection[0][0]>0 && 
			intersection[1][1]-intersection[1][0]>0 ;
		}
		
		void _find_rectangle_box()
		{
			rectangle_box[0][0] = min(min(min(v1[0],v2[0]),v3[0]),v4[0]) - m_epsilon;
			rectangle_box[0][1] = max(max(max(v1[0],v2[0]),v3[0]),v4[0]) + m_epsilon;
			rectangle_box[1][0] = min(min(min(v1[1],v2[1]),v3[1]),v4[1]) - m_epsilon;
			rectangle_box[1][1] = max(max(max(v1[1],v2[1]),v3[1]),v4[1]) + m_epsilon;
		}
		
		void _compute_norm(const Real *v2,const Real *v1,Real *nor1,Real *nor2)
		{
			nor1[0] = 1.0;
			nor1[1] = 0.0;
			Real v[2] = {0.0,0.0};
			v[0] = v2[0]-v1[0];
			v[1] = v2[1]-v1[1];
			Real mod2 = sqrt(v[0]*v[0]+v[1]*v[1]);
			nor2[0] = v[0]/mod2;
			nor2[1] = v[1]/mod2;
			if(v[0]==0.0)
			{
				nor1[0] = 1.0;
				nor1[1] = 0.0;
				return;
			}
			if(v[1]==0.0)
			{
				nor1[0] = 0.0;
				nor1[1] = 1.0;
				return;
			}
			nor1[1] = -(v[0]*nor1[0])/v[1];
			Real mod1 = sqrt(nor1[0]*nor1[0]+nor1[1]*nor1[1]);
			nor1[0] /= mod1;
			nor1[1] /= mod1;
		}
		
		void _compute_vertices(Real *v1,Real *v2,Real *v3,Real *v4,Real *nor1,Real *nor2,const Real *cm,const Real length,const Real width,const Real angle)
		{
			nor1[0] = 0.0;
			nor1[1] = 0.0;
			nor2[0] = 0.0;
			nor2[1] = 0.0;
			Real cm1[2] = {cm[0],cm[1]};
			cm1[0] += length*cos(angle);
			cm1[1] += length*sin(angle);
			_compute_norm(cm1,cm,nor1,nor2);
			v1[0] = cm[0] + nor1[0]*width/2.0;
			v1[1] = cm[1] + nor1[1]*width/2.0;
			v2[0] = cm[0] - nor1[0]*width/2.0;
			v2[1] = cm[1] - nor1[1]*width/2.0;			
			v3[0] = cm1[0] + nor1[0]*width/2.0;
			v3[1] = cm1[1] + nor1[1]*width/2.0;
			v4[0] = cm1[0] - nor1[0]*width/2.0;
			v4[1] = cm1[1] - nor1[1]*width/2.0;
		}
		
		Real _compute_char_func_box(const Real *v1,const Real *v2,const Real *v3,const Real *v4,const Real *nor1,const Real *nor2,const Real *pos,const Real length,const Real width) const
		{
			Real dist = 0.0;
			Real a[2] = {v1[0]-pos[0],v1[1]-pos[1]};
			Real b[2] = {v2[0]-pos[0],v2[1]-pos[1]};
			Real c[2] = {v2[0]-pos[0],v2[1]-pos[1]};
			Real d[2] = {v4[0]-pos[0],v4[1]-pos[1]};
			Real distA = fabs(a[0]*nor1[0]+a[1]*nor1[1]);
			Real distB = fabs(b[0]*nor1[0]+b[1]*nor1[1]);
			Real distC = fabs(c[0]*nor2[0]+c[1]*nor2[1]);
			Real distD = fabs(d[0]*nor2[0]+d[1]*nor2[1]);
			if(distA<=width && distB<=width)
			{
				dist = min(distA,distB);
			}
			else
			{
				dist = -min(distA,distB);
			}
			if(distC>length || distD>length)
			{
				dist = -1e30;
			}
			return dist;
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
						
						Real dist = _compute_char_func_box(v1,v2,v3,v4,nor1,nor2,p,m_length,m_width);
						
						b(ix,iy).tmp = max(I2D_LinkedBodies::mollified_heaviside(-dist,m_epsilon),b(ix,iy).tmp);
					}
			}
		}
	};	
	
	struct FillBlocks
	{
		const Real radius, epsilon;
		Real sphere_position[2];
		Real sphere_box[2][2];
		
		void _find_sphere_box()
		{
			sphere_box[0][0] = sphere_position[0] - radius;
			sphere_box[0][1] = sphere_position[0] + radius;
			sphere_box[1][0] = sphere_position[1] - radius;
			sphere_box[1][1] = sphere_position[1] + radius;
		}
		
		FillBlocks(Real radius, Real epsilon,  Real sphere_position[2]): 
		radius(radius), epsilon(epsilon)
		{
			this->sphere_position[0] = sphere_position[0];
			this->sphere_position[1] = sphere_position[1];
			
			_find_sphere_box();
		}
		
		FillBlocks(const FillBlocks& c): 
		radius(c.radius), epsilon(c.epsilon)
		{
			sphere_position[0] = c.sphere_position[0];
			sphere_position[1] = c.sphere_position[1];
			
			_find_sphere_box();
		}
		
		bool _is_touching(const BlockInfo& info) const
		{
			Real min_pos[2], max_pos[2];
			
			info.pos(min_pos, 0,0);
			info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1);
			
			Real intersection[2][2] = {
				max(min_pos[0], sphere_box[0][0]), min(max_pos[0], sphere_box[0][1]),
				max(min_pos[1], sphere_box[1][0]), min(max_pos[1], sphere_box[1][1]),
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
						
						const Real r = sqrt(pow(p[0]-sphere_position[0],2) + 
											pow(p[1]-sphere_position[1],2) );
						
						b(ix, iy).tmp = max(I2D_LinkedBodies::mollified_heaviside(r-radius, epsilon), b(ix, iy).tmp);
					}
			}		
		}
	};
}

Real I2D_LinkedBodies::mollified_heaviside(const Real x, const Real eps)
{
	const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));
	
	return 0.5+0.5*cos(alpha);
}

void I2D_LinkedBodies::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	Real CM[2] = {position[0], position[1]};
	Real cm[2] = {position[0], position[1]};
		
	vector<Real> anglesPlusOne;
	
	for(unsigned int i=0; i<angles.size(); i++)
		anglesPlusOne.push_back(angles[i]);
	
	anglesPlusOne.push_back(0.0);
	
	for(int i=0;i<(int)anglesPlusOne.size();i++){
		if(anglesPlusOne[i]<0.0){anglesPlusOne[i]+=360.0;}
		anglesPlusOne[i]*=(M_PI/180.0); 
	}	
	//for(int i=0;i<anglesPlusOne.size();i++){ anglesPlusOne[i]*=(M_PI/180.0); }
	
	LinksStuff::CleanTmp clean;
	block_processing.process(vInfo, coll, clean);
	
	for(int k=0;k<(int)anglesPlusOne.size();k++)
	{		
		LinksStuff::FillBlocks fill(width/2.0, epsilon, cm);
		block_processing.process(vInfo, coll, fill);
		cm[0] += length*cos(anglesPlusOne[k]);
		cm[1] += length*sin(anglesPlusOne[k]);
	}
	
	cm[0] = CM[0];
	cm[1] = CM[1];
	
	for(int k=0;k<(int)anglesPlusOne.size()-1;k++)
	{
		LinksStuff::FillBlocksRectangle fillRectangle(length, width, epsilon, cm, anglesPlusOne[k]);
		block_processing.process(vInfo, coll, fillRectangle);
	
		cm[0] += length*cos(anglesPlusOne[k]);
		cm[1] += length*sin(anglesPlusOne[k]);
	}	
}

