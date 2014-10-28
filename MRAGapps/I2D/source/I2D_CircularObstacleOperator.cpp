/*
 *  I2D_SphereObstacleOperator.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include <limits>
#include "I2D_CircularObstacleOperator.h"

Real I2D_CircularObstacleOperator::_get_smooth_radius(const Real R, const Real eps)
{
	const double D = 2*R;
	const double a = M_PI*D*D/4;
	const double num = sqrt(4*a*M_PI - eps*eps*(M_PI*M_PI-8));
	const double denom = 2*M_PI;
	
	return num/denom;
}

struct FillBlocks
{
	const Real smooth_radius, smoothing_length, support_radius;
	Real sphere_position[2];
	Real sphere_box[2][2];
	
	void _find_sphere_box()
	{
		sphere_box[0][0] = sphere_position[0] - support_radius;
		sphere_box[0][1] = sphere_position[0] + support_radius;
		sphere_box[1][0] = sphere_position[1] - support_radius;
		sphere_box[1][1] = sphere_position[1] + support_radius;
	}
	
	FillBlocks(Real support_radius, Real smooth_radius, Real smoothing_length,  Real sphere_position[3]): 
	support_radius(support_radius), smooth_radius(smooth_radius), smoothing_length(smoothing_length)
	{
		this->sphere_position[0] = sphere_position[0];
		this->sphere_position[1] = sphere_position[1];
		
		_find_sphere_box();
	}
	
	FillBlocks(const FillBlocks& c): 
	support_radius(c.support_radius), smooth_radius(c.smooth_radius), smoothing_length(c.smoothing_length)
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
					
					b(ix, iy).tmp = I2D_CircularObstacleOperator::mollified_heaviside(r-smooth_radius, smoothing_length);
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

Real I2D_CircularObstacleOperator::mollified_heaviside(const Real x, const Real eps)
{
	const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));
	
	return 0.5+0.5*cos(alpha);
}
/*
void I2D_CircularObstacleOperator::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	FillBlocks fill(smooth_radius+smoothing_length, smooth_radius, smoothing_length, position);
	block_processing.process(vInfo, coll, fill);
}
*/

namespace myass
{
double getHeavisideFDMH1(double * temp, const double h)
{
  const double dist = temp[1];
  if(dist >= +h)
    return 1;

  if(dist <= -h) return 0;
  assert(std::abs(dist)<=h);

  // compute first primitive of H(x): I(x) = int_0^x H(y) dy                                                                                                                        \
                                                                                                                                                                                     
  Real IplusX = temp[2];
  Real IminuX = temp[0];
  Real IplusY = temp[5];
  Real IminuY = temp[3];

  // set it to zero outside the cylinder                                                                                                                                             

  IplusX = IplusX < 0 ? 0 : IplusX;
  IminuX = IminuX < 0 ? 0 : IminuX;
  IplusY = IplusY < 0 ? 0 : IplusY;
  IminuY = IminuY < 0 ? 0 : IminuY;

  assert(IplusX>=0);
  assert(IminuX>=0);
  assert(IplusY>=0);
  assert(IminuY>=0);

  // gradI                                                                                                                                                                          \
                                                                                                                                                                                     
  const Real gradIX = 0.5/h * (IplusX - IminuX);
  const Real gradIY = 0.5/h * (IplusY - IminuY);

  // gradU                                                                                                                                                                          \
                                                                                                                                                                                     
  const Real gradUX = 0.5/h * (temp[2]-temp[0]);
  const Real gradUY = 0.5/h * (temp[5]-temp[3]);

  const Real denom = gradUX*gradUX+gradUY*gradUY;
  const Real numer = gradIX*gradUX + gradIY*gradUY;
  const Real H = denom==0 ?  numer : numer/denom;

  assert(H>=0 && H<=1);

  return H;
}
};

double sample(const double x, const double y, const double h, const double position[2], const double radius)
{
  const int stencil = 6;
  double temp[stencil];
  for(int i=0; i<stencil; i++)
    temp[i]=0;

  int idy=0;
    for(int idx=-1;idx<2;idx++)
    {
      const double my_x = x+idx*h;
      const double my_y = y+idy*h;

      const Real dist = sqrt(pow(my_x-position[0],2) +
			     pow(my_y-position[1],2) ) - radius;

      temp[idx+1] -= dist;
    }

    int idx=0;
    for(int idy=-1;idy<2;idy++)
      {
	if (idy==0) continue;

	const double my_x = x+idx*h;
	const double my_y = y+idy*h;

	const Real dist = sqrt(pow(my_x-position[0],2) +
			       pow(my_y-position[1],2) ) - radius;

	temp[idy+4] -= dist;
      }

    return myass::getHeavisideFDMH1(temp, h);
}

struct FillBlocksTowers
{
  const Real smooth_radius, smoothing_length, support_radius;
  Real sphere_position[2];
  Real sphere_box[2][2];

  void _find_sphere_box()
  {
    sphere_box[0][0] = sphere_position[0] - support_radius;
    sphere_box[0][1] = sphere_position[0] + support_radius;
    sphere_box[1][0] = sphere_position[1] - support_radius;
    sphere_box[1][1] = sphere_position[1] + support_radius;
  }

  FillBlocksTowers(Real support_radius, Real smooth_radius, Real smoothing_length,  Real sphere_position[2]):
    support_radius(support_radius), smooth_radius(smooth_radius), smoothing_length(smoothing_length)
  {
    this->sphere_position[0] = sphere_position[0];
    this->sphere_position[1] = sphere_position[1];

    _find_sphere_box();
  }

  FillBlocksTowers(const FillBlocksTowers& c):
    support_radius(c.support_radius), smooth_radius(c.smooth_radius), smoothing_length(c.smoothing_length)
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
    //    bool bEmpty = true;

    if(_is_touching(info))
      {
        for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
          for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
            {
              double p[2];
              info.pos(p, ix, iy);
	      //	      const Real r = sqrt(pow(p[0]-sphere_position[0],2) +
	      //				  pow(p[1]-sphere_position[1],2) );

	      //b(ix, iy).tmp = I2D_CircularObstacleOperator::mollified_heaviside(r-smooth_radius, smoothing_length);

	      b(ix,iy).tmp = sample(p[0], p[1], info.h[0], sphere_position, smooth_radius);
            }

	//        bEmpty = false;
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

void I2D_CircularObstacleOperator::characteristic_function()
{
  vector<BlockInfo> vInfo = grid.getBlocksInfo();
  const BlockCollection<B>& coll = grid.getBlockCollection();

  if(!SHARP)
    {
      FillBlocks fill(smooth_radius+smoothing_length, smooth_radius, smoothing_length, position);                                                                                  
      block_processing.process(vInfo, coll, fill);
    }
  else
    {
      FillBlocksTowers fill(smooth_radius+smoothing_length, smooth_radius, smoothing_length, position);
      block_processing.process(vInfo, coll, fill);
    }
}
