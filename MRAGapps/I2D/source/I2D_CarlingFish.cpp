/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include "I2D_CarlingFish.h"
#include "I2D_Interpolation1D.h"
#include "I2D_Frenet2D.h"
#include "I2D_VectorBlockLab.h"
#include "I2D_GradOfVector.h"
#include <limits>
#include <omp.h>

I2D_CarlingFish::Fish::Fish(double xm, double ym, double _D, double _T, double phase, double angle_rad, double angleInSpace_rad, double eps, const int LMAX, const bool isSharp):
angle(angle_rad), angleInSpace(angleInSpace_rad), D(_D), xm(xm), ym(ym), phase(phase), angular_velocity(0), vx(0), vy(0), vdefx(0), vdefy(0), N(0), SIZEX(0), SIZEY(0), MAPSIZEX(0), MAPSIZEY(0),
m(1.0), J(1.0), WH(0), SB(0), WT(0), ST(0), T(_T), EXTENSION(4.0), EPS(0), DS( 0), traslX(0.0), traslY(0.0), angular_velocityInSpace(0.0), TOrig(_T), LMAX(LMAX), SHARP(isSharp)
{
	// Grid sizes adapt to max MRAG resolution
	const int RATIOXYFISH = 4;
	int MAXRESFISH = ceil( D * (double)(FluidBlock2D::sizeX*pow(2.0,LMAX)) );
	MAXRESFISH += 2*RATIOXYFISH - (MAXRESFISH % (2*RATIOXYFISH));
	MAPSIZEX = MAXRESFISH;
	MAPSIZEY = MAPSIZEX;
	
	// Adaptive
	N = 1 + MAXRESFISH;
	SIZEX = N;
	SIZEY = 1 + (int)(((double)SIZEX-1.0)/(double)RATIOXYFISH);

	// Fixed model
	//N = 2401;
	//SIZEX = 601;
	//SIZEY = 121;
	//MAPSIZEX = 512;
	//MAPSIZEY = 512;
	//N = 2401;
	//SIZEX = 2401;
	//SIZEY = 481;
	//MAPSIZEX = 2048;
	//MAPSIZEY = 2048;

	printf("----------------------------------------\n");
	printf("FISH DISCRETIZATION SUMMARY\n");
	printf("\tMAXRESFISH=%d\n",MAXRESFISH);
	printf("\tN=%d\n",N);
	printf("\tSIZEX=%d\n",SIZEX);
	printf("\tSIZEY=%d\n",SIZEY);
	printf("\tMAPSIZEX=%d\n",MAPSIZEX);
	printf("\tMAPSIZEY=%d\n",MAPSIZEY);
	printf("----------------------------------------\n");

	if(MAPSIZEX%2!=0){ printf("MAPSIZEX not even!!\n"); abort(); }
	if(MAPSIZEY%2!=0){ printf("MAPSIZEY not even!!\n"); abort(); }

	EPS = eps/_D;
	WH = 0.04*(1.0-EXTENSION*EPS);
	SB = WH;
	WT =0.01*(1.0-EXTENSION*EPS);
	ST = 0.95*(1.0-EXTENSION*EPS);
	DS = 1.0/((double)N-1.0);

	X = new double[N];	
	Y = new double[N];	
	S = new double[N];	
	W = new double[N];	
	NORX = new double[N];	
	NORY = new double[N];	
	VX = new double[N];	
	VY = new double[N];	
	VNORX = new double[N];	
	VNORY = new double[N];		

	memset(X,0,N*sizeof(double));
	memset(Y,0,N*sizeof(double));
	memset(S,0,N*sizeof(double));
	memset(W,0,N*sizeof(double));
	memset(NORX,0,N*sizeof(double));
	memset(NORY,0,N*sizeof(double));
	memset(VX,0,N*sizeof(double));
	memset(VY,0,N*sizeof(double));
	memset(VNORX,0,N*sizeof(double));
	memset(VNORY,0,N*sizeof(double));

	SHAPEX = new double[2*N];	
	SHAPEY = new double[2*N];	

	memset(SHAPEX,0,2*N*sizeof(double));
	memset(SHAPEY,0,2*N*sizeof(double));

	assert(SIZEX > 4);	
	const int size = SIZEX*SIZEY;
	dataX = new double[size];	
	dataY = new double[size];	
	dataVX = new double[size];	
	dataVY = new double[size];
	dataDist = new double[size];
	dataSDF = new double[size];

	memset(dataX,0,size*sizeof(double));
	memset(dataY,0,size*sizeof(double));
	memset(dataVX,0,size*sizeof(double));
	memset(dataVY,0,size*sizeof(double));
	memset(dataDist,0,size*sizeof(double));
	memset(dataSDF,0,size*sizeof(double));

	const int sizeMap = MAPSIZEX*MAPSIZEY;
	CHI = new double[sizeMap];
	CHI2= new double[sizeMap];
	VDEFX = new double[sizeMap];
	VDEFY = new double[sizeMap];
	memset(CHI,0,sizeMap*sizeof(double));
	memset(CHI2,0,sizeMap*sizeof(double));
	memset(VDEFX,0,sizeMap*sizeof(double));
	memset(VDEFY,0,sizeMap*sizeof(double));

	for(int i = 0; i < N; i++ )
	{
		S[i] = DS*(double)i;
		W[i] = _getWidth(S[i]);
	}

	fish_xmin[0] = 0.0;
	fish_xmin[1] = 0.0;
	fish_xmax[0] = 0.0;
	fish_xmax[1] = 0.0;

	_set_ow_mapping();

	_computeMinStartMaxEnd(W, minStart, maxEnd);
}

I2D_CarlingFish::Fish::~Fish()
{
	delete [] X;
	delete [] Y;
	delete [] S;
	delete [] W;
	delete [] NORX;	
	delete [] NORY;	
	delete [] VX;	
	delete [] VY;	
	delete [] VNORX;	
	delete [] VNORY;	
	delete [] SHAPEX;	
	delete [] SHAPEY;		
	delete [] dataX;
	delete [] dataY;
	delete [] dataVX;
	delete [] dataVY;
	delete [] dataDist;	
	delete [] CHI;
	delete [] VDEFX;
	delete [] VDEFY;
	delete [] dataSDF;
	delete [] CHI2;
}

int I2D_CarlingFish::Fish::_dd2l(const int & ix, const int & iy) const
{
	assert(ix>=0 && ix<SIZEX);
	assert(iy>=0 && iy<SIZEY);	
	return min(max(0,ix+SIZEX*iy),SIZEX*SIZEY-1);
}

int I2D_CarlingFish::Fish::_ud2l(const int & ix, const int & iy) const
{
	assert(ix>=0 && ix<MAPSIZEX);
	assert(iy>=0 && iy<MAPSIZEY);
	return min(max(0,ix+MAPSIZEX*iy),MAPSIZEX*MAPSIZEY-1);
}

void I2D_CarlingFish::Fish::clearUniformGrids()
{
	const int sizeMap = MAPSIZEX*MAPSIZEY;
	memset(CHI,0,sizeMap*sizeof(double));
	memset(VDEFX,0,sizeMap*sizeof(double));
	memset(VDEFY,0,sizeMap*sizeof(double));
	memset(CHI2,0,sizeMap*sizeof(double));
}

double I2D_CarlingFish::Fish::mollified_heaviside(const double dist, const double eps)
{
	//Positive outside/negative inside
	const double alpha = M_PI*min(1., max(0., (dist+0.5*eps)/eps));			
	return 0.5+0.5*cos(alpha);
}

double I2D_CarlingFish::Fish::_bsp2(double x) const
{
  const double t = fabs(x);
  if (t>1) return 0;
  return (1-t);
}

double I2D_CarlingFish::Fish::_bsp4(double x) const
{
	const double t = fabs(x);	
	if (t>2) return 0;	
	if (t>1) return pow(2-t,3)/6;	
	return (1 + 3*(1-t)*(1 + (1-t)*(1 - (1-t))))/6;
}

template<typename R>
void I2D_CarlingFish::Fish::_w2o(const R xw[2], R xo[2]) const
{
	xo[0] = W2O[0][0]*xw[0] + W2O[0][1]*xw[1] + W2O[0][2];
	xo[1] = W2O[1][0]*xw[0] + W2O[1][1]*xw[1] + W2O[1][2];
}

template<typename R>
void I2D_CarlingFish::Fish::_o2w(const R xo[2], R xw[2]) const
{
	xw[0] = O2W[0][0]*xo[0] + O2W[0][1]*xo[1] + O2W[0][2];
	xw[1] = O2W[1][0]*xo[0] + O2W[1][1]*xo[1] + O2W[1][2];
}

void I2D_CarlingFish::Fish::_check_ow_mapping() const
{
	double A[3][3];

	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
		{
			A[i][j] = 0;

			for(int d=0; d<3; d++)
				A[i][j] += O2W[i][d]* W2O[d][j];
		}

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if (i==j)
				assert(fabs(A[i][j]-1)<1e-5);
			else
				assert(fabs(A[i][j]-0)<1e-5);
		}
	}
}

void I2D_CarlingFish::Fish::_set_ow_mapping()
{
	const double tx = xm;
	const double ty = ym;
	const double ca = cos(angle);
	const double sa = sin(angle);
	const double s = D;

	W2O[0][0] = +1/s*ca;
	W2O[0][1] = +1/s*sa;
	W2O[0][2] = +1/s*(-tx*ca - ty*sa);
	W2O[1][0] = -1/s*sa;
	W2O[1][1] = +1/s*ca;
	W2O[1][2] = +1/s*(+tx*sa - ty*ca);
	W2O[2][0] = 0;
	W2O[2][1] = 0;
	W2O[2][2] = 1;

	O2W[0][0] = +s*ca;
	O2W[0][1] = -s*sa;
	O2W[0][2] = tx;
	O2W[1][0] = s*sa;
	O2W[1][1] = s*ca;
	O2W[1][2] = ty;
	O2W[2][0] = 0;
	O2W[2][1] = 0;
	O2W[2][2] = 1;

	_check_ow_mapping();
}


void I2D_CarlingFish::Fish::_rotateAboutTheOriginSingle(const double & theta, double & x, double & y) const
{
	const double a00 = cos(theta);
	const double a01 = -sin(theta);
	const double a10 = sin(theta);
	const double a11 = cos(theta);

	const double xx = x;
	const double yy = y;

	x = a00*xx + a01*yy;
	y = a10*xx + a11*yy;
}

double I2D_CarlingFish::Fish::_sample(const double x, const double y) const
{
	const int ixbase = max((int)(floor(x)-1),0);
	const int iybase = max((int)(floor(y)-1),0);

	const double xr = x - floor(x);
	const double yr = y - floor(y);

	const double wx[4] = {_bsp4(xr+1), _bsp4(xr+0), _bsp4(xr-1), _bsp4(xr-2)};
	const double wy[4] = {_bsp4(yr+1), _bsp4(yr+0), _bsp4(yr-1), _bsp4(yr-2)};

	double s = 0;

	for(int ly=0; ly<4; ly++)
		for(int lx=0; lx<4; lx++)
		{
			const int gx = min(MAPSIZEX-1, max(0, ixbase + lx));
			const int gy = min(MAPSIZEY-1, max(0, iybase + ly));

			s += wx[lx]*wy[ly]*CHI[_ud2l(gx,gy)];
		}

	return s;
}

void I2D_CarlingFish::Fish::_sample(const double x, const double y, double & Xs, double & defvelx, double & defvely ) const
{
	Xs = 0.0;
	defvelx = 0.0;
	defvely = 0.0;

	const int ixbase = max((int)(floor(x)-1),0);
	const int iybase = max((int)(floor(y)-1),0);

	const double xr = x - floor(x);
	const double yr = y - floor(y);

	const double wx[4] = {_bsp4(xr+1), _bsp4(xr+0), _bsp4(xr-1), _bsp4(xr-2)};
	const double wy[4] = {_bsp4(yr+1), _bsp4(yr+0), _bsp4(yr-1), _bsp4(yr-2)};

	for(int ly=0; ly<4; ly++)
		for(int lx=0; lx<4; lx++)
		{
			const int gx = min(MAPSIZEX-1, max(0, ixbase + lx));
			const int gy = min(MAPSIZEY-1, max(0, iybase + ly));

			Xs += wx[lx]*wy[ly]*CHI[_ud2l(gx,gy)];
			defvelx += wx[lx]*wy[ly]*VDEFX[_ud2l(gx,gy)];
			defvely += wx[lx]*wy[ly]*VDEFY[_ud2l(gx,gy)];
		}
}

double getHeavisideFDMH1(double * temp, const double h)
{
  const double dist = temp[1];
  if(dist >= +h) 
      return 1;

  if(dist <= -h) return 0;
  assert(std::abs(dist)<=h);

  // compute first primitive of H(x): I(x) = int_0^x H(y) dy                                                                                                                                                                                                              
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
  
  // gradI                                                                                                                                                                                                                                                                
  const Real gradIX = 0.5/h * (IplusX - IminuX);
  const Real gradIY = 0.5/h * (IplusY - IminuY);

  // gradU                                                                                                                                                                                                                                                                
  const Real gradUX = 0.5/h * (temp[2]-temp[0]);
  const Real gradUY = 0.5/h * (temp[5]-temp[3]);

  const Real denom = gradUX*gradUX+gradUY*gradUY;
  const Real numer = gradIX*gradUX + gradIY*gradUY;
  const Real H = denom==0 ?  numer : numer/denom;
  
  assert(H>=0 && H<=1);
  
return H;
}

double I2D_CarlingFish::Fish::_sample(const double x, const double y, const double h) const
{
  const int stencil = 6;
  double temp[stencil];
  for(int i=0; i<stencil; i++)
    temp[i]=0;

  const int idy=0;
  for(int idx=-1;idx<2;idx++)
    {
      const double my_x = x+(double)idx;
      const double my_y = y+(double)idy;

      const int ixbase = max((int)floor(my_x),0);
      const int iybase = max((int)floor(my_y),0);

      const double wx_chi[2] = {_bsp2(my_x-floor(my_x)), _bsp2(my_x-floor(my_x)-1)};
      const double wy_chi[2] = {_bsp2(my_y-floor(my_y)), _bsp2(my_y-floor(my_y)-1)};

      for(int ly=0; ly<2; ly++)
        for(int lx=0; lx<2; lx++)
          {
            const int gx = min(MAPSIZEX-1, max(0, ixbase + lx));
            const int gy = min(MAPSIZEY-1, max(0, iybase + ly));

            temp[idx+1] -= wx_chi[lx]*wy_chi[ly]*CHI2[_ud2l(gx,gy)];
	    assert(wx_chi[lx]*wy_chi[ly]>0);
          }
    }

  const int idx=0;
  for(int idy=-1;idy<2;idy++)
    {
      if (idy==0) continue;

      const double my_x = x+(double)idx;
      const double my_y = y+(double)idy;

      const int ixbase = max((int)floor(my_x),0);
      const int iybase = max((int)floor(my_y),0);

      const double wx_chi[2] = {_bsp2(my_x-floor(my_x)), _bsp2(my_x-floor(my_x)-1)};
      const double wy_chi[2] = {_bsp2(my_y-floor(my_y)), _bsp2(my_y-floor(my_y)-1)};

      for(int ly=0; ly<2; ly++)
        for(int lx=0; lx<2; lx++)
          {
            const int gx = min(MAPSIZEX-1, max(0, ixbase + lx));
            const int gy = min(MAPSIZEY-1, max(0, iybase + ly));

            temp[4+idy] -= wx_chi[lx]*wy_chi[ly]*CHI2[_ud2l(gx,gy)];
	    assert( wx_chi[lx]*wy_chi[ly]>0);
          }
    }

  return getHeavisideFDMH1(temp, h/D);
}

void I2D_CarlingFish::Fish::_sample(const double x, const double y, double & Xs, double & defvelx, double & defvely, const double h) const
{
  Xs = 0.0;
  defvelx = 0.0;
  defvely = 0.0;

  const int ixbase = max((int)(floor(x)-1),0);
  const int iybase = max((int)(floor(y)-1),0);

  const double xr = x - floor(x);
  const double yr = y - floor(y);

  const double wx_chi[2] = {_bsp2(xr), _bsp2(xr-1)};
  const double wy_chi[2] = {_bsp2(yr), _bsp2(yr-1)};

  const double wx[4] = {_bsp4(xr+1), _bsp4(xr+0), _bsp4(xr-1), _bsp4(xr-2)};                                                                                                                                                                                        
  const double wy[4] = {_bsp4(yr+1), _bsp4(yr+0), _bsp4(yr-1), _bsp4(yr-2)};

  for(int ly=0; ly<4; ly++)
    for(int lx=0; lx<4; lx++)
      {
	const int gx = min(MAPSIZEX-1, max(0, ixbase + lx));
	const int gy = min(MAPSIZEY-1, max(0, iybase + ly));

	defvelx += wx[lx]*wy[ly]*VDEFX[_ud2l(gx,gy)];
	defvely += wx[lx]*wy[ly]*VDEFY[_ud2l(gx,gy)];
      }

  const int stencil = 6;
  double temp[stencil];
  for(int i=0; i<stencil; i++)
      temp[i]=0;

  const int idy=0;
  for(int idx=-1;idx<2;idx++)
    {
      const double my_x = x+(double)idx;
      const double my_y = y+(double)idy;
      
      const int ixbase = max((int)floor(my_x),0);
      const int iybase = max((int)floor(my_y),0);

      const double wx_chi[2] = {_bsp2(my_x-floor(my_x)), _bsp2(my_x-floor(my_x)-1)};
      const double wy_chi[2] = {_bsp2(my_y-floor(my_y)), _bsp2(my_y-floor(my_y)-1)};

      for(int ly=0; ly<2; ly++)
	for(int lx=0; lx<2; lx++)
	  {
	    const int gx = min(MAPSIZEX-1, max(0, ixbase + lx));
	    const int gy = min(MAPSIZEY-1, max(0, iybase + ly));

	    temp[idx+1] -= wx_chi[lx]*wy_chi[ly]*CHI2[_ud2l(gx,gy)];
	  }
    }
  
  const int idx=0;
    for(int idy=-1;idy<2;idy++)                                                                                                                                                                                                                                         
    {
      if (idy==0) continue;

      const double my_x = x+(double)idx;
      const double my_y = y+(double)idy;

      const int ixbase = max((int)floor(my_x),0);
      const int iybase = max((int)floor(my_y),0);

      const double wx_chi[2] = {_bsp2(my_x-floor(my_x)), _bsp2(my_x-floor(my_x)-1)};
      const double wy_chi[2] = {_bsp2(my_y-floor(my_y)), _bsp2(my_y-floor(my_y)-1)};

      for(int ly=0; ly<2; ly++)
        for(int lx=0; lx<2; lx++)
          {
            const int gx = min(MAPSIZEX-1, max(0, ixbase + lx));
            const int gy = min(MAPSIZEY-1, max(0, iybase + ly));

            temp[4+idy] -= wx_chi[lx]*wy_chi[ly]*CHI2[_ud2l(gx,gy)];
          }
    }  

   Xs = getHeavisideFDMH1(temp, h/D);
}

double I2D_CarlingFish::Fish::sample(const double x_, const double y_, const double h_) const
{	
	//x_,y_ are in world coordinates
	const double xworld[2] = {x_, y_};

	double xobject[2];
	_w2o(xworld, xobject);
	xobject[0] += 0.5 + traslX;
	xobject[1] += 0.5 + traslY;
	//xobject[0] = xworld[0];
	//xobject[1] = xworld[1];

	const double myh = 1./(MAPSIZEX-1);
	if (xobject[0]< myh || xobject[0] > 1.0 ) return 0;
	if (xobject[1]< myh || xobject[1] > 1.0 ) return 0;
	return _sample(xobject[0]*(MAPSIZEX-1), xobject[1]*(MAPSIZEY-1), h_);
}

void I2D_CarlingFish::Fish::sample(const double x_, const double y_, double & Xs, double & defvelx, double & defvely, const double h_) const
{	
	Xs = 0.0;
	defvelx = 0.0;
	defvely = 0.0;

	//x_,y_ are in world coordinates
	const double xworld[2] = {x_, y_};

	double xobject[2];
	_w2o(xworld, xobject);
	xobject[0] += 0.5 + traslX;
	xobject[1] += 0.5 + traslY;
	//xobject[0] = xworld[0];
	//xobject[1] = xworld[1];

	const double myh = 1./(MAPSIZEX-1);
	if (xobject[0]< myh || xobject[0] > 1.0 ) return;
	if (xobject[1]< myh || xobject[1] > 1.0 ) return;
	_sample(xobject[0]*(MAPSIZEX-1), xobject[1]*(MAPSIZEY-1), Xs, defvelx, defvely, h_);

	//SIMPLE BACK ROTATION!
	_rotateAboutTheOriginSingle(angle, defvelx, defvely);
	defvelx *= D;
	defvely *= D;
}

double I2D_CarlingFish::Fish::sample(const double x_, const double y_) const
{
  //x_,y_ are in world coordinates                                                                                                                                                                                                                                                                                      
  const double xworld[2] = {x_, y_};

  double xobject[2];
  _w2o(xworld, xobject);
  xobject[0] += 0.5 + traslX;
  xobject[1] += 0.5 + traslY;
  //xobject[0] = xworld[0];                                                                                                                                                                                                                                                                                             
  //xobject[1] = xworld[1];                                                                                                                                                                                                                                                                                             

  const double myh = 1./(MAPSIZEX-1);
  if (xobject[0]< myh || xobject[0] > 1.0 ) return 0;
  if (xobject[1]< myh || xobject[1] > 1.0 ) return 0;
  return _sample(xobject[0]*(MAPSIZEX-1), xobject[1]*(MAPSIZEY-1));
}

void I2D_CarlingFish::Fish::sample(const double x_, const double y_, double & Xs, double & defvelx, double & defvely) const
{
  Xs = 0.0;
  defvelx = 0.0;
  defvely = 0.0;

  //x_,y_ are in world coordinates                                                                                                                                                                                                                                                                                      
  const double xworld[2] = {x_, y_};

  double xobject[2];
  _w2o(xworld, xobject);
  xobject[0] += 0.5 + traslX;
  xobject[1] += 0.5 + traslY;
  //xobject[0] = xworld[0];                                                                                                                                                                                                                                                                                             
  //xobject[1] = xworld[1];                                                                                                                                                                                                                                                                                             

  const double myh = 1./(MAPSIZEX-1);
  if (xobject[0]< myh || xobject[0] > 1.0 ) return;
  if (xobject[1]< myh || xobject[1] > 1.0 ) return;
  _sample(xobject[0]*(MAPSIZEX-1), xobject[1]*(MAPSIZEY-1), Xs, defvelx, defvely);

  //SIMPLE BACK ROTATION!                                                                                                                                                                                                                                                                                               
  _rotateAboutTheOriginSingle(angle, defvelx, defvely);
  defvelx *= D;
  defvely *= D;
}

double I2D_CarlingFish::Fish::_getWidth(const double & ss) const
{
	const double s = ss - (EXTENSION/2.0)*EPS;
	double width = 0.0;

	if( (s>=0.0) && (s<SB) )
		width = sqrt( 2.0*WH*s-s*s );

	if( (s>=SB) && (s<ST) )
		width = WH - (WH-WT)*((s-SB)/(ST-SB));
	//width = WH - (WH-WT)*((s-SB)/(ST-SB))*((s-SB)/(ST-SB));	

	if( (s>=ST) && (s<=(1.0-EXTENSION*EPS)) )
	{
		double end = 1.0-EXTENSION*EPS;
		width = WT*(end-s)/(end-ST);		
	}

	return width;
}

void I2D_CarlingFish::Fish::_cubicInterpolation(double t0, double t1, double t, double k0, double k1, double & k, double & dkdt)
{
	k = 0.0;
	dkdt = 0.0;

	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double d = 0.0;

	const double dT = t1-t0;
	const double dt = t-t0;

	// CUBIC
	d = k0;
	c = 0.0;
	b = 3.0*(k1-k0)/(dT*dT);
	a = -2.0*(k1-k0)/(dT*dT*dT);
	k = a*dt*dt*dt + b*dt*dt + c*dt + d;
	dkdt = 3.0*a*dt*dt + 2.0*b*dt + c;
}

void I2D_CarlingFish::Fish::updateInSpace(const double t)
{	
	// Head reference frame quantities
	double hCMX = 0.0;
	double hCMY = 0.0;
	double hVCMX = 0.0;
	double hVCMY = 0.0;

	// Center of mass reference frame quantities
	double cmL = 0.0;
	double cmII = 0.0;
	double omega = 0.0;

	// Frame of reference shift
	const double centerx = 0.5;
	const double centery = 0.5;
	double yyy = 0.0;
	double xxx = 0.0;

	// Calculations in the head reference frame
	double rumpUpFactor = 1.0;
	double dummy = 0.0;
	if( t <= T )
	{
		rumpUpFactor = sin((2.0*M_PI)/(4.0*T)*t);
		//_cubicInterpolation(0.0, T, t, 0.0, 1.0, rumpUpFactor, dummy);
	}
	_calculateCenterlineHead(t, rumpUpFactor, T, DS, S, X, Y, N);
	_calculateNormalsAndVelocities(t, T, DS, S, X, Y, VX, VY, NORX, NORY, VNORX, VNORY, N);
	_createShapeBoundary(X, Y, W, NORX, NORY, SHAPEX, SHAPEY, N);
	_fillDefGrid();	

	// Go into center of mass frame of reference
	_getCenterOfMassFull(hCMX, hCMY); xxx = hCMX; yyy = hCMY;
	if(xxx<0.0){ printf("cazzo xxx negative!!\n"); abort(); }

	_getVelCenterOfMassFull(hVCMX, hVCMY);
	_centerlineCenterOfMassFrameTransform(hCMX, hCMY, hVCMX, hVCMY);
	_defGridCenterOfMassFrameTransform(hCMX, hCMY, hVCMX, hVCMY);

	// Correction for the rotational impulse in the center of mass reference frame
	_getAngularMomentumFull(cmL);
	_getScalarMomentOfInertiaFull(cmII);
	omega = -cmL/cmII;
	_correctCenterlineForRotationalImpulse(omega);
	_correctDefGridForRotationalImpulse(omega);

#ifndef NDEBUG
	// Check, the angular momentum here should be almost zero
	double hL = 0.0;
	_getCenterOfMassFull(hCMX, hCMY);		
	_getVelCenterOfMassFull(hVCMX, hVCMY);	
	_getAngularMomentumFull(hL);	
	if( (fabs(hL) > 1e-10) || (fabs(hCMX) > 1e-10) || (fabs(hCMY) > 1e-10) || (fabs(hVCMX) > 1e-10) || (fabs(hVCMY) > 1e-10) )
	{
		printf("(Everything should be zero!) cmx=%e, cmy=%e\n", hCMX, hCMY);
		printf("(Everything should be zero!) vxcm=%e, vycm=%e\n", hVCMX, hVCMY);
		printf("(Everything should be zero!) L=%e\n", hL);
		abort();
	}	
#endif

	// Rotate for the internal angle
	_rotateAboutTheOrigin(angleInSpace, X, Y, N);
	_rotateDefGridAboutTheOrigin(angleInSpace);

	// OPTIONAL
	_rotateAboutTheOrigin(angleInSpace, VX, VY, N);
	_rotateAboutTheOrigin(angleInSpace, NORX, NORY, N);
	_createShapeBoundary(X, Y, W, NORX, NORY, SHAPEX, SHAPEY, N);

	//FILE * ppFile = NULL;
	//ppFile = fopen("shapecontour", "w");
	//for(int i=0; i<2*N; i++){ fprintf(ppFile,"%f %f\n",SHAPEX[i],SHAPEY[i]); }
	//fclose(ppFile);

	// Compute bounding fish box
	double height = 0.0;
	double width = 0.0;
	_bboxShapeBoundary(SHAPEX, SHAPEY, N, width, height);
	const double spaceX = (1.0-width)/2.0;
	const double spaceY = (1.0-height)/2.0;
	xxx = (xxx<=0.0)?(1.0-spaceX+xxx):xxx+spaceX;
	yyy = (yyy<=0.0)?(1.0-spaceY+yyy):yyy+spaceY;

	// Fish should fit in domain [0,1]
	_rigidTranslation(xxx,yyy);
	traslX = xxx-centerx;
	traslY = yyy-centery;
	printf("traslY=%e, yyy=%e\n", traslY, yyy);
	angular_velocityInSpace = omega;

	// Compute bounding fish box
	_bbox();

	// Reconstruct shape on uniform fine mesh
	bilinearInterpolation();
}

void I2D_CarlingFish::Fish::_calculateCenterlineHead(const double t, const double factor, const double T, const double ds, const double * ss, double * rX, double * rY, const int nn)
{
	double dxx[nn-1];  
	double dyy[nn-1];

	// Calculate y using Carling's formula
	const double offset = 0.0;//-(EXTENSION/2.0)*EPS;
	const double ll = 1.0;
	for(int i=0; i<nn; i++ ){ rY[i] = factor*0.125*ll*((ss[i]+offset)/ll+0.03125)/1.03125*sin(2.0*M_PI*((ss[i]+offset)/ll-t/T) + phase); }

	// Calculate dyy
	for(int i=0; i<nn-1; i++){ dyy[i] = rY[i+1]-rY[i]; }

	// Calculate dxx
	for(int i=0; i<nn-1; i++)
	{
#ifndef NDEBUG
		if( (ds*ds - dyy[i]*dyy[i])<= 0.0 ){ std::cout << ds << "  " << dyy[i] << " orca troia lo sapevo" << std::endl; }
#endif
		dxx[i] = sqrt( ds*ds - dyy[i]*dyy[i] );
	}

	// Calculate x
	rX[0] = 0.0;
	for(int i=1; i<nn; i++){ rX[i] = rX[i-1] + dxx[i-1]; }

#ifndef NDEBUG	
	//Check
	double length = 0.0;
	for(int i=1; i<nn; i++){ length += sqrt( (rX[i]-rX[i-1])*(rX[i]-rX[i-1]) + (rY[i]-rY[i-1])*(rY[i]-rY[i-1])  ); }
	printf("length: %e, %e", length, ll);
#endif
}

void I2D_CarlingFish::Fish::_calculateNormalsAndVelocities(const double t, const double TT, const double ds, const double * ss, const double * rX, const double * rY,
		double * vX, double * vY, double * norX, double * norY, double * vNorX, double * vNorY, const int nn)
{
	const double a00 = cos(M_PI/2.0);
	const double a01 = -sin(M_PI/2.0);
	const double a10 = sin(M_PI/2.0);
	const double a11 = cos(M_PI/2.0);

	double rXplus[nn];
	double rYplus[nn];
	double tanX[nn];
	double tanY[nn];
	double tanXPlus[nn];
	double tanYPlus[nn];
	double norXPlus[nn];
	double norYPlus[nn];

	// Calculate advanced position centerline
	const double dt = 1e-6;
	const double tPlus = t + dt;
	double rumpUpFactor = 1.0;
	double dummy = 0.0;
	if( t <= TT )
	{
		rumpUpFactor = sin((2.0*M_PI)/(4.0*TT)*tPlus);
		//_cubicInterpolation(0.0, TT, tPlus, 0.0, 1.0, rumpUpFactor, dummy);
	}
	_calculateCenterlineHead(tPlus, rumpUpFactor, TT, ds, ss, rXplus, rYplus, nn);

	// Calculate vX and Vy
	for(int i=0; i<nn; i++)
	{
		vX[i] = (rXplus[i]-rX[i])/dt;
		vY[i] = (rYplus[i]-rY[i])/dt;
	}

	// Calculate tanX and tanY
	for(int i=0; i<nn-1; i++)
	{
		const double tX = rX[i+1]-rX[i];
		const double tY = rY[i+1]-rY[i];
		const double mod = sqrt(tX*tX+tY*tY);
		tanX[i] = tX/mod;
		tanY[i] = tY/mod;
	}
	tanX[nn-1] = tanX[nn-2];
	tanY[nn-1] = tanY[nn-2];

	// Calculate norX and norY (rotation tanX and tanY of 90 degrees)
	for(int i=0; i<nn; i++)
	{
		norX[i] = a00*tanX[i] + a01*tanY[i];
		norY[i] = a10*tanX[i] + a11*tanY[i];
	}

	// Calculate tanXPlus and tanYPlus
	for(int i=0; i<nn-1; i++)
	{
		const double tXPlus = rXplus[i+1]-rXplus[i];
		const double tYPlus = rYplus[i+1]-rYplus[i];
		const double mod = sqrt(tXPlus*tXPlus+tYPlus*tYPlus);
		tanXPlus[i] = tXPlus/mod;
		tanYPlus[i] = tYPlus/mod;
	}
	tanXPlus[nn-1] = tanXPlus[nn-2];
	tanYPlus[nn-1] = tanYPlus[nn-2];

	// Calculate norXPlus and norYPlus
	for(int i=0; i<nn; i++)
	{
		norXPlus[i] = a00*tanXPlus[i] + a01*tanYPlus[i];
		norYPlus[i] = a10*tanXPlus[i] + a11*tanYPlus[i];
	}

	// Calculate vNorX and vNorY
	for(int i=0; i<nn; i++)
	{
		vNorX[i] = (norXPlus[i]-norX[i])/dt;
		vNorY[i] = (norYPlus[i]-norY[i])/dt;
	}
}

void I2D_CarlingFish::Fish::_createShapeBoundary( const double * rX, const double * rY, const double * width, const double * norX, const double * norY, double * shapeX, double * shapeY, const int n)
{
	for(int i=0; i<n; i++)
	{
		shapeX[i] = rX[i] + norX[i]*width[i];
		shapeY[i] = rY[i] + norY[i]*width[i];
	}
	for(int i=0; i<n; i++)
	{
		shapeX[2*n-1-i] = rX[i] - norX[i]*width[i];
		shapeY[2*n-1-i] = rY[i] - norY[i]*width[i];
	}
}

void I2D_CarlingFish::Fish::_getDefGridCoord(const int & ix, const int & iy, const double & dg, const int & coeff, const double * rX, const double * rY, const double * norX, const double * norY, double & xCoord, double & yCoord) const
{
	xCoord = rX[coeff*ix] + norX[coeff*ix]*iy*dg;
	yCoord = rY[coeff*ix] + norY[coeff*ix]*iy*dg;
}

void I2D_CarlingFish::Fish::_getDefGridVel(const int & ix, const int & iy, const double & dg, const int & coeff, const double * vX, const double * vY, const double * vNorX, const double * vNorY, double & vxCoord, double & vyCoord) const
{
	vxCoord = vX[coeff*ix] + vNorX[coeff*ix]*iy*dg;
	vyCoord = vY[coeff*ix] + vNorY[coeff*ix]*iy*dg;
}

void I2D_CarlingFish::Fish::_computeMinStartMaxEnd(const double * width, int & _minStart, int & _maxEnd)
{
	_minStart = 0;
	for(int i=0; i<N;i++){ if(width[i] > 0.0){ _minStart = (int)i-1; break; } }
	assert(_minStart>=0 && _minStart<N);

	_maxEnd = 0;
	for(int i=N-1; i>0;i--){ if(width[i] > 0.0){ maxEnd = (int)i+1; break; } }
	assert(_maxEnd>=0 && _maxEnd<N);
	assert(_maxEnd>=_minStart);

	_minStart = (int)max((int)_minStart,(int)0);
	_maxEnd = (int)min((int)_maxEnd,(int)N);
}

void I2D_CarlingFish::Fish::_getUniformGridProperties()
{
	const double H = 1.0/(MAPSIZEX-1);

	double lcum = 0.0;
	double IIcum = 0.0;
	double M = 0.0;
	double meanVx = 0.0;
	double meanVy = 0.0;
	double meanX = 0.0;
	double meanY = 0.0;

	// Compute actual center of mass and mean velocity
#pragma omp parallel for reduction(+:M,meanVx,meanVy,meanX,meanY) schedule(static,1)
	for(int ix=0; ix<MAPSIZEX; ix++)
		for(int iy=0; iy<MAPSIZEY; iy++)
		{
			const double xx = (double)ix*H;
			const double yy = (double)iy*H;
			const int idx = _ud2l(ix,iy);
			const double Xs = CHI[idx];
			const double vxx = VDEFX[idx];
			const double vyy = VDEFY[idx];
			M += Xs;
			meanVx += Xs*vxx;
			meanVy += Xs*vyy;
			meanX += Xs*xx;
			meanY += Xs*yy;
		}
	const double corrVx = meanVx/M;
	const double corrVy = meanVy/M;
	const double xxcm = meanX/M;
	const double yycm = meanY/M;

	// Compute actual angular momentum and scalar momentum of inertia
	lcum = IIcum = M = meanVx = meanVy = meanX = meanY = 0.0;
#pragma omp parallel for reduction(+:lcum,IIcum) schedule(static,1)
	for(int ix=0; ix<MAPSIZEX; ix++)
		for(int iy=0; iy<MAPSIZEY; iy++)
		{
			const double xx = (double)ix*H;
			const double yy = (double)iy*H;
			const int idx = _ud2l(ix,iy);
			VDEFX[idx] -= corrVx;
			VDEFY[idx] -= corrVy;

			const double Xs = CHI[idx];
			const double vxx = VDEFX[idx];
			const double vyy = VDEFY[idx];
			lcum += Xs*( (xx-xxcm)*vyy - (yy-yycm)*vxx );
			IIcum += Xs*( (xx-xxcm)*(xx-xxcm) + (yy-yycm)*(yy-yycm) );
		}
	const double omega = -lcum/IIcum;

	// Perform the necessary corrections
	lcum = IIcum = M = meanVx = meanVy = meanX = meanY = 0.0;
#pragma omp parallel for
	for(int ix=0; ix<MAPSIZEX; ix++)
		for(int iy=0; iy<MAPSIZEY; iy++)
		{
			const double xx = (double)ix*H;
			const double yy = (double)iy*H;
			const int idx = _ud2l(ix,iy);

			double vRotX = 0.0;
			double vRotY = 0.0;
			_getRotationalVelocityAboutTheOrigin(omega, (xx-xxcm), (yy-yycm), vRotX, vRotY);
			VDEFX[idx] += vRotX;
			VDEFY[idx] += vRotY;
		}

#ifndef NDEBUG
	lcum = IIcum = M = meanVx = meanVy = meanX = meanY = 0.0;
#pragma omp parallel for reduction(+:M,meanVx,meanVy,lcum,IIcum) schedule(static,1)
	for(int ix=0; ix<MAPSIZEX; ix++)
		for(int iy=0; iy<MAPSIZEY; iy++)
		{
			const double xx = (double)ix*H;
			const double yy = (double)iy*H;
			const int idx = _ud2l(ix,iy);

			const double Xs = CHI[idx];
			const double vxx = VDEFX[idx];
			const double vyy = VDEFY[idx];

			M += Xs;
			meanVx += Xs*vxx;
			meanVy += Xs*vyy;
			lcum += Xs*( (xx-xxcm)*vyy - (yy-yycm)*vxx );
			IIcum += Xs*( (xx-xxcm)*(xx-xxcm) + (yy-yycm)*(yy-yycm) );
		}
	const double meanVxAfter = meanVx/M;
	const double meanVyAfter = meanVy/M;
	const double omegaAfter = -lcum/IIcum;
	const double offsetX = (0.5+traslX);
	const double offsetY = (0.5+traslY);

	//printf("meanVxAfter=%e\n",meanVxAfter);
	//printf("meanVyAfter=%e\n",meanVyAfter);
	//printf("offsetX=%e\n",(offsetX-xxcm));
	//printf("offsetY=%e\n",(offsetY-yycm));
	//printf("omegaAfter=%e\n",omegaAfter);

	if( (fabs(meanVxAfter) > 1e-10) || (fabs(meanVyAfter) > 1e-10) || (fabs(omegaAfter) > 1e-10) )
	{
		printf("translational or rotational velocities fucked up in uniform shape representation!\n"); abort();
	}
#endif
}

void I2D_CarlingFish::Fish::_getDefGridCharFunc(const int & ix, const int & iy, const int & spany, const double & dg, const int & coeff, const double * rX, const double * rY, const double * width, double & sdf) const
{	
	const double x = dataX[_dd2l(ix,iy)];
	const double y = dataY[_dd2l(ix,iy)];
	const double signInOne = 1.0;
	const double signInRight = -(width[coeff*ix] - fabs((double)(iy-spany)*dg))/fabs(width[coeff*ix] - fabs((double)(iy-spany)*dg));
	const double signIn = (width[coeff*ix]==0.0)?signInOne:signInRight;
	double minDist = std::numeric_limits<double>::max();

	const int range = ceil((EXTENSION/2.0)*EPS/DS);
	const int start = max(coeff*ix - range, minStart);
	const int end = min(coeff*ix + range, maxEnd );

	if( (iy-spany) >= 0)
	{	
		for(int i=start; i<=end; i++)
		{
			const double d = sqrt( (x-SHAPEX[i])*(x-SHAPEX[i]) + (y-SHAPEY[i])*(y-SHAPEY[i]) );
			minDist = min(d,minDist);
		}	
	}
	else
	{
		for(int i=start; i<=end; i++)
		{
			const double d = sqrt( (x-SHAPEX[2*N-1-i])*(x-SHAPEX[2*N-1-i]) + (y-SHAPEY[2*N-1-i])*(y-SHAPEY[2*N-1-i]) );
			minDist = min(d,minDist);
		}
	}

	sdf = signIn*minDist;
}

void I2D_CarlingFish::Fish::_fillDefGrid()
{	
	// Set variables
	const int coeff = floor((double)(N-1)/(double)(SIZEX-1));
	const double dg = (double)(coeff)*DS;

	// Check consistency
	if( coeff*(SIZEX-1) != (N-1) ){ printf("non multiple!\n"); abort(); }

	const int spanY = (int)((double)(SIZEY-1)/2.0);
	if( SIZEY != 2*spanY+1 ){ printf("Error span!\n"); abort(); }

#pragma omp parallel for 
	for(int ix = 0; ix<SIZEX; ix++)
		for(int iy = -spanY; iy<=spanY; iy++)
		{
			double xCoord = 0.0;
			double yCoord = 0.0;								
			_getDefGridCoord(ix, iy, dg, coeff, X, Y, NORX, NORY, xCoord, yCoord);
			dataX[ _dd2l(ix,iy+spanY) ] = xCoord;
			dataY[ _dd2l(ix,iy+spanY) ] = yCoord;

			double vxCoord = 0.0;
			double vyCoord = 0.0;			
			_getDefGridVel(ix, iy, dg, coeff, VX, VY, VNORX, VNORY, vxCoord, vyCoord);
			dataVX[ _dd2l(ix,iy+spanY) ] = vxCoord;
			dataVY[ _dd2l(ix,iy+spanY) ] = vyCoord;

			double sdf = 0.0;
			_getDefGridCharFunc(ix, iy+spanY, spanY, dg, coeff, X, Y, W, sdf);
			dataDist[ _dd2l(ix,iy+spanY) ] = sdf;
			dataSDF [ _dd2l(ix,iy+spanY) ] = sdf;
		}

#pragma omp parallel for
	for(int ix = 0; ix<SIZEX; ix++)
		for(int iy = 0; iy<SIZEY; iy++)
		{
		  dataDist[ _dd2l(ix,iy) ] = mollified_heaviside(dataDist[ _dd2l(ix,iy) ],EPS);
		}

}

void I2D_CarlingFish::Fish::_rigidTranslation(const Real x, const Real y )
{
#pragma omp parallel for
	for(int ix = 0; ix<SIZEX; ix++)
		for(int iy = 0; iy<SIZEY; iy++)
		{
			const int idx = _dd2l(ix,iy);
			dataX[idx] += x;
			dataY[idx] += y;			
		}
}

void I2D_CarlingFish::Fish::_cross(const double * v1, const double * v2, double * v3) const
{ 
	v3[0] = v1[0]*v2[1];
	v3[1] = - v1[1]*v2[0];
}


double I2D_CarlingFish::Fish::_dot(const double * v1, const double * v2) const
{
	return v1[0]*v2[0] + v1[1]*v2[1];
}

void I2D_CarlingFish::Fish::_boundingBoxTriangle(const double * a, const double * b, const double * c, const double H, int & ixMax, int & ixMin, int & iyMax, int & iyMin) const
{
	const double xMax = max(c[0],max(a[0],b[0]));
	const double xMin = min(c[0],min(a[0],b[0]));
	const double yMax = max(c[1],max(a[1],b[1]));
	const double yMin = min(c[1],min(a[1],b[1]));

	ixMax = min(MAPSIZEX-1, max(0,(int)ceil ( xMax/H ) ) ); 
	ixMin = min(MAPSIZEX-1, max(0,(int)floor( xMin/H ) ) ); 
	iyMax = min(MAPSIZEY-1, max(0,(int)ceil ( yMax/H ) ) ); 
	iyMin = min(MAPSIZEY-1, max(0,(int)floor( yMin/H ) ) );
}

bool I2D_CarlingFish::Fish::_pointInTriangle(const double * p, const double * a, const double * b, const double * c) const
{ 
	// Compute vectors        
	const double v0[2] = {c[0]-a[0],c[1]-a[1]};
	const double v1[2] = {b[0]-a[0],b[1]-a[1]};
	const double v2[2] = {p[0]-a[0],p[1]-a[1]};

	// Compute dot products
	const double dot00 = _dot(v0, v0);
	const double dot01 = _dot(v0, v1);
	const double dot02 = _dot(v0, v2);
	const double dot11 = _dot(v1, v1);
	const double dot12 = _dot(v1, v2);

	// Compute barycentric coordinates
	const double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
	const double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	const double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	return (u>0.0) && (v>0.0) && (u+v<1.0);
}

bool I2D_CarlingFish::Fish::_bilinearTriangle(const double * p, const double * a, const double * b, const double * c, const double * va, const double * vb, const double * vc, double * interp) const
{	
	// Compute vectors        
	const double v0[2] = {c[0]-a[0],c[1]-a[1]};
	const double v1[2] = {b[0]-a[0],b[1]-a[1]};
	const double v2[2] = {p[0]-a[0],p[1]-a[1]};

	// Compute dot products
	const double dot00 = _dot(v0, v0);
	const double dot01 = _dot(v0, v1);
	const double dot02 = _dot(v0, v2);
	const double dot11 = _dot(v1, v1);
	const double dot12 = _dot(v1, v2);

	// Compute barycentric coordinates
	const double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
	const double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	const double t = (dot00 * dot12 - dot01 * dot02) * invDenom;

	interp[0] = va[0] + t*(vb[0]-va[0]) + u*(vc[0]-va[0]);
	interp[1] = va[1] + t*(vb[1]-va[1]) + u*(vc[1]-va[1]);
	interp[2] = va[2] + t*(vb[2]-va[2]) + u*(vc[2]-va[2]);

	// Check if point is in triangle
	return (u>=0.0) && (t>=0.0) && (u+t<=1.0);
}

void I2D_CarlingFish::Fish::_getCenterOfMassFull(double & xCoord, double & yCoord) const
{
	double M = 0.0;
	double xCum = 0.0;
	double yCum = 0.0;
#pragma omp parallel for reduction(+:xCum,yCum,M) schedule(static,1)
	for(int ix=0; ix<SIZEX; ix++)
		for(int iy=0; iy<SIZEY; iy++)
		{
			const int idx = _dd2l(ix,iy);
			const double Xs = dataDist[idx];
			xCum += Xs*dataX[idx];
			yCum += Xs*dataY[idx];
			M += Xs;
		}
	xCoord = xCum/M;
	yCoord = yCum/M;
}

void I2D_CarlingFish::Fish::_getVelCenterOfMassFull(double & vxCoord, double & vyCoord) const
{
	double M = 0.0;
	double vxCum = 0.0;
	double vyCum = 0.0;
#pragma omp parallel for reduction(+:vxCum,vyCum,M) schedule(static,1)
	for(int ix=0; ix<SIZEX; ix++)
		for(int iy=0; iy<SIZEY; iy++)
		{
			const int idx = _dd2l(ix,iy);
			const double Xs = dataDist[idx];
			vxCum += Xs*dataVX[idx];
			vyCum += Xs*dataVY[idx];
			M += Xs;
		}
	vxCoord = vxCum/M;
	vyCoord = vyCum/M;
}

void I2D_CarlingFish::Fish::_getAngularMomentumFull(double & L) const
{
	double lcum = 0.0;
#pragma omp parallel for reduction(+:lcum) schedule(static,1)
	for(int ix=0; ix<SIZEX; ix++)
		for(int iy=0; iy<SIZEY; iy++)
		{
			const int idx = _dd2l(ix,iy);
			const double Xs = dataDist[idx];
			lcum += Xs*( dataX[idx]*dataVY[idx] - dataY[idx]*dataVX[idx] );
		}
	L = lcum*DS*DS;
}

void I2D_CarlingFish::Fish::_getScalarMomentOfInertiaFull(double & II) const
{
	double IIcum = 0.0;
#pragma omp parallel for reduction(+:IIcum) schedule(static,1)
	for(int ix=0; ix<SIZEX; ix++)
		for(int iy=0; iy<SIZEY; iy++)
		{			
			const int idx = _dd2l(ix,iy);
			const double Xs = dataDist[idx];
			IIcum += Xs*( dataX[idx]*dataX[idx] + dataY[idx]*dataY[idx] );
		}
	II = IIcum*DS*DS;
}

void I2D_CarlingFish::Fish::_centerlineCenterOfMassFrameTransform(const double & cmX, const double & cmY, const double & vcmX, const double & vcmY)
{
	for(int i=0; i<N; i++)
	{
		X[i] -= cmX;
		Y[i] -= cmY;
		VX[i] -= vcmX;
		VY[i] -= vcmY;
	}
}

void I2D_CarlingFish::Fish::_defGridCenterOfMassFrameTransform(const double & cmX, const double & cmY, const double & vcmX, const double & vcmY)
{
#pragma omp parallel for
	for(int ix=0; ix<SIZEX; ix++)
		for(int iy=0; iy<SIZEY; iy++)
		{
			const int idx = _dd2l(ix,iy);
			dataX[idx] -= cmX;
			dataY[idx] -= cmY;
			dataVX[idx] -= vcmX;
			dataVY[idx] -= vcmY;
		}
}

void I2D_CarlingFish::Fish::_getRotationalVelocityAboutTheOrigin(const double & omega, const double & rX, const double & rY, double & vRotX, double & vRotY) const
{
	vRotX = -omega*rY;
	vRotY = omega*rX;
}

void I2D_CarlingFish::Fish::_correctCenterlineForRotationalImpulse(const double & omega)
{
	double vRotX = 0.0;
	double vRotY = 0.0;
	for(int i=0; i<N; i++)
	{
		_getRotationalVelocityAboutTheOrigin(omega, X[i], Y[i], vRotX, vRotY);
		VX[i] += vRotX;
		VY[i] += vRotY;
	}
}

void I2D_CarlingFish::Fish::_correctDefGridForRotationalImpulse(const double & omega)
{
#pragma omp parallel for
	for(int ix=0; ix<SIZEX; ix++)
		for(int iy=0; iy<SIZEY; iy++)
		{
			const int idx = _dd2l(ix,iy);
			double vRotX = 0.0;
			double vRotY = 0.0;
			_getRotationalVelocityAboutTheOrigin(omega, dataX[idx], dataY[idx], vRotX, vRotY);
			dataVX[idx] += vRotX;
			dataVY[idx] += vRotY;
		}
}

void I2D_CarlingFish::Fish::_rotateAboutTheOrigin(const double & theta, double * x, double * y, const int n)
{
	const double a00 = cos(theta);
	const double a01 = -sin(theta);
	const double a10 = sin(theta);
	const double a11 = cos(theta);
	for(int i=0; i<n; i++)
	{
		const double xx = x[i];
		const double yy = y[i];
		x[i] = a00*xx + a01*yy;
		y[i] = a10*xx + a11*yy;
	}
}

void I2D_CarlingFish::Fish::_rotateDefGridAboutTheOrigin(const double & theta)
{
	const double a00 = cos(theta);
	const double a01 = -sin(theta);
	const double a10 = sin(theta);
	const double a11 = cos(theta);
#pragma omp parallel for
	for(int ix=0; ix<SIZEX; ix++)
		for(int iy=0; iy<SIZEY; iy++)
		{
			const int idx = _dd2l(ix,iy);
			const double xx = dataX[idx];
			const double yy = dataY[idx];
			const double vx = dataVX[idx];
			const double vy = dataVY[idx];
			dataX[idx] = a00*xx + a01*yy;
			dataY[idx] = a10*xx + a11*yy;
			dataVX[idx] = a00*vx + a01*vy;
			dataVY[idx] = a10*vx + a11*vy;
		}
}

void I2D_CarlingFish::Fish::bilinearInterpolation()
{
	const double H = 1.0/(MAPSIZEX-1);
	const int sX = 0;
	const int sY = 0;
	const int eX = SIZEX-1;
	const int eY = SIZEY-1;

#pragma omp parallel for
	for(int ixD=sX; ixD<eX; ixD++)
		for(int iyD=sY; iyD<eY; iyD++)
		{		
			// First triangle
			const int ixx1 = ixD;
			const int iyy1 = iyD;
			const int linIdx1 = _dd2l(ixx1,iyy1);
			const double A[2]  = { dataX[ linIdx1 ], dataY[ linIdx1 ] };
			const double va[3] = { dataVX[ linIdx1 ], dataVY[ linIdx1 ], dataDist[ linIdx1 ] };
			const double va2[3] = { dataVX[ linIdx1 ], dataVY[ linIdx1 ], dataSDF[ linIdx1 ] };

			const int ixx2 = ixD;
			const int iyy2 = iyD+1;
			const int linIdx2 = _dd2l(ixx2,iyy2);
			const double B[2]  = { dataX[ linIdx2 ], dataY[ linIdx2 ] };
			const double vb[3] = { dataVX[ linIdx2 ], dataVY[ linIdx2 ], dataDist[ linIdx2 ] };
                        const double vb2[3] = { dataVX[ linIdx2 ], dataVY[ linIdx2 ], dataSDF[ linIdx2 ] };

			const int ixx3 = ixD+1;
			const int iyy3 = iyD;
			const int linIdx3 = _dd2l(ixx3,iyy3);
			const double C[2] =  { dataX[ linIdx3 ], dataY[ linIdx3 ] };
			const double vc[3] = { dataVX[ linIdx3 ], dataVY[ linIdx3 ], dataDist[ linIdx3 ] };
			const double vc2[3] = { dataVX[ linIdx3 ], dataVY[ linIdx3 ], dataSDF[ linIdx3 ] };

			int ixMax = 0.0;
			int ixMin = 0.0;
			int iyMax = 0.0;
			int iyMin = 0.0;
			_boundingBoxTriangle(A, B, C, H, ixMax, ixMin, iyMax, iyMin);
			for(int iy=iyMin; iy<iyMax; iy++)
				for(int ix=ixMin; ix<ixMax; ix++)
				{
					const double P[2] = { (double)ix*H, (double)iy*H };	
					double interp[3] = {0.0, 0.0, 0.0};
					const bool isIn = _bilinearTriangle(P, A, B, C, va, vb, vc, interp);
					double interp2[3] = {0.0, 0.0, 0.0};
					const bool isIn2 =  _bilinearTriangle(P, A, B, C, va2, vb2, vc2, interp2);
					if(isIn)
					{
						const int lIdx = _ud2l(ix,iy);
						VDEFX[ lIdx ] = interp[0];
						VDEFY[ lIdx ] = interp[1];
						CHI[ lIdx ] = interp[2];
					}
					const int lIdx = _ud2l(ix,iy);
						CHI2[ lIdx ] = interp2[2];					
				}

			// Second triangle
			const int ixx21 = ixD+1;
			const int iyy21 = iyD;
			const int linIdx21 = _dd2l(ixx21,iyy21);
			const double AA[2] =  { dataX[ linIdx21 ], dataY[ linIdx21 ] };
			const double vaa[3] = { dataVX[ linIdx21 ], dataVY[ linIdx21 ], dataDist[ linIdx21 ] };
                        const double vaa2[3] = { dataVX[ linIdx21 ], dataVY[ linIdx21 ], dataSDF[ linIdx21 ] };

			const int ixx22 = ixD+1;
			const int iyy22 = iyD+1;
			const int linIdx22 = _dd2l(ixx22,iyy22);
			const double BB[2] =  { dataX[ linIdx22 ], dataY[ linIdx22 ] };
			const double vbb[3] = { dataVX[ linIdx22 ], dataVY[ linIdx22 ], dataDist[ linIdx22 ] };
			const double vbb2[3] = { dataVX[ linIdx22 ], dataVY[ linIdx22 ], dataSDF[ linIdx22 ] };

			const int ixx23 = ixD;
			const int iyy23 = iyD+1;
			const int linIdx23 = _dd2l(ixx23,iyy23);
			const double CC[2] =  { dataX[ linIdx23 ], dataY[ linIdx23 ] };
			const double vcc[3] = { dataVX[ linIdx23 ], dataVY[ linIdx23 ], dataDist[ linIdx23 ] };
			const double vcc2[3] = { dataVX[ linIdx23 ], dataVY[ linIdx23 ], dataSDF[ linIdx23 ] };

			int ixMax2 = 0.0;
			int ixMin2 = 0.0;
			int iyMax2 = 0.0;
			int iyMin2 = 0.0;
			_boundingBoxTriangle(AA, BB, CC, H, ixMax2, ixMin2, iyMax2, iyMin2);
			for(int iy=iyMin2; iy<iyMax2; iy++)
				for(int ix=ixMin2; ix<ixMax2; ix++)
				{
					const double P[2] = { (double)ix*H, (double)iy*H };
					double interp[3] = {0.0, 0.0, 0.0};
					const bool isIn = _bilinearTriangle(P, AA, BB, CC, vaa, vbb, vcc, interp);					
					double interp2[3] = {0.0, 0.0, 0.0};
					const bool isIn2 = _bilinearTriangle(P, AA, BB, CC, vaa2, vbb2, vcc2, interp2);
					if(isIn)
					{
						const int lIdx = _ud2l(ix,iy);
						VDEFX[ lIdx ] = interp[0];
						VDEFY[ lIdx ] = interp[1];
						CHI[ lIdx ] = interp[2];
					}
					const int lIdx = _ud2l(ix,iy);
					CHI2[ lIdx ] = interp2[2];
				}
		}
}

void I2D_CarlingFish::Fish::update_all(double dt)
{
	xm += vx*dt;
	ym += vy*dt;
	angle += angular_velocity*dt;
	angleInSpace += angular_velocityInSpace*dt;

	_set_ow_mapping();
}

void I2D_CarlingFish::Fish::restart(FILE * f)
{
	float val;

	fscanf(f, "xm: %e\n", &val);
	xm = val;
	printf("CarlingFish::restart(): xm is %e\n", xm);

	fscanf(f, "ym: %e\n", &val);
	ym = val;
	printf("CarlingFish::restart(): ym is %e\n", ym);

	fscanf(f, "angle: %e\n", &val);
	angle = val;
	printf("CarlingFish::restart(): angle is %e\n", angle);

	fscanf(f, "angleInSpace: %e\n", &val);
	angleInSpace = val;
	printf("CarlingFish::restart(): angleInSpace is %e\n", angleInSpace);

	// REMEMBER TO DO THE MAPPING OTHERWISE YOU ARE FUCKED!
	_set_ow_mapping();
}

void I2D_CarlingFish::Fish::save(FILE * f) const
{
	fprintf(f, "xm: %20.20e\n", xm);
	fprintf(f, "ym: %20.20e\n", ym);
	fprintf(f, "angle: %20.20e\n", angle);
	fprintf(f, "angleInSpace: %20.20e\n", angleInSpace);
}

void I2D_CarlingFish::Fish::_bbox(void)
{
	double minX = std::numeric_limits<double>::max();
	double maxX = -std::numeric_limits<double>::max();
	double minY = std::numeric_limits<double>::max();
	double maxY = -std::numeric_limits<double>::max();

	for(int ix=0; ix<SIZEX; ix++)
	{
		const int iy = 0;
		const int linIdx = _dd2l(ix,iy);
		minX = min(minX, dataX[linIdx]);
		maxX = max(maxX, dataX[linIdx]);
		minY = min(minY, dataY[linIdx]);
		maxY = max(maxY, dataY[linIdx]);		
	}

	for(int ix=0; ix<SIZEX; ix++)
	{
		const int iy = SIZEY-1;
		const int linIdx = _dd2l(ix,iy);
		minX = min(minX, dataX[linIdx]);
		maxX = max(maxX, dataX[linIdx]);
		minY = min(minY, dataY[linIdx]);
		maxY = max(maxY, dataY[linIdx]);
	}

	for(int iy=0; iy<SIZEY; iy++)
	{
		const int ix = 0;
		const int linIdx = _dd2l(ix,iy);
		minX = min(minX, dataX[linIdx]);
		maxX = max(maxX, dataX[linIdx]);
		minY = min(minY, dataY[linIdx]);
		maxY = max(maxY, dataY[linIdx]);		
	}

	for(int iy=0; iy<SIZEY; iy++)
	{
		const int ix = SIZEX-1;
		const int linIdx = _dd2l(ix,iy);
		minX = min(minX, dataX[linIdx]);
		maxX = max(maxX, dataX[linIdx]);
		minY = min(minY, dataY[linIdx]);
		maxY = max(maxY, dataY[linIdx]);		
	}	

	const double offsetX = (0.5+traslX);
	const double offsetY = (0.5+traslY);

	const Real v1[2] = {minX-offsetX,minY-offsetY};
	const Real v2[2] = {minX-offsetX,maxY-offsetY};
	const Real v3[2] = {maxX-offsetX,minY-offsetY};
	const Real v4[2] = {maxX-offsetX,maxY-offsetY};

	Real xw1[2], xw2[2], xw3[2], xw4[2];

	_o2w(v1, xw1);
	_o2w(v2, xw2);
	_o2w(v3, xw3);
	_o2w(v4, xw4);

	fish_xmin[0] = min((Real)min((Real)min((Real)xw1[0],(Real)xw2[0]),(Real)xw3[0]),(Real)xw4[0]);
	fish_xmin[1] = min((Real)min((Real)min((Real)xw1[1],(Real)xw2[1]),(Real)xw3[1]),(Real)xw4[1]);
	fish_xmax[0] = max((Real)max((Real)max((Real)xw1[0],(Real)xw2[0]),(Real)xw3[0]),(Real)xw4[0]);
	fish_xmax[1] = max((Real)max((Real)max((Real)xw1[1],(Real)xw2[1]),(Real)xw3[1]),(Real)xw4[1]);

	//fish_xmin[0] = 0.0;
	//fish_xmin[1] = 0.0;
	//fish_xmax[0] = 1.0;
	//fish_xmax[1] = 1.0;

	assert(fish_xmin[0]<=fish_xmax[0]);
	assert(fish_xmin[1]<=fish_xmax[1]);
}

void I2D_CarlingFish::Fish::_bboxShapeBoundary(const double * shapeX, const double * shapeY, const int n, double & width, double & height)
{
	double minX = std::numeric_limits<double>::max();
	double maxX = -std::numeric_limits<double>::max();
	double minY = std::numeric_limits<double>::max();
	double maxY = -std::numeric_limits<double>::max();

	for(int i=0; i<2*n; i++)
	{
		minX = min(minX, shapeX[i]);
		maxX = max(maxX, shapeX[i]);
		minY = min(minY, shapeY[i]);
		maxY = max(maxY, shapeY[i]);
	}

	//FILE * ppFile = NULL;
	//ppFile = fopen("shapecontour", "w");
	//for(int i=0; i<2*n; i++){ fprintf(ppFile,"%f %f\n",shapeX[i],shapeY[i]); }
	//fclose(ppFile);
	//exit(0);

	assert(minX<=maxX);
	assert(minY<=maxY);

	width = maxX - minX;
	height = maxY - minY;

	assert(width<=1.0000001);
	assert(height<=1.0000001);
}

void I2D_CarlingFish::Fish::bbox(const double eps, double xmin[2], double xmax[2]) const
{
	xmin[0] = fish_xmin[0];
	xmin[1] = fish_xmin[1];
	xmax[0] = fish_xmax[0];
	xmax[1] = fish_xmax[1];	
}

namespace FloatingFish
{
struct FillBlocks
{
	double eps;
	I2D_CarlingFish::Fish * shape;

	FillBlocks(double eps, I2D_CarlingFish::Fish *_shape): eps(eps) { shape = _shape; }

	FillBlocks(const FillBlocks& c): eps(c.eps) { shape = c.shape; }

	static bool _is_touching(double eps, const I2D_CarlingFish::Fish * fish, const BlockInfo& info)
	{
		double min_pos[2], max_pos[2];

		info.pos(min_pos, 0,0);
		info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1);

		double bbox[2][2];
		fish->bbox(eps, bbox[0], bbox[1]);

		double intersection[2][2] = {
				max(min_pos[0], bbox[0][0]), min(max_pos[0], bbox[1][0]),
				max(min_pos[1], bbox[0][1]), min(max_pos[1], bbox[1][1]),
		};

		return
				intersection[0][1]-intersection[0][0]>0 &&
				intersection[1][1]-intersection[1][0]>0 ;
	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{
		bool bEmpty = true;

		if(_is_touching(eps, shape, info))
		{
			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
				{
					double p[2];
					info.pos(p, ix, iy);
					b(ix,iy).tmp = max( (Real)shape->sample(p[0], p[1]), b(ix,iy).tmp );
				}

			bEmpty = false;
		}
	}
};

struct FillBlocksTowers
{
  double eps;
  I2D_CarlingFish::Fish * shape;

  FillBlocksTowers(double eps, I2D_CarlingFish::Fish *_shape): eps(eps) { shape = _shape; }

  FillBlocksTowers(const FillBlocksTowers& c): eps(c.eps) { shape = c.shape; }

  inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
  {
    bool bEmpty = true;

    if(FillBlocks::_is_touching(eps, shape, info))
      {
	for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
	  for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
	    {
	      double p[2];
	      info.pos(p, ix, iy);
	      b(ix,iy).tmp = max( (Real)shape->sample(p[0], p[1], info.h[0]), b(ix,iy).tmp );
	    }

	bEmpty = false;
      }
  }
};

struct DivergenceOp: I2D_GradOfVector_4thOrder
{
	Real t;
	double eps;
	I2D_CarlingFish::Fish * shape;
	int stencil_start[3], stencil_end[3];

	DivergenceOp(double eps, I2D_CarlingFish::Fish *_shape): eps(eps), t(0)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;

		shape = _shape;
	}

	DivergenceOp(const DivergenceOp & copy): eps(eps), t(0)
	{
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -2;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;

		shape = copy.shape;
	}

	struct TmpAdd { static inline void stream(FluidElement2D& out, Real in) { out.tmp += in; } };

	template<typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const
	{
		if(FillBlocks::_is_touching(eps, shape, info))
		{
			// Make my little copy :)
			FluidElement2D * const e = &out(0,0);
			static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
			Real tmpCopy[n];
			for(int i=0; i<n; i++)
				tmpCopy[i] = e[i].tmp;

			// Compute divergence on tmp (no masking)
			_dfdx<TmpAdd, 0 >(lab, info, out);
			_dfdy<TmpAdd, 1 >(lab, info, out);

			// Perform masking
			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
				{
					const int idx = iy*FluidBlock2D::sizeX+ix;
					double p[2];
					info.pos(p,ix,iy);
					const double Xs = shape->SHARP? shape->sample(p[0], p[1], info.h[0]) : shape->sample(p[0], p[1]);
					out(ix,iy).tmp = tmpCopy[idx] + Xs*out(ix,iy).tmp;
				}
		}
	}
};

struct ComputeDivergenceDefVel
{
	double eps, xcm, ycm, vxcorr, vycorr, omegacorr;
	I2D_CarlingFish::Fish * shape;

	ComputeDivergenceDefVel(double xcm, double ycm, double vxcorr, double vycorr, double omegacorr, double eps, I2D_CarlingFish::Fish *_shape): xcm(xcm), ycm(ycm), vxcorr(vxcorr), vycorr(vycorr), omegacorr(omegacorr), eps(eps)
	{
		shape = _shape;
	}

	ComputeDivergenceDefVel(const ComputeDivergenceDefVel& c): xcm(c.xcm), ycm(c.ycm), vxcorr(c.vxcorr), vycorr(c.vycorr), omegacorr(c.omegacorr), eps(c.eps)
	{
		shape = c.shape;
	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{
		//clear all vel component
		{
			FluidElement2D * const e = &b(0,0);

			static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
			for(int i=0; i<n; i++)
			{
				e[i].u[0] = 0;
				e[i].u[1] = 0;
			}
		}

		if(FillBlocks::_is_touching(eps, shape, info))
		{
			const double xm = xcm;
			const double ym = ycm;
			const double av = -omegacorr;
			const double vx = -vxcorr;
			const double vy = -vycorr;

			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
				{
					double p[2];
					info.pos(p, ix, iy);

					double Xs = 0.0;
					double vdefx = 0.0;
					double vdefy = 0.0;

					if(shape->SHARP)
					  shape->sample(p[0], p[1], Xs, vdefx, vdefy, info.h[0]);
					else
					  shape->sample(p[0], p[1], Xs, vdefx, vdefy);

					b(ix,iy).u[0] = -av*(p[1]-ym) + vx + vdefx;
					b(ix,iy).u[1] =  av*(p[0]-xm) + vy + vdefy;
				}
		}
	}
};

struct ComputeCenterOfMass
{
	Real Uinf[2];
	double eps;
	I2D_CarlingFish::Fish * shape;
	map<int, vector<double> >& b2sum;
	map<int, bool>& b2nonempty;

	ComputeCenterOfMass(double eps, I2D_CarlingFish::Fish *_shape, map<int, vector<double> >& b2sum, Real Uinf[2], map<int, bool>& b2nonempty): eps(eps), b2sum(b2sum), b2nonempty(b2nonempty)
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];

		shape = _shape;
	}

	ComputeCenterOfMass(const ComputeCenterOfMass& c): eps(c.eps), b2sum(c.b2sum), b2nonempty(c.b2nonempty)
	{
		Uinf[0] = c.Uinf[0];
		Uinf[1] = c.Uinf[1];

		shape = c.shape;
	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{
		bool bNonEmpty = false;

		if(FillBlocks::_is_touching(eps, shape, info))
		{
			double mass = 0;
			double xbar = 0;
			double ybar = 0;

			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
				{
					double p[2];
					info.pos(p, ix, iy);

					const double Xs = shape->SHARP? shape->sample(p[0], p[1], info.h[0]) : shape->sample(p[0], p[1]);
					bNonEmpty |= Xs>0;

					mass += Xs;
					xbar += Xs*p[0];
					ybar += Xs*p[1];
				}

			assert(b2sum.find(info.blockID) != b2sum.end());
			assert(b2nonempty.find(info.blockID) != b2nonempty.end());

			b2sum[info.blockID][0] = mass*info.h[0]*info.h[0];
			b2sum[info.blockID][1] = xbar*info.h[0]*info.h[0];
			b2sum[info.blockID][2] = ybar*info.h[0]*info.h[0];

			b2nonempty[info.blockID] = bNonEmpty;
		}
	}
};

struct ComputeAllCorrections
{
	Real Uinf[2];
	double eps;
	double xcm, ycm;
	I2D_CarlingFish::Fish * shape;
	map<int, vector<double> >& b2sum;
	map<int, bool>& b2nonempty;

	ComputeAllCorrections(double xcm, double ycm, double eps, I2D_CarlingFish::Fish *_shape, map<int, vector<double> >& b2sum, Real Uinf[2], map<int, bool>& b2nonempty): xcm(xcm), ycm(ycm), eps(eps), b2sum(b2sum), b2nonempty(b2nonempty)
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];

		shape = _shape;
	}

	ComputeAllCorrections(const ComputeAllCorrections& c): xcm(c.xcm), ycm(c.ycm), eps(c.eps), b2sum(c.b2sum), b2nonempty(c.b2nonempty)
	{
		Uinf[0] = c.Uinf[0];
		Uinf[1] = c.Uinf[1];

		shape = c.shape;
	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{
		bool bNonEmpty = false;

		if(FillBlocks::_is_touching(eps, shape, info))
		{
			double mass = 0.0;
			double J = 0.0;
			double vxbar = 0.0;
			double vybar = 0.0;
			double omegabar = 0.0;

			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
				{
					double p[2];
					info.pos(p, ix, iy);

					double Xs = 0.0;
					double vdefx = 0.0;
					double vdefy = 0.0;
					if (shape->SHARP)
					  shape->sample(p[0], p[1], Xs, vdefx, vdefy, info.h[0]);
					else
					  shape->sample(p[0], p[1], Xs, vdefx, vdefy);

					bNonEmpty |= Xs>0;

					mass += Xs;
					J += Xs*((p[0]-xcm)*(p[0]-xcm) + (p[1]-ycm)*(p[1]-ycm));
					vxbar += Xs*vdefx;
					vybar += Xs*vdefy;
					omegabar += Xs*(vdefy*(p[0]-xcm)-vdefx*(p[1]-ycm));
				}

			assert(b2sum.find(info.blockID) != b2sum.end());
			assert(b2nonempty.find(info.blockID) != b2nonempty.end());

			b2sum[info.blockID][0] = mass*info.h[0]*info.h[0];
			b2sum[info.blockID][1] = J*info.h[0]*info.h[0];
			b2sum[info.blockID][2] = vxbar*info.h[0]*info.h[0];
			b2sum[info.blockID][3] = vybar*info.h[0]*info.h[0];
			b2sum[info.blockID][4] = omegabar*info.h[0]*info.h[0];

			b2nonempty[info.blockID] = bNonEmpty;
		}
	}
};

struct ComputeAll
{
	Real Uinf[2];
	double eps, xcm, ycm;
	I2D_CarlingFish::Fish * shape;
	map<int, vector<double> >& b2sum;
	map<int, bool>& b2nonempty;

	ComputeAll(double xcm, double ycm, double eps, I2D_CarlingFish::Fish *_shape, map<int, vector<double> >& b2sum, Real Uinf[2], map<int, bool>& b2nonempty): xcm(xcm), ycm(ycm), eps(eps), b2sum(b2sum), b2nonempty(b2nonempty)
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];

		shape = _shape;
	}

	ComputeAll(const ComputeAll& c): xcm(c.xcm), ycm(c.ycm), eps(c.eps), b2sum(c.b2sum), b2nonempty(c.b2nonempty)
	{
		Uinf[0] = c.Uinf[0];
		Uinf[1] = c.Uinf[1];

		shape = c.shape;
	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{
		bool bNonEmpty = false;

		if(FillBlocks::_is_touching(eps, shape, info))
		{
			double mass = 0;
			double J = 0;
			double vxbar = 0;
			double vybar = 0;
			double omegabar = 0;

			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
				{
					double p[2];
					info.pos(p, ix, iy);

					const double Xs =shape->SHARP?  shape->sample(p[0], p[1], info.h[0]) : shape->sample(p[0], p[1]);
					bNonEmpty |= Xs>0;

					mass += Xs;
					J += Xs*((p[0]-xcm)*(p[0]-xcm) + (p[1]-ycm)*(p[1]-ycm));
					vxbar += Xs*(b(ix, iy).u[0]+Uinf[0]);
					vybar += Xs*(b(ix, iy).u[1]+Uinf[1]);
					omegabar += Xs*((b(ix, iy).u[1]+Uinf[1])*(p[0]-xcm)-(b(ix, iy).u[0]+Uinf[0])*(p[1]-ycm));
				}

			assert(b2sum.find(info.blockID) != b2sum.end());
			assert(b2nonempty.find(info.blockID) != b2nonempty.end());

			b2sum[info.blockID][0] = mass*info.h[0]*info.h[0];
			b2sum[info.blockID][1] = J*info.h[0]*info.h[0];
			b2sum[info.blockID][2] = vxbar*info.h[0]*info.h[0];
			b2sum[info.blockID][3] = vybar*info.h[0]*info.h[0];
			b2sum[info.blockID][4] = omegabar*info.h[0]*info.h[0];
			b2nonempty[info.blockID] = bNonEmpty;
		}
	}
};


struct FillVelblocks
{
	double eps, xcm, ycm, vxcorr, vycorr, omegacorr;
	I2D_CarlingFish::Fish * shape;
	vector<pair< BlockInfo, VelocityBlock *> >& workitems;

	FillVelblocks(vector<pair< BlockInfo, VelocityBlock *> >& workitems, double xcm, double ycm, double vxcorr, double vycorr, double omegacorr, double eps, I2D_CarlingFish::Fish *_shape):
		workitems(workitems), xcm(xcm), ycm(ycm), vxcorr(vxcorr), vycorr(vycorr), omegacorr(omegacorr), eps(eps)
	{
		shape = _shape;
	}

	FillVelblocks(const FillVelblocks& c): workitems(c.workitems), xcm(c.xcm), ycm(c.ycm), vxcorr(c.vxcorr), vycorr(c.vycorr), omegacorr(c.omegacorr), eps(c.eps)
	{
		shape = c.shape;
	}

	inline void operator()(blocked_range<int> range) const
	{
		for(int iblock=range.begin(); iblock<range.end(); iblock++)
		{
			BlockInfo info = workitems[iblock].first;
			VelocityBlock * u_desired = workitems[iblock].second;

			const double xm = xcm;
			const double ym = ycm;
			const double av = shape->angular_velocity - omegacorr;
			const double vx = shape->vx - vxcorr;
			const double vy = shape->vy - vycorr;

			if(FillBlocks::_is_touching(eps, shape, info))
				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						double p[2];
						info.pos(p, ix, iy);

						double Xs = 0.0;
						double vdefx = 0.0;
						double vdefy = 0.0;
						if (shape->SHARP)
						  shape->sample(p[0], p[1], Xs, vdefx, vdefy, info.h[0]);
						else
						  shape->sample(p[0], p[1], Xs, vdefx, vdefy);

						if (Xs > 0)
						{
							u_desired->u[0][iy][ix] = - av*(p[1]-ym) + vx + vdefx;
							u_desired->u[1][iy][ix] = + av*(p[0]-xm) + vy + vdefy;
						}
						else
						{
							u_desired->u[0][iy][ix] = 0.0;
							u_desired->u[1][iy][ix] = 0.0;
						}
					}
		}
	}
};
}


I2D_CarlingFish::I2D_CarlingFish(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real T, const Real phase, const Real angle, const Real angleInSpace, const Real eps, const Real Uinf[2],
		I2D_PenalizationOperator& penalization, const int LMAX, const int ID, RL::RL_TabularPolicy ** policy, const int seed):
		I2D_FloatingObstacleOperator(parser, grid, D, eps, Uinf, penalization, ID, policy, seed), xm_corr(0.0), ym_corr(0.0), vx_corr(0.0), vy_corr(0.0), omega_corr(0.0)
{
	this->Uinf[0] = Uinf[0];
	this->Uinf[1] = Uinf[1];

	const bool isSharp = parser("-sharp").asBool();

	shape = new Fish(_xm, _ym, D, T, phase, angle, angleInSpace, eps, LMAX, isSharp);
	
	shape->updateInSpace(0);
}


I2D_CarlingFish::~I2D_CarlingFish()
{
	assert(shape!=NULL);
	if(shape!=NULL)
	{
		delete shape;
		shape = NULL;
	}

	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

void I2D_CarlingFish::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	if(!(shape->SHARP))
	    {
	      FloatingFish::FillBlocks fill(eps,shape);
	      block_processing.process(vInfo, coll, fill);
	    }
	  else
	    {
	    FloatingFish::FillBlocksTowers fill(eps,shape); 
	    block_processing.process(vInfo, coll, fill);
	    }
}

void I2D_CarlingFish::restart(const double t, string filename)
{
	// Restart shape
	FILE * ppFile = NULL;
	assert(filename!=std::string());
	ppFile = fopen(filename.c_str(), "r");
	assert(ppFile!=NULL);
	shape->restart(ppFile);

	if(policy!=NULL)
		if((*policy)!=NULL)
			(*policy)->restart();

	// Cloase file
	fclose(ppFile);
}

void I2D_CarlingFish::save(const double t, string filename)
{
	FILE * ppFile = NULL;
	assert(filename!=std::string());

	// If string is set I open the corresponding file
	ppFile = fopen(filename.c_str(), "w");
	assert(ppFile!=NULL);

	// Actual restart
	shape->save(ppFile);

	// Cloase file
	fclose(ppFile);
}

void I2D_CarlingFish::refresh(const double t, string filename)
{
	// Open file stream
	ifstream filestream;
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		filestream.open("restart_I2D_CarlingFish.txt");
		if (!filestream.good())
		{
			cout << "ooops: file not found. Exiting now." << endl;
			exit(-1);
		}
	}
	else
	{
		filestream.open(filename.c_str());
		if (!filestream.good())
		{
			cout << "ooops: file not found. Exiting now." << endl;
			exit(-1);
		}
	}

	// Open file
	int c = 0;
	string line;
	while( getline(filestream, line) ) c++;
	filestream.close();

	FILE * f = NULL;
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		f = fopen("restart_I2D_CarlingFish.txt", "r");
		assert(f!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		f = fopen(filename.c_str(), "r");
		assert(f!=NULL);
	}

	//data stored in dataVect until t=t_restart
	int N=0;
	vector < vector<float> > dataVect;
	for(int i=0; i<c; i++)
	{
		vector<float> row;
		float variable = 0.0;
		fscanf(f,"%f",&variable);
		if (variable <= t)
		{
			const Real time = variable;
			row.push_back(time);
			fscanf(f,"%f",&variable);
			const Real xm = variable;
			row.push_back(xm);
			fscanf(f,"%f",&variable);
			const Real ym = variable;
			row.push_back(ym);
			fscanf(f,"%f",&variable);
			const Real vx = variable;
			row.push_back(vx);
			fscanf(f,"%f",&variable);
			const Real vy = variable;
			row.push_back(vy);
			fscanf(f,"%f",&variable);
			const Real angle = variable;
			row.push_back(angle);
			fscanf(f,"%f",&variable);
			const Real angular_velocity = variable;
			row.push_back(angular_velocity);
			fscanf(f,"%f",&variable);
			const Real angle_in_space = variable;
			row.push_back(angle_in_space);
			fscanf(f,"%f",&variable);
			const Real m = variable;
			row.push_back(m);
			N=i;
		}
		else
		{
			break;
		}
		dataVect.push_back(row);
	}
	fclose(f);
	f = NULL;

	//dataVect is copied in a new "shape_001" file
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		f = fopen("restart_I2D_CarlingFish.txt", "w");
		assert(f!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		f = fopen(filename.c_str(), "w");
		assert(f!=NULL);
	}
	for (int i=0; i<N; i++)
	{
		fprintf(f, "%e %e %e %e %e %e %e %e %e\n", dataVect[i][0],dataVect[i][1],dataVect[i][2],dataVect[i][3],dataVect[i][4],dataVect[i][5],dataVect[i][6],dataVect[i][7],dataVect[i][8]);
	}
	fclose(f);
}

void I2D_CarlingFish::create(const double t)
{
	// Here fish computations!!
	shape->clearUniformGrids();
	shape->updateInSpace(t);

	const int NQUANTITIES = 5;

	double mass = 0.0;
	double J = 0.0;
	double xbar = 0.0;
	double ybar = 0.0;
	double vxbar = 0.0;
	double vybar = 0.0;
	double omegabar = 0.0;

	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	map<int, vector<double> > integrals;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		integrals[it->blockID] = vector<double>(NQUANTITIES);

	map<int, bool> nonempty;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		nonempty[it->blockID] = false;

	// Compute actual center of mass given interpolated Xs
	mass = J = xbar = ybar = vxbar = vybar = omegabar = 0.0;
	FloatingFish::ComputeCenterOfMass computeCenterOfMass(eps, shape, integrals, Uinf, nonempty);
	block_processing.process(vInfo, coll, computeCenterOfMass);
	for(map<int, vector< double> >::const_iterator it= integrals.begin(); it!=integrals.end(); ++it)
	{
		mass += (it->second)[0];
		xbar += (it->second)[1];
		ybar += (it->second)[2];
	}
	xbar /= mass;
	ybar /= mass;
	xm_corr = xbar;
	ym_corr = ybar;

	// Compute all corrections terms for deformation velocity field
	mass = J = xbar = ybar = vxbar = vybar = omegabar = 0.0;
	FloatingFish::ComputeAllCorrections computeAllCorrections(xm_corr, ym_corr, eps, shape, integrals, Uinf, nonempty);
	block_processing.process(vInfo, coll, computeAllCorrections);
	for(map<int, vector< double> >::const_iterator it= integrals.begin(); it!=integrals.end(); ++it)
	{
		mass += (it->second)[0];
		J += (it->second)[1];
		vxbar += (it->second)[2];
		vybar += (it->second)[3];
		omegabar += (it->second)[4];
	}
	vxbar /= mass;
	vybar /= mass;
	omegabar /= J;
	vx_corr = vxbar;
	vy_corr = vybar;
	omega_corr = omegabar;
	//printf("vxCorrection_before=%e\n",vx_corr);
	//printf("vyCorrection_before=%e\n",vy_corr);
	//printf("omegaCorrection_before=%e\n",omega_corr);

	// Interpolate corrected deformation velocities on computational grid (no masking with Xs)
	FloatingFish::ComputeDivergenceDefVel computeDivergenceDefVel(xm_corr, ym_corr, vx_corr, vy_corr, omega_corr, eps, shape);
	block_processing.process(vInfo, coll, computeDivergenceDefVel);

	// Compute the divergence of the raw def vel field and mask it with Xs
	BoundaryInfo& binfo=grid.getBoundaryInfo();
	FloatingFish::DivergenceOp div(eps, shape);
	block_processing.process< I2D_VectorBlockLab< Streamer_Velocity, 2 >::Lab >(vInfo, coll, binfo, div);
}

void I2D_CarlingFish::computeDesiredVelocity(const double t)
{
	const int NQUANTITIES = 5;

	double mass = 0.0;
	double J = 0.0;
	double xbar = 0.0;
	double ybar = 0.0;
	double vxbar = 0.0;
	double vybar = 0.0;
	double omegabar = 0.0;

	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	map<int, vector<double> > integrals;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		integrals[it->blockID] = vector<double>(NQUANTITIES);

	map<int, bool> nonempty;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		nonempty[it->blockID] = false;

	// Compute all from the flow
	mass = J = xbar = ybar = vxbar = vybar = omegabar = 0.0;
	FloatingFish::ComputeAll computeAll(xm_corr, ym_corr, eps, shape, integrals, Uinf, nonempty);
	block_processing.process(vInfo, coll, computeAll);
	for(map<int, vector< double> >::const_iterator it= integrals.begin(); it!=integrals.end(); ++it)
	{
		mass += (it->second)[0];
		J += (it->second)[1];
		vxbar += (it->second)[2];
		vybar += (it->second)[3];
		omegabar += (it->second)[4];			
	}
	vxbar /= mass;
	vybar /= mass;
	omegabar /= J;

	// Set the right vx, vy and angular velocity
	shape->vx = vxbar;
	shape->vy = vybar;
	shape->angular_velocity = omegabar;
	shape->m = mass;
	shape->J = J;

	// Prepare desired velocity blocks
	for	(map<int, const VelocityBlock *>::iterator it = desired_velocity.begin(); it!= desired_velocity.end(); it++)
	{
		assert(it->second != NULL);
		VelocityBlock::deallocate(it->second);
	}
	desired_velocity.clear();
	vector<pair< BlockInfo, VelocityBlock *> > velblocks;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
	{
		if(nonempty[it->blockID] == true)
		{
			VelocityBlock * velblock = VelocityBlock::allocate(1);
			desired_velocity[it->blockID] = velblock;
			velblocks.push_back(pair< BlockInfo, VelocityBlock *>(*it, velblock));
		}
	}

	// Set desired velocities
	FloatingFish::FillVelblocks fillvelblocks(velblocks, xm_corr, ym_corr, vx_corr, vy_corr, omega_corr, eps, shape);
	tbb::parallel_for(blocked_range<int>(0, velblocks.size()), fillvelblocks, auto_partitioner());
}

void I2D_CarlingFish::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	// If string is set I open the corresponding file
	assert(filename!=std::string());
	FILE * ppFile = NULL;
	ppFile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
	assert(ppFile!=NULL);
	fprintf(ppFile, "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", t, shape->xm, shape->ym, shape->vx, shape->vy, shape->angle, shape->angular_velocity, shape->angleInSpace, shape->m);
	fflush(ppFile);
	fclose(ppFile);

	// Cumulative diagnostics
	string cumulative = filename + "_cum";
	ppFile = fopen(cumulative.c_str(),"a");
	assert(ppFile!=NULL);
	fprintf(ppFile, "%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", t, shape->xm, shape->ym, shape->vx, shape->vy, shape->angle, shape->angular_velocity, shape->angleInSpace, shape->m);
	fflush(ppFile);
	fclose(ppFile);

	shape->update_all(dt);
}

void I2D_CarlingFish::getInfo(vector< vector<double> > & output) const
{
	output.clear();

	// Gather my info
	const double myX[2] = {shape->xm,shape->ym};
	const double myV[2] = {shape->vx,shape->vy};
	const double mySize = shape->D;

	vector<double> tmp;
	tmp.push_back(myX[0]);
	tmp.push_back(myX[1]);
	tmp.push_back(myV[0]);
	tmp.push_back(myV[1]);
	tmp.push_back(mySize);

	output.push_back(tmp);
}

