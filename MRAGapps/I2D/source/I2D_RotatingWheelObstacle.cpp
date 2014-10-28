/*
 *  I2D_RotatingWheelObstacle.cpp
 *  I2D_ROCKS
 *
 *  Created by Diego Rossinelli on 1/11/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>

using namespace std;

#include "I2D_RotatingWheelObstacle.h"

namespace RotatingWheels
{
	class Wheel
	{
		Real O2W[3][3], W2O[3][3];
		
	public:
		double angle, xm, ym, width, angular_velocity;
		double J, xg[2], area;
	private:
		
		Real * data;
		
		const int sizeX;
		const int N;
		
		void _compute_rotational_info()
		{
			const double mydx = width/(sizeX-1);
			
			{				
				double s = 0;
				
				for(int i=0; i<N; i++) 
					s += data[i];
				
				area = s*mydx*mydx;
			}
			
			{
				double sx = 0, sy = 0, m=0;
				
				//baricentro
				for(int iy=0; iy<sizeX; iy++)
				{
					const Real yw = ym-width/2 + iy*mydx;
					
					for(int ix=0; ix<sizeX; ix++)
					{
						const double xw = xm-width/2 + ix*mydx;
						
						const Real Xs = data[ix + sizeX*iy];
						
						sx += Xs*xw;
						sy += Xs*yw;
						m += Xs;
					}
				}
				
				xg[0] = sx/m;
				xg[1] = sy/m;
			}
			
			{
				double s = 0;
				double m = 0;
				
				//momento di inerzia
				for(int iy=0; iy<sizeX; iy++)
				{
					const Real yw = ym-width/2 + iy*mydx;
					
					for(int ix=0; ix<sizeX; ix++)
					{
						const double xw = xm-width/2 + ix*mydx;
						
						const Real Xs = data[ix + sizeX*iy];
						assert(Xs>=0 && Xs<=1);
						
						s += Xs*(pow(xw-xm, 2) + pow(yw-ym, 2));
						m += Xs;
					}
				}
				
				J = s*mydx*mydx;
			}
			
		}
		
		void _lowpass(const int iterations)
		{
			Real * tmp = new Real[N];
			
			for(int iter=0; iter<iterations; iter++)
			{
				for(int i=0; i<N; i++)
					tmp[i] = data[i];
				
				//x-sweep
				for(int iy=0; iy<sizeX; iy++)
					for(int ix=1; ix<sizeX-1; ix++)
					{
						assert((ix > 0 && ix < sizeX-1));
						data[ix + sizeX*iy] = (tmp[ix-1+ sizeX*iy]+ tmp[ix+1+ sizeX*iy])*0.25f + tmp[ix+ sizeX*iy]*0.5f;
					}
				
				for(int iy=0; iy<sizeX; iy++)
					for(int ix=0; ix<sizeX; ix++)
					{
						if (ix > 0 && ix < sizeX-1) 
						{
							ix = sizeX-2;
							continue;
						}
						
						assert(!(ix > 0 && ix < sizeX-1));
						
						const Real myleft = (ix==0)? tmp[ix + sizeX*iy] : tmp[ix-1 + sizeX*iy];
						const Real myright = (ix==sizeX-1)? tmp[ix+ sizeX*iy] : tmp[ix+1 + sizeX*iy];
						
						data[ix + sizeX*iy] = (myleft + myright)*0.25f + tmp[ix+ sizeX*iy]*0.5f;
					}
				
				for(int i=0; i<N; i++)
					tmp[i] = data[i];
				
				//y-sweep
				for(int iy=1; iy<sizeX-1; iy++)
					for(int ix=0; ix<sizeX; ix++)
					{
						assert((iy > 0 && iy < sizeX-1));
						data[ix + sizeX*iy] = (tmp[ix + sizeX*(iy-1)]+ tmp[ix+ sizeX*(iy+1)])*0.25f + tmp[ix+ sizeX*iy]*0.5f;
					}
				
				for(int iy=0; iy<sizeX; iy++)
					for(int ix=0; ix<sizeX; ix++)
					{
						if (iy > 0 && iy < sizeX-1) 
						{
							ix = sizeX-1;
							iy = sizeX-2;
							continue;
						}
						
						assert(!(iy > 0 && iy < sizeX-1));
						
						const Real mydown = (iy==0)? tmp[ix + sizeX*iy] : tmp[ix + sizeX*(iy-1)];
						const Real myup = (iy==sizeX-1)? tmp[ix+ sizeX*iy] : tmp[ix + sizeX*(iy+1)];
						
						data[ix + sizeX*iy] = (mydown + myup)*0.25f + tmp[ix+ sizeX*iy]*0.5f;
					}
			}
			
			delete [] tmp;
		}
		
		Real _bsp4(Real x)  const
		{
			const Real t = fabs(x);
			
			if (t>2) return 0;
			
			if (t>1) return pow(2-t,3)/6;
			
			return (1 + 3*(1-t)*(1 + (1-t)*(1 - (1-t))))/6;
		}
		
		template<typename R>
		void _w2o(const R xw[2], R xo[2]) const
		{
			xo[0] = W2O[0][0]*xw[0] + W2O[0][1]*xw[1] + W2O[0][2];
			xo[1] = W2O[1][0]*xw[0] + W2O[1][1]*xw[1] + W2O[1][2];
		}
		
		template<typename R>
		void _o2w(const R xo[2], R xw[2]) const
		{
			xw[0] = O2W[0][0]*xo[0] + O2W[0][1]*xo[1] + O2W[0][2];
			xw[1] = O2W[1][0]*xo[0] + O2W[1][1]*xo[1] + O2W[1][2];
		}
		
		void _check_ow_mapping() const
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
		
		void _set_ow_mapping()
		{
			const double tx = xm;
			const double ty = ym;
			const double ca = cos(angle);
			const double sa = sin(angle);
			const double s = width;
			
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
		
		Real _sample(const Real x, const Real y) const
		{
			const int ixbase = max((int)(floor(x)-1),0);
			const int iybase = max((int)(floor(y)-1),0);
			
			const Real xr = x - floor(x);
			const Real yr = y - floor(y);
			
			const Real wx[4] = {_bsp4(xr+1), _bsp4(xr+0), _bsp4(xr-1), _bsp4(xr-2)};
			const Real wy[4] = {_bsp4(yr+1), _bsp4(yr+0), _bsp4(yr-1), _bsp4(yr-2)};
			
			Real s = 0;
			
			for(int ly=0; ly<4; ly++)
				for(int lx=0; lx<4; lx++)
				{
					const int gx = min(sizeX-1, max(0, ixbase + lx));
					const int gy = min(sizeX-1, max(0, iybase + ly));
					
					s += wx[lx]*wy[ly]*data[gx + sizeX*gy];
				}
			
			return s;
		}
		
		bool _check_compact_support() const
		{
			for(int iy=0; iy<sizeX; iy++)
				for(int ix=0; ix<sizeX; ix++)
				{
					if (ix > 0 && ix < sizeX-1 && iy > 0 && iy < sizeX-1) 
					{
						ix = sizeX-2;
						
						continue;
					}
					
					assert(data[ix + sizeX*iy]>=0);
					
					if (data[ix + sizeX*iy] > 0) return false;
				}
			
			return true;
		}
		
	public:
		
		const Real * get_dataptr() { return data; } 
		
		Wheel(Real xm, Real ym, Real width, Real angle_rad, const int sizeX=512): 
		sizeX(sizeX), data(NULL), N(sizeX*sizeX), angle(angle_rad), width(width), 
		xm(xm), ym(ym), J(0), angular_velocity(0), area(0)
		{
			assert(sizeX > 4);
			
			data = new Real[sizeX*sizeX];
			
			const Real xcenter = sizeX*0.5;
			
			string sWHEELTYPE = "cookie";
			
			{
				//FILE * f = fopen("wheels.type", "r");
				ifstream input("wheels.type");
				
				if (input.good())
				{
					input >> sWHEELTYPE;
					
					cout << "Wheel::Wheel(): type is: "<< sWHEELTYPE << endl;
				}
			}
			
			if (sWHEELTYPE == "cookie")
			{
				for(int iy=0; iy<sizeX; iy++)
					for(int ix=0; ix<sizeX; ix++)
					{
						const Real x = (ix-xcenter);
						const Real y = (iy-xcenter);
						
						const Real r = sqrt(x*x + y*y);
						
						const Real alpha = atan2(y, x);
						const Real lambda = pow(sin(1.5*alpha), 4);
						
						const Real r_min = sizeX*0.5*0.2;//0.2;
						const Real r_max = sizeX*0.5*0.75;
						
						data[ix + sizeX*iy] = (r<=r_min + lambda*(r_max-r_min));
					}
			}
			else if (sWHEELTYPE == "dentellato")
			{
				for(int iy=0; iy<sizeX; iy++)
					for(int ix=0; ix<sizeX; ix++)
					{
						const double x = (ix-xcenter);
						const double y = (iy-xcenter);
						
						const double r = sqrt(x*x + y*y);
						const double R = sizeX*0.5*0.85;
						
						const double alpha_ = atan2(y, x) + 2*M_PI;
						
						const double alphafan = 2*M_PI/3;
						const double alpha = alpha_ - floor(alpha_/alphafan)*alphafan;
						
						const double x0 = R*0.2*cos(alphafan);
						const double y0 = R*0.2*sin(alphafan);
						
						const double x1 = R;
						const double y1 = 0;
						
						const double invIvI =  1./sqrt(pow(y1-y0,2) + pow(x1-x0,2));
						const double vx = -(y1-y0)*invIvI;
						const double vy = +(x1-x0)*invIvI;
						
						const double t = -((double) (vx * x1) - (double) (vy * y0) - (double) (vx * x0) + 
										   sqrt((double) (-2 * vx * x1 * vy * y0 + 
														  2 * vy * y0 * vx * x0 - vx * vx * y0 * y0 + 
														  4 * vx * vx * R * R + 
														  2 * vy * vy * x0 * x1 - vy * vy * x1 * x1 - vy * vy * x0 * x0 + 
														  4 * vy * vy * R * R))) / (double) (vx * vx + vy * vy) / 0.2e1;
						
						const double xc = (x0+x1)*0.5 + vx*t;
						const double yc = (y0+y1)*0.5 + vy*t;
						
						data[ix + sizeX*iy] = pow(xc - r*cos(alpha),2) + pow(yc - r*sin(alpha), 2) < R*R;// < 0);
					}
			}
			else if (sWHEELTYPE == "attempt2")
			{
				for(int iy=0; iy<sizeX; iy++)
					for(int ix=0; ix<sizeX; ix++)
					{
						const float x_ = (ix-xcenter);
						const float y_ = (iy-xcenter);
						
						const float r = sqrt(x_*x_ + y_*y_);
						const float R = sizeX*0.5*1.0;
						
						const float alpha_ = atan2(y_, x_) + 2*M_PI;
						
						const double alphafan = 2*M_PI/3;
						const float alpha = alpha_ - floor(alpha_/alphafan)*alphafan;
						
						const double k = 0.22;
						const double x0 = R*k*cos(alphafan);
						const double y0 = R*k*sin(alphafan);
						
						const double x1 = 1.3*k*R;
						const double y1 = 0;
						
						const double x2 = cos(atan2((double)(y0),(double)(R)))*R;
						
						const double x = r*cos(alpha);
						const double y = r*sin(alpha);
						
						const double f = y0;//pow(x-x0,2)*(y0/(5*(x1-x0)*(x2-x0))) + (x-x0)*(-y0/(5*(x1-x0))) + y0;
						
						
						bool bExtra;
						
						{
							const double xg = x1*0.65+x2*0.35;
							const double yg = 1/(R*2)*pow(xg-x1,2) + (y0/(x2-x1) - (x2-x1)/(R*2))*(xg-x1);
							const double fprime = 1/R*(xg-x1)+ (y0/(x2-x1) - (x2-x1)/(R*2));
							const double angle = M_PI/2 - atan(fprime);
							assert(angle>=0);
							assert(angle<=M_PI/2);
							assert(sin(alpha)>0);
							const double rfit = (y0-yg)*0.69/(1+sin(alpha));
							const double xfit = xg + cos(atan(fprime)+M_PI/2)*rfit;
							const double yfit = yg + sin(atan(fprime)+M_PI/2)*rfit;
							assert(xfit<xg);
							assert(yfit<y0);
							assert(yfit>yg);
							assert(xg>x1);
							assert(xg<x2);
							//assert(fabs(y0-yfit - rfit)<1e-6);
							const double check_dist = sqrt(pow(xg-xfit,2)+pow(yg-yfit,2));
							assert(fabs(check_dist - rfit)<1e-8);
							
							const double rinfo = sqrt(pow(x-xfit,2) + pow(y-yfit,2));
							const double ainfo = atan2(y-yfit, x-xfit);
							
							const bool bCool = ainfo<-angle || ainfo>M_PI/2;
							const bool bCheck = rinfo<rfit;
							
							bExtra = bCool || bCheck;
						}
						
						const bool b1 = (x>=0 && x<=x0 && y <= x*y0/x0);
						const bool b2 = (x>=x0 && x<=x1 && y>=y1 && y<=f);
						const bool b3 = (x>=x1 && x<=x2 && y<=f && y>= 1/(R*2)*pow(x-x1,2) + (y0/(x2-x1) - (x2-x1)/(R*2))*(x-x1)); 
						data[ix + sizeX*iy] = (b1 || b2 || b3) && bExtra;//pow(xc - r*cos(alpha),2) + pow(yc - r*sin(alpha), 2) < R*R;// < 0);
					}
			}
			else 
			{
				printf("Wheel::Wheel: incorrect wheel type: %s  -> aborting\n", sWHEELTYPE.c_str());
				abort();
			}
			
			_lowpass(20);
			assert(_check_compact_support()); //exit(0);
			
			_set_ow_mapping();
			_compute_rotational_info();
		}
		
		const Wheel& operator=(const Wheel& c)
		{
			delete [] data;
			
			memcpy(this, &c, sizeof(Wheel));
			
			data = new Real[N];
			
			memcpy(data, c.data, N*sizeof(Real));
			
			return *this;
		}
		
		void update_angle(double torque, double dt)
		{
			angular_velocity = (torque/J);
			angle += angular_velocity*dt;
			
			_set_ow_mapping();
			
			{
				FILE * f = fopen("rot.data", "a");
				fprintf(f, "%e %e %e %e %e %e\n", angle, angular_velocity, torque, (torque/J)*(dt/area), torque*dt, 1./J/area);
				fclose(f);
			}
		}
		
		void restart(FILE * f)
		{
			//	FILE *f = fopen("rotatingwheel.restart", "r");
			float val;
			fscanf(f, "angle: %e\n", &val);
			//	fclose(f);
			
			angle = val;
			_set_ow_mapping();
			
			printf("Wheel::restart(): angle is %e\n", angle);
		}
		
		void save(FILE * f) const
		{
			//FILE * f = fopen("rotatingwheel.restart", "w");
			fprintf(f, "angle: %20.20e\n", angle);
			//fclose(f);
		}
		
		Real sample(const Real x_, const Real y_) const
		{
			//x_,y_ are in world coordinates
			Real xworld[2] = {x_, y_};
			
			Real xobject[2];
			_w2o(xworld, xobject);
			
			const Real myh = 1./(sizeX-1);
			
			if (xobject[0]< -0.5-myh || xobject[0] > 0.5+myh ) return 0;
			if (xobject[1]< -0.5-myh || xobject[1] > 0.5+myh ) return 0;
			
			return _sample((xobject[0]+0.5)*(sizeX-1), (xobject[1]+0.5)*(sizeX-1));
		}
		
		void bbox(const Real eps, Real xmin[2], Real xmax[2]) const
		{
			assert(eps>=0);
			
			for(int i=0; i<4; i++)
			{
				const int dx = i%2;
				const int dy = i/2;
				
				Real xo[2] = {
					-0.5 + dx,
					-0.5 + dy
				};
				
				Real xw[2];
				_o2w(xo, xw);
				
				if (i==0)
				{
					xmin[0] = xmax[0] = xw[0];
					xmin[1] = xmax[1] = xw[1];
				}
				else
				{
					xmin[0] = min(xw[0], xmin[0]);
					xmin[1] = min(xw[1], xmin[1]);
					xmax[0] = max(xw[0], xmax[0]);
					xmax[1] = max(xw[1], xmax[1]);
				}
			}
			
			xmin[0] -= 2*eps;
			xmin[1] -= 2*eps;
			xmax[0] += 2*eps;
			xmax[1] += 2*eps;
			
			assert(xmin[0]<=xmax[0]);
			assert(xmin[1]<=xmax[1]);
		}
		
	public:
		
		static vector<Wheel *> wheels;
		
		static void restart()
		{
			FILE *f = fopen("rotatingwheel.restart", "r");
			for(unsigned int i=0; i<wheels.size(); i++)
				wheels[i]->restart(f);
			fclose(f);
		}
		
		static void save()
		{
			FILE * f = fopen("rotatingwheel.restart", "w");
			for(unsigned int i=0; i<wheels.size(); i++)
				wheels[i]->save(f);
			fclose(f);
		}
		
		static void update_angle(vector<double> omegabar, double dt, double t)
		{
			assert(omegabar.size() == wheels.size());
			
			for(unsigned int i=0; i<wheels.size(); i++)
				wheels[i]->update_angle(omegabar[i], dt);
			
			//-------------------------------------- the following part is a request from petros
			ofstream output("wheels.output", t==0? ios::out : ios::app);
			
			output.precision(20);
			output << scientific << t << "\t"; 
			for(unsigned int i=0; i<wheels.size(); i++)
				output << scientific << wheels[i]->angle << "\t";
			
			output << endl;
			//-------------------------------------- the preceding part is a request from petros
		}
	};	
	
	vector<Wheel *> Wheel::wheels;
	
	struct FillBlocks
	{
		Real eps;
		
		FillBlocks(Real eps): eps(eps) {}
		
		FillBlocks(const FillBlocks& c): eps(c.eps){}
		
		static bool _is_touching(Real eps, const Wheel& wheel, const BlockInfo& info) 
		{
			Real min_pos[2], max_pos[2];
			
			info.pos(min_pos, 0,0);
			info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1);
			
			Real bbox[2][2];
			wheel.bbox(eps, bbox[0], bbox[1]);
			
			Real intersection[2][2] = {
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
			
			const int NWHEELS = Wheel::wheels.size();
			
			for(int iwheel = 0; iwheel<NWHEELS; iwheel++)
			{
				const Wheel& mywheel = *Wheel::wheels[iwheel];
				
				if(_is_touching(eps, mywheel, info))
				{
					for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
						for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
						{
							Real p[2];
							info.pos(p, ix, iy);
							
							b(ix, iy).tmp = mywheel.sample(p[0], p[1]);
						}
					
					bEmpty = false;
				}					
			}
			
			if (bEmpty) 
			{
				FluidElement2D * const e = &b(0,0);
				static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
				for(int i=0; i<n; i++) 
					e[i].tmp = 0;
			}
		}
	};
	
	struct ComputeOmegaBar
	{
		Real Uinf[2];
		Real eps;
		map<int, vector<double> >& b2sum;
		map<int, bool>& b2nonempty;
		
		ComputeOmegaBar(Real eps, map<int, vector<double> >& b2sum, Real Uinf[2], map<int, bool>& b2nonempty): eps(eps), b2sum(b2sum), b2nonempty(b2nonempty)
		{
			this->Uinf[0] = Uinf[0];
			this->Uinf[1] = Uinf[1];
		}
		
		ComputeOmegaBar(const ComputeOmegaBar& c): eps(c.eps), b2sum(c.b2sum), b2nonempty(c.b2nonempty)
		{
			Uinf[0] = c.Uinf[0];
			Uinf[1] = c.Uinf[1];
		}
		
		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{	
			const int NWHEELS = Wheel::wheels.size();
			
			bool bNonEmpty = false;
			
			for(int iwheel = 0; iwheel<NWHEELS; iwheel++)
			{
				const Wheel& mywheel = *Wheel::wheels[iwheel];
				
				if(FillBlocks::_is_touching(eps, mywheel, info))
				{
					double s = 0;
					
					for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
						for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
						{
							Real p[2];
							info.pos(p, ix, iy);
							
							const Real Xs = mywheel.sample(p[0], p[1]);
							bNonEmpty |= Xs>0;
							
							s += Xs*((b(ix, iy).u[1]+Uinf[1])*(p[0]-mywheel.xm)-(b(ix, iy).u[0]+Uinf[0])*(p[1]-mywheel.ym));
						}
					
					assert(b2sum.find(info.blockID) != b2sum.end());
					assert(b2nonempty.find(info.blockID) != b2nonempty.end());
					
					b2sum[info.blockID][iwheel] = s*info.h[0]*info.h[0];
					b2nonempty[info.blockID] = bNonEmpty;
				}
			}
		}
	};
	
	struct FillVelblocks
	{
		vector<pair< BlockInfo, VelocityBlock *> >& workitems;
		
		FillVelblocks(vector<pair< BlockInfo, VelocityBlock *> >& workitems): workitems(workitems) {}
		
		FillVelblocks(const FillVelblocks& c): workitems(c.workitems) {}
		
		inline void operator()(blocked_range<int> range) const
		{	
			for(int iblock=range.begin(); iblock<range.end(); iblock++)
			{
				BlockInfo info = workitems[iblock].first;
				VelocityBlock * u_desired = workitems[iblock].second;
				
				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						Real p[2];
						info.pos(p, ix, iy);
						
						u_desired->u[0][iy][ix] = 0;// - Uinf[0];
						u_desired->u[1][iy][ix] = 0;// - Uinf[1];
					}	
				
				const int NWHEELS = Wheel::wheels.size();
				
				for(int iwheel = 0; iwheel<NWHEELS; iwheel++)
				{
					const Wheel& mywheel = *Wheel::wheels[iwheel];
					
					const Real xm = mywheel.xm;
					const Real ym =	mywheel.ym;
					const Real av = mywheel.angular_velocity;
					
					if(FillBlocks::_is_touching(0, mywheel, info))
						for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
							for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
							{
								Real p[2];
								info.pos(p, ix, iy);
								
								const Real Xs = mywheel.sample(p[0], p[1]);
								
								if (Xs > 0)
								{
									u_desired->u[0][iy][ix] += - av*(p[1]-ym);// - Uinf[0];
									u_desired->u[1][iy][ix] += + av*(p[0]-xm);// - Uinf[1];
								}
							}				
				}
				
			}
		}
	};
	
	struct NonEmpty
	{
		Real eps;
		map<int, bool>& b2nonempty;
		
		NonEmpty(Real eps, map<int, bool>& b2nonempty): eps(eps), b2nonempty(b2nonempty) { }
		NonEmpty(const NonEmpty& c): eps(c.eps), b2nonempty(c.b2nonempty) { }
		
		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{	
			const int NWHEELS = Wheel::wheels.size();
			
			bool bNonEmpty = false;
			
			for(int iwheel = 0; iwheel<NWHEELS; iwheel++)
			{
				const Wheel& mywheel = *Wheel::wheels[iwheel];
				
				if(FillBlocks::_is_touching(eps, mywheel, info))
				{
					for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
						for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
						{
							Real p[2];
							info.pos(p, ix, iy);
							
							const Real Xs = mywheel.sample(p[0], p[1]);
							bNonEmpty |= Xs>0;							
						}
					
					assert(b2nonempty.find(info.blockID) != b2nonempty.end());
					
					b2nonempty[info.blockID] = bNonEmpty;
				}
			}
		}
	};	
}

I2D_RotatingWheelObstacle::I2D_RotatingWheelObstacle(ArgumentParser & parser, Grid<W,B>& grid, const Real smoothing_length, const Real radius, const Real Uinf[2], I2D_PenalizationOperator& penalization):
I2D_FloatingObstacleOperator(parser, grid, 2*radius, smoothing_length, Uinf, penalization),
grid(grid),  smoothing_length(smoothing_length), radius(radius), penalization(penalization)
{
	this->Uinf[0] = Uinf[0];
	this->Uinf[1] = Uinf[1];
	
	ifstream input("wheels.input");

	if (!input.good())
	{
		const double a = 2*M_PI/3;
		const double d = 1.9*1.5*radius; 
		
		RotatingWheels::Wheel::wheels.push_back(new RotatingWheels::Wheel(0.5, 0.5, 2*radius, -M_PI*0.36));
		RotatingWheels::Wheel::wheels.push_back(new RotatingWheels::Wheel(0.5 + d*cos(a),   0.5 + d*sin(a), 2*radius, -M_PI*-1.3));
		RotatingWheels::Wheel::wheels.push_back(new RotatingWheels::Wheel(0.5 + d*cos(2*a), 0.5 + d*sin(2*a), 2*radius, -M_PI/2*2.34));	
		
		return;
	}
	
	//-------------------------------------- the following part is a request from petros
	
	cout << "Reading wheels information from file ..." ;
	
	vector<double> v;
	
	while (input.good())
	{
		double f;
		input >> f;
		
		v.push_back(f);
	}
	
	assert(v.size() % 3 == 0);	
	
	const int NWHEELS = v.size() / 3;

	cout << "done: " << endl;
	
	for(int i=0; i<NWHEELS; i++)
	{
		const double safety = 4*radius;

		const double x = min(1-safety, max(safety, v[3*i]));
		const double y = min(1-safety, max(safety, v[3*i+1]));
		const double alpha = v[3*i+2]*M_PI/180;
		
		cout << "New wheel (x, y, alpha) = " << x << ", " << y << ", " << alpha << " degrees." << endl;
		
		RotatingWheels::Wheel::wheels.push_back(new RotatingWheels::Wheel(x, y, 2*radius, alpha));
	}
	
	//-------------------------------------- the preceding part is a request from petros
}

void I2D_RotatingWheelObstacle::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	for(vector<RotatingWheels::Wheel*>::const_iterator it = RotatingWheels::Wheel::wheels.begin(); it!= RotatingWheels::Wheel::wheels.end(); ++it)
		assert(*it != NULL);
	
	RotatingWheels::FillBlocks fill(smoothing_length);
	block_processing.process(vInfo, coll, fill);
}

void I2D_RotatingWheelObstacle::restart(const double t, string)
{
	RotatingWheels::Wheel::restart();
}

void I2D_RotatingWheelObstacle::save(const double t, string)
{
	RotatingWheels::Wheel::save();
}

void I2D_RotatingWheelObstacle::computeDesiredVelocity(const double t)
{
	characteristic_function();
	
	const int NWHEELS = RotatingWheels::Wheel::wheels.size();
	
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	for	(map<int, const VelocityBlock *>::iterator it = desired_velocity.begin(); it!= desired_velocity.end(); it++)
	{
		assert(it->second != NULL);
		VelocityBlock::deallocate(it->second);
	}
	
	map<int, bool> nonempty;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		nonempty[it->blockID] = false;
	
	RotatingWheels::NonEmpty compute(smoothing_length, nonempty);
	block_processing.process(vInfo, coll, compute);
	
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
#ifndef NDEBUG
		else 
		{
			const BlockCollection<FluidBlock2D>& coll = grid.getBlockCollection();
			FluidBlock2D& b = coll[it->blockID];
			
			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					assert(b(ix, iy).tmp == 0);
		}
#endif
	}
	
	RotatingWheels::FillVelblocks fillvelblocks(velblocks);
	tbb::parallel_for(blocked_range<int>(0, velblocks.size()), fillvelblocks, auto_partitioner());
	
	penalization.set_desired_velocity(&desired_velocity);
}

void I2D_RotatingWheelObstacle::update(const double dt, const double t, string, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	const int NWHEELS = RotatingWheels::Wheel::wheels.size();
	
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	map<int, vector<double> > integrals;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		integrals[it->blockID] = vector<double>(NWHEELS);
	
	map<int, bool> nonempty;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		nonempty[it->blockID] = false;
	
	RotatingWheels::ComputeOmegaBar compute(smoothing_length, integrals, Uinf, nonempty);
	block_processing.process(vInfo, coll, compute);
	
	vector<double> omegabar(NWHEELS);
	for(map<int, vector< double> >::const_iterator it= integrals.begin(); it!=integrals.end(); ++it)
		for(int i=0; i< NWHEELS; i++)
			omegabar[i] += (it->second)[i];
	
	RotatingWheels::Wheel::update_angle(omegabar, dt, t);
}
