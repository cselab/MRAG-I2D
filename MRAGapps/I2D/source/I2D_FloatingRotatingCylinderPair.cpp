//
//  I2D_FloatingRotatingCylinderPair.cpp
//  I2D_ROCKS
//
//  Created by Wim van Rees on 02/12/13.
//
//

#include "I2D_FloatingRotatingCylinderPair.h"
#include "I2D_PenalizationOperator.h"

#include <cmath>

I2D_FloatingRotatingCylinderPair::RotatingCylinderPair::RotatingCylinderPair(Real xm, Real ym, Real D, Real width, Real angle_rad, Real GammaLeft, Real GammaRight, Real burgersTime):xm(xm),ym(ym),D(D),width(width),angle(angle_rad),GammaLeft(GammaLeft),GammaRight(GammaRight),burgersTime(burgersTime),vx(0),vy(0),vdefx(0),vdefy(0),rho(1.0)
{
    const double rad = 0.5*D;
    const double area = M_PI*rad*rad; // area of single cylinder
    const double mass = rho*area; // mass of single cylinder
    const double momInertDisk = 0.5*mass*rad*rad; // again for single guy
    
    m = 2*mass; // total mass
	J = 2.0 * (momInertDisk + mass*std::pow(0.5*width,2.)); // parallel axis theorem
    
    std::cout << "MASS = " << m << " INERTIA = " << J << std::endl;
    _setTemporalPreFactor(0.0);
}

Real I2D_FloatingRotatingCylinderPair::RotatingCylinderPair::_mollified_heaviside(const double dist, const double eps) const
{
	//Positive outside/negative inside
	const double alpha = M_PI*min(1., max(0., (dist+0.5*eps)/eps));
	return 0.5+0.5*cos(alpha);
}

void I2D_FloatingRotatingCylinderPair::RotatingCylinderPair::update_all(double dt,double t)
{
	xm += vx*dt;
	ym += vy*dt;
	angle += angular_velocity*dt;
    _setTemporalPreFactor(t);
}

void I2D_FloatingRotatingCylinderPair::RotatingCylinderPair::restart(FILE * f)
{
	float val;
    
	fscanf(f, "xm: %e\n", &val);
	xm = val;
	printf("Cylinder::restart(): xm is %e\n", xm);
    
	fscanf(f, "ym: %e\n", &val);
	ym = val;
	printf("Cylinder::restart(): ym is %e\n", ym);
    
	fscanf(f, "angle: %e\n", &val);
	angle = val;
	printf("Cylinder::restart(): angle is %e\n", angle);
    
	fscanf(f, "vx: %e\n", &val);
	vx = val;
	printf("Cylinder::restart(): vx is %e\n", vx);
    
	fscanf(f, "vy: %e\n", &val);
	vy = val;
	printf("Cylinder::restart(): vy is %e\n", vy);
    
	fscanf(f, "angular_velocity: %e\n", &val);
	angular_velocity = val;
	printf("Cylinder::restart(): angular_velocity is %e\n", angular_velocity);
	
    fscanf(f, "timePrefac: %e\n", &val);
	timePrefac = val;
	printf("Cylinder::restart(): timePrefac is %e\n", timePrefac);
}

void I2D_FloatingRotatingCylinderPair::RotatingCylinderPair::save(FILE * f) const
{
	fprintf(f, "xm: %20.20e\n", xm);
	fprintf(f, "ym: %20.20e\n", ym);
	fprintf(f, "angle: %20.20e\n", angle);
	fprintf(f, "vx: %20.20e\n", vx);
	fprintf(f, "vy: %20.20e\n", vy);
	fprintf(f, "angular_velocity: %20.20e\n", angular_velocity);
	fprintf(f, "timePrefac: %20.20e\n", timePrefac);
}

Real I2D_FloatingRotatingCylinderPair::RotatingCylinderPair::_getBurgersFactor(double tRel)
{
  const Real tol = 1.0e-5;

  Real x0 = tRel;
  Real err = 2.0*tol;
  int iter=0;
  //now we iterate until we satisfy x0 = x - t*u0(x0)
  // in our case x = tRel (input) and t is burgersTime (fixed)
  while(err>tol or iter>1e5)
    {
      const Real newx0 = tRel - burgersTime*std::sin(2.0*M_PI*x0);
      err = std::abs(newx0-x0);
      x0 = newx0;
      iter++;
    }

  return std::sin(2.0*M_PI*x0);//between -1 and +1
}

void I2D_FloatingRotatingCylinderPair::RotatingCylinderPair::_setTemporalPreFactor(double t)
{
    const double g = GammaLeft == 0 ? std::abs(GammaRight) : (GammaRight==0 ? std::abs(GammaLeft) : std::min(std::abs(GammaLeft),std::abs(GammaRight)) ); // want to take the minimum but should still work if one of them is zero
    assert(g>0);
    const double Omega = g / (0.5 * M_PI * D * D); // angular velocity in rad/s
    const double tPerRevolution = 2.0*M_PI/Omega; // time for one revolution (the slowest one)
    const double nRevs = 20.0;

    //timePrefac = 1.0; // impulsively started
    //timePrefac = t < nRevs*tPerRevolution ? std::sin(0.5*M_PI*t/(nRevs*tPerRevolution)) : 1.0; // ramp it up over nRevs revolutions
    //    timePrefac = std::sin(2.0*M_PI*t/(nRevs*tPerRevolution)); // oscillate it with period of nRevs revolutions
    timePrefac = _getBurgersFactor(t/(nRevs*tPerRevolution)); // oscillate it with burgers wave, period is nRevs revolutions
}

Real I2D_FloatingRotatingCylinderPair::RotatingCylinderPair::sample(const Real x_, const Real y_, const Real eps) const
{
    const double CM_left[2] = {
        xm - 0.5*width*std::cos(angle),
        ym - 0.5*width*std::sin(angle)
    };
    
    const double CM_right[2] = {
        xm + 0.5*width*std::cos(angle),
        ym + 0.5*width*std::sin(angle)
    };
    
    const double dist_left = std::sqrt( (x_-CM_left[0])*(x_-CM_left[0]) + (y_-CM_left[1])*(y_-CM_left[1]) ) - 0.5*D;
    
    const double dist_right = std::sqrt( (x_-CM_right[0])*(x_-CM_right[0]) + (y_-CM_right[1])*(y_-CM_right[1]) ) - 0.5*D;
    
	return _mollified_heaviside(std::min(dist_left,dist_right), eps);
}

void I2D_FloatingRotatingCylinderPair::RotatingCylinderPair::bbox(const Real eps, Real xmin[2], Real xmax[2]) const
{
	assert(eps>=0);
	
    const Real rad = 0.5*D;
    
    const Real xdist[4] = {
        xm + 0.5*width*std::cos(angle) + rad,
        xm + 0.5*width*std::cos(angle) - rad,
        xm - 0.5*width*std::cos(angle) + rad,
        xm - 0.5*width*std::cos(angle) - rad
    };
    
    const Real ydist[4] = {
        ym + 0.5*width*std::cos(angle) + rad,
        ym + 0.5*width*std::cos(angle) - rad,
        ym - 0.5*width*std::cos(angle) + rad,
        ym - 0.5*width*std::cos(angle) - rad
    };

    xmin[0] = xdist[0];
    xmin[1] = ydist[0];
    xmax[0] = xdist[0];
    xmax[1] = ydist[0];
    
    for(int i=0; i<4;++i)
    {
        xmin[0] = std::min(xmin[0], xdist[i]);
        xmin[1] = std::min(xmin[1], ydist[i]);
        xmax[0] = std::max(xmax[0], xdist[i]);
        xmax[1] = std::max(xmax[1], ydist[i]);
    }

	xmin[0] -= 2*eps;
	xmin[1] -= 2*eps;
	xmax[0] += 2*eps;
	xmax[1] += 2*eps;
	
	assert(xmin[0]<=xmax[0]);
	assert(xmin[1]<=xmax[1]);
}


namespace FloatingRotatingCylinderPair
{
	struct FillBlocks
	{
		Real eps;
		I2D_FloatingRotatingCylinderPair::RotatingCylinderPair * shape;
		
		FillBlocks(Real eps, I2D_FloatingRotatingCylinderPair::RotatingCylinderPair *_shape): eps(eps) { shape = _shape; }
		
		FillBlocks(const FillBlocks& c): eps(c.eps) { shape = c.shape; }
		
		static bool _is_touching(Real eps, const I2D_FloatingRotatingCylinderPair::RotatingCylinderPair * shape, const BlockInfo& info)
		{
			Real min_pos[2], max_pos[2];
			
			info.pos(min_pos, 0,0);
			info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1);
			
			Real bbox[2][2];
			shape->bbox(eps, bbox[0], bbox[1]);
			
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
			
			if(_is_touching(eps, shape, info))
			{
				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
                for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
                {
                    Real p[2];
                    info.pos(p, ix, iy);
                    
                    b(ix,iy).tmp = max( shape->sample(p[0], p[1], eps), b(ix,iy).tmp );
                }
				
				bEmpty = false;
			}
		}
	};
	
	struct ComputeAll
	{
		Real Uinf[2];
		double eps;
		I2D_FloatingRotatingCylinderPair::RotatingCylinderPair * shape;
		map<int, vector<double> >& b2sum;
		map<int, bool>& b2nonempty;
		
		ComputeAll(double eps, I2D_FloatingRotatingCylinderPair::RotatingCylinderPair *_shape, map<int, vector<double> >& b2sum, Real Uinf[2], map<int, bool>& b2nonempty): eps(eps), b2sum(b2sum), b2nonempty(b2nonempty)
		{
			this->Uinf[0] = Uinf[0];
			this->Uinf[1] = Uinf[1];
			
			shape = _shape;
		}
		
		ComputeAll(const ComputeAll& c): eps(c.eps), b2sum(c.b2sum), b2nonempty(c.b2nonempty)
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
                    
                    const double Xs = shape->sample(p[0], p[1], eps);
                    bNonEmpty |= Xs>0;
                    
                    mass += Xs;
                    J += Xs*((p[0]-shape->xm)*(p[0]-shape->xm) + (p[1]-shape->ym)*(p[1]-shape->ym));
                    vxbar += Xs*(b(ix, iy).u[0]+Uinf[0]);
                    vybar += Xs*(b(ix, iy).u[1]+Uinf[1]);
                    omegabar += Xs*((b(ix, iy).u[1]+Uinf[1])*(p[0]-shape->xm)-(b(ix, iy).u[0]+Uinf[0])*(p[1]-shape->ym));
                }
				
				assert(b2sum.find(info.blockID) != b2sum.end());
				assert(b2nonempty.find(info.blockID) != b2nonempty.end());
				
				b2sum[info.blockID][0] = mass*info.h[0]*info.h[0]*shape->rho;
				b2sum[info.blockID][1] = J*info.h[0]*info.h[0]*shape->rho;
				b2sum[info.blockID][2] = vxbar*info.h[0]*info.h[0]*shape->rho;
				b2sum[info.blockID][3] = vybar*info.h[0]*info.h[0]*shape->rho;
				b2sum[info.blockID][4] = omegabar*info.h[0]*info.h[0]*shape->rho;
				b2nonempty[info.blockID] = bNonEmpty;
			}
		}
	};
	
	struct FillVelblocks
	{
		double eps;
		I2D_FloatingRotatingCylinderPair::RotatingCylinderPair * shape;
		vector<pair< BlockInfo, VelocityBlock *> >& workitems;
		
		FillVelblocks(vector<pair< BlockInfo, VelocityBlock *> >& workitems, double eps, I2D_FloatingRotatingCylinderPair::RotatingCylinderPair *_shape):
		workitems(workitems), eps(eps)
		{
			shape = _shape;
		}
		
		FillVelblocks(const FillVelblocks& c): workitems(c.workitems), eps(c.eps)
		{
			shape = c.shape;
		}
		
		inline void operator()(blocked_range<int> range) const
		{
			for(int iblock=range.begin(); iblock<range.end(); iblock++)
			{
				BlockInfo info = workitems[iblock].first;
				VelocityBlock * u_desired = workitems[iblock].second;
				
				const Real xm = shape->xm;
				const Real ym =	shape->ym;
				const Real av = shape->angular_velocity;
				const Real vx = shape->vx;
				const Real vy = shape->vy;
                // v_theta = Omega * r (Omega is angular velocity in rad/s)
                // Gamma = int omega dx dy // Gamma is circulation (imposed here), omega is the vorticity = 2 * Omega, so:
                // Gamma = 2*Omega*pi*R^2 = 0.5*pi*Omega*D^2
                // Omega = Gamma / (0.5 * pi * D^2)
                const Real fac = 2.0/(M_PI*shape->D*shape->D);
                const Real av_left = shape->timePrefac*shape->GammaLeft * fac;
                const Real av_right = shape->timePrefac*shape->GammaRight * fac;
				
                const double CM_left[2] = {
                    xm + 0.5*shape->width*std::cos(shape->angle),
                    ym + 0.5*shape->width*std::sin(shape->angle)
                };
                
                const double CM_right[2] = {
                    xm - 0.5*shape->width*std::cos(shape->angle),
                    ym - 0.5*shape->width*std::sin(shape->angle)
                };
                
				if(FillBlocks::_is_touching(eps, shape, info))
                for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
                for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
                {
                    Real p[2];
                    info.pos(p, ix, iy);
                    
                    const Real Xs = shape->sample(p[0], p[1], eps);
                    
                    if (Xs > 0)
                    {
                        // need to figure out which one we have
                        const double dist_left = std::sqrt( (p[0]-CM_left[0])*(p[0]-CM_left[0]) + (p[1]-CM_left[1])*(p[1]-CM_left[1]) );
                        
                        const double dist_right = std::sqrt( (p[0]-CM_right[0])*(p[0]-CM_right[0]) + (p[1]-CM_right[1])*(p[1]-CM_right[1]) );
                        
                        const bool isLeft = dist_left < dist_right;
                        const double av_here = isLeft ? av_left : av_right;
                        const double CM_here[2] = {
                            isLeft ? CM_left[0] : CM_right[0],
                            isLeft ? CM_left[1] : CM_right[1]
                        };
                        
                        const double vdefx = - av_here * (p[1] - CM_here[1]);
                        const double vdefy = + av_here * (p[0] - CM_here[0]);
                        
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



I2D_FloatingRotatingCylinderPair::I2D_FloatingRotatingCylinderPair(ArgumentParser & parser, Grid<W,B> & grid, const Real _xm, const Real _ym, const Real D, const Real angle, const Real width, const Real GammaLeft, const Real GammaRight, const Real eps, const Real Uinf[2],I2D_PenalizationOperator& penalization, const int LMAX, const int ID, RL::RL_TabularPolicy ** policy, const int seed):
    I2D_FloatingObstacleOperator(parser, grid, D, eps, Uinf, penalization, ID, policy, seed)
{
	this->Uinf[0] = Uinf[0];
	this->Uinf[1] = Uinf[1];
	parser.unset_strict_mode();
	// some explanation on the next line
	// default of argumentparser is to return 0 so if we argument is not passed the circulation will vary like a sine
	// then whatever we pass to the burgers equation should be between 0 and 1.0/(2*pi) because that's the break time of the wave
	// so to make it easier we normalize by 1/(2pi) here so that we can pass an argument between 0 (sine) and <1 (1 = broken wave)
	const Real burgersTime = parser("-burgersfactor").asDouble()/(2.0*M_PI);
	assert(burgersTime>=0 && burgersTime<1.0);
	shape = new RotatingCylinderPair(_xm, _ym, D, width, angle, GammaLeft,GammaRight,burgersTime);
//	const bool isSharp = parser("-sharp").asBool();
    
}

I2D_FloatingRotatingCylinderPair::~I2D_FloatingRotatingCylinderPair()
{
	assert(shape!=NULL);
	if(shape!=NULL)
	{
		delete shape;
		shape = NULL;
	}
}

void I2D_FloatingRotatingCylinderPair::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	FloatingRotatingCylinderPair::FillBlocks fill(eps,shape);
	block_processing.process(vInfo, coll, fill);
}

void I2D_FloatingRotatingCylinderPair::restart(const double t, string filename)
{
	FILE * ppFile = NULL;
	
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("restart_I2D_FloatingRotatingCylinderPair.txt", "r");
		assert(ppFile!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		ppFile = fopen(filename.c_str(), "r");
		assert(ppFile!=NULL);
	}
	
	// Actual restart
	shape->restart(ppFile);
	
	// Cloase file
	fclose(ppFile);
}

void I2D_FloatingRotatingCylinderPair::save(const double t, string filename)
{
	FILE * ppFile = NULL;
	
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("restart_I2D_FloatingRotatingCylinderPair.txt", "w");
		assert(ppFile!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		ppFile = fopen(filename.c_str(), "w");
		assert(ppFile!=NULL);
	}
	
	// Actual restart
	shape->save(ppFile);
	
	// Cloase file
	fclose(ppFile);
}

void I2D_FloatingRotatingCylinderPair::refresh(const double t, string filename)
{
    // wim: I will not implemented this unless i figure out i need it for something
    std::cout << " REFRESH NOT IMPLEMENTED " << std::endl;
    abort();
}

void I2D_FloatingRotatingCylinderPair::computeDesiredVelocity(const double t)
{
	const int NQUANTITIES = 5;
	
	double mass = 0.0;
	double J = 0.0;
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
	
	// Compute all
	FloatingRotatingCylinderPair::ComputeAll computeAll(eps, shape, integrals, Uinf, nonempty);
	block_processing.process(vInfo, coll, computeAll);
	
	// Sum up contrubutions
	for(map<int, vector< double> >::const_iterator it= integrals.begin(); it!=integrals.end(); ++it)
	{
		mass += (it->second)[0];
		J += (it->second)[1];
		vxbar += (it->second)[2];
		vybar += (it->second)[3];
		omegabar += (it->second)[4];
	}
	
	// Normalization (divide by mass or moment of inertia - in the most complicated cases J varies and must be recomputed on the fly)
	vxbar /= mass;
	vybar /= mass;
	omegabar /= J;
	
	// Set the right vx, vy and angular velocity for each single object
	shape->vx = vxbar;
	shape->vy = vybar;
	shape->angular_velocity = omegabar;
	
    // wim: set mass and J as well
    shape->m = mass;
    shape->J = J;
    
	printf("\n\n");
	printf("vx=%e\n", shape->vx);
	printf("vy=%e\n", shape->vy);
	printf("omegabar=%e\n", shape->angular_velocity);
	printf("mass=%e\n", mass);
	printf("J=%e\n", J);
	printf("\n\n");
	
	// Set desired velocities
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
	
	FloatingRotatingCylinderPair::FillVelblocks fillvelblocks(velblocks, eps, shape);
	tbb::parallel_for(blocked_range<int>(0, velblocks.size()), fillvelblocks, auto_partitioner());
}

void I2D_FloatingRotatingCylinderPair::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
    shape->update_all(dt,t);

	FILE * ppFile = NULL;
    
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("I2D_FloatingRotatingCylinderPair.txt", t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		ppFile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);
	}
    
	// Write update data
	fprintf(ppFile, "%e %e %e %e %e %e %e %e %e %e %e %e\n", t, shape->xm, shape->ym, shape->vx, shape->vy, shape->angle, shape->angular_velocity, shape->J, shape->m, shape->rho, shape->timePrefac*shape->GammaLeft, shape->timePrefac*shape->GammaRight);
    
	// Cloase file
	fclose(ppFile);
}
