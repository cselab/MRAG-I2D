/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include "I2D_FloatingCylinder.h"
#include <iostream>
#include <fstream>

I2D_FloatingCylinder::Cylinder::Cylinder(): 
angle(0), D(0.1), xm(0.5), ym(0.5), angular_velocity(0), vx(0), vy(0), vdefx(0), vdefy(0), rho(1.0)
{
	m = rho*M_PI*(D/2.0)*(D/2.0);
	J = 0.5*m*(D/2.0)*(D/2.0);
}

I2D_FloatingCylinder::Cylinder::Cylinder(Real xm, Real ym, Real D, Real angle_rad): 
angle(angle_rad), D(D), xm(xm), ym(ym), angular_velocity(0), vx(0), vy(0), vdefx(0), vdefy(0), rho(1.0)
{
	m = rho*M_PI*(D/2.0)*(D/2.0);
	J = 0.5*m*(D/2.0)*(D/2.0);
}

Real I2D_FloatingCylinder::Cylinder::_mollified_heaviside(const double dist, const double eps) const
{
	//Positive outside/negative inside
	const double alpha = M_PI*min(1., max(0., (dist+0.5*eps)/eps));			
	return 0.5+0.5*cos(alpha);
}

void I2D_FloatingCylinder::Cylinder::update_all(double dt,double t)
{
	xm += vx*dt;
	ym += vy*dt;			
	angle += angular_velocity*dt;
}

void I2D_FloatingCylinder::Cylinder::restart(FILE * f)
{
	float val;

	fscanf(f, "xm: %e\n", &val);
	xm = val;
	printf("Cylinder::restart(): xm is %e\n", xm);

	fscanf(f, "ym: %e\n", &val);
	ym = val;
	printf("Cylinder::restart(): ym is %e\n", ym);

	fscanf(f, "vx: %e\n", &val);
	vx = val;
	printf("Cylinder::restart(): vx is %e\n", vx);

	fscanf(f, "vy: %e\n", &val);
	vy = val;
	printf("Cylinder::restart(): vy is %e\n", vy);

	fscanf(f, "angular_velocity: %e\n", &val);
	angular_velocity = val;
	printf("Cylinder::restart(): angular_velocity is %e\n", angular_velocity);

	fscanf(f, "angle: %e\n", &val);
	angle = val;
	printf("Cylinder::restart(): angle is %e\n", angle);
}

void I2D_FloatingCylinder::Cylinder::save(FILE * f) const
{
	fprintf(f, "xm: %20.20e\n", xm);
	fprintf(f, "ym: %20.20e\n", ym);
	fprintf(f, "vx: %20.20e\n", vx);
	fprintf(f, "vy: %20.20e\n", vy);
	fprintf(f, "angular_velocity: %20.20e\n", angular_velocity);
	fprintf(f, "angle: %20.20e\n", angle);
}

Real I2D_FloatingCylinder::Cylinder::sample(const Real x_, const Real y_, const Real eps) const
{
	const double dist = sqrt( (x_-xm)*(x_-xm) + (y_-ym)*(y_-ym) ) - D/2.0;	
	return _mollified_heaviside(dist, eps);
}

void I2D_FloatingCylinder::Cylinder::bbox(const Real eps, Real xmin[2], Real xmax[2]) const
{
	assert(eps>=0);
	
	xmin[0] = xm - D/2.0;
	xmin[1] = ym - D/2.0;
	xmax[0] = xm + D/2.0;
	xmax[1] = ym + D/2.0;
	
	xmin[0] -= 2*eps;
	xmin[1] -= 2*eps;
	xmax[0] += 2*eps;
	xmax[1] += 2*eps;
	
	assert(xmin[0]<=xmax[0]);
	assert(xmin[1]<=xmax[1]);
}


void I2D_FloatingCylinder::Cylinder::restart()
{
	//FILE *f = fopen("rotatingwheel.restart", "r");
	//for(int i=0; i<shapes.size(); i++)
	//	shapes[i]->restart(f);
	//fclose(f);
}

void I2D_FloatingCylinder::Cylinder::save()
{
	//FILE * f = fopen("rotatingwheel.restart", "w");
	//for(int i=0; i<shapes.size(); i++)
	//	shapes[i]->save(f);
	//fclose(f);
}


namespace FloatingCylinder
{
	struct FillBlocks
	{
		Real eps;
		I2D_FloatingCylinder::Cylinder * shape;
		
		FillBlocks(Real eps, I2D_FloatingCylinder::Cylinder *_shape): eps(eps) { shape = _shape; }
		
		FillBlocks(const FillBlocks& c): eps(c.eps) { shape = c.shape; }
		
		static bool _is_touching(Real eps, const I2D_FloatingCylinder::Cylinder * wheel, const BlockInfo& info) 
		{
			Real min_pos[2], max_pos[2];
			
			info.pos(min_pos, 0,0);
			info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1);
			
			Real bbox[2][2];
			wheel->bbox(eps, bbox[0], bbox[1]);
			
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
		I2D_FloatingCylinder::Cylinder * shape;
		map<int, vector<double> >& b2sum;
		map<int, bool>& b2nonempty;
		
		ComputeAll(double eps, I2D_FloatingCylinder::Cylinder *_shape, map<int, vector<double> >& b2sum, Real Uinf[2], map<int, bool>& b2nonempty): eps(eps), b2sum(b2sum), b2nonempty(b2nonempty)
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
		I2D_FloatingCylinder::Cylinder * shape;
		vector<pair< BlockInfo, VelocityBlock *> >& workitems;
		
		FillVelblocks(vector<pair< BlockInfo, VelocityBlock *> >& workitems, double eps, I2D_FloatingCylinder::Cylinder *_shape):
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
				const Real vdefx = shape->vdefx;
				const Real vdefy = shape->vdefy;
				
				if(FillBlocks::_is_touching(eps, shape, info))
					for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
						for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
						{
							Real p[2];
							info.pos(p, ix, iy);
							
							const Real Xs = shape->sample(p[0], p[1], eps);
							
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



I2D_FloatingCylinder::I2D_FloatingCylinder(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization):
I2D_FloatingObstacleOperator(parser, grid, D, eps, Uinf, penalization)
{
	this->Uinf[0] = Uinf[0];
	this->Uinf[1] = Uinf[1];
	
	shape = new Cylinder(_xm, _ym, D, 0.0);
}

I2D_FloatingCylinder::~I2D_FloatingCylinder()
{
	assert(shape!=NULL);
	if(shape!=NULL)
	{
		delete shape;
		shape = NULL;
	}
}

void I2D_FloatingCylinder::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	FloatingCylinder::FillBlocks fill(eps,shape);
	block_processing.process(vInfo, coll, fill);
}

void I2D_FloatingCylinder::restart(const double t, string filename)
{
	FILE * ppFile = NULL;
	
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("restart_I2D_FloatingCylinder.txt", "r");
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

void I2D_FloatingCylinder::save(const double t, string filename)
{
	FILE * ppFile = NULL;
	
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("restart_I2D_FloatingCylinder.txt", "w");
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

void I2D_FloatingCylinder::refresh(const double t, string filename)
{
	// Open file stream
	ifstream filestream;
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		filestream.open("restart_I2D_FloatingCylinder.txt");
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
		f = fopen("restart_I2D_FloatingCylinder.txt", "r");
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
			const Real J = variable;
			row.push_back(J);
			fscanf(f,"%f",&variable);
			const Real m = variable;
			row.push_back(m);
			fscanf(f,"%f",&variable);
			const Real rho = variable;
			row.push_back(rho);
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
		f = fopen("restart_I2D_FloatingCylinder.txt", "w");
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
		fprintf(f, "%e %e %e %e %e %e %e %e %e %e \n", dataVect[i][0],dataVect[i][1],dataVect[i][2],dataVect[i][3],dataVect[i][4],dataVect[i][5],dataVect[i][6],dataVect[i][7],dataVect[i][8],dataVect[i][9]);
	}
	fclose(f);
}

void I2D_FloatingCylinder::computeDesiredVelocity(const double t)
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
	FloatingCylinder::ComputeAll computeAll(eps, shape, integrals, Uinf, nonempty);
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
	
	FloatingCylinder::FillVelblocks fillvelblocks(velblocks, eps, shape);
	tbb::parallel_for(blocked_range<int>(0, velblocks.size()), fillvelblocks, auto_partitioner());
}

void I2D_FloatingCylinder::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	shape->update_all(dt,t);
	
	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("update_I2D_FloatingCylinder.txt", t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);		
	}
	else
	{
		// If string is set I open the corresponding file
		ppFile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);
	}

	// Write update data
	fprintf(ppFile, "%e %e %e %e %e %e %e %e %e %e\n", t, shape->xm, shape->ym, shape->vx, shape->vy, shape->angle, shape->angular_velocity, shape->J, shape->m, shape->rho);

	// Cloase file
	fclose(ppFile);
}
