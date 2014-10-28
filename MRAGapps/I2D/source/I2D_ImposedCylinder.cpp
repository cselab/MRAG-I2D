/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#include "I2D_ImposedCylinder.h"

I2D_ImposedCylinder::Cylinder::Cylinder(): 
angle(0), D(0.1), xm(0.5), ym(0.5), angular_velocity(0), vx(0), vy(0), rho(1.0)
{
	m = rho*M_PI*(D/2.0)*(D/2.0);
	J = 0.5*m*(D/2.0)*(D/2.0);
}

I2D_ImposedCylinder::Cylinder::Cylinder(Real xm, Real ym, Real D, Real angle_rad): 
																angle(angle_rad), D(D), xm(xm), ym(ym), angular_velocity(0), vx(0), vy(0), rho(1.0)
{
	m = rho*M_PI*(D/2.0)*(D/2.0);
	J = 0.5*m*(D/2.0)*(D/2.0);
}

Real I2D_ImposedCylinder::Cylinder::_mollified_heaviside(const double dist, const double eps) const
{
	//Positive outside/negative inside
	const double alpha = M_PI*min(1., max(0., (dist+0.5*eps)/eps));			
	return 0.5+0.5*cos(alpha); // outside: 0, inside: 1, between: mollified heaviside
}

/*
const I2D_ImposedCylinder::Cylinder::Cylinder& I2D_ImposedCylinder::Cylinder::operator=(const Cylinder& c)
{
	memcpy(this, &c, sizeof(Cylinder));

	return *this;
}
*/

void I2D_ImposedCylinder::Cylinder::update_all(double dt,double t)
{
	xm += vx*dt;
	ym += vy*dt;			
	angle += angular_velocity*dt;
}

void I2D_ImposedCylinder::Cylinder::restart(FILE * f)
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
	printf("Cylinder::restart(): angular_velocity is %e\n", vx);

	fscanf(f, "angle: %e\n", &val);
	angle = val;
	printf("Cylinder::restart(): angle is %e\n", vy);
}

void I2D_ImposedCylinder::Cylinder::save(FILE * f) const
{
	fprintf(f, "xm: %20.20e\n", xm);
	fprintf(f, "ym: %20.20e\n", ym);
	fprintf(f, "vx: %20.20e\n", vx);
	fprintf(f, "vy: %20.20e\n", vy);
	fprintf(f, "angular_velocity: %20.20e\n", angular_velocity);
	fprintf(f, "angle: %20.20e\n", angle);
}

Real I2D_ImposedCylinder::Cylinder::sample(const Real x_, const Real y_, const Real eps) const
{
	const double dist = sqrt( (x_-xm)*(x_-xm) + (y_-ym)*(y_-ym) ) - D/2.0;	// negative = inside cylinder, positive = outside cylinder
	return _mollified_heaviside(dist, eps);
}

void I2D_ImposedCylinder::Cylinder::bbox(const Real eps, Real xmin[2], Real xmax[2]) const // calculate coordinates of a box around the cylinder
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


void I2D_ImposedCylinder::Cylinder::restart()
{
	//FILE *f = fopen("rotatingwheel.restart", "r");
	//for(int i=0; i<shapes.size(); i++)
	//	shapes[i]->restart(f);
	//fclose(f);
}

void I2D_ImposedCylinder::Cylinder::save()
{
	//FILE * f = fopen("rotatingwheel.restart", "w");
	//for(int i=0; i<shapes.size(); i++)
	//	shapes[i]->save(f);
	//fclose(f);
}


namespace ImposedCylinder
{
struct FillBlocks
{
	Real eps;
	I2D_ImposedCylinder::Cylinder * shape;

	FillBlocks(Real eps, I2D_ImposedCylinder::Cylinder *_shape): eps(eps) { shape = _shape; }

	FillBlocks(const FillBlocks& c): eps(c.eps) { shape = c.shape; }

	static bool _is_touching(Real eps, const I2D_ImposedCylinder::Cylinder * wheel, const BlockInfo& info) // true if cylinder coincides with block
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
					info.pos(p, ix, iy); // returns physical coordinates of each point within the box

					b(ix,iy).tmp = max( shape->sample(p[0], p[1], eps), b(ix,iy).tmp );
				}

			bEmpty = false;
		}
	}
};


struct GetNonEmpty
{
	Real eps;
	I2D_ImposedCylinder::Cylinder * shape;
	map<int, bool>& b2nonempty;

	GetNonEmpty(Real eps, I2D_ImposedCylinder::Cylinder *_shape, map<int, bool>& b2nonempty): eps(eps), b2nonempty(b2nonempty)
	{
		shape = _shape;
	}

	GetNonEmpty(const GetNonEmpty& c): eps(c.eps), b2nonempty(c.b2nonempty)
	{
		shape = c.shape;
	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{
		bool bNonEmpty = false;

		if(FillBlocks::_is_touching(eps, shape, info))
		{
			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
				{
					Real p[2];
					info.pos(p, ix, iy);

					const Real Xs = shape->sample(p[0], p[1], eps);
					bNonEmpty |= Xs>0;
				}

			assert(b2nonempty.find(info.blockID) != b2nonempty.end());

			b2nonempty[info.blockID] = bNonEmpty;
		}
	}
};

struct FillVelblocks
{
	double eps;
	I2D_ImposedCylinder::Cylinder * shape;
	vector<pair< BlockInfo, VelocityBlock *> >& workitems;

	FillVelblocks(vector<pair< BlockInfo, VelocityBlock *> >& workitems, double eps, I2D_ImposedCylinder::Cylinder *_shape):
		workitems(workitems), eps(eps)
	{
		shape = _shape;
	}

	FillVelblocks(const FillVelblocks& c): workitems(c.workitems), eps(c.eps)
	{
		shape = c.shape;
	}

	void operator()(blocked_range<int> range) const
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

			if(FillBlocks::_is_touching(eps, shape, info))
				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						Real p[2];
						info.pos(p, ix, iy);

						const Real Xs = shape->sample(p[0], p[1], eps);

						u_desired->u[0][iy][ix] = (Xs>0)?(- av*(p[1]-ym) + vx):0.0;
						u_desired->u[1][iy][ix] = (Xs>0)?(+ av*(p[0]-xm) + vy):0.0;
					}
		}
	}
};
}


I2D_ImposedCylinder::I2D_ImposedCylinder(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real _vx, const Real _vy, const Real _D, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization):
																I2D_FloatingObstacleOperator(parser, grid, D, eps, Uinf, penalization), vx_imposed(_vx), vy_imposed(_vy)
{
	this->Uinf[0] = Uinf[0];
	this->Uinf[1] = Uinf[1];

	shape = new Cylinder(_xm, _ym, _D, 0.0);
}

I2D_ImposedCylinder::~I2D_ImposedCylinder()
{
	assert(shape!=NULL);
	if(shape!=NULL)
	{
		delete shape;
		shape = NULL;
	}
}

void I2D_ImposedCylinder::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	ImposedCylinder::FillBlocks fill(eps,shape);
	block_processing.process(vInfo, coll, fill);
}

void I2D_ImposedCylinder::restart(const double t, string filename)
{
	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("restart_I2D_ImposedCylinder.txt", "r");
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

void I2D_ImposedCylinder::save(const double t, string filename)
{
	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("restart_I2D_ImposedCylinder.txt", "w");
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

void I2D_ImposedCylinder::refresh(const double t, string filename)
{
	// Open file stream
	ifstream filestream;
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		filestream.open("restart_I2D_ImposedCylinder.txt");
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
		f = fopen("restart_I2D_ImposedCylinder.txt", "r");
		assert(f!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		f = fopen(filename.c_str(), "r");
		assert(f!=NULL);
	}

	//data stored in dataVect until t=t_restart
	// TODO: make it automatic!
	vector < vector<float> > dataVect;
	for(int i=0; i<c; i++)
	{
		vector<float> row;
		float variable = 0.0;
		fscanf(f,"%f",&variable);
		const Real dimensionlessTime = variable;
		fscanf(f,"%f",&variable);
		const Real time = variable;
		if(time <= t)
		{
			row.push_back(dimensionlessTime);
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
			fscanf(f,"%f",&variable);
			const Real cD = variable;
			row.push_back(cD);
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
		f = fopen("restart_I2D_ImposedCylinder.txt", "w");
		assert(f!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		f = fopen(filename.c_str(), "w");
		assert(f!=NULL);
	}

	// Print stuff before t<t_restart
	for(unsigned int i=0; i<dataVect.size(); i++)
	{
		for(unsigned int j=0; j<dataVect[i].size(); j++)
			fprintf(f, "%e ", dataVect[i][j]);

		fprintf(f, "\n");
	}

	fclose(f);
}

void I2D_ImposedCylinder::_setMotionPattern(const Real t)
{
	shape->vx = vx_imposed;
	shape->vy = vy_imposed;
	shape->angular_velocity = 0.0;
}

void I2D_ImposedCylinder::computeDesiredVelocity(const double t)
{
	_setMotionPattern(t);

	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	map<int, vector<double> > integrals;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		integrals[it->blockID] = vector<double>(1);

	map<int, bool> nonempty;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		nonempty[it->blockID] = false;

	// Get non-empty blocks
	ImposedCylinder::GetNonEmpty getNonEmpty(eps, shape, nonempty);
	block_processing.process(vInfo, coll, getNonEmpty);

	// Set desired velocities
	for	(map<int, const VelocityBlock *>::iterator it = desired_velocity.begin(); it!= desired_velocity.end(); it++)
	{
		assert(it->second != NULL);
		VelocityBlock::deallocate(it->second);
		it->second = NULL;
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

	ImposedCylinder::FillVelblocks fillvelblocks(velblocks, eps, shape);
	tbb::parallel_for(blocked_range<int>(0, velblocks.size()), fillvelblocks, auto_partitioner());
}

void I2D_ImposedCylinder::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	shape->update_all(dt,t);

	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("update_I2D_ImposedCylinder.txt", t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);		
	}
	else
	{
		// If string is set I open the corresponding file
		ppFile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
		assert(ppFile!=NULL);
	}

	// Write update data
	fprintf(ppFile, "%e %e %e %e %e %e %e %e %e %e %e %e\n", this->dimT, t, shape->xm, shape->ym, shape->vx, shape->vy, shape->angle, shape->angular_velocity, shape->J, shape->m, shape->rho, this->Cd);
	fflush(ppFile);

	// Cloase file
	fclose(ppFile);
}
