/**
 * @file I2D_PassiveTracer.cpp
 * @author Mattia Gazzola
 * @date Mar 1, 2012
 */

#include <iostream>
#include <fstream>
#include "I2D_PassiveTracer.h"
#include "I2D_ParticleBlockLab.h"
#include "I2D_TracerAdvection_RK.h"

// Constructor
// @param grid
// @param _xm
// @param _ym
// @param D
// @param eps
// @param Uinf
// @param penalization
I2D_PassiveTracer::I2D_PassiveTracer(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization): I2D_FloatingObstacleOperator(parser, grid, D, eps, Uinf, penalization)
{
	this->Uinf[0] = Uinf[0];
	this->Uinf[1] = Uinf[1];

	particles.clear();
	block2particles.clear();

	vector<Real> pos;
	pos.push_back(_xm);
	pos.push_back(_ym);
	particles[0] = pos;
}

// Destructor
I2D_PassiveTracer::~I2D_PassiveTracer()
{
}

// Update position and write pertinent info to file
// @param dt
// @param t
// @param filename
void I2D_PassiveTracer::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();
	BoundaryInfo& binfo = grid.getBoundaryInfo();

	// Create tree
	vector<I3> leaves;
	for(int i=0; i<vInfo.size(); i++)
	{
		I3 node(vInfo[i].index[0],vInfo[i].index[1],vInfo[i].level);
		leaves.push_back(node); //-----> input to tree
		block2particles[node] = vector<long int>();
	}

	// Instantiate octree and perform neighbors search
	const int LMAX =  grid.getCurrentMaxLevel();
	double origin[2] = {0.0,0.0};
	const double width = 1.0;
	QuadTree tree(LMAX+1,width,origin);
	tree.split(leaves);

	assert(particles.size()>0);
	for(map<long int, vector<Real> >::iterator it=particles.begin(); it!=particles.end(); ++it)
	{
		assert(it->second.size()==2);
		const double p[2] = {it->second[0],it->second[1]};
		I3 myblock = tree.locateCellIndex(p);

#ifndef NDEBUG
		const Real dx = 1.0/pow(2.0, myblock.i[2]); /// dx of the block
		Real startBlock[2] = {myblock.i[0]*dx, myblock.i[1]*dx}; /// actual location of the start of the block
		Real endBlock[2] = {startBlock[0]+dx, startBlock[1]+dx}; /// actual location of the end limits of the block
		const bool isInBlock = ( p[0]>=startBlock[0] && p[0]<endBlock[0] && p[1]>=startBlock[1] && p[1]<endBlock[1] );
		assert(isInBlock);
#endif

		(block2particles[myblock]).push_back(it->first);
	}

	// Update passive tracer positions
	I2D_TracerAdvection_RK advect(particles, block2particles, dt, Uinf);
	block_processing.process<I2D_ParticleBlockLab>(vInfo, coll, binfo, advect);

	const Real x = particles[0][0];
	const Real y = particles[0][1];

	/// File I/O
	FILE * ppFile = NULL;
	assert(filename!=std::string());
	ppFile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
	assert(ppFile != NULL);
	fprintf(ppFile, "%e %e %e\n", t, x, y);
	fclose(ppFile);
}

// Save restart to file
// @param t
// @param filename
void I2D_PassiveTracer::save(const double t, string filename)
{
	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		ppFile = fopen("restart_I2D_PassiveTracer.txt", "w");
		assert(ppFile != NULL);
	}
	else
	{
		ppFile = fopen(filename.c_str(), "w");
		assert(ppFile != NULL);
	}

	const Real x = particles[0][0];
	const Real y = particles[0][1];
	fprintf(ppFile, "x: %20.20e\n", x);
	fprintf(ppFile, "y: %20.20e\n", y);
	fclose(ppFile);
}


// Read restart conditions
// @param t
// @param filename
void I2D_PassiveTracer::restart(const double t, string filename)
{
	FILE * ppFile = NULL;
	float val;

	if(filename==std::string()) /// string is not set
	{
		ppFile = fopen("restart_I2D_PassiveTracer.txt", "r");
		assert(ppFile!=NULL);
	}
	else
	{
		ppFile = fopen(filename.c_str(), "r");
		assert(ppFile!=NULL);
	}

	particles.clear();
	block2particles.clear();

	vector<Real> pos;
	fscanf(ppFile, "x: %e\n", &val);
	pos.push_back(val);
	fscanf(ppFile, "y: %e\n", &val);
	pos.push_back(val);
	particles[0] = pos;
	printf("PassiveTracer restart (x, y) is: (%e, %e)\n", particles[0][0], particles[0][1]);
}


// Backup to last restart position?
// @param t
// @param filename
void I2D_PassiveTracer::refresh(const double t, string filename)
{
	ifstream filestream;
	if(filename == std::string())
	{
		filestream.open("restart_I2D_PassiveTracer.txt");
		if (!filestream.good())
		{
			cout << "File not found. Exiting now." << endl;
			exit(-1);
		}
	}
	else
	{
		filestream.open(filename.c_str());
		if (!filestream.good())
		{
			cout << "File not found. Exiting now." << endl;
			exit(-1);
		}
	}

	/// Open file and count lines
	int c = 0;
	string line;
	while (getline(filestream, line)) c++;
	filestream.close();

	FILE * f = NULL;
	if (filename == std::string())
	{
		f = fopen("restart_I2D_PassiveTracer.txt", "r");
		assert(f != NULL);
	}
	else
	{
		f = fopen(filename.c_str(), "r");
		assert(f != NULL);
	}

	// data stored in dataVect until t = t_restart
	int N = 0;
	vector < vector<float> > dataVect;
	for (int i = 0; i < c; i++)
	{
		vector<float> row;
		float variable = 0.0;
		fscanf(f, "%f", &variable);
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
			N = i;
		}
		else
		{
			break;
		}
		dataVect.push_back(row);
	}
	fclose(f);
	f = NULL;

	// dataVect is copied in a new "shape" file
	if(filename==std::string()) {
		f = fopen("restart_I2D_PassiveTracer.txt", "w");
		assert(f!=NULL);
	}
	else {
		f = fopen(filename.c_str(), "w");
		assert(f != NULL);
	}
	for (int i = 0; i < N; i++)
		fprintf(f, "%e %e %e \n", dataVect[i][0], dataVect[i][1], dataVect[i][2]);

	fclose(f);
}


