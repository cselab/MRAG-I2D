/*
 *  I2D_FlowPastFixedObstacle.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 *  PATER NOSTER, qui es in caelis, sanctificetur nomen tuum. 
 *	Adveniat regnum tuum. 
 *	Fiat voluntas tua, sicut in caelo et in terra. 
 *	Panem nostrum quotidianum da nobis hodie, 
 *	et dimitte nobis debita nostra sicut et nos dimittimus debitoribus nostris. 
 *	Et ne nos inducas in tentationem, 
 *	sed libera nos a malo. 
 *	Amen.
 *
 */
#include <dirent.h>
#include <time.h>
#include "I2D_FTLE.h"
#include "I2D_AdvectionOperator_Particles.h"
#include "I2D_VelocitySolver_Mani.h"

#include "I2D_CircularObstacleOperator.h"
#include "I2D_RectangularObstacleOperator.h"
#include "I2D_WingObstacleOperator.h"
#include "I2D_LinkedBodies.h"
#include "I2D_EllipticalObstacleOperator.h"
#include "I2D_TracerAdvection_RK.h"

#ifdef _I2D_MPI_
#include "I2D_VelocitySolverMPI_Mani.h"
#endif

static const int maxParticleStencil[2][3] = {
		-4, -4, 0,
		+5, +5, +1
};

namespace FTLEStuff
{

struct FTLE
{
	double T;
	map< I3, unsigned int > & code;
	map<long int, vector<Real> > & particles;

	FTLE(double T, map<long int, vector<Real> > & particles, map< I3, unsigned int > & code): T(T), particles(particles), code(code)
	{
	}

	FTLE(const FTLE& c): T(c.T), particles(c.particles), code(c.code)
	{
	}

	inline Real _computeFTLE(const Real J[2][2], const Real T) const
	{
		//maple code to compute lambda max
		const Real t1 = J[1][0];
		const Real t2 = t1 * t1;
		const Real t3 = J[1][1];
		const Real t4 = t3 * t3;
		const Real t5 = J[0][0];
		const Real t6 = t5 * t5;
		const Real t7 = J[0][1];
		const Real t8 = t7 * t7;
		const Real t9 = t2 * t2;
		const Real t16 = t4 * t4;
		const Real t21 = t6 * t6;
		const Real t24 = t8 * t8;
		const Real t29 = t9 + 0.2e1 * t4 * t2 + 0.2e1 * t6 * t2 - 0.2e1 * t2 * t8 + t16 - 0.2e1 * t4 * t6 + 0.2e1 * t8 * t4 + t21 + 0.2e1 * t6 * t8 + t24 + 0.8e1 * t5 * t1 * t7 * t3;
		const Real t29bis = (t29>=0)?t29:0.0;
		const Real t30 = sqrt(t29bis);
		const Real lambda_max = t2 / 0.2e1 + t4 / 0.2e1 + t6 / 0.2e1 + t8 / 0.2e1 + t30 / 0.2e1;
		const Real ftle = log(sqrt(lambda_max))/fabs(T);
		assert(ftle==ftle);
		return ftle;
	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{
		const unsigned int ppb = (B::sizeX+2)*(B::sizeY+2);
		const double factor = 1.0 / (2.0*info.h[0] );
		I3 node(info.index[0],info.index[1],info.level);
		const unsigned int offset = ppb*code[node];

		for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
			{
				Real J[2][2];

				Real eastX = 0.0;
				Real eastY = 0.0;
				Real westX = 0.0;
				Real westY = 0.0;
				Real northX = 0.0;
				Real northY = 0.0;
				Real southX = 0.0;
				Real southY = 0.0;

				// east
				{
					const unsigned int ixx = ix + 1;
					const unsigned int iyy = iy;
					const unsigned int idx = offset + (ixx+1) + (iyy+1)*(B::sizeX+2);
					map<long int, vector<Real> >::iterator it = particles.find(idx);
					assert(it!=particles.end());
					eastX = it->second[0];
					eastY = it->second[1];
					assert(eastX==eastX);
					assert(eastY==eastY);
				}

				// west
				{
					const unsigned int ixx = ix - 1;
					const unsigned int iyy = iy;
					const unsigned int idx = offset + (ixx+1) + (iyy+1)*(B::sizeX+2);
					map<long int, vector<Real> >::iterator it = particles.find(idx);
					assert(it!=particles.end());
					westX = it->second[0];
					westY = it->second[1];
					assert(westX==westX);
					assert(westY==westY);
				}

				// north
				{
					const unsigned int ixx = ix;
					const unsigned int iyy = iy + 1;
					const unsigned int idx = offset + (ixx+1) + (iyy+1)*(B::sizeX+2);
					map<long int, vector<Real> >::iterator it = particles.find(idx);
					assert(it!=particles.end());
					northX = it->second[0];
					northY = it->second[1];
					assert(northX==northX);
					assert(northY==northY);
				}

				// south
				{
					const unsigned int ixx = ix;
					const unsigned int iyy = iy - 1;
					const unsigned int idx = offset + (ixx+1) + (iyy+1)*(B::sizeX+2);
					map<long int, vector<Real> >::iterator it = particles.find(idx);
					assert(it!=particles.end());
					southX = it->second[0];
					southY = it->second[1];
					assert(southX==southX);
					assert(southY==southY);
				}

				J[0][0] = factor*(eastX  - westX);
				J[0][1] = factor*(northX - southX);
				J[1][0] = factor*(eastY  - westY);
				J[1][1] = factor*(northY - southY);

				b(ix,iy).u[0] = _computeFTLE(J,T);
			}
	}
};


}

I2D_FTLE::I2D_FTLE(const int argc, const char ** argv): parser(argc, argv), t_completed(-1.0)
{
	printf("////////////////////////////////////////////////////////////\n");
	printf("////////////            FTLE         ///////////////\n");
	printf("////////////////////////////////////////////////////////////\n");

	// Parse input
	parser.set_strict_mode();
	PATH = parser("-path").asString();
	JUMP = parser("-jump").asInt();
	LMAX = parser("-lmax").asInt();
	RTOL = parser("-rtol").asDouble();
	bUNIFORM = parser("-uniform").asBool();
	Uinf[0] = parser("-uinfx").asDouble();
	Uinf[1] = parser("-uinfy").asDouble();
	DTFTLE = parser("-dtFTLE").asDouble();
	TFTLE = parser("-TFTLE").asDouble();
	TSTARTFTLE = parser("-tStartFTLE").asDouble();
	TENDFTLE = parser("-tEndFTLE").asDouble();
	bRESTART = parser("-restart").asBool();
	parser.save_options();

	// Instantiate grid
	grid = new Grid<W,B>(8,8,1);
	assert(grid != NULL);
	assert(JUMP >= 1);
	assert(LMAX >= 0);
	assert(RTOL > 0);
	assert(DTFTLE > 0);
	assert(TFTLE > 0);
	assert(TSTARTFTLE > 0);
	assert(TENDFTLE > 0);

	// Set this stuff even if we dont really need it
	refiner = new Refiner_BlackList(JUMP,LMAX);
	compressor = new Compressor(JUMP);
	grid->setRefiner(refiner);
	grid->setCompressor(compressor);

	if(bRESTART)
		_restart();
}

I2D_FTLE::~I2D_FTLE()
{
	if(grid!=NULL){ delete grid; grid=NULL; }
	if(refiner!=NULL){ delete refiner; refiner=NULL; }
	if(compressor!=NULL){ delete compressor; compressor=NULL; }
}

void I2D_FTLE::_getdirContent(string dir, vector<string> &files)
{
	struct dirent **filelist = {0};
	int fcount = -1;

	fcount = scandir(dir.c_str(), &filelist, 0, alphasort);

	if(fcount < 0)
		perror(dir.c_str());

	for(int i = 0; i < fcount; i++)
	{
		//printf("%02d: %s\n", i, filelist[i]->d_name);
		files.push_back(filelist[i]->d_name);
		free(filelist[i]);
	}

	free(filelist);
}

void I2D_FTLE::_getdirRestarts(vector<string> &restarts)
{
	restarts.clear();

	vector<string> files;
	vector<string> filesmrg;

	_getdirContent(PATH, files);

	for(int i=0;i<files.size();i++)
	{
		const size_t found1 = files[i].find(".mrg");
		const size_t found2 = files[i].find("_");
		if(found1!=string::npos && found2!=string::npos)
			filesmrg.push_back(files[i].c_str());
	}

	for(int i=0;i<filesmrg.size();i++)
	{
		const size_t found1 = filesmrg[i].find(".mrg");
		restarts.push_back(filesmrg[i].substr(0,found1));
	}
}

void I2D_FTLE::_mapping(map< string, double > &mapping)
{
	mapping.clear();

	vector<string> restarts;
	_getdirRestarts(restarts);

	for(unsigned int i=0; i<restarts.size(); ++i)
		mapping[restarts[i]] = _loadStatus(restarts[i]);

	//for(map< string, double >::iterator it=mapping.begin(); it!=mapping.end(); ++it)
	//	printf("%s ---> %e\n",it->first.c_str(), it->second);
}

double I2D_FTLE::_loadStatus(string restart)
{
	string numbered_status;
	numbered_status = PATH + "/" + restart + ".status";
	FILE * f = fopen(numbered_status.c_str(), "r");
	assert(f != NULL);
	float val = -1;
	fscanf(f, "time: %e\n", &val);
	assert(val>=0);
	const double time = val;
	fclose(f);
	return time;
}

void I2D_FTLE::_loadGrid(string restart)
{
	string numbered_restart;
	numbered_restart = restart;
	IO_Binary<W,B> serializer;
	serializer.Read(*grid, numbered_restart.c_str());
}

void I2D_FTLE::_createDtsSet(const map< string, double > & subset, vector<Real> & dts)
{
	dts.clear();
	unsigned int counter = 0;
	for(map< string, double >::const_iterator it=subset.begin(); it!=subset.end(); ++it)
	{
		const double t1 = it->second;
		it++;
		const double t2 = it->second;
		it--;
		dts.push_back(t2-t1);
		counter++;

		if( counter==(subset.size()-1) )
			break;
	}
}

void I2D_FTLE::_createParticleSet(string restart, map<long int, vector<Real> > & particles, map< I3, unsigned int > & code)
{
	particles.clear();
	code.clear();

	_loadGrid( restart );
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const unsigned int ppb = (B::sizeX+2)*(B::sizeY+2);
	for(unsigned int i=0; i<vInfo.size(); i++)
	{
		I3 node(vInfo[i].index[0],vInfo[i].index[1],vInfo[i].level);
		code[node] = i;

		const unsigned int offset = ppb*i;
		for(int iy=-1; iy<B::sizeY+1; iy++)
			for(int ix=-1; ix<B::sizeX+1; ix++)
			{
				const unsigned int idx = (ix+1) + (iy+1)*(B::sizeX+2);
				Real p[2];
				vInfo[i].pos(p,ix,iy);
				assert(p[0]==p[0]);
				assert(p[1]==p[1]);
				vector<Real> pos;
				pos.push_back(p[0]);
				pos.push_back(p[1]);
				particles[offset+idx] = pos;
			}
	}
}

void I2D_FTLE::_extractSubset(double tStart, double tEnd, const map< string, double > &mapping, map< string, double > &subset)
{
	subset.clear();
	for(map< string, double >::const_iterator it=mapping.begin(); it!=mapping.end(); ++it)
		if( it->second>=tStart && it->second<tEnd )
			subset[it->first] = it->second;
}

void I2D_FTLE::_FTLE(string initialField, double T, map<long int, vector<Real> > & particles, vector<Real> & dts, map< I3, unsigned int > & code)
{
	_loadGrid( initialField );
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const BlockCollection<B>& coll = grid->getBlockCollection();
	FTLEStuff::FTLE ftle(T,particles,code);
	block_processing.process(vInfo,coll,ftle);
}

double I2D_FTLE::_advectParticles(map<long int, vector<Real> > & particles, vector<Real> & dts, map< string, double > &subset)
{
	unsigned int counter = 0;
	double T = 0.0;
	for(map< string, double >::const_iterator it=subset.begin(); it!=subset.end(); ++it)
	{
		vector<I3> leaves;
		map<I3, vector<long int> > block2particles;

		// Load grid
		profiler.push_start("LOADGRID");
		_loadGrid(PATH + "/" + it->first);
		profiler.pop_stop();

		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		const BlockCollection<B>& coll = grid->getBlockCollection();
		BoundaryInfo& binfo = grid->getBoundaryInfo();

		// Get leaves
		profiler.push_start("GETLEAVES");
		for(int i=0; i<vInfo.size(); i++)
		{
			I3 node(vInfo[i].index[0],vInfo[i].index[1],vInfo[i].level);
			leaves.push_back(node);
			block2particles[node] = vector<long int>();
		}
		profiler.pop_stop();

		// Instantiate octree and perform neighbors search
		profiler.push_start("CREATETREE");
		const int LMAX =  grid->getCurrentMaxLevel();
		double origin[2] = {0.0,0.0};
		const double width = 1.0;
		QuadTree tree(LMAX+1,width,origin);
		tree.split(leaves);
		profiler.pop_stop();

		// Assign particles to blocks
		// THIS PART SHOULD BE PARALLELIZED!
		profiler.push_start("LOCATECELL");
		assert(particles.size()>0);
		for(map<long int, vector<Real> >::const_iterator it=particles.begin(); it!=particles.end(); ++it)
		{
			assert(it->second.size()==2);
			const double p[2] = {it->second[0],it->second[1]};
			if( p[0]>=0.0 && p[0]<1.0 && p[1]>=0.0 && p[1]<1.0 )
			{

				I3 myblock = tree.locateCellIndex(p);
				(block2particles[myblock]).push_back(it->first);


#ifndef NDEBUG
				const Real dx = 1.0/pow(2.0, myblock.i[2]); /// dx of the block
				Real startBlock[2] = {myblock.i[0]*dx, myblock.i[1]*dx}; /// actual location of the start of the block
				Real endBlock[2] = {startBlock[0]+dx, startBlock[1]+dx}; /// actual location of the end limits of the block
				const bool isInBlock = ( p[0]>=startBlock[0] && p[0]<endBlock[0] && p[1]>=startBlock[1] && p[1]<endBlock[1] );
				if(!isInBlock)
				{
					printf("p[0]=%e, p[1]=%e\n",p[0],p[1]);
					printf("startBlock[0]=%e, startBlock[1]=%e\n",startBlock[0],startBlock[1]);
					printf("endBlock[0]=%e, endBlock[1]=%e\n",endBlock[0],endBlock[1]);
				}
				assert(isInBlock);
#endif
			}
		}
		profiler.pop_stop();

		// Update passive tracer positions
		const double dt = dts[counter];
		I2D_TracerAdvection_RK advect(particles, block2particles, dt, Uinf);
		profiler.push_start("ADVECT");
		block_processing.process<I2D_ParticleBlockLab>(vInfo, coll, binfo, advect);
		profiler.pop_stop();
		T += dt;

		counter++;
		if( counter==(subset.size()-1) )
			break;
	}

	return T;
}

void I2D_FTLE::_refine()
{
	if (bUNIFORM) return;

	set<int> boundary_blocks;

	while(true)
	{
		((Refiner_BlackList*)refiner)->set_blacklist(&boundary_blocks);
		const int refinements = Science::AutomaticRefinement<0,3>(*grid, fwt_wuvx, RTOL, LMAX, 1, NULL, (void (*)(Grid<W,B>&))NULL, &boundary_blocks);
		if (refinements == 0) break;
	}
}

void I2D_FTLE::_dump(string filename)
{
	IO_VTKNative<W,B, 4,0> vtkdumper;
	vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
}

void I2D_FTLE::_save(string filename)
{
	printf("****SERIALIZING****\n");
	IO_Binary<W,B> serializer;
	serializer.Write(*grid, filename.c_str());
	printf("****SERIALIZING DONE****\n");
}

void I2D_FTLE::_computeFwdFTLE(double tStart, double tEnd, int idx, map< string, double > & mapping)
{
	map< string, double > subset;
	map< I3, unsigned int > code;
	map<long int, vector<Real> > particles;
	vector<Real> dts;

	_extractSubset(tStart, tEnd, mapping, subset);
	_createDtsSet(subset, dts);
	_createParticleSet(PATH + "/" + (subset.begin())->first, particles, code);

	// Compute preview FTLE
	const Real T1 = _advectParticles(particles,dts,subset);
	_FTLE( PATH + "/" + (subset.begin())->first, T1, particles, dts, code);

	// Refine obtained field
	_refine();
	string ftleTmp("ftleTmp");
	_save(ftleTmp);
	sleep(2);

	// Compute final FTLE
	_createParticleSet(ftleTmp, particles, code);
	const Real T2 = _advectParticles(particles,dts,subset);
	assert(T1==T2);
	sleep(2);
	_FTLE( ftleTmp, T2, particles, dts, code);

	// Print out obtained field
	char buf[500];
	sprintf(buf, "fwdFTLE_%07d", (int)idx);
	string dumpFile(buf);
	_dump(dumpFile);
}

void I2D_FTLE::_save()
{
	FILE * f = fopen("restart.status", "w");
	if (f != NULL)
	{
		fprintf(f, "time: %20.20e\n", t_completed);
		fclose(f);
	}
}

void I2D_FTLE::_restart()
{
	//read status
	{
		FILE * f = fopen("restart.status", "r");
		assert(f != NULL);
		float val = -1;
		fscanf(f, "time: %e\n", &val);
		assert(val>=0);
		t_completed=val;
		fclose(f);
	}

	printf("DESERIALIZATION: time completed is %f\n", t_completed);
}

void I2D_FTLE::run()
{
	// Load data file names
	map< string, double > mapping;
	_mapping(mapping);

	const int FTLE_ITER = floor((TENDFTLE - TSTARTFTLE) / DTFTLE);
	const int idxOffset = floor(TSTARTFTLE / DTFTLE);

	for (int i = 0; i < FTLE_ITER; i++)
	{
		const Real tStart = TSTARTFTLE + (Real) i * DTFTLE;
		const Real tEnd = min(TENDFTLE, tStart+TFTLE);

		if(tStart > t_completed)
		{
			printf("...computing fwd FTLE from time=%f to time=%f (t_completed=%e)\n",tStart,tEnd,t_completed);

			_computeFwdFTLE(tStart, tEnd, i+idxOffset, mapping);

			t_completed = tStart + DTFTLE/2.0;
			_save();

			profiler.printSummary();
		}
	}

	printf("oki I m done!\n");
}
