/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#include "I2D_FloatingObstacleVector.h"
#include "I2D_Clear.h"

namespace FloatingObstacleVectorStuff
{

struct Mass
{
	Real density;
	map< int, Real>& local_mass;

	Mass(Real density, map< int, Real>& local_mass): density(density), local_mass(local_mass)
	{
	}

	Mass(const Mass& c): density(c.density), local_mass(c.local_mass)
	{
	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& out) const
	{
		const int n = FluidBlock2D::sizeX*FluidBlock2D::sizeY;

		FluidElement2D * const b = &out(0,0);
		Real total = 0;
		for(int i=0; i<n; i++)
			total += b[i].tmp;

		total *= density*info.h[0]*info.h[0];

		map< int, Real>::iterator it = local_mass.find(info.blockID);
		assert(it != local_mass.end());
		it->second = total;
	}
};

}


I2D_FloatingObstacleVector::I2D_FloatingObstacleVector(ArgumentParser & parser, Grid<W,B>& grid, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, map< string, vector<I2D_FloatingObstacleOperator *> > _data, Real charLength, Real charVel):
												I2D_FloatingObstacleOperator(parser, grid, 1, eps, Uinf, penalization), data(_data), arrayVelBlock(NULL), charLength(charLength), charVel(charVel)
{
	for( map< string, vector<I2D_FloatingObstacleOperator *> >::iterator it = data.begin(); it!=data.end(); ++it)
		for( vector<I2D_FloatingObstacleOperator *>::iterator it2 = it->second.begin(); it2!=it->second.end(); ++it2)
			agents.push_back(*it2);
}

I2D_FloatingObstacleVector::~I2D_FloatingObstacleVector()
{	
	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		if((*it)!=NULL)
		{
			delete (*it);
			(*it) = NULL;
		}

	agents.clear();
	data.clear();

	if(policy!=NULL)
		if( (*policy)!=NULL )
		{
			delete (*policy);
			(*policy) = NULL;
		}
}

Real I2D_FloatingObstacleVector::getD() const
{
	printf("The call to this method is for an eterogeneous collection of obstacle, not implemented yet!\n");
	abort();
	return 1.0;
}

Real I2D_FloatingObstacleVector::getDrag() const
{
	Real SumcD = 0.0;
	for( vector<I2D_FloatingObstacleOperator *>::const_iterator it = agents.begin(); it!=agents.end(); ++it)
		SumcD += (*it)->getDrag();

	return SumcD;
}

void I2D_FloatingObstacleVector::characteristic_function()
{
	I2D_Clear cleaner;
	cleaner.clearTmp(grid);

	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->characteristic_function();
}

void I2D_FloatingObstacleVector::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
	if(filename!="")
	{
		printf("What?? I2D_FloatingObstacleVector does not print foffuckffake! \n");
		abort();
	}

	unsigned int counter = 0;
	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
	{
		char buf[500];
		sprintf(buf, "shape_%04d", counter);
		string f(buf);
		(*it)->update(dt,t,f,&data);
		counter++;
	}
}

void I2D_FloatingObstacleVector::create(const double t)
{
	I2D_Clear cleaner;
	cleaner.clearTmp(grid);

	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->create(t);
}

void I2D_FloatingObstacleVector::computeDesiredVelocity(const double t)
{
	// Compute desired velocities single objects
	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->computeDesiredVelocity(t);

	// Clean up obstacle vector desired velocity (TODO: blocks could be allocated/deallocated sync with refinement/compression)
	if(arrayVelBlock!=NULL)
		VelocityBlock::deallocate(arrayVelBlock);

	desired_velocity.clear();

	// Instantiate global map
	map<int , VelocityBlock *> gmap;

	// Pack all obstacles' desired velocities into a global one
	vector< map<int , const VelocityBlock *> > vectorDesiredVelocities;
	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		vectorDesiredVelocities.push_back( (*it)->getDesiredVelocity() );

	// Create set of ids of all involved blocks
	set<int> idx;
	for(vector< map<int , const VelocityBlock *> >::iterator ito = vectorDesiredVelocities.begin(); ito!=vectorDesiredVelocities.end(); ito++)
		for(map<int, const VelocityBlock *>::iterator it = (*ito).begin(); it!=(*ito).end(); it++)
			idx.insert(it->first);

	// Allocate all blocks
	int counter = 0;
	arrayVelBlock = VelocityBlock::allocate(idx.size());
	memset(arrayVelBlock,0,idx.size()*sizeof(VelocityBlock));
	for(set<int>::iterator it = idx.begin(); it!=idx.end(); it++, counter++ )
		gmap[*it] = &arrayVelBlock[counter];

	// Cumulate values
	for(vector< map<int , const VelocityBlock *> >::iterator ito = vectorDesiredVelocities.begin(); ito!=vectorDesiredVelocities.end(); ito++)
		for(map<int, const VelocityBlock *>::iterator it = (*ito).begin(); it!=(*ito).end(); it++)
		{
			VelocityBlock * dest = gmap[it->first];
			const VelocityBlock * src = it->second;

			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
				{
					dest->u[0][iy][ix] += src->u[0][iy][ix];
					dest->u[1][iy][ix] += src->u[1][iy][ix];
				}
		}

	// Register penalization and global desired velocity
	for	(map<int, VelocityBlock *>::iterator it = gmap.begin(); it!= gmap.end(); it++)
		desired_velocity[it->first] = it->second;

	penalization.set_desired_velocity(&desired_velocity);
}

void I2D_FloatingObstacleVector::save(const double t, string filename)
{
	if(filename=="")
	{
		printf("What?? I2D_FloatingObstacleVector does not have the step_id foffuckffake!\n");
		abort();
	}

	// Write the restart files
	{
		unsigned int counter = 0;
		for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		{
			char buf[500];
			sprintf(buf, "restart.shape%04d", counter);
			string f(buf);
			(*it)->save(t,f);
			counter++;
		}
	}

	// Write the numbered restart files for safety
	{
		unsigned int counter = 0;
		for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		{
			char buf[500];
			sprintf(buf, "restart_%s.shape%04d", filename.c_str(),counter);
			string f(buf);
			(*it)->save(t,f);
			counter++;
		}
	}
}

void I2D_FloatingObstacleVector::restart(const double t, string filename)
{
	if(filename!="")
	{
		printf("What?? I2D_FloatingObstacleVector does not restart foffuckffake! \n");
		abort();
	}

	unsigned int counter = 0;
	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
	{
		char buf[500];
		sprintf(buf, "restart.shape%04d", counter);
		string f(buf);
		(*it)->restart(t,f);
		counter++;
	}
}

void I2D_FloatingObstacleVector::refresh(const double t, string filename)
{
	if(filename!="")
	{
		printf("What?? I2D_FloatingObstacleVector does not refresh foffuckffake! \n");
		abort();
	}

	unsigned int counter = 0;
	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
	{
		char buf[500];
		sprintf(buf, "shape_%04d", counter);
		string f(buf);
		(*it)->refresh(t,f);
		counter++;
	}
}

void I2D_FloatingObstacleVector::computeDragAndStuff(const Real time, const Real charLength, const Real charVel)
{
	I2D_Clear cleaner;

	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
	{
		cleaner.clearTmp(grid);
		(*it)->characteristic_function();
		(*it)->computeDragAndStuff(time, this->charLength, this->charVel);
	}
}

std::vector<Real> I2D_FloatingObstacleVector::getMass()
{
	std::vector<Real> result;

	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
	{
		std::vector<Real> tmp = (*it)->getMass();
		result.push_back(tmp[0]);
	}

	I2D_Clear cleaner;
	cleaner.clearTmp(grid);

	characteristic_function();

	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	map<int, Real> local_data;
	for(vector<BlockInfo>::iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		local_data[it->blockID] = 0;

	// Compute mass
	const Real density = 1;
	FloatingObstacleVectorStuff::Mass mass(density,local_data);
	block_processing.process(vInfo, coll, mass);

	Real totalmass = 0;
	for(map< int, Real>::iterator it=local_data.begin(); it!=local_data.end(); it++)
		totalmass += it->second;

	result.push_back(totalmass);

	return result;
}

bool I2D_FloatingObstacleVector::choose(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data )
{
	bool valid = true;
	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		valid *= (*it)->choose(t,&data);

	return valid;
}

void I2D_FloatingObstacleVector::learn(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data, string name)
{
	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->stopTest(t,&data);

	unsigned int counter = 0;
	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
	{
		char ending[100];
		sprintf(ending, "shape_%04d", counter);
		string dummy(ending);
		string name("learning_" + dummy);
		(*it)->learn(t,&data,name);
		counter++;
	}
}

void I2D_FloatingObstacleVector::reward(const double t, map< string, vector<I2D_FloatingObstacleOperator *> > * _data )
{
	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
		(*it)->reward(t,&data);
}

void I2D_FloatingObstacleVector::savePolicy(string name)
{
	unsigned int counter = 0;
	for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
	{
		char ending[100];
		sprintf(ending, "shape_%04d", counter);
		string dummy(ending);
		string name("saveQ_" + dummy);
		(*it)->savePolicy(name);
		counter++;
	}
}

void I2D_FloatingObstacleVector::restartPolicy(string name)
{
	unsigned int counter = 0;
		for( vector<I2D_FloatingObstacleOperator *>::iterator it = agents.begin(); it!=agents.end(); ++it)
	{
		char ending[100];
		sprintf(ending, "shape_%04d", counter);
		string dummy(ending);
		string name("saveQ_" + dummy);
		(*it)->restartPolicy(name);
		counter++;
	}
}

