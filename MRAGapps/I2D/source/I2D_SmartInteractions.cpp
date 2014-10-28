/*
 * IF2DSmartInteractions.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: mgazzola
 */

#include "I2D_SmartInteractions.h"
#include "I2D_FloatingObstacleVector.h"
#include "I2D_CarlingFishMorph.h"
#include "RL_QLearning.h"
#include "rng.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

namespace I2D_SmartInteractions_stuff
{

class FitnessSpeed
{
	bool started;
	Real xStart, yStart;
	I2D_FloatingObstacleOperator * floatingObstacle;

public:
	FitnessSpeed() : started(false), xStart(0.0), yStart(0.0) {}

	void fitness(I2D_FloatingObstacleOperator * floatingObstacle, double time, double startTime, double endTime)
	{
		assert(startTime<=endTime);

		I2D_FloatingObstacleVector * vec = static_cast<I2D_FloatingObstacleVector*>(floatingObstacle);

		if(time>=startTime && !started)
		{
			string object("I2D_CarlingFishMorph");
			vector<I2D_FloatingObstacleOperator *> &agents = (vec->data)[object];
			assert(agents.size()==1);
			vector<I2D_FloatingObstacleOperator *>::iterator it=agents.begin();
			I2D_CarlingFishMorph * b = static_cast<I2D_CarlingFishMorph*>(*it);
			xStart = b->shape->xm;
			yStart = b->shape->ym;
			started = true;
		}

		if(time>=startTime && time<=endTime)
		{
			string object("I2D_CarlingFishMorph");
			vector<I2D_FloatingObstacleOperator *> &agents = (vec->data)[object];
			assert(agents.size()==1);
			vector<I2D_FloatingObstacleOperator *>::iterator it=agents.begin();
			I2D_CarlingFishMorph * b = static_cast<I2D_CarlingFishMorph*>(*it);
			const Real xx = b->shape->xm;
			const Real yy = b->shape->ym;
			const Real dist = sqrt( (xStart-xx)*(xStart-xx) + (yStart-yy)*(yStart-yy) );
			const Real speed = dist/(time-startTime);
			const Real sign = ((xx-xStart)>0.0)?1:-1;

			// Dump fitness
			FILE * fitnessFile = fopen("fitness","w");
			assert(fitnessFile!=NULL);
			fprintf(fitnessFile,"%e", sign*speed);
			fclose(fitnessFile);
		}
	}
};

}

I2D_SmartInteractions::I2D_SmartInteractions(const int argc, const char ** argv): I2D_FluidMediatedInteractions(argc,argv), efficiency(NULL), policy(NULL), TSTARTFTLE(0), TENDFTLE(0), TSTARTEFF(0), TENDEFF(0)
{
	// Reset stuff set by studycases
	_dispose();
	//assert(floatingObstacle==obstacle);
	//if(floatingObstacle!=NULL && floatingObstacle==obstacle){ delete floatingObstacle; floatingObstacle = NULL; obstacle = NULL;}

	// Extra parsing
	std::string sFACTORY = parser("-factory").asString();
	TSTARTFTLE = parser("-tStartFTLE").asDouble();
	TENDFTLE = parser("-tEndFTLE").asDouble();
	TSTARTEFF = parser("-tStartEff").asDouble();
	TENDEFF = parser("-tEndEff").asDouble();

	assert(sFACTORY != "");
	assert(TSTARTFTLE >= 0);
	assert(TENDFTLE >= 0);

	_prepareAgents();

	FILE * ppFile = fopen("header.txt","w");
	fprintf(ppFile, "LCFL=%e\n", LCFL);
	fprintf(ppFile, "LMAX=%d\n", LMAX);
	fprintf(ppFile, "RE=%e\n", RE);
	fprintf(ppFile, "charLength=%e\n",charLength);
	fprintf(ppFile, "charVel=%e\n",charVel);
	fprintf(ppFile, "nu=%10.10e\n", nu);
	fflush(ppFile);
	fclose(ppFile);

	// Initial conditions
	if(floatingObstacle!=NULL && obstacle!=NULL)
	{
		if(bRESTART)
		{
			_restart();
			_dump("restartedcondition");
		}
		else
		{
			_ic(*grid);
			_refine(true);
			_compress(true);
			_dump("initialcondition");
		}
	}

	// This must be deleted and reallocated in refresh!! Why? Look at the second input!!
	//efficiency = new I2D_ComputeEfficiency(*grid, floatingObstacle, nu);
}


I2D_SmartInteractions::~I2D_SmartInteractions()
{
	_dispose();
}

void I2D_SmartInteractions::_dispose()
{
	//if(efficiency!=NULL){ delete efficiency; efficiency=NULL; }

	if(floatingObstacle!=NULL)
	{
		delete floatingObstacle;
		floatingObstacle = NULL;
		obstacle = NULL;
	}

	assert(floatingObstacle==obstacle);
	if(floatingObstacle!=NULL && floatingObstacle==obstacle){ delete floatingObstacle; floatingObstacle = NULL; obstacle = NULL;}
	assert(policy==NULL);

	if(policy!=NULL)
	{
		delete policy;
		policy = NULL;
	}
}

void I2D_SmartInteractions::_prepareAgents()
{
	if( parser("-obstacle").asString() == "heterogeneous" )
	{
		// Pack objects in shapesArray using factory object
		map< string, vector<I2D_FloatingObstacleOperator *> > shapesMap;
		I2D_ObjectFactory factory(*grid, charLength, XPOS, YPOS, epsilon, Uinf, *penalization, LMAX);
		factory.create(parser,shapesMap,true,&policy);

		// Compute characteristic velocity (the highest characteristic velocity among the shapes and also Uinf is considered)
		charVel = factory.getModulusMaxVel();
		assert(charVel>0.0);

		// Compute viscosity given characteristic velocity
		nu = charVel*charLength/RE;
		diffusion->set_viscosity(nu);

		// Allocate heterogeneous collection of objects
		floatingObstacle = new I2D_FloatingObstacleVector(parser, *grid, epsilon, Uinf, *penalization, shapesMap, charLength, charVel);
	}
	else
	{
		printf("Study case not defined!\n");
		abort();
	}

	// Set obstacle pointer to floating obstacle to inherit
	obstacle = floatingObstacle;

	assert(obstacle==floatingObstacle);
	assert(floatingObstacle!=NULL);

	floatingObstacle->restartPolicy();
}

void I2D_SmartInteractions::refresh()
{
	_dispose();

	_prepareAgents();

	// Initial conditions
	_ic(*grid);
	_refine(true);
	_compress(true);
	_dump("initialcondition");

	// Reset time, etc
	t = 0.0;

	//efficiency = new I2D_ComputeEfficiency(*grid, floatingObstacle, nu);

	_save();

	printf("Simulation successfully refreshed!!\n");
}

Real I2D_SmartInteractions::_initial_dt(int nsteps)
{
	const Real dt_min = 1e-5;
	const Real dt_max = 1e-2;

	Real time = 0.0;
	for(unsigned int i=0; i<nsteps; i++)
		time += dt_min + (Real)(i)/(Real)(nsteps)*(dt_max - dt_min);

	Real time_current = 0.0;
	Real dt_current = 0.0;
	for(unsigned int i=0; i<nsteps; i++)
	{
		dt_current = dt_min + (Real)(i)/(Real)(nsteps)*(dt_max - dt_min);
		time_current += dt_current;
		if(time_current>t)
			break;
	}

	if( t<=time )
		return dt_current;
	else
		return 1000.0;

	return 1e-3;
}

void I2D_SmartInteractions::run()
{
	assert(obstacle==floatingObstacle);
	assert(floatingObstacle!=NULL);

	const tbb::tick_count start_instant = tbb::tick_count::now();

	vector<Real> infoObstacle;

	const Real mindt = 9.99999e-6;

	bool refreshed = false;

	//I2D_SmartInteractions_stuff::FitnessSpeed fit;

	while(true)
	{
		printf("\n\n\n\n------------------ STEP %d ------------------\n",(int)step_id);
		profiler.push_start("REF");
		_refine(false);
		profiler.pop_stop();

		for(int i=0; i<ADAPTFREQ; i++)
		{
			// Choose an action
			const bool valid = floatingObstacle->choose(t);
			//if(!valid){ printf("NON VALID: REFRESH!\n"); floatingObstacle->savePolicy(); refresh(); break; } // ----> ONLY FOR LEARNING

			// Create shape here and load div def into tmp
			profiler.push_start("SHAPE");
			floatingObstacle->create(t);
			profiler.pop_stop();

			// Reconstruct velocity field from vorticity and potential
			profiler.push_start("VEL");
			velsolver->compute_velocity();
			if(bUSEPOTENTIAL)
			{
				potsolver->compute_velocity();
			}
			profiler.pop_stop();

			double tnext, tnext_dump, tend;
			_tnext(tnext, tnext_dump, tend);
			const Real dt = (Real)tnext - t;
			//if(dt<mindt){ printf("TIMESTEP TO SMALL (dt=%e): REFRESH!\n",dt); floatingObstacle->savePolicy(); refresh(); break; } // ----> ONLY FOR LEARNING
			//if(dt<1e-7){ exit(0); } // ----> ONLY FOR OPTIMIZATION

			profiler.push_start("DESVEL");
			floatingObstacle->computeDesiredVelocity(t);
			profiler.pop_stop();
			profiler.push_start("PEN");
			floatingObstacle->characteristic_function();
			penalization->perform_timestep(dt);
			profiler.pop_stop();

			profiler.push_start("DIFF");
			diffusion->perform_timestep(dt);
			profiler.pop_stop();

			profiler.push_start("ADV");
			advection->perform_timestep(dt);
			profiler.pop_stop();

			//profiler.push_start("EFF");
			//efficiency->compute(floatingObstacle,t,dt,TSTARTEFF,TENDEFF);
			//profiler.pop_stop();

			floatingObstacle->update(dt,t);
			floatingObstacle->reward(t);
			floatingObstacle->learn(t);

			//fit.fitness(floatingObstacle, t, 5.0, 6.0);

			t = tnext;
			step_id++;

			if(tnext >= tnext_dump || step_id <= 5)
			{
				printf("DUMPING...\n");
				_dump();
				printf("DONE WITH DUMPING\n");
			}

			if (t >= tend)
			{
				printf("T=TEND=%2.2f reached! (t=%2.2f)\n", TEND, t);
				exit(0);
			}

			{
				const tbb::tick_count now = tbb::tick_count::now();
				const Real nondim_factor = sqrt(Uinf[0]*Uinf[0]+Uinf[1]*Uinf[1])*2/D;
				const Real Tcurr = t*nondim_factor;
				const Real wallclock = (now-start_instant).seconds();
				const int nblocks = grid->getBlocksInfo().size();

				FILE * f = fopen("more-perfmon.txt", step_id == 1 ? "w" : "a");

				if (f!=NULL)
				{
					if (step_id == 0)
						fprintf(f,"stepid\twall-clock[s]\tT[-]\tblocks\tADV-rhs\tDIFF-rhs\n");

					fprintf(f,"%d\t%e\t%e\t%d\t%d\t%d\n", (int)step_id, wallclock, Tcurr, nblocks, advection->get_nofrhs(), diffusion->get_nofrhs());

					fclose(f);
				};
			}

			if(t>=TSTARTFTLE && t<=TENDFTLE)
				_save();
		}

		profiler.push_start("COMPRESS");
		_compress(false);
		profiler.pop_stop();

		//profiler.printSummary();

		// Save should be after compressing otherwise one compressing stage would be lost!
		if(step_id % SAVEFREQ == 0)
		{
			printf("SAVING...\n");
			_save();
			printf("DONE SAVING\n");
		}
	}
}

