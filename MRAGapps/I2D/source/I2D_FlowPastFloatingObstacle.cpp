/*
 *  I2D_FlowPastFloatingObstacle.cpp
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include "I2D_FlowPastFloatingObstacle.h"
#include "I2D_AdvectionOperator_Particles.h"
#include "I2D_VelocitySolver_Mani.h"

#include "I2D_CircularObstacleOperator.h"
#include "I2D_RectangularObstacleOperator.h"
#include "I2D_WingObstacleOperator.h"
#include "I2D_RotatingWheelObstacle.h"
#include "I2D_LinkedBodies.h"
#include "I2D_EllipticalObstacleOperator.h"
#include "I2D_FloatingCylinder.h"
#include "I2D_FloatingObstacleVector.h"
#include "I2D_ImposedCylinder.h"
#include "I2D_CarlingFish.h"

#ifdef _I2D_MPI_
#include "I2D_VelocitySolverMPI_Mani.h"
#endif

static const int maxParticleStencil[2][3] = {
		-3, -3, 0,
		+4, +4, +1
};

I2D_FlowPastFloatingObstacle::I2D_FlowPastFloatingObstacle(const int argc, const char ** argv): 
						I2D_FlowPastFixedObstacle(argc,argv), floatingObstacle(NULL)
{
	// Reset stuff set by I2D_FlowPastFixedObstacle	
	if(obstacle!=NULL){ delete obstacle; obstacle = NULL;}

	const Real modv = sqrt(Uinf[0]*Uinf[0]+Uinf[1]*Uinf[1]);

	// Allocate shapes
	if (sOBSTACLE == "floatingCyl")
	{
		map< string, vector<I2D_FloatingObstacleOperator *> > shapesMap;
		shapesMap[sOBSTACLE] = vector<I2D_FloatingObstacleOperator *>();
		shapesMap[sOBSTACLE].push_back( new I2D_FloatingCylinder(parser, *grid, XPOS, YPOS, D, epsilon, Uinf, *penalization) );
		floatingObstacle = new I2D_FloatingObstacleVector(parser, *grid, epsilon, Uinf, *penalization, shapesMap, D, modv); // ObstacleVector enables single object drag evaluation
	}
	else if (sOBSTACLE == "floatingCylArray")
	{
		// Pack objects here
		map< string, vector<I2D_FloatingObstacleOperator *> > shapesMap;
		shapesMap[sOBSTACLE] = vector<I2D_FloatingObstacleOperator *>();
		shapesMap[sOBSTACLE].push_back( new I2D_FloatingCylinder(parser, *grid, 0.3, 0.5, D, epsilon, Uinf, *penalization) );
		shapesMap[sOBSTACLE].push_back( new I2D_FloatingCylinder(parser, *grid, 0.7, 0.5, D, epsilon, Uinf, *penalization) );
		floatingObstacle = new I2D_FloatingObstacleVector(parser, *grid, epsilon, Uinf, *penalization, shapesMap, D, modv); // ObstacleVector enables single object drag evaluation
	}
	else if (sOBSTACLE == "wheels")
	{
		floatingObstacle = new I2D_RotatingWheelObstacle(parser, *grid, epsilon, D/2, Uinf, *penalization);
	}

	if (floatingObstacle != NULL) 
		obstacle = floatingObstacle;

	// Note Mattia:	As in I2D_FlowPastFixedObstacle, this check must be commented out because
	//				the constructor of this class, is used by other study cases with different
	//				obstacles. The same way this class uses the constructor of I2D_FlowPastFixedObstacle
	//if (obstacle == NULL)
	//{
	//	cout << "oops, obstacle not assigned. Aborting now." << endl;
	//	abort();
	//}

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
}

I2D_FlowPastFloatingObstacle::~I2D_FlowPastFloatingObstacle()
{
	if(floatingObstacle!=NULL)
	{
		delete floatingObstacle;
		floatingObstacle = NULL;
		obstacle = NULL;
	}
}

void I2D_FlowPastFloatingObstacle::_restart()
{
	//read status
	{
		FILE * f = fopen("restart.status", "r");
		assert(f != NULL);
		float val = -1;
		fscanf(f, "time: %e\n", &val);
		assert(val>=0);
		t=val;
		int step_id_fake = -1;
		fscanf(f, "stepid: %d\n", &step_id_fake);
		step_id = step_id_fake;
		assert(step_id >= 0);
		fclose(f);
	}

	printf("DESERIALIZATION: time is %f and step id is %d\n", t, (int)step_id);

	//read grid
	IO_Binary<W,B> serializer;
	serializer.Read(*grid, "restart");

	floatingObstacle->restart(t);
	floatingObstacle->refresh(t);
	floatingObstacle->create(t);
	floatingObstacle->computeDesiredVelocity(t);
}

void I2D_FlowPastFloatingObstacle::_save()
{
	obstacle->characteristic_function();

	printf("****SERIALIZING****\n");

	//write status
	{
		FILE * f = fopen("restart.status", "w");
		if (f != NULL)
		{
			fprintf(f, "time: %20.20e\n", t);
			fprintf(f, "stepid: %d\n", (int)step_id);
			fclose(f);
		}

		printf( "time: %20.20e\n", t);
		printf( "stepid: %d\n", (int)step_id);
	}

	//write numbered status (extra safety measure)
	{
		string numbered_status;
		char buf[500];
		sprintf(buf, "restart_%07d.status", (int)step_id);
		numbered_status = string(buf);

		FILE * f = fopen(numbered_status.c_str(), "w");
		if (f != NULL)
		{
			fprintf(f, "time: %20.20e\n", t);
			fprintf(f, "stepid: %d\n", (int)step_id);
			fclose(f);
		}
	}

	string numbered_filename;

	{
		char buf[500];
		sprintf(buf, "restart_%07d", (int)step_id);
		numbered_filename = string(buf);
	}

	//write grid
	IO_Binary<W,B> serializer;

	serializer.Write(*grid, "restart");
	serializer.Write(*grid, numbered_filename.c_str()); //akamon no shit now

	printf("****SERIALIZING DONE****\n");

	{
		FILE * f = fopen("history.txt", step_id == 0? "w" : "a");
		if (f!= NULL)
		{
			fprintf(f, "%10.10f %d\n", t, (int)step_id);
			fclose(f);
		}
	}

	string step_id_string;
	{
		char buf[500];
		sprintf(buf, "%07d", (int)step_id);
		step_id_string = string(buf);
	}

	floatingObstacle->save(t,step_id_string);

	floatingObstacle->savePolicy();
}

void I2D_FlowPastFloatingObstacle::run()
{
	printf("/////////////////////////////////////////////////////////////////////////\n");
	printf("////////////          FLOW PAST FLOATING OBSTACLE         ///////////////\n");
	printf("/////////////////////////////////////////////////////////////////////////\n");

	if( Uinf[0]==0.0 && Uinf[1]==0.0 )
	{
		printf("I2D_FlowPastFloatingObstacle::run(): Uinf = 0!\n");
		abort();
	}

	const tbb::tick_count start_instant = tbb::tick_count::now();

	vector<Real> infoObstacle;

	while(true)
	{		
		printf("REFINING..\n");
		profiler.push_start("REF");
		_refine(false);
		profiler.pop_stop();
		printf("DONE WITH REFINEMENT\n");

		for(int i=0; i<ADAPTFREQ; i++)
		{
			printf("INIT STEP\n");

			profiler.push_start("VEL");
			velsolver->compute_velocity();
			profiler.pop_stop();
			printf("DONE WITH COMPUTE VELOCITY\n");

			double tnext, tnext_dump, tend;
			_tnext(tnext, tnext_dump, tend);
			const Real dt = (Real)tnext - t;
			printf("DONE WITH TNEXT\n");			

			profiler.push_start("DESVEL");
			floatingObstacle->computeDesiredVelocity(t);
			profiler.pop_stop();
			profiler.push_start("PEN");
			floatingObstacle->characteristic_function();
			penalization->perform_timestep(dt);
			floatingObstacle->characteristic_function();
			floatingObstacle->getObstacleInfo(infoObstacle);
			profiler.pop_stop();
			printf("DONE WITH PENALIZATION\n");

			profiler.push_start("DIFF");
			diffusion->perform_timestep(dt);
			profiler.pop_stop();
			printf("DONE WITH DIFFUSION\n");

			profiler.push_start("ADV");
			advection->perform_timestep(dt);
			profiler.pop_stop();
			printf("DONE WITH ADVECTION\n");

			floatingObstacle->update(dt,t);
			t = tnext;
			step_id++;

			if(tnext >= tnext_dump || step_id < 20)
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

			printf("END TIME STEP (%d over %d)\n", i, ADAPTFREQ);
		}

		profiler.printSummary();

		if(step_id % SAVEFREQ == 0)
		{
			printf("SAVING...\n");
			_save();
			printf("DONE SAVING\n");

		}

		printf("COMPRESS..\n");
		profiler.push_start("COMPRESS");
		_compress(false);
		profiler.pop_stop();
		printf("DONE WITH COMPRESS\n");

	}
}

