/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include "I2D_FluidMediatedInteractions.h"
#include "I2D_AdvectionOperator_Particles.h"
#include "I2D_VelocitySolver_Mani.h"
#include "I2D_PotentialSolver_Mattia.h"

#include "I2D_FloatingObstacleVector.h"
#include "I2D_StefanFishMorph.h"
#include "I2D_CStartLarva.h"
#include <limits>

//#define _FTLE_

#ifdef _I2D_MPI_
#include "I2D_VelocitySolverMPI_Mani.h"
#endif

static const int maxParticleStencil[2][3] = {
		-3, -3, 0,
		+4, +4, +1
};

class FitnessAcceleration
{
  bool started;
  Real xStart, yStart;
  I2D_FloatingObstacleOperator * floatingObstacle;

public:
  FitnessAcceleration() : started(false), xStart(0.0), yStart(0.0) {}

  void fitness(I2D_FloatingObstacleOperator * floatingObstacle, double time, double startTime, double endTime)
  {
    I2D_FloatingObstacleVector * vec = static_cast<I2D_FloatingObstacleVector*>(floatingObstacle);

    string object("I2D_CStartLarva");
    vector<I2D_FloatingObstacleOperator *> &agents = (vec->data)[object];
    assert(agents.size()==1);
    vector<I2D_FloatingObstacleOperator *>::iterator it=agents.begin();
    I2D_CStartLarva * b = static_cast<I2D_CStartLarva*>(*it);

    assert(startTime<=endTime);

    if(time>=startTime && !started)
      {
	xStart = b->shape->xm;
	yStart = b->shape->ym;
	started = true;
      }

    if(time>=startTime && time<=endTime)
      {
	const Real xx = b->shape->xm;
	const Real yy = b->shape->ym;
	const Real dist = sqrt( (xStart-xx)*(xStart-xx) + (yStart-yy)*(yStart-yy) );
	const Real speed = dist;///(time-startTime);

	// Dump fitness                                                                                                                                                                                                                                                                                                                                     
	FILE * fitnessFile = fopen("fitness","w");
	assert(fitnessFile!=NULL);
	fprintf(fitnessFile,"%10.10e", speed);
	fclose(fitnessFile);
      }

    if(time>=endTime)
      exit(0);
  }
};

I2D_FluidMediatedInteractions::I2D_FluidMediatedInteractions(const int argc, const char ** argv): I2D_FlowPastFloatingObstacle(argc,argv)
{
	bUSEPOTENTIAL = parser("-fmm-potential").asBool();

	if(bUSEPOTENTIAL)
	  {
	    if(sFMMSOLVER == "velocity")
	      potsolver = new I2D_PotentialSolver_Mattia(*grid, parser);
	    else
	      {
		printf("No potential solver MPI!!\n");
		abort();
	      }
	  }

	// Extra parsing
	std::string sFACTORY = parser("-factory").asString();
	assert(sFACTORY != "");

	// Set overall system's characteristic length (Every shape must have its own characteristic length <= D )
	charLength = D;

	// Reset stuff set by studycases
	assert(floatingObstacle==obstacle);
	if(floatingObstacle!=NULL && floatingObstacle==obstacle){ delete floatingObstacle; floatingObstacle = NULL; obstacle = NULL;}

	if( parser("-obstacle").asString() == "heterogeneous" )
	{
		// Pack objects in shapesArray using factory object
		map< string, vector<I2D_FloatingObstacleOperator *> > shapesArray;

		I2D_ObjectFactory factory(*grid, charLength, XPOS, YPOS, epsilon, Uinf, *penalization, LMAX);
		factory.create(parser,shapesArray);

		// Compute characteristic velocity (the highest characteristic velocity among the shapes and also Uinf is considered)
		charVel = factory.getModulusMaxVel();
		assert(charVel>0.0);

		// Compute viscosity given characteristic velocity
		nu = charVel*charLength/RE;
		diffusion->set_viscosity(nu);

		// Allocate heterogeneous collection of objects
		floatingObstacle = new I2D_FloatingObstacleVector(parser,*grid, epsilon, Uinf, *penalization, shapesArray, charLength, charVel);
	}
	else
	{
		printf("Study case not defined!\n");
		abort();
	}

	parser.unset_strict_mode();

	// use vorticity killing?
	bUSEKILLVORT = parser("-usekillvort").asBool();

	if(bUSEKILLVORT)
	{
		// initialize vorticity killing stuff
		KILLVORT = parser("-killvort").asInt();
		killVort = new I2D_KillVortRightBoundaryOperator(*grid, KILLVORT);
	}
	else
	{
		// vorticity killing stuff not used --> defang it
		KILLVORT = 0;
		killVort = NULL;
	}

	// coupled with optimization?
	bUSEOPTIMIZER = parser("-useoptimizer").asBool();
	if(bUSEOPTIMIZER)
	{
		TBOUND = parser("-tbound").asDouble();
		// as long as t<TBOUND the fitness file contains a huge number
		FILE * fitnessFile = fopen("fitness","w");
		assert(fitnessFile!=NULL);
		fprintf(fitnessFile,"%e",std::numeric_limits<Real>::max());
		fclose(fitnessFile);
	}
	else
		TBOUND = 0.0;


	FILE * ppFile = fopen("header.txt","w");
	fprintf(ppFile, "LCFL=%e\n", LCFL);
	fprintf(ppFile, "LMAX=%d\n", LMAX);
	fprintf(ppFile, "RE=%e\n", RE);
	fprintf(ppFile, "charLength=%e\n",charLength);
	fprintf(ppFile, "charVel=%e\n",charVel);
	fprintf(ppFile, "nu=%10.10e\n", nu);
	fflush(ppFile);
	fclose(ppFile);

	// Set obstacle pointer to floating obstacle to inherit
	obstacle = floatingObstacle;

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
			_ic(*grid); // initial condition: set omega, vx, vy and t to 0
			_refine(true);
			_compress(true);
			_dump("initialcondition");
		}
	}
}

I2D_FluidMediatedInteractions::~I2D_FluidMediatedInteractions()
{
	if(floatingObstacle!=NULL)
	{
		delete floatingObstacle;
		floatingObstacle = NULL;
		obstacle = NULL;
	}

	if(killVort!=NULL)
	{
		delete killVort;
		killVort = NULL;
	}
}

void I2D_FluidMediatedInteractions::run()
{
	const tbb::tick_count start_instant = tbb::tick_count::now();

	vector<Real> infoObstacle;

	// cD integration stuff
	Real cD_t = 0.0;
	Real SumcD = 0.0;
	Real dimless_f = 0.0;
	Real T_dimless = 0.0;
	Real dT_dimless = 0.0;
	Real interval = 0.0;
	Real fitness = std::numeric_limits<Real>::max();

	// TODO: or implement or remove repulsion
	Real repulsion = 1.0;

	//	FitnessAcceleration fitness_accel;

	while(true)
	{
		printf("\n\n\n\n------------------ STEP %d ------------------\n",(int)step_id);
		printf("REFINING..\n");
		profiler.push_start("REF");
		_refine(false); // don't use IC for refinement
		profiler.pop_stop();
		printf("DONE WITH REFINEMENT\n");

		for(int i=0; i<ADAPTFREQ; i++)
		{
			printf("INIT STEP\n");

			// Create shape here and load div def into tmp
			profiler.push_start("SHAPE");
			floatingObstacle->create(t);
			profiler.pop_stop();
			printf("DONE WITH SHAPE CREATION\n");

			// TEMPORARY
			//const Real dt = 1e-2;
			//floatingObstacle->update(dt,t);
			//_dump();
			//step_id++;
			//t += dt;
			//if(t>2.0) exit(0);
			// TEMPORARY

			// Kill vorticity at right boundary
			if(bUSEKILLVORT)
			{
				killVort->killVorticity();
				printf("DONE WITH KILLING VORTICITY AT RIGHT BOUNDARY\n");
			}

			// Reconstruct velocity field from vorticity and potential
			profiler.push_start("VEL");
			velsolver->compute_velocity();
			printf("DONE WITH VELOCITY FROM VORTICITY\n");
			if(bUSEPOTENTIAL)
			{
				potsolver->compute_velocity();
				printf("DONE WITH VELOCITY FROM POTENTIAL\n");
			}
			profiler.pop_stop();

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
			profiler.pop_stop();
			printf("DONE WITH PENALIZATION\n");

			floatingObstacle->computeDragAndStuff(t);
			printf("DONE WITH DIAGNOSTICS\n");

			// integrate drag coefficient over time
			if(bUSEOPTIMIZER)
			{
				dimless_f = sqrt(Uinf[0]*Uinf[0]+Uinf[1]*Uinf[1])*2.0/D;
				T_dimless = t*dimless_f;
				// TODO put next line into if statement as soon as cD is not written to file any more
				cD_t = floatingObstacle->getDrag();
				if (T_dimless>TBOUND) // avoid initial transient (T_dimless <= TBOUND)
				{
					dT_dimless = dt*dimless_f;
					SumcD += cD_t*dT_dimless;
					// compute fitness function and write it to file
					interval = T_dimless-TBOUND;
					fitness = (1.0/interval)*(SumcD/repulsion);
					FILE * fitnessFile = fopen("fitness","w");
					assert(fitnessFile!=NULL);
					fprintf(fitnessFile,"%e", fitness);
					fclose(fitnessFile);
				}

				// temporary: write total drag coefficient to file
				FILE * DragFile = fopen("cD_tot.txt", t == 0.0 ? "w" : "a");
				assert(DragFile!=NULL);
				fprintf(DragFile,"%e\t%e\n", T_dimless, cD_t);
				fclose(DragFile);
			}

			profiler.push_start("DIFF");
			diffusion->perform_timestep(dt);
			profiler.pop_stop();
			printf("DONE WITH DIFFUSION\n");

			profiler.push_start("ADV");
			advection->perform_timestep(dt);
			profiler.pop_stop();
			printf("DONE WITH ADVECTION\n");

			floatingObstacle->update(dt,t);

			//                        fitness_accel.fitness(floatingObstacle, t, 0.0, TEND);

                        t = tnext;
                        step_id++;

			if(step_id%SAVEFREQ==0)
			{
				printf("DUMPING...\n");
				_dump(); 
				printf("DONE WITH DUMPING\n");
			}

#ifdef _FTLE_
			if(tnext >= tnext_dump)
			{
				printf("FTLE SAVING...\n");
				_save();
				printf("DONE FTLE SAVING\n");
			}
#endif

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

		printf("COMPRESS..\n");
		profiler.push_start("COMPRESS");
		_compress(false);
		profiler.pop_stop();
		printf("DONE WITH COMPRESS\n");

		profiler.printSummary();

#ifndef _FTLE_
		// Save should be after compressing otherwise one compressing stage would be lost!
		if(step_id % SAVEFREQ == 0)
		{
			printf("SAVING...\n");
			_save();
			printf("DONE SAVING\n");
		}
#endif
	}
}

