//
//  I2D_RotatingCylinderPair.cpp
//  I2D_ROCKS
//
//  Created by Wim van Rees on 05/12/13.
//
//

#include "I2D_RotatingCylinderPair.h"

#include "I2D_AdvectionOperator_Particles.h"
#include "I2D_VelocitySolver_Mani.h"

#include "I2D_FloatingObstacleVector.h"
#include "I2D_FloatingRotatingCylinderPair.h"

#include <limits>

static const int maxParticleStencil[2][3] = {
    -3, -3, 0,
    +4, +4, +1
};


I2D_RotatingCylinderPair::I2D_RotatingCylinderPair(const int argc, const char ** argv): I2D_FlowPastFloatingObstacle(argc,argv)
{
    
    // parse and create the object
    assert(parser("-obstacle").asString() == "rotcyl");
    
    const Real angle = 0.0;
    const Real width = parser("-width").asDouble() * D; // normalize with D
    const Real g1 = parser("-g1").asDouble();
    const Real g2 = -g1;
    std::cout << "WDITH = " << width << std::endl;
    // set scales
    LENGTHSCALE = D;
    CIRCULATION = g1;
    OMEGA = CIRCULATION / (0.5 * M_PI * D * D);
    TIMESCALE = 2.0 * M_PI / OMEGA ; // time for one revolution
    
    // we normalize lambda, so kill and redo penalization operator
    delete penalization;penalization=NULL;
    const double LT_ref = LAMBDA *  0.00098696;
    LAMBDA = LT_ref/TIMESCALE;
    penalization = new I2D_PenalizationOperator(*grid, LAMBDA, Uinf, bRESTART);

    std::cout << "SOME INFO ON THE SCALES" << std::endl;
    std::cout << "LENGTHSCALE : " << D << std::endl;
    std::cout << "CIRCULATION : " << CIRCULATION << std::endl;
    std::cout << "OMEGA : " << OMEGA << std::endl;
    std::cout << "TIMESCALE : " << TIMESCALE << std::endl;
    std::cout << "LAMBDA : " << LAMBDA << std::endl;
    
    // create the obstacle
    // fuck we can only do it through obstacle vector
    // ok fuck it we fake one
    {
      map< string, vector<I2D_FloatingObstacleOperator *> > shapesArray;
        
      I2D_FloatingRotatingCylinderPair * object = new I2D_FloatingRotatingCylinderPair(parser, *grid, XPOS, YPOS, D, angle, width, g1, g2, epsilon, Uinf, *penalization, LMAX);

      shapesArray["I2D_FloatingRotatingCylinderPair"].push_back(object);
      floatingObstacle = new I2D_FloatingObstacleVector(parser,*grid, epsilon, Uinf, *penalization, shapesArray, charLength, charVel);
    }
    
    // set diffusion
    nu = CIRCULATION/RE;
	diffusion->set_viscosity(nu);
        
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
	fprintf(ppFile, "charLength=%e\n",LENGTHSCALE);
	fprintf(ppFile, "charVel=%e\n",LENGTHSCALE/TIMESCALE);
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

I2D_RotatingCylinderPair::~I2D_RotatingCylinderPair()
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

void I2D_RotatingCylinderPair::_tnext(double &tnext, double &tnext_dump, double &tend)
{
	//0. setup
	//1. prepare candidate tnexts
	//2. pick the smallest one, report who wins
    
	//0.
	vector<double> tnext_candidates;
	vector<string> tnext_names;
    
	const double max_dx = (1./B::sizeX)*pow(0.5,grid->getCurrentMinLevel());
	const double min_dx = (1./B::sizeX)*pow(0.5,grid->getCurrentMaxLevel());
	const double max_vel = advection->compute_maxvel();
    
	//1.
	tend = TEND/TIMESCALE;
    
    {
        const int nsteps = RAMP;
        const Real dt_max = LCFL/(2.0*OMEGA); // final dt = LCFL / max_omega, where max_omega = 2.0*OMEGA (rad/s)
        const Real dt_min = 1e-2*dt_max;

        const double dtramp = (step_id >= nsteps) ? 1000.0 : dt_min + (Real)(step_id)/(Real)(nsteps)*(dt_max - dt_min);
        tnext_candidates.push_back(t + dtramp);
        tnext_names.push_back("RAMP");
    }
    
	tnext_candidates.push_back(t + max_dx/max_vel*CFL*0.999);
	tnext_names.push_back("CFL");
    
	if (LAMBDADT>0)
	{
		tnext_candidates.push_back(t + LAMBDADT/LAMBDA);
		tnext_names.push_back("LAMBDA");
	}
    
	if (bPARTICLES)
	{
		tnext_candidates.push_back(t + advection->estimate_largest_dt());
		tnext_names.push_back("LCFL");
	}
    
	if (bFMMSKIP)
	{
		tnext_candidates.push_back(t + 0.95*B::sizeX*min_dx/max_vel);
		tnext_names.push_back("CFL-FMMSKIP");
	}
    
	const bool bDiffusionLTS = true;
	if (bDiffusionLTS)
	{
		tnext_candidates.push_back(t + diffusion->estimate_largest_dt());
		tnext_names.push_back("LTS-Fc");
	}
	else
	{
		tnext_candidates.push_back(t + diffusion->estimate_smallest_dt());
		tnext_names.push_back("GTS-Fc");
	}
    
	{
//		const double delta = 1./(DUMPFREQ*nondim_factor);
          const Real delta = TIMESCALE / DUMPFREQ ; // wim: changed this to DUMPFREQ dumps per non-dim time. dunno what was the shit above.. wtf mattia
		tnext_dump = (DUMPFREQ!=0)?(1.+floor(t/delta))*delta:10000.0;
	}
    
	//2.
	{
		tnext = tnext_candidates[0];
		int imin = 0;
        
		for(int i=1; i<(int)tnext_candidates.size(); i++)
        if (tnext_candidates[i] < tnext)
        {
            imin = i;
            tnext = tnext_candidates[i];
        }
        
		const string dtBound = tnext_names[imin] + " bound";
        
		printf("####################### t is %e, dt is %e, %s,\t{", t, tnext - t, dtBound.c_str());
        
		for(int i=0; i<(int)tnext_candidates.size(); i++)
        printf("%6s:%2.2e ", tnext_names[i].c_str(), tnext_candidates[i] - t);
		printf("}\n");
        
		FILE * f = fopen("report.txt", step_id == 0 ? "w" : "a");
		assert(f!=NULL);
        
		fprintf(f, "####################### t is %e, dt is %e and is bound by: %s", t, tnext - t, dtBound.c_str());
		fprintf(f, "stepid=%d\tT=%e\tDT=%e\t", (int)step_id, tnext/TIMESCALE, (tnext - t)/TIMESCALE);
		for(int i=0; i<(int)tnext_candidates.size(); i++)
        fprintf(f, "%s: %2.2e,\t", tnext_names[i].c_str(), tnext_candidates[i] - t);
        
		fprintf(f, "GTS-Fc: %2.2e,\t", pow(min_dx,2)/(8.0*nu));
        
		fprintf(f, "\n");
		fclose(f);
	}
}


void I2D_RotatingCylinderPair::run()
{
	const tbb::tick_count start_instant = tbb::tick_count::now();
    
	vector<Real> infoObstacle;
    
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
            
			// Reconstruct velocity field from vorticity and potential
			profiler.push_start("VEL");
			velsolver->compute_velocity();
			printf("DONE WITH VELOCITY FROM VORTICITY\n");
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
            
			if(bUSEOPTIMIZER)
			{
                // not yet implemented
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
            
            //                        fitness_accel.fitness(floatingObstacle, t, 3.0);
            
			if(tnext >= tnext_dump)
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
            
            t = tnext;
            step_id++;
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

