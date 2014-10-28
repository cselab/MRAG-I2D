/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include "I2D_MRAGOptimisation.h"
#include "I2D_AdvectionOperator_Particles.h"
#include "I2D_VelocitySolver_Mani.h"

#include "I2D_CircularObstacleOperator.h"
#include "I2D_RectangularObstacleOperator.h"
#include "I2D_WingObstacleOperator.h"
#include "I2D_LinkedBodies.h"
#include "I2D_Ingo.h"
#include "I2D_WingAngleOfAttack.h"
#include "I2D_EllipticalObstacleOperator.h"
#include "I2D_KillVortRightBoundaryOperator.h"

#ifdef _I2D_MPI_
#include "I2D_VelocitySolverMPI_Mani.h"
#endif

static const int maxParticleStencil[2][3] = {
	-3, -3, 0,
	+4, +4, +1
};

I2D_MRAGOptimisation::I2D_MRAGOptimisation(const int argc, const char ** argv): 
I2D_FlowPastFixedObstacle(argc,argv)
{
	if(obstacle!=NULL)
	{
		delete obstacle;
		obstacle = NULL;
	}
				
	KILLVORT = parser("-killvort").asInt();
	sCTRL = parser("-ctrl").asString();
	
	parser.save_options();
	
	assert(sCTRL != "");
	
	Real pos[2] = {0.5,0.5};
	if(XPOS>0 && XPOS<1.0){ pos[0] = XPOS; }
	if(YPOS>0 && YPOS<1.0){ pos[1] = YPOS; }
	
	
	killVort = new I2D_KillVortRightBoundaryOperator(*grid, KILLVORT);
	
	// Open ctrl file and retrieve information
	FILE * ppFile;
	ppFile = fopen(sCTRL.c_str(),"r");
	if(ppFile==NULL){ std::cout << "could not open ctrl file " << sCTRL << "!" << std::endl;}
	char workDir[10000];
	char jobName[10000];
	unsigned int dim = 0;
	unsigned int ID = 0;
	float dummy = 0.0;
	std::vector<Real> parameterCMA;
	fscanf(ppFile,"%s",workDir);
	fscanf(ppFile,"%s",jobName);
	fscanf(ppFile,"%d",&dim);
	for(unsigned int i=0; i<dim-1; i++)
	{
		fscanf(ppFile,"%f",&dummy);
		parameterCMA.push_back((Real)dummy);
	}
	fscanf(ppFile,"%d",&ID);
	fclose(ppFile);
	
	if (sOBSTACLE == "links")
	{
		vector<Real> angles;
		angles.push_back(parameterCMA[0]);
		angles.push_back(parameterCMA[1]);
		angles.push_back(parameterCMA[2]);
		angles.push_back(parameterCMA[3]);
		angles.push_back(parameterCMA[4]);
		angles.push_back(parameterCMA[5]);
		angles.push_back(parameterCMA[6]);
		angles.push_back(parameterCMA[7]);
		angles.push_back(parameterCMA[8]);
		angles.push_back(parameterCMA[9]);
		angles.push_back(parameterCMA[10]);
		angles.push_back(parameterCMA[11]);
		angles.push_back(parameterCMA[12]);
		Real length = 0.03;
		Real width = length/6.0;
		obstacle = new I2D_LinkedBodies(*grid,length,width,pos,angles,epsilon);
	}
	else if (sOBSTACLE == "ingo")
	{
		vector<Real> angles;
		angles.push_back(parameterCMA[0]);
		angles.push_back(parameterCMA[1]);
		angles.push_back(parameterCMA[2]);
		angles.push_back(parameterCMA[3]);
		angles.push_back(parameterCMA[4]);
		Real length = 0.03;
		Real width = length/6.0;
		obstacle = new I2D_Ingo(*grid,length,width,pos,angles,epsilon);
	}
	else if (sOBSTACLE == "angleAttack")
	{
		vector<Real> angles;
		angles.push_back(parameterCMA[0]);
		angles.push_back(parameterCMA[1]);
		angles.push_back(parameterCMA[2]);
		angles.push_back(parameterCMA[3]);
		angles.push_back(parameterCMA[4]);
		angles.push_back(parameterCMA[5]);
		angles.push_back(parameterCMA[6]);
		angles.push_back(parameterCMA[7]);
		angles.push_back(parameterCMA[8]);
		angles.push_back(parameterCMA[9]);
		angles.push_back(parameterCMA[10]);
		angles.push_back(parameterCMA[11]);
		angles.push_back(parameterCMA[12]);
		Real length = 0.03;
		Real width = length/6.0;
		obstacle = new I2D_WingAngleOfAttack(*grid,length,width,pos,angles,epsilon);
	}	
	else if (sOBSTACLE == "NACA")
	{
		obstacle = new I2D_WingObstacleOperator(*grid,parser,pos,epsilon);
	}
	else
	{
		printf("ooops, obstacle not specified. Aborting. \n");
		abort();
	}
	
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

void I2D_MRAGOptimisation::run()
{
	const tbb::tick_count start_instant = tbb::tick_count::now();
	
	Real SumcD = 0.0;
	
	Real SumcL = 0.0;
	
	Real DragAverage, LiftAverage ;
	
	int step_id_fitness = 0;
	
	Real T =0.0;
	
	vector<Real> infoObstacle;
		
	Real cor[2] = {0.0,0.0};
	
	while(true)
	{		
		printf("REFINING..\n");
		profiler.push_start("REF");
		_refine(false);
		profiler.pop_stop();
		printf("DONE WITH REFINEMENT\n");
		
		for(int i=0; i<ADAPTFREQ; i++)
		{
			T= 2*max(fabs(Uinf[0]), fabs(Uinf[1]))*t/D;
			printf("INIT STEP\n");
			
			if(bPARTICLES==1){
				killVort->killVorticity();
				printf("DONE WITH KILLING VORTICITY AT RIGHT BOUNDARY\n");
			}

			profiler.push_start("VEL");
			velsolver->compute_velocity();
			profiler.pop_stop();
			printf("DONE WITH COMPUTE VELOCITY\n");
			
			double tnext, tnext_dump, tend;
			_tnext(tnext, tnext_dump, tend);
			const Real dt = (Real)tnext - t;
			printf("DONE WITH TNEXT\n");			

			profiler.push_start("PEN");
			obstacle->characteristic_function();
			penalization->perform_timestep(dt);
			obstacle->characteristic_function();
			obstacle->getObstacleInfo(infoObstacle);
			profiler.pop_stop();
			printf("DONE WITH PENALIZATION\n");

			if(infoObstacle.size()!=0){ cor[0] = infoObstacle[0]; cor[1] = infoObstacle[1]; }
			penalization->compute_dragandstuff(t,D,cor,"diag");
			if(T >= 30.0){
				SumcD += penalization->getDrag();
				SumcL += penalization->getLift();
				step_id_fitness++;
			}
			printf("DONE WITH DRAG\n");

			profiler.push_start("DIFF");
			diffusion->perform_timestep(dt);
			profiler.pop_stop();
			printf("DONE WITH DIFFUSION\n");

			profiler.push_start("ADV");
			advection->perform_timestep(dt);
			profiler.pop_stop();
			printf("DONE WITH ADVECTION\n");
			
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
				DragAverage = SumcD/(Real)(step_id_fitness);
				LiftAverage = SumcL/(Real)(step_id_fitness);
				FILE * fitnessFile = fopen("fitness","w");
				assert(fitnessFile!=NULL);
				fprintf(fitnessFile,"%e", -(LiftAverage/DragAverage) );
				fclose(fitnessFile);
			
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
