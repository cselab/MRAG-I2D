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

#include "I2D_FlowPastFixedObstacle.h"
#include "I2D_AdvectionOperator_Particles.h"
#include "I2D_VelocitySolver_Mani.h"
#include "I2D_VelocitySolver_Wim.h"

#include "I2D_CircularObstacleOperator.h"
#include "I2D_RectangularObstacleOperator.h"
#include "I2D_WingObstacleOperator.h"
#include "I2D_LinkedBodies.h"
#include "I2D_EllipticalObstacleOperator.h"

#ifdef _I2D_MPI_
#include "I2D_VelocitySolverMPI_Mani.h"
#endif

static const int maxParticleStencil[2][3] = {
		-3, -3, 0,
		+4, +4, +1
};

I2D_FlowPastFixedObstacle::I2D_FlowPastFixedObstacle(const int argc, const char ** argv): 
												parser(argc, argv), t(0), step_id(0),
												velsolver(NULL),obstacle(NULL), penalization(NULL), advection(NULL), diffusion(NULL)
{
	printf("////////////////////////////////////////////////////////////\n");
	printf("////////////            AVE MARIA         ///////////////\n");
	printf("////////////////////////////////////////////////////////////\n");

	bRESTART = parser("-restart").asBool();
	sOBSTACLE = parser("-obstacle").asString();
	sRIGID_INLET_TYPE = parser("-rio").asString();
	bREFINEOMEGAONLY = parser("-refine-omega-only").asBool();
	bFMMSKIP = parser("-fmm-skip").asBool();
	LAMBDADT = parser("-lambdadt").asDouble();
	XPOS = parser("-xpos").asDouble();
	YPOS = parser("-ypos").asDouble();

	if (sOBSTACLE == "")
		sOBSTACLE = "cyl";

	parser.set_strict_mode();

	BPD = parser("-bpd").asInt();
	JUMP = parser("-jump").asInt();
	LMAX = parser("-lmax").asInt();
	TEND = parser("-tend").asDouble();
	RAMP = parser("-ramp").asInt();
	ADAPTFREQ = parser("-adaptfreq").asInt();
	DUMPFREQ = parser("-dumpfreq").asDouble();
	SAVEFREQ = parser("-savefreq").asInt();
	RE = parser("-re").asDouble();
	CFL = parser("-cfl").asDouble();
	LCFL = parser("-lcfl").asDouble();
	RTOL = parser("-rtol").asDouble();
	CTOL = parser("-ctol").asDouble(); 
	LAMBDA = parser("-lambda").asDouble();
	D = parser("-D").asDouble();
	bPARTICLES = parser("-particles").asBool();
	bUNIFORM = parser("-uniform").asBool();
	sFMMSOLVER = parser("-fmm").asString();
	MOLLFACTOR = parser("-mollfactor").asInt();
	Uinf[0] = parser("-uinfx").asDouble();
	Uinf[1] = parser("-uinfy").asDouble();
	FC = parser("-fc").asDouble();
	const bool HILBERT = parser("-hilbert").asBool();

	parser.save_options();

	assert(TEND >= 0.0);
	assert(BPD > 1);
	assert(ADAPTFREQ > 0);
	assert(DUMPFREQ >= 0);
	assert(JUMP >= 1);
	assert(LMAX >= 0);
	assert(RE > 0);
	assert(CFL > 0 && CFL<1);
	assert(LCFL > 0 && LCFL<1);
	assert(RTOL > 0);
	assert(CTOL > 0);
	assert(LAMBDA > 0);
	assert(sFMMSOLVER != "");
	assert(sRIGID_INLET_TYPE != "");
	assert(MOLLFACTOR > 0);
	assert(FC>0.0 && FC<1.0);

#ifdef _I2D_MPI_
	if (sFMMSOLVER == "mpi-velocity")
		velsolver = new I2D_VelocitySolverMPI_Mani(argc, argv);
#endif

	const Real h_min = 1./FluidBlock2D::sizeX*pow(0.5, LMAX);

	const Real modV = sqrt(Uinf[0]*Uinf[0]+Uinf[1]*Uinf[1]);

	//Note by Mattia:	All flow past cases MUST have a defined Uinf. Here the case of Uinf=0 is allowed only in the constructor.
	//					The reason of this is that the constructor of I2D_FlowPastFixedObstacle and I2D_FlowPastFloatingObstacle is
	//					used also in the more general I2D_FluidMediatedInteractions, which may have Uinf=0.
	// 					Therefore I put a check in the I2D_FlowPastFixedObstacle::run() and I2D_FlowPastFloatingObstacle::run() to make sure
	//					that at the moment of running Uinf!=0.
	//					As here Uinf could be zero, in order to set a viscosity!=0 and avoid the assert in the diffusion operator,
	//					I just define the characteristic velocity as below, therefore if Uinf!=0 ----> nu = charVel*D/RE, otherwise nu = D/RE
	const Real charVel = (modV==0.0)?1.0:modV;

	Real pos[2] = {0.5,0.5};
	if(XPOS>0 && XPOS<1.0){ pos[0] = XPOS; }
	if(YPOS>0 && YPOS<1.0){ pos[1] = YPOS; }

	nu = charVel*D/RE;

	if (HILBERT)
	{
		if (bPARTICLES)
			grid = new Grid_Hilbert2D<W,B>(BPD,BPD,1, maxParticleStencil);
		else
			grid = new Grid_Hilbert2D<W,B>(BPD,BPD,1);
	}
	else
	{
		if (bPARTICLES)
			grid = new Grid<W,B>(BPD,BPD,1, maxParticleStencil);
		else
			grid = new Grid<W,B>(BPD,BPD,1);
	}


	assert(grid != NULL);

	refiner = new Refiner_BlackList(JUMP, LMAX);
	compressor = new Compressor(JUMP);
	grid->setRefiner(refiner);
	grid->setCompressor(compressor);

	//const Real h_spaceconv = 1./FluidBlock2D::sizeX*pow(0.5, 4);
	//epsilon = (Real)MOLLFACTOR*sqrt(2.)*h_spaceconv;

	epsilon = (Real)MOLLFACTOR*sqrt(2.)*h_min;

	penalization = new I2D_PenalizationOperator(*grid, LAMBDA, Uinf, bRESTART);

	if (bPARTICLES)
		advection =new I2D_AdvectionOperator_Particles(*grid, CFL, LCFL);
	else
		advection =new I2D_AdvectionOperator(*grid, CFL);

	advection->set_Uinfinity(Uinf);


	if(sFMMSOLVER == "velocity")
		velsolver = new I2D_VelocitySolver_Mani(*grid, parser);
	else if(sFMMSOLVER == "velocity-wim")
		velsolver = new I2D_VelocitySolver_Wim(*grid, parser);
	else 
	{
#ifdef _I2D_MPI_
		if (sFMMSOLVER == "mpi-velocity")
			((I2D_VelocitySolverMPI_Mani *) velsolver)->set_grid(*grid);
#endif
		if(velsolver == NULL)
		{
			printf("VELOCITY SOLVER CANNOT BE NULL AT THIS POINT. aborting...");
			abort();
		}
	}

	diffusion = new I2D_DiffusionOperator_4thOrder(*grid, nu, FC);

	if (sOBSTACLE=="cyl")
	  {
	const bool isSharp = parser("-sharp").asBool();
	  obstacle = new I2D_CircularObstacleOperator(*grid,D/2,pos,epsilon, isSharp);
	  }
	else if (sOBSTACLE=="plate")
		obstacle = new I2D_RectangularObstacleOperator(*grid,1/90., pos, D, epsilon);
	else if (sOBSTACLE=="square")
		obstacle = new I2D_RectangularObstacleOperator(*grid,1., pos, D, epsilon);
	else if (sOBSTACLE == "NACA")
		obstacle = new I2D_WingObstacleOperator(*grid,parser,pos,epsilon);
	else if (sOBSTACLE == "links")
	{
		vector<Real> angles;
		angles.push_back(90.0);
		//angles.push_back(-45.0);
		//angles.push_back(20.0);
		//angles.push_back(-10.0);
		//angles.push_back(15.0);
		//angles.push_back(-25.0);
		//angles.push_back(30.0);
		//angles.push_back(15.0);
		Real length = D;
		Real width = length/5.0;
		obstacle = new I2D_LinkedBodies(*grid,length,width,pos,angles,epsilon);
	}
	else if (sOBSTACLE == "elly")
	{
		Real aspectRatio = parser("-elly_aspectRatio").asDouble(); 
		Real angle = parser("-elly_angle").asDouble();
		obstacle = new I2D_EllipticalObstacleOperator(*grid,D/2.0,aspectRatio,angle,pos,epsilon);
	}

	if(obstacle!=NULL)
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

I2D_FlowPastFixedObstacle::~I2D_FlowPastFixedObstacle()
{
	if(velsolver!=NULL){ delete velsolver; velsolver=NULL; }
	if(obstacle!=NULL){ delete obstacle; obstacle=NULL; }
	if(penalization!=NULL){ delete penalization; penalization=NULL; }
	if(advection!=NULL){ delete advection; advection=NULL; }
	if(diffusion!=NULL){ delete diffusion; diffusion=NULL; }

	if(grid!=NULL){ delete grid; grid=NULL; }
	if(refiner!=NULL){ delete refiner; refiner=NULL; }
	if(compressor!=NULL){ delete compressor; compressor=NULL; }

#ifdef _MRAG_GLUT_VIZ
	if(viewer!=NULL){ delete viewer; viewer=NULL; }
#endif	
}

void I2D_FlowPastFixedObstacle::_restart()
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
}

void I2D_FlowPastFixedObstacle::_save()
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
}

set<int> I2D_FlowPastFixedObstacle::_getBoundaryBlockIDs()
{
	set<int> result;
	vector<BlockInfo> vInfo = grid->getBlocksInfo();

	if (sRIGID_INLET_TYPE == "rigid_frame") //RIGID AT EACH BOUNDARY
		for(vector<BlockInfo>::iterator it = vInfo.begin(); it != vInfo.end(); it++)
		{
			const bool bX = it->index[0] == 0 || it->index[0] == pow(2.0, it->level)-1;
			const bool bY = it->index[1] == 0 || it->index[1] == pow(2.0, it->level)-1;

			if (bX || bY) result.insert(it->blockID);
		}
	else if (sRIGID_INLET_TYPE == "rigid_inlet_only") //RIGID ONLY AT THE INLET
		for(vector<BlockInfo>::iterator it = vInfo.begin(); it != vInfo.end(); it++)
		{
			const bool bX = it->index[0] == 0;

			if (bX) result.insert(it->blockID);
		}
	else if (sRIGID_INLET_TYPE == "free_outlet_only") //RIGID AT EACH BOUNDARY EXCEPT FOR THE OUTLET
		for(vector<BlockInfo>::iterator it = vInfo.begin(); it != vInfo.end(); it++)
		{
			const bool bX = it->index[0] == 0;
			const bool bY = it->index[1] == 0 || it->index[1] == pow(2.0, it->level)-1;

			if (bX || bY) result.insert(it->blockID);
		}
	else if (sRIGID_INLET_TYPE == "free_frame") //NON RIGID AT EACH BOUNDARY
		result.clear();
	else 
	{
		printf("I2D_FlowPastFixedObstacle::_getBoundaryBlockIDs: the chosen sRIGID_INLET_TYPE (=%s) is not supported. Aborting\n", sRIGID_INLET_TYPE.c_str());
		abort();
	}


	return result;
}

void I2D_FlowPastFixedObstacle::_dump(string filename)
{
	obstacle->characteristic_function();

	IO_VTKNative<W,B, 4,0> vtkdumper;
	vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
}

void I2D_FlowPastFixedObstacle::_dump()
{
	char buf[500];
	sprintf(buf, "avemaria_%07d", (int)step_id);
	_dump(buf);
}

void I2D_FlowPastFixedObstacle::_ic(Grid<W,B>& grid)
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();

	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		B& b = grid.getBlockCollection()[vInfo[i].blockID];

		for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
			{
				b(ix, iy).omega = 0;
				b(ix, iy).u[0] = 0;
				b(ix, iy).u[1] = 0;
				b(ix, iy).tmp = 0;
			}
	}

	obstacle->characteristic_function();
}

void I2D_FlowPastFixedObstacle::_refine(bool bUseIC)
{
	if (bUNIFORM) return;

	set<int> boundary_blocks = _getBoundaryBlockIDs();

	// For initial condition refine given Xs only
	if (bUseIC)
	{
		while(true)
		{
			((Refiner_BlackList*)refiner)->set_blacklist(&boundary_blocks);
			const int refinements = Science::AutomaticRefinement<0,0>(*grid, fwt_obstacle, RTOL, LMAX, 1, NULL, (void (*)(Grid<W,B>&))NULL, &boundary_blocks);
			_ic(*grid);
			if (refinements == 0) break;
		}
	}

	// For the rest refine given omega AND velocity
	if (!bREFINEOMEGAONLY)
	{
		while(true)
		{
			((Refiner_BlackList*)refiner)->set_blacklist(&boundary_blocks);
			const int refinements = Science::AutomaticRefinement<0,2>(*grid, fwt_wuv, RTOL, LMAX, 1, NULL, (void (*)(Grid<W,B>&))NULL, &boundary_blocks);
			if (refinements == 0) break;
		}
	}
	else
	{
		while(true)
		{
			((Refiner_BlackList*)refiner)->set_blacklist(&boundary_blocks);
			const int refinements = Science::AutomaticRefinement<0,0>(*grid, fwt_omega, RTOL, LMAX, 1, NULL, (void (*)(Grid<W,B>&))NULL, &boundary_blocks);
			if (refinements == 0) break;
		}
	}
}

void I2D_FlowPastFixedObstacle::_compress(bool bUseIC)
{
	if (bUNIFORM) return;

	set<int> boundary_blocks = _getBoundaryBlockIDs();

	// For initial condition compress given Xs only
	if (bUseIC)
	{
		obstacle->characteristic_function();
		Science::AutomaticCompression<0,0>(*grid, fwt_obstacle, CTOL, 1, NULL, (void (*)(Grid<W,B>&))NULL);
		return;
	}

	// For the rest refine given omega AND velocity
	if (!bREFINEOMEGAONLY)
		Science::AutomaticCompression<0,2>(*grid, fwt_wuv, CTOL, 1, NULL, (void (*)(Grid<W,B>&))NULL);
	else
		Science::AutomaticCompression<0,0>(*grid, fwt_omega, CTOL, 1, NULL, (void (*)(Grid<W,B>&))NULL);
}

double I2D_FlowPastFixedObstacle::_initial_dt(int nsteps)
{
	if ( step_id>=nsteps ) return 1000.0;

	const double dt_min = 1e-7;//1e-6;
	const double dt_max = 1e-4;//1e-3;
	return dt_min + (double)(step_id)/(double)(nsteps)*(dt_max - dt_min);
	//return 1e-3; //For convergence the best is ramp = 1 step with dt 1e-3
}

void I2D_FlowPastFixedObstacle::_tnext(double &tnext, double &tnext_dump, double &tend)
{
	//0. setup
	//1. prepare candidate tnexts
	//2. pick the smallest one, report who wins

	//0.
	vector<double> tnext_candidates;
	vector<string> tnext_names;

	const double moduinf = sqrt(Uinf[0]*Uinf[0]+Uinf[1]*Uinf[1]);
	const double nondim_factor = (moduinf==0.0)?1.0:(moduinf*2.0/D);
	const double max_dx = (1./B::sizeX)*pow(0.5,grid->getCurrentMinLevel());
	const double min_dx = (1./B::sizeX)*pow(0.5,grid->getCurrentMaxLevel());
	const double max_vel = advection->compute_maxvel();

	//1.
	tend = TEND/nondim_factor;

	tnext_candidates.push_back(t + _initial_dt(RAMP));
	tnext_names.push_back("RAMP");

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
		const double delta = 1./(DUMPFREQ*nondim_factor);
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
		fprintf(f, "stepid=%d\tT=%e\tDT=%e\t", (int)step_id, tnext*nondim_factor, (tnext - t)*nondim_factor);
		for(int i=0; i<(int)tnext_candidates.size(); i++)
			fprintf(f, "%s: %2.2e,\t", tnext_names[i].c_str(), tnext_candidates[i] - t);

		fprintf(f, "GTS-Fc: %2.2e,\t", FC*pow(min_dx,2)/(8.0*nu));

		fprintf(f, "\n");
		fclose(f);
	}
}

void I2D_FlowPastFixedObstacle::run()
{
	if( Uinf[0]==0.0 && Uinf[1]==0.0 )
	{
		printf("I2D_FlowPastFixedObstacle::run(): Uinf = 0!\n");
		abort();
	}

	if(obstacle == NULL)
	{
		printf("I2D_FlowPastFixedObstacle::run(): obstacle = NULL!\n");
		abort();
	}

	printf("\n\nI AM NEW!\n\n");

	const tbb::tick_count start_instant = tbb::tick_count::now();

	vector<Real> infoObstacle;

	Real cor[2] = {0.0,0.0};

	const double moduinf = sqrt(Uinf[0]*Uinf[0]+Uinf[1]*Uinf[1]);
	const double nondim_factor = (moduinf==0.0)?1.0:(moduinf*2.0/D);

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

			profiler.push_start("PEN");
			obstacle->characteristic_function();
			penalization->perform_timestep(dt);
			obstacle->characteristic_function();
			obstacle->getObstacleInfo(infoObstacle);
			profiler.pop_stop();
			printf("DONE WITH PENALIZATION\n");

			if(infoObstacle.size()!=0){ cor[0] = infoObstacle[0]; cor[1] = infoObstacle[1]; }
			penalization->compute_dragandstuff(t,D,cor,"diag");
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

			if(step_id%SAVEFREQ==0)
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
