/*
 *  I2D_Axisymmetrization.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Roman Schaerer on 30/03/11.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 *
 */

#include "I2D_Axisymmetrization.h"
#include "I2D_AdvectionOperator_Particles.h"
#include "I2D_VelocitySolver_Mani.h"
#include <cmath>

#ifdef _I2D_MPI_
#include "I2D_VelocitySolverMPI_Mani.h"
#endif

#include <sys/stat.h>

static const int maxParticleStencil[2][3] = {
	-3, -3, 0,
	+4, +4, +1
};

class Diagnostics
{
	Real m_J0;

public:
	Diagnostics () {
		m_J0 = -1;
	}

	void compute (Grid<W,B>& grid, Real t)
    {
        std::vector<BlockInfo> vInfo = grid.getBlocksInfo();

		Real J20 = 0, J02 = 0, J11 = 0;
		Real Gamma = 0;
		Real max_omega = 0;
        Real enstrophy = 0;
        
		for(int i=0; i<vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			B& b = grid.getBlockCollection()[vInfo[i].blockID];
			Real dxdy = info.h[0] * info.h[1];
			for(int iy=0; iy<B::sizeY; iy++)
            {
				for(int ix=0; ix<B::sizeX; ix++)
                {
					Real pos[2];
					info.pos(pos,ix,iy);
					J20 += b(ix,iy).omega * (pos[0]-0.5) * (pos[0]-0.5) * dxdy;
					J02 += b(ix,iy).omega * (pos[1]-0.5) * (pos[1]-0.5) * dxdy;
					J11 += b(ix,iy).omega * (pos[0]-0.5) * (pos[1]-0.5) * dxdy;
					Gamma += b(ix,iy).omega * dxdy;
					max_omega = std::max(max_omega, b(ix,iy).omega);
                    enstrophy += b(ix,iy).omega * b(ix,iy).omega * dxdy;
				}
			}
		}

		const Real J = J02 + J20;

		if (m_J0 == -1 && J > 0) {
			m_J0 = J;
			return;
		}

		const Real D = J20 - J02;
		const Real R = sqrt (D*D+4*J11*J11);

		const Real nu_eff = m_J0 < 0 ? 0.0 : std::abs((J-m_J0)/(4*Gamma*t));

		const Real lambda_eff = sqrt((J+R)/(J-R));
		const Real Re_eff = std::abs(Gamma/nu_eff);

		fstream file;

		file.open ("diagnostics.txt",  std::fstream::app | std::fstream::out);
		if (file.fail ()) {
			std::cerr << "\nFailed to open diagnostics.txt file. Aborting now.\n";
			abort ();
		}
		file << setprecision (6) << scientific << t << "\t" << max_omega << "\t" << Gamma << "\t" \
				<< nu_eff << "\t" << Re_eff << "\t" << lambda_eff << "\t" << enstrophy << "\n";

		file.close ();
		if (file.fail ()) {
			std::cerr << "\nFailed to close measurements.txt file. Aborting now.\n";
			abort ();
		}
	}
};

class InitialCondition
{
protected:
    const Real scaleFactL;
    const Real scaleFactT;
public:
    InitialCondition():scaleFactL(8.0), scaleFactT(1.0){}
	virtual Real getOmega (Real, Real) = 0;
	virtual Real getMaxOmega () = 0;
};

class EllipticalVortexCase1 : public InitialCondition
{
	Real m_q, m_lambda, m_r0, m_AR;
	Real m_cx, m_cy;
    
public:
	EllipticalVortexCase1 ()
    {
        // values from the pk1997 paper
		m_q = 2.56085 ; // dimensionless
		m_lambda = 20.0 * scaleFactT; // lambda -> 1/T
		m_r0 = 0.8 / scaleFactL; // r0 -> L
        m_AR = 2.0 ; // aspect ratio
		m_cx = m_cy = 0.5; // already in computational coordinates
	}

	Real getOmega (Real x, Real y)
    {
		const Real r = std::sqrt (m_AR*m_AR*(x-m_cx)*(x-m_cx)+(y-m_cy)*(y-m_cy));
        const Real z = r/m_r0;
        return (z <= 1.0  ? m_lambda * (1.0 - std::exp (-(m_q/z) * std::exp(1./(z-1)) ) ) : 0.0 );
	}

	Real getMaxOmega ()
    {
		return m_lambda * scaleFactT;
	}
};

class EllipticalVortexCase2 : public InitialCondition
{
    Real m_lambda, m_r0, m_AR;
	Real m_cx, m_cy;

public:

	EllipticalVortexCase2 ()
    {
		m_lambda = 20.0 * scaleFactT ; // lambda -> 1/T
		m_r0 = 0.8 / scaleFactL; // r0 -> L
        m_AR = 2.0 ; // aspect ratio
		m_cx = m_cy = 0.5; // already in computational coordinates
	}

	Real getOmega (Real x,Real y)
    {
		const Real r = std::sqrt(m_AR*m_AR*(x-m_cx)*(x-m_cx)+(y-m_cy)*(y-m_cy));
        const Real z = r/m_r0;
        return (z <= 1.0  ? m_lambda * (1.-std::pow (z, 4.0) ) : 0.0 );
	}

	Real getMaxOmega ()
    {
		return m_lambda * scaleFactT;
	}
};


I2D_Axisymmetrization::I2D_Axisymmetrization (const int argc, const char ** argv):
I2D_Test(), parser(argc, argv), t(0), step_id(0),
velsolver(NULL), advection(NULL), initialCondition(NULL)
{
	printf("////////////////////////////////////////////////////////////\n");
	printf("////////////         Axisymmetrization       ///////////////\n");
	printf("////////////////////////////////////////////////////////////\n");

	bRESTART = parser("-restart").asBool();
	sRIGID_INLET_TYPE = parser("-rio").asString();
	bREFINEOMEGAONLY = parser("-refine-omega-only").asBool();
	bFMMSKIP = parser("-fmm-skip").asBool();
	sIC = parser ("-ic").asString();

	parser.set_strict_mode();

	BPD = parser("-bpd").asInt();
	JUMP = parser("-jump").asInt();
	LMAX = parser("-lmax").asInt();
	TEND = parser("-tend").asDouble();
	ADAPTFREQ = parser("-adaptfreq").asInt();
	DUMPFREQ = parser("-dumpfreq").asDouble();
	SAVEFREQ = parser("-savefreq").asInt();
	CFL = parser("-cfl").asDouble();
	LCFL = parser("-lcfl").asDouble();
	RTOL = parser("-rtol").asDouble();
	CTOL = parser("-ctol").asDouble();
	bPARTICLES = parser("-particles").asBool();
	bUNIFORM = parser("-uniform").asBool();
	sFMMSOLVER = parser("-fmm").asString();

	parser.save_options();

	assert(TEND >= 0.0);
	assert(BPD > 1);
	assert(ADAPTFREQ > 0);
	assert(DUMPFREQ > 0);
	assert(JUMP >= 1);
	assert(LMAX >= 0);
	assert(CFL > 0 && CFL<1);
	assert(LCFL > 0 && LCFL<1);
	assert(RTOL > 0);
	assert(CTOL > 0);
	assert(sFMMSOLVER != "");
	assert(sIC != "");

#ifdef _I2D_MPI_
	if (sFMMSOLVER == "mpi-velocity")
		velsolver = new I2D_VelocitySolverMPI_Mani(argc, argv);
#endif

	const Real h_min = 1./FluidBlock2D::sizeX*pow(0.5, LMAX);

	if (sIC == "case1") {
		initialCondition = new EllipticalVortexCase1;
	} else if (sIC == "case2") {
		initialCondition = new EllipticalVortexCase2;
	} else {
		std::cout << "\nNo initial vorticity field specified. Aborting...\n";
		abort ();
	}

	nondim_factor_time = 4*M_PI/(initialCondition->getMaxOmega());

	if (bPARTICLES)
		grid = new Grid<W,B>(BPD,BPD,1, maxParticleStencil);
	else
		grid = new Grid<W,B>(BPD,BPD,1);

	assert(grid != NULL);

	refiner = new Refiner_BlackList(JUMP, LMAX);
	compressor = new Compressor(JUMP);
	grid->setRefiner(refiner);
	grid->setCompressor(compressor);

	if (bPARTICLES)
		advection =new I2D_AdvectionOperator_Particles(*grid, CFL, LCFL);
	else
		advection =new I2D_AdvectionOperator(*grid, CFL);

	Real Uinf[2];
	Uinf [0]= Uinf[1] = 0;
	advection->set_Uinfinity(Uinf);


	if(sFMMSOLVER == "velocity")
		velsolver = new I2D_VelocitySolver_Mani(*grid, parser);
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

void I2D_Axisymmetrization::_restart()
{
	//read status
	{
		FILE * f = fopen("restart.status", "r");
		assert(f != NULL);
		float val = -1;
		fscanf(f, "time: %e\n", &val);
		assert(val>=0);
		t=val;
		step_id = -1;
		fscanf(f, "stepid: %d\n", &step_id);
		assert(step_id >= 0);
		fclose(f);
	}

	printf("DESERIALIZATION: time is %f and step id is %d\n", t, step_id);

	//read grid
	IO_Binary<W,B> serializer;
	serializer.Read(*grid, "restart");
}

void I2D_Axisymmetrization::_save()
{
	printf("****SERIALIZING****\n");

	//write status
	{
		FILE * f = fopen("restart.status", "w");
		if (f != NULL)
		{
			fprintf(f, "time: %20.20e\n", t);
			fprintf(f, "stepid: %d\n", step_id);
			fclose(f);
		}

		printf( "time: %20.20e\n", t);
		printf( "stepid: %d\n", step_id);
	}

	string numbered_filename;

	{
		char buf[500];
		sprintf(buf, "restart_%07d", step_id);
		numbered_filename = string(buf);
	}

	//write grid
	IO_Binary<W,B> serializer;

	serializer.Write(*grid, "restart");
	serializer.Write(*grid, numbered_filename.c_str());

	printf("****SERIALIZING DONE****\n");

	{
		FILE * f = fopen("history.txt", step_id == 0? "w" : "a");
		if (f!= NULL)
		{
			fprintf(f, "%10.10f %d\n", t, step_id);
			fclose(f);
		}
	}
}

set<int> I2D_Axisymmetrization::_getBoundaryBlockIDs()
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
		printf("I2D_Axisymmetrization::_getBoundaryBlockIDs: the chosen sRIGID_INLET_TYPE (=%s) is not supported. Aborting\n", sRIGID_INLET_TYPE.c_str());
		abort();
	}
    
    
	return result;
}

void I2D_Axisymmetrization::_dump(string filename)
{
	IO_VTKNative<W,B, 4,0> vtkdumper;
	vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
}

void I2D_Axisymmetrization::_dump()
{
	char buf[500];
	sprintf(buf, "avemaria_%07d", step_id);
	_dump(buf);
}



void I2D_Axisymmetrization::_ic(Grid<W,B>& grid)
{
	assert (initialCondition != NULL);

	vector<BlockInfo> vInfo = grid.getBlocksInfo();

	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		B& b = grid.getBlockCollection()[vInfo[i].blockID];
		for(int iy=0; iy<B::sizeY; iy++)
			for(int ix=0; ix<B::sizeX; ix++)
			{
				Real point[2];
				info.pos(point,ix,iy);
				b(ix,iy).omega = initialCondition->getOmega (point[0],point[1]);
				b(ix,iy).u[0] = b(ix,iy).u[1] = b(ix,iy).tmp = 0.0;
			}
	}
}

void I2D_Axisymmetrization::_refine(bool bUseIC)
{
	if (bUNIFORM) return;

	set<int> boundary_blocks = _getBoundaryBlockIDs();

	if (bUseIC)
	{
		//fill with Xs
		while(true)
		{
			((Refiner_BlackList*)refiner)->set_blacklist(&boundary_blocks);

			const int refinements = Science::AutomaticRefinement<0,0>(*grid, fwt_omega, RTOL, LMAX, 1, NULL, (void (*)(Grid<W,B>&))NULL, &boundary_blocks);

			_ic(*grid);

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

	//velocity
	if (!bREFINEOMEGAONLY)
		while(true)
		{
			((Refiner_BlackList*)refiner)->set_blacklist(&boundary_blocks);

			const int refinements = Science::AutomaticRefinement<0,1>(*grid, fwt_velocity, RTOL, LMAX, 1, NULL, (void (*)(Grid<W,B>&))NULL, &boundary_blocks);

     		if (refinements == 0) break;
		}
}

void I2D_Axisymmetrization::_compress(bool bUseIC)
{
	if (bUNIFORM) return;

	set<int> boundary_blocks = _getBoundaryBlockIDs();

	//omega AND velocity
	if (!bREFINEOMEGAONLY)
		Science::AutomaticCompression<0,3>(*grid, fwt_wuvx, CTOL, 1, NULL, (void (*)(Grid<W,B>&))NULL);
	else
		Science::AutomaticCompression<0,0>(*grid, fwt_omega, CTOL, 1, NULL, (void (*)(Grid<W,B>&))NULL);

}

Real I2D_Axisymmetrization::_initial_dt(int nsteps)
{
	if (step_id>nsteps) return 1000.0;
	const Real dt_min = 1e-6;
	const Real dt_max = 1e-2;
	return dt_min + (Real)(step_id)/(Real)(nsteps)*(dt_max - dt_min);
}

void I2D_Axisymmetrization::_tnext(Real &tnext, Real& tnext_dump, Real& tend)
{
	//0. setup
	//1. prepare candidate tnexts
	//2. pick the smallest one, report who wins

	//0.
	vector<Real> tnext_candidates;
	vector<string> tnext_names;

	const Real max_dx = (1./B::sizeX)*pow(0.5,grid->getCurrentMinLevel());
	const Real min_dx = (1./B::sizeX)*pow(0.5,grid->getCurrentMaxLevel());
	const Real max_vel = advection->compute_maxvel();

	//1.
	tend = TEND/nondim_factor_time;

	tnext_candidates.push_back(t + max_dx/max_vel*CFL*0.999);
	tnext_names.push_back("CFL");

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

	{
	  //		const Real delta = 1./(DUMPFREQ*nondim_factor_time);
	  const Real delta = nondim_factor_time / DUMPFREQ ; // wim: changed this to DUMPFREQ dumps per non-dim time. dunno what was the shit above.. wtf mattia
		tnext_dump = (1.+floor(t/delta))*delta;
	}

	//2.
	{
		tnext = tnext_candidates[0];
		int imin = 0;

		for(int i=1; i<tnext_candidates.size(); i++)
			if (tnext_candidates[i] < tnext)
			{
				imin = i;
				tnext = tnext_candidates[i];
			}

		const string dtBound = tnext_names[imin] + " bound";

		printf("####################### t is %e, dt is %e, %s,\t{", t, tnext - t, dtBound.c_str());

		for(int i=0; i<tnext_candidates.size(); i++)
			printf("%6s:%2.2e ", tnext_names[i].c_str(), tnext_candidates[i] - t);
		printf("}\n");

		FILE * f = fopen("report.txt", step_id == 0 ? "w" : "a");
		assert(f!=NULL);

		fprintf(f, "####################### t is %e, dt is %e and is bound by: %s", t, tnext - t, dtBound.c_str());
		fprintf(f, "stepid=%d\tT=%e\tDT=%e\t", step_id, tnext*nondim_factor_time, (tnext - t)*nondim_factor_time);
		for(int i=0; i<tnext_candidates.size(); i++)
			fprintf(f, "%s: %2.2e,\t", tnext_names[i].c_str(), tnext_candidates[i] - t);

		fprintf(f, "\n");
		fclose(f);
	}
}

void I2D_Axisymmetrization::run()
{
	const tbb::tick_count start_instant = tbb::tick_count::now();

	vector<Real> infoObstacle;

	Real cor[2] = {0.0,0.0};

	Diagnostics diagnostics;

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

			Real tnext, tnext_dump, tend;
			_tnext(tnext, tnext_dump, tend);
			const Real dt = tnext - t;
			printf("DONE WITH TNEXT\n");

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
				printf("T=TEND=%2.2f reached! (t=%2.2f)\n", TEND, t);
				exit(0);
			}

			{
				const tbb::tick_count now = tbb::tick_count::now();
				const Real Tcurr = t*nondim_factor_time;
				const Real wallclock = (now-start_instant).seconds();
				const int nblocks = grid->getBlocksInfo().size();

				FILE * f = fopen("more-perfmon.txt", step_id == 1 ? "w" : "a");

				if (f!=NULL)
				{
					if (step_id == 0)
						fprintf(f,"stepid\twall-clock[s]\tT[-]\tblocks\tADV-rhs\tDIFF-rhs\n");

					fprintf(f,"%d\t%e\t%e\t%d\t%d\n", step_id, wallclock, Tcurr, nblocks, advection->get_nofrhs());

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

		diagnostics.compute (*grid,t);


		printf("COMPRESS..\n");
		profiler.push_start("COMPRESS");
		_compress(false);
		profiler.pop_stop();
		printf("DONE WITH COMPRESS\n");

	}
}
