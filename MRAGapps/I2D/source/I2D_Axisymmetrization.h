/*
 *  I2D_FlowPastFixedObstacle.h
 *  IncompressibleFluids2D
 *
 *  Created by Roman Schaerer on 30/03/11.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 *
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"
#include "I2D_VelocityOperator.h"
#include "I2D_PenalizationOperator.h"
#include "I2D_AdvectionOperator.h"
#include "I2D_DiffusionOperator.h"

class InitialCondition;

class I2D_Axisymmetrization : public I2D_Test
{
protected:
	//"constants" of the sim
	int BPD, JUMP, LMAX, ADAPTFREQ, SAVEFREQ, RAMP;
	Real DUMPFREQ, RE, CFL, LCFL, RTOL, CTOL, LAMBDA, TEND, Uinf[2], LAMBDADT, XPOS, YPOS, epsilon;
	bool bPARTICLES, bUNIFORM, bCORRECTION, bRESTART, bREFINEOMEGAONLY, bFMMSKIP;
	string sFMMSOLVER, sIC, sRIGID_INLET_TYPE;

	//state of the sim
	double t;
	int step_id;
	Real nondim_factor_time;

	ArgumentParser parser;

	Grid<W,B> * grid;

	Refiner * refiner;
	Compressor * compressor;

	BlockFWT<W, B, vorticity_projector, false, 1> fwt_omega;
	BlockFWT<W, B, velocity_projector, false, 2> fwt_velocity;
	BlockFWT<W, B, vorticityANDvelocityANDchi_projector, false, 4> fwt_wuvx;

	set<int> _getBoundaryBlockIDs();

	void _dump(string filename);
	void _dump();

	void _ic(Grid<W,B>& grid);
	void _restart();
	void _save();

	Real _initial_dt(int nsteps);
	void _tnext(Real &tnext, Real& tnext_dump, Real& tend);

	void _refine(bool bUseIC);
	void _compress(bool bUseIC);

	I2D_VelocityOperator * velsolver;
	I2D_AdvectionOperator * advection;

	InitialCondition* initialCondition;
	Profiler profiler;

#ifdef _MRAG_GLUT_VIZ
	GridViewer * viewer;
#endif

public:

	I2D_Axisymmetrization (const int argc, const char ** argv);
	~I2D_Axisymmetrization (){}

	void run();
	void paint(){}
};

