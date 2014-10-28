/*
 *  I2D_FlowPastFixedObstacle.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 *  
 *   
 *	PATER NOSTER, qui es in caelis, sanctificetur nomen tuum. 
 *	Adveniat regnum tuum. 
 *	Fiat voluntas tua, sicut in caelo et in terra. 
 *	Panem nostrum quotidianum da nobis hodie, 
 *	et dimitte nobis debita nostra sicut et nos dimittimus debitoribus nostris. 
 *	Et ne nos inducas in tentationem, 
 *	sed libera nos a malo. 
 *	Amen.
 *  
 *  
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

class I2D_FlowPastFixedObstacle: public I2D_Test
{
protected:
	//"constants" of the sim
	int BPD, JUMP, LMAX, ADAPTFREQ, SAVEFREQ, RAMP, MOLLFACTOR;
	Real DUMPFREQ, RE, CFL, LCFL, RTOL, CTOL, LAMBDA, D, TEND, Uinf[2], nu, LAMBDADT, XPOS, YPOS, epsilon, FC;
	bool bPARTICLES, bUNIFORM, bCORRECTION, bRESTART, bREFINEOMEGAONLY, bFMMSKIP;
	string sFMMSOLVER, sOBSTACLE, sRIGID_INLET_TYPE;
	
	//state of the sim
	double t;
	long unsigned int step_id;
	
	ArgumentParser parser;
	
	Grid<W,B> * grid;
	
	Refiner * refiner;
	Compressor * compressor;
	
	BlockFWT<W, B, vorticity_projector, false, 1> fwt_omega;
	BlockFWT<W, B, velocity_projector, false, 2> fwt_velocity;
	BlockFWT<W, B, obstacle_projector, false, 1> fwt_obstacle;
	BlockFWT<W, B, vorticityANDvelocityANDchi_projector, false, 4> fwt_wuvx;
	BlockFWT<W, B, vorticityANDvelocity_projector, false, 3> fwt_wuv;
	
	set<int> _getBoundaryBlockIDs();
	
	void _dump(string filename);
	void _dump();
	
	void _ic(Grid<W,B>& grid);
	void _restart();
	void _save();
	
	virtual Real _initial_dt(int nsteps);
	virtual void _tnext(double &tnext, double& tnext_dump, double& tend);
	
	void _refine(bool bUseIC);
	void _compress(bool bUseIC);
	
	I2D_VelocityOperator * velsolver;
	I2D_ObstacleOperator * obstacle;
	I2D_PenalizationOperator * penalization;
	I2D_AdvectionOperator * advection;
	I2D_DiffusionOperator * diffusion;
	
	Profiler profiler;
	
#ifdef _MRAG_GLUT_VIZ
	GridViewer * viewer;
#endif
	
public:
	
	I2D_FlowPastFixedObstacle(const int argc, const char ** argv);
	~I2D_FlowPastFixedObstacle();
	
	virtual void run();
	void paint(){}
};
