/*
 *  I2D_FTLE.h
 *  IncompressibleFluids2D
 *
 *  Created by Mattia Gazzola on 20/01/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
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

class I2D_FTLE: public I2D_Test
{
protected:
	//"constants" of the sim
	string PATH;
	Real RTOL, DTFTLE, TFTLE, TSTARTFTLE, TENDFTLE, Uinf[2];
	int JUMP, LMAX;
	bool bUNIFORM, bRESTART;

	//state of the sim
	Real t_completed;

	ArgumentParser parser;

	Profiler profiler;

	Grid<W,B> * grid;

	BlockProcessing block_processing;

	Refiner * refiner;
	Compressor * compressor;

	BlockFWT<W, B, vorticity_projector, false, 1> fwt_omega;
	BlockFWT<W, B, velocity_projector, false, 2> fwt_velocity;
	BlockFWT<W, B, obstacle_projector, false, 1> fwt_obstacle;
	BlockFWT<W, B, vorticityANDvelocityANDchi_projector, false, 4> fwt_wuvx;
	BlockFWT<W, B, vorticityANDvelocity_projector, false, 3> fwt_wuv;

	void _getdirContent(string dir, vector<string> &files);
	void _getdirRestarts(vector<string> &restarts);
	void _mapping(map< string, double > &mapping);
	void _extractSubsed(double tStart, double tEnd, const map< string, double > &mapping, map< string, double > &subset);
	double _loadStatus(string restart);
	void _loadGrid(string restart);
	void _dump(string filename);
	void _computeFwdFTLE(double tStart, double tEnd, int idx, map< string, double > & mapping);
	void _save(string filename);
	void _save();
	void _restart();
	void _refine();
	double _advectParticles(map<long int, vector<Real> > & particles, vector<Real> & dts, map< string, double > &subset);
	void _FTLE(string initialField, double T, map<long int, vector<Real> > & particles, vector<Real> & dts, map< I3, unsigned int > & code);
	void _extractSubset(double tStart, double tEnd, const map< string, double > &mapping, map< string, double > &subset);
	void _createParticleSet(string restart, map<long int, vector<Real> > & particles, map< I3, unsigned int > & code);
	void _createDtsSet(const map< string, double > & subset, vector<Real> & dts);

public:

	I2D_FTLE(const int argc, const char ** argv);
	~I2D_FTLE();

	void run();
	void paint(){}
};
