/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"
#include "I2D_FloatingObstacleOperator.h"

class I2D_CarlingFish: public I2D_FloatingObstacleOperator
{
protected:
	double xm_corr, ym_corr, vx_corr, vy_corr, omega_corr;

public:
	
	class Fish  
	{
	protected:
		int minStart, maxEnd;
		
		double O2W[3][3], W2O[3][3];
		double fish_xmin[2], fish_xmax[2];
		
		double * dataX;
		double * dataY;
		double * dataVX;
		double * dataVY;
		double * dataDist;
		double * X;
		double * Y;
		double * S;
		double * W;
		double * NORX;
		double * NORY;
		double * VX;
		double * VY;
		double * VNORX;
		double * VNORY;
		double * SHAPEX;
		double * SHAPEY;
		double * CHI, *CHI2;
		double * VDEFX;
		double * VDEFY;
		double * dataSDF;
		
		template<typename R> void _w2o(const R xw[2], R xo[2]) const;
		template<typename R> void _o2w(const R xo[2], R xw[2]) const;
		int _dd2l(const int & ix, const int & iy) const;
		int _ud2l(const int & ix, const int & iy) const;
		double _bsp4(double x) const;
		double _bsp2(double x) const;
		void _check_ow_mapping() const;
		void _set_ow_mapping();
		double _sample(const double x, const double y, const double h) const;
		void _sample(const double x, const double y, double & Xs, double & defvelx, double & defvely, const double h) const;
                double _sample(const double x, const double y) const;
                void _sample(const double x, const double y, double & Xs, double & defvelx, double & defvely) const;
		virtual double _getWidth(const double & ss) const;
		void _calculateCenterlineHead(const double t, const double factor, const double T, const double ds, const double * ss, double * rX, double * rY, const int nn);
		void _calculateNormalsAndVelocities(const double t, const double T, const double ds, const double * ss, const double * rX, const double * rY, double * vX, double * vY, double * norX, double * norY, double * vNorX, double * vNorY, const int nn);
		void _createShapeBoundary( const double * rX, const double * rY, const double * width, const double * norX, const double * norY, double * shapeX, double * shapeY, const int n);
		void _getDefGridCoord(const int & ix, const int & iy, const double & dg, const int & coeff, const double * rX, const double * rY, const double * norX, const double * norY, double & xCoord, double & yCoord) const;
		void _getDefGridVel(const int & ix, const int & iy, const double & dg, const int & coeff, const double * vX, const double * vY, const double * vNorX, const double * vNorY, double & vxCoord, double & vyCoord) const;
		void _getDefGridCharFunc(const int & ix, const int & iy, const int & spany, const double & dg, const int & coeff, const double * rX, const double * rY, const double * width, double & sdf) const;
		void _fillDefGrid();
		void _rigidTranslation(const Real x, const Real y );
		void _getCenterOfMassFull(double & xCoord, double & yCoord) const;
		void _getVelCenterOfMassFull(double & vxCoord, double & vyCoord) const;
		void _getAngularMomentumFull(double & L) const;
		void _getScalarMomentOfInertiaFull(double & II) const;
		void _centerlineCenterOfMassFrameTransform(const double & cmX, const double & cmY, const double & vcmX, const double & vcmY);
		void _defGridCenterOfMassFrameTransform(const double & cmX, const double & cmY, const double & vcmX, const double & vcmY);
		void _getRotationalVelocityAboutTheOrigin(const double & omega, const double & rX, const double & rY, double & vRotX, double & vRotY) const;
		void _correctCenterlineForRotationalImpulse(const double & omega);
		void _correctDefGridForRotationalImpulse(const double & omega);
		void _rotateAboutTheOrigin(const double & theta, double * x, double * y, const int n);
		void _rotateDefGridAboutTheOrigin(const double & theta);
		void _bbox();
		void _bboxShapeBoundary(const double * shapeX, const double * shapeY, const int n, double & width, double & height);
		void _computeMinStartMaxEnd(const double * width, int & _minStart, int & _maxEnd);
		void _getUniformGridProperties();
		void _rotateAboutTheOriginSingle(const double & theta, double & x, double & y) const;
		void _cubicInterpolation(double t0, double t1, double t, double k0, double k1, double & k, double & dkdt);

		void _cross(const double * v1, const double * v2, double * v3) const;
		double _dot(const double * v1, const double * v2) const;
		void _boundingBoxTriangle(const double * a, const double * b, const double * c, const double H, int & ixMax, int & ixMin, int & iyMax, int & iyMin) const;
		bool _pointInTriangle(const double * p, const double * a, const double * b, const double * c) const;
		bool _bilinearTriangle(const double * p, const double * a, const double * b, const double * c, const double * va, const double * vb, const double * vc, double * interp) const;
		
						
	public:		
		const double TOrig;
		int SIZEX, SIZEY, MAPSIZEX, MAPSIZEY, N;
		double WH, SB, WT, ST, T, DS, EXTENSION, EPS;
		double angle, angleInSpace, xm, ym, D, angular_velocity, angular_velocityInSpace, vx, vy, m, J, vdefx, vdefy, phase, tau, traslX, traslY;
		int LMAX;
		bool SHARP;

		Fish(double xm, double ym, double _D, double _T, double phase, double angle_rad, double angleInSpace_rad, double eps, const int LMAX, const bool isSharp=false);
		virtual ~Fish();
		
		const double * get_dataXptr() { return dataX; }
		const double * get_dataYptr() { return dataY; } 		
		const double * get_dataVXptr() { return dataVX; } 
		const double * get_dataVYptr() { return dataVY; }
		const double * get_dataDistptr() { return dataDist; }		
		void update_all(double dt);
		virtual void restart(FILE * f);
		virtual void save(FILE * f) const;
		double sample(const double x_, const double y_, const double h_) const;
		void sample(const double x_, const double y_, double & Xs, double & defvelx, double & defvely, const double h_) const;
                double sample(const double x_, const double y_) const;
                void sample(const double x_, const double y_, double & Xs, double & defvelx, double & defvely) const;
		void bbox(const double eps, double xmin[2], double xmax[2]) const;
		virtual void updateInSpace(const double t);
		void bilinearInterpolation();		
		static double mollified_heaviside(const double dist, const double eps);
		void clearUniformGrids();
		virtual void clear(){};
	};
	
	Fish * shape;

	I2D_CarlingFish(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real T, const Real phase, const Real angle, const Real angleInSpace, const Real eps, const Real Uinf[2],
			I2D_PenalizationOperator& penalization, const int LMAX, const int ID = 0, RL::RL_TabularPolicy ** policy = NULL, const int seed = 0);
	virtual ~I2D_CarlingFish();
	
	void characteristic_function();		
	Real getD() const {return D;}
	
	void create(const double t);
	virtual void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	void computeDesiredVelocity(const double t);

	virtual void save(const double t, string filename = std::string());
	virtual void restart(const double t, string filename = std::string());
	virtual void refresh(const double t, string filename = std::string());
	Real getModulusMaxVel(){ return shape->D/shape->T; }
	virtual void getInfo(vector< vector<double> > & output) const;
};
