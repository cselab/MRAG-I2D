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

class I2D_ImposedCylinder: public I2D_FloatingObstacleOperator
{	
public:
	
	class Cylinder  
	{
	private:
		Real _mollified_heaviside(const double dist, const double eps) const;
		
	public:		
		double angle, xm, ym, D, angular_velocity, vx, vy, m, J, rho;		
		
		Cylinder();
		Cylinder(Real xm, Real ym, Real D, Real angle_rad);		
		
		//const Cylinder& operator=(const Cylinder& c);
		void update_all(double dt, double t);
		void restart(FILE * f);
		void save(FILE * f) const;
		Real sample(const Real x_, const Real y_, const Real eps) const;
		void bbox(const Real eps, Real xmin[2], Real xmax[2]) const;
		void restart();
		void save();
	};
	
	Real vx_imposed, vy_imposed;
	Cylinder * shape;
	
	I2D_ImposedCylinder(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real _vx, const Real _vy, const Real _D, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization);
	~I2D_ImposedCylinder();
	
	void characteristic_function();		
	Real getD() const {return D;}
	
	void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	void computeDesiredVelocity(const double t);
	void create(const double t){}

	void save(const double t, string filename = std::string());
	void restart(const double t, string filename = std::string());
	void refresh(const double t, string filename = std::string());
	
protected:
	virtual void _setMotionPattern(const Real t);
};
