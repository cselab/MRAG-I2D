/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#pragma once

#include <iostream>

class I2D_Frenet2D
{
public:
	// Constructors
	I2D_Frenet2D(){};
	~I2D_Frenet2D(){};
	
	// Methods
	void solve(const double ds, const double * k, double * rX, double * rY, double * tanX, double * tanY, double * norX, double * norY, const unsigned int n)
	{
		double epsilonX = 1.0;
		double epsilonY = 0.0;
		double nuX = 0.0;
		double nuY = 1.0;
		
		double rX_old = 0.0;
		double rY_old = 0.0;
		double eX_old = epsilonX;
		double eY_old = epsilonY;
		double nuX_old = nuX;
		double nuY_old = nuY;
		
		double eX_new = 0.0;
		double eY_new = 0.0;
		double nuX_new = 0.0;
		double nuY_new = 0.0;
		double rX_new = 0.0;
		double rY_new = 0.0;
		
		double d_eX_ds = 0.0;
		double d_eY_ds = 0.0;
		double d_nuX_ds = 0.0;
		double d_nuY_ds = 0.0;
		double d_rX_ds = 0.0;
		double d_rY_ds = 0.0;
		
		double d1 = 0.0;
		double d2 = 0.0;
		
		for(unsigned int i=0; i<n; i++)
		{
			// Store coordinates
			rX[i] = rX_old;
			rY[i] = rY_old;
			tanX[i] = eX_old;
			tanY[i] = eY_old;
			norX[i] = nuX_old;
			norY[i] = nuY_old;
			
			// Calcualte derivatives
			d_eX_ds = k[i]*nuX_old;
			d_eY_ds = k[i]*nuY_old;
			d_nuX_ds = -k[i]*eX_old;
			d_nuY_ds = -k[i]*eY_old;
			d_rX_ds = eX_old;
			d_rY_ds = eY_old;
			
			// Update with euler scheme
			eX_new = eX_old + d_eX_ds*ds;
			eY_new = eY_old + d_eY_ds*ds;
			nuX_new = nuX_old + d_nuX_ds*ds;
			nuY_new = nuY_old + d_nuY_ds*ds;
			rX_new = rX_old + d_rX_ds*ds;
			rY_new = rY_old + d_rY_ds*ds;
			
			// Normalize nu and e
			d1 = sqrt(eX_new*eX_new + eY_new*eY_new);
			d2 = sqrt(nuX_new*nuX_new + nuY_new*nuY_new);
			
			// Swap variables
			eX_old = eX_new/d1;
			eY_old = eY_new/d1;
			nuX_old = nuX_new/d2;
			nuY_old = nuY_new/d2;
			rX_old = rX_new;
			rY_old = rY_new;
		}
	}
	
	void solveVelProfile(const double ds, const double * kt, const double * kt_dt, const double * tanX, const double * tanY, const double * norX, const double * norY,
						 double * vX, double * vY, double * vTanX, double * vTanY, double * vNorX, double * vNorY, const unsigned int n)
	{
		double rVelX_old = 0.0;
		double rVelY_old = 0.0;
		double eVelX_old = 0.0;
		double eVelY_old = 0.0;
		double nuVelX_old = 0.0;
		double nuVelY_old = 0.0;
		
		double d_eVelX_ds = 0.0;
		double d_eVelY_ds = 0.0;
		double d_nuVelX_ds = 0.0;
		double d_nuVelY_ds = 0.0;
		double d_rVelX_ds = 0.0;
		double d_rVelY_ds = 0.0;
		
		double eVelX_new = 0.0;
		double eVelY_new = 0.0;
		double nuVelX_new = 0.0;
		double nuVelY_new = 0.0;
		double rVelX_new = 0.0;
		double rVelY_new = 0.0;
		
		for(unsigned int i=0; i<n; i++)
		{
			// Store coordiantes
			vX[i] = rVelX_old;
			vY[i] = rVelY_old;
			vTanX[i] = eVelX_old;
			vTanY[i] = eVelY_old;
			vNorX[i] = nuVelX_old;
			vNorY[i] = nuVelY_old;
			
			// Calcualte derivatives
			d_eVelX_ds = kt_dt[i]*norX[i] + kt[i]*nuVelX_old;
			d_eVelY_ds = kt_dt[i]*norY[i] + kt[i]*nuVelY_old;
			d_nuVelX_ds = -kt_dt[i]*tanX[i] - kt[i]*eVelX_old;
			d_nuVelY_ds = -kt_dt[i]*tanY[i] - kt[i]*eVelY_old;
			d_rVelX_ds = eVelX_old;
			d_rVelY_ds = eVelY_old;
			
			// Update with euler scheme
			eVelX_new = eVelX_old + d_eVelX_ds*ds;
			eVelY_new = eVelY_old + d_eVelY_ds*ds;
			nuVelX_new = nuVelX_old + d_nuVelX_ds*ds;
			nuVelY_new = nuVelY_old + d_nuVelY_ds*ds;
			rVelX_new = rVelX_old + d_rVelX_ds*ds;
			rVelY_new = rVelY_old + d_rVelY_ds*ds;
			
			// Swap variables
			eVelX_old = eVelX_new;
			eVelY_old = eVelY_new;
			nuVelX_old = nuVelX_new;
			nuVelY_old = nuVelY_new;
			rVelX_old = rVelX_new;
			rVelY_old = rVelY_new;
		}
	}
};

