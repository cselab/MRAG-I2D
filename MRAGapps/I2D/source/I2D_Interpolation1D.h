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

class I2D_Interpolation1D
{
public:
	// Constructors
	I2D_Interpolation1D(){};
	~I2D_Interpolation1D(){};
	
	// Methods
	void naturalCubicSpline(const double * x, const double * y, const unsigned int n, const double * xx, double * yy, const unsigned int nn)
	{
		
		double y2[n];
		double u[n-1];
		double p,qn,sig,un,h,b,a;
		
		y2[0] = 0.0;
		u[0] = 0.0;
		for(unsigned int i=1; i<n-1; i++)
		{
			sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]); 
			p=sig*y2[i-1]+2.0; 
			y2[i]=(sig-1.0)/p; 
			u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]); 
			u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p; 
		}
		
		qn = 0.0;
		un = 0.0;
		y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0); 
		
		
		for(unsigned int k=n-2; k>0; k--)
			y2[k]=y2[k]*y2[k+1]+u[k]; 
		
		for(unsigned int j=0; j<nn; j++)
		{
			unsigned int klo = 0;
			unsigned int khi = n-1;
			unsigned int k = 0; 
			while(khi-klo>1)
			{ 
				k=(khi+klo)>>1; 
				if( x[k]>xx[j])
					khi=k; 
				else
					klo=k; 
			} 
			
			
			h = x[khi]-x[klo]; 
			if(h==0.0){ std::cout << "Interpolation points must be distinct!" << std::endl; abort(); }
			a = (x[khi]-xx[j])/h; 
			b = (xx[j]-x[klo])/h;
			yy[j] = a*y[klo]+b*y[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0; 
		}
	}
	
	void naturalCubicSpline(const double * x, const double * y, const unsigned int n, const double * xx, double * yy, const unsigned int nn, double offset)
	{
		{
			
			double y2[n];
			double u[n-1];
			double p,qn,sig,un,h,b,a;
			
			y2[0] = 0.0;
			u[0] = 0.0;
			for(unsigned int i=1; i<n-1; i++)
			{
				sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]); 
				p=sig*y2[i-1]+2.0; 
				y2[i]=(sig-1.0)/p; 
				u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]); 
				u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p; 
			}
			
			qn = 0.0;
			un = 0.0;
			y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0); 
			
			
			for(unsigned int k=n-2; k>0; k--)
				y2[k]=y2[k]*y2[k+1]+u[k]; 
			
			for(unsigned int j=0; j<nn; j++)
			{
				unsigned int klo = 0;
				unsigned int khi = n-1;
				unsigned int k = 0; 
				while(khi-klo>1)
				{ 
					k=(khi+klo)>>1; 
					if( x[k]>(xx[j]+offset))
						khi=k; 
					else
						klo=k; 
				} 
				
				
				h = x[khi]-x[klo]; 
				if(h==0.0){ std::cout << "Interpolation points must be distinct!" << std::endl; abort(); }
				a = (x[khi]-(xx[j]+offset))/h; 
				b = ((xx[j]+offset)-x[klo])/h;
				yy[j] = a*y[klo]+b*y[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0; 
			}
		}
	}
};



