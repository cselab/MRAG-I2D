/*
 * I2D_CStartLarva.cpp
 *
 *  Created on: Dec 8, 2011
 *      Author: mgazzola
 */

#include "I2D_StefanFishMorph.h"
#include "I2D_Interpolation1D.h"
#include "I2D_Frenet2D.h"
#include "I2D_VectorBlockLab.h"
#include "I2D_GradOfVector.h"
#include <limits>
#include <omp.h>

I2D_StefanFishMorph::StefanFishMorph::StefanFishMorph(double xm, double ym, double _D, double _T, double phase, double tau, double angle_rad, vector<double> WIDTH, vector<double> BASELINE, vector<double> CURVATURE, double angleInSpace_rad, double eps, const int LMAX, const bool isSharp):
  StefanFish(xm, ym, _D, _T, phase, tau, angle_rad, BASELINE, CURVATURE, angleInSpace_rad, eps, LMAX, isSharp)
{
	this->BASELINE = BASELINE;
	this->CURVATURE = CURVATURE;
	this->WIDTH = WIDTH;

	// Checks
	const double dx = 1.0/(double)(FluidBlock2D::sizeX*pow(2.0,LMAX));
	if(this->WIDTH.size()<3){ printf("WIDTH vector too small!\n"); abort(); }
	if(this->WIDTH[1] < dx){ printf("WIDTH[1] < dx!\n"); abort(); }
	if(this->WIDTH[this->WIDTH.size()-1] < dx){ printf("WIDTH[1] < dx!\n"); abort(); }

	printf("--------------SUMMARY STEFAN FISH MORPH--------------\n");
	printf("xm=%f\n",xm);
	printf("ym=%f\n",ym);
	printf("D=%f\n",_D);
	printf("T=%f\n",T);
	printf("phase=%f\n",phase);
	printf("tau=%f\n",tau);
	printf("angle_rad=%f\n",angle_rad);
	printf("angleInSpace_rad=%f\n",angleInSpace_rad);
	printf("eps=%f\n",eps);
	printf("LMAX=%d\n",LMAX);
	printf("WIDTH\n");
	for(unsigned int i=0; i<this->WIDTH.size();i++){ printf("%f, ", this->WIDTH[i]); }
	printf("BASELINE\n");
	for(unsigned int i=0; i<this->BASELINE.size();i++){ printf("%f, ", this->BASELINE[i]); }
	printf("\nCURVATURE\n");
	for(unsigned int i=0; i<this->CURVATURE.size();i++){ printf("%f, ", this->CURVATURE[i]); }
	printf("--------------SUMMARY STEFAN FISH MORPH--------------\n");

	if( this->BASELINE.size()!=this->CURVATURE.size() || this->BASELINE.size()!=6){ printf("Something wrong with BASELINE or CURVATURE!\n"); abort(); }

	this->WH = 0.0;
	this->SB = 0.0;
	this->WT = 0.0;
	this->ST = 0.0;

	for(int i = 0; i < N; i++ )
	{
		S[i] = DS*(double)i;
		W[i] = _getWidth(S[i]);
		if(W[i]<0.0){ printf("X-CROSSING IN THE WIDTH!\n"); abort(); }
	}
}

double I2D_StefanFishMorph::StefanFishMorph::_coxDeBoorRecursion(const int j, const int d, const double * t, const double u) const
{
	const double thres = 1e-6;
	double fac1 = 0.0;
	double fac2 = 0.0;

	if( fabs((double)d-1.0)<thres )
		return (t[j]<=u && u<t[j+1])?1:0;

	if(fabs(t[j+d-1]-t[j])<thres)
		fac1 = 0.0;
	else
		fac1 = (u-t[j])/(t[j+d-1]-t[j]);

	if(fabs(t[j+d]-t[j+1])<thres)
		fac2 = 0.0;
	else
		fac2 = (t[j+d]-u)/(t[j+d]-t[j+1]);

	const double bOld = _coxDeBoorRecursion(j,d-1,t,u);
	const double bOldRight = _coxDeBoorRecursion(j+1,d-1,t,u);

	return fac1*bOld+fac2*bOldRight;
}

double I2D_StefanFishMorph::StefanFishMorph::_getWidth(const double & ss) const
{
	// Very inefficient implementation, but the filling of the W vector
	// is done once at the beginning and that s it, so no big deal!

	// Constants
	const int Nw = 6;
	const int Nx = 2*N;
	const int n = Nw - 1;
	const int ndegree = 3;
	const int d = ndegree + 1;
	const int m = n + d + 1;
	double t[m];

	// Allocate control points
	double Pwx[Nw];
	double Pwy[Nw];

	// Allocate x and y vectors and set them to zero
	double xr[Nx];
	double yr[Nx];
	for(int i=0; i<Nx; i++){ xr[i] = 0.0; yr[i] = 0.0; }

	// Define control points, first the ends
	Pwx[0] = 0.0;
	Pwy[0] = 0.0;
	Pwx[Nw-1] = 1.0;
	Pwy[Nw-1] = 0.0;

	// Now the interior
	int counter = 0;
	const double dxw = 1.0/((double)Nw-3.0);
	for(int i=1; i<Nw-1; i++)
	{
		Pwx[i] = (i-1)*dxw;
		Pwy[i] = WIDTH[counter];
		counter++;
	}

	// Build t vector
	for(int i=0; i<m; i++)
	{
		if(i<d)
			t[i] = 0.0;
		else if(i>n)
			t[i] = n - d + 2.0;
		else
			t[i] = i - d + 1.0;
	}

	for(int i=0; i<m; i++)
		t[i] /= t[m-1];

	// Build the spline
	const double du = t[m-1]/((double)Nx-1.0);
	for(int i=0; i<Nx; i++)
	{
		const double u = (double)i*du;
		for(int j=0; j<n+1; j++)
		{
			xr[i] += Pwx[j]*_coxDeBoorRecursion(j,d,t,u);
			yr[i] += Pwy[j]*_coxDeBoorRecursion(j,d,t,u);
		}
	}

	xr[Nx-1] = Pwx[Nw-1];
	yr[Nx-1] = Pwy[Nw-1];

	// Scale down to account for eps..
	for(int i=0; i<Nx; i++)
	{
		xr[i] *= 1.0-EXTENSION*EPS;
		yr[i] *= 1.0-EXTENSION*EPS;
	}

	// Generate profile via linear interpolation
	const double end = 1.0-EXTENSION*EPS;
	const double s = ss - (EXTENSION/2.0)*EPS;
	double width = 0.0;
	if( (s>=0.0) && (s<=end) )
	{
		int idx = 0;
		for(int i=0; i<Nx; i++)
		{
			if(xr[i]>s)
			{
				idx = max(0,i-1);
				break;
			}
		}
		const int idxPlus = min(idx+1,Nx-1);
		width = yr[idx] + (yr[idxPlus]-yr[idx])*(s-xr[idx])/(xr[idxPlus]-xr[idx]);
	}

	return width;
}

I2D_StefanFishMorph::I2D_StefanFishMorph(ArgumentParser & parser, Grid<W,B>& grid, const Real _xm, const Real _ym, const Real D, const Real T, const Real phase, const Real tau, const Real angle, vector<double> WIDTH, vector<double> BASELINE, vector<double> CURVATURE, const Real angleInSpace, const Real eps, const Real Uinf[2], I2D_PenalizationOperator& penalization, const int LMAX, const int ID):
		I2D_StefanFish(parser, grid, _xm, _ym, D, T, phase, tau, angle, BASELINE, CURVATURE, angleInSpace, eps, Uinf, penalization, LMAX, ID)
{
	assert(shape!=NULL);
	if(shape!=NULL)
	{
		delete shape;
		shape = NULL;
	}

	const bool isSharp = parser("-sharp").asBool();

	shape = new StefanFishMorph(_xm, _ym, D, T, phase, tau, angle, WIDTH, BASELINE, CURVATURE, angleInSpace, eps, LMAX, isSharp);
	shape->updateInSpace(0.0);
	//shape->clear();
}

