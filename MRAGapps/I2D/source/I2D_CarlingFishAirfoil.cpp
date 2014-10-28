/*
 * I2D_CarlingFishAirfoil.cpp
 *
 *  Created on: Feb 16, 2012
 *      Author: mgazzola
 */

#include "I2D_CarlingFishAirfoil.h"

I2D_CarlingFishAirfoil::CarlingFishAirfoil::CarlingFishAirfoil(string _naca, double xm, double ym, double _D, double _T, double phase, double angle_rad, double angleInSpace_rad, double eps, const int LMAX):
Fish(xm,ym,_D,_T,phase,angle_rad,angleInSpace_rad,eps,LMAX), naca(_naca)
{
	assert(naca.size()==4);
	printf("--------------SUMMARY CARLING AIRFOIL--------------\n");
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
	printf("NACA=%s\n",naca.c_str());
	printf("--------------SUMMARY CARLING AIRFOIL--------------\n");

	for(int i = 0; i < N; i++ )
	{
		S[i] = DS*(double)i;
		W[i] = _getWidth(S[i]);
		if(W[i]<0.0){ printf("X-CROSSING IN THE WIDTH!\n"); abort(); }
	}
}

double I2D_CarlingFishAirfoil::CarlingFishAirfoil::_getWidth(const double & ss) const
{
	const double end = 1.0-EXTENSION*EPS;
	const double s = ss - (EXTENSION/2.0)*EPS;
	double width = 0.0;
	if( (s>=0.0) && (s<=end) )
	{
		char buf[2] = {naca[2],naca[3]};
		int twoDigits = atoi(buf);
		const double t = (double)twoDigits/100.0;
		const double c = end;
		const double xc = s/c;
		const double sqrtXC = sqrt(xc);
		const double xc2 = xc*xc;
		const double xc3 = xc*xc*xc;
		const double xc4 = xc*xc*xc*xc;
		width = (t*c/0.2)*( 0.2969*sqrtXC - 0.1260*xc - 0.3516*xc2 + 0.2843*xc3 - 0.1015*xc4 );
	}

	return width;
}

I2D_CarlingFishAirfoil::I2D_CarlingFishAirfoil(ArgumentParser & parser, Grid<W,B>& grid, string _naca, const Real _xm, const Real _ym, const Real D, const Real T, const Real phase, const Real angle, const Real angleInSpace, const Real eps, const Real Uinf[2],
		I2D_PenalizationOperator& penalization, const int LMAX, const int ID, RL::RL_TabularPolicy ** policy):
		I2D_CarlingFish(parser, grid, _xm, _ym, D, T, phase, angle, angleInSpace, eps, Uinf, penalization, LMAX, ID, policy)
{
	assert(shape!=NULL);
	if(shape!=NULL)
	{
		delete shape;
		shape = NULL;
	}

	shape = new CarlingFishAirfoil(_naca, _xm, _ym, D, T, phase, angle, angleInSpace, eps, LMAX);
	shape->updateInSpace(0.0);
}

I2D_CarlingFishAirfoil::~I2D_CarlingFishAirfoil()
{
}

