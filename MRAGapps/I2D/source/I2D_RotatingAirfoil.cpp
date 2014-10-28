/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */
#include "I2D_RotatingAirfoil.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;



// ------------------- MEMBER FUNCTIONS OF CLASS "NACA4gen"--------------------------------

Real I2D_RotatingAirfoil::NACA4gen::getNoseDistance(Real xx[2]) const
{
	assert(xx[0]>=-2.0*epsilon && xx[0]<0.05*leadingEdgeRadius);

	Real dist = sqrt( (leadingEdgeCenter[0]-xx[0])*(leadingEdgeCenter[0]-xx[0]) + (leadingEdgeCenter[1]-xx[1])*(leadingEdgeCenter[1]-xx[1]) ) - leadingEdgeRadius;

	return dist;
}


Real I2D_RotatingAirfoil::NACA4gen::getYUpperProfile(Real xx[2], Real & dist) const
{
	// nPanels = number of base points --> prerequisite: nPanels == xu.size == yu.size
	assert( xu.size() == yu.size() );
	assert( xu.size() == nPanels );
	assert(xx[0]>=0.0 && xx[0]<=chord);

	const Real dx = chord/((Real)nPanels-1.0); // interval length
	const int midIdx = floor(xx[0]/dx); // get base point before sought x-coordinate x[0]
	// base points to iterate over
	const int startIdx = ((midIdx-maxDiffX)<0) ? 0 : (midIdx-maxDiffX); // symmetrical profiles: maxDiff = 10
	const int endIdx = ((midIdx+maxDiffX)>(int)(nPanels))?nPanels:(midIdx+maxDiffX);

	assert(startIdx>=0);
	assert(endIdx>=0);
	assert(startIdx<=xu.size());
	assert(endIdx<=xu.size());
	assert(startIdx<endIdx);

	int lowIdx = 0;
	int highIdx = 0;
	for(int i=startIdx; i<endIdx; i++)
	{
		if(xu[i] > xx[0])
		{
			highIdx = i;
			lowIdx = i-1;
			break;
		}
	}

	assert(lowIdx>=0);
	assert(highIdx>=0);
	assert(lowIdx<=xu.size());
	assert(highIdx<=xu.size());
	assert(lowIdx<highIdx);

	// interpolate y value of the upper airfoil shape at position xx[0]
	const Real v = yu[lowIdx] + (yu[highIdx] - yu[lowIdx])/(xu[highIdx] - xu[lowIdx]) * (xx[0] - xu[lowIdx]);

	// base points to iterate over
	const int startIdxEps = ((lowIdx-maxDiffEpsilon)<0) ? 0 : (lowIdx-maxDiffEpsilon);
	const int endIdxEps = ((highIdx+maxDiffEpsilon)>(int)(nPanels))?nPanels:(highIdx+maxDiffEpsilon);

	assert(startIdxEps>=0);
	assert(endIdxEps>=0);
	assert(startIdxEps<=xu.size());
	assert(endIdxEps<=xu.size());
	assert(startIdxEps<endIdxEps);
	Real myDist = std::numeric_limits<Real>::max(); // ensure assignment in for loop (myDist for sure smaller than max value of Real)
	// calculate minimal distance between upper shape and point xx
	for(int i=startIdxEps; i<endIdxEps; i++)
	{
		const Real dxx = (xu[i]-xx[0]);
		const Real dyy = (yu[i]-xx[1]);
		myDist = min(myDist,(Real)sqrt(dxx*dxx+dyy*dyy));
	}
	dist = myDist;
	assert( dist!=std::numeric_limits<Real>::max() );

	return v;
}

Real  I2D_RotatingAirfoil::NACA4gen::getYLowerProfile(Real xx[2], Real & dist) const
{
	assert( xl.size() == yl.size() );
	assert( xl.size() == nPanels );
	assert(xx[0]>=0.0 && xx[0]<=chord);

	const Real dx = chord/((Real)nPanels-1.0);
	const int midIdx = floor(xx[0]/dx);
	const int startIdx = ((midIdx-maxDiffX)<0) ? 0 : (midIdx-maxDiffX);
	const int endIdx = ((midIdx+maxDiffX)>(int)(nPanels))?nPanels:(midIdx+maxDiffX);

	assert(startIdx>=0);
	assert(endIdx>=0);
	assert(startIdx<=xl.size());
	assert(endIdx<=xl.size());
	assert(startIdx<endIdx);

	int lowIdx = 0;
	int highIdx = 0;
	for(int i=startIdx; i<endIdx; i++)
	{
		if(xl[i] > xx[0])
		{
			highIdx = i;
			lowIdx = i-1;
			break;
		}
	}

	assert(lowIdx>=0);
	assert(highIdx>=0);
	assert(lowIdx<=xl.size());
	assert(highIdx<=xl.size());
	assert(lowIdx<highIdx);

	const Real v = yl[lowIdx] + (yl[highIdx] - yl[lowIdx])/(xl[highIdx] - xl[lowIdx]) * (xx[0] - xl[lowIdx]);

	const int startIdxEps = ((lowIdx-maxDiffEpsilon)<0) ? 0 : (lowIdx-maxDiffEpsilon);
	const int endIdxEps = ((highIdx+maxDiffEpsilon)>(int)(nPanels))?nPanels:(highIdx+maxDiffEpsilon);

	assert(startIdxEps>=0);
	assert(endIdxEps>=0);
	assert(startIdxEps<=xl.size());
	assert(endIdxEps<=xl.size());
	assert(startIdxEps<endIdxEps);
	Real myDist = std::numeric_limits<Real>::max();
	for(int i=startIdxEps; i<endIdxEps; i++)
	{
		const Real dxx = (xl[i]-xx[0]);
		const Real dyy = (yl[i]-xx[1]);
		myDist = min(myDist,(Real)sqrt(dxx*dxx+dyy*dyy));
	}
	dist = myDist;
	assert( dist!=std::numeric_limits<Real>::max() );

	return v;
}

Real I2D_RotatingAirfoil:: NACA4gen::sdf(Real xo[2]) const
{
	// xo not within [0,chord]
	if (xo[0]<-2.0*epsilon){ return sqrt(pow(xo[0],2) + pow(xo[1], 2)); } // return distance to origin

	if (xo[0]>chord){ return sqrt(pow(xo[0]-chord,2) + pow(xo[1], 2)); } // return distance to trailing edge

	if (xo[0]>=-2.0*epsilon && xo[0]<0.05*leadingEdgeRadius){ return getNoseDistance(xo); } // return distance to leading edge

	// xo within [0,chord]
	Real distU = 0.0;
	Real distL = 0.0;
	const Real high = this->getYUpperProfile(xo,distU);
	const Real low = this->getYLowerProfile(xo,distL);

	assert(high >= low);

	const Real s = (Real)(xo[1]<low || xo[1]>high)*2.0 - 1.0; // s=-1, if xo inside body, s=1 if xo outside body
	const Real D = min(distU,distL); // positive inside body, negative outside body

	return s*D;
}



I2D_RotatingAirfoil::NACA4gen::NACA4gen(Real chord_in,unsigned int nPanels_in,bool isFiniteTE,Real epsilon_in, const int _d1, const int _d2, const int _d3d4)
:d1(_d1), d2(_d2), d3d4(_d3d4)
{
	x.clear();
	yt.clear();
	xu.clear();
	yu.clear();
	xl.clear();
	yl.clear();
	yC.clear();

	/*
    char cstrD1;
	char cstrD2;
	char cstrD3D4[3];
	string strD1;
	string strD2;
	string strD3D4;

    cstrD1 = wingType[0];
	cstrD2 = wingType[1];
    cstrD3D4[0] = wingType[2];
	cstrD3D4[1] = wingType[3];

	strD1 = cstrD1;
	strD2 = cstrD2;
	strD3D4 = cstrD3D4;

	int d1 = atoi(strD1.c_str());
	int d2 = atoi(strD2.c_str());
	int d3d4 = atoi(strD3D4.c_str());
	 */

	/*
    assert(d1>=0 && d1<=9);
    assert(d2>=0 && d2<=9);
	assert(d3d4>=0 && d3d4<=99);
	 */

	printf("\nNACA-%d%d%d\n",d1,d2,d3d4);

	nPanels = nPanels_in;
	chord = chord_in;
	epsilon = epsilon_in;

	const Real t = (Real)d3d4/100.0;
	const Real m = (Real)d1/100.0;
	Real p = (Real)d2/10.0;

	// coefficients for symmetrical 4-digit NACA airfoil
	const Real a0 = 0.2969;
	const Real a1 = -0.1260;
	const Real a2 = -0.3516;
	const Real a3 = 0.2843;
	Real a4 = 0.0;

	if (isFiniteTE)
		a4=-0.1015; // For finite thick trailing edge
	else
		a4=-0.1036;  // For zero thick trailing edge

	const Real dx = 1.0/((Real)nPanels-1.0);
	for(unsigned int i=0; i<nPanels; i++)
	{
		x.push_back( (Real)i*dx );
		yt.push_back( (t/0.2)*(a0*sqrt(x[i])+a1*x[i]+a2*x[i]*x[i]+a3*x[i]*x[i]*x[i]+a4*x[i]*x[i]*x[i]*x[i]) );
	}

	Real diffX = 0.0;
	if(p==0.0) // symmetrical airfoil
	{
		for(unsigned int i=0; i<nPanels; i++)
		{
			xu.push_back( x[i] );
			yu.push_back( yt[i] );

			xl.push_back( x[i] );
			yl.push_back( -yt[i] );

			yC.push_back(0.0);
		}
	}
	else // non-symmetrical airfoil
	{
		for(unsigned int i=0; i<nPanels; i++)
		{
			Real zc = 0.0;
			Real dyc_dx = 0.0;

			if(x[i]<=p)
			{
				zc=(m/(p*p))*((Real)2.0*p*x[i]-x[i]*x[i]);
				dyc_dx=(m/(p*p))*((Real)2.0*p-(Real)2.0*x[i]);
			}
			else
			{
				zc=(m/((1.0-p)*(1.0-p)))*((1.0-2.0*p)+2.0*p*x[i]-x[i]*x[i]);
				dyc_dx=(m/((1.0-p)*(1.0-p)))*(2.0*p-2.0*x[i]);
			}

			const Real theta=atan(dyc_dx);

			yC.push_back(zc);

			xu.push_back( x[i] - yt[i]*sin(theta) );
			diffX = max(diffX,(Real)fabs(x[i]-xu[i]));
			yu.push_back( zc + yt[i]*cos(theta) );

			xl.push_back( x[i] + yt[i]*sin(theta) );
			diffX = max(diffX,(Real)fabs(x[i]-xl[i]));
			yl.push_back( zc - yt[i]*cos(theta) );
		}
	}

	assert(nPanels==x.size());
	assert(nPanels==yt.size());
	assert(nPanels==xu.size());
	assert(nPanels==yu.size());
	assert(nPanels==xl.size());
	assert(nPanels==yl.size());
	assert(nPanels==yC.size());

	maxDiffX = ceil(diffX/dx)+10; // symmetrical profiles: maxDiffX = 10; ceil(x) returns the smallest int value not smaller than x
	maxDiffEpsilon = 2.0*ceil(epsilon/dx);

	for(unsigned int i=0; i<nPanels; i++)
	{
		x[i] *= chord;
		yt[i] *= chord;
		xu[i] *= chord;
		yu[i] *= chord;
		xl[i] *= chord;
		yl[i] *= chord;
		yC[i] *= chord;
	}

	// assure assignment of one of the boxes (values are > 0 but smaller than max of Real for sure)

	maxBoxX = 0.0;
	minBoxX = std::numeric_limits<Real>::max(); // maximum value for a variable of type Real
	maxBoxY = 0.0;
	minBoxY = std::numeric_limits<Real>::max();

	leadingEdgeRadius = 1.1019*chord*t*t; // general formula

	//std::cout << x[0] << std::endl;
	//std::cout << x[1] << std::endl;
	//std::cout << yC[0] << std::endl;
	//std::cout << yC[1] << std::endl;

	const Real noseAngle = atan2((yC[1]-yC[0]),(x[1]-x[0]));
	leadingEdgeCenter[0] = leadingEdgeRadius*cos(noseAngle);
	leadingEdgeCenter[1] = leadingEdgeRadius*sin(noseAngle);

	//std::cout << noseAngle << std::endl;
	//std::cout << leadingEdgeCenter[0] << std::endl;
	//std::cout << leadingEdgeCenter[1] << std::endl;

	for(unsigned int i=0; i<nPanels; i++)
	{
		maxBoxX = std::max((Real)maxBoxX,(Real)xu[i]);
		maxBoxX = std::max((Real)maxBoxX,(Real)xl[i]);

		minBoxX = std::min((Real)minBoxX,(Real)xu[i]);
		minBoxX = std::min((Real)minBoxX,(Real)xl[i]);

		maxBoxY = std::max((Real)maxBoxY,(Real)yu[i]);
		minBoxY = std::min((Real)minBoxY,(Real)yl[i]);
	}

	//printf("maxBoxX %f\n",maxBoxX);
	//printf("minBoxX %f\n",minBoxX);
	//printf("maxBoxY %f\n",maxBoxY);
	//printf("minBoxY %f\n",minBoxY);
}

I2D_RotatingAirfoil::NACA4gen::~NACA4gen(){}

// -------------------END OF MEMBER FUNCTIONS OF CLASS "NACA4gen"--------------------------------


// -------------------MEMBER FUNCTIONS OF CLASS "DISCRETISED WING" --------------------------------


void I2D_RotatingAirfoil::DiscretizedWing::check_mapping() const // check if O2W and W2O are inverse matrices
{
	Real A[3][3];

	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
		{
			A[i][j] = 0;

			for(int d=0; d<3; d++)
				A[i][j] += O2W[i][d]* W2O[d][j];
		}

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if (i==j)
				assert(fabs(A[i][j]-1)<1e-5);
			else
				assert(fabs(A[i][j]-0)<1e-5);
		}
	}
}

void I2D_RotatingAirfoil::DiscretizedWing::_w2o(const Real xw[2], Real xo[2]) const
{
	xo[0] = W2O[0][0]*xw[0] + W2O[0][1]*xw[1] + W2O[0][2];
	xo[1] = W2O[1][0]*xw[0] + W2O[1][1]*xw[1] + W2O[1][2];
}

void I2D_RotatingAirfoil::DiscretizedWing::_o2w(const Real xo[2], Real xw[2]) const
{
	xw[0] = O2W[0][0]*xo[0] + O2W[0][1]*xo[1] + O2W[0][2];
	xw[1] = O2W[1][0]*xo[0] + O2W[1][1]*xo[1] + O2W[1][2];
}

I2D_RotatingAirfoil::DiscretizedWing::DiscretizedWing
       (Real _tx,  Real _ty,  Real width, Real epsilon, Real _angle,
    	Real _self_angle,const int d1, const int d2, const int d3d4, const  Real l) :

		naca(width,5000,false,epsilon,d1,d2,d3d4),
		scaling(1.0),
		epsilon(epsilon),
		width(width),
		tx(_tx), ty(_ty), s(1.0), angular_velocity(0), rho(1.0), m(0.0), J(0.0), l(0.0),
        vx(0.0), vy(0.0), self_angle(_self_angle), angle(_angle)

{
	ca = cos(self_angle*M_PI/180.);
	sa = sin(self_angle*M_PI/180.);

	printf("self_angle=%e\n",self_angle);
	printf("angle=%e\n",angle);


    //Inverse Rotation Matrix
	W2O[0][0] = +1/s*ca;
	W2O[0][1] = +1/s*sa;
	W2O[0][2] = +1/s*(-tx*ca - ty*sa);
	W2O[1][0] = -1/s*sa;
	W2O[1][1] = +1/s*ca;
	W2O[1][2] = +1/s*(+tx*sa - ty*ca);
	W2O[2][0] = 0;
	W2O[2][1] = 0;
	W2O[2][2] = 1;

	//Rotation Matrix
	O2W[0][0] = +s*ca;
	O2W[0][1] = -s*sa;
	O2W[0][2] = tx;
	O2W[1][0] = s*sa;
	O2W[1][1] = s*ca;
	O2W[1][2] = ty;
	O2W[2][0] = 0;
	O2W[2][1] = 0;
	O2W[2][2] = 1;

	check_mapping();
	printf("checked map in Discretized wing constructor \n");
	printf("tx=%e\n", tx);
	printf("ty=%e\n", ty);


}

void I2D_RotatingAirfoil::DiscretizedWing::update_all (const double dt, const double t, const double perist, const double _D, const double _theta)
{

	angle=angle+_theta;
    vx=-perist*_D/2.0*sin(angle*M_PI/180.);
    vy= perist*_D/2.0*cos(angle*M_PI/180.);

    tx=vx*dt;
    ty=vy*dt;

	printf("timestep dt=%e\n", dt);
	printf("new global angle of airfoil =%e\n\n", angle);
	printf("new local angle of airfoil =%e\n\n", self_angle);

	ca = cos(self_angle*M_PI/180.);
	sa = sin(self_angle*M_PI/180.);

	W2O[0][0] = +1/s*ca;
	W2O[0][1] = +1/s*sa;
	W2O[0][2] = +1/s*(-tx*ca - ty*sa);
	W2O[1][0] = -1/s*sa;
	W2O[1][1] = +1/s*ca;
	W2O[1][2] = +1/s*(+tx*sa - ty*ca);
	W2O[2][0] = 0;
	W2O[2][1] = 0;
	W2O[2][2] = 1;

	O2W[0][0] = +s*ca;
	O2W[0][1] = -s*sa;
	O2W[0][2] = tx;
	O2W[1][0] = s*sa;
	O2W[1][1] = s*ca;
	O2W[1][2] = ty;
	O2W[2][0] = 0;
	O2W[2][1] = 0;
	O2W[2][2] = 1;
	check_mapping();
	printf("checked map in Discretized wing update all \n");

}

void I2D_RotatingAirfoil::DiscretizedWing::getAerodynamicCenter(Real cor[2]) const // calculates position and orientation of airfoil
{
	cor[0] = 0.25*width;
	cor[1] = 0.0;

	Real corw[2] = {0.0,0.0};
	_o2w(cor, corw); // rotate cor according to angle

	cor[0] = corw[0];
	cor[1] = corw[1];
}

Real I2D_RotatingAirfoil::DiscretizedWing::sdf(Real xw[2]) const
{
	Real xo[2];
	_w2o(xw, xo);
	return scaling*naca.sdf(xo);
}

void I2D_RotatingAirfoil::DiscretizedWing::bbox(const Real eps, Real xmin[2], Real xmax[2]) const
{
	assert(eps>0);

	Real v1[2] = {naca.getMinX(),naca.getMinY()};
	Real v2[2] = {naca.getMinX(),naca.getMaxY()};
	Real v3[2] = {naca.getMaxX(),naca.getMinY()};
	Real v4[2] = {naca.getMaxX(),naca.getMaxY()};

	Real xw1[2], xw2[2], xw3[2], xw4[2];

	_o2w(v1, xw1);
	_o2w(v2, xw2);
	_o2w(v3, xw3);
	_o2w(v4, xw4);


	xmin[0] = min((Real)min((Real)min((Real)xw1[0],(Real)xw2[0]),(Real)xw3[0]),(Real)xw4[0]); // xmin
	xmin[1] = min((Real)min((Real)min((Real)xw1[1],(Real)xw2[1]),(Real)xw3[1]),(Real)xw4[1]); // ymin
	xmax[0] = max((Real)max((Real)max((Real)xw1[0],(Real)xw2[0]),(Real)xw3[0]),(Real)xw4[0]); // xmax
	xmax[1] = max((Real)max((Real)max((Real)xw1[1],(Real)xw2[1]),(Real)xw3[1]),(Real)xw4[1]); // ymax

	xmin[0] -= 2*eps;
	xmin[1] -= 2*eps;
	xmax[0] += 2*eps;
	xmax[1] += 2*eps;

	assert(xmin[0]<=xmax[0]);
	assert(xmin[1]<=xmax[1]);
}

void I2D_RotatingAirfoil::DiscretizedWing::save(FILE * f) const
{
	fprintf(f, "xm: %20.20e\n", tx);
	fprintf(f, "ym: %20.20e\n", ty);
	fprintf(f, "angular_velocity: %20.20e\n", angular_velocity);
	fprintf(f, "angle: %20.20e\n", angle);
}

void I2D_RotatingAirfoil::DiscretizedWing::restart(FILE * f)
{
	float val;

	fscanf(f, "xm: %e\n", &val);
	tx = val;
	printf("DiscretizedWing::restart(): xm is %e\n", tx);

	fscanf(f, "ym: %e\n", &val);
	ty = val;
	printf("DiscretizedWing::restart(): ym is %e\n", ty);

	fscanf(f, "angular_velocity: %e\n", &val);
	angular_velocity = val;
	printf("DiscretizedWing::restart(): angular_velocity is %e\n", angular_velocity);

	fscanf(f, "angle: %e\n", &val);
	angle = val;
	printf("DiscretizedWing::restart(): angle is %e\n", angle);
}

Real I2D_RotatingAirfoil::DiscretizedWing::_mollified_heaviside(const Real x, const Real eps) const
{
	const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));
	return 0.5+0.5*cos(alpha);
}

Real I2D_RotatingAirfoil::DiscretizedWing::sample(Real p[2], const Real eps) const // calculate characteristic function at point (p[0],p[1])
{
	Real dist = sdf(p);
	return _mollified_heaviside(dist,eps);
}




// -------------------END OF MEMBER FUNCTIONS OF CLASS "DISCRETISED WING"---------------------------------



// -------------------MEMBER FUNCTIONS OF CLASS "TOO MANY AIRFOILS"---------------------------------------


I2D_RotatingAirfoil::ToomanyAirfoils::ToomanyAirfoils
       (const Real _xm , const Real _ym, const Real _D,
		const Real _c  , const int d1,   const int d2,
		const int d3d4 , const Real _l,  const Real epsilon):

		theta(0.0), xm(_xm), ym(_ym), D(_D),rho(1.0),
		m(0.0),J(0.0),angular_velocity(0.0),epsilon(epsilon)

{
	wing_1=new I2D_RotatingAirfoil::DiscretizedWing(_xm+_D/2,_ym,_c,epsilon,0.0,-120., d1,d2,d3d4,_l);
	wing_2=new I2D_RotatingAirfoil::DiscretizedWing(_xm-(cos(M_PI/3)*_D/2.0),_ym+(sin(M_PI/3)*_D/2.0),_c,epsilon,120.0,0.0,d1,d2,d3d4,_l);
	wing_3=new I2D_RotatingAirfoil::DiscretizedWing(_xm-(sin(M_PI/6)*_D/2.0),_ym-(cos(M_PI/6)*_D/2.0),_c,epsilon,240.0,120.0,d1,d2,d3d4,_l);
}

Real I2D_RotatingAirfoil::ToomanyAirfoils::sample(Real p[2], const Real eps) const // calculate characteristic function at point (p[0],p[1])
{
	const Real dist1 = wing_1->sdf(p);
	const Real FirstWing=_mollified_heaviside(dist1,eps);
	Real myMax= FirstWing;

	const Real dist2 = wing_2->sdf(p);
	const Real SecondWing= _mollified_heaviside(dist2,eps);
	myMax=max(myMax, SecondWing);


	const Real dist3 = wing_3->sdf(p);
	const Real ThirdWing= _mollified_heaviside(dist3,eps);
	myMax=max(myMax,ThirdWing);

	return myMax;
}

Real I2D_RotatingAirfoil::ToomanyAirfoils::_mollified_heaviside(const Real x, const Real eps) const
{
	const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));
	return 0.5+0.5*cos(alpha);
}

void I2D_RotatingAirfoil::ToomanyAirfoils::bbox(const Real eps, Real xmin[2], Real xmax[2]) const
{
	    xmin[0] = xm - ((D/2.0)+wing_1->width);
		xmax[0] = xm + ((D/2.0)+wing_1->width);
		xmin[1] = ym - ((D/2.0)+wing_1->width);
		xmax[1] = ym + ((D/2.0)+wing_1->width);

		Real v1[2] = { xmin[0], xmin[1] }; /// bottom left corner
		Real v2[2] = { xmax[0], xmax[1] }; /// upper right corner
		Real v3[2] = { xmin[0], xmax[1] }; /// upper left corner
		Real v4[2] = { xmax[0], xmin[1] }; /// bottom right corner

		rotate(v1, theta); /// rotate the box by body orientation
		rotate(v2, theta);
		rotate(v3, theta);
		rotate(v4, theta);

		xmin[0] = min((Real)min((Real)min(v1[0], v2[0]), v3[0]), v4[0]); /// min x
		xmax[0] = max((Real)max((Real)max(v1[0], v2[0]), v3[0]), v4[0]); /// max x
		xmin[1] = min((Real)min((Real)min(v1[1], v2[1]), v3[1]), v4[1]); /// min y
		xmax[1] = max((Real)max((Real)max(v1[1], v2[1]), v3[1]), v4[1]); /// max y

		assert(eps >= 0);

		xmin[0] -= 2*eps; /// add smoothing length
		xmin[1] -= 2*eps;
		xmax[0] += 2*eps;
		xmax[1] += 2*eps;

		assert(xmin[0] <= xmax[0]);
		assert(xmin[1] <= xmax[1]);
}

void I2D_RotatingAirfoil::ToomanyAirfoils::rotate(Real v[2], Real _theta) const
{
	const Real _th=_theta*M_PI/180.;
	const Real a00 = cos(_th);
	const Real a01 = -sin(_th);
	const Real a10 = sin(_th);
	const Real a11 = cos(_th);

	const Real xv = v[0];
	const Real yv = v[1];

	v[0] = a00*(xv-xm) + a01*(yv-ym);
	v[1] = a10*(xv-xm) + a11*(yv-ym);

	v[0] += xm;
	v[1] += ym;
}

void I2D_RotatingAirfoil::ToomanyAirfoils::getAerodynamicCenter(Real cor[2]) const
{
	wing_1->getAerodynamicCenter(cor);
    wing_2->getAerodynamicCenter(cor);
	wing_3->getAerodynamicCenter(cor);
}

void I2D_RotatingAirfoil::ToomanyAirfoils::update_all(double dt, double t)
{
	printf("angular_velocity=%e\n\n",angular_velocity);
	theta +=angular_velocity*dt;
	printf("theta=%e\n\n",theta);

	wing_1->update_all(dt,t,angular_velocity,D,theta);
	wing_2->update_all(dt,t,angular_velocity,D,theta);
	wing_3->update_all(dt,t,angular_velocity,D,theta);

}

void I2D_RotatingAirfoil::ToomanyAirfoils::save(FILE * f) const
{
	fprintf(f, "xm: %20.20e\n", xm);
	fprintf(f, "ym: %20.20e\n", ym);
	fprintf(f, "angular_velocity: %20.20e\n", angular_velocity);
	fprintf(f, "angle: %20.20e\n", theta);
}

void I2D_RotatingAirfoil::ToomanyAirfoils::restart(FILE * f)// TO DO: change this
{
	float val;

		fscanf(f, "xm: %e\n", &val);
		xm = val;
		printf("Wheel::restart(): xm is %e\n", xm);

		fscanf(f, "ym: %e\n", &val);
		ym = val;
		printf("Wheel::restart(): ym is %e\n", ym);

		fscanf(f, "angular_velocity: %e\n", &val);
		angular_velocity = val;
		printf("Wheel::restart(): angular_velocity is %e\n", angular_velocity);

		fscanf(f, "angle: %e\n", &val);
		theta = val;
		printf("Wheel::restart(): angle is %e\n", theta);



}

// -------------------END OF MEMBER FUNCTIONS OF CLASS "TOO MANY AIRFOILS"--------------------------------


// --------------------------------------------NAMESPACE STUFF--------------------------------------------


namespace AirfoilStuff
{

struct FillBlocks
{
	Real eps;
	I2D_RotatingAirfoil::ToomanyAirfoils * pterygio;

	FillBlocks(Real eps,I2D_RotatingAirfoil::ToomanyAirfoils * _pterygio ): eps(eps)
	{
		pterygio=_pterygio;
	}

	FillBlocks(const FillBlocks& c): eps(c.eps)
	{
		pterygio=c.pterygio;
	}

	static bool _is_touching(Real eps, const I2D_RotatingAirfoil::ToomanyAirfoils * wing, const BlockInfo& info)
	{

		Real wing_min[2], wing_max[2];
		wing->bbox(eps, wing_min, wing_max); // store shapes' xmin, ymin in wing_min and shapes' xmax, ymax in wing_max

		Real min_pos[2], max_pos[2];
		info.pos(min_pos, 0,0);
		info.pos(max_pos, FluidBlock2D::sizeX-1, FluidBlock2D::sizeY-1);

		Real intersection[2][2] = {
				max(min_pos[0], wing_min[0]), min(max_pos[0],  wing_max[0]),
				max(min_pos[1], wing_min[1]), min(max_pos[1],  wing_max[1])
		};

		return
				intersection[0][1]-intersection[0][0]>0 &&
				intersection[1][1]-intersection[1][0]>0 ;

	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{
		bool bEmpty = true;

		if(_is_touching(eps, pterygio, info))
		{
			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
				{
					Real p[2];
					info.pos(p, ix, iy);

					b(ix,iy).tmp = max( pterygio->sample(p, eps), b(ix,iy).tmp );
				}

			bEmpty = false;
		}
	}


};


struct ComputeAll
{
	Real Uinf[2];
	double eps;
	I2D_RotatingAirfoil::ToomanyAirfoils * pterygio;
	map<int, vector<double> >& b2sum;
	map<int, bool>& b2nonempty;

	ComputeAll(double eps,I2D_RotatingAirfoil::ToomanyAirfoils * _pterygio, map<int, vector<double> >& b2sum, Real Uinf[2], map<int, bool>& b2nonempty): eps(eps), b2sum(b2sum), b2nonempty(b2nonempty)
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];

		pterygio=_pterygio;
	}

	ComputeAll(const ComputeAll& c): eps(c.eps), b2sum(c.b2sum), b2nonempty(c.b2nonempty)
	{
		Uinf[0] = c.Uinf[0];
		Uinf[1] = c.Uinf[1];
		pterygio=c.pterygio;
	}

	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{
		bool bNonEmpty = false;

		if(FillBlocks::_is_touching(eps, pterygio, info))
		{
			double omegabar = 0;
			double J=0;
			double mass=0;

			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
				{
					Real p[2];
					info.pos(p, ix, iy);

					const double Xs = pterygio->sample(p, eps);
					bNonEmpty |= Xs>0;

					mass += Xs;
					J += Xs*(pow((p[0]-pterygio->xm),2) + pow((p[0]-pterygio->ym),2));
					omegabar += Xs*((b(ix, iy).u[1]+Uinf[1])*(p[0]-pterygio->xm)-(b(ix, iy).u[0]+Uinf[0])*(p[1]-pterygio->ym));

				}

			assert(b2sum.find(info.blockID) != b2sum.end());
			assert(b2nonempty.find(info.blockID) != b2nonempty.end());

			b2sum[info.blockID][0] = mass*info.h[0]*info.h[0]*pterygio->rho;
			b2sum[info.blockID][1] = J*info.h[0]*info.h[0]*pterygio->rho;
			b2sum[info.blockID][2] = omegabar*info.h[0]*info.h[0]*pterygio->rho;
			b2nonempty[info.blockID] = bNonEmpty;
		}
	}

};

struct FillVelblocks
{
	double eps;
	vector<pair< BlockInfo, VelocityBlock *> >& workitems;
	I2D_RotatingAirfoil::ToomanyAirfoils * pterygio;

	FillVelblocks(vector<pair< BlockInfo, VelocityBlock *> >& workitems, double eps, I2D_RotatingAirfoil::ToomanyAirfoils * _pterygio):workitems(workitems), eps(eps)
	{
		pterygio=_pterygio;
	}

	FillVelblocks(const FillVelblocks& c): workitems(c.workitems), eps(c.eps)
	{
		pterygio=c.pterygio;
	}

	inline void operator()(blocked_range<int> range) const
	{
		for(int iblock=range.begin(); iblock<range.end(); iblock++)
		{
			BlockInfo info = workitems[iblock].first;
			VelocityBlock * u_desired = workitems[iblock].second;

			const Real xm = pterygio->xm;
			const Real ym =	pterygio->ym;
			const Real av = pterygio->angular_velocity;

			if(FillBlocks::_is_touching(eps, pterygio, info))
				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						Real p[2];
						info.pos(p, ix, iy);

						const Real Xs = pterygio->sample(p, eps);

						if (Xs > 0)
						{
							u_desired->u[0][iy][ix] = - av*(p[1]-ym) ;
							u_desired->u[1][iy][ix] = + av*(p[0]-xm) ;
						}
						else
						{
							u_desired->u[0][iy][ix] = 0.0;
							u_desired->u[1][iy][ix] = 0.0;
						}
					}
		}

	}
};
}
// ---------------------------------  END OF NAMESPACE STUFF    -------------------------------------------



// ------------------- MEMBER FUNCTIONS OF CLASS "I2D_RotatingAirfoil"--------------------------------

I2D_RotatingAirfoil::I2D_RotatingAirfoil
        (    Grid<W,B>& grid, ArgumentParser& parser, const Real _xm,
		const Real _ym, const Real _D, const Real _c, const int d1,    const int d2,
		const int d3d4, const Real _l, const Real _eps, const Real Uinf[2],
		I2D_PenalizationOperator& penalization):

I2D_FloatingObstacleOperator(parser, grid, D, eps, Uinf, penalization)

{
	this->Uinf[0] = Uinf[0];
	this->Uinf[1] = Uinf[1];
    airfoil=new I2D_RotatingAirfoil::ToomanyAirfoils(_xm ,_ym, _D,_c, d1, d2, d3d4,_l,eps);
}

I2D_RotatingAirfoil::~I2D_RotatingAirfoil()
{}


void I2D_RotatingAirfoil::characteristic_function()
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	AirfoilStuff::FillBlocks fill(eps,airfoil);
	block_processing.process(vInfo, coll, fill);

}

void I2D_RotatingAirfoil::getObstacleInfo(vector<Real> & infoObstacle)
{
	infoObstacle.clear();

	Real cor[2] = {0.0,0.0};
	airfoil->getAerodynamicCenter(cor);

	infoObstacle.push_back(cor[0]);
	infoObstacle.push_back(cor[1]);
}

void I2D_RotatingAirfoil::restart(const double t, string filename)
{
	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("restart_I2D_RotatingAirfoil.txt", "r");
		assert(ppFile!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		ppFile = fopen(filename.c_str(), "r");
		assert(ppFile!=NULL);
	}

	// Actual restart
	airfoil->restart(ppFile);

	// Close file
	fclose(ppFile);

}

void I2D_RotatingAirfoil::save(const double t, string filename)
{

	FILE * ppFile = NULL;

	if(filename==std::string())
	{
		// If the string is not set I write my own file
		ppFile = fopen("restart_I2D_RotatingAirfoil.txt", "w");
		assert(ppFile!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		ppFile = fopen(filename.c_str(), "w");
		assert(ppFile!=NULL);
	}

	// Actual save
		airfoil->save(ppFile);

	// Close file
	fclose(ppFile);

}

void I2D_RotatingAirfoil::refresh(const double t, string filename)
{
	// Open file stream
	ifstream filestream;
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		filestream.open("restart_I2D_RotatingAirfoil.txt");
		if (!filestream.good())
		{
			cout << "ooops: file not found. Exiting now." << endl;
			exit(-1);
		}
	}
	else
	{
		filestream.open(filename.c_str());
		if (!filestream.good())
		{
			cout << "ooops: file not found. Exiting now." << endl;
			exit(-1);
		}
	}

	// Open file
	int c = 0;
	string line;
	while( getline(filestream, line) ) c++;
	filestream.close();

	FILE * f = NULL;
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		f = fopen("restart_I2D_RotatingAirfoil.txt", "r");
		assert(f!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		f = fopen(filename.c_str(), "r");
		assert(f!=NULL);
	}

	//data stored in dataVect until t=t_restart
	int N=0;
	vector < vector<float> > dataVect;
	for(int i=0; i<c; i++)
	{
		vector<float> row;
		float variable = 0.0;
		fscanf(f,"%f",&variable);
		const Real dimensionlessTime = variable;
		fscanf(f,"%f",&variable);
		const Real time = variable;
		if(time <= t)
		{
			row.push_back(dimensionlessTime);
			row.push_back(time);
			fscanf(f,"%f",&variable);
			const Real xm = variable;
			row.push_back(xm);
			fscanf(f,"%f",&variable);
			const Real ym = variable;
			row.push_back(ym);
			fscanf(f,"%f",&variable);
			const Real vx = variable;
			row.push_back(vx);
			fscanf(f,"%f",&variable);
			const Real vy = variable;
			row.push_back(vy);
			fscanf(f,"%f",&variable);
			const Real angle = variable;
			row.push_back(angle);
			fscanf(f,"%f",&variable);
			const Real angular_velocity = variable;
			row.push_back(angular_velocity);
			fscanf(f,"%f",&variable);
			const Real J = variable;
			row.push_back(J);
			fscanf(f,"%f",&variable);
			const Real m = variable;
			row.push_back(m);
			fscanf(f,"%f",&variable);
			const Real rho = variable;
			row.push_back(rho);
			fscanf(f,"%f",&variable);
			const Real cD = variable;
			row.push_back(cD);
			N=i;
		}
		else
		{
			break;
		}
		dataVect.push_back(row);
	}
	fclose(f);
	f = NULL;

	//dataVect is copied in a new "shape_001" file
	if(filename==std::string())
	{
		// If the string is not set I write my own file
		f = fopen("restart_I2D_RotatingAirfoil.txt", "w");
		assert(f!=NULL);
	}
	else
	{
		// If string is set I open the corresponding file
		f = fopen(filename.c_str(), "w");
		assert(f!=NULL);
	}
	for (int i=0; i<N; i++)
	{
		fprintf(f, "%e %e %e %e %e %e %e %e %e %e \n", dataVect[i][0],dataVect[i][1],dataVect[i][2],dataVect[i][3],dataVect[i][4],dataVect[i][5],dataVect[i][6],dataVect[i][7],dataVect[i][8],dataVect[i][9],dataVect[i][10],dataVect[i][11]);
	}
	fclose(f);

}


void I2D_RotatingAirfoil::computeDesiredVelocity(const double t)

{
	const int NQUANTITIES = 3;

	double mass=0;
	double J=0;
	double omegabar=0;

	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	const BlockCollection<B>& coll = grid.getBlockCollection();

	map<int, vector<double> > integrals;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
	{
		integrals[it->blockID]= vector<double>(NQUANTITIES);
	}


	map<int, bool> nonempty;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
	{
		nonempty[it->blockID] = false;
	}

// -------------------------Compute all--------------------------------------------------------

	AirfoilStuff::ComputeAll computeAll(eps, airfoil, integrals, Uinf, nonempty);
	block_processing.process(vInfo, coll, computeAll);

// -------------------------Sum up contributions----------------------------------------------------------


	for(map<int, vector<double> > ::const_iterator it=integrals.begin(); it!=integrals.end(); ++it)
	{
			mass += (it->second)[0];
			J+= (it->second)[1];
			omegabar += (it->second)[2];
	}

// Normalization (divide by mass or moment of inertia - in the most complicated cases J varies and must be recomputed on the fly)

	  omegabar /= J;

// Set the right angular velocity for each single object
	airfoil->angular_velocity = omegabar;

		printf("\n\n");
		printf("mass=%e\n", mass);
		printf("J=%e\n", J);
		printf("anglular_velocity=%e\n\n", airfoil->angular_velocity);
		printf("\n\n");

// -----------------------Set desired velocities-----------------------------------------------
	for	(map<int, const VelocityBlock *>::iterator it = desired_velocity.begin(); it!= desired_velocity.end(); it++)
	{
		assert(it->second != NULL);
		VelocityBlock::deallocate(it->second);
	}

	desired_velocity.clear();

	vector<pair< BlockInfo, VelocityBlock *> > velblocks;
	for(vector<BlockInfo>::const_iterator it=vInfo.begin(); it!=vInfo.end(); it++)
	{
		if(nonempty[it->blockID] == true)
		{
			VelocityBlock * velblock = VelocityBlock::allocate(1);
			desired_velocity[it->blockID] = velblock;
			velblocks.push_back(pair< BlockInfo, VelocityBlock *>(*it, velblock));
		}
	}

		AirfoilStuff::FillVelblocks fillvelblocks(velblocks, eps, airfoil);
		tbb::parallel_for(blocked_range<int>(0, velblocks.size()), fillvelblocks, auto_partitioner());
}


void I2D_RotatingAirfoil::update(const double dt, const double t, string filename, map< string, vector<I2D_FloatingObstacleOperator *> > * _data)
{
				airfoil->update_all(dt,t);

				FILE * ppFile = NULL;

				if(filename==std::string())
				{
					// If the string is not set I write my own file
					ppFile = fopen("update_I2D_RotatingAirfoil.txt", t == 0.0 ? "w" : "a");
					assert(ppFile!=NULL);
				}
				else
				{
					// If string is set I open the corresponding file
					ppFile = fopen(filename.c_str(), t == 0.0 ? "w" : "a");
					assert(ppFile!=NULL);
				}

				// Write update data
				fprintf(ppFile, "%e %e %e %e %e %e %e %e %e %e %e %e\n", this->dimT, t, airfoil->xm, airfoil->ym, airfoil->theta, airfoil->angular_velocity, airfoil->J, airfoil->m, airfoil->rho, this->Cd);
				fflush(ppFile);

				// Close file
				fclose(ppFile);

}
// -------------------END OF MEMBER FUNCTIONS OF CLASS "I2D_RotatingAirfoil"--------------------------------
