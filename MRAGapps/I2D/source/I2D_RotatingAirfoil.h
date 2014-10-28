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
#include "I2D_PenalizationOperator.h"

class I2D_RotatingAirfoil: public I2D_FloatingObstacleOperator
{
public:

	class NACA4gen
	{
			private:
				unsigned int nPanels;
				int maxDiffX, maxDiffEpsilon;
				vector<Real> x,yt,xu,yu,xl,yl,yC;
				// maxBoxX: stores physical coordinate of maximal base point in x; the others: analogue
				Real chord, maxBoxX, minBoxX, maxBoxY, minBoxY,leadingEdgeRadius,epsilon;
				Real leadingEdgeCenter[2];
				const int d1,d2,d3d4; // NACA digits: NACA d1-d2-d3d4

			public:
				Real getMaxX(void) const { return maxBoxX; }
				Real getMinX(void) const { return minBoxX; }
				Real getMaxY(void) const { return maxBoxY; }
				Real getMinY(void) const { return minBoxY; }
				Real getLeadingEdgeRadius(void){ return leadingEdgeRadius; }
				Real getNoseDistance(Real xx[2]) const;
				Real getYUpperProfile(Real xx[2], Real & dist) const;
				Real getYLowerProfile(Real xx[2], Real & dist) const;
				Real sdf(Real xo[2]) const;
				NACA4gen(Real chord_in,unsigned int nPanels_in,bool isFiniteTE,Real epsilon_in, const int _d1, const int _d2, const int _d3d4);
				~NACA4gen();
	};


	class DiscretizedWing
	{
	        private:
		        Real O2W[3][3], W2O[3][3]; // W2O = inv(O2W)
		        Real ca, sa;



	        public:
		        NACA4gen naca;
		        const Real scaling, width, epsilon, s;
	        	Real tx, ty, angle, self_angle,angular_velocity, vx, vy, rho, m, J, l;

	        	DiscretizedWing (Real _tx,      Real _ty,         Real width,   Real epsilon,
	        			         Real _angle,   Real _self_angle, const int d1, const int d2,
	        			         const int d3d4, const Real l);


		        void check_mapping() const;
		        void _w2o(const Real xw[2], Real xo[2]) const;
				void _o2w(const Real xo[2], Real xw[2]) const;
	        	void update_all(double dt, double t, const double peristrofiki, const double _D, const double _theta);
		        void getAerodynamicCenter(Real cor[2]) const; // calculates position and orientation of airfoil
				Real sdf(Real xw[2]) const;
				void bbox(const Real eps, Real xmin[2], Real xmax[2]) const;
				void restart(FILE * f);
		        void save(FILE * f) const;
		        Real _mollified_heaviside(const Real x, const Real eps) const;
		        Real sample(Real p[2], const Real eps) const; // calculate characteristic function at point (p[0],p[1])


	};


	class ToomanyAirfoils
	{
	    public:

		     DiscretizedWing *wing_1, *wing_2, *wing_3;
		     Real xm, ym,theta, D, angular_velocity, rho, m, J, epsilon;
		     ToomanyAirfoils(   Real _xm ,     Real _ym,    const Real _D,  const Real _c,
		    		            const int d1,  const int d2,const int d3d4, const Real _l,
		    		            const Real epsilon);

		     Real _mollified_heaviside(const Real x, const Real eps) const;
		     Real sample(Real p[2], const Real eps) const;
		     void getAerodynamicCenter(Real cor[2]) const; // calculates position and orientation of airfoil
		     void bbox(const Real eps, Real xmin[2], Real xmax[2]) const;
		     void rotate(Real v[2], Real _theta) const;
		     void restart(FILE * f);
		     void save(FILE * f) const;
		     void update_all(double dt, double t);
	};


    ToomanyAirfoils * airfoil;
	I2D_RotatingAirfoil(Grid<W,B>& grid, ArgumentParser& parser, const Real _xm,
			const Real _ym, const Real _D, const Real _c,  const int d1, const int d2,
			const int d3d4, const Real _l, const Real _eps,const Real Uinf[2],
			I2D_PenalizationOperator& penalization);

	~I2D_RotatingAirfoil();

	void characteristic_function();
	Real getD() const {return D ;}

	void getObstacleInfo(vector<Real> & infoObstacle);
	void computeDesiredVelocity(const double t);
	void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	void create(const double t){}

	void save(const double t, string filename = std::string());
	void restart(const double t, string filename = std::string());
	void refresh(const double t, string filename = std::string());

};
