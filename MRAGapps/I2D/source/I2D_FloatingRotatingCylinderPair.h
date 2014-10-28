//
//  I2D_FloatingRotatingCylinderPair.h
//  I2D_ROCKS
//
//  Created by Wim van Rees on 02/12/13.
//
//

#include "I2D_FloatingObstacleOperator.h"

#ifndef __I2D_ROCKS__I2D_FloatingRotatingCylinderPair__
#define __I2D_ROCKS__I2D_FloatingRotatingCylinderPair__

class I2D_FloatingRotatingCylinderPair : public I2D_FloatingObstacleOperator
{
protected:

public:
    
	class RotatingCylinderPair
	{
        private:
		Real _mollified_heaviside(const double dist, const double eps) const;
		
        void _setTemporalPreFactor(double t);
	Real _getBurgersFactor(double tRel); // timeRel assumes a periodicity of 1
	
        public:
		double angle, xm, ym, D, width, angular_velocity, vx, vy, m, J, rho, vdefx, vdefy, Uinf[2];
		double GammaLeft, GammaRight; // strength
		double timePrefac;
		double burgersTime; // this is the time with which the sine wave has evolved according to burgers equation
        
		RotatingCylinderPair(Real xm, Real ym, Real D, Real width, Real angle_rad, Real GammaLeft, Real GammaRight, Real burgersTime=0.0);
		
		void update_all(double dt,double t);
		void restart(FILE * f);
		void save(FILE * f) const;
		Real sample(const Real x_, const Real y_, const Real eps) const;
		void bbox(const Real eps, Real xmin[2], Real xmax[2]) const;
		void restart();
		void save();
	};
	
	RotatingCylinderPair * shape;
    
    I2D_FloatingRotatingCylinderPair(ArgumentParser & parser, Grid<W,B> & grid, const Real _xm, const Real _ym, const Real D, const Real angle, const Real width, const Real GammaLeft, const Real GammaRight, const Real eps, const Real Uinf[2],I2D_PenalizationOperator& penalization, const int LMAX, const int ID = 0, RL::RL_TabularPolicy ** policy = NULL, const int seed = 0);
    
	virtual ~I2D_FloatingRotatingCylinderPair();
	
	void characteristic_function();
	Real getD() const {return D;}
	
	virtual void update(const double dt, const double t, string filename = std::string(), map< string, vector<I2D_FloatingObstacleOperator *> > * _data = NULL);
	void computeDesiredVelocity(const double t);
	void create(const double t){}
    
	virtual void save(const double t, string filename = std::string());
	virtual void restart(const double t, string filename = std::string());
	virtual void refresh(const double t, string filename = std::string());
	Real getModulusMaxVel(){
        // take total circulation
        const Real gamma = std::abs(shape->GammaLeft) + std::abs(shape->GammaRight);
        // Re will be this velocity * D / nu
        // we want Re = Gamma / nu
        // so we return Gamma / D
        return gamma/shape->D;
    }
};

#endif /* defined(__I2D_ROCKS__I2D_FloatingRotatingCylinderPair__) */
