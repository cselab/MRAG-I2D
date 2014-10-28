//
//  I2D_TestMultipole.h
//  I2D_ROCKS
//
//  Created by Wim van Rees on 4/1/13.
//
//

#ifndef __I2D_ROCKS__I2D_TestMultipole__
#define __I2D_ROCKS__I2D_TestMultipole__

#include <cmath>
#include <limits>

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_VelocityOperator.h"

struct LambOseenVortex
{
private:
    Real Gamma;
    Real org[2];
    
public:
    LambOseenVortex(const Real _Gamma, const Real _org[2]):Gamma(_Gamma)
    {
        org[0]=_org[0];
        org[1]=_org[1];
    }
    
    void gimmeVort(const Real x[2], const Real nut, Real & w)
    {
        const Real relx = x[0]-org[0];
        const Real rely = x[1]-org[1];
        const Real rsq = relx*relx + rely*rely;

        const Real fac = Gamma/(4.0*M_PI*nut);
        w = fac*std::exp(-rsq/(4.0*nut));
    }
    
    void gimmeVel(const Real x[2], const Real nut, Real u[2])
    {
        const Real relx = x[0]-org[0];
        const Real rely = x[1]-org[1];
        const Real rsq = relx*relx + rely*rely;
        const Real r = std::sqrt(rsq);
        
        if(rsq<std::numeric_limits<Real>::epsilon())
        {
            u[0]=u[1]=0;
            return;
        }
        const Real fac = Gamma/(2.0*M_PI*r);
        const Real utheta = fac*(1.0 - std::exp(-rsq/(4.0*nut)));
        
        u[0] = -utheta*rely/r;
        u[1] =  utheta*relx/r;
    }
};

struct ExactSolutionBlock
{
    Real u_exact[_BLOCKSIZE_][_BLOCKSIZE_][2];
};


class I2D_TestMultipole: public I2D_Test
{
	ArgumentParser parser;
	
	Grid<W,B> * grid;
	Refiner * refiner;
	Compressor * compressor;
	
	BlockFWT<W, B, vorticity_projector, false, 1> fwt_omega;
	BlockFWT<W, B, velocity_projector, false, 2> fwt_velocity;
    
	set<int> _getBoundaryBlockIDs();
	void _dump(string filename);
	void _ic_omega(Grid<W,B>& grid);
	void _computeError();
	void _refine(bool bUseIC);
	void _compress(bool bUseIC);
    
	void _restart();
	void _save();
    
    void _storeExactSolution();
    std::vector<ExactSolutionBlock> exactsol;
    
	Real t;
	int step_id;
    bool bRestart;
    
	I2D_VelocityOperator * poisson_solver;
    
	Profiler profiler;
    
    Real LambOseenGamma,LambOseenNut, LambOseenOrg[2];
public:
	
	I2D_TestMultipole(const int argc, const char ** argv);
	
	void run();
	void paint();
};


#endif /* defined(__I2D_ROCKS__I2D_TestMultipole__) */
