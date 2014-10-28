//
//  I2D_VelocitySolver_Wim.h
//  I2D_ROCKS
//
//  Created by Wim van Rees on 4/2/13.
//
//

#ifndef __I2D_ROCKS__I2D_VelocitySolver_Wim__
#define __I2D_ROCKS__I2D_VelocitySolver_Wim__

#pragma once

#include <limits>

#include "I2D_VelocityOperator.h"
#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_FMMTypes.h"

class I2D_VelocitySolver_Wim: public I2D_VelocityOperator
{
protected:
    
	Grid<W,B>* grid_ptr;
	BlockProcessing block_processing;
    
	Real theta;
	Real tolParticle;
	Real scaling_factor;

public:
    
	I2D_VelocitySolver_Wim(Grid<W,B>& grid, ArgumentParser& parser_): 
    grid_ptr(&grid)
	{

		theta = parser_("-fmm-theta").asDouble();
        const int lmax = parser_("-lmax").asDouble();
        
		const double min_dV = pow(pow(0.5,lmax)/B::sizeX, 2);
		scaling_factor = max(1., 1e-4/min_dV);
		tolParticle = 2*numeric_limits<Real>::epsilon();
        
		printf("scaling_factor=%f tolParticle=%f\n", scaling_factor, tolParticle);
        assert(theta>=0 && theta<1);
	}
    
	virtual void compute_velocity();
};


#endif /* defined(__I2D_ROCKS__I2D_VelocitySolver_Wim__) */
