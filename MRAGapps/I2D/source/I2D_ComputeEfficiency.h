/*
 * I2D_ComputeEfficiency.h
 *
 *  Created on: Mar 8, 2012
 *      Author: mgazzola
 */

#ifndef I2D_COMPUTEEFFICIENCY_H_
#define I2D_COMPUTEEFFICIENCY_H_

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_FloatingObstacleOperator.h"

class I2D_ComputeEfficiency
{
protected:
	double viscosity, density;
	Grid<W,B>& grid;
	BlockProcessing block_processing;
	I2D_FloatingObstacleOperator * floatingObstacle;

public:
	I2D_ComputeEfficiency(Grid<W,B>& grid, I2D_FloatingObstacleOperator * floatingObstacle, double viscosity, double density = 1.0): grid(grid), floatingObstacle(floatingObstacle), viscosity(viscosity), density(density)
	{
		assert(viscosity>0);
		assert(density>0);
		assert(floatingObstacle!=NULL);
	}

	virtual ~I2D_ComputeEfficiency(){};

	void compute(I2D_FloatingObstacleOperator * floatingObstacle, double t, double dt, double startTime, double endTime);
};

#endif /* I2D_COMPUTEEFFICIENCY_H_ */
