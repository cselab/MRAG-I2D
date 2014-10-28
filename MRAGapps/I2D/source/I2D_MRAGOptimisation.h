/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Chloe Mimeau on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"
#include "I2D_VelocityOperator.h"
#include "I2D_PenalizationOperator.h"
#include "I2D_AdvectionOperator.h"
#include "I2D_DiffusionOperator.h"
#include "I2D_KillVortRightBoundaryOperator.h"
#include "I2D_FlowPastFixedObstacle.h"


class I2D_MRAGOptimisation: public I2D_FlowPastFixedObstacle
{
	int KILLVORT;
	string sCTRL;
	
	I2D_KillVortRightBoundaryOperator * killVort;
	
public:
	
	I2D_MRAGOptimisation(const int argc, const char ** argv);
	
	void run();
	void paint(){}
};
