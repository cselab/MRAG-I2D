/*
 *  I2D_FlowPastFloatingObstacle.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_FloatingObstacleOperator.h"
#include "I2D_VelocityOperator.h"
#include "I2D_PenalizationOperator.h"
#include "I2D_AdvectionOperator.h"
#include "I2D_DiffusionOperator.h"
#include "I2D_FlowPastFixedObstacle.h"

class I2D_FlowPastFloatingObstacle: public I2D_FlowPastFixedObstacle
{	
protected:
	I2D_FloatingObstacleOperator * floatingObstacle;
	
	void _restart();
	void _save();
	void _refresh(const Real t);
	
public:	

	I2D_FlowPastFloatingObstacle(const int argc, const char ** argv);
	~I2D_FlowPastFloatingObstacle();
	
	virtual void run();
};
