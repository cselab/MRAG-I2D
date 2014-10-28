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
#include "I2D_FloatingObstacleOperator.h"
#include "I2D_VelocityOperator.h"
#include "I2D_PenalizationOperator.h"
#include "I2D_AdvectionOperator.h"
#include "I2D_DiffusionOperator.h"
#include "I2D_FlowPastFloatingObstacle.h"
#include "I2D_ObjectFactory.h"
#include "I2D_KillVortRightBoundaryOperator.h"

class I2D_FluidMediatedInteractions: public I2D_FlowPastFloatingObstacle
{	
protected:
	bool bUSEPOTENTIAL;
	Real charLength, charVel;
	I2D_ObjectFactory * factory;
	I2D_VelocityOperator * potsolver;
	bool bUSEKILLVORT;
	int KILLVORT;
	I2D_KillVortRightBoundaryOperator * killVort;
	bool bUSEOPTIMIZER;
	Real TBOUND;
	
public:	
	I2D_FluidMediatedInteractions(const int argc, const char ** argv);
	~I2D_FluidMediatedInteractions();
	virtual void run();
};
