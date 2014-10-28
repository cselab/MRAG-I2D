//
//  I2D_RotatingCylinderPair.h
//  I2D_ROCKS
//
//  Created by Wim van Rees on 05/12/13.
//
//

#ifndef __I2D_ROCKS__I2D_RotatingCylinderPair__
#define __I2D_ROCKS__I2D_RotatingCylinderPair__

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

class I2D_RotatingCylinderPair: public I2D_FlowPastFloatingObstacle
{
    protected:
	Real charLength, charVel;
	I2D_ObjectFactory * factory;
	bool bUSEKILLVORT;
	int KILLVORT;
	I2D_KillVortRightBoundaryOperator * killVort;
	bool bUSEOPTIMIZER;
	Real TBOUND;

	Real TIMESCALE;
    Real LENGTHSCALE;
    Real CIRCULATION;
    Real OMEGA; // angular velocity rad/s
    public:
	I2D_RotatingCylinderPair(const int argc, const char ** argv);
	~I2D_RotatingCylinderPair();
	virtual void run();
	void _tnext(double &tnext, double& tnext_dump, double& tend); // overriding FlowPastFixedObstacle::_tnext
};


#endif /* defined(__I2D_ROCKS__I2D_RotatingCylinderPair__) */
