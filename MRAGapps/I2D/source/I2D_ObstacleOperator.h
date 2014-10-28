/*
 *  I2D_ObstacleOperator.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#pragma once

class I2D_ObstacleOperator
{
public:	
	virtual void characteristic_function() = 0;
	virtual Real getD() const = 0;
	virtual void getObstacleInfo(vector<Real> & infoObstacle){}
	virtual ~I2D_ObstacleOperator(){}
};
