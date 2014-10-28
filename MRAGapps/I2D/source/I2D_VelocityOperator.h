/*
 *  I2D_VelocityOperator.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#pragma once

class I2D_VelocityOperator
{
public:
	virtual void compute_velocity() = 0;
};