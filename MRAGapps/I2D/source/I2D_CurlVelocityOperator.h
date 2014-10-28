/*
 *  I2D_CurlVelocityOperator.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"

class I2D_CurlVelocityOperator
{
protected:
	Grid<W,B>& grid;
	BlockProcessing block_processing;
	
public:
	
	I2D_CurlVelocityOperator(Grid<W,B>& grid): grid(grid){}
	
	virtual void perform()=0;
};


class I2D_CurlVelocityOperator_2ndOrder: public I2D_CurlVelocityOperator
{	
public:
	
	I2D_CurlVelocityOperator_2ndOrder(Grid<W,B>& grid): I2D_CurlVelocityOperator(grid){}
	
	void perform();
};


class I2D_CurlVelocityOperator_4thOrder: public I2D_CurlVelocityOperator
{	
public:
	
	I2D_CurlVelocityOperator_4thOrder(Grid<W,B>& grid): I2D_CurlVelocityOperator(grid){}
	
	void perform();
};