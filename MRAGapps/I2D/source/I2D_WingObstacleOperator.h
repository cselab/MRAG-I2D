/*
 *  untitled.h
 *  I2D_ROCKS
 *
 *  Created by Mattia Gazzola on 4/1/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#pragma once

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_ObstacleOperator.h"

#include <map>

struct ChiBlock
{
	Real chi[_BLOCKSIZE_][_BLOCKSIZE_];
};

class I2D_WingObstacleOperator: public I2D_ObstacleOperator
{
	map<int, ChiBlock *> cached_blocks;
	int nBleeds, nTrailingBleeds;
	bool plainAirfoil;
	Real smoothing_length;
	Grid<W,B>& grid;
	BlockProcessing block_processing;
	
public:
	I2D_WingObstacleOperator(Grid<W,B>& grid, ArgumentParser& parser, Real cm[2], const Real smoothing_length);
	
	void getObstacleInfo(vector<Real> & infoObstacle);
	
	void characteristic_function();
	
	Real getD() const {return 2./16;}
};