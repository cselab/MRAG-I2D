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

class I2D_Clear
	{		
		BlockProcessing block_processing;
	public:
		void clearTmp(Grid<W,B> & grid);
		void clearVel(Grid<W,B> & grid);
	};
