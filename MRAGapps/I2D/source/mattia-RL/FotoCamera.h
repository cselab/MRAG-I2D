/*
 *  FotoCamera.h
 *  Swimmers
 *
 *  Created by menahel on 7/2/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#pragma once

#include "RL_Environment.h"

#ifdef _RL_VIZ

#include <stdlib.h>
#include <stdio.h>
#include "Screenshot.h"

class FotoCamera
{
public:
	// Constructor
	FotoCamera(){};
	~FotoCamera(){};
	
	// Methods
	void shoot(void)
	{
		int frequency = 1;
		static int iFrameCounterAbsolute = 0;
		static int iFrameCounter = 0;
		char buf[300];
		if( iFrameCounterAbsolute % frequency == 0 )
		{
			iFrameCounter++;
			//if(iFrameCounter < 7000)
			//{
			sprintf(buf, "img%07d.tga", iFrameCounter);
			printf("Writing to %s\n", buf);
			gltWriteTGA(buf);
			//}
		}
		iFrameCounterAbsolute++;
	}
};

#endif
