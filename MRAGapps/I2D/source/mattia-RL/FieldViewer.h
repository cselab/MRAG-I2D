/*
 *  FieldViewer.h
 *  Exercise 9
 *
 *  Created by Diego Rossinelli on 4/26/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "RL_Environment.h"
#include <assert.h>

#ifdef _RL_VIZ
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <vector>
using namespace std;

//DO NOT LOOK AT ME - i barely exist

namespace RL
{

class Viewer
{
protected:
	int size[2];
	float * texture_data;
	GLuint texture_handle;
	
public:
	
	virtual void view() = 0;
	
	Viewer(): texture_data(NULL), texture_handle(0)
	{
		size[0] = size[1] = 0;
	}

	void allocate_texture(int sizeX, int sizeY)
	{
		assert(texture_handle == 0);
		
		texture_data = new float[sizeX*sizeY];
		memset(texture_data, 0, sizeof(float)*sizeX*sizeY);
		
		size[0] = sizeX;
		size[1] = sizeY;
		
		glGenTextures(1, &texture_handle);
		assert(texture_handle != 0);
		
		glBindTexture(GL_TEXTURE_2D, texture_handle);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, size[0], size[1], 0, GL_RED, GL_FLOAT, 0);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
	
	void upload_texture()
	{
		assert(texture_handle != 0);
		
		glBindTexture(GL_TEXTURE_2D, texture_handle);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size[0], size[1], GL_LUMINANCE, GL_FLOAT, texture_data);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
};

class FieldViewer
{
public:
	vector<Viewer *> views;
	
	void register_view(Viewer * v)
	{
		assert(v!=NULL);
		views.push_back(v);
	}
	
	void view()
	{
		if (views.size() == 0) return;
		
		for(int i=0; i<(int)views.size(); i++)
			views[i]->view();
	}
};

}

#endif
