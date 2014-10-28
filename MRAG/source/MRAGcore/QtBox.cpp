/*
 *  QtBox.cpp
 *  Multipole
 *
 *  Created by Mattia Gazzola on 9/27/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 *
 *  THE QUADRANT STRUCTURE ASSUMES THE FOLLOWING INDEXING!!
 *
 *			-------------
 *			|     |     |
 *			|  2  |  3  |
 *			|     |     |
 *			-------------
 *			|     |     |
 *			|  0  |  1  |
 *			|     |     |
 *			-------------
 *
 */

#include <cstdio>
#include "QtBox.h"
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include "assert.h"

QtBox::QtBox(const unsigned int _maxLevel, const double _width, const double _minVertex[2], const int _invLevel, QtBox * _parent, const unsigned int _idxX, const unsigned int _idxY):
MAX_N_LEVELS(_maxLevel), width(_width), invLevel(_invLevel), parent(_parent), idxX(_idxX), idxY(_idxY), bEmpty(true)
{
	//printf("idxX=%d, idxY=%e\n",idxX,idxY);

	assert(MAX_N_LEVELS>=1);

	ROOT_LEVEL = MAX_N_LEVELS-1;
	MAX_VAL = pow(2,ROOT_LEVEL);

	level = _invLevelToLevel(invLevel);

	minVertex[0] = _minVertex[0];
	minVertex[1] = _minVertex[1];

	xLocCode = _getLocationalCode(minVertex[0]);
	yLocCode = _getLocationalCode(minVertex[1]);

	for(int i=0;i<4;i++)
		children[i] = NULL;
}

QtBox::~QtBox()
{
	if(!bEmpty)
		for(int i=0;i<4;i++)
			if(children[i]!=NULL)
			{
				delete children[i];
				children[i] = NULL;
			}
			else
			{
				printf("Ma come cazzo e' possibile!\n");
				abort();
			}
}

I3 QtBox::getIndices()
{
	return I3(idxX,idxY,invLevel);
}

I3 QtBox::getIndicesTranslated(QtBox * cell)
{
	{
		const int maxlevel = max(cell->invLevel,invLevel);
		assert(maxlevel>1);
		const int shift = 1 << invLevel;
		const int deltaL = invLevel - cell->invLevel;

		int xg, yg, xt, yt;

		if (deltaL >= 0)
		{
			xg = idxX >> deltaL;
			yg = idxY >> deltaL;

			xt = cell->idxX;
			yt = cell->idxY;
		}
		else
		{
			xg = idxX ;
			yg = idxY ;

			xt = cell->idxX >> -deltaL;
			yt = cell->idxY >> -deltaL;
		}

		const int dx = xg - xt;
		const int dy = yg - yt;

		int fx = 0;
		int fy = 0;

		if( cell->invLevel == 1 || invLevel==1 )
		{
			const double dxCell = 1.0 / (double)(1 << cell->invLevel);
			const double dxGhost = 1.0 / (double)(1 << invLevel);

			const double minCellX = (double)cell->idxX*dxCell;
			const double maxCellX = (double)(cell->idxX+1)*dxCell;
			const double minCellY = (double)cell->idxY*dxCell;
			const double maxCellY = (double)(cell->idxY+1)*dxCell;

			const double minGhostX = (double)idxX*dxGhost;
			const double maxGhostX = (double)(idxX+1)*dxGhost;
			const double minGhostY = (double)idxY*dxGhost;
			const double maxGhostY = (double)(idxY+1)*dxGhost;

			assert(minCellX>=0.0 && minCellX<=1.0);
			assert(maxCellX>=0.0 && maxCellX<=1.0);
			assert(minCellY>=0.0 && minCellY<=1.0);
			assert(maxCellY>=0.0 && maxCellY<=1.0);
			assert(minGhostX>=0.0 && minGhostX<=1.0);
			assert(maxGhostX>=0.0 && maxGhostX<=1.0);
			assert(minGhostY>=0.0 && minGhostY<=1.0);
			assert(maxGhostY>=0.0 && maxGhostY<=1.0);

			const double minX = min(minCellX,minGhostX);
			const double minY = min(minCellY,minGhostY);

			const double maxX = max(maxCellX,maxGhostX);
			const double maxY = max(maxCellY,maxGhostY);

			if(!((maxX-minX)==(dxCell+dxGhost)))
				fx = dx <= -1 ? +1 : dx >= 1 ? -1 : 0;
			else
				fx = dx < -1 ? +1 : dx > 1 ? -1 : 0;

			if(!((maxY-minY)==(dxCell+dxGhost)))
				fy = dy <= -1 ? +1 : dy >= 1 ? -1 : 0;
			else
				fy = dy < -1 ? +1 : dy > 1 ? -1 : 0;
		}
		else
		{
			fx = dx < -1 ? +1 : dx > 1 ? -1 : 0;
			fy = dy < -1 ? +1 : dy > 1 ? -1 : 0;
		}

		const int newX = idxX + fx*shift;
		const int newY = idxY + fy*shift;

		return I3(newX, newY, invLevel);
	}

	printf("If you are here something is very wrong!\n");
	abort();
	return I3(0,0,0);
}

void QtBox::split(std::vector<QtBox*> & leaves)
{
	// Reached max depth, set box to leaf and return!
	if(level<=0)
	{
		bEmpty = true;
		leaves.push_back(this);
		return;
	}

	// Constant levels and indices
	const unsigned int pix = idxX;
	const unsigned int piy = idxY;

	if(_splitCondition())
	{
		// Create the children boxes and recursively split again
		bEmpty = false;
		double childMinVertex[2] = {0.0,0.0};
		double childWidth;

		for(unsigned int i=0; i<2; i++)
			for(unsigned int j=0; j<2; j++)
			{
				const unsigned int coordIdx[2] = {i,j};
				_getChildBoxDimensions(minVertex,width,childMinVertex,childWidth,coordIdx);
				const unsigned int idx = i+2*j;
				children[idx] = new QtBox(MAX_N_LEVELS,childWidth, childMinVertex,invLevel+1,this,2*pix+i,2*piy+j);
				children[idx]->split(leaves);
			}
	}
	else
	{
		bEmpty = true;
		leaves.push_back(this);
		return;
	}
}

void QtBox::split(I3 input, std::vector<QtBox*> & leaves)
{
	// Reached max depth, set box to leaf and return!
	if(level<=0)
	{
		bEmpty = true;
		leaves.push_back(this);
		return;
	}

	// Constant levels and indices
	const unsigned int ix = input.i[0];
	const unsigned int iy = input.i[1];
	const unsigned int l = input.i[2];
	const unsigned int pix = idxX;
	const unsigned int piy = idxY;
	const unsigned int pl = invLevel;

	// Reached target depth, set box to leaf and return!
	if(pl==l)
	{
		bEmpty = true;
		leaves.push_back(this);
		return;
	}

	{
		if(bEmpty)
		{
			// Create the children boxes
			bEmpty = false;
			double childMinVertex[2] = {0.0,0.0};
			double childWidth;
			for(unsigned int i=0; i<2; i++)
				for(unsigned int j=0; j<2; j++)
				{
					const unsigned int coordIdx[2] = {i,j};
					_getChildBoxDimensions(minVertex,width,childMinVertex,childWidth,coordIdx);
					children[i+2*j] = new QtBox(MAX_N_LEVELS,childWidth, childMinVertex,invLevel+1,this,2*pix+i,2*piy+j);
				}
		}

		unsigned int selectX = 1;
		unsigned int selectY = 1;
		for(unsigned int i=0; i<2; i++)
			for(unsigned int j=0; j<2; j++)
			{
				const unsigned int newix = 2*pix + i;
				const unsigned int newiy = 2*piy + j;
				const unsigned int case1 = pow( 2,(l-(pl+1)) ) * newix;
				const unsigned int case2 = pow( 2,(l-(pl+1)) ) * newiy;
				selectX *= (ix>=case1);
				selectY *= (iy>=case2);
			}
		children[selectX+2*selectY]->split(input,leaves);
	}
}

bool QtBox::isInside(const double p[2])
{
	const double x = p[0];
	const double y = p[1];
	const double minx = minVertex[0];
	const double maxx = minx + width;
	const double miny = minVertex[1];
	const double maxy = miny + width;

	return (x>=minx) && (x<maxx) && (y>=miny) && (y<maxy);
}

bool QtBox::isMyLeftNeighbor(const QtBox * leftNeigh)
{
	double intpart;
	const double myMinx = minVertex[0];
	const double leftNeighMaxx = modf( (leftNeigh->minVertex[0]+leftNeigh->width), &intpart);
	const bool touchLeftSide = (myMinx==leftNeighMaxx);

	// y-direction
	const double myMiny = minVertex[1];
	const double myMaxy = minVertex[1] + width;
	const double myLength = fabs(myMaxy-myMiny);
	const double leftNeighMiny = leftNeigh->minVertex[1];
	const double leftNeighMaxy = leftNeigh->minVertex[1]+leftNeigh->width;
	const double leftNeigLength = fabs(leftNeighMaxy-leftNeighMiny);
	const double length = myLength + leftNeigLength;
	const double maxy = max(myMaxy,leftNeighMaxy);
	const double miny = min(myMiny,leftNeighMiny);
	const double maxDist = fabs(maxy-miny);
	const bool withinHeight = (maxDist<length);

	if( !(touchLeftSide && withinHeight) ){ printf("Not LEFT neighbor\n"); exit(0); }

	return touchLeftSide && withinHeight;
}

bool QtBox::isMyRightNeighbor(const QtBox * rightNeigh)
{
	double intpart;
	const double myMaxx = modf( (minVertex[0]+width), &intpart);
	const double rightNeighMinx = rightNeigh->minVertex[0];
	const bool touchRightSide = (myMaxx==rightNeighMinx);

	// y-direction
	const double myMiny = minVertex[1];
	const double myMaxy = minVertex[1] + width;
	const double myLength = fabs(myMaxy-myMiny);
	const double rightNeighMiny = rightNeigh->minVertex[1];
	const double rightNeighMaxy = rightNeigh->minVertex[1]+rightNeigh->width;
	const double rightNeigLength = fabs(rightNeighMaxy-rightNeighMiny);
	const double length = myLength + rightNeigLength;
	const double maxy = max(myMaxy,rightNeighMaxy);
	const double miny = min(myMiny,rightNeighMiny);
	const double maxDist = fabs(maxy-miny);
	const bool withinHeight = (maxDist<length);

	if( !(touchRightSide && withinHeight) ){ printf("Not RIGHT neighbor\n"); exit(0); }

	return touchRightSide && withinHeight;
}

bool QtBox::isMyBottomNeighbor(const QtBox * bottomNeigh)
{
	double intpart;
	const double myMiny = minVertex[1];
	const double bottomNeighMaxy = modf( (bottomNeigh->minVertex[1]+bottomNeigh->width), &intpart);
	const bool touchBottomSide = (myMiny==bottomNeighMaxy);

	// x-direction
	const double myMinx = minVertex[0];
	const double myMaxx = minVertex[0] + width;
	const double myLength = fabs(myMaxx-myMinx);
	const double bottomNeighMinx = bottomNeigh->minVertex[0];
	const double bottomNeighMaxx = bottomNeigh->minVertex[0]+bottomNeigh->width;
	const double bottomNeigLength = fabs(bottomNeighMaxx-bottomNeighMinx);
	const double length = myLength + bottomNeigLength;
	const double maxx = max(myMaxx,bottomNeighMaxx);
	const double minx = min(myMinx,bottomNeighMinx);
	const double maxDist = fabs(maxx-minx);
	const bool withinWidth = (maxDist<length);

	if( !(touchBottomSide && withinWidth) ){ printf("Not BOTTOM neighbor\n"); exit(0); }

	return touchBottomSide && withinWidth;
}

bool QtBox::isMyTopNeighbor(const QtBox * topNeigh)
{
	double intpart;
	const double myMaxy = modf( (minVertex[1]+width), &intpart);
	const double topNeighMiny = topNeigh->minVertex[1];
	const bool touchTopSide = (myMaxy==topNeighMiny);

	// x-direction
	const double myMinx = minVertex[0];
	const double myMaxx = minVertex[0] + width;
	const double myLength = fabs(myMaxx-myMinx);
	const double topNeighMinx = topNeigh->minVertex[0];
	const double topNeighMaxx = topNeigh->minVertex[0]+topNeigh->width;
	const double topNeigLength = fabs(topNeighMaxx-topNeighMinx);
	const double length = myLength + topNeigLength;
	const double maxx = max(myMaxx,topNeighMaxx);
	const double minx = min(myMinx,topNeighMinx);
	const double maxDist = fabs(maxx-minx);
	const bool withinWidth = (maxDist<length);

	if( !(touchTopSide && withinWidth) ){ printf("Not TOP neighbor\n"); exit(0); }

	return touchTopSide && withinWidth;
}

bool QtBox::isTopLeftCornerNeighbor(const QtBox * cornerNeigh)
{
	// y-direction
	const double myMiny = minVertex[1];
	const double myMaxy = minVertex[1] + width;
	const bool periodicY = ((minVertex[1]+width)==1.0);
	const double offsetY = periodicY?1.0:0.0;
	const double myLengthY = fabs(myMaxy-myMiny);
	const double cornerNeighMiny = periodicY?(cornerNeigh->minVertex[1]+offsetY):(cornerNeigh->minVertex[1]);
	const double cornerNeighMaxy = periodicY?(cornerNeigh->minVertex[1]+cornerNeigh->width+offsetY):(cornerNeigh->minVertex[1]+cornerNeigh->width);
	const double cornerNeigLengthY = fabs(cornerNeighMaxy-cornerNeighMiny);
	const double lengthY = myLengthY + cornerNeigLengthY;
	const double maxy = max(myMaxy,cornerNeighMaxy);
	const double miny = min(myMiny,cornerNeighMiny);
	const double maxDistY = fabs(maxy-miny);
	const bool withinHeight = (maxDistY==lengthY) && (cornerNeighMiny>myMiny);

	// x-direction
	const double myMinx = minVertex[0];
	const double myMaxx = minVertex[0] + width;
	const bool periodicX = (minVertex[0]==0.0);
	const double offsetX = periodicX?-1.0:0.0;
	const double myLengthX = fabs(myMaxx-myMinx);
	const double cornerNeighMinx = periodicX?(cornerNeigh->minVertex[0]+offsetX):(cornerNeigh->minVertex[0]);
	const double cornerNeighMaxx = periodicX?(cornerNeigh->minVertex[0]+cornerNeigh->width+offsetX):(cornerNeigh->minVertex[0]+cornerNeigh->width);
	const double cornerNeigLengthX = fabs(cornerNeighMaxx-cornerNeighMinx);
	const double lengthX = myLengthX + cornerNeigLengthX;
	const double maxx = max(myMaxx,cornerNeighMaxx);
	const double minx = min(myMinx,cornerNeighMinx);
	const double maxDistX = fabs(maxx-minx);
	const bool withinWidth = (maxDistX==lengthX) && (cornerNeighMinx<myMinx);

	if( !(withinHeight && withinWidth) ){ printf("Not TOP-LEFT CORNER neighbor\n"); exit(0); }

	return withinHeight && withinWidth;
}

bool QtBox::isTopRightCornerNeighbor(const QtBox * cornerNeigh)
{
	// y-direction
	const double myMiny = minVertex[1];
	const double myMaxy = minVertex[1] + width;
	const bool periodicY = ((minVertex[1]+width)==1.0);
	const double offsetY = periodicY?1.0:0.0;
	const double myLengthY = fabs(myMaxy-myMiny);
	const double cornerNeighMiny = periodicY?(cornerNeigh->minVertex[1]+offsetY):(cornerNeigh->minVertex[1]);
	const double cornerNeighMaxy = periodicY?(cornerNeigh->minVertex[1]+cornerNeigh->width+offsetY):(cornerNeigh->minVertex[1]+cornerNeigh->width);
	const double cornerNeigLengthY = fabs(cornerNeighMaxy-cornerNeighMiny);
	const double lengthY = myLengthY + cornerNeigLengthY;
	const double maxy = max(myMaxy,cornerNeighMaxy);
	const double miny = min(myMiny,cornerNeighMiny);
	const double maxDistY = fabs(maxy-miny);
	const bool withinHeight = (maxDistY==lengthY) && (cornerNeighMiny>myMiny);

	// x-direction
	const double myMinx = minVertex[0];
	const double myMaxx = minVertex[0] + width;
	const bool periodicX = (minVertex[0]+width==1.0);
	const double offsetX = periodicX?1.0:0.0;
	const double myLengthX = fabs(myMaxx-myMinx);
	const double cornerNeighMinx = periodicX?(cornerNeigh->minVertex[0]+offsetX):(cornerNeigh->minVertex[0]);
	const double cornerNeighMaxx = periodicX?(cornerNeigh->minVertex[0]+cornerNeigh->width+offsetX):(cornerNeigh->minVertex[0]+cornerNeigh->width);
	const double cornerNeigLengthX = fabs(cornerNeighMaxx-cornerNeighMinx);
	const double lengthX = myLengthX + cornerNeigLengthX;
	const double maxx = max(myMaxx,cornerNeighMaxx);
	const double minx = min(myMinx,cornerNeighMinx);
	const double maxDistX = fabs(maxx-minx);
	const bool withinWidth = (maxDistX==lengthX) && (cornerNeighMinx>myMinx);

	if( !(withinHeight && withinWidth) ){ printf("Not TOP-RIGHT CORNER neighbor\n"); exit(0); }

	return withinHeight && withinWidth;
}


bool QtBox::isBottomLeftCornerNeighbor(const QtBox * cornerNeigh)
{
	// y-direction
	const double myMiny = minVertex[1];
	const double myMaxy = minVertex[1] + width;
	const bool periodicY = (minVertex[1]==0.0);
	const double offsetY = periodicY?-1.0:0.0;
	const double myLengthY = fabs(myMaxy-myMiny);
	const double cornerNeighMiny = periodicY?(cornerNeigh->minVertex[1]+offsetY):(cornerNeigh->minVertex[1]);
	const double cornerNeighMaxy = periodicY?(cornerNeigh->minVertex[1]+cornerNeigh->width+offsetY):(cornerNeigh->minVertex[1]+cornerNeigh->width);
	const double cornerNeigLengthY = fabs(cornerNeighMaxy-cornerNeighMiny);
	const double lengthY = myLengthY + cornerNeigLengthY;
	const double maxy = max(myMaxy,cornerNeighMaxy);
	const double miny = min(myMiny,cornerNeighMiny);
	const double maxDistY = fabs(maxy-miny);
	const bool withinHeight = (maxDistY==lengthY) && (cornerNeighMiny<myMiny);

	// x-direction
	const double myMinx = minVertex[0];
	const double myMaxx = minVertex[0] + width;
	const bool periodicX = (minVertex[0]==0.0);
	const double offsetX = periodicX?-1.0:0.0;
	const double myLengthX = fabs(myMaxx-myMinx);
	const double cornerNeighMinx = periodicX?(cornerNeigh->minVertex[0]+offsetX):(cornerNeigh->minVertex[0]);
	const double cornerNeighMaxx = periodicX?(cornerNeigh->minVertex[0]+cornerNeigh->width+offsetX):(cornerNeigh->minVertex[0]+cornerNeigh->width);
	const double cornerNeigLengthX = fabs(cornerNeighMaxx-cornerNeighMinx);
	const double lengthX = myLengthX + cornerNeigLengthX;
	const double maxx = max(myMaxx,cornerNeighMaxx);
	const double minx = min(myMinx,cornerNeighMinx);
	const double maxDistX = fabs(maxx-minx);
	const bool withinWidth = (maxDistX==lengthX) && (cornerNeighMinx<myMinx);

	if( !(withinHeight && withinWidth) ){ printf("Not BOTTOM-LEFT CORNER neighbor\n"); exit(0); }

	return withinHeight && withinWidth;
}

bool QtBox::isBottomRightCornerNeighbor(const QtBox * cornerNeigh)
{
	// y-direction
	const double myMiny = minVertex[1];
	const double myMaxy = minVertex[1] + width;
	const bool periodicY = (minVertex[1]==0.0);
	const double offsetY = periodicY?-1.0:0.0;
	const double myLengthY = fabs(myMaxy-myMiny);
	const double cornerNeighMiny = periodicY?(cornerNeigh->minVertex[1]+offsetY):(cornerNeigh->minVertex[1]);
	const double cornerNeighMaxy = periodicY?(cornerNeigh->minVertex[1]+cornerNeigh->width+offsetY):(cornerNeigh->minVertex[1]+cornerNeigh->width);
	const double cornerNeigLengthY = fabs(cornerNeighMaxy-cornerNeighMiny);
	const double lengthY = myLengthY + cornerNeigLengthY;
	const double maxy = max(myMaxy,cornerNeighMaxy);
	const double miny = min(myMiny,cornerNeighMiny);
	const double maxDistY = fabs(maxy-miny);
	const bool withinHeight = (maxDistY==lengthY) && (cornerNeighMiny<myMiny);

	// x-direction
	const double myMinx = minVertex[0];
	const double myMaxx = minVertex[0] + width;
	const bool periodicX = (minVertex[0]+width==1.0);
	const double offsetX = periodicX?1.0:0.0;
	const double myLengthX = fabs(myMaxx-myMinx);
	const double cornerNeighMinx = periodicX?(cornerNeigh->minVertex[0]+offsetX):(cornerNeigh->minVertex[0]);
	const double cornerNeighMaxx = periodicX?(cornerNeigh->minVertex[0]+cornerNeigh->width+offsetX):(cornerNeigh->minVertex[0]+cornerNeigh->width);
	const double cornerNeigLengthX = fabs(cornerNeighMaxx-cornerNeighMinx);
	const double lengthX = myLengthX + cornerNeigLengthX;
	const double maxx = max(myMaxx,cornerNeighMaxx);
	const double minx = min(myMinx,cornerNeighMinx);
	const double maxDistX = fabs(maxx-minx);
	const bool withinWidth = (maxDistX==lengthX) && (cornerNeighMinx>myMinx);

	if( !(withinHeight && withinWidth) ){ printf("Not BOTTOM-RIGHT CORNER neighbor\n"); exit(0); }

	return withinHeight && withinWidth;
}

unsigned int QtBox::spotQuadrant(const double p[2])
{
	const unsigned int idxX = floor((p[0]-minVertex[0])/(width/2.0));
	const unsigned int idxY = floor((p[1]-minVertex[1])/(width/2.0));
	return idxX + 2*idxY;
}

void QtBox::twoChildren(unsigned int q0, unsigned int q1, std::vector< QtBox * > & two)
{
	if(bEmpty)
	{
		two.push_back(this);
		return;
	}
	else
	{
		children[q0]->twoChildren(q0,q1,two);
		children[q1]->twoChildren(q0,q1,two);
	}
}

void QtBox::oneChildren(unsigned int q0, std::vector< QtBox * > & one)
{
	if(bEmpty)
	{
		one.push_back(this);
		return;
	}
	else
	{
		children[q0]->oneChildren(q0,one);
	}
}

void QtBox::leftmostChildren(std::vector< QtBox * > & two){ twoChildren(2,0,two); }
void QtBox::rightmostChildren(std::vector< QtBox * > & two){ twoChildren(3,1,two); }
void QtBox::bottommostChildren(std::vector< QtBox * > & two){ twoChildren(0,1,two); }
void QtBox::topmostChildren(std::vector< QtBox * > & two){ twoChildren(2,3,two); }

void QtBox::bottomLeftChildren(std::vector< QtBox * > & one){ oneChildren(0,one); }
void QtBox::bottomRightChildren(std::vector< QtBox * > & one){ oneChildren(1,one); }
void QtBox::topLeftChildren(std::vector< QtBox * > & one){ oneChildren(2,one); }
void QtBox::topRightChildren(std::vector< QtBox * > & one){ oneChildren(3,one); }

void QtBox::print()
{
	const double minx = minVertex[0];
	const double maxx = minx + width;
	const double miny = minVertex[1];
	const double maxy = miny + width;
	printf("minx=%e, maxx=%e\n",minx,maxx);
	printf("miny=%e, maxy=%e\n",miny,maxy);
}

bool QtBox::_splitCondition(void)
{
	// count number of particles per box
	double threshold = 0.2;
	if(level==ROOT_LEVEL)
		threshold = 0.05;

	double rng = ((double) rand() / RAND_MAX);
	//return (rng > threshold);
	return true;
}

void QtBox::_getChildBoxDimensions(const double parentMinVertex[2], double parentWidth, double childMinVertex[2], double & childWidth, const unsigned int coordIdx[2])
{
	assert(coordIdx[0]>=0 && coordIdx[0]<=1);
	assert(coordIdx[1]>=0 && coordIdx[1]<=1);

	childWidth = 0.5*parentWidth;
	childMinVertex[0] = parentMinVertex[0] + coordIdx[0]*childWidth;
	childMinVertex[1] = parentMinVertex[1] + coordIdx[1]*childWidth;
}

unsigned int QtBox::_invLevelToLevel(unsigned int _invLevel)
{
	return abs(_invLevel - ROOT_LEVEL);
}

unsigned int QtBox::_getLocationalCode(double x)
{
	return (unsigned int) (x * MAX_VAL);
}




