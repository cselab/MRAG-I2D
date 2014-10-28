/*
 *  QtBox.h
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

#ifndef QTBOX_H_
#define QTBOX_H_

#include <vector>
#include <iostream>
#include <cmath>
#include "assert.h"
#include "MRAGCommon.h"

struct QtBox
{
public:
	const unsigned int MAX_N_LEVELS;	// Number of possible levels in the quadtree
	unsigned int ROOT_LEVEL;			// Level of root cell (QT_N_LEVELS - 1)
	double MAX_VAL;						// For converting positions to locational codes
	const unsigned int invLevel;		// inverse level of the box (root=0...)
	const unsigned int idxX;			// space index x
	const unsigned int idxY;			// space index y
	unsigned int level;					// level of the box (root=ROOT_LEVEL, ..., leaves=0)
	const double width;					// box width
	QtBox * parent;						// Parent box
	QtBox * children[4];				// Children boxes contiguous in memory
	bool bEmpty; 						// am I a leaf node?
	double minVertex[2];		 		// minimum vertex of the box
	unsigned int xLocCode;				// X location code
	unsigned int yLocCode; 				// Y location code

	QtBox(const unsigned int _maxLevel, const double _width, const double _minVertex[2], const int _invLevel = 0, QtBox * _parent = NULL, const unsigned int _idxX = 0, const unsigned int _idxY = 0);
	~QtBox();

	I3 getIndices();
	I3 getIndicesTranslated(QtBox * cell);
	void split(std::vector<QtBox*> & leaves);
	void split(I3 input, std::vector<QtBox*> & leaves);
	bool isInside(const double p[2]);
	bool isMyLeftNeighbor(const QtBox * leftNeigh);
	bool isMyRightNeighbor(const QtBox * rightNeigh);
	bool isMyBottomNeighbor(const QtBox * bottomNeigh);
	bool isMyTopNeighbor(const QtBox * topNeigh);
	bool isTopLeftCornerNeighbor(const QtBox * cornerNeigh);
	bool isTopRightCornerNeighbor(const QtBox * cornerNeigh);
	bool isBottomLeftCornerNeighbor(const QtBox * cornerNeigh);
	bool isBottomRightCornerNeighbor(const QtBox * cornerNeigh);
	unsigned int spotQuadrant(const double p[2]);
	void twoChildren(unsigned int q0, unsigned int q1, std::vector< QtBox * > & two);
	void oneChildren(unsigned int q0, std::vector< QtBox * > & one);
	void leftmostChildren(std::vector< QtBox * > & two);
	void rightmostChildren(std::vector< QtBox * > & two);
	void bottommostChildren(std::vector< QtBox * > & two);
	void topmostChildren(std::vector< QtBox * > & two);
	void bottomLeftChildren(std::vector< QtBox * > & one);
	void bottomRightChildren(std::vector< QtBox * > & one);
	void topLeftChildren(std::vector< QtBox * > & one);
	void topRightChildren(std::vector< QtBox * > & one);
	void print();

protected:
	bool _splitCondition(void);
	void _getChildBoxDimensions(const double parentMinVertex[2], double parentWidth, double childMinVertex[2], double & childWidth, const unsigned int coordIdx[2]);
	unsigned int _invLevelToLevel(unsigned int _invLevel);
	unsigned int _getLocationalCode(double x);
};

#endif /* QTBOX_H_ */
