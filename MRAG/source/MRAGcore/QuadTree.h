/*
 * QuadTree.h
 *
 *  Created on: Jan 30, 2012
 *      Author: mgazzola
 *
 *  Point location, region location and neighbor searches were stripped out
 *  from the paper "Simple and efficient traversal methods for quadtrees and
 *  octrees" by Frisken and Perry.
 *
 *  THEIR CODE WAS THE MOST FAULTY EVER FUCKING MOTHERFUCKERS!!
 */



#ifndef QUADTREE_H_
#define QUADTREE_H_

#include "QtBox.h"
#include <map>
#include "MRAGCommon.h"

class QuadTree
{
protected:
	const unsigned int MAX_N_LEVELS;	// Number of possible levels in the quadtree
	unsigned int ROOT_LEVEL;			// Level of root cell (QT_N_LEVELS - 1)
	double MAX_VAL;						// For converting positions to locational codes
	QtBox * root;
	std::vector<QtBox*> leaves;
	std::map<I3,QtBox*> mapping;

public:
	QuadTree(const unsigned int _maxLevel, const double _width, const double _minVertex[2]);
	virtual ~QuadTree();

	//-------------------------------------------------------------------------------
	// Split the domain in sub-cells according to a predefined criterium
	//-------------------------------------------------------------------------------
	void split();
	void split(std::vector<I3> leavesInput);

	//-------------------------------------------------------------------------------
	// Macro to traverse a quadtree from a specified cell (typically the root cell)
	// to a leaf cell by following the x and y locational codes, xLocCode and
	// yLocCode. Upon entering, cell is the specified cell and nextLevel is one less
	// than the level of the specified cell. Upon termination, cell is the leaf cell
	// and nextLevel is one less than the level of the leaf cell.
	//-------------------------------------------------------------------------------
	// BUG: in their implementation the pointer was not modified
	//-------------------------------------------------------------------------------
	void traverse(QtBox ** cell,unsigned int & nextLevel, const unsigned int xLocCode, const unsigned int yLocCode);

	// Macro to traverse a quadtree from a specified cell to an offspring cell by
	// following the x and y locational codes, xLocCode and yLocCode. The offpring
	// cell is either at a specified level or is a leaf cell if a leaf cell is
	// reached before the specified level. Upon entering, cell is the specified
	// cell and nextLevel is one less than the level of the specified cell. Upon
	// termination, cell is the offspring cell and nextLevel is one less than the
	// level of the offspring cell.
	//-------------------------------------------------------------------------------
	// BUG: in their implementation the pointer was not modified
	//-------------------------------------------------------------------------------
	void traverseToLevel(QtBox ** cell, unsigned int & nextLevel, const unsigned int xLocCode, const unsigned int yLocCode, const unsigned int level);

	//-------------------------------------------------------------------------------
	// Macro for traversing a quadtree to a common ancestor of a specified cell
	// and its neighbor, whose x or y locational code differs from the cell's
	// corresponding x or y locational code by binaryDiff (determined by XOR'ing the
	// appropriate pair of x or y locational codes). Upon entering, cell is the
	// specified cell and cellLevel is the cell's level. Upon termination, cell is
	// the common ancestor and cellLevel is the common ancestor's level.
	//-------------------------------------------------------------------------------
	void getCommonAncestor(QtBox ** cell, unsigned int & cellLevel, const unsigned int binaryDiff);

	//-------------------------------------------------------------------------------
	// Locate the leaf cell containing the specified point p, where p lies in
	// [0,1)x[0,1).
	//-------------------------------------------------------------------------------
	QtBox * locateCell(const double p[2]);
	I3 locateCellIndex(const double p[2]);

	//-------------------------------------------------------------------------------
	// Locate the smallest cell that entirely contains a rectangular region defined
	// by its bottom-left vertex v0 and its top-right vertex v1, where v0 and v1
	// lie in [0,1)x[0,1).
	//-------------------------------------------------------------------------------
	QtBox * locateRegion(const double v0[2], const double v1[2]);

	//-------------------------------------------------------------------------------
	// Locate the left edge neighbor of the same size or larger than a specified
	// cell. Periodic tree is assumed, therefore the returned pointer cannot be NULL.
	// To each box a bool is associated: true if the box is a ghost, false if the box
	// is adjacent within the unit square doamain
	//-------------------------------------------------------------------------------
	std::pair< QtBox *, bool > locateLeftNeighborGELevelPeriodic(QtBox * cell);

	//-------------------------------------------------------------------------------
	// Locate the right edge neighbor of the same size or larger than a specified
	// cell. Periodic tree is assumed, therefore the returned pointer cannot be NULL.
	// To each box a bool is associated: true if the box is a ghost, false if the box
	// is adjacent within the unit square doamain
	//-------------------------------------------------------------------------------
	std::pair< QtBox *, bool > locateRightNeighborGELevelPeriodic(QtBox * cell);

	//-------------------------------------------------------------------------------
	// Locate the bottom edge neighbor of the same size or larger than a specified
	// cell. Periodic tree is assumed, therefore the returned pointer cannot be NULL.
	// To each box a bool is associated: true if the box is a ghost, false if the box
	// is adjacent within the unit square doamain
	//-------------------------------------------------------------------------------
	std::pair< QtBox *, bool > locateBottomNeighborGELevelPeriodic(QtBox * cell);

	//-------------------------------------------------------------------------------
	// Locate the top edge neighbor of the same size or larger than a specified
	// cell. Periodic tree is assumed, therefore the returned pointer cannot be NULL.
	// To each box a bool is associated: true if the box is a ghost, false if the box
	// is adjacent within the unit square doamain
	//-------------------------------------------------------------------------------
	std::pair< QtBox *, bool > locateTopNeighborGELevelPeriodic(QtBox * cell);

	//-------------------------------------------------------------------------------
	// Locate all the left edge neighbor of a specified
	// cell. Periodic tree is assumed, therefore the returned pointer cannot be NULL.
	// To each box a bool is associated: true if the box is a ghost, false if the box
	// is adjacent within the unit square doamain
	//-------------------------------------------------------------------------------
	std::vector< std::pair< QtBox *, bool > > locateLeftNeighborsPeriodic(QtBox * cell);

	//-------------------------------------------------------------------------------
	// Locate all the right edge neighbors of a specified
	// cell. Periodic tree is assumed, therefore the returned pointer cannot be NULL.
	// To each box a bool is associated: true if the box is a ghost, false if the box
	// is adjacent within the unit square doamain
	//-------------------------------------------------------------------------------
	std::vector< std::pair< QtBox *, bool > > locateRightNeighborsPeriodic(QtBox * cell);

	//-------------------------------------------------------------------------------
	// Locate all the bottom edge neighbors of a specified
	// cell. Periodic tree is assumed, therefore the returned pointer cannot be NULL.
	// To each box a bool is associated: true if the box is a ghost, false if the box
	// is adjacent within the unit square doamain
	//-------------------------------------------------------------------------------
	std::vector< std::pair< QtBox *, bool > > locateBottomNeighborsPeriodic(QtBox * cell);

	//-------------------------------------------------------------------------------
	// Locate all the top edge neighbors of a specified
	// cell. Periodic tree is assumed, therefore the returned pointer cannot be NULL.
	// To each box a bool is associated: true if the box is a ghost, false if the box
	// is adjacent within the unit square doamain
	//-------------------------------------------------------------------------------
	std::vector< std::pair< QtBox *, bool > > locateTopNeighborsPeriodic(QtBox * cell);

	//-------------------------------------------------------------------------------
	// Locate all neighbors of a specified
	// cell. Periodic tree is assumed, therefore the returned pointer cannot be NULL.
	// To each box a bool is associated: true if the box is a ghost, false if the box
	// is adjacent within the unit square doamain
	//-------------------------------------------------------------------------------
	std::vector< std::pair< QtBox *, bool > > locateAllNeighborsPeriodic(QtBox * cell);

	//-------------------------------------------------------------------------------
	// Determined the neighborhood as required by MRAG.
	// The first entry in the main pair is the vector containing all adjacent blocks.
	// The second entry in the main pair is a vector of pairs containing all
	// non existing (outside square unit domain) ghosts. The first entry of this second
	// pair is the non-existing ghost indices, while the second is the non-existing
	// ghost indices, re-mapped periodically
	//-------------------------------------------------------------------------------
	std::pair< std::vector<I3>, std::vector<std::pair<I3, I3> > > locateNeighborhoodMRAG(QtBox * cell);

	//-------------------------------------------------------------------------------
	// Get MRAG neighbors
	//-------------------------------------------------------------------------------
	std::vector< std::pair< std::vector<I3>, std::vector<std::pair<I3, I3> > > > getNeighborsMRAG(std::vector<I3> const & input);
};

#endif /* QUADTREE_H_ */
