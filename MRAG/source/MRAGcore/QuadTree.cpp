/*
 * QuadTree.cpp
 *
 *  Created on: Jan 30, 2012
 *      Author: mgazzola
 */
#include "QuadTree.h"

QuadTree::QuadTree(const unsigned int _maxLevel, const double _width, const double _minVertex[2]): MAX_N_LEVELS(_maxLevel)
{
	assert(MAX_N_LEVELS>=1);
	ROOT_LEVEL = MAX_N_LEVELS-1;
	MAX_VAL = pow(2,ROOT_LEVEL);

	root = new QtBox(MAX_N_LEVELS,_width,_minVertex);
}

QuadTree::~QuadTree()
{
	if(root!=NULL)
		delete root;
}

void QuadTree::split()
{
	leaves.clear();
	root->split(leaves);

	//for(std::vector<QtBox*>::iterator it=leaves.begin(); it<leaves.end(); it++)
	//	mapping[I3((*it)->idxX,(*it)->idxY,(*it)->invLevel)] = (*it);
}

void QuadTree::split(std::vector<I3> input)
{
	leaves.clear();
	for(std::vector<I3>::iterator it=input.begin(); it<input.end(); it++)
		root->split((*it),leaves);

	//for(std::vector<QtBox*>::iterator it=leaves.begin(); it<leaves.end(); it++)
	//	mapping[I3((*it)->idxX,(*it)->idxY,(*it)->invLevel)] = (*it);
}

void QuadTree::traverse(QtBox ** cell,unsigned int & nextLevel, const unsigned int xLocCode, const unsigned int yLocCode)
{
	//-------------------------------------------------------------------------------
	// BUG: in their code wrong indexing
	//-------------------------------------------------------------------------------
	while( !(*cell)->bEmpty )
	{
		const unsigned int childBranchBit = 1 << (nextLevel);
		const unsigned int childIndexX = (((xLocCode) & childBranchBit) >> (nextLevel));
		const unsigned int childIndexY = (((yLocCode) & childBranchBit) >> (nextLevel));
		const unsigned int childIndex = childIndexX + 2*childIndexY;
		--nextLevel;
		(*cell) = (*cell)->children[childIndex];
	}
}

void QuadTree::traverseToLevel(QtBox ** cell, unsigned int & nextLevel, const unsigned int xLocCode, const unsigned int yLocCode, const unsigned int level)
{
	//-------------------------------------------------------------------------------
	// BUG1: in their code wrong indexing
	// BUG2: in their code the break was at the end of the loop
	//-------------------------------------------------------------------------------
	unsigned int n = (nextLevel) - (level) + 1;
	while (n--)
	{
		if ( (*cell)->bEmpty ) break;
		unsigned int childBranchBit = 1 << (nextLevel);
		const unsigned int childIndexX = (((xLocCode) & childBranchBit) >> (nextLevel));
		const unsigned int childIndexY = (((yLocCode) & childBranchBit) >> (nextLevel));
		const unsigned int childIndex = childIndexX + 2*childIndexY;
		--nextLevel;
		(*cell) = (*cell)->children[childIndex];
	}
}

void QuadTree::getCommonAncestor(QtBox ** cell, unsigned int & cellLevel, const unsigned int binaryDiff)
{
	//-------------------------------------------------------------------------------
	// BUG: in their code if parent==NULL it doesnt break!\
	// Who da fuck are this motherfcukers?!
	//-------------------------------------------------------------------------------
	while( (binaryDiff) & (1 << (cellLevel)) )
	{
		QtBox * parent = (*cell)->parent;
		if( parent == NULL ) break;
		(*cell) = (*cell)->parent;
		(cellLevel)++;
	}
}

QtBox * QuadTree::locateCell(const double p[2])
{
	//----Determine the x and y locational codes of the point's position. Refer
	//----to [King2001] for more efficient methods for converting floating point
	//----numbers to integers.
	const unsigned int xLocCode = (unsigned int) (p[0] * MAX_VAL);
	const unsigned int yLocCode = (unsigned int) (p[1] * MAX_VAL);
	//----Follow the branching patterns of the locational codes from the root cell
	//----to locate the leaf cell containing p
	QtBox * cell = root;
	unsigned int nextLevel = ROOT_LEVEL - 1;
	traverse(&cell,nextLevel,xLocCode,yLocCode);
	return cell;
}

I3 QuadTree::locateCellIndex(const double p[2])
{
	QtBox * cell = locateCell(p);
	return I3(cell->idxX, cell->idxY, cell->invLevel);
}

QtBox * QuadTree::locateRegion(const double v0[2], const double v1[2])
{
	//----Determine the x and y locational codes of the region boundaries. Refer
	//----to [King2001] for more efficient methods for converting floating point
	//----numbers to integers.
	const unsigned int x0LocCode = (unsigned int) (v0[0] * MAX_VAL);
	const unsigned int y0LocCode = (unsigned int) (v0[1] * MAX_VAL);
	const unsigned int x1LocCode = (unsigned int) (v1[0] * MAX_VAL);
	const unsigned int y1LocCode = (unsigned int) (v1[1] * MAX_VAL);
	//----Determine the XOR'ed pairs of locational codes of the region boundaries
	const unsigned int xDiff = x0LocCode ^ x1LocCode;
	const unsigned int yDiff = y0LocCode ^ y1LocCode;
	//----Determine the level of the smallest possible cell entirely containing
	//----the region
	QtBox *cell = root;
	unsigned int level = ROOT_LEVEL;
	unsigned int minLevel = ROOT_LEVEL;
	while (!(xDiff & (1 << level)) && level) level--;
	while (!(yDiff & (1 << minLevel)) && (minLevel > level)) minLevel--;
	minLevel++;

	//----Follow the branching patterns of the locational codes of v0 from the
	//----root cell to the smallest cell entirely containing the region
	level = ROOT_LEVEL-1;
	traverseToLevel(&cell,level,x0LocCode,y0LocCode,minLevel);

	//-------------------------------------------------------------------------------
	// BUG1: in their code it was always returning the SECOND smallest box
	// Added the following piece of code
	//-------------------------------------------------------------------------------
	if(!cell->bEmpty)
	{
		const unsigned int idxV0 = cell->spotQuadrant(v0);
		const unsigned int idxV1 = cell->spotQuadrant(v1);
		if( idxV0 == idxV1 )
			cell = cell->children[idxV0];
	}

	return(cell);
}

std::pair< QtBox *, bool > QuadTree::locateLeftNeighborGELevelPeriodic(QtBox * cell)
{
	//----Get cell's x and y locational codes and the x locational code of the
	//----cell's smallest possible left neighbor
	const bool leftEdge = (cell->xLocCode==0);
	const unsigned int xLocCodePeriodic = (unsigned int)(1.0 * MAX_VAL);
	const unsigned int xLocCode = leftEdge?xLocCodePeriodic:cell->xLocCode;
	const unsigned int yLocCode = cell->yLocCode;
	const unsigned int xLeftLocCode = xLocCode - 0x00000001;

	//----Determine the smallest common ancestor of the cell and the cell's
	//----smallest possible left neighbor
	unsigned int cellLevel = cell->level;
	unsigned int nextLevel = cell->level;
	const unsigned int diff = xLocCode ^ xLeftLocCode;
	QtBox *pCell = cell;
	getCommonAncestor(&pCell,nextLevel,diff);

	//----Start from the smallest common ancestor and follow the branching
	//----patterns of the locational codes downward to the smallest left
	//----neighbor of size greater than or equal to cell
	nextLevel--;
	traverseToLevel(&pCell,nextLevel,xLeftLocCode,yLocCode,cellLevel);
	std::pair< QtBox *, bool > tag(pCell,leftEdge);
	return(tag);
}

std::pair< QtBox *, bool > QuadTree::locateRightNeighborGELevelPeriodic(QtBox * cell)
{
	//----Get cell's x and y locational codes and the x locational code of the
	//----cell's right neighbors
	const unsigned int binaryCellSize = 1 << cell->level;
	const bool rightEdge = ((cell->xLocCode + binaryCellSize) >= (1 << ROOT_LEVEL));
	const unsigned int xLocCodePeriodic = (unsigned int)(0.0 * MAX_VAL);
	const unsigned int xLocCode = rightEdge?xLocCodePeriodic:cell->xLocCode;
	const unsigned int yLocCode = cell->yLocCode;
	const unsigned int xRightLocCode = rightEdge?xLocCode:(xLocCode+binaryCellSize);

	//----Determine the smallest common ancestor of the cell and the cell's
	//----right neighbors
	const unsigned int cellLevel = cell->level;
	unsigned int nextLevel = rightEdge?ROOT_LEVEL:cell->level;
	QtBox *pCell = rightEdge?root:cell;
	if(!rightEdge)
	{
		const unsigned int diff = xLocCode ^ xRightLocCode;
		getCommonAncestor(&pCell,nextLevel,diff);
	}

	//----Start from the smallest common ancestor and follow the branching
	//----patterns of the locational codes downward to the smallest right
	//----neighbor of size greater than or equal to cell
	nextLevel--;
	traverseToLevel(&pCell,nextLevel,xRightLocCode,yLocCode,cellLevel);
	std::pair< QtBox *, bool > tag(pCell,rightEdge);
	return(tag);
}

std::pair< QtBox *, bool > QuadTree::locateBottomNeighborGELevelPeriodic(QtBox * cell)
{
	//----Get cell's x and y locational codes and the y locational code of the
	//----cell's smallest possible bottom neighbor
	const bool bottomEdge = (cell->yLocCode==0);
	const unsigned int xLocCode = cell->xLocCode;
	const unsigned int yLocCodePeriodic = (unsigned int)(1.0 * MAX_VAL);
	const unsigned int yLocCode = bottomEdge?yLocCodePeriodic:cell->yLocCode;
	const unsigned int yBottomLocCode = yLocCode - 0x00000001;

	//----Determine the smallest common ancestor of the cell and the cell's
	//----smallest possible bottom neighbor
	unsigned int cellLevel = cell->level;
	unsigned int nextLevel = cell->level;
	const unsigned int diff = yLocCode ^ yBottomLocCode;
	QtBox *pCell = cell;
	getCommonAncestor(&pCell,nextLevel,diff);

	//----Start from the smallest common ancestor and follow the branching
	//----patterns of the locational codes downward to the smallest bottom
	//----neighbor of size greater than or equal to cell
	nextLevel--;
	traverseToLevel(&pCell,nextLevel,xLocCode,yBottomLocCode,cellLevel);
	std::pair< QtBox *, bool > tag(pCell,bottomEdge);
	return(tag);
}

std::pair< QtBox *, bool > QuadTree::locateTopNeighborGELevelPeriodic(QtBox * cell)
{
	//----Get cell's x and y locational codes and the x locational code of the
	//----cell's top neighbors
	const unsigned int binaryCellSize = 1 << cell->level;
	const bool topEdge = ((cell->yLocCode + binaryCellSize) >= (1 << ROOT_LEVEL));
	const unsigned int xLocCode = cell->xLocCode;
	const unsigned int yLocCodePeriodic = (unsigned int)(0.0 * MAX_VAL);
	const unsigned int yLocCode = topEdge?yLocCodePeriodic:cell->yLocCode;
	const unsigned int yTopLocCode = topEdge?yLocCode:(yLocCode+binaryCellSize);

	//----Determine the smallest common ancestor of the cell and the cell's
	//----top neighbors
	const unsigned int cellLevel = cell->level;
	unsigned int nextLevel = topEdge?ROOT_LEVEL:cell->level;
	QtBox *pCell = topEdge?root:cell;
	if(!topEdge)
	{
		const unsigned int diff = yLocCode ^ yTopLocCode;
		getCommonAncestor(&pCell,nextLevel,diff);
	}

	//----Start from the smallest common ancestor and follow the branching
	//----patterns of the locational codes downward to the smallest top
	//----neighbor of size greater than or equal to cell
	nextLevel--;
	traverseToLevel(&pCell,nextLevel,xLocCode,yTopLocCode,cellLevel);
	std::pair< QtBox *, bool > tag(pCell,topEdge);
	return(tag);
}


std::vector< std::pair< QtBox *, bool > > QuadTree::locateLeftNeighborsPeriodic(QtBox * cell)
{
	std::vector< QtBox * > cellChildren;
	std::pair< QtBox *, bool > b = locateLeftNeighborGELevelPeriodic(cell);
	b.first->rightmostChildren(cellChildren);
	std::vector< std::pair< QtBox *, bool > > tagged;
	for(std::vector< QtBox * >::iterator it=cellChildren.begin(); it<cellChildren.end(); it++)
		tagged.push_back( std::pair< QtBox*,bool>((*it), b.second) );

	return tagged;
}

std::vector< std::pair< QtBox *, bool > > QuadTree::locateRightNeighborsPeriodic(QtBox * cell)
{
	std::vector< QtBox * > cellChildren;
	std::pair< QtBox *, bool > b = locateRightNeighborGELevelPeriodic(cell);
	b.first->leftmostChildren(cellChildren);
	std::vector< std::pair< QtBox *, bool > > tagged;
	for(std::vector< QtBox * >::iterator it=cellChildren.begin(); it<cellChildren.end(); it++)
		tagged.push_back( std::pair< QtBox*,bool>((*it), b.second) );

	return tagged;
}

std::vector< std::pair< QtBox *, bool > > QuadTree::locateBottomNeighborsPeriodic(QtBox * cell)
{
	std::vector< QtBox * > cellChildren;
	std::pair< QtBox *, bool > b = locateBottomNeighborGELevelPeriodic(cell);
	b.first->topmostChildren(cellChildren);
	std::vector< std::pair< QtBox *, bool > > tagged;
	for(std::vector< QtBox * >::iterator it=cellChildren.begin(); it<cellChildren.end(); it++)
		tagged.push_back( std::pair< QtBox*,bool>((*it), b.second) );

	return tagged;
}

std::vector< std::pair< QtBox *, bool > > QuadTree::locateTopNeighborsPeriodic(QtBox * cell)
{
	std::vector< QtBox * > cellChildren;
	std::pair< QtBox *, bool > b = locateTopNeighborGELevelPeriodic(cell);
	b.first->bottommostChildren(cellChildren);
	std::vector< std::pair< QtBox *, bool > > tagged;
	for(std::vector< QtBox * >::iterator it=cellChildren.begin(); it<cellChildren.end(); it++)
		tagged.push_back( std::pair< QtBox*,bool>((*it), b.second) );

	return tagged;
}

std::vector< std::pair< QtBox *, bool > > QuadTree::locateAllNeighborsPeriodic(QtBox * cell)
{
	std::vector< std::vector< std::pair< QtBox *, bool > > > total;

	// Collect left neighbors
	std::pair< QtBox *, bool > left = locateLeftNeighborGELevelPeriodic(cell);
	{
		std::vector< QtBox * > leftChildren;
		left.first->rightmostChildren(leftChildren);
		std::vector< std::pair< QtBox *, bool > > tagged;
		for(std::vector< QtBox * >::iterator it=leftChildren.begin(); it<leftChildren.end(); it++)
			tagged.push_back( std::pair< QtBox*,bool>((*it), left.second) );
		total.push_back(tagged);
#ifndef NDEBUG
		for(unsigned int j=0; j<leftChildren.size(); j++){ cell->isMyLeftNeighbor(leftChildren[j]); }
#endif
	}

	// Collect right neighbors
	std::pair< QtBox *, bool > right = locateRightNeighborGELevelPeriodic(cell);
	{
		std::vector< QtBox * > rightChildren;
		right.first->leftmostChildren(rightChildren);
		std::vector< std::pair< QtBox *, bool > > tagged;
		for(std::vector< QtBox * >::iterator it=rightChildren.begin(); it<rightChildren.end(); it++)
			tagged.push_back( std::pair< QtBox*,bool>((*it), right.second) );
		total.push_back(tagged);
#ifndef NDEBUG
		for(unsigned int j=0; j<rightChildren.size(); j++){ cell->isMyRightNeighbor(rightChildren[j]); }
#endif
	}

	// Collect bottom neighbors
	std::pair< QtBox *, bool > bottom = locateBottomNeighborGELevelPeriodic(cell);
	{
		std::vector< QtBox * > bottomChildren;
		bottom.first->topmostChildren(bottomChildren);
		std::vector< std::pair< QtBox *, bool > > tagged;
		for(std::vector< QtBox * >::iterator it=bottomChildren.begin(); it<bottomChildren.end(); it++)
			tagged.push_back( std::pair< QtBox*,bool>((*it), bottom.second) );
		total.push_back(tagged);
#ifndef NDEBUG
		for(unsigned int j=0; j<bottomChildren.size(); j++){ cell->isMyBottomNeighbor(bottomChildren[j]); }
#endif
	}

	// Collect top neighbors
	std::pair< QtBox *, bool > top = locateTopNeighborGELevelPeriodic(cell);
	{
		std::vector< QtBox * > topChildren;
		top.first->bottommostChildren(topChildren);
		std::vector< std::pair< QtBox *, bool > > tagged;
		for(std::vector< QtBox * >::iterator it=topChildren.begin(); it<topChildren.end(); it++)
			tagged.push_back( std::pair< QtBox*,bool>((*it), top.second) );
		total.push_back(tagged);
#ifndef NDEBUG
		for(unsigned int j=0; j<topChildren.size(); j++){ cell->isMyTopNeighbor(topChildren[j]); }
#endif
	}


	// Top-left corner (it gives you a corner ONLY if it is not a lateral block!!)
	{
		std::pair< QtBox *, bool > lt = locateTopNeighborGELevelPeriodic(left.first);
		std::vector< QtBox * > first;
		lt.first->bottomRightChildren(first);

		std::pair< QtBox *, bool > tl = locateLeftNeighborGELevelPeriodic(top.first);
		std::vector< QtBox * > second;
		tl.first->bottomRightChildren(second);

		assert( first.size()==1 && second.size() );
		if( first[0]==second[0] )
		{
			std::vector< std::pair< QtBox *, bool > > tagged;
			tagged.push_back(std::pair< QtBox*,bool>(first[0], (lt.second || tl.second) ));
			total.push_back(tagged);
#ifndef NDEBUG
			for(unsigned int j=0; j<first.size(); j++){ cell->isTopLeftCornerNeighbor(first[j]); }
#endif
		}
	}

	// Bottom-left corner (it gives you a corner ONLY if it is not a lateral block!!)
	{
		std::pair< QtBox *, bool > lb = locateBottomNeighborGELevelPeriodic(left.first);
		std::vector< QtBox * > first;
		lb.first->topRightChildren(first);

		std::pair< QtBox *, bool > bl = locateLeftNeighborGELevelPeriodic(bottom.first);
		std::vector< QtBox * > second;
		bl.first->topRightChildren(second);

		assert( first.size()==1 && second.size() );
		if( first[0]==second[0] )
		{
			std::vector< std::pair< QtBox *, bool > > tagged;
			tagged.push_back(std::pair< QtBox*,bool>(first[0], (lb.second || bl.second) ));
			total.push_back(tagged);
#ifndef NDEBUG
			for(unsigned int j=0; j<first.size(); j++){ cell->isBottomLeftCornerNeighbor(first[j]); }
#endif
		}
	}

	// Top-right corner (it gives you a corner ONLY if it is not a lateral block!!)
	{
		std::pair< QtBox *, bool > ru = locateTopNeighborGELevelPeriodic(right.first);
		std::vector< QtBox * > first;
		ru.first->bottomLeftChildren(first);

		std::pair< QtBox *, bool > ur = locateRightNeighborGELevelPeriodic(top.first);
		std::vector< QtBox * > second;
		ur.first->bottomLeftChildren(second);

		assert( first.size()==1 && second.size() );
		if( first[0]==second[0] )
		{
			std::vector< std::pair< QtBox *, bool > > tagged;
			tagged.push_back(std::pair< QtBox*,bool>(first[0], (ru.second || ur.second) ));
			total.push_back(tagged);
#ifndef NDEBUG
			for(unsigned int j=0; j<first.size(); j++){ cell->isTopRightCornerNeighbor(first[j]); }
#endif
		}
	}

	// Bottom-right corner (it gives you a corner ONLY if it is not a lateral block!!)
	{
		std::pair< QtBox *, bool > rb = locateBottomNeighborGELevelPeriodic(right.first);
		std::vector< QtBox * > first;
		rb.first->topLeftChildren(first);

		std::pair< QtBox *, bool > br = locateRightNeighborGELevelPeriodic(bottom.first);
		std::vector< QtBox * > second;
		br.first->topLeftChildren(second);

		assert( first.size()==1 && second.size() );
		if( first[0]==second[0] )
		{
			std::vector< std::pair< QtBox *, bool > > tagged;
			tagged.push_back(std::pair< QtBox*,bool>(first[0], (rb.second || br.second) ));
			total.push_back(tagged);
#ifndef NDEBUG
			for(unsigned int j=0; j<first.size(); j++){ cell->isBottomRightCornerNeighbor(first[j]); }
#endif
		}
	}

	// Put all neighbors together
	std::vector< std::pair< QtBox *, bool > > totLinear;
	for(std::vector< std::vector< std::pair< QtBox *, bool > > >::iterator it=total.begin(); it<total.end(); it++)
		for(std::vector< std::pair< QtBox *, bool > >::iterator it2=(*it).begin(); it2<(*it).end(); it2++)
			totLinear.push_back(*it2);

	return totLinear;
}

std::pair< std::vector<I3>, std::vector<std::pair<I3, I3> > > QuadTree::locateNeighborhoodMRAG(QtBox * cell)
{
	std::vector< std::pair<I3, I3> > ghosts;
	std::vector<I3> adjacent;
	std::vector< std::pair< QtBox *, bool > > tagged = locateAllNeighborsPeriodic(cell);

	for(std::vector< std::pair< QtBox *, bool > >::iterator it=tagged.begin(); it<tagged.end(); it++)
	{
		if( (*it).second )
			ghosts.push_back( std::pair<I3, I3>( (*it).first->getIndicesTranslated(cell), (*it).first->getIndices() ) ); // Ghosts
		else
			adjacent.push_back( (*it).first->getIndices() ); // Adjacent
	}

	std::pair< std::vector<I3>, std::vector<std::pair<I3, I3> > > outputMRAG(adjacent,ghosts);

	return outputMRAG;
}

std::vector< std::pair< std::vector<I3>, std::vector<std::pair<I3, I3> > > > QuadTree::getNeighborsMRAG(std::vector<I3> const & input)
{
	std::vector< std::pair< std::vector<I3>, std::vector<std::pair<I3, I3> > > > output;

	// Generate corresponding tree
	split(input);
	assert(input.size()==leaves.size());

#ifndef NDEBUG
	for(unsigned int i=0; i<input.size(); i++)
	{
		assert(input[i].i[0]==leaves[i]->idxX);
		assert(input[i].i[1]==leaves[i]->idxY);
		assert(input[i].i[2]==leaves[i]->invLevel);
	}
#endif

	// Prepare output
	for(std::vector<QtBox*>::iterator it=leaves.begin(); it<leaves.end(); it++)
		output.push_back( locateNeighborhoodMRAG( (*it) ) );

	return output;
}








