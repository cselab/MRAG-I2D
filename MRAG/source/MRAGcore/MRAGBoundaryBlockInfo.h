 /*
 *  BoundaryBlockInfo.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/24/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once
#undef min
#undef max
#include <vector>
#undef min
#undef max
#include <map>
#undef min
#undef max
#include <set>
#undef min
#undef max

using namespace std;

#include "MRAGCommon.h"
#include "MRAGHuffmanEncoder.h"
#include "MRAGGridNode.h"

namespace MRAG
{
	struct PointIndex
	{
		int blockID;
		/*unsigned short*/  int index;
		
		PointIndex& operator=(const PointIndex culo)
		{
			blockID = culo.blockID;
			index = culo.index;
			
			return *this;
		}
		
		//why should i implement this??
		PointIndex& operator=(pair<PointIndex, int> culo)
		{
			blockID = culo.first.blockID;
			index = culo.first.index;
			return *this;
		}
		
		bool operator<(const PointIndex culo) const
		{
			return blockID<culo.blockID || (blockID==culo.blockID && index<culo.index);
		}
	};
	
	struct IndexWP
	{
		/*unsigned short*/ int point_index;
		unsigned char weights_index[3];
	};
	
	struct BlockOfGhosts
	{
		unsigned int start;
		unsigned int nGhosts;
	};
	
	struct FaceInfoHeader
	{
		int start[3], end[3];
		unsigned int start_data;
	};
	
	struct BoundaryBlockHeader
	{
		int memsize;
		FaceInfoHeader faces[6];
		int points;
		int weights;
		unsigned int start_weights;
	};
	
	struct BoundaryInfoHeader
	{
		int nBlocks;
		int memsize;
	};
	

	struct BoundaryInfoBlock
	{	
		enum BBIState {BBIState_Initialized, BBIState_Locked, BBIState_Unlocked};
	private:
		
		BBIState state;
		int nLocks;
		int block_size[3];
		vector<PointIndex> indexPool;
		
		const GridNode& node;
		vector<GridNode*> neighbors;
		vector<GridNode> myneighbors;
		int diff_res,  nof_neighbors;
		
		typedef vector<IndexWP> ReconstructionInfo;
		vector< ReconstructionInfo > ghosts;
	
		HuffmanEncoder<unsigned short> encodedInstructionSizes;
		Encoder<unsigned char> encodedInstructionItemsWs;
		Encoder<unsigned short> encodedInstructionItemsPts;
		Encoder<unsigned char> vBlockID_encodedPointIndices3D;
		vector< pair<int, int> > vBlockID_Points;
		bool bCompressed;
		
		void _compress();
		void _decompress();
		void _discard_decompression();
		
		int _get_diff_res()
		{
			//assert(node != NULL);
			
			int min_lev = node.level;
			int max_lev = node.level;
			assert(nof_neighbors == neighbors.size());
#ifndef NDEBUG
			for(int i=0; i<nof_neighbors; i++)
			{
				GridNode& a = myneighbors[i];
				GridNode& b = *neighbors[i];
				assert(a.isEmpty == false);
				assert(a.index[0] == b.index[0]);
				assert(a.index[1] == b.index[1]);
				assert(a.index[2] == b.index[2]);
				assert(a.level == b.level);
				assert(a.blockID == b.blockID);
				assert(a.isEmpty == b.isEmpty);
			}
#endif
			
			for(vector<GridNode>::iterator it = myneighbors.begin(); it!= myneighbors.end(); ++it)
			{
				//assert(it != NULL);
				assert(it->level >= 0);
				
				min_lev = min(it->level, min_lev);
				max_lev = max(it->level, max_lev);
				
				assert(it->isEmpty == false);
			}
			
			assert(node.isEmpty == false);
			//assert(diff_res == max_lev - min_lev);
			
			return max_lev - min_lev;
		}
	public:
		
		BlockOfGhosts boundary[27];
		vector<int> dependentBlockIDs;
		vector<double> weightsPool;
		
		void lock()
		{
			if (state == BBIState_Unlocked)
				
			if(bCompressed)
				_decompress();
			
			nLocks++;
			
			state = BBIState_Locked;
		}
		
		inline const vector<PointIndex>& getIndexPool() const
		{
			assert(state == BBIState_Initialized || state == BBIState_Locked);
			return indexPool;
		}
		
		const vector< ReconstructionInfo >& getGhosts() const
		{
			assert(state == BBIState_Initialized || state == BBIState_Locked);
			return ghosts;
		}
		
		void release()
		{
			assert(nLocks>0);
			assert(state == BBIState_Initialized || state == BBIState_Locked);
			
			nLocks--;
			
			if(nLocks==0)
			{
				const double oldMB = getMemorySize();
				if (state == BBIState_Initialized  && oldMB>128. && !bCompressed)
				{
					_compress();
					_discard_decompression();
					
					bCompressed = true;
					
			//		const double newMB = getMemorySize();
			//		const int this_diff_res = _get_diff_res();
			//		assert(diff_res ==this_diff_res);
			//		printf("Compression factor %f, (%f MB -> %f MB) diff_res=%d\n", oldMB/newMB, oldMB, newMB, diff_res);
				}
				else if (bCompressed)
					_discard_decompression();
				else 
				{
			//		const double newMB = getMemorySize();
			//		const int this_diff_res = _get_diff_res();
			//		assert(diff_res ==this_diff_res);
			//		printf("Compression factor %f, (%f MB -> %f MB) diff_res=%d\n", oldMB/newMB, oldMB, newMB, diff_res);
				}

				
			}
			
			state = BBIState_Unlocked;
		}
		
		void * createBBPack();
		
		float getMemorySize() const;
		
		BoundaryInfoBlock(const int block_size_[3], const GridNode& b, const vector<GridNode*>& neighbors): 
			node(b), neighbors(neighbors),
			indexPool(), weightsPool(), ghosts(), 
			dependentBlockIDs(), state(BBIState_Initialized), nLocks(1),
			encodedInstructionSizes(), encodedInstructionItemsWs(), encodedInstructionItemsPts(),
			vBlockID_encodedPointIndices3D(), vBlockID_Points(), bCompressed(false)
		{
			nof_neighbors = neighbors.size();
			
			
			for(vector<GridNode*>::const_iterator it=neighbors.begin(); it!=neighbors.end(); ++it)
				myneighbors.push_back(**it);
			
			diff_res = _get_diff_res();
			
			block_size[0] = block_size_[0];
			block_size[1] = block_size_[1];
			block_size[2] = block_size_[2];
		}
		
		bool valid_content()
		{
			if (node.blockID == -1) return false;
			
		/*	for(int i=0; i<nof_neighbors; i++)
			{
				GridNode& a = myneighbors[i];
				GridNode& b = *neighbors[i];
				assert(a.isEmpty == false);
				assert(a.index[0] == b.index[0]);
				assert(a.index[1] == b.index[1]);
				assert(a.index[2] == b.index[2]);
				assert(a.level == b.level);
				assert(a.blockID == b.blockID);
				assert(a.isEmpty == b.isEmpty);
			}*/
			
			for(vector<GridNode*>::iterator it = neighbors.begin(); it!= neighbors.end(); ++it)
			{
				if ((*it)->blockID == -1) return false;
			}
			
			return true;
		}
		
	/*	~BoundaryInfoBlock()
		{
            indexPool.clear();
			weightsPool.clear();
			dependentBlockIDs.clear();
			for(int i=0; i<ghosts.size(); i++) 
				ghosts[i].clear();
			
			ghosts.clear();
		}*/
		
		void check() 
		{
			const int diff_res = _get_diff_res();
			float newMB = getMemorySize();
			printf("%f MB diff_res=%d\n", newMB, diff_res);
		}
		
		template<typename WaveletsType, typename BlockType> friend	class MRAG_BBInfoCreator;
	};
		
	struct BoundaryInfo
	{
		char stencil_start[3], stencil_end[3];
		map<int, BoundaryInfoBlock*> boundaryInfoOfBlock;
        std::map<int, std::vector<BlockInfo> > neighbors_vInfo;//wim added
		
		//set<int> invalidBlocks;
		
		float getMemorySize() const;
		void clear();
		void erase(int blockID, bool bCritical=true);
		BoundaryInfo():boundaryInfoOfBlock(){}
		~BoundaryInfo() { clear(); }
		
		void check()
		{
			for(map<int, BoundaryInfoBlock*>::iterator it= boundaryInfoOfBlock.begin(); it != boundaryInfoOfBlock.end(); ++it)
			{
				assert(it->second!=NULL);
				
				it->second->check();
			}
		}
		
		set<int> block_to_erase()
		{
			set<int> r;
			
			for(map<int, BoundaryInfoBlock*>::iterator itBBI=boundaryInfoOfBlock.begin(); itBBI!= boundaryInfoOfBlock.end(); ++itBBI)
			{
				if (!itBBI->second->valid_content())
					r.insert(itBBI->first);
			}
			
			return r;
		}
	};
}
