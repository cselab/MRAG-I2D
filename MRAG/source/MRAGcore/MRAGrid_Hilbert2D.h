/*
 *  MRAGrid.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/23/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once
#undef min
#undef max
#include <vector>
#include <map>
#undef min
#undef max
#include <set>
using namespace std;

#include "MRAGcore/MRAGrid.h"

namespace MRAG
{
	template <typename WaveletType, typename BlockType>
	class Grid_Hilbert2D : public Grid<WaveletType, BlockType>
	{	
		class HilbertIndexer2D
		{
		public:
			struct CartesianIndex { int x,y; };
			
			CartesianIndex decode(const int hilbertcode, const int level) const
			{
				int ix = 0, iy = 0;
				
				int rule0 = 0, rule1 = 1, rule2 = 3, rule3 = 2;
				
				for(int l=level-1; l>=0; l--)
				{			
					const int hilbert = (hilbertcode >> 2*l) & 0x3;
					
					const int peano = 0*(hilbert == rule0) +  1*(hilbert == rule1) +  2*(hilbert == rule2) +  3*(hilbert == rule3);
					
					ix |= (peano & 1) << l;
					iy |= ((peano & 2)>>1) << l;
					
					const int h0 = (hilbert == 0);
					const int h12 = (hilbert == 1 || hilbert == 2);
					const int h3 = (hilbert == 3);
					
					const int lut0[] = {0,3,2,1};
					const int lut3[] = {2,1,0,3};
					
					const int newrule0 = h0*lut0[rule0] + h12*rule0 +  h3*lut3[rule0];
					const int newrule1 = h0*lut0[rule1] + h12*rule1 +  h3*lut3[rule1];
					const int newrule2 = h0*lut0[rule2] + h12*rule2 +  h3*lut3[rule2];
					const int newrule3 = h0*lut0[rule3] + h12*rule3 +  h3*lut3[rule3];
					
					rule0 = newrule0;
					rule1 = newrule1;
					rule2 = newrule2;
					rule3 = newrule3;
				}
				
				CartesianIndex result;
				result.x = ix;
				result.y = iy;
				
				return result;
			}
			
			int encode(const int ix, const int iy, const int level) const
			{
				assert(ix >= 0);
				assert(ix < (1<<level));
				assert(iy >= 0);
				assert(iy < (1<<level));
				
				int encoding = 0;
				
				int rule0 = 0, rule1 = 1, rule2 = 3, rule3 = 2;
				
				for(int l=level-1; l>=0; l--)
				{
					const int peano = (ix>>l & 1) + ((iy>>l & 1) << 1);
					
					const int hilbert = rule0*(peano==0) + rule1*(peano==1) + rule2*(peano==2) + rule3*(peano==3);
					
					encoding |= (hilbert << 2*l);
					
					const int h0 = (hilbert == 0);
					const int h12 = (hilbert == 1 || hilbert == 2);
					const int h3 = (hilbert == 3);
					
					const int lut0[] = {0,3,2,1};
					const int lut3[] = {2,1,0,3};
					
					const int newrule0 = h0*lut0[rule0] + h12*rule0 +  h3*lut3[rule0];
					const int newrule1 = h0*lut0[rule1] + h12*rule1 +  h3*lut3[rule1];
					const int newrule2 = h0*lut0[rule2] + h12*rule2 +  h3*lut3[rule2];
					const int newrule3 = h0*lut0[rule3] + h12*rule3 +  h3*lut3[rule3];
					
					rule0 = newrule0;
					rule1 = newrule1;
					rule2 = newrule2;
					rule3 = newrule3;
				}
				
				return encoding;
			}	
		};
		
		vector<BlockInfo> _sort(const vector<BlockInfo>& blocksinfo) const
		{			
			const int N = blocksinfo.size();
			
			vector< pair<int, int > > tosort(N);
			
			{
				HilbertIndexer2D indexer;
				const int mylevel = this->getCurrentMaxLevel();
				
#pragma omp parallel for
				for(int i=0; i<N; i++)
				{
					const int ix = blocksinfo[i].index[0] << (mylevel-blocksinfo[i].level);
					const int iy = blocksinfo[i].index[1] << (mylevel-blocksinfo[i].level);
					
					tosort[i] = pair<int, int> (indexer.encode(ix, iy, mylevel), i);
				}
			}
			
			sort(tosort.begin(), tosort.end());
			
			vector<BlockInfo> result(N);
			
			for(int i=0; i<N; i++)
				result[i] = blocksinfo[tosort[i].second];
			
			return result;
		}
		
		static bool sort_pred(const pair<int, int> left, const pair<int, int>  right)
		{
			return left.second < right.second;
		}
		
		
	public:
		
		Grid_Hilbert2D(int nBlocksX, int nBlocksY=1, int nBlocksZ=1, BlockCollection<BlockType>* collection = NULL, bool bVerbose=true):
			Grid<WaveletType, BlockType>(nBlocksX, nBlocksY, nBlocksZ, collection, bVerbose)	{}
		
        Grid_Hilbert2D(int nBlocksX, int nBlocksY, int nBlocksZ, const int maxStencil[2][3], BlockCollection<BlockType>* collection = NULL, bool bVerbose=true):
		Grid<WaveletType, BlockType>(nBlocksX, nBlocksY, nBlocksZ, maxStencil, collection, bVerbose)	{}
		
		vector<BlockInfo> getBlocksInfo() const
		{			
			return _sort( Grid<WaveletType, BlockType>::getBlocksInfo() );
		}
		
		virtual void _refresh(bool bUpdateLazyData = false)
		{
			Grid<WaveletType, BlockType>::_refresh(bUpdateLazyData);
			
			for(int level=0; level<this->m_blockAtLevel.size(); level++)
			{
				vector<BlockInfo> tmp = this->m_blockAtLevel[level];
				this->m_blockAtLevel[level] = _sort(tmp);
			}
		}
	};
}


