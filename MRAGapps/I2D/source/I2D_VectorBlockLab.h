/*
 *  I2D_VectorBlockLab.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <xmmintrin.h>

template<typename Streamer, int components>
struct I2D_VectorBlockLab
{
public:
	template<typename BlockType>
	class Lab
	{
	protected:
		typedef typename BlockType::ElementType ElementType;
		
		enum eBlockLab_State {Prepared, ReadyToLoad, Loaded, Uninitialized};
		eBlockLab_State m_state;
		
		Real * m_weightPool, (*m_valuePool)[components], * m_sourceData[components];
		Real target_time;
		int m_weightPoolSize, m_valuePoolSize;
		int m_stencilStart[2], m_stencilEnd[2];
		int size[2];
		int row_size, base_offset;
		
		const BlockCollection<BlockType>* m_refCollection;
		const BoundaryInfo* m_refBoundaryInfo;
		
		template<typename T> inline _MRAG_BLOCKLAB_ALLOCATOR<T> allocator() const 
		{ 
			return _MRAG_BLOCKLAB_ALLOCATOR<T>();
		};
		
		template<typename T>
		void _release_vector(T *& t, int n)
		{ 
			if (t != NULL)
				allocator<T>().deallocate(t, n); 
			
			t = NULL; 
		}
		
	public:
		
		Lab():
		m_state(Uninitialized),
		m_weightPool(NULL), m_weightPoolSize(0),
		m_valuePool(NULL), m_valuePoolSize(0),
		row_size(0),  base_offset(0),
		m_refCollection(NULL), m_refBoundaryInfo(NULL)
		{
			m_stencilStart[0] = m_stencilStart[1] = 0;
			m_stencilEnd[0] = m_stencilEnd[1] = 0;
			
			for(int c=0; c<components; c++)
				m_sourceData[c] = NULL;
		}
		
		~Lab()
		{
			for(int c=0; c<components; c++)
				_release_vector(m_sourceData[c], size[0]*size[1]);
			
			_release_vector(m_weightPool, m_weightPoolSize);
			_release_vector(m_valuePool, m_valuePoolSize);
		}
		
		void prepare(const BlockCollection<BlockType>& collection, const BoundaryInfo& boundaryInfo, const int stencil_start[3], const int stencil_end[3])
		{
			m_refCollection = &collection;
			m_refBoundaryInfo = &boundaryInfo;
			
			assert(stencil_start[0]>= boundaryInfo.stencil_start[0]);
			assert(stencil_start[1]>= boundaryInfo.stencil_start[1]);
			assert(stencil_start[2]>= boundaryInfo.stencil_start[2]);
			assert(stencil_start[2] == 0);
			assert(stencil_end[0]<= boundaryInfo.stencil_end[0]);
			assert(stencil_end[1]<= boundaryInfo.stencil_end[1]);
			assert(stencil_end[2]<= boundaryInfo.stencil_end[2]);
			assert(stencil_end[2] == 1);
			
			m_stencilStart[0] = boundaryInfo.stencil_start[0];
			m_stencilStart[1] = boundaryInfo.stencil_start[1];
			
			m_stencilEnd[0] = boundaryInfo.stencil_end[0];
			m_stencilEnd[1] = boundaryInfo.stencil_end[1];
			
			assert(m_stencilStart[0]<=m_stencilEnd[0]);
			assert(m_stencilStart[1]<=m_stencilEnd[1]);
			
			if (m_sourceData[0] == NULL || 
				size[0]!= BlockType::sizeX + m_stencilEnd[0] - m_stencilStart[0] -1 ||
				size[1]!= BlockType::sizeY + m_stencilEnd[1] - m_stencilStart[1] -1 )			{
				
				if (m_sourceData != NULL)
					for(int c=0; c<components; c++)
						_release_vector(m_sourceData[c], size[0]*size[1]);
				
				size[0] = BlockType::sizeX + m_stencilEnd[0] - m_stencilStart[0] -1;
				size[1] = BlockType::sizeY + m_stencilEnd[1] - m_stencilStart[1] -1;
				
				row_size = size[0];
				base_offset = - m_stencilStart[0] - m_stencilStart[1]*row_size;
				
				for(int c=0; c<components; c++)
				{
					m_sourceData[c] = allocator<Real>().allocate(size[0]*size[1]);
					assert((((unsigned long int)m_sourceData[c]) & 0xf) ==  0);
				}
			}
			
			m_state = Prepared;
		}
		
		template<typename Processing>
		void inspect(Processing& p)
		{
			assert(m_state == Prepared || m_state==Loaded);
			
			target_time = p.t;
			
			m_state = ReadyToLoad;
		}
		
		void load(const BlockInfo& info)
		{
			const BlockCollection<BlockType>& collection = *m_refCollection;
			const BoundaryInfo& boundaryInfo = *m_refBoundaryInfo;
			
			//0. couple of checks
			//1. load the block into the cache
			//2. create the point and the weight pools
			//3. compute the ghosts, put them into the cache
			
			//0.
			assert(m_state == ReadyToLoad || m_state == Loaded);
			assert(m_sourceData[0] != NULL);
			
			const int nX = BlockType::sizeX;
			const int nY = BlockType::sizeY;
			const int nZ = BlockType::sizeZ;
			
			BoundaryInfoBlock& bbinfo = *boundaryInfo.boundaryInfoOfBlock.find(info.blockID)->second;
			bbinfo.lock();
			
			const vector< vector<IndexWP> >& ghosts = bbinfo.getGhosts();
			
			//1.
			{
				collection.lock(info.blockID);
				
				BlockType& block = collection[info.blockID];
				
				for(int iz=0; iz<nZ; iz++)
					for(int iy=0; iy<nY; iy++)
					{
						Real * ptrDestination[components];
						
						for(int c=0; c<components; c++)
							ptrDestination[c] = m_sourceData[c] + base_offset + iy*row_size;
						
						ElementType * ptrSource = &block(0, iy, iz);
						
						for(int ix=0; ix<nX; ix++)
						{
							Real output[components];
							Streamer::template operate<components>(ptrSource[ix], target_time, output);
							
							for(int c=0; c<components; c++)
								ptrDestination[c][ix] = output[c];
						}
					}
				
				collection.release(info.blockID);
			}
			
			//2.
			{
				const int nValuePoolSize = bbinfo.getIndexPool().size();
				
				if (nValuePoolSize > m_valuePoolSize)
				{
					_release_vector(m_valuePool, m_valuePoolSize);
					
					m_valuePool = allocator<Real[components]>().allocate(nValuePoolSize);
					
					m_valuePoolSize = nValuePoolSize;
				}
				
				collection.lock(bbinfo.dependentBlockIDs);
				
				/* here we are going to do this:
				 for(int i=0; i<nValuePoolSize; i++)
				 m_valuePool[i] = collection[bbinfo.getIndexPool()[i].blockID][bbinfo.getIndexPool()[i].index];
				 which is simply rewritten as:
				 */
				
				{  
					const PointIndex  * p = &(bbinfo.getIndexPool()[0]);
					Real (* e)[components] = &(m_valuePool[0]);
					BlockType * b = NULL;
					int oldBlockId = -1;
					
					for(int i=0;i<nValuePoolSize;++i){
						if(oldBlockId!=(*p).blockID){
							oldBlockId = (*p).blockID;
							b = &(collection[oldBlockId]);
						}
						
						Streamer::template operate<components>((*b)[(*p).index], target_time, *e);
						//(*b)[(*p).index].template time_rec<component>(target_time);
						++e; ++p;
					}             
				}
				
				collection.release(bbinfo.dependentBlockIDs);
				
				const int nWeightPoolSize = bbinfo.weightsPool.size();
				if (nWeightPoolSize > m_weightPoolSize)
				{
					_release_vector(m_weightPool, m_weightPoolSize);
					
					m_weightPoolSize = nWeightPoolSize;
					m_weightPool = allocator<Real>().allocate(nWeightPoolSize);
				}
				
				for(int i=0; i<nWeightPoolSize; i++)
					m_weightPool[i] = bbinfo.weightsPool[i];
			}
			
			//3.
			{
				for(int icode2D=0; icode2D<9; icode2D++)
				{
					if (icode2D == 1*1 + 3*1) continue;
					
					const int code[3] = { icode2D%3-1, (icode2D/3)%3-1, 0};
					
					const int s[3] = { 
						code[0]<1? (code[0]<0 ? m_stencilStart[0]:0 ) : nX, 
						code[1]<1? (code[1]<0 ? m_stencilStart[1]:0 ) : nY, 
						0 };
					
					const int e[3] = {
						code[0]<1? (code[0]<0 ? 0:nX ) : nX+m_stencilEnd[0]-1, 
						code[1]<1? (code[1]<0 ? 0:nY ) : nY+m_stencilEnd[1]-1, 
						1 };
					
					
					assert(bbinfo.boundary[icode2D + 9*1].nGhosts == (e[2]-s[2])*(e[1]-s[1])*(e[0]-s[0]));
					
					int currentghost = bbinfo.boundary[icode2D + 9*1].start;
					
					for(int iy=s[1]; iy<e[1]; iy++)
						for(int ix=s[0]; ix<e[0]; ix++, currentghost++)
						{
							const vector<IndexWP>& vWP= ghosts[currentghost];
							
							const vector<IndexWP>::const_iterator start = vWP.begin();
							const vector<IndexWP>::const_iterator end = vWP.end();
							
							Real ghost[components];
							for(int c=0; c<components; c++)
								ghost[c] = 0;
							
							for(vector<IndexWP>::const_iterator it=start; it!= end; it++)
							{
								const Real w = (m_weightPool[it->weights_index[0]]*m_weightPool[it->weights_index[1]]);
								
								for(int c=0; c<components; c++)
									ghost[c] += w*m_valuePool[it->point_index][c];
							}
							
							for(int c=0; c<components; c++)
								m_sourceData[c][base_offset + ix + iy*row_size] = ghost[c];
						}
				}
			}
			
			bbinfo.release();
			
			m_state = Loaded;
		}
		
		template<int component>
		Real get(int ix, int iy=0) const
		{
			assert(m_state == Loaded);
			
			assert(ix-m_stencilStart[0]>=0 && ix-m_stencilStart[0]<size[0]);
			assert(iy-m_stencilStart[1]>=0 && iy-m_stencilStart[1]<size[1]);
			
			return m_sourceData[component][base_offset + ix + iy*row_size];
		}
		
	private:
		
		//forbidden
		Lab(const Lab&):
		m_state(Uninitialized),
		m_sourceData(NULL),
		m_weightPool(NULL), m_weightPoolSize(0),
		m_valuePool(NULL), m_valuePoolSize(0) {abort();}
		
		Lab& operator=(const Lab&){abort(); return *this;}	
	};
};