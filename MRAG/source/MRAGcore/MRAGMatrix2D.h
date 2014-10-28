/*
 *  Matrix3D.h
 *  ReactionDiffusion
 *
 *  Created by Diego Rossinelli on 10/19/06.
 *  Copyright 2006 ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <iostream>
using namespace std;

namespace MRAG
{
	//lo metto qui perche' mi fa schifo
#define _FORBID_COPIES(x) \
x(const x& culo){abort();}\
x& operator=(const x& culatone){abort(); return *this;}
	
	template <class DataType>  class Matrix2D
	{
	private:
		DataType * m_pData;
		unsigned int m_vSize[2];
		unsigned int m_nElements;
		bool m_bDataOwner;
	public:
		_FORBID_COPIES(Matrix2D)
		
		Matrix2D(unsigned int nSizeX, unsigned int nSizeY, DataType * data=NULL);
		~Matrix2D() { if (m_bDataOwner) delete [] m_pData; }
		
		DataType & Access(unsigned int ix, unsigned int iy) const;
		DataType & LinAccess(unsigned int i) const;
		inline DataType & operator()(unsigned int ix, unsigned int iy) const
		{
			assert(ix<m_vSize[0]);
			assert(iy<m_vSize[1]);
			
			return m_pData[iy*m_vSize[0] + ix];
		}
		
		unsigned int * getSize();
		unsigned int getNumberOfElements() const;
		DataType * getData();
		
		Matrix2D& operator= (DataType val);
		Matrix2D& operator*= (DataType val);
	};
	
	template <class DataType>  Matrix2D<DataType>::Matrix2D(unsigned int nSizeX, unsigned int nSizeY, DataType * data):
	m_pData(NULL),
	m_nElements(0),
	m_bDataOwner(data==NULL)
	{
		m_vSize[0] = nSizeX;
		m_vSize[1] = nSizeY;
		
		m_nElements = nSizeX*nSizeY;
		
		if (m_bDataOwner)
			m_pData = new DataType[m_nElements];
		else
			m_pData = data;
	}
	
	template <class DataType> inline DataType & Matrix2D<DataType>::Access(unsigned int ix, unsigned int iy) const
	{	
#ifndef NDEBUG
		assert(ix<m_vSize[0]);
		assert(iy<m_vSize[1]);
#endif
		return m_pData[ iy*m_vSize[0] + ix];
	}
	
	template <class DataType> inline DataType & Matrix2D<DataType>::LinAccess(unsigned int i) const
	{
		assert(i<m_nElements);
		return m_pData[i];
	}
	
	template <class DataType> inline unsigned int  Matrix2D<DataType>::getNumberOfElements() const
	{
		return m_nElements;
	}
	
	template <class DataType> inline unsigned int *  Matrix2D<DataType>::getSize() 
	{
		return m_vSize;
	}
	
	template <class DataType> inline Matrix2D<DataType>&  Matrix2D<DataType>::operator= (DataType val)
	{
		for(unsigned int i=0;i<m_nElements;i++) m_pData[i] = val;
		return *this;
	}
	
	template <class DataType> inline Matrix2D<DataType>&  Matrix2D<DataType>::operator*= (DataType val)
	{
		for(unsigned int i=0;i<m_nElements;i++) m_pData[i] *= val;
		return *this;
	}
	
	template <class DataType> DataType * Matrix2D<DataType>::getData ()
	{
		return m_pData;
	}
	
	/*template <class DataType> DataType Matrix2D<DataType>::getData ()
	 {
	 return m_pData;
	 }*/
	
}

