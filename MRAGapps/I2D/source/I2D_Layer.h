/*
 *  Layer.h
 *  VM2D
 *
 *  Created by Diego Rossinelli on 2/9/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once 
#include <math.h>
#include <string>
#include <vector>
using namespace std;
#include <assert.h>

#ifdef _WIN32
#define ALIGN_ATTRIBUTE __declspec( align( 16 ) )
//
#else
#define ALIGN_ATTRIBUTE __attribute__((aligned(16))) 
#endif

template<int _sizeX, int _sizeY, int _nDim=2>
struct ALIGN_ATTRIBUTE I2D_Layer
{
	static const int sizeX = _sizeX;
	static const int sizeY = _sizeY;
	static const int nDim = _nDim;

	Real data[nDim][sizeY][sizeX];

	template<int iChannel>
	inline  I2D_Layer<_sizeX, _sizeY, 1>& getSubLayer() // DO IT AT YOUR RISK!!!!!!!
	{
		return *(I2D_Layer<_sizeX, _sizeY, 1> *)(&data[iChannel]);
	}

	inline Real& operator ()(int ix = 0, int iy=0, int dim = 0)
	{
#ifndef NDEBUG
		assert(ix>=0 && ix<sizeX);
		assert(iy>=0 && iy<sizeY);
#endif
		
		return data[dim][iy][ix];
	}

	inline Real read(int ix = 0, int iy=0, int dim = 0) const
	{
#ifndef NDEBUG
		assert(ix>=0 && ix<sizeX);
		assert(iy>=0 && iy<sizeY);
#endif
		return data[dim][iy][ix];
	}

	inline Real read_bczero(int ix = 0, int iy = 0, int dim = 0) const
	{
		if(ix < 0 || ix >= sizeX || iy < 0 || iy >= sizeY ) return 0;
		else return data[dim][iy][ix];
	}

	const I2D_Layer& operator=(const Real val)
	{
		for(int idim = 0; idim<nDim; idim++)
		for(int iy = 0; iy<sizeY; iy++)
		for(int ix = 0; ix<sizeX; ix++)
			data[idim][iy][ix] = val;
		
		return *this;
	}

	const I2D_Layer& operator=(const I2D_Layer& l)
	{
		for(int idim = 0; idim<nDim; idim++)
			for(int iy = 0; iy<sizeY; iy++)
				for(int ix = 0; ix<sizeX; ix++)
					data[idim][iy][ix] = l.data[idim][iy][ix];
		
		return *this;
	}

	template<int dim>
	void clear(Real val)
	{
		for(int iy = 0; iy<sizeY; iy++)
		for(int ix = 0; ix<sizeX; ix++)
			data[dim][iy][ix] = val;
	}

	const vector<double> operator -(const I2D_Layer& layer)
	{
		vector<double> result;
		
		//compute linf distance
		{
			double LInf_diff = 0;
			for(int idim = 0; idim<nDim; idim++)
				for(int iy = 0; iy<sizeY; iy++)
					for(int ix = 0; ix<sizeX; ix++)
						LInf_diff = max(LInf_diff, (double)fabs(data[idim][iy][ix] - layer.data[idim][iy][ix]));
			
			result.push_back(LInf_diff);
		}
		
		//compute linf distance
		{
			double L2error = 0;
			for(int idim = 0; idim<nDim; idim++)
				for(int iy = 0; iy<sizeY; iy++)
					for(int ix = 0; ix<sizeX; ix++)
						L2error += pow(data[idim][iy][ix] - layer.data[idim][iy][ix], 2);
			
			result.push_back(sqrt((double)L2error/(sizeY*sizeX)));
		}
		
		return result;
	}

	const vector<double> computeNorms(int start_d=0, int end_d=nDim) const
	{
		vector<double> result;
		
		//compute linf distance
		{
			for(int idim = start_d; idim<end_d; idim++)
			{
				double LInf_diff = 0;
				for(int iy = 0; iy<sizeY; iy++)
					for(int ix = 0; ix<sizeX; ix++)
						LInf_diff = max(LInf_diff, (double)fabs(data[idim][iy][ix]));
				
				result.push_back(LInf_diff);
			}
			
			
		}
		
		//compute linf distance
		{
			for(int idim = start_d; idim<end_d; idim++)
			{
				double L2error = 0;
				
				for(int iy = 0; iy<sizeY; iy++)
					for(int ix = 0; ix<sizeX; ix++)
						L2error += pow(data[idim][iy][ix], 2);
				
				result.push_back(sqrtf(L2error/(sizeY*sizeX)));
			}
		}
		
		return result;
	}

	template<int iDim>
	Real * getPlane() 
	{
		return (Real*)&data[iDim][0][0];
	}

	static  Real getH0() {return 1.0/sizeX;}
	static  Real getH1() {return 1.0/sizeX;}

	void MatlabDelCazzo(string sFileName) const
	{
		FILE * f = fopen(sFileName.data(), "w");
		assert(f!=NULL);
		
		for (int iy=0; iy<sizeY; iy++)
			for (int ix=0; ix<sizeX; ix++)
			{
				for(int idim = 0; idim<nDim; idim++)
					fprintf(f, "%e\t", data[idim][iy][ix]);
				
				fprintf(f, "\n");
			}
		
		fclose(f);
	}

	void serialize(string sFileName) const
	{
		vector<double> norms = computeNorms();
		
		FILE * f = fopen(sFileName.data(),"wb");
		fwrite(&norms[0], sizeof(double), 1, f);
		fwrite(&norms[1], sizeof(double), 1, f);
		fwrite(data, sizeof(Real), nDim*sizeX*sizeY, f);
		
		printf("serialized,  norm are : Linf norm  = %e, L2 norm = %e\n", norms[0], norms[1]);
		
		fclose(f);
	}

	void deserialize(string sFileName) 
	{
		FILE * f = fopen(sFileName.data(),"rb");
		if (f==NULL)
		{
			printf("\nCould not open file to deserialize...\n"); 
			abort();
		}
		
		double n1=0, n2=0;
		fread(&n1, sizeof(double), 1, f);
		fread(&n2, sizeof(double), 1, f);
		
		fread(data, sizeof(Real), _nDim*_sizeX*_sizeY, f);
		fclose(f);
		
		/*
		cout << "closed";
		vector<double>norms = computeNorms();
		
		//printf("deserialized, checking norm errors (che palle): Linf norm error = %e, L2 norm error= %e\n", fabs(n1-norms[0]), fabs(n2-norms[1]));
		
		cout << "checking norms";	   
		if(n1 != norms[0]) abort();
		if(n2 != norms[1]) abort();
		cout << " fine\n";
		*/
	}
};


#include <xmmintrin.h>
template <typename LayerT> 
LayerT * allocate()
{
	void * data = _mm_malloc(sizeof(LayerT), 16);
	return new (data) LayerT;
}

