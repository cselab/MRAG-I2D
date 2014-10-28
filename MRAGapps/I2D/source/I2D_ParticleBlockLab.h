/*
 *  I2D_ParticleBlockLab.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "I2D_VectorBlockLab.h"
#include "I2D_CoreParticles.h"

struct Streamer_OmegaAndVelocity
{
	template <int components>
	inline static void operate(FluidElement2D& input, Real time, Real output[components]){ abort(); }
};

template <> inline void Streamer_OmegaAndVelocity::operate<3>(FluidElement2D& input, Real time, Real output[3])
{
	output[0] = input.omega;
	output[1] = input.u[0];
	output[2] = input.u[1];
}

template<typename BlockType>
class I2D_ParticleBlockLab: public I2D_VectorBlockLab<Streamer_OmegaAndVelocity, 3>::Lab<BlockType>
{
public:
	I2D_CoreParticles pcore;
	
	void load(const BlockInfo& info)
	{
		I2D_VectorBlockLab<Streamer_OmegaAndVelocity, 3>::Lab<BlockType>::load(info);
		
		//kill vorticity at the inlet
		if(info.index[0] == 0)
		{
			for(int iy=this->m_stencilStart[1]; iy<this->m_stencilEnd[1]; iy++)
				for(int ix=this->m_stencilStart[0]; ix<0; ix++)
					this->m_sourceData[0][this->base_offset + ix + iy*this->row_size] = 0;
		}
	}
};