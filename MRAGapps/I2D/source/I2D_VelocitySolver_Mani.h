/*
 *  I2D_VelocitySolver_Mani.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 *	This operator takes omega[0-2] and compute the velocity
 *	in u[0-2]. It allocates and uses external memory. Tested!
 *
 *	IN:		omega[0-2]
 *	OUT:	velocity[0-2]
 *
 */

#pragma once

#include <limits>

#include "I2D_VelocityOperator.h"
#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_FMMTypes.h"

#include "I2D_CoreFMM_AggressiveVel.h"
#include "I2D_CoreFMM_SSE.h"
#include "I2D_CoreFMM_Check.h"

class I2D_VelocitySolver_Mani: public I2D_VelocityOperator
{
protected:

	template<typename streamer, int stage>
	struct ThresholdParticles
	{
		const Real tolParticle;
		const Real scaling_factor;
		VelocitySourceParticle * destptr;
		map<int, SourceParticlesInfo>& blockid2info;

		ThresholdParticles(Real _tolParticle, Real _scaling_factor, map<int, SourceParticlesInfo>& blockid2info): tolParticle(_tolParticle), scaling_factor(_scaling_factor), blockid2info(blockid2info), destptr(NULL){}
		ThresholdParticles(Real _tolParticle, Real _scaling_factor, const ThresholdParticles& c): tolParticle(_tolParticle), scaling_factor(_scaling_factor), blockid2info(c.blockid2info), destptr(c.destptr){}

		inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
		{

			map<int, SourceParticlesInfo>::iterator it = blockid2info.find(info.blockID);
			assert(it != blockid2info.end());
			SourceParticlesInfo& result = it->second;

			if (stage == 0)
			{
				int n = 0;

				FluidElement2D * const e = &b(0,0);
				static const int BS = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
				for(int i=0; i<BS; i++)
					n += (int)(fabs(streamer::stream(e[i])) > tolParticle);

				result.nsource_particles = n;

				//for(int i=0; i<BS; i++)
				//{
				//	e[i].u[0] = 0;
				//	e[i].u[1] = 0;
				//}
			}
			else if (stage == 1)
			{
				const Real dV = pow(info.h[0],2);
				const Real prefac = scaling_factor*dV;

				VelocitySourceParticle * sources = destptr + result.start;

				FluidElement2D * const e = &b(0,0);

				for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
					for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
					{
						const int idx = iy*FluidBlock2D::sizeX+ix;
						const Real variable = streamer::stream(e[idx]);

						if(fabs(variable) > tolParticle)
						{
							Real p[2];
							info.pos(p,ix,iy);

							sources->x[0] = p[0];
							sources->x[1] = p[1];
							sources->w[0] = prefac*variable;

							sources++;
						}
					}

				assert(sources == destptr + result.start + result.nsource_particles);
			}
			else abort();
		}
	};

	Grid<W,B>* grid_ptr;
	BlockProcessing block_processing;
	I2D_CoreFMM_AggressiveVel * coreFMM;
	bool bSKIPBLOCKS;

	Real theta;
	Real tolParticle;
	Real scaling_factor;

	int nsource_particles;
	int lmax;	

	map<int, SourceParticlesInfo> blockid2info;
	vector<BlockInfo> vDest;
	VelocitySourceParticle * srcparticles;
	VelocityBlock * my_velBlocks;

	virtual void _count_sourceparticles();
	virtual void _collect_sourceparticles();
	void _compute();
	virtual void _updateBlocks();
	void _cleanup();

	void _setup(ArgumentParser& parser)
	{
		parser.unset_strict_mode();

		if (parser("-core-fmm").asString() == "sse")
			coreFMM = new I2D_CoreFMM_SSE();
		else if (parser("-core-fmm").asString() == "check")
			coreFMM = new I2D_CoreFMM_Check ();
		else 
			coreFMM = new I2D_CoreFMM_AggressiveVel();

		bSKIPBLOCKS = parser("-fmm-skip").asBool(); 

		parser.set_strict_mode();

		lmax = parser("-lmax").asDouble();
		const double min_dV = pow(pow(0.5,lmax)/B::sizeX, 2);
		scaling_factor = max(1., 1e-4/min_dV);
		tolParticle = 2*numeric_limits<Real>::epsilon();
		theta = parser("-fmm-theta").asDouble();
		printf("scaling_factor=%f tolParticle=%f\n", scaling_factor, tolParticle);

		assert(theta>0);
		assert(coreFMM != NULL);
	}


	I2D_VelocitySolver_Mani(const int argc, const char ** argv):
		srcparticles(NULL), nsource_particles(0), my_velBlocks(NULL), grid_ptr(NULL)
	{
		ArgumentParser parser(argc, argv);

		_setup(parser);
	}

	void set_grid_ptr(Grid<W,B>* ptr)
	{
		grid_ptr = ptr;
	}

public:

	I2D_VelocitySolver_Mani(Grid<W,B>& grid, ArgumentParser& parser_): 
		srcparticles(NULL), nsource_particles(0), my_velBlocks(NULL), grid_ptr(&grid)
	{
		_setup(parser_);
	}

	virtual void compute_velocity();
};
