//
//  I2D_VelocitySolver_Wim.cpp
//  I2D_ROCKS
//
//  Created by Wim van Rees on 4/2/13.
//
//
#include <tbb/parallel_sort.h>
#include <utility>

#include "I2D_VelocitySolver_Wim.h"
#include "I2D_Clear.h"

extern  double _THETA;

#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGEnvironment.h"
#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGProfiler.h"

#include "mani-fmm2d/VortexExpansions.h"
#include "mani-fmm2d/hcfmm_box.h"
#include "mani-fmm2d/hcfmm_boxBuilder_serial.h"
#ifndef _MRAG_TBB
#include "mani-fmm2d/hcfmm_evaluator_serial.h"
#else
#include "mani-fmm2d/hcfmm_evaluator_tbb.h"
#endif


struct myExpansions: _VortexExpansions<VelocitySourceParticle, _ORDER_>
{
	void evaluateExpansions (Real *location, VelocityRHS *out_RHS, const Real theta, Real & errorBound)
	{
		const ExpansionsValueType rp = ExpansionsValueType(location[0],location[1])-ExpansionsValueType(this->Center[0],this->Center[1]);
        
		std::complex<double> csum = std::complex<double>(0,0);
		std::complex<double> prod = std::complex<double>(1,0);
		
#pragma unroll
		for (int n=0;n<_ORDER_;++n)
		{
			csum+=(prod*(std::complex<double>)this->values[0][n]);
			prod/=(std::complex<double>)rp;
		}
		
		csum*= -ExpansionsValueType(0,1)/rp;
		
		out_RHS->x[0]+=csum.real();
		out_RHS->x[1]-=csum.imag();
        
        errorBound = 1.0/std::abs(rp) * std::pow(theta,_ORDER_)/(1-theta);// prefactor A multiplication done outside
	}
};

struct Core_VelocitySolver
{
    typedef VelocitySourceParticle tParticle;
    typedef myExpansions tExpansions;

    typedef HCFMM::Box<tExpansions,_FMM_MAX_LEVEL_> tBox;
    typedef HCFMM::boxBuilder_serial<tExpansions,_FMM_MAX_LEVEL_> tBoxBuilder;
    
    Grid<W,B>& grid;
    Profiler * pProfiler;
    
    int _createSourceParticles(tParticle *& sources, const Real tolParticle, const Real scaling_factor) const
    {
        assert(sources == NULL);
        
        static const int BS = FluidBlock2D::sizeZ*FluidBlock2D::sizeY*FluidBlock2D::sizeX;
        
        vector<BlockInfo> vInfo = grid.getBlocksInfo();
        const BlockCollection<B>& coll = grid.getBlockCollection();
        
        int nof_sources = -1;
        
        //find the amount of sources
        {
            int n = 0;
            for(vector<BlockInfo>::iterator itBlock = vInfo.begin(); itBlock!=vInfo.end(); ++itBlock)
            {
                FluidElement2D * const e = &coll[itBlock->blockID](0,0,0);
                for(int i=0; i<BS; i++) n += (int)(fabs(e[i].omega) > tolParticle);
            }
            
            nof_sources = n;
        }
        
        //create an array with the sources
        if (nof_sources == 0) return 0;
        sources = new tParticle[nof_sources];
        assert(sources != NULL);
        
        {
            int c = 0;
            for(vector<BlockInfo>::iterator itBlock = vInfo.begin(); itBlock!=vInfo.end(); ++itBlock)
            {
                const Real dV = std::pow(itBlock->h[0],2);
                const Real prefac = scaling_factor*dV;
                
                B& b = coll[itBlock->blockID];
                for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
                    for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
                        if (fabs(b(ix,iy).omega) > tolParticle)
                        {
                            Real p[2];
                            itBlock->pos(p, ix, iy);
                            
                            sources[c].x[0] = p[0];
                            sources[c].x[1] = p[1];
                            sources[c].w[0] = prefac*b(ix,iy).omega;
                            
                            c++;
                        }
            }
        }
        
        return nof_sources;
    }
        
    struct VelocityEvaluatorWim
    {
        tBox * const rootNode;
        const Real inv_scaling;
        const Real theta;
        
        VelocityEvaluatorWim(tBox * rootNode, const Real inv_scaling, const Real theta):
        rootNode(rootNode), inv_scaling(inv_scaling), theta(theta)
        {
        }
        
        VelocityEvaluatorWim(const VelocityEvaluatorWim& c):
        rootNode(c.rootNode), inv_scaling(c.inv_scaling), theta(c.theta)
        {
        }
        
        bool isClosePoint(Real srcBoxCtr[2],Real srcBoxRad,const Real trgPos[2]) const
        {
            const Real p2b_dist = std::sqrt((srcBoxCtr[0]-trgPos[0])*(srcBoxCtr[0]-trgPos[0]) +
                                            (srcBoxCtr[1]-trgPos[1])*(srcBoxCtr[1]-trgPos[1]));
            
            const Real denom = std::max(p2b_dist, std::numeric_limits<Real>::epsilon());
            return (srcBoxRad>=theta*denom);
        }
        
        bool isIntersectingPoint(Real srcBoxCtr[2],Real srcBoxRad, const Real trgBoxCtr[2]) const
        {
            const Real intersection[2] = {
                min(trgBoxCtr[0], srcBoxCtr[0]+srcBoxRad) - max(trgBoxCtr[0],srcBoxCtr[0]-srcBoxRad),
                min(trgBoxCtr[1], srcBoxCtr[1]+srcBoxRad) - max(trgBoxCtr[1],srcBoxCtr[1]-srcBoxRad)
            };
            
            return intersection[0]>=0 && intersection[1]>=0;
        }
        
        bool isCloseBox(Real srcBoxCtr[2],Real srcBoxRad, const Real trgBoxCtr[2],const Real halfTrgBoxWidth) const
        {
            
            // minimum distance given that target box is a cube
            const Real x0 = trgBoxCtr[0]-halfTrgBoxWidth;
            const Real x1 = trgBoxCtr[0]+halfTrgBoxWidth;
            const Real y0 = trgBoxCtr[1]-halfTrgBoxWidth;
            const Real y1 = trgBoxCtr[1]+halfTrgBoxWidth;
            const Real x = srcBoxCtr[0];
            const Real y = srcBoxCtr[1];
            
            const Real dx = x < x0 ? (x0-x) : (x > x1 ? (x-x1) : 0); // can replace this with min/max
            const Real dy = y < y0 ? (y0-y) : (y > y1 ? (y-y1) : 0);
            const Real min_p2bDist = std::sqrt(dx*dx+dy*dy);
            
            const Real denom = std::max(min_p2bDist, std::numeric_limits<Real>::epsilon());
            return (srcBoxRad>=theta*denom);
        }
        
        bool isIntersectingBox(Real srcBoxCtr[2],Real srcBoxRad, const Real trgBoxCtr[2],const Real halfTrgBoxWidth) const
        {
            const Real intersection[3] = {
                std::min(trgBoxCtr[0]+halfTrgBoxWidth, srcBoxCtr[0]+srcBoxRad) - std::max(trgBoxCtr[0]-halfTrgBoxWidth,srcBoxCtr[0]-srcBoxRad),
                std::min(trgBoxCtr[1]+halfTrgBoxWidth, srcBoxCtr[1]+srcBoxRad) - std::max(trgBoxCtr[1]-halfTrgBoxWidth,srcBoxCtr[1]-srcBoxRad),
            };
            
            return intersection[0]>=0 && intersection[1]>=0;
        }
        
        std::pair<Real,Real> evaluateDirect(tBox * srcBox, const Real trgPos[2]) const
        {
            Real u=0,v=0;
            const tParticle * const p = srcBox->vparticles;
            for(int i=0;i<srcBox->nParticles;i++)
            {
                const Real distX = trgPos[0]-p[i].x[0];
                const Real distY = trgPos[1]-p[i].x[1];
                const Real distSq = distX*distX + distY*distY;
                if(distSq>std::numeric_limits<Real>::epsilon())
                {
                    const Real invDistSq = 1./distSq;
                    u-=invDistSq*distY*p[i].w[0];
                    v+=invDistSq*distX*p[i].w[0];
                }
            }
            return std::make_pair(u,v);
        }
        
//            Real evaluateDirectAlgebraic(tBox * srcBox, const Real trgPos[3], const Real eps_in=-1) const
//            {
//                Real phi=0;
//                const tParticle * const p = srcBox->vparticles;
//                const Real epsSq = eps_in > 0 ? eps_in*eps_in : EpsilonSmoothKernel::epsSq;
//                const Real threeHalfepsSq = 1.5*epsSq;
//                
//                for(int i=0;i<srcBox->nParticles;i++)
//                {
//                    const Real distSq=( (p[i].x[0]-trgPos[0])*(p[i].x[0]-trgPos[0]) +
//                                       (p[i].x[1]-trgPos[1])*(p[i].x[1]-trgPos[1]) +
//                                       (p[i].x[2]-trgPos[2])*(p[i].x[2]-trgPos[2]));
//                    phi+=p[i].w[0]*(distSq + threeHalfepsSq)/((distSq + epsSq)*std::sqrt(distSq + epsSq));
//                }
//                return phi;
//            }
        
        //LOOPC-STYLE:
        template<typename tBlockType>
        void operator()(const BlockInfo& info, tBlockType & b) const
        {
            //clean potential field:
            {
                FluidElement2D * const e = &b(0,0,0);
                
                static const int n = FluidBlock2D::sizeZ*FluidBlock2D::sizeY*FluidBlock2D::sizeX;
                for(int i=0; i<n; i++)
                {
                    e[i].u[0]=e[i].u[1]=0;
                    e[i].tmp = 0; // we will store the error here
                }
            }
            
            // traverse per particle
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    Real target_pos[2];
                    info.pos(target_pos, ix, iy);
                    
                    HCFMM::BoxIterator<tBox,tbb::scalable_allocator> srcBox(rootNode);
                    while(srcBox!=NULL && (srcBox->nParticles>0))
                    {
                        bool canRemove=true;
                        
                        const bool isClose = isClosePoint(srcBox->expansions.Center,srcBox->expansions.Radius,target_pos);
                        const bool isIntersecting = isIntersectingPoint(srcBox->expansions.Center,srcBox->expansions.Radius,target_pos);
                        if(isClose || isIntersecting)
                        {
                            if(not srcBox->isleaf)
                                canRemove=false;
                            else
                            {
                                std::pair<Real,Real> vel = evaluateDirect(srcBox,target_pos);
                                b(ix, iy).u[0]+=vel.first;
                                b(ix, iy).u[1]+=vel.second;
                            }
                        }
                        else
                        {
                            
                            VelocityRHS rhs;
                            Real errorBound = 0.0;
                            srcBox->expansions.evaluateExpansions(target_pos, &rhs, theta, errorBound);
                            b(ix, iy).u[0] += rhs.x[0];
                            b(ix, iy).u[1] += rhs.x[1];
                            b(ix, iy).tmp += errorBound*srcBox->TotalMass/(2.0*M_PI)*inv_scaling;
                        }
                        
                        if(canRemove)
                            srcBox.advanceRemove();
                        else
                            srcBox++;
                    }
                }

            //multiply by scaling factor
            {
                FluidElement2D * const e = &b(0,0,0);
                
                static const int n = FluidBlock2D::sizeZ*FluidBlock2D::sizeY*FluidBlock2D::sizeX;
                
                const Real scale = 1./(2.0*M_PI)*inv_scaling;
                
                for(int i=0; i<n; i++)
                {
                    e[i].u[0] *= scale;
                    e[i].u[1] *= scale;
                }
            }            
        }
    };
    
public:
	
	Core_VelocitySolver(Grid<W,B>& grid):
	grid(grid)
	{}
	
	int execute(BlockProcessing& block_processing, double tolParticle, double in_scaling_factor, const Real theta)
	{
        Profiler profiler;
        
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		const BlockCollection<B>& coll = grid.getBlockCollection();
		
		VelocitySourceParticle * sources = NULL;
		const int nof_sources = _createSourceParticles(sources,tolParticle, in_scaling_factor);
        
		if (nof_sources == 0)
		{
			printf("no particles: NEED TO CLEAR\n");
//			ClearPsi clear;
//			block_processing.process(vInfo, coll ,clear);
			
			return 0;
		}
		
		assert(nof_sources > 0);
		assert(sources != NULL);
		
        
        tBox * rootBox=new tBox;
        
        profiler.push_start("tree");
        tBoxBuilder::buildBoxes(sources, nof_sources, rootBox);
        profiler.pop_stop();
        
        profiler.push_start("expansions");
        tBoxBuilder::generateExpansions(rootBox);
        profiler.pop_stop();
        
        profiler.push_start("evaluations");
		VelocityEvaluatorWim evaluator(rootBox, 1./in_scaling_factor,theta);
		block_processing.process(vInfo, coll ,evaluator);
        profiler.pop_stop();
        
        profiler.printSummary();
        
		delete [] sources;sources=NULL;
		delete rootBox;rootBox=NULL;
        
		return nof_sources;
	}
    
};


void I2D_VelocitySolver_Wim::compute_velocity()
{
	I2D_Clear cleaner;
	cleaner.clearVel(*grid_ptr);
    
	vector<BlockInfo> vInfo = grid_ptr->getBlocksInfo();
	
    Core_VelocitySolver core(*grid_ptr);
    const long long nsources = core.execute(block_processing, tolParticle, scaling_factor,theta);
    
    printf("nsources: %lld\n",nsources);    
}
