#include "I2D_SolverOperator.h"
#include "I2D_ScalarBlockLab.h"
#include "I2D_DivGradOfScalar.h"
#include "I2D_GradOfVector.h"
#include "I2D_VectorBlockLab.h"

struct Solver_DiffusionRHS_4thOrder: I2D_DivGradOfScalar_4thOrder
{
	Real t;
	Real viscosity;
	int stencil_start[3], stencil_end[3];
	
	Solver_DiffusionRHS_4thOrder(Real viscosity): viscosity(viscosity), t(0)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	Solver_DiffusionRHS_4thOrder(const Solver_DiffusionRHS_4thOrder & copy): viscosity(copy.viscosity), t(0)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	struct AddToTmp { static inline void stream(Real& out, Real in) { out += in; } };
	
	template<typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const 
	{
		//clear all tmps
		Real * const e = &out.external_data[0][0];
		
		static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
		
		for(int i=0; i<n; i++) 
			e[i] = 0;
		
		_ddpsiddx_ptr< AddToTmp >(lab, info, e);
		_ddpsiddy_ptr< AddToTmp >(lab, info, e);
		
		for(int i=0; i<n; i++)
			e[i] *= viscosity;
	}
};

struct RHSUpwind5thOrder
{
	Real t;
	Real Uinf[2];
	
	int stencil_start[3], stencil_end[3];
	
	RHSUpwind5thOrder(Real Uinf[2]): t(0)
	{
		stencil_start[0] = stencil_start[1] = -3;
		stencil_end[0] = stencil_end[1] = +4;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
		
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
	}
	
	RHSUpwind5thOrder(const RHSUpwind5thOrder & copy): t(0)
	{
		stencil_start[0] = stencil_start[1] = -3;
		stencil_end[0] = stencil_end[1] = +4;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
		
		Uinf[0] = copy.Uinf[0];
		Uinf[1] = copy.Uinf[1];
	}
	
	template<typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const 
	{
		Real dwdx[2] = {0.0,0.0};
		Real dwdy[2] = {0.0,0.0};
		Real u[2] = {0.0,0.0};
		const Real factor = 1./(60.0*info.h[0]);
		
		
		for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
			for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
			{
				dwdx[0] = -2.0*lab(ix-3, iy) + 15.0*lab(ix-2, iy) - 60*lab(ix-1, iy) + 20.0*lab(ix, iy) + 30.0*lab(ix+1, iy) - 3.0*lab(ix+2, iy);
				dwdx[1] =  2.0*lab(ix+3, iy) - 15.0*lab(ix+2, iy) + 60*lab(ix+1, iy) - 20.0*lab(ix, iy) - 30.0*lab(ix-1, iy) + 3.0*lab(ix-2, iy);
				
				dwdy[0] = -2.0*lab(ix, iy-3) + 15.0*lab(ix, iy-2) - 60*lab(ix, iy-1) + 20.0*lab(ix, iy) + 30.0*lab(ix, iy+1) - 3.0*lab(ix, iy+2);
				dwdy[1] =  2.0*lab(ix, iy+3) - 15.0*lab(ix, iy+2) + 60*lab(ix, iy+1) - 20.0*lab(ix, iy) - 30.0*lab(ix, iy-1) + 3.0*lab(ix, iy-2);
				
				u[0] = Uinf[0] + out(ix, iy).u[0];
				u[1] = Uinf[1] + out(ix, iy).u[1];
				
				out.external_data[iy][ix] -= factor*(max(u[0], (Real)0)*dwdx[0] + min(u[0], (Real)0)*dwdx[1] + 
													max(u[1], (Real)0)*dwdy[0] + min(u[1], (Real)0)*dwdy[1]);
			}		
		
		//correct for the x-inlet
		if(info.index[0] == 0)
		{
			//overwrite the rhs with the correct values for ix=0 and ix=1 and ix=2
			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
				for(int ix=0; ix<3; ix++)
				{
					dwdx[0] = -2.0*lab(ix-3, iy) + 15.0*lab(ix-2, iy) - 60*lab(ix-1, iy) + 20.0*lab(ix, iy) + 30.0*lab(ix+1, iy) - 3.0*lab(ix+2, iy);
					dwdx[1] =  2.0*lab(ix+3, iy) - 15.0*lab(ix+2, iy) + 60*lab(ix+1, iy) - 20.0*lab(ix, iy) - 30.0*lab(ix-1, iy) + 3.0*lab(ix-2, iy);
					
					dwdy[0] = -2.0*lab(ix, iy-3) + 15.0*lab(ix, iy-2) - 60*lab(ix, iy-1) + 20.0*lab(ix, iy) + 30.0*lab(ix, iy+1) - 3.0*lab(ix, iy+2);
					dwdy[1] =  2.0*lab(ix, iy+3) - 15.0*lab(ix, iy+2) + 60*lab(ix, iy+1) - 20.0*lab(ix, iy) - 30.0*lab(ix, iy-1) + 3.0*lab(ix, iy-2);
					
					u[0] = Uinf[0] + out(ix, iy).u[0];
					u[1] = Uinf[1] + out(ix, iy).u[1];
					
					out.external_data[iy][ix] += factor*(max(u[0], (Real)0)*dwdx[0] + min(u[0], (Real)0)*dwdx[1] + 
														 max(u[1], (Real)0)*dwdy[0] + min(u[1], (Real)0)*dwdy[1]);
					
					dwdx[0] = 20.0*lab(ix, iy) + 30.0*lab(ix+1, iy) - 3.0*lab(ix+2, iy);
					dwdx[1] = 2.0*lab(ix+3, iy) - 15.0*lab(ix+2, iy) + 60*lab(ix+1, iy) - 20.0*lab(ix, iy);
					
					dwdy[0] = -2.0*lab(ix, iy-3) + 15.0*lab(ix, iy-2) - 60*lab(ix, iy-1) + 20.0*lab(ix, iy) + 30.0*lab(ix, iy+1) - 3.0*lab(ix, iy+2);
					dwdy[1] =  2.0*lab(ix, iy+3) - 15.0*lab(ix, iy+2) + 60*lab(ix, iy+1) - 20.0*lab(ix, iy) - 30.0*lab(ix, iy-1) + 3.0*lab(ix, iy-2);
					
					out.external_data[iy][ix] -= factor*(max(u[0], (Real)0)*dwdx[0] + min(u[0], (Real)0)*dwdx[1] + 
														 max(u[1], (Real)0)*dwdy[0] + min(u[1], (Real)0)*dwdy[1]);
				}		 
		}		
	}
};

struct Solver_Penalization
{
	Real lambda;
	Real Uinf[2];
	
	Solver_Penalization(Real lambda, const Real Uinf[2]): lambda(lambda)
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
	}
	
	Solver_Penalization(const Solver_Penalization& c): lambda(c.lambda)
	{
		Uinf[0] = c.Uinf[0];
		Uinf[1] = c.Uinf[1];
	}
	
	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{		
		FluidElement2D * const e = &b(0,0);
		
		Real lambdaXs = 0.0;
		Real u_current[2] = {0.0,0.0};
		
		static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
		for(int i=0; i<n; i++)			
		{
			lambdaXs  = lambda*e[i].tmp;
		
			u_current[0] = e[i].u[0];
			u_current[1] = e[i].u[1];

			e[i].u[0] = -lambdaXs*(u_current[0] + Uinf[0]);
			e[i].u[1] = -lambdaXs*(u_current[1] + Uinf[1]);
			//e[i].u[0] = -lambdaXs*(Uinf[0]);
			//e[i].u[1] = -lambdaXs*(Uinf[1]);
		}
	}
};

struct Solver_CurlDiff_4thOrder: I2D_GradOfVector_4thOrder
{
	Real t;
	int stencil_start[3], stencil_end[3];
	
	Solver_CurlDiff_4thOrder(): t(0)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	Solver_CurlDiff_4thOrder(const Solver_CurlDiff_4thOrder & copy): t(0)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	struct AddToTmp { static inline void stream(Real& out, Real in) { out += in; } };
	
	struct SubToTmp { static inline void stream(Real& out, Real in) { out -= in; } };
	
	template<typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const 
	{
		Real * const e = &out.external_data[0][0];
		
		_dfdx_ptr<AddToTmp, 1 >(lab, info, e);
		_dfdy_ptr<SubToTmp, 0 >(lab, info, e);
	}
};

double I2D_SolverOperator_4thOrder::estimate_largest_dt() const
{
	double max_dx = (1./B::sizeX)*pow(0.5,grid.getCurrentMaxLevel());
	
	return pow(max_dx,2)/(viscosity * 4);
}

void I2D_SolverOperator_4thOrder::perform_timestep(double dt)
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	BoundaryInfo& binfo=grid.getBoundaryInfo();  
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	const Real lambda = 1e4;
	Real Uinf[2] = {0.1,0.0};
	
	Solver_DiffusionRHS_4thOrder diff(viscosity);
	RHSUpwind5thOrder advection(Uinf);
	Solver_Penalization penalization(lambda, Uinf);
	Solver_CurlDiff_4thOrder curl_diff;
		
	block_processing.process< I2D_ScalarBlockLab<Streamer_Omega>::Lab >(vInfo, coll, binfo, diff);
	block_processing.process< I2D_ScalarBlockLab<Streamer_Omega>::Lab >(vInfo, coll, binfo, advection);
	//block_processing.process(vInfo, coll, penalization);
	//block_processing.process< I2D_VectorBlockLab< Streamer_Velocity, 2 >::Lab >(vInfo, coll, binfo, curl_diff);	
	UpdateScalarRK2_Simple<1> stepA(dt);
	block_processing.process(vInfo, coll, stepA);
	
	//velsolver->compute_velocity();
	
	block_processing.process< I2D_ScalarBlockLab<Streamer_Omega>::Lab >(vInfo, coll, binfo, diff);
	block_processing.process< I2D_ScalarBlockLab<Streamer_Omega>::Lab >(vInfo, coll, binfo, advection);
	//block_processing.process(vInfo, coll, penalization);
	//block_processing.process< I2D_VectorBlockLab< Streamer_Velocity, 2 >::Lab >(vInfo, coll, binfo, curl_diff);
	UpdateScalarRK2_Simple<2> stepB(dt);
	block_processing.process(vInfo, coll, stepB);
	
	//velsolver->compute_velocity();
}


