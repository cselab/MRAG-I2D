/*
 *  I2D_PenalizationOperator.cpp
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include "I2D_PenalizationOperator.h"
#include "I2D_VectorBlockLab.h"
#include "I2D_GradOfVector.h"

struct Penalization
{
	Real dt, lambda;
	Real Uinf[2];
	
	Penalization(Real dt, Real lambda, const Real Uinf[2]): lambda(lambda), dt(dt)
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
	}
	
	Penalization(const Penalization& c): lambda(c.lambda), dt(c.dt)
	{
		Uinf[0] = c.Uinf[0];
		Uinf[1] = c.Uinf[1];
	}
	
	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{		
		FluidElement2D * const e = &b(0,0);
		Real * const ext = &b.external_data[0][0];
		
		const Real lamdt = lambda*dt;
		
		static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
		for(int i=0; i<n; i++)			
		{
			const double lamdtX  = lamdt*e[i].tmp;
			
			const double u_current[2] = {
				e[i].u[0],
				e[i].u[1],
			};
			
			e[i].tmp = u_current[0];
			ext[i] = u_current[1];
			
			e[i].u[0] = (u_current[0] - Uinf[0]*lamdtX)/(1 + lamdtX) - u_current[0];
			e[i].u[1] = (u_current[1] - Uinf[1]*lamdtX)/(1 + lamdtX) - u_current[1];
		}
	}
};

struct PenalizationCustomizedVelocity
{
	Real dt, lambda;
	Real Uinf[2];
	const map<int, const VelocityBlock *>& customized_velocity;
	
	PenalizationCustomizedVelocity(Real dt, Real lambda, const map<int, const VelocityBlock *>& customized_velocity, const Real Uinf[2]):  lambda(lambda), dt(dt), customized_velocity(customized_velocity)
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
	}
	
	PenalizationCustomizedVelocity(const PenalizationCustomizedVelocity& c): lambda(c.lambda), dt(c.dt), customized_velocity(c.customized_velocity)
	{
		this->Uinf[0] = c.Uinf[0];
		this->Uinf[1] = c.Uinf[1];
	}
	
	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{		
		FluidElement2D * const e = &b(0,0);
		Real * const ext = &b.external_data[0][0];
		
		const Real lamdt = lambda*dt;
						
		const map<int, const VelocityBlock *>::const_iterator it = customized_velocity.find(info.blockID);

		if (it!=customized_velocity.end())
		{
			const VelocityBlock * u_desired = it->second;
			assert(u_desired != NULL);

			for(int iy=0; iy<FluidBlock2D::sizeY; iy++)		
				for(int ix=0; ix<FluidBlock2D::sizeX; ix++)		
				{
					const int i = ix + FluidBlock2D::sizeX*iy;
					
					const Real lamdtX  = lamdt*e[i].tmp;
					
					const Real u_current[2] = {
						e[i].u[0],
						e[i].u[1],
					};
					
					e[i].tmp = u_current[0];
					ext[i] = u_current[1];
					
					e[i].u[0] = (u_current[0] + (u_desired->u[0][iy][ix]-Uinf[0]) * lamdtX)/(1 + lamdtX) - u_current[0];
					e[i].u[1] = (u_current[1] + (u_desired->u[1][iy][ix]-Uinf[1]) * lamdtX)/(1 + lamdtX) - u_current[1];
				}
		}
		else 
		{
			const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;

			for(int i=0; i<n; i++)			
			{
				assert(e[i].tmp == 0);
				
				e[i].tmp = e[i].u[0];
				ext[i] = e[i].u[1];
				
				e[i].u[0] = 0;
				e[i].u[1] = 0;
			}
		}
		
	}
};

struct CurlDiff_4thOrder: I2D_GradOfVector_4thOrder
{
	Real t;
	int stencil_start[3], stencil_end[3];
	
	CurlDiff_4thOrder(): t(0)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	CurlDiff_4thOrder(const CurlDiff_4thOrder & copy): t(0)
	{
		stencil_start[0] = stencil_start[1] = -2;
		stencil_end[0] = stencil_end[1] = +3;
		stencil_start[2] = 0;
		stencil_end[2] = 1;
	}
	
	struct AddToOmega
	{ static inline void stream(FluidElement2D& out, Real in) { out.omega += in; } };
	
	struct SubToOmega
	{ static inline void stream(FluidElement2D& out, Real in) { out.omega -= in; } };
	
	template<typename Lab>
	inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock2D& out) const 
	{
		_dfdx<AddToOmega, 1 >(lab, info, out);
		_dfdy<SubToOmega, 0 >(lab, info, out);
	}
};

struct UpdateVelocities
{
	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{		
		FluidElement2D * const e = &b(0,0);
		const Real * const ext = &b.external_data[0][0];
		
		static const int n = FluidBlock2D::sizeY*FluidBlock2D::sizeX;
		for(int i=0; i<n; i++)			
		{
			e[i].u[0] += e[i].tmp;
			e[i].u[1] += ext[i];
		}
	}
};

void I2D_PenalizationOperator::perform_timestep(Real dt)
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	BoundaryInfo& binfo=grid.getBoundaryInfo();  
	const BlockCollection<B>& coll = grid.getBlockCollection();
	
	if (desired_velocities == NULL)
	{
		Penalization penalization(dt, lambda, Uinf);
		block_processing.process(vInfo, coll, penalization);
	}
	else 
	{
		PenalizationCustomizedVelocity penalization(dt, lambda, *desired_velocities, Uinf);
		block_processing.process(vInfo, coll, penalization);
		
		desired_velocities = NULL;
	}

	CurlDiff_4thOrder curl_diff;
	block_processing.process< I2D_VectorBlockLab< Streamer_Velocity, 2 >::Lab >(vInfo, coll, binfo, curl_diff);
	
	UpdateVelocities update;
	block_processing.process(vInfo, coll, update);
}


struct Diagnostic
{
	Real force[2];
	Real torque[2];
	Real area;
	Real gridpoints;
	Real circulation;	
	
	Diagnostic():area(0), gridpoints(0), circulation(0)
	{
		force[0] = force[1] = 0;
		torque[0] = torque[1] = 0;
	}
	
	const Diagnostic& operator += (const Diagnostic& c)
	{
		force[0] += c.force[0];
		force[1] += c.force[1];
		torque[0] += c.torque[0];
		torque[1] += c.torque[1];
		area += c.area;
		gridpoints += c.gridpoints;
		circulation += c.circulation;
		
		return *this;
	}
};

struct ComputeDiagnostics
{
	Real lambda, Uinf[2], cor[2]; // cor is Center Of Rotation
	map<int, Diagnostic>& diagnostics;
	
	ComputeDiagnostics(map< int, Diagnostic>& diagnostics, Real lambda, const Real Uinf[2], const Real cor[2]):diagnostics(diagnostics), lambda(lambda)
	{
		this->Uinf[0] = Uinf[0];
		this->Uinf[1] = Uinf[1];
		this->cor[0] = cor[0];
		this->cor[1] = cor[1];
	}
	
	ComputeDiagnostics(const ComputeDiagnostics& c): diagnostics(c.diagnostics), lambda(c.lambda)
	{
		Uinf[0] = c.Uinf[0];
		Uinf[1] = c.Uinf[1];
		cor[0] = c.cor[0];
		cor[1] = c.cor[1];
	}
	
	inline void operator()(const BlockInfo& info, FluidBlock2D& b) const
	{
		map< int, Diagnostic>::iterator it = diagnostics.find(info.blockID);
		assert(it != diagnostics.end());
		Diagnostic& diag = it->second;
		
		diag.force[0] = diag.force[1]  = 0.0;
		diag.torque[0] = diag.torque[1]  = 0.0;
		diag.area = 0.0;
		diag.gridpoints = 0.0;
		diag.circulation = 0.0;
		
		for(int iy=0; iy<FluidBlock2D::sizeY; iy++)
			for(int ix=0; ix<FluidBlock2D::sizeX; ix++)
			{
				Real p[2];
				info.pos(p,ix,iy);
				
				const Real Xs = b(ix,iy).tmp;
				diag.force[0] += b(ix,iy).u[0]*Xs;
				diag.force[1] += b(ix,iy).u[1]*Xs;
				
				diag.torque[0] +=  (p[0]-cor[0]) * (b(ix,iy).u[1]+Uinf[1]) * Xs;
				diag.torque[1] += -(p[1]-cor[1]) * (b(ix,iy).u[0]+Uinf[0]) * Xs;
				
				diag.area += Xs;
				diag.gridpoints += Xs;				
				diag.circulation += b(ix,iy).omega;
			}
		
		const Real dA = pow(info.h[0], 2);
		
		diag.force[0] *= dA*lambda;
		diag.force[1] *= dA*lambda;
		
		diag.torque[0] *= dA*lambda;
		diag.torque[1] *= dA*lambda;
		
		diag.area *= dA;		
		diag.circulation *= dA;
		
		diag.force[0] += diag.area*lambda*Uinf[0];
		diag.force[1] += diag.area*lambda*Uinf[1];
	}
};

void I2D_PenalizationOperator::compute_dragandstuff(Real time, const Real D, const Real cor[2], string filename)
{
	const Real maxu = max(fabs(Uinf[0]), fabs(Uinf[1]));
	const Real U_infinity = (maxu==0.0)?1:maxu;
	
	map<int, Diagnostic> diagnostics;
	
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	for(vector<BlockInfo>::iterator it=vInfo.begin(); it!=vInfo.end(); it++)
		diagnostics[it->blockID] = Diagnostic();
	
	ComputeDiagnostics get_diag(diagnostics, lambda, Uinf, cor);
	block_processing.process(vInfo, grid.getBlockCollection(), get_diag);
	
	Diagnostic global;
	for(map< int, Diagnostic>::iterator it=diagnostics.begin(); it!=diagnostics.end(); it++)
		global += it->second;
	
	const Real cD = 2*global.force[0]/(pow(U_infinity, 2)*D);
	const Real cL = 2*global.force[1]/(pow(U_infinity, 2)*D);
	
	this->Cd = cD;
	this->Cl = cL;
	
	const Real rho = 1.0;
	const Real q = 0.5*rho*U_infinity*U_infinity;
	const Real Cm = -global.torque[0]/(q*D*D);
	
	const Real T = 2*U_infinity*time/D;
	
	const Real Area = global.area;
	const Real GridPoints = global.gridpoints;
	const Real Circulation = fabs(global.circulation);
	
	const int totalNumPoints = vInfo.size()*B::sizeX*B::sizeX;
	
	FILE * f = fopen(filename.data(), bAppendToFile ? "a" : "w");
	assert(f!=NULL);
	if (!bAppendToFile)
		fprintf(f, "T\t\tcD\t\tcL\t\tCm\t\tcirculation\t\tArea\t\tpoints\t\tN-obstcl-1D\n");
	
	fprintf(f, "%e\t%e\t%e\t%e\t%e\t%e\t%d\t%f\n", T, cD, cL, Cm, Circulation, Area, totalNumPoints, pow((float)GridPoints, (float)(1./2.)));
	
	fclose(f);
	
	printf("%e\t%e\t%e\t%e\t%e\t%e\t%d\t%f\n", T, cD, cL, Cm, Circulation, Area, totalNumPoints, pow((float)GridPoints, (float)(1./2.)));
	
	bAppendToFile =  true;
}

