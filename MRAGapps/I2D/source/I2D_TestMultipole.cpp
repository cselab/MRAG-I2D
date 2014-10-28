//
//  I2D_TestMultipole.cpp
//  I2D_ROCKS
//
//  Created by Wim van Rees on 4/1/13.
//
//

#include "I2D_TestMultipole.h"
#include "I2D_VelocitySolver_Mani.h"
#include "I2D_VelocitySolver_Wim.h"

#include <limits>


I2D_TestMultipole::I2D_TestMultipole(const int argc, const char ** argv): parser(argc, argv),t(0), step_id(0)
{
	
	printf("////////////////////////////////////////////////////////////\n");
	printf("////////////       MULTIPOLE TEST            ///////////////\n");
	printf("////////////////////////////////////////////////////////////\n");
    
    bRestart = parser("-restart").asBool();
    
	const int bpd = parser("-bpd").asInt();
	//assert(bpd > 1);
	
	const string fmm_name = parser("-fmm").asString();
	
	grid = new Grid<W,B>(bpd,bpd,1);
	assert(grid != NULL);
	
	const int res_jump = max(1, parser("-jump").asInt());
	const int lmax = max(1, parser("-lmax").asInt());
	
	refiner = new Refiner_BlackList(res_jump, lmax);
	compressor = new Compressor(res_jump);
	grid->setRefiner(refiner);
	grid->setCompressor(compressor);
    
	if(fmm_name == "velocity")
		poisson_solver = new I2D_VelocitySolver_Mani(*grid, parser);
	else if(fmm_name == "velocity-wim")
		poisson_solver = new I2D_VelocitySolver_Wim(*grid, parser);
	else
	{
        printf("No valid poisson solver specified! Aborting\n");
        abort();
	}
    
    LambOseenGamma = 10.0;
    LambOseenNut = 1e-3;
    LambOseenOrg[0] = 0.5;
    LambOseenOrg[1] = 0.5;
    
    if(bRestart)
    {
        _restart();
        _storeExactSolution();
        _dump("multipole_rst");
    }
    else
    {
        _ic_omega(*grid);
        _refine(true);
        _compress(true);
        _dump("multipole_ic");
    }
}


set<int> I2D_TestMultipole::_getBoundaryBlockIDs()
{
	set<int> result;
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	for(vector<BlockInfo>::iterator it = vInfo.begin(); it != vInfo.end(); it++)
	{
		const bool bX = it->index[0] == 0 || it->index[0] == pow(2, it->level)-1;
		const bool bY = it->index[1] == 0 || it->index[1] == pow(2, it->level)-1;
		
		if (bX || bY ) result.insert(it->blockID);
	}
	
	return result;
}

void I2D_TestMultipole::_refine(bool bUseIC)
{
	if (!parser("-uniform").asBool())
	{
		while(true)
		{
            //			set<int> boundary_blocks = _getBoundaryBlockIDs();
            //			((Refiner_BlackList*)refiner)->set_blacklist(&boundary_blocks);
            //			const int refinements = Science::AutomaticRefinement<0,2>(*grid, fwt_omega, parser("-rtol").asDouble(), parser("-lmax").asInt(), 1, NULL, initial_condition, &boundary_blocks);
			const int refinements = Science::AutomaticRefinement<0,0>(*grid, fwt_omega, parser("-rtol").asDouble(), parser("-lmax").asInt(), 1, NULL, (void (*)(Grid<W,B>&))NULL);
			
            if(bUseIC)
                _ic_omega(*grid);
            
			if (refinements == 0) break;
		}
	}
}

void I2D_TestMultipole::_compress(bool bUseIC)
{
	if (!parser("-uniform").asBool())
		Science::AutomaticCompression<0,0>(*grid, fwt_omega, parser("-ctol").asDouble(), 1, NULL, (void (*)(Grid<W,B>&))NULL);
    
    if(bUseIC)
        _ic_omega(*grid);
}

void I2D_TestMultipole::paint()
{
}


void I2D_TestMultipole::_ic_omega(Grid<W,B>& grid)
{
    LambOseenVortex lovortex(LambOseenGamma, LambOseenOrg);
    
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		B& b = grid.getBlockCollection()[info.blockID];
		
        for(int iy=0; iy<B::sizeY; iy++)
            for(int ix=0; ix<B::sizeX; ix++)
            {
                Real p[2];
                info.pos(p, ix, iy);
                
                Real w=0;
                lovortex.gimmeVort(p,LambOseenNut,w);
                
                // initialize field
                b(ix,iy).omega = w;
                b(ix,iy).u[0]=b(ix,iy).u[1]=0;
                b(ix,iy).tmp=0;
            }
	}
}

void I2D_TestMultipole::_computeError()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	LambOseenVortex lovortex(LambOseenGamma, LambOseenOrg);
    
    Real L1=0.;
    Real L2=0.;
    Real LI=0.;
    Real L1bnd=0.;
    Real L2bnd=0.;
    Real LIbnd=0.;

    Real L1_ex=0.;
    Real L2_ex=0.;
    Real LI_ex=0.;
    
    assert(exactsol.size() == vInfo.size());
    
	const bool isreference = parser("-fmm").asString() == "potential";
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		B& b = grid->getBlockCollection()[info.blockID];
		const Real dA = info.h[0]*info.h[1];
		
        for(int iy=0; iy<B::sizeY; iy++)
            for(int ix=0; ix<B::sizeX; ix++)
            {
                Real p[2];
                info.pos(p, ix, iy);

                const Real err = std::sqrt( std::pow(b(ix,iy).u[0] - exactsol[i].u_exact[iy][ix][0],2) + std::pow(b(ix,iy).u[1] - exactsol[i].u_exact[iy][ix][1], 2) );
                
                const Real errBnd = b(ix,iy).tmp;

                L1 += err*dA;
                L2 += err*err*dA;
                LI = std::max(LI,err);
                
                L1bnd += errBnd*dA;
                L2bnd += errBnd*errBnd*dA;
                LIbnd = std::max(errBnd,LIbnd);
                
                Real exactVel[2];
                lovortex.gimmeVel(p,LambOseenNut,exactVel);

                const Real err_ex = std::sqrt( std::pow(b(ix,iy).u[0] - exactVel[0],2) + std::pow(b(ix,iy).u[1] - exactVel[1], 2) );
                
                L1_ex += err_ex*dA;
                L2_ex += err_ex*err_ex*dA;
                LI_ex = std::max(LI_ex,err_ex);
            }
	}
    L2 = std::sqrt(L2);
    L2bnd = std::sqrt(L2bnd);
    L2_ex = std::sqrt(L2_ex);
    
    const int bpd = parser("-bpd").asInt();
    const int lmax =  grid->getCurrentMaxLevel();
    int res;
    if (parser("-uniform").asBool())
        res = bpd*FluidBlock2D::sizeX;
    else
        res = FluidBlock2D::sizeX*pow(2,lmax);
    const double nofgridpoints = grid->getBlocksInfo().size()*(double)(B::sizeX*B::sizeY*B::sizeZ);
    
    printf("========= VELOCITY ERRORS AT RES %d, ORDER %d, THETA %.1f ========\n", res,_ORDER_,parser("-fmm-theta").asDouble());
    printf("L1: %e \t %e \t %e\n",L1,L1bnd,L1_ex);
    printf("L2: %e \t %e \t %e\n",L2,L2bnd,L2_ex);
    printf("LI: %e \t %e \t %e\n",LI,LIbnd,LI_ex);
    printf("========= END VELOCITY ===========\n");
}

void I2D_TestMultipole::run()
{
    const string fmm_name = parser("-fmm").asString();
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    _dump("before_potential");
    
    profiler.push_start("POTENTIAL");
    poisson_solver->compute_velocity();
    profiler.pop_stop();;
    
    profiler.push_start("ERROR");
	if(bRestart) _computeError();
    profiler.pop_stop();
    
    if (!bRestart) _save();
    
    _dump("after_potential");
    
    {
		profiler.printSummary();
        
		FILE * f = fopen("perf.txt", "w");
		assert(f!=NULL);
		profiler.printSummaryToFile(f);
		fclose(f);
	}
    
    std::cout << "DONE: ABORTING" << std::endl;
    abort();
}


void I2D_TestMultipole::_dump(string filename)
{
	IO_VTKNative<W,B, 4,0> vtkdumper;
	vtkdumper.Write(*grid, grid->getBoundaryInfo(), filename);
}

void I2D_TestMultipole::_restart()
{
	//read status
	{
		FILE * f = fopen("restart.status", "r");
		assert(f != NULL);
		float val = -1;
		fscanf(f, "time: %e\n", &val);
		assert(val>=0);
		t=val;
		step_id = -1;
		fscanf(f, "stepid: %d\n", &step_id);
		assert(step_id >= 0);
		fclose(f);
	}
	
	printf("DESERIALIZATION: time is %f and step id is %d\n", t, step_id);
	
	//read grid
	IO_Binary<W,B> serializer;
	serializer.Read(*grid, "restart");
}

void I2D_TestMultipole::_save()
{
	printf("****SERIALIZING****\n");
	//write status
	{
		FILE * f = fopen("restart.status", "w");
		assert(f != NULL);
		fprintf(f, "time: %20.20e\n", t);
		fprintf(f, "stepid: %d\n", step_id);
		fclose(f);
		
		printf( "time: %20.20e\n", t);
		printf( "stepid: %d\n", step_id);
	}
	
	//write grid
    IO_Binary<W,B> serializer;
	serializer.Write(*grid, "restart");
	printf("****SERIALIZING DONE****\n");
    
}

void I2D_TestMultipole::_storeExactSolution()
{
    std::vector<BlockInfo> vInfo = grid->getBlocksInfo();
    exactsol.clear();
    exactsol.resize(vInfo.size());
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		B& b = grid->getBlockCollection()[info.blockID];
		
        for(int iy=0; iy<B::sizeY; iy++)
            for(int ix=0; ix<B::sizeX; ix++)
            {
                exactsol[i].u_exact[iy][ix][0] = b(ix,iy).u[0];
                exactsol[i].u_exact[iy][ix][1] = b(ix,iy).u[1];
            }
    }
}
