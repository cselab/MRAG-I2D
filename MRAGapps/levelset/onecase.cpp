/*
 * onecase.cpp — minimal SERIAL driver for the MRAG-I2D wavelet
 * multiresolution core (no TBB, no GLUT, no I2D vortex solver).
 *
 * What it demonstrates: the wavelet-detail refinement criterion that IS the
 * grid in ppm/c. We put a level-set (signed-distance) circle on a coarse
 * block grid, then let Science::AutomaticRefinement drive a fast wavelet
 * transform per block and refine wherever the detail coefficient exceeds a
 * threshold. The result is a grid whose fine cells track the interface.
 *
 * Build serially with -D_MRAG_SERIAL (see the Makefile here).
 */

#include <cmath>
#include <cstdio>
#include <vector>
#include <string>

#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGEnvironment.h"
#include "MRAGcore/MRAGBlock.h"
#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGRefiner.h"
#include "MRAGcore/MRAGCompressor.h"
#include "MRAGcore/MRAGBlockLab.h"
#include "MRAGcore/MRAGBlockFWT.h"
#include "MRAGcore/MRAGWavelets_AverageInterp5thOrder.h"
#include "MRAGscience/MRAGScienceCore.h"
#include "MRAGio/MRAG_IO_XDMF.h"

using namespace MRAG;

typedef Wavelets_AverageInterp5thOrder W;
typedef Block<Real, 8, 8, 1>        B;

// LeVeque single-vortex deformation benchmark initial shape:
// a circle of radius R0 centered at (XC,YC) in the unit square.
static const double R0 = 0.15;
static double XC = 0.5;
static double YC = 0.75;

// Smoothed-indicator "color" function: ~1 inside the circle, ~0 outside, with
// a sharp transition AT the interface. This is the quantity the LeVeque
// deformation benchmark advects, and it makes the wavelet detail coefficients
// localize on the zero contour (the field is flat away from the interface, so
// the refinement criterion only fires in the thin band around r=R0).
static inline double phi0(double x, double y)
{
	const double d   = std::sqrt((x - XC) * (x - XC) + (y - YC) * (y - YC)) - R0;
	const double eps = 0.02; // transition half-width (a few finest cells)
	return 0.5 * (1.0 - std::tanh(d / eps));
}

// fill the (possibly just-refined) grid with the initial level set
static void _ic(Grid<W, B>& grid)
{
	std::vector<BlockInfo> vInfo = grid.getBlocksInfo();
	for (int i = 0; i < (int)vInfo.size(); i++)
	{
		BlockInfo& info  = vInfo[i];
		B&         block = grid.getBlockCollection()[info.blockID];

		for (int iy = 0; iy < B::sizeY; iy++)
			for (int ix = 0; ix < B::sizeX; ix++)
			{
				double x[2];
				info.pos(x, ix, iy);
				block(ix, iy) = phi0(x[0], x[1]);
			}
	}
}

int main(int argc, char** argv)
{
	const int    nBlocks  = 4;
	const double tol      = 1e-3;
	const int    maxLevel = 5;
	const int    nsteps   = 20;

	Environment::setup();

	std::vector<IO_XDMF<W, B>::Frame> frames;
	for (int step = 0; step < nsteps; step++)
	{
		const double t = (double)step / nsteps;
		XC = 0.5 + 0.2 * std::cos(2.0 * M_PI * t);
		YC = 0.5 + 0.2 * std::sin(2.0 * M_PI * t);

		Grid<W, B>  grid(nBlocks, nBlocks, 1);
		Refiner     refiner;
		Compressor  compressor;
		grid.setRefiner(&refiner);
		grid.setCompressor(&compressor);
		BlockFWT<W, B> blockfwt;

		_ic(grid);
		Science::AutomaticRefinement<0, 0>(grid, blockfwt, tol, maxLevel, -1, NULL, _ic);

		char name[64];
		sprintf(name, "phi.%04d", step);
		const char* names[] = { "phi" };
		IO_XDMF<W, B> xdmf;
		frames.push_back(xdmf.Write(grid, name, t, step, names));
	}

	IO_XDMF<W, B>::WriteTemporalMaster("phi", frames);
	printf("done: %d frames; open phi.xdmf2 in ParaView\n", nsteps);
	return 0;
}
