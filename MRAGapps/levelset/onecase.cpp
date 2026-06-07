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
#include <map>

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

using namespace MRAG;

typedef Wavelets_AverageInterp5thOrder W;
typedef Block<Real, 32, 32, 1>        B;

// LeVeque single-vortex deformation benchmark initial shape:
// a circle of radius R0 centered at (XC,YC) in the unit square.
static const double R0 = 0.15;
static const double XC = 0.5;
static const double YC = 0.75;

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

// report how many blocks live at each refinement level
static void _report(Grid<W, B>& grid, const char* tag)
{
	std::vector<BlockInfo> vInfo = grid.getBlocksInfo();
	std::map<int, int> perLevel;
	for (int i = 0; i < (int)vInfo.size(); i++)
		perLevel[vInfo[i].level]++;

	printf("[%s] %zu blocks total (%zu cells)\n",
	       tag, vInfo.size(), vInfo.size() * (size_t)(B::sizeX * B::sizeY));
	for (std::map<int, int>::const_iterator it = perLevel.begin(); it != perLevel.end(); ++it)
		printf("    level %d : %d blocks\n", it->first, it->second);
}

// dump the leaf cells as an ASCII point cloud (x y phi level) for plotting
static void _dump(Grid<W, B>& grid, const char* fname)
{
	FILE* f = fopen(fname, "w");
	if (!f) { printf("could not open %s\n", fname); return; }

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
				fprintf(f, "%g\t%g\t%g\t%d\n", x[0], x[1], (double)block(ix, iy), info.level);
			}
	}
	fclose(f);
	printf("wrote %s\n", fname);
}

int main(int argc, char** argv)
{
	const int    nBlocks  = 4;      // 4x4 coarse blocks
	const double tol      = 1e-3;   // wavelet-detail refinement threshold
	const int    maxLevel = 5;      // cap on refinement depth

	// no-op under -D_MRAG_SERIAL; sets the TBB worker count otherwise
	Environment::setup();

	Grid<W, B>  grid(nBlocks, nBlocks, 1);
	Refiner     refiner;
	Compressor  compressor;
	grid.setRefiner(&refiner);
	grid.setCompressor(&compressor);

	BlockFWT<W, B> blockfwt;

	_ic(grid);
	_report(grid, "initial");

	// wavelet-detail-driven adaptive refinement around the interface.
	// AutomaticRefinement runs a per-block fast wavelet transform and refines
	// any block whose detail coefficient exceeds `tol`, looping until the whole
	// grid is below tolerance (or maxLevel is reached). The grid ends up just
	// fine enough to resolve the interface band to tolerance -- exactly the
	// multiresolution behaviour behind ppm/c.
	Science::AutomaticRefinement<0, 0>(grid, blockfwt, tol, maxLevel, -1, NULL, _ic);
	_report(grid, "refined");
	_dump(grid, "phi_refined.txt");

	// NOTE: Science::AutomaticCompression is intentionally not called here.
	// It does not terminate for this case: grid.compress() refuses to collapse
	// blocks that would violate the 2:1 grading constraint and returns early
	// with nCollapses==0, yet the loop sees a non-zero phantom count and never
	// breaks. The refinement above is the demonstration of the wavelet
	// multiresolution criterion; round-trip compression would need the
	// levelset-specific AutomaticCompressionForLevelsets path.

	printf("done.\n");
	return 0;
}
