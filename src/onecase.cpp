#include <cmath>
#include <cstdio>
#include <vector>
#include <string>

#include "MRAGCommon.h"
#include "MRAGEnvironment.h"
#include "MRAGBlock.h"
#include "MRAGrid.h"
#include "MRAGRefiner.h"
#include "MRAGCompressor.h"
#include "MRAGBlockLab.h"
#include "MRAGBlockFWT.h"
#include "MRAGWavelets_AverageInterp5thOrder.h"
#include "MRAGScienceCore.h"
#include "MRAG_IO_XDMF.h"

using namespace MRAG;

typedef Wavelets_AverageInterp5thOrder W;
typedef Block<Real, 8, 8, 1> B;
static const double R0 = 0.15;
static double XC = 0.5;
static double YC = 0.75;
static inline double phi0(double x, double y) {
  const double d = std::sqrt((x - XC) * (x - XC) + (y - YC) * (y - YC)) - R0;
  const double eps = 0.02; // transition half-width (a few finest cells)
  return 0.5 * (1.0 - std::tanh(d / eps));
}
static void _ic(Grid<W, B> &grid) {
  std::vector<BlockInfo> vInfo = grid.getBlocksInfo();
  for (int i = 0; i < (int)vInfo.size(); i++) {
    BlockInfo &info = vInfo[i];
    B &block = grid.getBlockCollection()[info.blockID];

    for (int iy = 0; iy < B::sizeY; iy++)
      for (int ix = 0; ix < B::sizeX; ix++) {
        double x[2];
        info.pos(x, ix, iy);
        block(ix, iy) = phi0(x[0], x[1]);
      }
  }
}
int main(int argc, char **argv) {
  const int nBlocks = 4;
  const double tol = 1e-3;
  const int maxLevel = 5;
  const int nsteps = 20;

  Environment::setup();

  Grid<W, B> grid(nBlocks, nBlocks, 1);
  Refiner refiner;
  Compressor compressor;
  grid.setRefiner(&refiner);
  grid.setCompressor(&compressor);
  BlockFWT<W, B> blockfwt;

  std::vector<IO_XDMF<W, B>::Frame> frames;
  for (int step = 0; step < nsteps; step++) {
    const double t = (double)step / nsteps;
    XC = 0.5 + 0.2 * std::cos(2.0 * M_PI * t);
    YC = 0.5 + 0.2 * std::sin(2.0 * M_PI * t);

    _ic(grid);
    Science::AutomaticRefinement<0, 0>(grid, blockfwt, tol, maxLevel, -1, _ic);
    Science::AutomaticCompression<0, 0>(grid, blockfwt, tol);
    printf("  step %d: grid has %zu blocks after refine+compress\n", step,
           grid.getBlocksInfo().size());

    char name[64];
    sprintf(name, "phi.%04d", step);
    const char *names[] = {"phi"};
    IO_XDMF<W, B> xdmf;
    frames.push_back(xdmf.Write(grid, name, t, step, names));
  }
  IO_XDMF<W, B>::WriteTemporalMaster("phi", frames);
  printf("done: %d frames; open phi.xdmf2 in ParaView\n", nsteps);
  return 0;
}
