#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>

#include "MRAGCommon.h"
#include "MRAGEnvironment.h"
#include "MRAGrid.h"
#include "MRAGRefiner.h"
#include "MRAGCompressor.h"
#include "MRAGBlockFWT.h"
#include "MRAGWavelets_AverageInterp5thOrder.h"
#include "MRAGSimpleLevelsetBlock.h"
#include "MRAGScienceCore.h"
#include "MRAG_IO_XDMF.h"
#include "MRAG_IO_Blocks.h"

using namespace MRAG;

struct LevelsetPoint {
  Real phi;
  LevelsetPoint(Real v = 0) : phi(v) {}
  void operator+=(LevelsetPoint t) { phi += t.phi; }
  operator Real() { return phi; }
  Real levelset() const { return phi; }
  Real giveMe(int, Real) const { return phi; }
};
LevelsetPoint operator*(const LevelsetPoint &p, Real v) {
  return LevelsetPoint(p.phi * v);
}
template <typename T, int i> inline Real ls_projector_impl(const LevelsetPoint &t) {
  return (Real)t.levelset();
}
make_projector(ls_projector, ls_projector_impl)

typedef Wavelets_AverageInterp5thOrder W;
typedef SimpleLevelsetBlock<LevelsetPoint, 4, 8, 8, 1> B;

static const double R0 = 0.15;
static double XC = 0.5;
static double YC = 0.75;
static inline double phi0(double x, double y) {
  return std::sqrt((x - XC) * (x - XC) + (y - YC) * (y - YC)) - R0;
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
        block(ix, iy) = LevelsetPoint(phi0(x[0], x[1]));
      }
  }
}
int main(int argc, char **argv) {
  const int nBlocks = 4;
  const double tol = 1e-7;
  const int maxLevel = 5;
  const int nsteps = 20;

  Environment::setup(argc > 1 ? atoi(argv[1]) : -1);

  Grid<W, B> grid(nBlocks, nBlocks, 1);
  Refiner refiner;
  Compressor compressor;
  grid.setRefiner(&refiner);
  grid.setCompressor(&compressor);
  BlockFWT<W, B, ls_projector> blockfwt;

  std::vector<IO_XDMF<W, B>::Frame> frames;
  std::vector<IO_Blocks<W, B>::Frame> bframes;
  for (int step = 0; step < nsteps; step++) {
    const double t = (double)step / nsteps;
    XC = 0.5 + 0.2 * std::cos(2.0 * M_PI * t);
    YC = 0.5 + 0.2 * std::sin(2.0 * M_PI * t);

    _ic(grid);
    Science::AutomaticRefinementForLevelsets(grid, blockfwt, tol, maxLevel, true,
                                             -1, NULL, _ic);
    Science::AutomaticCompressionForLevelsets(grid, blockfwt, tol, true, 1);
    printf("  step %d: grid has %zu blocks\n", step,
           grid.getBlocksInfo().size());

    char name[64];
    sprintf(name, "phi.%04d", step);
    const char *names[] = {"phi"};
    IO_XDMF<W, B> xdmf;
    frames.push_back(xdmf.Write(grid, name, t, step, names));

    char bname[64];
    sprintf(bname, "blocks.%04d", step);
    bframes.push_back(IO_Blocks<W, B>::Write(grid, bname, t));
  }
  IO_XDMF<W, B>::WriteTemporalMaster("phi", frames);
  IO_Blocks<W, B>::WriteTemporalMaster("blocks", bframes);
  return 0;
}
