#include <cmath>
#include <cstdio>
#include <vector>
#include <string>

#include "MRAGCommon.h"
#include "MRAGEnvironment.h"
#include "MRAGBlock.h"
#include "MRAGrid.h"
#include "MRAGRefiner_SpaceExtension.h"
#include "MRAGCompressor.h"
#include "MRAGBlockLab.h"
#include "MRAGBlockFWT.h"
#include "MRAGWavelets_Interp4thOrder.h"
#include "MRAGScienceCore.h"
#include "MRAGBlockProcessing_TBB.h"
#include "MRAG_IO_XDMF.h"
#include "MRAG_IO_Blocks.h"

using namespace MRAG;

struct PS {
  Real w, tmp, u[2];
  PS() : w(0), tmp(0) { u[0] = u[1] = 0; }
  void operator+=(PS t) { w += t.w; tmp += t.tmp; u[0] += t.u[0]; u[1] += t.u[1]; }
  operator Real() { return w; }
  Real giveMe(int, Real) const { return w; }
};
PS operator*(const PS &p, Real v) {
  PS t;
  t.w = p.w * v; t.tmp = p.tmp * v; t.u[0] = p.u[0] * v; t.u[1] = p.u[1] * v;
  return t;
}
template <typename T, int i> inline Real scalar_projector_impl(const T &t) {
  return (Real)t.w;
}
make_projector(WARP_projector, scalar_projector_impl)

static const int blockSize = 8;
static const int nBlocks = 4;
static const int maxLevel = 6;
static const int resJump = 1;
static const int maxStencil[2][3] = {-3, -3, 0, +4, +4, 1};
static const double ctol = 5e-5;
static const double rtol = ctol / 10;

typedef Wavelets_Interp4thOrder W;
struct WARPBlock : public Block<PS, blockSize, blockSize, 1> {
  typedef PS ElementType;
  Real maximumDT;
  WARPBlock(ElementType e = ElementType())
      : Block<PS, blockSize, blockSize, 1>(e) {}
};
typedef WARPBlock B;
typedef Multithreading::BlockProcessing_TBB<B> BlockProcessing;

#include "blockprocessingWARP.h"

static void _ic(Grid<W, B> &grid) {
  std::vector<BlockInfo> vInfo = grid.getBlocksInfo();
  for (int i = 0; i < (int)vInfo.size(); i++) {
    BlockInfo &info = vInfo[i];
    B &block = grid.getBlockCollection().lock(info.blockID);
    for (int iy = 0; iy < B::sizeY; iy++)
      for (int ix = 0; ix < B::sizeX; ix++) {
        float x[3];
        info.pos(x, ix, iy, 0);
        PS &p = block(ix, iy, 0);
        const double r1 = std::sqrt(0.25 * (x[0] - 0.5) * (x[0] - 0.5) +
                                    (x[1] - 0.5) * (x[1] - 0.5));
        const double eps = 0.1;
        const double alpha =
            M_PI * std::min(1.0, std::max(0.0, (r1 - 0.15 + 0.5 * eps) / eps));
        p.w = 0.5 + 0.5 * std::cos(alpha);
        const double px = M_PI * x[0], py = M_PI * x[1];
        p.u[0] = -2.0 * std::sin(px) * std::sin(px) * std::sin(py) * std::cos(py);
        p.u[1] = 2.0 * std::sin(py) * std::sin(py) * std::sin(px) * std::cos(px);
      }
    grid.getBlockCollection().release(info.blockID);
  }
}

int main(int argc, char **argv) {
  const int nsteps = 20;
  const double frameDt = 0.05;

  Environment::setup(argc > 1 ? atoi(argv[1]) : -1);

  Grid<W, B> grid(nBlocks, nBlocks, 1, maxStencil);
  Refiner_SpaceExtension refiner(resJump);
  Compressor compressor(resJump);
  grid.setRefiner(&refiner);
  grid.setCompressor(&compressor);
  BlockFWT<W, B, WARP_projector> blockfwt;
  BlockProcessing blockProcessing;

  _ic(grid);
  Science::AutomaticRefinement<0, 0>(grid, blockfwt, rtol, maxLevel, -1, _ic);
  Science::AutomaticCompression<0, 0>(grid, blockfwt, ctol, 1, _ic);

  double t = 0;
  std::vector<IO_XDMF<W, B>::Frame> frames;
  std::vector<IO_Blocks<W, B>::Frame> bframes;
  for (int step = 0; step < nsteps; step++) {
    Science::AutomaticRefinement<0, 0>(grid, blockfwt, rtol, maxLevel, -1);

    std::vector<BlockInfo> vInfo = grid.getBlocksInfo();
    BoundaryInfo &bInfo = grid.getBoundaryInfo();
    const BlockCollection<B> &coll = grid.getBlockCollection();

    double tnext = t + frameDt;
    while (t < tnext) {
      CLFParticle<blockSize>::FindMaximumDT findmaximumDT;
      BlockProcessing::process(vInfo, coll, findmaximumDT);
      Real dt = HUGE_VALF;
      for (int i = 0; i < (int)vInfo.size(); i++) {
        B &b = coll.lock(vInfo[i].blockID);
        dt = min(dt, b.maximumDT);
        coll.release(vInfo[i].blockID);
      }
      dt = min(dt, (Real)(tnext - t));

      CLFParticle<blockSize>::AdvectParticlesRK3 advect(dt);
      blockProcessing.process<CLFParticle<blockSize>::BlockLabWARP>(
          vInfo, coll, bInfo, advect);
      CLFParticle<blockSize>::UpdateOmega updateomega;
      BlockProcessing::process(vInfo, coll, updateomega);
      t += dt;
    }

    Science::AutomaticCompression<0, 0>(grid, blockfwt, ctol, 1);

    char name[64];
    sprintf(name, "omega.%04d", step);
    const char *names[] = {"omega"};
    IO_XDMF<W, B> xdmf;
    frames.push_back(xdmf.Write(grid, name, t, step, names));

    char bname[64];
    sprintf(bname, "blocks.%04d", step);
    bframes.push_back(IO_Blocks<W, B>::Write(grid, bname, t));

    printf("step %d t=%.4f blocks=%zu\n", step, t,
           grid.getBlocksInfo().size());
  }
  IO_XDMF<W, B>::WriteTemporalMaster("omega", frames);
  IO_Blocks<W, B>::WriteTemporalMaster("blocks", bframes);
  return 0;
}
