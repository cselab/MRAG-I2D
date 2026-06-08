/*
 *  mainWARP.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/27/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include <stdlib.h>
#include <iostream>
#include "../MRAG/MRAGcore/MRAGCommon.h"
#include "../MRAG/MRAGcore/MRAGEnvironment.h"

#if _MRAG_OS == _MRAG_OS_APPLE
#if defined(_MRAG_GLUT_VIZ) && defined(__APPLE__)
#include "GLUT/glut.h"
#elif defined(_MRAG_GLUT_VIZ)
#include <stdlib.h>
#include "GL/glew.h"
#include "GL/glut.h"
#endif
#endif

#undef min
#undef max

#include "../MRAG/MRAGcore/MRAGWavelets_AverageInterp5thOrder.h"
#include "../MRAG/MRAGcore/MRAGWavelets_Interp4thOrder.h"
#include "../MRAG/MRAGcore/MRAGWavelets_AverageInterp3rdOrder.h"
#include "../MRAG/MRAGcore/MRAGWavelets_Interp2ndOrder.h"
#include "../MRAG/MRAGcore/MRAGWavelets_Haar.h"
#include "../MRAG/MRAGcore/MRAGBlock.h"
#include "../MRAG/MRAGcore/MRAGrid.h"
#include "../MRAG/MRAGcore/MRAGRefiner.h"
#include "../MRAG/MRAGcore/MRAGCompressor.h"
#include "../MRAG/MRAGcore/MRAGBlockLab.h"
#include "../MRAG/MRAGcore/MRAGBlockFWT.h"

#ifdef _MRAG_GLUT_VIZ
#include "../MRAG/MRAGvisual/GridViewer.h"
#endif

#include "../MRAG/MRAGscience/MRAGScienceCore.h"
#include "../MRAG/MRAGscience/MRAGAutomaticRefiner.h"
#include "../MRAG/MRAGscience/MRAGSimpleLevelsetBlock.h"
#include "../MRAG/MRAGscience/MRAGSpaceTimeSorter.h"
#include "../MRAG/MRAGscience/MRAGRefiner_SpaceExtension.h"

#include "../MRAG/MRAGmultithreading/MRAGBlockProcessing_SingleCPU.h"
#include "../MRAG/MRAGmultithreading/MRAGBlockProcessing_TBB.h"

using namespace MRAG;

struct PS {
  float w, tmp;
  float u[2];

  PS() : w(0), tmp(0) { u[0] = u[1] = 0; }

  void operator+=(PS t) {
    w += t.w;
    tmp += t.tmp;
    u[0] += t.u[0];
    u[1] += t.u[1];
  }

  operator Real() { return (Real)w; }
};

PS operator*(const PS &p, Real v) {
  PS t;

  t.w = p.w * v;
  t.tmp = p.tmp * v;
  t.u[0] = p.u[0] * v;
  t.u[1] = p.u[1] * v;

  return t;
}

template <typename T, int i> inline Real scalar_projector_impl(const T &t) {
  return (Real)t.w;
}

make_projector(WARP_projector, scalar_projector_impl)

    static const int blockSize = 16;
static const int blocksPerDimension = 4;
static const int maxLevel = 12;
static const int resJump = 2;

double t = 0;
double dt = 0;

static const int maxStencil[2][3] = {-3, -3, 0, +4, +4, 1};

const double compression_tolerance = 5e-5;
const double refinement_tolerance = compression_tolerance / 10;

typedef Wavelets_Interp4thOrder W;
struct WARPBlock : public Block<PS, blockSize, blockSize, 1> {
  typedef PS ElementType;
  Real maximumDT;
  WARPBlock(ElementType e = ElementType())
      : Block<PS, blockSize, blockSize, 1>(e) {}
};

typedef WARPBlock B;
typedef Multithreading::BlockProcessing_TBB<B> BlockProcessing;

#ifdef _MRAG_GLUT_VIZ
GridViewer viewer(!W::bIsCellCentered, true);
RGBA convertToRGBA(PS &p) {
  const double c = p.w;
  const double R = max(0., min(1., 2 * (c - 0.5)));
  const double G = max(0., min(1., 1 - 2 * fabs(c - 0.5)));
  const double B = max(0., min(1., 1 - 2 * c));
  RGBA color(R, G, B, 0);
  return color;
}
#endif

#include "blockprocessingWARP.h"

class DemoWARP {
  Grid<W, B> *grid;
  BlockProcessing blockProcessing;
  Refiner_SpaceExtension refiner;
  Compressor compressor;
  BlockFWT<W, B, WARP_projector> blockfwt;
  SpaceTimeSorter stSorter;
  Real time;

  static const Real HeavySide(double x0, double x) {
    const Real eps = 0.1;
    const Real alpha = M_PI * min(1., max(0., (x - x0 + 0.5 * eps) / eps));

    return 0.5 + 0.5 * cos(alpha);
  }

  static void _getVelocity(float x[3], float u[2], Real time) {
    const Real factor = time > 3 ? -1 : 1;
    const double p[3] = {2 * M_PI * (x[0] - 0.00), 2 * M_PI * (x[1] - 0.00),
                         2 * M_PI * (x[2] - 0.00)};
    u[0] = -factor * 2.0 * pow(sin(p[0] * 0.5), 2) * sin(p[1] * 0.5) *
           cos(p[1] * 0.5);
    u[1] = 2.0 * factor * pow(sin(p[1] * 0.5), 2) * sin(p[0] * 0.5) *
           cos(p[0] * 0.5);
    if (time > 6)
      exit(0);
  }

  static void _ic_func4(float x[3], float &w, float u[2], Real time) {
    const float r1 = sqrt(pow((x[0] - 0.5) / 2, 2) + pow(x[1] - 0.5, 2));
    const double d1 = r1 - 0.15;

    w = HeavySide(0, d1);
    // u[0] = -(x[1]-0.5);
    // u[1] = +(x[0]-0.5);
    {
      _getVelocity(x, u, time);
    }
  }

  static void _ic(Grid<W, B> &grid) {
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    for (int i = 0; i < vInfo.size(); i++) {
      BlockInfo &info = vInfo[i];
      B &block = grid.getBlockCollection().lock(info.blockID);

      for (int iy = 0; iy < B::sizeY; iy++)
        for (int ix = 0; ix < B::sizeX; ix++) {
          float x[3];
          info.pos(x, ix, iy, 0);

          PS &point = block(ix, iy, 0);
          _ic_func4(x, point.w, point.u, 0);
        }

      grid.getBlockCollection().release(info.blockID);
    }
  }

public:
  DemoWARP() : refiner(resJump), compressor(resJump) {
    grid =
        new Grid<W, B>(blocksPerDimension, blocksPerDimension, 1, maxStencil);
    grid->setCompressor(&compressor);
    grid->setRefiner(&refiner);
    stSorter.connect(*grid);
    _ic(*grid);
    Science::AutomaticRefinement<0, 0>(*grid, blockfwt, refinement_tolerance,
                                       maxLevel, 1, _ic);
    _ic(*grid);
    Science::AutomaticCompression<0, 0>(*grid, blockfwt, compression_tolerance,
                                        1, _ic);
  }

  void simulation_init() { printf("WARP init\n"); }

  void simulation_run(int nsteps) {
    //	printf("WARP run\n");
    Science::AutomaticRefinement<0, 0>(*grid, blockfwt, refinement_tolerance,
                                       maxLevel, -1);

    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    BoundaryInfo &sBoundaryInfo = grid->getBoundaryInfo();
    const BlockCollection<B> &collection = grid->getBlockCollection();

    CLFParticle<B::sizeX>::FindMaximumDT findmaximumDT;
    BlockProcessing::process(vInfo, collection, findmaximumDT);

    Real dt = HUGE_VALF;
    for (int i = 0; i < vInfo.size(); i++) {
      B &b = collection.lock(vInfo[i].blockID);
      dt = min(dt, b.maximumDT);
      collection.release(vInfo[i].blockID);
    }

    printf("T=%f, minimum time step is %f\n", time, dt);

    for (int istep = 0; istep < 10; istep++) {
      CLFParticle<B::sizeX>::AdvectParticlesRK3 advectparticles(dt);
      blockProcessing.process<CLFParticle<B::sizeX>::BlockLabWARP>(
          vInfo, collection, sBoundaryInfo, advectparticles);

      CLFParticle<B::sizeX>::UpdateOmega updateomega;
      BlockProcessing::process(vInfo, collection, updateomega);

      time += dt;
    }

    Science::AutomaticCompression<0, 0>(*grid, blockfwt, compression_tolerance, 1);
  }

  void simulation_render(bool bDrawTextures) {
#ifdef _MRAG_GLUT_VIZ
    viewer.drawContent(*grid, grid->getBlockCollection());
    viewer.drawSketch(*grid, false);
#endif
  }
};

#ifdef _MRAG_GLUT_VIZ
#ifndef __APPLE__
#include "GL/glut.h"
#else
#include "GLUT/glut.h"
#endif
// #include "screenshot.h"
extern "C" GLint gltWriteTGA(const char *szFileName);
#endif

const int nMaxSteps = 100000;
int iCurrStep = 0;
DemoWARP *demo = NULL;
#ifdef _MRAG_GLUT_VIZ

static void display(void) {
  // glClearColor(1,0,1,0);
  glClear(GL_COLOR_BUFFER_BIT);

  demo->simulation_render(true);

  const bool bStop = (iCurrStep++ >= nMaxSteps);
  if (bStop)
    exit(0);
  glFinish();

  static int iFrameCounter = 0;
  char buf[300];
  sprintf(buf, "asdasd%05d.tga", iFrameCounter++);
  printf("Writing to %s\n", buf);
  gltWriteTGA(buf);
  glutSwapBuffers();
  demo->simulation_run(1);
}

static void idle(void) { glutPostRedisplay(); }
#endif

int main(int argc, char **argv) {
  srand(3290);

#ifdef _MRAG_GLUT_VIZ
  {
    glutInit(&argc, argv);
    glutInitWindowSize(800, 800);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_STENCIL | GLUT_RGBA | GLUT_DOUBLE);

    glutCreateWindow("WARP");

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(-0.2, 1.2, -0.2, 1.2, -1, 1);
    glMatrixMode(GL_MODELVIEW);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glEnable(GL_TEXTURE_2D);

    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    glutDisplayFunc(display);
    glutIdleFunc(idle);
  }
#endif
  Environment::setup();
  demo = new DemoWARP;
  demo->simulation_init();

#ifdef _MRAG_GLUT_VIZ
  glutMainLoop();
#else
  const int nsteps = 4;
  for (int i = 0; i < nsteps; i++)
    demo->simulation_run(1);
#endif
}
