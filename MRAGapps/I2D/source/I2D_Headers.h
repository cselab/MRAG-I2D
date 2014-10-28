/*
 *  I2D_Headers.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include <xmmintrin.h>

#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGEnvironment.h"

#ifdef __APPLE__
#ifdef _MRAG_GLUT_VIZ
#include "GLUT/glut.h"
#endif
#endif

#undef min
#undef max

#include "MRAGcore/MRAGWavelets_AverageInterp5thOrder.h"
#include "MRAGcore/MRAGWavelets_Interp4thOrder.h"
#include "MRAGcore/MRAGWavelets_AverageInterp3rdOrder.h"
#include "MRAGcore/MRAGWavelets_Interp2ndOrder.h"
#include "MRAGcore/MRAGWavelets_Haar.h"
#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGrid_Hilbert2D.h"
#include "MRAGcore/MRAGRefiner.h"
#include "MRAGcore/MRAGRefiner_Greedy.h"
#include "MRAGcore/MRAGRefiner_BlackList.h"
#include "MRAGcore/MRAGCompressor.h"
#include "MRAGcore/MRAGBlockLab.h"
#include "MRAGcore/MRAGBlockFWT.h"
#include "MRAGcore/MRAGBlock.h"
#include "MRAGcore/MRAGProfiler.h"
#include "MRAGcore/QuadTree.h"

#ifdef _MRAG_GLUT_VIZ
#include "MRAGvisual/GridViewer.h"
#endif

#include "MRAGscience/MRAGScienceCore.h"
#include "MRAGscience/MRAGAutomaticRefiner.h"
#include "MRAGscience/MRAGSimpleLevelsetBlock.h"
#include "MRAGscience/MRAGSpaceTimeSorter.h"
#include "MRAGscience/candidate_SpaceTimeSorterRK2.h"
#include "MRAGscience/MRAGRefiner_SpaceExtension.h"

#include "MRAGmultithreading/MRAGBlockProcessing_SingleCPU.h"
#include "MRAGmultithreading/MRAGBlockProcessing_TBB.h"

#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "MRAGio/MRAG_IO_Binary.h"
#include "MRAGio/MRAG_IO_VTKNative.h"

