/*
 *  I2D_Headers.h
 *  IncompressibleFluids2D
 *
 *  Created by Diego Rossinelli on 11/05/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include <xmmintrin.h>

#include "MRAGCommon.h"
#include "MRAGEnvironment.h"

#ifdef __APPLE__
#ifdef _MRAG_GLUT_VIZ
#include "GLUT/glut.h"
#endif
#endif

#undef min
#undef max

#include "MRAGWavelets_AverageInterp5thOrder.h"
#include "MRAGWavelets_Interp4thOrder.h"
#include "MRAGWavelets_AverageInterp3rdOrder.h"
#include "MRAGWavelets_Interp2ndOrder.h"
#include "MRAGWavelets_Haar.h"
#include "MRAGrid.h"
#include "MRAGrid_Hilbert2D.h"
#include "MRAGRefiner.h"
#include "MRAGRefiner_Greedy.h"
#include "MRAGRefiner_BlackList.h"
#include "MRAGCompressor.h"
#include "MRAGBlockLab.h"
#include "MRAGBlockFWT.h"
#include "MRAGBlock.h"
#include "QuadTree.h"

#ifdef _MRAG_GLUT_VIZ
#include "GridViewer.h"
#endif

#include "MRAGScienceCore.h"
#include "MRAGAutomaticRefiner.h"
#include "MRAGSimpleLevelsetBlock.h"
#include "MRAGSpaceTimeSorter.h"
#include "candidate_SpaceTimeSorterRK2.h"
#include "MRAGRefiner_SpaceExtension.h"

#include "MRAGBlockProcessing_SingleCPU.h"
#include "MRAGBlockProcessing_TBB.h"

#include "MRAG_IO_ArgumentParser.h"
#include "MRAG_IO_Binary.h"
#include "MRAG_IO_VTKNative.h"
