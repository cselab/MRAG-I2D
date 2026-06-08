# Full I2D vortex solver (I2D_main), oneTBB, headless.
# SSE intrinsics are mapped to NEON via ./sse-compat (sse2neon) on Apple Silicon.
# VTK output is stubbed out unless _MRAG_VTK is defined (no VTK install needed).
#
#   make -f Makefile.i2d            -> build ./i2d
#   make -f Makefile.i2d run        -> build + run a short advection study

MRAG = ../MRAG
CXX  = g++-15

DEFS = -D_BLOCKSIZE_=32 -D_ORDER_=12 -D_FMM_MAX_LEVEL_=9 -D_SMAX=200 \
       -D_FMM_MIXEDPREC_KERNELS_ -D_FMM_NOPRECDIV_KERNELS_
CXXFLAGS = -std=gnu++17 -O2 -DNDEBUG -fopenmp -DGL_SILENCE_DEPRECATION $(DEFS)
INCLUDES = -Isse-compat -I/opt/homebrew/include -I$(MRAG) -I. -I./mani-fmm2d -I./mattia-RL \
           $(shell pkg-config --cflags tbb)

TBB_LIBS   = $(shell pkg-config --libs tbb) -ltbbmalloc
TBB_LIBDIR = $(shell pkg-config --variable=libdir tbb)

# top-level I2D sources, minus the alternate main()s
I2D_SRC = $(filter-out I2D_AveMaria.cpp I2D_Chloe.cpp I2D_Mattia.cpp, $(wildcard I2D_*.cpp))
FMM_SRC = mani-fmm2d/binomial.cpp mani-fmm2d/precompute.cpp
RL_SRC  = mattia-RL/RL_QLearning.cpp mattia-RL/RL_MultiTable.cpp mattia-RL/rng.cpp
CORE    = $(MRAG)/MRAGcore/MRAGWavelets_StaticData.cpp \
          $(MRAG)/MRAGcore/MRAGBoundaryBlockInfo.cpp \
          $(MRAG)/MRAGcore/MRAGProfiler.cpp \
          $(MRAG)/MRAGcore/QtBox.cpp \
          $(MRAG)/MRAGcore/QuadTree.cpp

i2d: $(I2D_SRC) $(FMM_SRC) $(RL_SRC) $(CORE)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(TBB_LIBS)

# short advection study (level-set / scalar advection on the adapted grid)
ADV = -study advection -nsteps 50 -nthreads 1 -bpd 8 -lmax 4 \
      -rtol 1e-4 -ctol 1e-5 -dumpfreq 5 -jump 1 -particles 0 -uniform 0

run: i2d
	DYLD_LIBRARY_PATH=$(TBB_LIBDIR) ./i2d $(ADV)

clean:
	rm -f i2d advection.xdmf2 advection_*.xdmf2 advection_*.raw

.PHONY: run clean
