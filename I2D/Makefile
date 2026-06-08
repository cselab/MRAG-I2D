SRC  = ../src
CXX  = g++-15

DEFS = -D_BLOCKSIZE_=32 -D_ORDER_=12 -D_FMM_MAX_LEVEL_=9 -D_SMAX=200 \
       -D_FMM_MIXEDPREC_KERNELS_ -D_FMM_NOPRECDIV_KERNELS_
CXXFLAGS = -std=gnu++17 -O2 -DNDEBUG -fopenmp -DGL_SILENCE_DEPRECATION $(DEFS)
INCLUDES = -I$(SRC)/sse-compat -I/opt/homebrew/include -I$(SRC) \
           $(shell pkg-config --cflags tbb)

TBB_LIBS   = $(shell pkg-config --libs tbb) -ltbbmalloc
TBB_LIBDIR = $(shell pkg-config --variable=libdir tbb)

# local entry point + the I2D solver library (no main()) from src/.
# The alternate main()s (I2D_AveMaria/Chloe/Mattia.cpp) also live here but are
# intentionally not compiled.
I2D_SRC = I2D_main.cpp $(wildcard $(SRC)/I2D_*.cpp)
FMM_SRC = $(SRC)/binomial.cpp $(SRC)/precompute.cpp
RL_SRC  = $(SRC)/RL_QLearning.cpp $(SRC)/RL_MultiTable.cpp $(SRC)/rng.cpp
CORE    = $(SRC)/MRAGWavelets_StaticData.cpp \
          $(SRC)/MRAGBoundaryBlockInfo.cpp \
          $(SRC)/QtBox.cpp \
          $(SRC)/QuadTree.cpp

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
