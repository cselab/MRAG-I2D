# Minimal MRAG-I2D wavelet core driver (oneTBB).
# Apple-Silicon clang / C++03-era code.
#
#   make        -> build ./onecase
#   make run    -> build + run

SRC_MRAG = ../../MRAG/source
CXX      = clang++

WARN     = -Wno-deprecated -Wno-unknown-pragmas -Wno-reorder -Wno-comment \
           -Wno-writable-strings -Wno-c++11-compat-deprecated-writable-strings
CPPFLAGS = -I$(SRC_MRAG) $(WARN)
CXXFLAGS = -O2 -g

# oneTBB via pkg-config. tbb.pc only lists -ltbb, so append -ltbbmalloc
# (needed for tbb::scalable_allocator).
TBB_FLAGS  = $(shell pkg-config --cflags tbb)
TBB_LIBS   = $(shell pkg-config --libs tbb) -ltbbmalloc
TBB_LIBDIR = $(shell pkg-config --variable=libdir tbb)

CORE_CPP = \
  $(SRC_MRAG)/MRAGcore/MRAGWavelets_StaticData.cpp \
  $(SRC_MRAG)/MRAGcore/MRAGBoundaryBlockInfo.cpp \
  $(SRC_MRAG)/MRAGcore/MRAGProfiler.cpp

onecase: onecase.cpp $(CORE_CPP)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TBB_FLAGS) -o $@ $^ $(TBB_LIBS)

run: onecase
	DYLD_LIBRARY_PATH=$(TBB_LIBDIR) ./onecase

clean:
	rm -f onecase phi_refined.txt
	rm -rf onecase.dSYM

.PHONY: run clean
