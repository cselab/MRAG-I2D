SRC_MRAG = ../../MRAG/source
TBB_FLAGS  = $(shell pkg-config --cflags tbb)
TBB_LIBS   = $(shell pkg-config --libs tbb) -ltbbmalloc
TBB_LIBDIR = $(shell pkg-config --variable=libdir tbb)
CORE_CPP = \
  $(SRC_MRAG)/MRAGcore/MRAGWavelets_StaticData.cpp \
  $(SRC_MRAG)/MRAGcore/MRAGBoundaryBlockInfo.cpp \
  $(SRC_MRAG)/MRAGcore/MRAGProfiler.cpp
onecase: onecase.cpp $(CORE_CPP)
	$(CXX) $(CXXFLAGS) -I$(SRC_MRAG) $(TBB_FLAGS) -o $@ $^ $(TBB_LIBS)
run: onecase
	DYLD_LIBRARY_PATH=$(TBB_LIBDIR) ./onecase
clean:
	rm -f onecase phi_refined.txt phi_refined.xdmf2 phi_refined.xyz.raw phi_refined.phi.raw
.PHONY: run clean
