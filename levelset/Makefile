SRC_MRAG = ../MRAG
TBB_FLAGS  = $(shell pkg-config --cflags tbb)
TBB_LIBS   = $(shell pkg-config --libs tbb) -ltbbmalloc
TBB_LIBDIR = $(shell pkg-config --variable=libdir tbb)
CORE_CPP = \
  $(SRC_MRAG)/MRAGcore/MRAGWavelets_StaticData.cpp \
  $(SRC_MRAG)/MRAGcore/MRAGBoundaryBlockInfo.cpp
onecase: onecase.cpp $(CORE_CPP)
	$(CXX) $(CXXFLAGS) -I$(SRC_MRAG) $(TBB_FLAGS) -o $@ $^ $(TBB_LIBS)
run: onecase
	DYLD_LIBRARY_PATH=$(TBB_LIBDIR) ./onecase
clean:
	rm -f onecase phi.xdmf2 phi.*.xdmf2 phi.*.xyz.raw phi.*.phi.raw
.PHONY: run clean
