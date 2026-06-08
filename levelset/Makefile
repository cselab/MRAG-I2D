SRC = ../src
TBB_FLAGS  = $(shell pkg-config --cflags tbb)
TBB_LIBS   = $(shell pkg-config --libs tbb) -ltbbmalloc
TBB_LIBDIR = $(shell pkg-config --variable=libdir tbb)
CORE_CPP = \
  $(SRC)/MRAGWavelets_StaticData.cpp \
  $(SRC)/MRAGBoundaryBlockInfo.cpp
onecase: onecase.cpp $(CORE_CPP)
	$(CXX) $(CXXFLAGS) -I$(SRC) $(TBB_FLAGS) -o $@ $^ $(TBB_LIBS)
run: onecase
	DYLD_LIBRARY_PATH=$(TBB_LIBDIR) ./onecase
clean:
	rm -f onecase phi.xdmf2 phi.*.xdmf2 phi.*.xyz.raw phi.*.phi.raw
.PHONY: run clean
