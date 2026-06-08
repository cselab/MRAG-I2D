TBB_FLAGS  = $(shell pkg-config --cflags tbb --libs tbb glu GL) -ltbbmalloc
TBB_LIBDIR = $(shell pkg-config --variable=libdir tbb)
CXX = g++-15
CXXFLAGS ?= -O3 --fast-math -DNDEBUG -fopenmp -D_MRAG_GLUT_VIZ -DGL_SILENCE_DEPRECATION
test: mainWARP.cpp screenshot.cpp ../MRAGcore/*.cpp
	$(CXX) $(CXXFLAGS) $(TBB_FLAGS) $^ -o test -I..
clean:
	rm -f test
