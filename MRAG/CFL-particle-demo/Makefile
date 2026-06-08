TBB_FLAGS  = $(shell pkg-config --cflags tbb)
TBB_LIBS   = $(shell pkg-config --libs tbb) -ltbbmalloc
TBB_LIBDIR = $(shell pkg-config --variable=libdir tbb)
CXXFLAGS ?= -Ofast -DNDEBUG -fopenmp -D_MRAG_GLUT_VIZ
test: mainWARP.cpp screenshot.cpp ../MRAGcore/*.cpp
	$(CXX) $(CXXFLAGS) $(TBB_FLAGS) $^ -o test -I../ -lGL -ltbb -ltbbmalloc  -lglut
clean:
	rm -f test
