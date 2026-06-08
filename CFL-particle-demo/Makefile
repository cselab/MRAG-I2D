SRC = ../src
TBB_FLAGS  = $(shell pkg-config --cflags --libs tbb) -ltbbmalloc
TBB_LIBDIR = $(shell pkg-config --variable=libdir tbb)
GL_FRAMEWORKS = -framework GLUT -framework OpenGL
CXX = g++-15
CXXFLAGS ?= -O3 --fast-math -DNDEBUG -fopenmp -D_MRAG_GLUT_VIZ -DGL_SILENCE_DEPRECATION
CORE = $(SRC)/MRAGWavelets_StaticData.cpp \
       $(SRC)/MRAGBoundaryBlockInfo.cpp \
       $(SRC)/QtBox.cpp \
       $(SRC)/QuadTree.cpp
main: $(SRC)/main.cpp $(SRC)/screenshot.cpp $(CORE)
	$(CXX) $(CXXFLAGS) $(TBB_FLAGS) $^ -o $@ -I$(SRC) $(GL_FRAMEWORKS)
clean:
	rm -f main
