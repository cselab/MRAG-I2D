TBB_FLAGS  = $(shell pkg-config --cflags --libs tbb) -ltbbmalloc
TBB_LIBDIR = $(shell pkg-config --variable=libdir tbb)
GL_FRAMEWORKS = -framework GLUT -framework OpenGL
CXX = g++-15
CXXFLAGS ?= -O3 --fast-math -DNDEBUG -fopenmp -D_MRAG_GLUT_VIZ -DGL_SILENCE_DEPRECATION
main: main.cpp screenshot.cpp ../MRAG/MRAGcore/*.cpp
	$(CXX) $(CXXFLAGS) $(TBB_FLAGS) $^ -o $@ -I../MRAG $(GL_FRAMEWORKS)
clean:
	rm -f main
