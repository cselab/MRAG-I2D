CXXFLAGS ?= -Ofast -DNDEBUG -fopenmp -D_MRAG_GLUT_VIZ

test: mainWARP.cpp screenshot.cpp ../MRAGcore/*.cpp
	$(CXX) $(CXXFLAGS) $^ -o test -I../ -lGL -ltbb -ltbbmalloc  -lglut

clean:
	rm -f test

.PHONY = clean
