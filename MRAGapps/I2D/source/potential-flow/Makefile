# Potential-flow (PF) app built on the MRAG core, oneTBB.
# Headless by default; build with VIZ=1 for the GLUT viewer.
#
#   make            -> build ./pf
#   make VIZ=1      -> build with the GLUT viewer (_RL_VIZ)
#   make run        -> build + run a sample study

MRAG = ../../../../MRAG/source
RL   = ../mattia-RL
CXX  = g++-15

WARN = -Wno-deprecated -Wno-unknown-pragmas -Wno-reorder -Wno-comment
CXXFLAGS = -std=gnu++17 -O2 -DNDEBUG -fopenmp -DGL_SILENCE_DEPRECATION $(WARN)
INCLUDES = -I$(MRAG) -I.. -I$(RL) $(shell pkg-config --cflags tbb)

# scalable_allocator lives in libtbbmalloc; tbb.pc only lists -ltbb
TBB_LIBS   = $(shell pkg-config --libs tbb) -ltbbmalloc
TBB_LIBDIR = $(shell pkg-config --variable=libdir tbb)
GL_LIBS    = -framework OpenGL -framework GLUT

ifdef VIZ
CXXFLAGS += -D_RL_VIZ
endif

PF_SRC = PF_main.cpp PF_Agent.cpp PF_AgentVector.cpp PF_Dipole.cpp \
         PF_DipoleEight.cpp PF_DipoleInline.cpp PF_DipoleInlineBravais.cpp \
         PF_ObjectFactory.cpp PF_Solver.cpp PF_SwarmCylinders.cpp
RL_SRC = $(RL)/RL_QLearning.cpp $(RL)/RL_MultiTable.cpp $(RL)/rng.cpp
CORE   = $(MRAG)/MRAGcore/MRAGWavelets_StaticData.cpp \
         $(MRAG)/MRAGcore/MRAGBoundaryBlockInfo.cpp \
         $(MRAG)/MRAGcore/MRAGProfiler.cpp \
         $(MRAG)/MRAGcore/QtBox.cpp \
         $(MRAG)/MRAGcore/QuadTree.cpp

pf: $(PF_SRC) $(RL_SRC) $(CORE)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(TBB_LIBS) $(GL_LIBS)

run: pf
	DYLD_LIBRARY_PATH=$(TBB_LIBDIR) ./pf -study LEARNING_DIPOLE

clean:
	rm -f pf

.PHONY: run clean
