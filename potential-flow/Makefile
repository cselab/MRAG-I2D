# Potential-flow (PF) app built on the MRAG core, oneTBB.
# Headless by default; build with VIZ=1 for the GLUT viewer.
#
#   make            -> build ./pf
#   make VIZ=1      -> build with the GLUT viewer (_RL_VIZ)
#   make run        -> build + run a sample study

SRC  = ../src
CXX  = g++-15

WARN = -Wno-deprecated -Wno-unknown-pragmas -Wno-reorder -Wno-comment
CXXFLAGS = -std=gnu++17 -O2 -DNDEBUG -fopenmp -DGL_SILENCE_DEPRECATION $(WARN) -include cstring -include cassert
INCLUDES = -I$(SRC) $(shell pkg-config --cflags tbb)

# scalable_allocator lives in libtbbmalloc; tbb.pc only lists -ltbb
TBB_LIBS   = $(shell pkg-config --libs tbb) -ltbbmalloc
TBB_LIBDIR = $(shell pkg-config --variable=libdir tbb)
GL_LIBS    = -framework OpenGL -framework GLUT

ifdef VIZ
CXXFLAGS += -D_RL_VIZ
endif

PF_SRC = $(wildcard $(SRC)/PF_*.cpp)
RL_SRC = $(SRC)/RL_QLearning.cpp $(SRC)/RL_MultiTable.cpp $(SRC)/rng.cpp
CORE   = $(SRC)/MRAGWavelets_StaticData.cpp \
         $(SRC)/MRAGBoundaryBlockInfo.cpp \
         $(SRC)/QtBox.cpp \
         $(SRC)/QuadTree.cpp

pf: $(PF_SRC) $(RL_SRC) $(CORE)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(TBB_LIBS) $(GL_LIBS)

# full LEARNING_DIPOLE study args (from launch/natalieLaunchPF.sh); reads ./factoryPF.
# This is a long RL training run (learningTime 100000) -- Ctrl-C to stop.
STUDY = -study LEARNING_DIPOLE -nthreads 1 -D 0.015 -xpos 0.5 -ypos 0.5 \
        -lr 0.01 -gamma 0.98 -greedyEps 0.01 -shared 1 -isControlled 1 -smooth 0 \
        -savefreq 100 -learnDump 1000 -learningTime 100000 -nLearningLevels 1 \
        -fitnessTime 0 -navg 1 -fitnessSaveFreq 100000 -fitnessBuffer 1.5 \
        -individualFitness true -exploitation true -factory factoryPF \
        -fitselection 1 -islabframe 1 -isUsingTracers 0 -restart 0

run: pf factoryPF
	DYLD_LIBRARY_PATH=$(TBB_LIBDIR) ./pf $(STUDY)

factoryPF:
	cp ../I2D/launch/factoryPF factoryPF

clean:
	rm -f pf

# remove run artifacts (keeps factoryPF)
clean-run:
	rm -f RestartTime.txt argumentparser.log data.txt fitness \
	      learning_PF_* saveQ_PF_*

.PHONY: run clean clean-run
