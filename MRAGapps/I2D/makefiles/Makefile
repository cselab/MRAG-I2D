hn ?= $(shell hostname)
username ?= $(shell whoami)

ifeq "$(hn)" "taylor.ethz.ch"
include make.angeliki
endif		

ifeq "$(hn)" "olli.inf.ethz.ch"
include make.olli
endif

ifeq "$(hn)" "natalie.ethz.ch"
include make.mattia
endif

ifeq "$(hn)" "Mattias-MacBook-Pro.local"
include make.mattia
endif

ifeq "$(hn)" "public-docking-nw-4-122.ethz.ch"
include make.mattia
endif

ifeq "$(hn)" "macross.inf.ethz.ch"
include make.mattia
endif

ifeq "$(hn)" "sequoia.inf.ethz.ch"
include make.giovanni
endif

ifeq "$(hn)" "jerome.zyxel.com"
include make.atchieu.jerome
endif

ifeq "$(hn)" "marcel.ethz.ch"
include make.atchieu.marcel
endif

ifeq "$(hn)" "vpn-global-111-dhcp.ethz.ch"
include make.rosso
endif

ifeq "$(hn)" "rosso.inf.ethz.ch"
include make.rosso
endif

ifeq "$(hn)" "xenon.inf.ethz.ch"
include make.roman
endif

ifneq "$(findstring brutus,$(hn))" ""
	ifeq "$(username)" "mgazzola"
		include make.brutus
	endif

	ifeq "$(username)" "alexiad"
		include make.alexiaBrutus
	endif

	ifeq "$(username)" "atchieu"
		include make.atchieu.brutus
	endif

	ifeq "$(username)" "alaskari"
		include make.angelikiBrutus
	endif

        ifeq "$(username)" "tosimon"
	        include make.simon
        endif

        ifeq "$(username)" "hbabak"
	        include make.brutus.hbabak
        endif

        ifeq "$(username)" "wvanrees"
	        include make.brutus.wvanrees
        endif

endif

ap ?= double
bs ?= 32
fmm-order ?= 12
fmm-lmax ?= 9
fmm-smax ?= 200
fmm-sse-noprecdiv ?=1
fmm-sse-mixedprec ?= 1

nthreads ?= 48

CPPFLAGS += -D_BLOCKSIZE_=$(bs) -D_MRAG_TBB_NTHREADS_HINT=$(nthreads) -D_MRAG_REAL_=$(ap)
CPPFLAGS += -D_ORDER_=$(fmm-order) -D_FMM_MAX_LEVEL_=$(fmm-lmax) -D_SMAX=$(fmm-smax) 

ifeq "$(config)" "production"
CPPFLAGS += -DNDEBUG
else
CPPFLAGS += -g
endif 

ifeq "$(fmm-sse-mixedprec)" "1"
CPPFLAGS += -D_FMM_MIXEDPREC_KERNELS_
endif 

ifeq "$(fmm-sse-noprecdiv)" "1"
CPPFLAGS += -D_FMM_NOPRECDIV_KERNELS_
endif

CPPFLAGS+= -DVTK_EXCLUDE_STRSTREAM_HEADERS

CPPFLAGS += -I../source/ -I../source/mani-fmm2d/ -I../source/mattia-RL/ -I../source/potential-flow/ -I../../../MRAG/source/

LIBS ?= -lvtkHybrid \
	-lvtkVolumeRendering \
	-lvtkRendering \
	-lvtkIO \
	-lvtkGenericFiltering \
	-lvtkGraphics \
	-lvtkImaging \
	-lvtkFiltering \
	-lvtkCommon \
	-lvtkftgl \
	-lvtkDICOMParser \
	-lvtksys 

LIBS += -lstdc++ -ltbb -ltbbmalloc

TESTOBJS = \
	I2D_main.o \
	MRAGBoundaryBlockInfo.o \
	MRAGProfiler.o \
	MRAGWavelets_StaticData.o \
	binomial.o \
	QtBox.o \
	QuadTree.o \
	RL_MultiTable.o \
	RL_QLearning.o \
	rng.o \
	I2D_TestDumping.o \
	I2D_CircularObstacleOperator.o \
	I2D_TestDiffusion.o \
	I2D_DiffusionOperator_4thOrder.o \
	I2D_TestAdvection.o \
	I2D_AdvectionOperator.o \
	I2D_AdvectionOperator_Particles.o \
	I2D_TestPenalizationAndOther.o \
	I2D_PenalizationOperator.o \
	I2D_DivOperator.o \
	I2D_VelocitySolver_Mani.o \
	I2D_CoreFMM_AggressiveVel.o \
	I2D_AggressiveDiego.o \
	I2D_CoreFMM_SSE.o \
	I2D_CoreFMM_Check.o \
	I2D_CoreFMM_PlanBuilder.o \
	I2D_CoreFMM_PlanBuilderWim.o \
	I2D_FMMTypes.o \
	I2D_TestPoissonEquation.o \
	I2D_CurlVelocityOperator_4thOrder.o \
	I2D_CurlVelocityOperator_2ndOrder.o \
	I2D_WingObstacleOperator.o \
	I2D_Axisymmetrization.o \
	I2D_TestPoissonEquationPotential.o \
	I2D_PotentialSolver_Mattia.o \
	I2D_TestMultipole.o \
	I2D_VelocitySolver_Wim.o \
	I2D_Clear.o

NOMINEPATRIS_OBJS = \
	MRAGBoundaryBlockInfo.o \
	MRAGProfiler.o \
	MRAGWavelets_StaticData.o \
	binomial.o \
	QtBox.o \
	QuadTree.o \
	RL_MultiTable.o \
	RL_QLearning.o \
	rng.o \
	I2D_CircularObstacleOperator.o \
	I2D_RotatingWheelObstacle.o \
	I2D_RectangularObstacleOperator.o \
	I2D_DiffusionOperator_4thOrder.o \
	I2D_AdvectionOperator.o \
	I2D_AdvectionOperator_Particles.o \
	I2D_PenalizationOperator.o \
	I2D_DivOperator.o \
	I2D_FlowPastFixedObstacle.o \
	I2D_VelocitySolver_Mani.o \
	I2D_VelocitySolver_Wim.o \
	I2D_CoreFMM_AggressiveVel.o \
	I2D_AggressiveDiego.o \
	I2D_CoreFMM_SSE.o \
	I2D_CoreFMM_Check.o \
	I2D_CoreFMM_PlanBuilder.o \
	I2D_CoreFMM_PlanBuilderWim.o \
	I2D_FMMTypes.o \
	I2D_AveMaria.o \
	I2D_WingObstacleOperator.o \
	I2D_LinkedBodies.o \
	I2D_CurlVelocityOperator_4thOrder.o \
	I2D_EllipticalObstacleOperator.o \
	I2D_FlowPastFloatingObstacle.o \
	I2D_FloatingCylinder.o \
	I2D_FloatingRotatingCylinderPair.o \
	I2D_FloatingObstacleOperator.o \
	I2D_FloatingObstacleVector.o \
	I2D_Clear.o \
	I2D_ImposedCylinder.o \
	I2D_CarlingFish.o \
	I2D_Axisymmetrization.o

CHLOE_OBJS = \
	MRAGBoundaryBlockInfo.o \
	MRAGProfiler.o \
	MRAGWavelets_StaticData.o \
	binomial.o \
	QtBox.o \
	QuadTree.o \
	RL_MultiTable.o \
	RL_QLearning.o \
	rng.o \
	I2D_CircularObstacleOperator.o \
	I2D_FlowPastFixedObstacle.o \
	I2D_KillVortRightBoundaryOperator.o \
	I2D_RectangularObstacleOperator.o \
	I2D_DiffusionOperator_4thOrder.o \
	I2D_AdvectionOperator.o \
	I2D_AdvectionOperator_Particles.o \
	I2D_PenalizationOperator.o \
	I2D_DivOperator.o \
	I2D_MRAGOptimisation.o \
	I2D_VelocitySolver_Mani.o \
	I2D_VelocitySolver_Wim.o \
	I2D_CoreFMM_AggressiveVel.o \
	I2D_AggressiveDiego.o \
	I2D_CoreFMM_SSE.o \
	I2D_CoreFMM_Check.o \
	I2D_CoreFMM_PlanBuilder.o \
	I2D_CoreFMM_PlanBuilderWim.o \
	I2D_FMMTypes.o \
	I2D_Chloe.o \
	I2D_WingObstacleOperator.o \
	I2D_Ingo.o \
	I2D_LinkedBodies.o \
	I2D_WingAngleOfAttack.o \
	I2D_CurlVelocityOperator_4thOrder.o \
	I2D_EllipticalObstacleOperator.o \
	I2D_ImposedCylinder.o \
	I2D_FloatingCylinder.o \
	I2D_FloatingRotatingCylinderPair.o \
	I2D_FloatingObstacleVector.o \
	I2D_Clear.o \
	I2D_CarlingFish.o \
	I2D_StefanFish.o \
	I2D_ThreeLinkFish.o \
	I2D_FlowPastFloatingObstacle.o \
	I2D_FluidMediatedInteractions.o \
	I2D_RotatingCylinderPair.o \
	I2D_ImposedCylinderOscillation.o \
	I2D_ObjectFactory.o \
	I2D_ImposedEllipse.o \
	I2D_ImposedEllipseOscillation.o \
	I2D_ImposedEllipseRotation.o \
	I2D_ImposedCylinderRotation.o \
	I2D_FloatingEllipse.o \
	I2D_RotatingWheelObstacle.o \
	I2D_FloatingObstacleOperator.o \
	I2D_ImposedCylinderStopAndGo.o \
	I2D_PotentialSolver_Mattia.o \
	I2D_CStartLarva.o \
	I2D_CarlingFishMorph.o \
	I2D_StefanFishMorph.o \
	I2D_FTLE.o \
	I2D_SmartInteractions.o \
	I2D_StefanFishSmart.o \
	I2D_CarlingFishAirfoil.o \
	I2D_StefanFishSmartInline.o \
	I2D_StefanFishSmartLattice.o \
	I2D_PassiveTracer.o \
	I2D_StefanFishSmartSpiral.o \
	I2D_ComputeEfficiency.o \
	I2D_ImposedAirfoil.o \
	I2D_PitchingAirfoil.o \

RL_OBJS = \
	MRAGProfiler.o \
	RL_main.o \
	RL_TestTabular.o \
	RL_MultiTable.o \
	RL_QLearning.o \
	rng.o \
	RL_Agent.o \
	RL_AgentVector.o \
	RL_ObjectFactory.o \
	RL_Food.o \
	RL_SmartyGlutton.o \

ALEXIA_OBJS = \
	MRAGProfiler.o \
	PF_main.o \
	RL_MultiTable.o \
	RL_QLearning.o \
	rng.o \
	PF_AgentVector.o \
	PF_Solver.o \
	PF_ObjectFactory.o \
	PF_Dipole.o \
	PF_Agent.o \
	PF_DipoleInline.o \
	PF_SwarmCylinders.o \
	PF_DipoleInlineBravais.o \
	PF_DipoleEight.o

ifeq "$(mpi)" "1"
	TESTOBJS += \
		I2D_VelocitySolverMPI_Mani.o

	NOMINEPATRIS_OBJS += \
		I2D_VelocitySolverMPI_Mani.o

	CHLOE_OBJS += \
		I2D_VelocitySolverMPI_Mani.o

	override CC=mpicxx

	CPPFLAGS += -D_I2D_MPI_
endif


VPATH := ../source/ ../source/mani-fmm2d/ ../source/mattia-RL/ ../source/potential-flow/ ../../../MRAG/source/MRAGcore
.DEFAULT_GOAL := avemaria


if2d-tests: $(TESTOBJS)
	$(CC)  $^ $(LIBS) -o $@
	@echo done

avemaria: $(NOMINEPATRIS_OBJS)
	$(CC)  $^ $(LIBS) -o $@
	@echo done

chloe: $(CHLOE_OBJS)
	$(CC) $(LIBS) $^ -o $@
	@echo done

alexia: $(ALEXIA_OBJS)
	$(CC) $(LIBS) $^ -o $@
	@echo done

rl: $(RL_OBJS)
	$(CC) $(LIBS) $^ -o $@
	@echo done

#I2D_CoreFMM_AggressiveVel.o: I2D_CoreFMM_AggressiveVel.cpp
#	echo SOURCES = $@
#	$(CC ) $(CPPFLAGS) -axTW -ip -c $^ -o $@

ifeq "$(CC)" "icc"
KERNELOPTFLAGS = -O3 -xHOST -ip -ansi-alias -fno-fnalias -fno-alias -falign-functions -diag-disable remark -wd68
else
KERNELOPTFLAGS = -O3 -fstrict-aliasing -msse3 -Wno-missing-braces -fdiagnostics-show-option
endif

I2D_AggressiveDiego.o: I2D_AggressiveDiego.cpp
	echo SOURCES = $@
	$(CC) $(CPPFLAGS) $(KERNELOPTFLAGS) -O3 -c $^ -o $@
	$(CC) $(CPPFLAGS) $(KERNELOPTFLAGS) -O3 -c $^ -S 

%.o: %.cpp
	$(CC)  $(CPPFLAGS) $(CPPFLAGSOPT) -c $^ -o $@

clean:
	rm -f *.o
	rm -f *.s
	rm -f avemaria if2d-tests rl chloe alexia
