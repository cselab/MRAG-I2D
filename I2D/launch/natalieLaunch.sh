#!/bin/bash

EXECNAME=$1

SETTINGS=
SETTINGS+=" -study SMART_INTERACTIONS"
SETTINGS+=" -bpd 16"
SETTINGS+=" -nthreads 2"

SETTINGS+=" -lambda 1e4"
SETTINGS+=" -tend 10.0"
SETTINGS+=" -re 550"
SETTINGS+=" -D 0.1"
SETTINGS+=" -xpos 0.7"
SETTINGS+=" -ypos 0.5"
SETTINGS+=" -uinfx 0.0"
SETTINGS+=" -uinfy 0.0"
SETTINGS+=" -mollfactor 4"

SETTINGS+=" -lcfl 0.01"
SETTINGS+=" -cfl 0.5"	
SETTINGS+=" -fc 0.5"

SETTINGS+=" -lr 0.1"
SETTINGS+=" -gamma 0.9"
SETTINGS+=" -greedyEps 0.1"
SETTINGS+=" -shared 1"

SETTINGS+=" -lmax 8"
SETTINGS+=" -jump 2"
SETTINGS+=" -rtol 1e-4"
SETTINGS+=" -ctol 1e-6"
SETTINGS+=" -uniform 0"
SETTINGS+=" -hilbert 0"

SETTINGS+=" -fmm velocity"
SETTINGS+=" -fmm-potential 0"
SETTINGS+=" -fmm-theta 0.8"
SETTINGS+=" -core-fmm sse"
SETTINGS+=" -refine-omega-only 0"
SETTINGS+=" -fmm-skip 0"

SETTINGS+=" -vtu 0"
SETTINGS+=" -particles 1"
SETTINGS+=" -rio free_frame"

SETTINGS+=" -ramp 50"
SETTINGS+=" -adaptfreq 5"
SETTINGS+=" -dumpfreq 30"
SETTINGS+=" -savefreq 1000"
SETTINGS+=" -tStartFTLE 0.0"
SETTINGS+=" -tEndFTLE 0.0"
SETTINGS+=" -tStartEff 0.0"
SETTINGS+=" -tEndEff 0.0"

SETTINGS+=" -obstacle heterogeneous"
SETTINGS+=" -factory factory"
		
RESTART=" -restart 0"

OPTIONS=${SETTINGS}${RESTART}

export OMP_NUM_THREADS=2
export LD_LIBRARY_PATH=/Users/mgazzola/mattia/cse/lib/TBB/tbb40_233oss/build/macos_intel64_gcc_cc4.2.1_os10.7.2_release/:/Developer/opt/intel/Compiler/11.1/084/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=/Users/mgazzola/mattia/cse/lib/TBB/tbb40_233oss/build/macos_intel64_gcc_cc4.2.1_os10.7.2_release/:/Developer/opt/intel/Compiler/11.1/084/lib/:$DYLD_LIBRARY_PATH
echo $OMP_NUM_THREADS
if [[ ${RESTART} == " -restart 0" ]]; then
rm -fr ../run
fi
mkdir ../run
cp ../makefiles/${EXECNAME} ../run/executable  
cp ctrl* factory ../run/
cd ../run
./executable ${OPTIONS}
#./executable ${OPTIONS} | tee lafiga
#gdb executable














