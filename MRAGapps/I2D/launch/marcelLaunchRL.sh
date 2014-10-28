#!/bin/bash

EXECNAME=$1

SETTINGS=
SETTINGS+=" -study RL"
SETTINGS+=" -nthreads 4"
SETTINGS+=" -D 0.005"
SETTINGS+=" -xpos 0.8"
SETTINGS+=" -ypos 0.5"

SETTINGS+=" -lr 0.001"
SETTINGS+=" -gamma 0.95"
SETTINGS+=" -greedyEps 0.01"
SETTINGS+=" -shared 1"
		
SETTINGS+=" -savefreq 100000"
SETTINGS+=" -learnDump 1000"
SETTINGS+=" -factory factoryRL"

OPTIONS=${SETTINGS}

export OMP_NUM_THREADS=4
export LD_LIBRARY_PATH=/Users/mgazzola/mattia/cse/lib/TBB/tbb40_233oss/build/macos_intel64_gcc_cc4.2.1_os10.7.2_release/:/Developer/opt/intel/Compiler/11.1/084/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=/Users/mgazzola/mattia/cse/lib/TBB/tbb40_233oss/build/macos_intel64_gcc_cc4.2.1_os10.7.2_release/:/Developer/opt/intel/Compiler/11.1/084/lib/:$DYLD_LIBRARY_PATH
echo $OMP_NUM_THREADS
rm -fr ../run
mkdir ../run
cp ../makefiles/${EXECNAME} ../run/executable  
cp ctrl* factoryRL ../run/ 
cd ../run
./executable ${OPTIONS}
#./executable ${OPTIONS} | tee lafiga
#gdb executable









