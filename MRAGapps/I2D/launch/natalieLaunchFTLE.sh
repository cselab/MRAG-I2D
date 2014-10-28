#!/bin/bash

EXECNAME=$1

SETTINGS=
SETTINGS+=" -study FTLE"
SETTINGS+=" -nthreads 2"

SETTINGS+=" -path /Users/mgazzola/mattia/cse/eclipseProjects/MRAGapps/IF2D_ROCKS/runFTLE/"
SETTINGS+=" -uinfx 0.0"
SETTINGS+=" -uinfy 0.0"

SETTINGS+=" -TFTLE 0.5"
SETTINGS+=" -tStartFTLE 1.0"
SETTINGS+=" -tEndFTLE 1.5"
SETTINGS+=" -dtFTLE 0.05"

SETTINGS+=" -lmax 6"
SETTINGS+=" -jump 2"
SETTINGS+=" -rtol 1e-5"
SETTINGS+=" -uniform 0"

RESTART=" -restart 0"
		
OPTIONS=${SETTINGS}${RESTART}

export OMP_NUM_THREADS=2
export LD_LIBRARY_PATH=/Users/mgazzola/mattia/cse/lib/TBB/tbb40_233oss/build/macos_intel64_gcc_cc4.2.1_os10.7.2_release/:/Developer/opt/intel/Compiler/11.1/084/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=/Users/mgazzola/mattia/cse/lib/TBB/tbb40_233oss/build/macos_intel64_gcc_cc4.2.1_os10.7.2_release/:/Developer/opt/intel/Compiler/11.1/084/lib/:$DYLD_LIBRARY_PATH
echo $OMP_NUM_THREADS
rm -fr ../run
mkdir ../run
cp ../makefiles/${EXECNAME} ../run/executable  
cp ctrl* factory ../run/
cd ../run
./executable ${OPTIONS}
#./executable ${OPTIONS} | tee lafiga
#gdb executable














