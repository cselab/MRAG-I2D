#!/bin/bash

EXECNAME=$1

SETTINGS=
SETTINGS+=" -study RL"
SETTINGS+=" -D 0.025"
SETTINGS+=" -xpos 0.7"
SETTINGS+=" -ypos 0.2"

SETTINGS+=" -lr 0.1"
SETTINGS+=" -gamma 0.9"
SETTINGS+=" -greedyEps 0.1"
SETTINGS+=" -shared 1"
		
SETTINGS+=" -savefreq 10000"
SETTINGS+=" -learnDump 10"
SETTINGS+=" -factory factoryRL"

RESTART=" -restart 0"

OPTIONS=${SETTINGS}${RESTART}

export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/Users/alexia/Documents/Software/Library/tbb40_20120201oss/build/macos_intel64_gcc_cc4.2.1_os10.6.8_release/:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=/Users/alexia/Documents/Software/Library/tbb40_20120201oss/build/macos_intel64_gcc_cc4.2.1_os10.6.8_release/:$DYLD_LIBRARY_PATH
echo $OMP_NUM_THREADS
rm -fr ../run
mkdir ../run
cp ../makefiles/${EXECNAME} ../run/executable  
cp ctrl* factoryRL ../run/
cd ../run
./executable ${OPTIONS}
#./executable ${OPTIONS} | tee lafiga
#gdb executable









