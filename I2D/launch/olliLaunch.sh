#!/bin/sh

EXECNAME=$1

SETTINGS=
SETTINGS+=" -study CYLINDER_INTERACTION"
SETTINGS+=" -bpd 2"
SETTINGS+=" -nthreads 2"

SETTINGS+=" -lambda 1e4"
SETTINGS+=" -tend 6"
SETTINGS+=" -re 100"
SETTINGS+=" -D 0.05"
SETTINGS+=" -xpos 0.5"
SETTINGS+=" -ypos 0.5"
SETTINGS+=" -mollfactor 2"

SETTINGS+=" -lcfl 0.1"
SETTINGS+=" -cfl 1.0"	

SETTINGS+=" -lmax 3"
SETTINGS+=" -jump 2"
SETTINGS+=" -rtol 1e-4"
SETTINGS+=" -ctol 1e-6"
SETTINGS+=" -uniform 0"
SETTINGS+=" -hilbert 0"

SETTINGS+=" -fmm velocity"
SETTINGS+=" -fmm-theta 0.8"
SETTINGS+=" -core-fmm sse"
SETTINGS+=" -refine-omega-only 1"
SETTINGS+=" -fmm-skip 0"

SETTINGS+=" -vtu 0"
SETTINGS+=" -particles 1"
SETTINGS+=" -rio free_frame"

SETTINGS+=" -ramp 100"
SETTINGS+=" -adaptfreq 1"
SETTINGS+=" -dumpfreq 30"
SETTINGS+=" -savefreq 1000"

SETTINGS+=" -obstacle cylTransport"
SETTINGS+=" -factory factory"

SETTINGS+=" -ctrl ctrl.413"
SETTINGS+=" -vxo 0.1"
SETTINGS+=" -vyo 0.0"
SETTINGS+=" -killvort 4"
SETTINGS+=" -NACA_angle 0"
SETTINGS+=" -NACA_nRearBleeds 0"
SETTINGS+=" -NACA_nTrailingBleeds 0"
SETTINGS+=" -NACA_plain 1"

RESTART=" -restart 0"

OPTIONS=${SETTINGS}${RESTART}

export OMP_NUM_THREADS=2
export DYLD_LIBRARY_PATH=/Users/mimeauc/tbb30_131oss/build/macos_ia32_gcc_cc4.2.1_os10.6.6_release/
echo $OMP_NUM_THREADS
mkdir ../run
//rm ../run/*
cp ${EXECNAME} ../run/exec  
cp ctrl* factory ../run/
cd ../run
./exec ${OPTIONS}


