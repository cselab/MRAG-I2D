#!/bin/bash

PROGNAME=$1
EXECNAME=$2
TIMES=$3
WCLOCK=$4
NTHREADS=$5

NPROCESSES=1
NPROCESSORS=$(($NPROCESSES*$NTHREADS))

BASEPATH="/cluster/scratch_xl/public/hbabak/"

SETTINGS=
SETTINGS+=" -study FLUID_MEDIATED_INTERACTIONS"
SETTINGS+=" -bpd 32"
SETTINGS+=" -nthreads ${NTHREADS}"

SETTINGS+=" -lambda 1e4"
SETTINGS+=" -tend 6.0"
SETTINGS+=" -re 7142.857143"
SETTINGS+=" -D 0.05"
SETTINGS+=" -xpos 0.5"
SETTINGS+=" -ypos 0.5"
SETTINGS+=" -uinfx 0.0"
SETTINGS+=" -uinfy 0.0"
SETTINGS+=" -mollfactor 4"
SETTINGS+=" -sharp 1"

SETTINGS+=" -lcfl 0.01"
SETTINGS+=" -cfl 0.5"
SETTINGS+=" -fc 0.5"

SETTINGS+=" -lmax 8"
SETTINGS+=" -jump 2"
SETTINGS+=" -rtol 1e-3"
SETTINGS+=" -ctol 1e-4"
SETTINGS+=" -uniform 0"
SETTINGS+=" -hilbert 0"

SETTINGS+=" -fmm velocity"
SETTINGS+=" -fmm-potential 1"
SETTINGS+=" -fmm-theta 0.8"
SETTINGS+=" -core-fmm sse"
SETTINGS+=" -refine-omega-only 0"
SETTINGS+=" -fmm-skip 0"

SETTINGS+=" -vtu 1"
SETTINGS+=" -particles 1"
SETTINGS+=" -rio free_frame"

SETTINGS+=" -ramp 100"
SETTINGS+=" -adaptfreq 50"
SETTINGS+=" -dumpfreq 100"
SETTINGS+=" -savefreq 200"

SETTINGS+=" -obstacle heterogeneous"
SETTINGS+=" -factory factory"

# wim: put some parameters here, required but not used for non-learning cases
SETTINGS+=" -lr 0.3"
SETTINGS+=" -gamma 0.4"
SETTINGS+=" -greedyEps 0.3"
SETTINGS+=" -shared 0"

#SETTINGS+=" -ctrl ctrl.413"
#SETTINGS+=" -vxo 0.1"
#SETTINGS+=" -vyo 0.0"
#SETTINGS+=" -killvort 4"
#SETTINGS+=" -NACA_angle 0"
#SETTINGS+=" -NACA_nRearBleeds 0"
#SETTINGS+=" -NACA_nTrailingBleeds 0"
#SETTINGS+=" -NACA_plain 1"

RESTART=" -restart 0"

OPTIONS=${SETTINGS}${RESTART}

if [[ "$(hostname)" == brutus* ]]; then
mkdir -p ${BASEPATH}${EXECNAME}
cp ${0} ${BASEPATH}${EXECNAME}/
cp ctrl.* factory ${BASEPATH}${EXECNAME}/
cp ../makefiles/${PROGNAME} ${BASEPATH}${EXECNAME}/
cd ${BASEPATH}${EXECNAME}
export LD_LIBRARY_PATH+=:/cluster/home/infk/hbabak/tbb41_20130613oss/build/linux_intel64_gcc_cc4.7.0_libc2.12_kernel2.6.32_release/:/cluster/home/infk/hbabak/LIB/myVTK/lib/vtk-5.2/
export OMP_NUM_THREADS=${NTHREADS}

echo "Submission 0..."
bsub -J ${EXECNAME} -n ${NTHREADS} -R span[ptile=${NTHREADS}] -W ${WCLOCK} -o out time ./${PROGNAME} ${OPTIONS}

RESTART=" -restart 1"
OPTIONS=${SETTINGS}${RESTART}
for (( c=1; c<=${TIMES}-1; c++ ))
do
    echo "Submission $c..."
    bsub -J ${EXECNAME} -n ${NTHREADS} -R span[ptile=${NTHREADS}] -w "ended(${EXECNAME})" -W ${WCLOCK} -o out time ./${PROGNAME} ${OPTIONS}
done
fi