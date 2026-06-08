#!/bin/bash

PROGNAME=$1
EXECNAME=$2
TIMES=$3
WCLOCK=$4
NTHREADS=$5

BASEPATH=/cluster/work/infk/mgazzola/fallingDisk/STUDY/vehicles/

SETTINGS=
SETTINGS+=" -study LEARNING_DIPOLE"
SETTINGS+=" -nthreads ${NTHREADS}"
SETTINGS+=" -D 0.005"
SETTINGS+=" -xpos 0.5"
SETTINGS+=" -ypos 0.5"

SETTINGS+=" -lr 0.01"
SETTINGS+=" -gamma 0.9"
SETTINGS+=" -greedyEps 0.05"
SETTINGS+=" -shared 1"

SETTINGS+=" -smooth 1"
SETTINGS+=" -isControlled 1"
		
SETTINGS+=" -savefreq 100"
SETTINGS+=" -factory factoryPF"

RESTART=" -restart 0"

OPTIONS=${SETTINGS}${RESTART}

if [[ "$(hostname)" == brutus* ]]; then
mkdir -p ${BASEPATH}${EXECNAME}
cp brutusLaunchPF.sh ${BASEPATH}${EXECNAME}/
cp ctrl.* factoryPF ${BASEPATH}${EXECNAME}/
cp ../makefiles/${PROGNAME} ${BASEPATH}${EXECNAME}/
cd ${BASEPATH}${EXECNAME}
export LD_LIBRARY_PATH+=:/cluster/home/infk/mgazzola/mattia/LIB/myVTK/lib/vtk-5.2/:/cluster/home/infk/mgazzola/mattia/LIB/tbb30_20100406oss/build/linux_intel64_gcc_cc4.1.2_libc2.5_kernel2.6.18_release/:
export OMP_NUM_THREADS=${NTHREADS}

echo "Submission 0..."
bsub -J ${EXECNAME} -n ${NTHREADS} -R span[ptile=${NTHREADS}] -W ${WCLOCK} -o out time ./${PROGNAME} ${OPTIONS}

for (( c=1; c<=${TIMES}-1; c++ ))
do
echo "Submission $c..."		
	bsub -J ${EXECNAME} -n ${NTHREADS} -R span[ptile=${NTHREADS}] -w "ended(${EXECNAME})" -W ${WCLOCK} -o out time ./${PROGNAME} ${OPTIONS}
done
fi