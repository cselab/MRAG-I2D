#!/bin/bash

PROGNAME=$1
TIMES=$2
WCLOCK=$3
NTHREADS=$4

BASEPATH=/cluster/work/scr3/dalmassg/
EXECNAME=rl
FACTORY=factoryRL_test1

SETTINGS=
SETTINGS+=" -study RL"
SETTINGS+=" -nthreads ${NTHREADS}"
SETTINGS+=" -D 0.025"
SETTINGS+=" -xpos 0.5"
SETTINGS+=" -ypos 0.5"

SETTINGS+=" -lr 0.1"
SETTINGS+=" -gamma 0.9"
SETTINGS+=" -greedyEps 0.1"
SETTINGS+=" -shared 1"
		
SETTINGS+=" -savefreq 100000"
SETTINGS+=" -learnDump 10"
SETTINGS+=" -factory ${FACTORY}"

RESTART=" -restart 0"

OPTIONS=${SETTINGS}${RESTART}

if [[ "$(hostname)" == brutus* ]]; then
mkdir -p ${BASEPATH}${PROGNAME}
cp giovanniBrutusLaunchRL.sh ${BASEPATH}${PROGNAME}/
cp ctrl.* ${FACTORY} ${BASEPATH}${PROGNAME}/
cp ../makefiles/${EXECNAME} ${BASEPATH}${PROGNAME}/executable
cd ${BASEPATH}${PROGNAME}
export LD_LIBRARY_PATH+=:/cluster/home/infk/mgazzola/mattia/LIB/myVTK/lib/vtk-5.2/:/cluster/home/infk/mgazzola/mattia/LIB/tbb22_012oss/build/linux_intel64_icc_cc4.1.2_libc2.5_kernel2.6.18_release/:
export OMP_NUM_THREADS=${NTHREADS}

echo "Submission 0..."
bsub -J ${PROGNAME} -n ${NTHREADS} -R select[model==Opteron6174] -W ${WCLOCK} -o out time ./executable ${OPTIONS}

RESTART=" -restart 1"
OPTIONS=${SETTINGS}${RESTART}
for (( c=1; c<=${TIMES}-1; c++ ))
do
echo "Submission $c..."		
	bsub -J ${PROGNAME} -n ${NTHREADS} -R select[model==Opteron6174] -w "ended(${PROGNAME})" -W ${WCLOCK} -o out time ./executable ${OPTIONS}
done
fi

