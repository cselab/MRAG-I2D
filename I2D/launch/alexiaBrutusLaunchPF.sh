#!/bin/bash

EXECNAME=$1
PROGNAME=$2
TIMES=$3
WCLOCK=$4
NTHREADS=$5

BASEPATH=/cluster/work/scr3/alexiad/

SETTINGS=
SETTINGS+=" -study LEARNING_DIPOLE"
SETTINGS+=" -nthreads ${NTHREADS}"
SETTINGS+=" -D 0.005"
SETTINGS+=" -xpos 0.5"
SETTINGS+=" -ypos 0.1"

#SETTINGS+=" -lr 0.05"
SETTINGS+=" -gamma 0.9"
#SETTINGS+=" -greedyEps 0.01"
SETTINGS+=" -shared 1"

SETTINGS+=" -smooth 0"
SETTINGS+=" -isControlled 1"
SETTINGS+=" -exploitation 1"
SETTINGS+=" -learningTime 1000"
SETTINGS+=" -fitnessTime 100"
SETTINGS+=" -fitnessSaveFreq 1000"
SETTINGS+=" -navg 3" # number of runs to average the fitness over after learning time
		
SETTINGS+=" -savefreq 1000"
SETTINGS+=" -factory factoryPF"

RESTART=" -restart 0"

SETTINGS2=
SETTINGS2+=" -lr 0.01"
SETTINGS2+=" -greedyEps 0.01"

OPTIONS=${SETTINGS}${RESTART}${SETTINGS2}

if [[ "$(hostname)" == brutus* ]]; then
mkdir -p ${BASEPATH}${PROGNAME}
cp alexiaBrutusLaunchPF.sh ${BASEPATH}${PROGNAME}/
cp ctrl.* factoryPF ${BASEPATH}${PROGNAME}/
cp ../makefiles/${EXECNAME} ${BASEPATH}${PROGNAME}/
cd ${BASEPATH}${PROGNAME}
export LD_LIBRARY_PATH+=:/cluster/home/mavt/alexiad/LIB/tbb30_20100406oss/build/linux_intel64_gcc_cc4.1.2_libc2.5_kernel2.6.18_release/:/cluster/home/mavt/alexiad/LIB/myVTK/lib/vtk-5.2/:
export OMP_NUM_THREADS=${NTHREADS}

echo "Submission 0..."
bsub -J ${PROGNAME} -n ${NTHREADS} -R span[ptile=${NTHREADS}] -W ${WCLOCK} -o out time ./${EXECNAME} ${OPTIONS}

for (( c=1; c<=${TIMES}-1; c++ ))
do
#if [ "c" = "1" ] && [ "${RESTART}" = " -restart 1" ]; then
#SETTINGS2=
#SETTINGS2+=" -lr 0.0"
#SETTINGS2+=" -greedyEps 0.01"
#OPTIONS=${SETTINGS}${RESTART}${SETTINGS2}
#fi
#if [ "c" = "2" ] && [ "${RESTART}" = " -restart 1" ]; then
#SETTINGS2=
#SETTINGS2+=" -lr 0.0"
#SETTINGS2+=" -greedyEps 0.0"
#OPTIONS=${SETTINGS}${RESTART}${SETTINGS2}
#fi
echo "Submission $c..."		
	bsub -J ${PROGNAME} -n ${NTHREADS} -R span[ptile=${NTHREADS}] -w "ended(${PROGNAME})" -W ${WCLOCK} -o out time ./${EXECNAME} ${OPTIONS}
done
fi

