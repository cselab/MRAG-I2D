#!/bin/bash

EXECNAME=$1

SETTINGS=
SETTINGS+=" -study LEARNING_DIPOLE"
SETTINGS+=" -nthreads 2"
SETTINGS+=" -D 0.005"
SETTINGS+=" -xpos 0.5"
SETTINGS+=" -ypos 0.0001"

#SETTINGS+=" -lr 0.05"
SETTINGS+=" -gamma 0.9"
#SETTINGS+=" -greedyEps 0.005"
SETTINGS+=" -shared 1"

SETTINGS+=" -smooth 0"
SETTINGS+=" -isControlled 1"
SETTINGS+=" -exploitation true"
SETTINGS+=" -learningTime 101"
SETTINGS+=" -nLearningLevels 1" # number of levels to progressively halve the lr, e.g. lr = 0.01 for first learning interval, lr = 0.005 for second, lr = 0.0025 for third, ...
SETTINGS+=" -fitnessTime 100"
SETTINGS+=" -fitnessSaveFreq 100"
SETTINGS+=" -fitnessBuffer 1.1" # take (FINTESSBUFFER - 1)% longer to get rid of transient after a refresh
SETTINGS+=" -navg 5" # number of runs to average the fitness over after learning time
SETTINGS+=" -individualFitness false" # compute fitness for each dipole
		
SETTINGS+=" -savefreq 10000"
SETTINGS+=" -factory factoryPF"

RESTART=" -restart 0"

SETTINGS2=
SETTINGS2+=" -lr 0.0"
SETTINGS2+=" -greedyEps 0.01"

OPTIONS=${SETTINGS}${RESTART}${SETTINGS2}

export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/Users/menahel/Alexia/tbb40_20120408oss/build/macos_intel64_gcc_cc4.2.1_os10.6.6_release/:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=/Users/menahel/Alexia/tbb40_20120408oss/build/macos_intel64_gcc_cc4.2.1_os10.6.6_release/:$DYLD_LIBRARY_PATH
echo $OMP_NUM_THREADS
#rm -fr ../run
#mkdir ../run
cp ../makefiles/${EXECNAME} ../run/executable  
cp ctrl* factoryPF ../run/
cd ../run
./executable ${OPTIONS}

if [ "${RESTART}" = " -restart 1" ]; then
SETTINGS2=
SETTINGS2+=" -lr 0.0"
SETTINGS2+=" -greedyEps 0.01"
OPTIONS=${SETTINGS}${RESTART}${SETTINGS2}
./executable ${OPTIONS}

SETTINGS2=
SETTINGS2+=" -lr 0.0"
SETTINGS2+=" -greedyEps 0.0"
OPTIONS=${SETTINGS}${RESTART}${SETTINGS2}
./executable ${OPTIONS}
fi
#echo "Ok, bye! "
#./executable ${OPTIONS} | tee lafiga
#gdb executable
