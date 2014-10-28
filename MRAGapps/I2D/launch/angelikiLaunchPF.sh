#!/bin/bash

EXECNAME=$1
NTHREADS=1

SETTINGS=
SETTINGS+=" -study LEARNING_DIPOLE"
SETTINGS+=" -nthreads ${NTHREADS}"
SETTINGS+=" -D 0.01"
SETTINGS+=" -xpos 0.5"
SETTINGS+=" -ypos 0.5"

SETTINGS+=" -lr 0.01"
SETTINGS+=" -gamma 0.98"
SETTINGS+=" -greedyEps 0.01"
SETTINGS+=" -shared 1"

SETTINGS+=" -smooth 0"
SETTINGS+=" -isControlled 1"
SETTINGS+=" -savefreq 100000" # policy save frequency (writes a ton of data if state space is large!!)

SETTINGS+=" -learningTime 100000" # amount of time specified to learn, total learning time is nLearningLevels*learningTime 
SETTINGS+=" -nLearningLevels 1" # number of levels to progressively halve the lr, e.g. lr = 0.01 for first learning interval, lr = 0.005 for second, lr = 0.0025 for third, ...
SETTINGS+=" -fitnessTime 10" # time to average fitness function over, CAREFULLY SELECT THIS!
SETTINGS+=" -navg 1" # number of runs to average the fitness over after learningTime*nLearningLevels
SETTINGS+=" -fitnessSaveFreq 100000" # change to very large if dumping less, then it only saves at the end of each averaging stage
SETTINGS+=" -fitnessBuffer 1.5" # take (FINTESSBUFFER - 1)% longer to get rid of transient after a refresh
SETTINGS+=" -individualFitness true" # compute fitness for each dipole
SETTINGS+=" -exploitation true" # if true, the learning rate is turned to zero after learning time (and possibly greedyEps, check code!)

SETTINGS+=" -learnDump 100000"
SETTINGS+=" -factory factoryPF"

RESTART=" -restart 0"

OPTIONS=${SETTINGS}${RESTART}

export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/Users/laskariangeliki/Documents/tbb40_297oss/lib/:$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=/Users/laskariangeliki/Documents/tbb40_297oss/lib/:$LD_LIBRARY_PATH

if [[ ${RESTART} == " -restart 0" ]]; then
rm -fr ../run
mkdir ../run
fi

if [[ ${RESTART} == " -restart 1" ]]; then
rm ../run/*.tga
rm ../run/*.txt
fi

cp marcelLaunchPF.sh ../run/
cp ../makefiles/${EXECNAME} ../run/executable  
cp factoryPF ../run/
# cp ctrl* ../run/
cd ../run

./executable ${OPTIONS}










