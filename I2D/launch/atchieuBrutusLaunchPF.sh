PROGNAME=$1
EXECNAME=$2
TIMES=$3
WCLOCK=$4
NTHREADS=$5

NPROCESSES=1
NPROCESSORS=$(($NPROCESSES*$NTHREADS))

BASEPATH="/cluster/scratch_xl/public/atchieu/"

echo "** Defining program settings"

SETTINGS=
SETTINGS+=" -study LEARNING_DIPOLE"
SETTINGS+=" -nthreads ${NTHREADS}"
SETTINGS+=" -D 0.0025"
SETTINGS+=" -xpos 0.5"
SETTINGS+=" -ypos 0.0"

SETTINGS+=" -lr 0.01" 
SETTINGS+=" -gamma 0.95"
SETTINGS+=" -greedyEps 0.01"
SETTINGS+=" -shared 1"

SETTINGS+=" -smooth 0"
SETTINGS+=" -isControlled 1"
SETTINGS+=" -exploitation true"
SETTINGS+=" -individualFitness true" # compute fitness for each dipole
SETTINGS+=" -learningTime 10000"
SETTINGS+=" -nLearningLevels 2" # number of levels to progressively halve the lr, e.g. lr = 0.01 for first learning interval, lr = 0.005 for second, lr = 0.0025 for third, ...
SETTINGS+=" -fitnessTime 200"
SETTINGS+=" -fitnessSaveFreq 1000000"
SETTINGS+=" -fitnessBuffer 1.5" # take (FINTESSBUFFER - 1)% longer to get rid of transient after a refresh
SETTINGS+=" -navg 100"

SETTINGS+=" -savefreq 2000000"
SETTINGS+=" -learnDump 10000"
SETTINGS+=" -factory factoryPF"

RESTART=" -restart 0"	# true = restart from last 

OPTIONS=${SETTINGS}${RESTART}

if [[ "$(hostname)" == brutus* ]]; then

echo "** Making directory structure in ${BASEPATH}${EXECNAME}"
mkdir -p ${BASEPATH}${EXECNAME}
cp atchieuBrutusLaunchPF.sh ${BASEPATH}${EXECNAME}/
cp factoryPF ${BASEPATH}${EXECNAME}/
cp ../makefiles/${PROGNAME} ${BASEPATH}${EXECNAME}/
cd ${BASEPATH}${EXECNAME}

echo "** Exporting shell variables"
export LD_LIBRARY_PATH+=:/cluster/home/mavt/atchieu/Library/myVTK/lib/vtk-5.2:
export LD_LIBRARY_PATH+=:/cluster/home/mavt/atchieu/Library/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release:
export OMP_NUM_THREADS=${NTHREADS}

echo "** Submission 0..."
bsub -J ${EXECNAME} -n ${NTHREADS} -R span[ptile=${NTHREADS}] -W ${WCLOCK} -o output.txt time ./${PROGNAME} ${OPTIONS}

RESTART=" -restart 1"
OPTIONS=${SETTINGS}${RESTART}
for (( c=1; c<=${TIMES}-1; c++ ))
do
echo "** Submission $c..."		
	bsub -J ${EXECNAME} -n ${NTHREADS} -R span[ptile=${NTHREADS}] -w "ended(${EXECNAME})" -W ${WCLOCK} -o out time ./${PROGNAME} ${OPTIONS}
done

fi


