PROGNAME=$1
EXECNAME=$2
TIMES=$3
WCLOCK=$4
NTHREADS=$5

NPROCESSES=1
NPROCESSORS=$(($NPROCESSES*$NTHREADS))

BASEPATH="/cluster/scratch_xl/public/mgazzola/"

SETTINGS=
SETTINGS+=" -study FTLE"
SETTINGS+=" -nthreads 48"

SETTINGS+=" -path /cluster/scratch_xl/public/mgazzola/wimMostEfficient3D_HD/"
SETTINGS+=" -uinfx 0.0"
SETTINGS+=" -uinfy 0.0"

SETTINGS+=" -TFTLE 2.0"
SETTINGS+=" -tStartFTLE 6.0"
SETTINGS+=" -tEndFTLE 10.0"
SETTINGS+=" -dtFTLE 0.1"

SETTINGS+=" -lmax 10"
SETTINGS+=" -jump 2"
SETTINGS+=" -rtol 1e-6"
SETTINGS+=" -uniform 0"

RESTART=" -restart 0"

OPTIONS=${SETTINGS}${RESTART}

if [[ "$(hostname)" == brutus* ]]; then
mkdir -p ${BASEPATH}${EXECNAME}
cp ftleBrutusLaunch.sh ${BASEPATH}${EXECNAME}/
cp ctrl.* factory ${BASEPATH}${EXECNAME}/
cp ../makefiles/${PROGNAME} ${BASEPATH}${EXECNAME}/
cd ${BASEPATH}${EXECNAME}
export LD_LIBRARY_PATH+=:/cluster/home/infk/mgazzola/mattia/LIB/myVTK/lib/vtk-5.2/:/cluster/home/infk/mgazzola/mattia/LIB/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:
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

cd -