PROGNAME=$1
EXECNAME=$2
NPROCESSES=$3
TIMES=$4
WCLOCK=$5

NPROCESSORS=$(($NPROCESSES*16))

BASEPATH="/cluster/work/infk/mimeauc/STUDY/"

SETTINGS=
SETTINGS+=" -bpd 8"
SETTINGS+=" -nthreads 16"

SETTINGS+=" -lambda 1e4"
SETTINGS+=" -lambdadt 10"
SETTINGS+=" -tend 200"
SETTINGS+=" -re 1000"
SETTINGS+=" -D 0.03"
SETTINGS+=" -xpos 0.25"
SETTINGS+=" -ypos 0.5"
SETTINGS+=" -mollfactor 2"

SETTINGS+=" -lcfl 0.1"
SETTINGS+=" -cfl 0.25"	

SETTINGS+=" -lmax 5"
SETTINGS+=" -jump 2"
SETTINGS+=" -rtol 1e-5"
SETTINGS+=" -ctol 1e-7"
SETTINGS+=" -uniform 0"

SETTINGS+=" -fmm velocity"
SETTINGS+=" -fmm-theta 0.8"
SETTINGS+=" -core-fmm sse"
SETTINGS+=" -fmm-skip 0"
SETTINGS+=" -refine-omega-only 0"

SETTINGS+=" -vtu 0"
SETTINGS+=" -particles 1"
SETTINGS+=" -rio 2"
SETTINGS+=" -killvort 5"

SETTINGS+=" -ramp 100"
SETTINGS+=" -adaptfreq 5"
SETTINGS+=" -dumpfreq 2"
SETTINGS+=" -savefreq 1000"	

SETTINGS+=" -obstacle links"

SETTINGS+=" -ctrl ctrl.413"

SETTINGS+=" -elly_angle 45"
SETTINGS+=" -elly_aspectRatio 0.5"
SETTINGS+=" -NACA_angle 0"
SETTINGS+=" -NACA_plain 1"
SETTINGS+=" -NACA_nRearBleeds 0"
SETTINGS+=" -NACA_nTrailingBleeds 1"

RESTART=" -restart 0"

OPTIONS=${SETTINGS}${RESTART}

if [[ "$(hostname)" == brutus* ]]; then
	mkdir -p ${BASEPATH}${EXECNAME}
	cp ${PROGNAME} ctrl.413 ${BASEPATH}${EXECNAME}/
	cd ${BASEPATH}${EXECNAME}

	export LD_LIBRARY_PATH+=:/cluster/work/infk/diegor/myVTK/lib/vtk-5.2/:/cluster/work/infk/diegor/tbb22_012oss/build/linux_intel64_icc_cc4.1.2_libc2.5_kernel2.6.18_release/:
	export OMP_NUM_THREADS=16
	echo "Submission 0..."
	
	if [[ $NPROCESSES == 1 ]]; then
	    bsub -J ${EXECNAME} -R 'select[model==Opteron8384]' -n 16 -W ${WCLOCK} -o out time numactl --interleave=all ./${PROGNAME} ${OPTIONS}
		echo "bsub -J ${EXECNAME} -R 'select[model==Opteron8384]' -n 16 -W ${WCLOCK} -o out time numactl --interleave=all ./${PROGNAME} ${OPTIONS}" > submissionLog
	fi
	if [[ $NPROCESSES > 1 ]]; then
	    bsub -J ${EXECNAME} -n $NPROCESSORS -R 'span[ptile=16]' -W ${WCLOCK} -o out mpirun -np $NPROCESSES -pernode time numactl --interleave=all ./${PROGNAME} ${OPTIONS}
		echo "bsub -J ${EXECNAME} -n $NPROCESSORS -R 'span[ptile=16]' -W ${WCLOCK} -o out mpirun -np $NPROCESSES -pernode time numactl --interleave=all ./${PROGNAME} ${OPTIONS}" > submissionLog
	fi
	
	
	RESTART=" -restart 1"
	OPTIONS=${SETTINGS}${RESTART}
	for (( c=1; c<=${TIMES}-1; c++ ))
	do
	  echo "Submission $c..."		
	  
	  if [[ $NPROCESSES == 1 ]]; then
	      bsub -J ${EXECNAME} -R 'select[model==Opteron8384]' -w "ended(${EXECNAME})" -n 16 -W ${WCLOCK} -o out time numactl --interleave=all ./${PROGNAME} ${OPTIONS}
		  echo "bsub -J ${EXECNAME} -R 'select[model==Opteron8384] -w \"ended(${EXECNAME})\" -n 16 -W ${WCLOCK} -o out time numactl --interleave=all ./${PROGNAME} ${OPTIONS}" > submissionLog
	  fi
	  if [[ $NPROCESSES > 1 ]]; then
	      bsub -J ${EXECNAME} -w "ended(${EXECNAME})" -n $NPROCESSORS -R 'span[ptile=16]' -W ${WCLOCK} -o out mpirun -np $NPROCESSES -pernode time numactl --interleave=all ./${PROGNAME} ${OPTIONS}
		  echo "bsub -J ${EXECNAME} -w \"ended(${EXECNAME})\" -n $NPROCESSORS -R 'span[ptile=16]' -W ${WCLOCK} -o out mpirun -np $NPROCESSES -pernode time numactl --interleave=all ./${PROGNAME} ${OPTIONS}" > submissionLog
	  fi
	done
fi
