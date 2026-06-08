PROGNAME=$1
EXECNAME=$2
NPROCESSES=$3
TIMES=$4
WCLOCK=$5

NPROCESSORS=$(($NPROCESSES*16))
MAINDIR=${EXECNAME}

EXEC_LOCATION="/cluster/home/infk/mgazzola/mattia/MRAGapps/IF2D_ROCKS/makefiles/"
BASEPATH="/cluster/work/infk/mgazzola/fallingDisk/"
		
SETTINGS=
SETTINGS+=" -study FLOW_PAST_FIXED_OBSTACLE"
SETTINGS+=" -bpd 16"
SETTINGS+=" -nthreads 16"

SETTINGS+=" -lambda 1e4"
SETTINGS+=" -tend 10"
SETTINGS+=" -re 9500"
SETTINGS+=" -D 0.1"
SETTINGS+=" -xpos 0.25"
SETTINGS+=" -ypos 0.5"
SETTINGS+=" -mollfactor 2"

SETTINGS+=" -lcfl 0.1"
SETTINGS+=" -cfl 0.75"	

SETTINGS+=" -lmax 9"
SETTINGS+=" -jump 2"
SETTINGS+=" -rtol 1e-4"
SETTINGS+=" -ctol 1e-6"
SETTINGS+=" -uniform 0"
SETTINGS+=" -hilbert 0"

SETTINGS+=" -fmm velocity"
SETTINGS+=" -fmm-theta 0.8"
SETTINGS+=" -core-fmm sse"
SETTINGS+=" -refine-omega-only 0"
SETTINGS+=" -fmm-skip 0"

SETTINGS+=" -vtu 0"
SETTINGS+=" -particles 1"
SETTINGS+=" -rio free_outlet_only"

SETTINGS+=" -ramp 100"
SETTINGS+=" -adaptfreq 5"
SETTINGS+=" -dumpfreq 2"
SETTINGS+=" -savefreq 5000"

SETTINGS+=" -obstacle NACA"

SETTINGS+=" -NACA_angle 354"
SETTINGS+=" -NACA_plain 1"		
SETTINGS+=" -NACA_nTrailingBleeds 0"

if [[ "$(hostname)" == brutus* ]]; then
	export LD_LIBRARY_PATH+=:/cluster/work/infk/diegor/myVTK/lib/vtk-5.2/:/cluster/work/infk/diegor/tbb22_012oss/build/linux_intel64_icc_cc4.1.2_libc2.5_kernel2.6.18_release/:
	export OMP_NUM_THREADS=16

	for (( a=0; a<=8; a++ ))
	do		
		NBLEEDS=" -NACA_nRearBleeds "$a
		EXECNAME=${MAINDIR}"_rear_bleeds."$a

		mkdir -p ${BASEPATH}${MAINDIR}/${EXECNAME}
		cp ${EXEC_LOCATION}${PROGNAME} ${BASEPATH}${MAINDIR}/${EXECNAME}/
		cd ${BASEPATH}${MAINDIR}/${EXECNAME}

		RESTART=" -restart 0"
		OPTIONS=${SETTINGS}${RESTART}${NBLEEDS}

		echo "Submission 0..."
		
		if [[ $NPROCESSES == 1 ]]; then
			bsub -J ${EXECNAME} -R 'select[model==Opteron8384]' -n 16 -W ${WCLOCK} -o out time numactl --interleave=all ./${PROGNAME} ${OPTIONS}
			echo "bsub -J ${EXECNAME} -R 'select[model==Opteron8384]' -n 16 -W ${WCLOCK} -o out time numactl --interleave=all ./${PROGNAME} ${OPTIONS}" > submissionLog
		fi
		if [[ $NPROCESSES > 1 ]]; then
			bsub -J ${EXECNAME} -n $NPROCESSORS -R 'select[model==Opteron8384]' -R 'span[ptile=16]' -W ${WCLOCK} -o out mpirun -np $NPROCESSES -pernode time numactl --interleave=all ./${PROGNAME} ${OPTIONS}
			echo "bsub -J ${EXECNAME} -n $NPROCESSORS -R 'select[model==Opteron8384]' -R 'span[ptile=16]' -W ${WCLOCK} -o out mpirun -np $NPROCESSES -pernode time numactl --interleave=all ./${PROGNAME} ${OPTIONS}" > submissionLog
		fi


		RESTART=" -restart 1"
		OPTIONS=${SETTINGS}${RESTART}${NBLEEDS}
		for (( c=1; c<=${TIMES}-1; c++ ))
		do
			echo "Submission $c..."		

			if [[ $NPROCESSES == 1 ]]; then
				bsub -J ${EXECNAME} -w "ended(${EXECNAME})" -R 'select[model==Opteron8384]' -n 16 -W ${WCLOCK} -o out time numactl --interleave=all ./${PROGNAME} ${OPTIONS}
				echo "bsub -J ${EXECNAME} -w \"ended(${EXECNAME})\" -R 'select[model==Opteron8384]' -n 16 -W ${WCLOCK} -o out time numactl --interleave=all ./${PROGNAME} ${OPTIONS}" > submissionLog
			fi
			if [[ $NPROCESSES > 1 ]]; then
				bsub -J ${EXECNAME} -w "ended(${EXECNAME})" -n $NPROCESSORS -R 'select[model==Opteron8384]' -R 'span[ptile=16]' -W ${WCLOCK} -o out mpirun -np $NPROCESSES -pernode time numactl --interleave=all ./${PROGNAME} ${OPTIONS}
				echo "bsub -J ${EXECNAME} -w \"ended(${EXECNAME})\" -n $NPROCESSORS -R 'select[model==Opteron8384]' -R 'span[ptile=16]' -W ${WCLOCK} -o out mpirun -np $NPROCESSES -pernode time numactl --interleave=all ./${PROGNAME} ${OPTIONS}" > submissionLog
			fi
		done
	done
fi



