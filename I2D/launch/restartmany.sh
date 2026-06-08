#!/bin/bash

CODEPATH=/cluster/home/infk/mgazzola/mattia/MRAGapps/IF2D_ROCKS/
REPOSITORY=/cluster/work/infk/mgazzola/fallingDisk/STUDY/vehicles/factories/cutoff_5/
BASEPATH=/cluster/work/infk/mgazzola/fallingDisk/STUDY/vehicles_parallel/


ARRAY=()

ARRAY[0]=solo_4actions_100
ARRAY[1]=solo_4actions_250
ARRAY[2]=solo_4actions_500
ARRAY[3]=solo_4actions_750
ARRAY[4]=solo_4actions_1000
ARRAY[5]=solo_4actions_1250
ARRAY[6]=solo_4actions_1500

ARRAY[7]=solo_3actions_100
ARRAY[8]=solo_3actions_250
ARRAY[9]=solo_3actions_500
ARRAY[10]=solo_3actions_750
ARRAY[11]=solo_3actions_1000
ARRAY[12]=solo_3actions_1250
ARRAY[13]=solo_3actions_1500

ARRAY[14]=food_3actions_100
ARRAY[15]=food_3actions_250
ARRAY[16]=food_3actions_500
ARRAY[17]=food_3actions_750
ARRAY[18]=food_3actions_1000
ARRAY[19]=food_3actions_1250
ARRAY[20]=food_3actions_1500

ARRAY[21]=food_4actions_100
ARRAY[22]=food_4actions_250
ARRAY[23]=food_4actions_500
ARRAY[24]=food_4actions_750
ARRAY[25]=food_4actions_1000
ARRAY[26]=food_4actions_1250
ARRAY[27]=food_4actions_1500

ARRAY[28]=food_3actions_100_noavoid
ARRAY[29]=food_3actions_250_noavoid
ARRAY[30]=food_3actions_500_noavoid
ARRAY[31]=food_3actions_750_noavoid
ARRAY[32]=food_3actions_1000_noavoid
ARRAY[33]=food_3actions_1250_noavoid
ARRAY[34]=food_3actions_1500_noavoid

ARRAY[35]=food_4actions_100_noavoid
ARRAY[36]=food_4actions_250_noavoid
ARRAY[37]=food_4actions_500_noavoid
ARRAY[38]=food_4actions_750_noavoid
ARRAY[39]=food_4actions_1000_noavoid
ARRAY[40]=food_4actions_1250_noavoid
ARRAY[41]=food_4actions_1500_noavoid

for (( c=14; c<=41; c++ ))
do
	NAME=${ARRAY[${c}]}
	
	NCORES=16
	
	#case "${c}" in
	#8) 	echo  "48 cores" NCORES=48 ;;
	#9) 	echo  "48 cores" NCORES=48 ;;
	#10) echo  "48 cores" NCORES=48 ;;
	#11) echo  "48 cores" NCORES=48 ;;
   	#12) echo  "48 cores" NCORES=48 ;;
   	#13) echo  "48 cores" NCORES=48 ;;
   	#*)  echo  "16 cores" ;;
	#esac
	
	mkdir -p ${BASEPATH}${NAME}
	cp ${REPOSITORY}${NAME}/factoryRL ${BASEPATH}${NAME}
	cp ${CODEPATH}/launch/brutusLaunchRL.sh ${BASEPATH}${NAME}
	cp ${CODEPATH}/makefiles/rl ${BASEPATH}${NAME}
	cd ${BASEPATH}${NAME}
	./brutusLaunchRL.sh rl ${NAME} 1 08:00 ${NCORES}	
done


