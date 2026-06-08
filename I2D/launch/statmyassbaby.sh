#!/bin/bash

BASEPATH=/Volumes/mattiaScratch/vehicles/FINAL/

ARRAY=()

ARRAY[0]=solo_4actions_100
ARRAY[1]=solo_4actions_250
ARRAY[2]=solo_4actions_500
ARRAY[3]=solo_4actions_750
ARRAY[4]=solo_4actions_1000
ARRAY[5]=solo_4actions_1250
ARRAY[6]=solo_4actions_1500
ARRAY[7]=solo_4actions_1750
ARRAY[8]=solo_4actions_2000

ARRAY[9]=solo_3actions_100
ARRAY[10]=solo_3actions_250
ARRAY[11]=solo_3actions_500
ARRAY[12]=solo_3actions_750
ARRAY[13]=solo_3actions_1000
ARRAY[14]=solo_3actions_1250
ARRAY[15]=solo_3actions_1500
ARRAY[16]=solo_3actions_1750
ARRAY[17]=solo_3actions_2000

ARRAY[18]=food_3actions_100
ARRAY[19]=food_3actions_250
ARRAY[20]=food_3actions_500
ARRAY[21]=food_3actions_750
ARRAY[22]=food_3actions_1000
ARRAY[23]=food_3actions_1250
ARRAY[24]=food_3actions_1500
ARRAY[25]=food_3actions_1750
ARRAY[26]=food_3actions_2000

ARRAY[27]=food_4actions_100
ARRAY[28]=food_4actions_250
ARRAY[29]=food_4actions_500
ARRAY[30]=food_4actions_750
ARRAY[31]=food_4actions_1000
ARRAY[32]=food_4actions_1250
ARRAY[33]=food_4actions_1500
ARRAY[34]=food_4actions_1750
ARRAY[35]=food_4actions_2000

ARRAY[36]=food_3actions_100_noavoid
ARRAY[37]=food_3actions_250_noavoid
ARRAY[38]=food_3actions_500_noavoid
ARRAY[39]=food_3actions_750_noavoid
ARRAY[40]=food_3actions_1000_noavoid
ARRAY[41]=food_3actions_1250_noavoid
ARRAY[42]=food_3actions_1500_noavoid
ARRAY[43]=food_3actions_1750_noavoid
ARRAY[44]=food_3actions_2000_noavoid

ARRAY[45]=food_4actions_100_noavoid
ARRAY[46]=food_4actions_250_noavoid
ARRAY[47]=food_4actions_500_noavoid
ARRAY[48]=food_4actions_750_noavoid
ARRAY[49]=food_4actions_1000_noavoid
ARRAY[50]=food_4actions_1250_noavoid
ARRAY[51]=food_4actions_1500_noavoid
ARRAY[52]=food_4actions_1750_noavoid
ARRAY[53]=food_4actions_2000_noavoid

for (( c=1; c<=8; c++ ))
do

NAME=${ARRAY[${c}]}
mkdir -p anal/${NAME}
cp ../matlab/clustermyass.m ../matlab/clusterEdi.m ../matlab/fromMatlabToPlot.m anal/${NAME}
screen -d -m -S ${NAME}
screen -X -S ${NAME} screen -t test 1
screen -X -S ${NAME} -p 1 stuff "cd anal/${NAME}; matlab -nodisplay -nodesktop -nojvm -nosplash -r \"clustermyass('${BASEPATH}/${NAME}')\" > lafiga
" 

done





