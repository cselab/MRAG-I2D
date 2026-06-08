#!/bin/bash

EXECNAME=$1

SETTINGS=
SETTINGS+=" -study CART"

SETTINGS+=" -lr 0.001"
SETTINGS+=" -gamma 0.99"
SETTINGS+=" -greedyEps 0.1"
		
SETTINGS+=" -savefreq 100000"

RESTART=" -restart 0"

OPTIONS=${SETTINGS}${RESTART}

export OMP_NUM_THREADS=2
echo $OMP_NUM_THREADS
//rm -fr ../run/
mkdir ../run
cp ../makefiles/${EXECNAME} ../run/executable  
cp ctrl* factory ../run/
cd ../run
./executable ${OPTIONS}









