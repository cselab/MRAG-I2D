#!/bin/bash

#${CODEPATH}/matlab/clustermyass.m ${CODEPATH}/matlab/clusterEdi.m ${CODEPATH}/matlab/fromMatlabToPlot.m

bsub -J prova -n 1 -W 00:10 "matlab -nodisplay -nodesktop -nosplash -r \"clustermyass('./')\" > lafiga"