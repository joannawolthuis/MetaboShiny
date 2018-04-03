#!/bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -l h_rt=00:45:00

/home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save /home/cog/jwolthuis/testscripts/peakcalling/xcms_call_peaks.R $1 $SGE_TASK_ID

#/home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save xcms_call_peaks.R 1