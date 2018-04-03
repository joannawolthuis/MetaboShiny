#!/bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -l h_rt=12:00:00
#$ -l h_vmem=30G

/home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save /home/cog/jwolthuis/testscripts/peakcalling/xcms_end.R $1

#/home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save xcms_call_peaks.R 1