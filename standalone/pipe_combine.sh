#!/bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -l h_rt=06:00:00
#$ -l h_vmem=40G

/home/cog/jwolthuis/R-3.3.3/bin/Rscript --no-save /hpc/cog_bioinf/ridder/users/jwolthuis/MassChecker/standalone/pipe_combine.R $1