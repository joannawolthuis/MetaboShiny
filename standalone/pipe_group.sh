#!/bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -l h_rt=00:30:00

/home/cog/jwolthuis/R-3.3.3/bin/Rscript --no-save /hpc/cog_bioinf/ridder/users/jwolthuis/MassChecker/standalone/pipe_group.R $SGE_TASK_ID "$1" "$2" "$3"