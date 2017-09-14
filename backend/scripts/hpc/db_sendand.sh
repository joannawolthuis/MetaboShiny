#!/bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -l h_rt=00:30:00

/home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save /hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/scripts/hpc/db_and.R $SGE_TASK_ID