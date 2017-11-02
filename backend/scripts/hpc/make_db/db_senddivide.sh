#!/bin/bash
#$ -N db_divide
#$ -l h_rt=03:00:00
#$ -l h_vmem=40G
#$ -M j.c.wolthuis-2@umcutrecht.nl
#$ -m beas

/home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save /hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/scripts/hpc/db_divide.R