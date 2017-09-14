#!/bin/bash
#$ -N db_conquer
#$ -l h_rt=12:00:00
#$ -l h_vmem=40G
#$ -M j.c.wolthuis-2@umcutrecht.nl
#$ -m beas

TMPDIR=/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/tmp /home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save /hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/scripts/hpc/db_conquer.R

