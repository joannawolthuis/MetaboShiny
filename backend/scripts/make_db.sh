#!/bin/bash
#$ -N makedb_cont
#$ -l h_rt=72:00:00
#$ -l h_vmem=40G
#$ -pe threaded 39
#$ -m beas
#$ -M j.c.wolthuis-2@umcutrecht.nl
#$ -m beas

TMPDIR=/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/tmp /home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save /hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/scripts/db.make.hpc.R
