#!/bin/bash
# -l h_rt=07:00:00
# -l h_vmem=40G
# -pe threaded 40
# -m beas
# -M j.c.wolthuis-2@umcutrecht.nl

/home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save /hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/scripts/db.make.hpc.R