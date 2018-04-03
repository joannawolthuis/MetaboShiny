#$ -S /bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl

nfiles="$(find /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/BrSp/MZXML -type f -name '*.raw' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"
arrayjob="qsub -N xcms_call -t 1-${nfiles_nw}:1 backend/scripts/hpc/run_raw.sh -cwd"

$arrayjob
