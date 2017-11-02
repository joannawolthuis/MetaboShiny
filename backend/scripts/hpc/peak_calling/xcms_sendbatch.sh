#$ -S /bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl

nfiles="$(find /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/RES-2017-04-26_DBS_Turkey_DSM/MZXML -type f -name '*.mzXML' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"
arrayjob="qsub -N xcms_call -t 1-${nfiles_nw}:1 backend/scripts/hpc/xcms_runcall.sh -cwd"

$arrayjob
