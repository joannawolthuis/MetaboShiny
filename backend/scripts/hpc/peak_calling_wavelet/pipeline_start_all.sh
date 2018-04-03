#$ -S /bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -m beas

# ----------------

HPCBASE="/hpc/cog_bioinf/ridder/users/jwolthuis"
PROJFOLDER="RES-2017-04-26_DBS_Turkey_DSM"
SCRIPTS="/home/cog/jwolthuis/testscripts/peakcalling2/"

PROJFOLDER="BrazilFarm1and2"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER

nfiles="$(find ${INFOLDER}/MZXML -type f -name '*.raw' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"
arrayjob="qsub -N raw -t 1-${nfiles_nw}:1 backend/scripts/hpc/run_raw.sh -cwd"
$arrayjob

# average
arrayjob="qsub -hold_jid raw -N average ${SCRIPTS}/hpc/run_average.sh -cwd"

# peak finding
nfiles="$(find {INFOLDER}/MZXML/results/averaged -type f -name '*.RData' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"
arrayjob="qsub -hold_jid average -N peaks -t 1-${nfiles_nw}:1 ${SCRIPTS}/hpc/run_peaks.sh -cwd"
$arrayjob
