#$ -S /bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -m beas

nfiles="$(find ${dbdir}/${dbname}_csv_ext -type f -name '*.csv' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"
arrayjob="qsub -N db_and -t 1-${nfiles_nw}:1 db_sendand.sh -cwd"

$arrayjob
