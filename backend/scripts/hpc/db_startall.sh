#$ -S /bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -m beas

# ---------

dividejob="qsub -N db_divide db_senddivide.sh -l h_rt=03:00:00 -l h_vmem=10G -cwd"

$dividejob

# --- when done... ---

sendbatchjob="qsub -N db_sendbatch db_sendbatch.sh -hold_jid db_divide -cwd"

$sendbatchjob

# ------

collectjob="qsub -N db_conquer -hold_jid db_and db_sendconquer.sh -cwd"

$collectjob

# --- remove residuals ---

# rm -r ${dbdir}/${dbname}_csv
# rm -r ${dbdir}/${dbname}_csv_ext

# ------

echo 'Submitted all jobs! ;)"