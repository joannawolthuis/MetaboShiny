#$ -S /bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -m beas

# ----------------

HPCBASE="/hpc/cog_bioinf/ridder/users/jwolthuis"
PROJFOLDER="RES-2017-04-26_DBS_Turkey_DSM"
SCRIPTS="/home/cog/jwolthuis/testscripts/peakcalling/"

PROJFOLDER="RES-2017-04-26_DBS_Turkey_DSM"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
OUTFOLDER=${INFOLDER}"/peakcalled"
qsub -m beas -M J.C.Wolthuis-2@umcutrecht.nl -N xcms_end_tr xcms_runend.sh ${OUTFOLDER}
 
PROJFOLDER="RES-2017-05-01_DBS_France_DSM"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
OUTFOLDER=${INFOLDER}"/peakcalled"
qsub -m beas -M J.C.Wolthuis-2@umcutrecht.nl -N xcms_end_fr xcms_runend.sh ${OUTFOLDER}

PROJFOLDER="RES-2017-05-02_DBS_France-Spain_DSM"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
OUTFOLDER=${INFOLDER}"/peakcalled"
qsub -m beas -M J.C.Wolthuis-2@umcutrecht.nl -N xcms_end_fr_sp xcms_runend.sh ${OUTFOLDER}

PROJFOLDER="RES-2017-05-03_DBS_Spain-Gent_DSM"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
OUTFOLDER=${INFOLDER}"/peakcalled"
qsub -m beas -M J.C.Wolthuis-2@umcutrecht.nl -N xcms_end_sp_be xcms_runend.sh ${OUTFOLDER}

PROJFOLDER="Project2016_018_IBD"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
OUTFOLDER=${INFOLDER}"/peakcalled"
qsub -m beas -M J.C.Wolthuis-2@umcutrecht.nl -N xcms_end_ibd xcms_runend.sh ${OUTFOLDER}

# ----------------------

INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
OUTFOLDER=${INFOLDER}"/peakcalled"

mkdir -p ${OUTFOLDER}

cd $SCRIPTS

# --- send out ---

nfiles="$(find ${INFOLDER} -type f -name '*.mzXML' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"

qsub -N xcms_call -t 1-${nfiles_nw}:1 xcms_runcall.sh ${INFOLDER} -cwd

# ------

qsub -N xcms_end xcms_runend.sh ${OUTFOLDER}

# ------
echo 'Submitted all jobs! \(. w.)/"
