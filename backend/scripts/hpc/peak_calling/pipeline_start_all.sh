#$ -S /bin/bash
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -m beas

# ----------------

HPCBASE="/hpc/cog_bioinf/ridder/users/jwolthuis"
PROJFOLDER="RES-2017-04-26_DBS_Turkey_DSM"
SCRIPTS="/home/cog/jwolthuis/testscripts/peakcalling2/"

PROJFOLDER="BrazilFarm1and2"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER
OUTFOLDER=${INFOLDER}"/peakcalled"
# hold only works if not multiple xcms_calls are running, otherwise plz give them seperate names...
qsub -m beas -M J.C.Wolthuis-2@umcutrecht.nl -hold_jid xcms_call -N xcms_end_br12 $SCRIPTS/xcms_runend.sh ${OUTFOLDER}
 
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

PROJFOLDER="RES-2017-10-31_DSM_DBS_Brazil_1"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
OUTFOLDER=${INFOLDER}"/peakcalled"
qsub -m beas -M J.C.Wolthuis-2@umcutrecht.nl -N xcms_end_br1_8 xcms_runend.sh ${OUTFOLDER}

mkdir -p ${OUTFOLDER}

cd $SCRIPTS

# --- send out ---
HPCBASE="/hpc/cog_bioinf/ridder/users/jwolthuis"

PROJFOLDER="BrazilFarm1and2"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/$PROJFOLDER"
nfiles="$(find ${INFOLDER} -type f -name '*.mzXML' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"
qsub -N xcms_call_br12 -t 1-${nfiles_nw}:1 xcms_runcall.sh ${INFOLDER} -cwd

PROJFOLDER="RES-2017-05-01_DBS_France_DSM"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
nfiles="$(find ${INFOLDER} -type f -name '*.mzXML' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"
qsub -N xcms_call_france -t 1-${nfiles_nw}:1 xcms_runcall.sh ${INFOLDER} -cwdx

PROJFOLDER="RES-2017-04-26_DBS_Turkey_DSM"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
nfiles="$(find ${INFOLDER} -type f -name '*.mzXML' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"
qsub -N xcms_call_turkey -t 1-${nfiles_nw}:1 xcms_runcall.sh ${INFOLDER} -cwd

HPCBASE="/hpc/cog_bioinf/ridder/users/jwolthuis"
PROJFOLDER="RES-2017-04-26_DBS_Turkey_DSM"
SCRIPTS="/home/cog/jwolthuis/testscripts/peakcalling2/"

PROJFOLDER="RES-2017-10-31_DSM_DBS_Brazil_1"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
nfiles="$(find ${INFOLDER} -type f -name '*.mzXML' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"
#nfiles_nw="$(echo -e 5 | tr -d '[:space:]')"

qsub -N xcms_call_brfull -t 1-${nfiles_nw}:1 xcms_runcall.sh ${INFOLDER} -cwd

PROJFOLDER="RES-2017-05-02_DBS_France-Spain_DSM"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
nfiles="$(find ${INFOLDER} -type f -name '*.mzXML' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"
qsub -N xcms_call_fr_sp -t 1-${nfiles_nw}:1 xcms_runcall.sh ${INFOLDER} -cwd

PROJFOLDER="RES-2017-05-03_DBS_Spain-Gent_DSM"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
nfiles="$(find ${INFOLDER} -type f -name '*.mzXML' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"
qsub -N xcms_call_sp_be -t 1-${nfiles_nw}:1 xcms_runcall.sh ${INFOLDER} -cwd

PROJFOLDER="Project2016_018_IBD"
INFOLDER=$HPCBASE"/Data/Metabolomics/DSM/"$PROJFOLDER"/MZXML"
nfiles="$(find ${INFOLDER} -type f -name '*.mzXML' | wc -l)"
nfiles_nw="$(echo -e ${nfiles} | tr -d '[:space:]')"
qsub -N xcms_call_ibd -t 1-${nfiles_nw}:1 xcms_runcall.sh ${INFOLDER} -cwd

# ------

qsub -N xcms_end xcms_runend.sh ${OUTFOLDER}

# ------
echo 'Submitted all jobs! \(. w.)/'

# --- cleanup / zipping ---

tar -czvf /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/Chicken_pilot/chick_data_cent.tar.gz /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/Chicken_pilot/data
tar -czvf /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/Chicken_pilot/chick_data_profile.tar.gz /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/Chicken_pilot/data_profile

tar -czvf /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/Project2016_018_IBD/ibd_data_profile.tar.gz /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/Project2016_018_IBD/data_profile
tar -czvf /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/Project2016_018_IBD/ibd_data_mzxml.tar.gz /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/Project2016_018_IBD/MZXML

tar -czvf /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/RES-2017-04-26_DBS_Turkey_DSM/turkey_raw.tar.gz /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/RES-2017-04-26_DBS_Turkey_DSM/RAW
tar -czvf /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/RES-2017-05-01_DBS_France_DSM/france_raw.tar.gz /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/RES-2017-05-01_DBS_France_DSM/RAW
tar -czvf /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/RES-2017-05-02_DBS_France-Spain_DSM/france_spain_raw.tar.gz /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/RES-2017-05-02_DBS_France-Spain_DSM/RAW
tar -czvf /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/RES-2017-05-03_DBS_Spain-Gent_DSM/spain_ghent_raw.tar.gz /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/RES-2017-05-03_DBS_Spain-Gent_DSM/RAW
tar -czvf /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/RES-2017-10-31_DSM_DBS_Brazil_1/brazil_1_raw.tar.gz /hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/RES-2017-10-31_DSM_DBS_Brazil_1/RAW

