#!/bin/bash
#####################################################################
#####################################################################
user=marcel
#####################################################################
#####################################################################

# scripts="./scripts"
# outdir="./results"
# inpdir="./data"

#scripts="/hpc/shared/dbg_mz/$user/Direct-Infusion-Pipeline_2.1/scripts"
#outdir="/hpc/shared/dbg_mz/$user/Direct-Infusion-Pipeline_2.1/results"
#inpdir="/hpc/shared/dbg_mz/$user/Direct-Infusion-Pipeline_2.1/data"

scripts="/hpc/cog_bioinf/ridder/users/jwolthuis/Pipelines/changed/Direct-Infusion-Pipeline_2.1/scripts"
projdir="/hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/RES-2017-10-31_DSM_DBS_Brazil_1"
inpdir="${projdir}/MZXML"
outdir="${projdir}/res_2.1_flexdb"

#thresh=5000
#resol=140000
#trim=0.1

thresh_pos=2000
thresh_neg=2000
dimsThresh=100
resol=140000
trim=0.1
nrepl=3 # 3 or 5!

# thresh2remove = 1*10^9 ============> set in averageTechReplicates.R
# thresh2remove = 5*10^8

echo "`pwd`"

it=0
find $inpdir -iname "*.mzXML" | while read mzXML;
do
echo "Processing file $mzXML"
it=$((it+1))
if [ $it == 1 ]; then
qsub -l h_rt=01:00:00 -l h_vmem=16G -N "breaks" $scripts/runGenerateBreaks.sh $mzXML $outdir $trim $resol $scripts $nrepl
fi

qsub -l h_rt=01:00:00 -l h_vmem=16G -N "dims" -hold_jid "breaks" $scripts/runDIMS.sh $mzXML $scripts $outdir $trim $dimsThresh $resol
done

qsub -l h_rt=01:00:00 -l h_vmem=16G -N "average" -hold_jid "dims" $scripts/runAverageTechReps.sh $scripts $outdir $nrepl

scanmode="negative"
qsub -l h_rt=01:00:00 -l h_vmem=16G -N "queueFinding_$scanmode" -hold_jid "average" $scripts/queuePeakFinding.sh $scripts $outdir $inpdir $thresh_neg $resol $scanmode
scanmode="positive"
qsub -l h_rt=01:00:00 -l h_vmem=16G -N "queueFinding_$scanmode" -hold_jid "average" $scripts/queuePeakFinding.sh $scripts $outdir $inpdir $thresh_pos $resol $scanmode
