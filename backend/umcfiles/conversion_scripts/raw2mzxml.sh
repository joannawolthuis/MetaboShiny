#!/bin/sh

#  raw2mzxml.sh
#  
#
#  Created by jwolthuis on 16/10/2017.
#

#READW=/Users/jwolthuis/Applications/ReAdW.2016010.msfilereader.exe
#INFOLDER="/Users/jwolthuis/Downloads/HPCdownloads/TURKEY/RAW"
#OUTFOLDER="/Users/jwolthuis/Downloads/HPCdownloads/TURKEY/MZXML"

#for INFILE in ${INFOLDER}/*.raw
#do
#filename=$(basename "$INFILE")
#OUTFILE=${OUTFOLDER}/${filename%.*}.mzXML
#echo "Converting $INFILE to $OUTFILE ..."
#env WINEPREFIX=/Users/jwolthuis/.wine32 wine $READW $INFILE $OUTFILE

#done

# ------ ver that gets from and sends back to hpc --------

# --- LOCAL ---
READW=/Users/jwolthuis/Applications/ReAdW.2016010.msfilereader.exe
PROCFOLDER="/Users/jwolthuis/PROCESSING_HPC"

# --- HPC ---
HPCBASE="/hpc/cog_bioinf/ridder/users/jwolthuis"
INFOLDER="Data/Metabolomics/DSM/Turkey/RAW"
OUTFOLDER="Data/Metabolomics/DSM/Turkey/MZXML"

mkdir -p ${PROCFOLDER}

# get files to temp location (will get asked for password)

scp hpctransfer:${HPCBASE}/${INFOLDER}/*.raw ${PROCFOLDER}

# process
for INFILE in ${PROCFOLDER}/*.raw
do
filename=$(basename "$INFILE")
OUTFILE=${PROCFOLDER}/${filename%.*}.mzXML
#############################################
if [ -f $OUTFILE ]; then                    #
echo "File already converted! Skipping ..." #
continue                                    #
fi                                          #
#############################################
echo "Converting $INFILE to $OUTFILE ..."
env WINEPREFIX=/Users/jwolthuis/.wine32 wine $READW $INFILE $OUTFILE

done

# send files back
scp ${PROCFOLDER}/*.mzXML hpctransfer:${HPCBASE}/${OUTFOLDER}

# delete residuals ?
