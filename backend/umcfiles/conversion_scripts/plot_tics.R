library(xcms)
library(pbapply)
library(parallel)

cl = makeCluster(4, type="FORK") 
##  NB: peak finding, peak group finding and fill missing done on HPC !

## Set directory where output files are going to be on new Windows PC
outdir <- "/Users/jwolthuis/PROCESSING_HPC/"   # project dir E:\Metabolomics\Project2016_03_biggen

# source("C:/Users/mraves/Metabolomics/AddOnFunctions.R")
# source("E:/Metabolomics/Mia_rescuedC/Metabolomics/AddOnFunctions.R") # Temporary
# source("Z:/Project2015_003_Toth/AddOnFunctions.R") # More complete in adducts
# source add on functions

## Path where sample directories are stored    
# Centroided data:
# filepath <- "Y:/Diagnostiek/Analyse_Data_Metabonomics/Q-Exact/BSP-2015-04-14/Bioinformatics Toth/mzXML_centroid" 

xmlfiles <- list.files(outdir, recursive=TRUE, full.names=TRUE, pattern="*.mzXML")

resol <- 140000

# theor.MZ <- read.table(file="C:/Users/mraves/Metabolomics/TheoreticalMZ_NegPos.txt", sep="\t", header=TRUE)
options(digits=16)

### plots and checks  ### include checks on "good" versus "bad" inputfiles?
# Ctrl <- xmlfiles.centroid[31]
# Ctrl <- xmlfiles[1]  # Ctrl2 <- xmlfiles.centroid[4]
# rawCtrl <- xcmsRaw(Ctrl, profstep=0.01) # NB 0.001 is too big for plotSurf
# rawCtrl  # check mz range:  69.3032-606.0467. Time range: 0.3-181.1 seconds (0-3 minutes)
# rawCtrl@scantime  # check RT range: 0.3 to 181.1 in steps of 0.6 seconds
# rawCtrl@polarity  # 155 positive scans, 150 negative
# Check out matrix of intensities
#allY <- rawMat(rawCtrl)  # matrix with time, mz, intensity for Pos and Neg
#allY.subset <- allY[1:100000, ]
#write.table(allY.subset, file=paste("Profile_C1A_022.txt", sep=""), sep="\t")
# plotChrom(rawCtrl, mz=c(172, 172.2), rt=c(0,max(rawCtrl@scantime)))  # 13C6 Phe in Pos
# plotChrom(rawCtrl, mz=c(188.08, 188.12), rt=c(0,max(rawCtrl@scantime)))  # 13C6 Tyr in Pos
# plotChrom(rawCtrl, mz=c(186.07, 186.10), rt=c(0,max(rawCtrl@scantime)))  # 13C6 Tyr in Neg

# plotTIC(rawCtrl, ident=FALSE, msident=FALSE) # ident=TRUE waits for mouse input; hit Esc

# save all TICs to file
TICdir <- paste(outdir, "/TICs", sep="")

dir.create(TICdir)

for(x in 1:length(xmlfiles)){
  rawF <- xcmsRaw(xmlfiles[x], profstep=0.1)
  png(filename=paste(TICdir, "/", sprintf("%04d", x), ".png", sep=""), 320, 240)
  plotTIC(rawF, ident=FALSE, msident=FALSE)
  dev.off()
  }