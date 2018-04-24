require(data.table)
require(MALDIquant)


args = commandArgs(trailingOnly=TRUE)

args = c("/hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/BrSp", 1)

# load file

scriptdir = normalizePath("/hpc/cog_bioinf/ridder/users/jwolthuis/MassChecker/scripts")
resdir = file.path(args[1], "results")

files <- list.files(path = file.path(resdir, "averaged"),
                    full.names = T)

current <- files[as.numeric(args[2])]

for(f in list.files(file.path(scriptdir), full.names = T)) source(f)
for(f in list.files(file.path(scriptdir, "AddOnFunctions"), full.names = T)) source(f)

load(current)

df <- matrix(data = averaged[[1]]@intensity,
             ncol=1)
rownames(df) <- averaged[[1]]@mass

values = list("mean" = NULL,
              "area" = NULL,
              "nr" = NULL,
              "min" = NULL,
              "max" = NULL,
              "qual" = NULL,
              "spikes" = 0)

int.factor = 1*10^5
factor = FALSE
resol = 140000
sampname = names(averaged)
plot = FALSE
thresh = 2000
scale = 1
scanmode = "pos"
peaks = searchMZRange(df[,1],values,int.factor,scale,resol,outdir,sampname,scanmode,FALSE,100,100,thresh)  
validPeaks = which(peaks$qual > 0)
mzs <- as.numeric(peaks$mean[validPeaks])
ints <- as.numeric(peaks$area[validPeaks])
qual <- as.numeric(peaks$qual[validPeaks])

peaks_gauss = MALDIquant::createMassPeaks(mass = as.numeric(peaks$mean[validPeaks]),
                                          intensity = peaks$area[validPeaks],
                                          metaData = list(qual=qual))