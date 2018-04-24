require(data.table)
require(MALDIquant)
require(MassSpecWavelet)
require(pastecs)
require(pbapply)


args = commandArgs(trailingOnly=TRUE)

#args = c("/hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/BrSp", 1)

#args = c("/Users/jwolthuis/Documents/umc/data/Data/BrSpIt/MZXML", 637)

# load file

scriptdir = normalizePath("/hpc/cog_bioinf/ridder/users/jwolthuis/MassChecker/scripts")

#scriptdir = normalizePath("/Users/jwolthuis/Google Drive/MassChecker/scripts")


resdir = file.path(args[1], "results")

files <- list.files(path = file.path(resdir, "averaged"),
                    full.names = T)

current <- files[as.numeric(args[2])]

cl = 0
#cl = parallel::makeCluster(3, "FORK")

fn = file.path(resdir, "peaks_wavelet", basename(current))
if(file.exists(fn)) return(NULL)

print(current)

load(current)

for(f in list.files(file.path(scriptdir), full.names = T)) source(f)
for(f in list.files(file.path(scriptdir, "AddOnFunctions"), full.names = T)) source(f)

sampname = gsub(basename(current), pattern = "_(pos|neg)\\.RData", replacement = "")

# ----------------

df <- matrix(data = averaged[[1]]@intensity,
             ncol=1)
rownames(df) <- averaged[[1]]@mass

scales <- seq(1,5,0.5)

wCoefs <- cwt(df, scales=scales, wavelet='mexh')
localMax <- getLocalMaximumCWT(wCoefs,
                               amp.Th = 100)

ridgeList <- getRidge(localMax)

majorPeakInfo <- identifyMajorPeaks(df, 
                                    ridgeList, 
                                    wCoefs, 
                                    SNR.Th = sig_noise, 
                                    peakScaleRange = 1,
                                    nearbyPeak = FALSE)
# ----------------
keep <- which(df[majorPeakInfo$peakCenterIndex] > 0)
# ----------------
peaks = MALDIquant::createMassPeaks(mass = as.numeric(rownames(df)[keep]),
                                    intensity = majorPeakInfo$peakValue[keep],
                                    snr = majorPeakInfo$peakSNR[keep],
                                    metaData = list(sample = sampname))
# ----------------

if(!dir.exists(file.path(resdir, "peaks_wavelet"))) dir.create(file.path(resdir, "peaks_wavelet"))

print(fn)

save(x=peaks, file=fn)


