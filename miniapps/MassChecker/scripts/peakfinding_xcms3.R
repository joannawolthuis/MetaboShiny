## ----load-libs, message = FALSE, results = "hide"--------------------------
library(xcms)
library(MassSpecWavelet)
library(msdata)
library(data.table)

register(SerialParam())
register(MulticoreParam())

## Loading a small subset of direct injection, single spectrum files

examplefiles <- list.files(system.file("fticr", package = "msdata"),
                     recursive = TRUE, full.names = TRUE)
myfiles <- list.files("/Users/jwolthuis/MountPoint/Data/Metabolomics/DSM/BrazilFarm1and2",
                      recursive = TRUE, full.names = TRUE, pattern = ".mzXML")

myfiles_1 <- readMSData(myfiles[1:2], 
                          msLevel. = 1, 
                          mode = "onDisk",
                          verbose = T, 
                          centroided=F)

examp_raw1 <- readMSData(examplefiles[1:2], msLevel. = 1, mode = "onDisk")
examp_raw2 <- readMSData(examplefiles[3:4], msLevel. = 1, mode = "onDisk")



intensity(myfiles_raw)
intensity(bin(myfiles_raw, binSize=2))
msnset.comb <- combineFeatures(msnset, grp, "sum")
msnset.nrm <- normalise(myfiles_raw, "sum")
msnset.nrm
quantify(myfiles_raw, method = "trap")

#> Error in .nextMethod(.Object, ...): object 'tan2009r1' not found
avg <- averageMSnSet(x)

## Perform the MSW peak detection on these:
p <- MSWParam(scales = c(1, 7), peakThr = 80000, ampTh = 0.005,
              SNR.method = "data.mean", winSize.noise = 500)

examppeaks1 <- findChromPeaks(examp_raw1, param = p)
examppeaks2 <- findChromPeaks(examp_raw2, param = p)

c(examppeaks1, examppeaks2)

examppeaks@msFeatureData

mytestpeaks <- mypeaks
tofix <- data.table(mytestpeaks@featureData@data,keep.rownames = T)
samples <- unique(gsub(tofix$rn, pattern = "\\..*$",replacement = ""))

mytestpeaks@featureData@data
mytestpeaks@featureData@varMetadata

## Now create the MzClustParam parameter object: we're assuming here that
## both samples are from the same sample group.
p <- MzClustParam(sampleGroups = c(1, 2))
mytestpeaks_grouped <- groupChromPeaks(mytestpeaks, param = p)

exampeaks_grouped <- groupChromPeaks(examppeaks, param = p)

## Get the definition of the features.
featureDefinitions(mytestpeaks_grouped)


mytestpeaks_filled <- fillChromPeaks(mytestpeaks_grouped, param = FillChromPeaksParam())

exampeaks_filled <- fillChromPeaks(exampeaks_grouped, param = FillChromPeaksParam())

chromPeaks(exampeaks_filled)
featureDefinitions(exampeaks_filled)
featureValues(exampeaks_filled)

## Extract the mz and intensity values
mz(exampeaks_filled, bySample = TRUE)
intensity(exampeaks_filled, bySample = TRUE)

## ----message = FALSE-------------------------------------------------------
## Subset to the first file.
first_file <- filterFile(ham_prep, file = 1)

## Extract 3 m/z values
calib_mz <- chromPeaks(first_file)[c(1, 4, 7), "mz"]
calib_mz <- calib_mz + 0.00001 * runif(1, 0, 0.4) * calib_mz + 0.0001


## ----message = FALSE-------------------------------------------------------
## Set-up the parameter class for the calibration
prm <- CalibrantMassParam(mz = calib_mz, method = "edgeshift",
                          mzabs = 0.0001, mzppm = 5)
first_file_calibrated <- calibrate(first_file, param = prm)


## ----calibration-result, fig = TRUE, fig.align = "center"------------------
diffs <- chromPeaks(first_file_calibrated)[, "mz"] -
  chromPeaks(first_file)[, "mz"]

plot(x = chromPeaks(first_file)[, "mz"], xlab = expression(m/z[raw]),
     y = diffs, ylab = expression(m/z[calibrated] - m/z[raw]))


## ----correspondence, message = FALSE, results = "hide"---------------------
## Using default settings but define sample group assignment
mzc_prm <- MzClustParam()
ham_prep <- groupChromPeaks(ham_prep, param = mzc_prm)


## --------------------------------------------------------------------------
ham_prep

## --------------------------------------------------------------------------
featureDefinitions(ham_prep)

## ----feature-FT01, fig = TRUE, fig.width = 6, fig.height = 4, fig.align = "center"----
## Get the peaks belonging to the first feature
pks <- chromPeaks(ham_prep)[featureDefinitions(ham_prep)$peakidx[[1]], ]

## Define the m/z range
mzr <- c(min(pks[, "mzmin"]) - 0.001, max(pks[, "mzmax"]) + 0.001)

## Subset the object to the m/z range
ham_prep_sub <- filterMz(ham_prep, mz = mzr)

## Extract the mz and intensity values
mzs <- mz(ham_prep_sub, bySample = TRUE)
ints <- intensity(ham_prep_sub, bySample = TRUE)

## Plot the data
plot(3, 3, pch = NA, xlim = range(mzs), ylim = range(ints), main = "FT01",
     xlab = "m/z", ylab = "intensity")
## Define colors
cols <- rep("#ff000080", length(mzs))
cols[ham_prep_sub$sample_group == "ham5"] <- "#0000ff80"
tmp <- mapply(mzs, ints, cols, FUN = function(x, y, col) {
  points(x, y, col = col, type = "l")
})


## --------------------------------------------------------------------------
feat_vals <- featureValues(ham_prep, value = "into")
head(feat_vals)


## ----fillpeaks, message = FALSE--------------------------------------------
ham_prep <- fillChromPeaks(ham_prep, param = FillChromPeaksParam())

head(featureValues(ham_prep, value = "into"))
