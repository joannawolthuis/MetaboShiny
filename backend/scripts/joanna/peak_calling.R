library(xcms)
library(msdata)
library(MassSpecWavelet)

# This documentation uses raw mzdata files from msdata as example data set. Assuming
# that msdata is installed, we locate the path of the package and extract the datafiles.

mzdatapath <- system.file("fticr", package = "msdata")
mzdatafiles <- list.files(mzdatapath, recursive = TRUE, full.names = TRUE)
cat("Starting xcmsDirect.Rnw")

# The xcmsSet-Constructor parses the given files and applies peakpicking using the
# MassSpecWavelet algorithm, leading to a xcmsSet object with 2 sampleclasses, ham4
# and ham5, and 5 samples, respectively

mzdatafiles

data.mean <- "data.mean"

xs <- xcmsSet(
   method="MSW",
   files=mzdatafiles,
   scales=c(1,4,9),
   nearbyPeak=T,
   verbose.columns = FALSE,
   winSize.noise=500,
   SNR.method="data.mean",
   snthr=10,
   BPPARAM = 3
   )

# calibrate can be used to correct the m/z values in a xcmsSet. It needs a xcmsSet and
# a list of m/z value which should be found in the object. To show this on a example a
# sample of ham4 is created and discalibrated a bit after getting some m/z:

mzdatapath <- "/Users/jwolthuis/MountPoint/jwolthuis/Data/Metabolomics/Pilot/Profile/Pig"
mzdatafiles <- list.files(mzdatapath, recursive = TRUE, full.names = TRUE,pattern=".mzXML$")
mzdatafiles[1]

xs_eur <- xcmsSet(
   method = "MSW",
   files = mzdatafiles,
   scales = c(1,4,9),
   nearbyPeak = T,
   verbose.columns = FALSE,
   winSize.noise = 500,
   SNR.method = "data.mean",
   snthr = 10)

masslist <- xs4@peaks
masslist

xs4@peaks[,"mz"] <- xs4@peaks[,"mz"] + 0.00001*runif(1,0,0.4)*xs4@peaks[,"mz"] + 0.0001

# The xcmsSet now can be calibrated again with the m/z from the masslist. The plot
# shows the reference masses with the distances to the found ones and the regression-line.


xs4c <- calibrate(xs4,
                     calibrants=masslist,
                     method="edgeshift",
                     mzabs=0.0001,
                     mzppm=5,
                     neighbours=3,
                     plotres=TRUE
                     )

# The method ”shift” adds a value to each m/z, ”linear” does a regression and edgeshift
# does a regression but uses a shift before the smallest and after the biggest m/z from the
# calibrants.
# These steps are necessary to create a usable input for mzClust. However, if you have
# already stored the data in a xcmsSet, you can skip the steps above.

# ---------------------------

# Now we can align xs with mzClust. The result is a clone of xs enhanced by the result of
# mzClust. For a description of the arguments mzClust takes, see helppage of the function.

xsg <- group(xs, method="mzClust")

groups(xsg)[1:10,]

peaks(xsg)[groupidx(xsg)[[1]]]

# In most cases not all samples are in one group. This can be the origin of serious problems
# in code, which is based on e.g. groupval. groupval sets missing peaks to NA. The solution
# is fillPeaks. It changes all NA values to random noise based on the raw data file.

groupval(xsg)[1,]

xsgf <- fillPeaks(xsg, method="MSW")

groupval(xsgf, "medret", "into")[1:10,]

reporttab <- diffreport(xsgf, "ham4", "ham5", "example", eicmax=4,
                        + h=480, w=640)

reporttab[1:4,]

  