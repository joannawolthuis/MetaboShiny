# --- get array num ---

args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
i <- as.numeric(args[2])

# --- load libraries ---

library(xcms)
f = "ABCDEFG_ASDASDAS_AASDASD_neg.cdf"
gsub(basename(f), pattern="_[neg|pos].*$", replacement="")
xcmsSet(snames = gsub(basename(f), pattern="_.*$", replacement=""))
# --- set defaults ---

scale <- seq(1,64,2)

outfolder = file.path(folder, "peakcalled")

files = list.files(folder, pattern="\\.mzXML", full.names = T)

f <- files[i]

# --- own function(s) ---

getSumSpec = function (object, ...) 
{
  sel <- profRange(object, ...)
  print(sel)
  scans <- list(length(sel$scanidx))
  uniquemz <- numeric()
  for (i in seq(along = sel$scanidx)) {
    currScan <- getScan(object, sel$scanidx[i], sel$mzrange)
    uniquemz <- unique(c(uniquemz, currScan[, "mz"]))
    scans[[i]] <- currScan
  }
  uniquemz <- sort(uniquemz)
  intmat <- matrix(nrow = length(uniquemz), ncol = length(sel$scanidx))
  for (i in seq(along = sel$scanidx)) {
    scan <- getScan(object, sel$scanidx[i], sel$mzrange)
    intmat[, i] <- approx(scan, xout = uniquemz)$y
  }
  # remove NA rows
  intmat[is.na(intmat)] <- 0
  # --------------
  points <- cbind(mz = uniquemz, intensity = rowSums(intmat))
  points
}

# load raw file, specify scan range here
xraw_pos <- xcmsRaw(f,profstep = 0.01,scanrange = c(10,150))

# merge profiles

pos_profile <- xraw_pos@env$profile

p <- pos_profile[,1]
for(i in 2:ncol(pos_profile)){
  vec <- pos_profile[,i]
  p = p + vec
}
pmat <- matrix(p)

# merge scans

scansum <- as.data.frame(getSumSpec(xraw_pos))
xr_sum <- xraw_pos
xr_sum@env$mz <- scansum$mz
xr_sum@env$intensity <- scansum$intensity
xr_sum@env$profile <- pmat
xr_sum@scanindex <- as.integer(0)
xr_sum@scantime <- xr_sum@scantime[1]
xr_sum@tic <- sum(scansum$intensity)
xr_sum@polarity <- as.factor("positive")
xr_sum@acquisitionNum <- as.integer(1)
xr_sum@scanrange <- 1

# write to cdf

cdf_pos = file.path(outfolder,
                    gsub(basename(f),
                         pattern="\\.mzXML",
                         replacement="_pos.cdf"))

write.cdf(xr_sum,file=cdf_pos)

# load in and call peaks
xset_pos <- xcmsSet(method="MSW",
                    files=cdf_pos,
                    scales=scale,
                    snthresh=3,
                    nearbyPeak=T)


if(!is.null(xset_pos)) npeaks <- nrow(xset_pos@peaks)

# --- report back ---
print(paste("Found", npeaks, "peaks in file", basename(f), "in positive mode!", sep=" "))

# --------------------------------------

npeaks = 0

# load raw file, specify scan range here
xraw_neg <- xcmsRaw(f,profstep = 0.01,scanrange = c(151,290))

# merge profiles

neg_profile <- xraw_neg@env$profile

p <- neg_profile[,1]
for(i in 2:ncol(neg_profile)){
  vec <- neg_profile[,i]
  p = p + vec
}
pmat <- matrix(p)

# merge scans

scansum <- as.data.frame(getSumSpec(xraw_neg))
xr_sum <- xraw_neg
xr_sum@env$mz <- scansum$mz
xr_sum@env$intensity <- scansum$intensity
xr_sum@env$profile <- pmat
xr_sum@scanindex <- as.integer(0)
xr_sum@scantime <- xr_sum@scantime[1]
xr_sum@tic <- sum(scansum$intensity)
xr_sum@polarity <- as.factor("negative")
xr_sum@acquisitionNum <- as.integer(1)
xr_sum@scanrange <- 1

# write to cdf

cdf_neg = file.path(outfolder,
                    gsub(basename(f),
                         pattern="\\.mzXML",
                         replacement="_neg.cdf"))

write.cdf(xr_sum,file=cdf_neg)

# load in and call peaks
xset_neg <- xcmsSet(method="MSW",
                    files=cdf_neg,
                    scales=scale,
                    snthresh=3,
                    nearbyPeak=T)


if(!is.null(xset_neg)) npeaks <- nrow(xset_neg@peaks)

# --- report back ---
print(paste("Found", npeaks, "peaks in file", basename(f), "in negative mode!", sep=" "))

# ---- save ----

save(xset_pos,
     xset_neg,
     file=file.path(outfolder, 
                    gsub(basename(f), 
                         pattern="\\.mzXML", 
                         replacement="_called.RData")
     ))