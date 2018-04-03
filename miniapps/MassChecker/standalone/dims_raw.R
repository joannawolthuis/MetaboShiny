dims_raw <- function(infile, #.raw file
                     outfile, #.mzxml name and location 
                     outdir,
                     thresh=100,
                     trim=0.1,
                     resol=140000, 
                     scriptdir, 
                     mzmin=70, 
                     mzmax=600){
  
  require(data.table)
  require(MALDIquant)
  require(gsubfn)
  
  # --------------------------]
  
  if(file.exists(filename)) return(NULL)
  
  # thresh=dimsThresh
  mzrange = c(mzmin, 
              mzmax)
  
  load(file.path(outdir, "breaks.RData"))
  
  cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/scantimes' -raw '$rawfile'")
  
  scantimes <- as.numeric(system(cmd,intern = T))
  
  trimLeft = round(scantimes[length(scantimes) * trim],
                   digits = 22)  
  trimRight = round(scantimes[length(scantimes) * (1-trim)], 
                    digits=22)
  
  bins <- rep(0, length(breaks.fwhm) - 1)
  
  posRes=NULL
  negRes=NULL
  
  # Generate a matrix
  cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/scanpolarity' -raw '$rawfile'")
  polarity <- system(cmd, intern = T)
  
  # Get time values for positive and negative scans
  posTimes <- scantimes[polarity == "positive"]
  negTimes <- scantimes[polarity == "negative"]
  
  # Select scans where sample is present
  posTimes <- posTimes[posTimes > trimLeft & posTimes < trimRight]
  negTimes <- negTimes[negTimes > trimLeft & negTimes < trimRight]
  
  # Generate an index with which to select values for each mode
  posInd <- which(scantimes %in% posTimes)
  negInd <- which(scantimes %in% negTimes)
  
  # get mzvals
  pos = lapply(posInd, FUN=function(scan){
    cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/printspectrum' -raw '$rawfile' -sn $scan")
    spec <- system(cmd,
                   intern = T)
    temp.list <- strsplit(spec, " ")
    mzvals = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[1]])))
    intensities = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[2]])))
    # --- mass spectrum obj ---
    data.table(mz = mzvals, 
               i = intensities)
  })
  neg = lapply(negInd, FUN=function(scan){
    cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/printspectrum' -raw '$rawfile' -sn $scan")
    spec <- system(cmd,intern = T)
    temp.list <- strsplit(spec, " ")
    mzvals = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[1]])))
    intensities = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[2]])))
    # --- mass spectrum obj ---
    data.table(mz = mzvals, 
               i = intensities)
  })
  
  # Separate each mode into its own matrix
  posY <- rbindlist(pos)
  negY <- rbindlist(neg)
  
  # Get index for binning intensity values
  ## This doesn't round the value for mz - is this an issue?
  
  setkey(posY, mz)
  setkey(negY, mz)
  
  posY <<- posY[mz %between% mzrange]
  negY <<- negY[mz %between% mzrange]
  
  posY[, bin:= cut(mz, breaks.fwhm, include.lowest=TRUE, right=TRUE, labels=FALSE)]
  negY[, bin:= cut(mz, breaks.fwhm, include.lowest=TRUE, right=TRUE, labels=FALSE)]
  
  # Empty the bins
  posBins <- bins
  negBins <- bins
  
  # Get the list of intensity values for each bin, and add the
  # intensity values which are in the same bin
  
  if (nrow(posY) > 0) {
    posY[, aggr := sum(i), by=bin]
    ap <- posY[!is.na(bin),c("bin", "aggr")]
    posBins[ap$bin] <- ap$aggr
  }
  if (nrow(negY) > 0) {
    negY[, aggr := sum(i), by=bin]
    an <- negY[!is.na(bin),c("bin", "aggr")]
    negBins[an$bin] <- an$aggr
  }
  
  # Zero any values that are below the threshold
  
  #posBins[posBins < thresh] <- 0
  #negBins[negBins < thresh] <- 0
  
  posRes = posBins
  negRes = negBins
  
  posRes = t(posRes)
  negRes = t(negRes)
  
  label=gsub(".+/(.+)\\.raw$", "\\1", rawfile, ignore.case=TRUE)
  
  # Add in file names as row names
  rownames(posRes) = label
  rownames(negRes) = label
  
  # Add 0.5 to the values in breaks.fwhm, and delete the last value
  a <- breaks.fwhm.avg[-length(breaks.fwhm.avg)]  # + 0.5*deltaMZ
  
  # Format as string and show precision of float to 2 digits
  b <- sprintf("%.5f",a)
  
  # Use this as the column names
  colnames(posRes) <- b
  colnames(negRes) <- b
  
  # omit rows with only zeros
  posResT <- t(posRes)
  negResT <- t(negRes)
  
  pos_spec = MALDIquant::createMassSpectrum(as.numeric(rownames(posResT)), 
                                            c(posResT), 
                                            metaData = list(sample=label))
  neg_spec = MALDIquant::createMassSpectrum(as.numeric(rownames(negResT)), 
                                            c(negResT), 
                                            metaData = list(sample=label))
  # --- return ---
  
  pklist <- list("pos"=pos_spec,"neg"=neg_spec)
  
  save(x=pklist, 
       file=outfile)
}


require(data.table)
require(MALDIquant)
require(gsubfn)
require(pbapply)

rawfile = "/Users/jwolthuis/Documents/umc/data/Data/BrSpIt/MZXML/RES_20171031_001.raw"
scriptdir = "/Users/jwolthuis/Google Drive/MassChecker/scripts/"

cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/scantimes' -raw '$rawfile'")

scantimes <- as.numeric(system(cmd,intern = T))

posRes=NULL
negRes=NULL

# Generate a matrix
cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/scanpolarity' -raw '$rawfile'")
polarity <- system(cmd, intern = T)

# Get time values for positive and negative scans
posTimes <- scantimes[polarity == "positive"]
negTimes <- scantimes[polarity == "negative"]

# Generate an index with which to select values for each mode
posInd <- which(scantimes %in% posTimes)
negInd <- which(scantimes %in% negTimes)

# get mzvals
pos = pbapply::pblapply(posInd, FUN=function(scan){
  cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/printspectrum' -raw '$rawfile' -sn $scan")
  spec <- system(cmd,
                 intern = T)
  temp.list <- strsplit(spec, " ")
  mzvals = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[1]])))
  intensities = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[2]])))
  # --- mass spectrum obj ---
  data.table(mz = mzvals, 
             i = intensities)
})

neg = pbapply::pblapply(negInd, FUN=function(scan){
  cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/printspectrum' -raw '$rawfile' -sn $scan")
  spec <- system(cmd,intern = T)
  temp.list <- strsplit(spec, " ")
  mzvals = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[1]])))
  intensities = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[2]])))
  # --- mass spectrum obj ---
  data.table(mz = mzvals, 
             i = intensities)
})

# Separate each mode into its own matrix
posY <- rbindlist(pos)
negY <- rbindlist(neg)
