dims <- function(xmlfile,outdir,thresh,trim,resol){
# thresh=dimsThresh
  
  trimLeft=NULL
  trimRight=NULL
  breaks.fwhm=NULL
  breaks.fwhm.avg=NULL
  bins=NULL
  posRes=NULL
  negRes=NULL
  
  x = NULL
  try({x = xcmsRaw(xmlfile)}, silent = TRUE)
  if (is.null(x)){
    return(NULL)
  }

  load(paste(outdir, "breaks.fwhm.RData", sep="/"))
  
  # Create empty placeholders for later use
  bins <- rep(0, length(breaks.fwhm) - 1)
  
  # Generate a matrix
  y <- rawMat(x)
  
  # Get time values for positive and negative scans
  posTimes <- x@scantime[x@polarity == "positive"]
  negTimes <- x@scantime[x@polarity == "negative"]
  
  # Select scans where sample is present
  posTimes <- posTimes[posTimes > trimLeft & posTimes < trimRight]
  negTimes <- negTimes[negTimes > trimLeft & negTimes < trimRight]
  
  
  # Generate an index with which to select values for each mode
  posInd <- which(y[,"time"] %in% posTimes)
  negInd <- which(y[,"time"] %in% negTimes)
  
  # Separate each mode into its own matrix
  posY <- y[posInd,]
  negY <- y[negInd,]
  
  # Get index for binning intensity values
  ## This doesn't round the value for mz - is this an issue?
  yp <- cut(posY[,"mz"], breaks.fwhm, include.lowest=TRUE, right=TRUE, labels=FALSE)
  yn <- cut(negY[,"mz"], breaks.fwhm, include.lowest=TRUE, right=TRUE, labels=FALSE)

#     Z <- seq(from=1, to=10, by=0.5)
#     cut(Z, breaks = 1:10, include.lowest=TRUE, right=TRUE, labels=FALSE)

  # Empty the bins
  posBins<-bins
  negBins<-bins
  
  # Get the list of intensity values for each bin, and add the
  # intensity values which are in the same bin
  if (nrow(posY) > 0) {
#       ap <- aggregate(posY[,"intensity"],list(yp),sum)
#       posBins[ap[,1]] <- posBins[ap[,1]] + ap[,2] / length(posTimes)
     ap <- aggregate(posY[,"intensity"],list(yp), FUN = function(x){if (is.na(mean(x[which(x>thresh)]))){
                                                                      0 
                                                                    } else {
                                                                      mean(x[which(x>thresh)])
                                                                    }})
     posBins[ap[,1]] <- ap[,2]
    
  }
  if (nrow(negY) > 0) {
#       an <- aggregate(negY[,"intensity"],list(yn),sum)
#       negBins[an[,1]] <- negBins[an[,1]] + an[,2] / length(negTimes)
     an <- aggregate(negY[,"intensity"],list(yn), FUN = function(x){if (is.na(mean(x[which(x>thresh)]))){
                                                                      0 
                                                                    } else {
                                                                      mean(x[which(x>thresh)])
                                                                    }})
     negBins[an[,1]] <- an[,2]
  }

  # Zero any values that are below the threshold
  posBins[posBins < thresh] <- 0
  negBins[negBins < thresh] <- 0
  
  posRes = cbind(posRes, posBins)
  negRes = cbind(negRes, negBins)
  # }
  
  #which(posRes[,3]!=0)
  posRes = t(posRes)
  negRes = t(negRes)
  
  label=gsub(".+/(.+)\\.mzXML$", "\\1", xmlfile, ignore.case=TRUE)
  
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
#  sumsp <- apply(posResT,1,sum)
#  posResT.nonzero <- posResT[(sumsp != 0), ] # <=============================================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  negResT <- t(negRes)
#  sums <- apply(negResT,1,sum)
#  negResT.nonzero <- negResT[(sums != 0), ] # <=============================================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  return(list("pos"=posResT,"neg"=negResT, "breaksFwhm"=breaks.fwhm))
}

dims_raw <- function(rawfile,
                     outdir,
                     thresh=100,
                     trim=0.1,
                     resol=140000, 
                     scriptdir, 
                     mzmin=70, 
                     mzmax=600){
  
  filename = file.path(outdir,
                       "spectra",
                       paste0(gsub(basename(rawfile),
                                   pattern="\\.raw",
                                   replacement=""),
                              ".RData"))
  
  if(file.exists(filename)) return(NULL)
  
  # thresh=dimsThresh
  mzrange = c(mzmin, 
              mzmax)
  
  #rawfile = "/Users/jwolthuis/Documents/umc/data/Data/BrazilAndSpain/MZXML/RES_20171031_016.raw"
  
  load(file.path(outdir, "breaks.RData"))
  
  cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/scantimes' -raw '$rawfile'")
  
  scantimes <- as.numeric(system(cmd,intern = T))
  
  trimLeft = round(scantimes[length(scantimes)*trim],
                   digits = 22)  
  trimRight = round(scantimes[length(scantimes)*(1-trim)], 
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
       file=filename)
}
