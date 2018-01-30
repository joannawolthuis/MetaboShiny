dims <- function(xmlfile,outdir,thresh,trim,resol){
# thresh=dimsThresh
  
  print(thresh)
  
  trimLeft=NULL
  trimRight=NULL
  breaks.fwhm=NULL
  breaks.fwhm.avg=NULL
  bins=NULL
  posRes=NULL
  negRes=NULL
  
  print("here1")
  x = NULL
  try({x = xcmsRaw(xmlfile)}, silent = TRUE)
  if (is.null(x)){
    return(NULL)
  }

  load(paste(outdir, "breaks.fwhm.RData", sep="/"))
  
  # Create empty placeholders for later use
  bins <- rep(0,length(breaks.fwhm)-1)
  
  # Generate a matrix
  y <- rawMat(x)
  
  print("here2")
  
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
  
  print("here3")
  
  # Get index for binning intensity values
  ## This doesn't round the value for mz - is this an issue?
  yp <- cut(posY[,"mz"], breaks.fwhm, include.lowest=TRUE, right=TRUE, labels=FALSE)
  yn <- cut(negY[,"mz"], breaks.fwhm, include.lowest=TRUE, right=TRUE, labels=FALSE)

#     Z <- seq(from=1, to=10, by=0.5)
#     cut(Z, breaks = 1:10, include.lowest=TRUE, right=TRUE, labels=FALSE)

  # Empty the bins
  posBins<-bins
  negBins<-bins
  
  print("here4")
  
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

  print("here5")
  
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
  
  print("here6")
  
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
  
  print("here7")
  
  
  # omit rows with only zeros
  posResT <- t(posRes)
#  sumsp <- apply(posResT,1,sum)
#  posResT.nonzero <- posResT[(sumsp != 0), ] # <=============================================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  negResT <- t(negRes)
#  sums <- apply(negResT,1,sum)
#  negResT.nonzero <- negResT[(sums != 0), ] # <=============================================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  return(list("pos"=posResT,"neg"=negResT, "breaksFwhm"=breaks.fwhm))
}
