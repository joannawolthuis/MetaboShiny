require(MassSpecWavelet)

cwt_J <- function (ms, scales = 1, wavelet = "mexh") 
{
  if (wavelet == "mexh") {
    psi_xval <- seq(-8, 8, length = 1024)
    psi <- (2/sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) * 
      exp(-psi_xval^2/2)
  }
  else if (is.matrix(wavelet)) {
    if (nrow(wavelet) == 2) {
      psi_xval <- wavelet[1, ]
      psi <- wavelet[2, ]
    }
    else if (ncol(wavelet) == 2) {
      psi_xval <- wavelet[, 1]
      psi <- wavelet[, 2]
    }
    else {
      stop("Unsupported wavelet format!")
    }
  }
  else {
    stop("Unsupported wavelet!")
  }
  oldLen <- length(ms)
  ms <- MassSpecWavelet:::extendNBase(ms, nLevel = NULL, base = 2)
  len <- length(ms)
  nbscales <- length(scales)
  wCoefs <- NULL
  psi_xval <- psi_xval - psi_xval[1]
  dxval <- psi_xval[2]
  xmax <- psi_xval[length(psi_xval)]
  for (i in 1:length(scales)) {
    scale.i <- scales[i]
    f <- rep(0, len)
    j <- 1 + floor((0:(scale.i * xmax))/(scale.i * dxval))
    if (length(j) == 1) 
      j <- c(1, 1)
    lenWave <- length(j)
    f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
    if (length(f) > len) 
      stop(paste("scale", scale.i, "is too large!"))
    wCoefs.i <- 1/sqrt(scale.i) * convolve(ms, f)
    wCoefs.i <- c(wCoefs.i[(len - floor(lenWave/2) + 1):len], 
                  wCoefs.i[1:(len - floor(lenWave/2))])
    wCoefs <- cbind(wCoefs, wCoefs.i)
  }
  if (length(scales) == 1) 
    wCoefs <- matrix(wCoefs, ncol = 1)
  colnames(wCoefs) <- scales
  wCoefs <- wCoefs[1:oldLen, , drop = FALSE]
  return(wCoefs)
}

getLocalMaximumCWT_J <- function (wCoefs, minWinSize = 5, amp.Th = 0) 
{
  localMax <- NULL
  scales <- as.numeric(colnames(wCoefs))
  for (i in 1:length(scales)) {
    scale.i <- scales[i]
    winSize.i <- round(scale.i, 0) * 2 + 1
    if (winSize.i < minWinSize) {
      winSize.i <- minWinSize
    }
    temp <- localMaximum(wCoefs[, i], winSize.i)
    localMax <- cbind(localMax, temp)
  }
  localMax[wCoefs < amp.Th] <- 0
  colnames(localMax) <- colnames(wCoefs)
  rownames(localMax) <- rownames(wCoefs)
  return(localMax)
}

getRidge_J <- function (localMax, iInit = ncol(localMax), step = -1, iFinal = 1, 
                        minWinSize = 5, gapTh = 3, skip = NULL){
  scales <- as.numeric(colnames(localMax))
  if (is.null(scales)) 
    scales <- 1:ncol(localMax)
  maxInd_curr <- which(localMax[, iInit] > 0)
  nMz <- nrow(localMax)
  if (is.null(skip)) {
    skip <- iInit + 1
  }
  if (ncol(localMax) > 1) {
    colInd <- seq(iInit + step, iFinal, step)
  }
  else {
    colInd <- 1
  }
  ridgeList <- as.list(maxInd_curr)
  print(ridgeList)
  names(ridgeList) <- maxInd_curr
  peakStatus <- as.list(rep(0, length(maxInd_curr)))
  names(peakStatus) <- maxInd_curr
  orphanRidgeList <- NULL
  orphanRidgeName <- NULL
  nLevel <- length(colInd)
  for (j in 1:nLevel) {
    col.j <- colInd[j]
    scale.j <- scales[col.j]
    if (colInd[j] == skip) {
      oldname <- names(ridgeList)
      ridgeList <- lapply(ridgeList, function(x) c(x, 
                                                   x[length(x)]))
      names(ridgeList) <- oldname
      next
    }
    if (length(maxInd_curr) == 0) {
      maxInd_curr <- which(localMax[, col.j] > 0)
      next
    }
    winSize.j <- round(scale.j,0) * 2 + 1
    if (winSize.j < minWinSize) {
      winSize.j <- minWinSize
    }
    selPeak.j <- NULL
    remove.j <- NULL
    for (k in 1:length(maxInd_curr)) {
      ind.k <- maxInd_curr[k]
      start.k <- ifelse(ind.k - winSize.j < 1, 1, ind.k - 
                          winSize.j)
      end.k <- ifelse(ind.k + winSize.j > nMz, nMz, ind.k + 
                        winSize.j)
      ind.curr <- which(localMax[start.k:end.k, col.j] > 
                          0) + start.k - 1
      if (length(ind.curr) == 0) {
        status.k <- peakStatus[[as.character(ind.k)]]
        if (is.null(status.k)) 
          status.k <- gapTh + 1
        if (status.k > gapTh & scale.j >= 2) {
          temp <- ridgeList[[as.character(ind.k)]]
          orphanRidgeList <- c(orphanRidgeList, list(temp[1:(length(temp) - 
                                                               status.k)]))
          orphanRidgeName <- c(orphanRidgeName, paste(col.j + 
                                                        status.k + 1, ind.k, sep = "_"))
          remove.j <- c(remove.j, as.character(ind.k))
          next
        }
        else {
          ind.curr <- ind.k
          peakStatus[[as.character(ind.k)]] <- status.k + 
            1
        }
      }
      else {
        peakStatus[[as.character(ind.k)]] <- 0
        if (length(ind.curr) >= 2) 
          ind.curr <- ind.curr[which.min(abs(ind.curr - 
                                               ind.k))]
      }
      ridgeList[[as.character(ind.k)]] <- c(ridgeList[[as.character(ind.k)]], 
                                            ind.curr)
      selPeak.j <- c(selPeak.j, ind.curr)
    }
    if (length(remove.j) > 0) {
      removeInd <- which(names(ridgeList) %in% remove.j)
      ridgeList <- ridgeList[-removeInd]
      peakStatus <- peakStatus[-removeInd]
    }
    dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])
    if (length(dupPeak.j) > 0) {
      removeInd <- NULL
      for (dupPeak.jk in dupPeak.j) {
        selInd <- which(selPeak.j == dupPeak.jk)
        selLen <- sapply(ridgeList[selInd], length)
        removeInd.jk <- which.max(selLen)
        removeInd <- c(removeInd, selInd[-removeInd.jk])
        orphanRidgeList <- c(orphanRidgeList, ridgeList[removeInd.jk])
        orphanRidgeName <- c(orphanRidgeName, paste(col.j, 
                                                    selPeak.j[removeInd.jk], sep = "_"))
      }
      selPeak.j <- selPeak.j[-removeInd]
      ridgeList <- ridgeList[-removeInd]
      peakStatus <- peakStatus[-removeInd]
    }
    names(ridgeList) <- selPeak.j
    names(peakStatus) <- selPeak.j
    if (scale.j >= 2) {
      maxInd_next <- which(localMax[, col.j] > 0)
      unSelPeak.j <- maxInd_next[!(maxInd_next %in% selPeak.j)]
      newPeak.j <- as.list(unSelPeak.j)
      names(newPeak.j) <- unSelPeak.j
      ridgeList <- c(ridgeList, newPeak.j)
      maxInd_curr <- c(selPeak.j, unSelPeak.j)
      newPeakStatus <- as.list(rep(0, length(newPeak.j)))
      names(newPeakStatus) <- newPeak.j
      peakStatus <- c(peakStatus, newPeakStatus)
    }
    else {
      maxInd_curr <- selPeak.j
    }
  }
  names(ridgeList) <- paste(1, names(ridgeList), sep = "_")
  names(orphanRidgeList) <- orphanRidgeName
  ridgeList <- c(ridgeList, orphanRidgeList)
  ridgeList <- lapply(ridgeList, rev)
  ridgeList <- ridgeList[!duplicated(names(ridgeList))]
  attr(ridgeList, "class") <- "ridgeList"
  attr(ridgeList, "scales") <- scales
  return(ridgeList)
}

tuneInPeakInfo_J <- function (ms, majorPeakInfo = NULL, peakIndex = NULL, peakScale = NULL, 
                              maxScale = 128, ...) 
{
  
  keep.i = NULL
  
  if (!is.null(majorPeakInfo)) {
    if (!all(c("peakIndex", "peakCenterIndex", "peakScale") %in% 
             names(majorPeakInfo))) 
      stop("Format of majorPeakInfo is incorret!")
    peakIndex <- majorPeakInfo$peakIndex
    peakScale <- majorPeakInfo$peakScale[names(peakIndex)]
    peakCenterIndex <- majorPeakInfo$peakCenterIndex[names(peakIndex)]
    peakValue <- majorPeakInfo$peakValue[names(peakIndex)]
    peakSNR <- majorPeakInfo$peakSNR[names(peakIndex)]
  }
  else {
    if (is.null(peakIndex) | is.null(peakScale)) 
      stop("majorPeakInfo or peakIndex and peakScale should be provided!")
    peakCenterIndex <- peakIndex
    peakSNR <- NULL
    peakValue <- NULL
  }
  peakName <- names(peakIndex)
  if (is.null(peakIndex)) 
    peakName <- as.character(peakIndex)
  peakCenterIndex.new <- NULL
  peakScale.new <- NULL
  peakValue.new <- NULL
  unProcessedInd <- NULL
  pb <- pbapply::startpb(min = 1, max = length(peakIndex))
  for (i in 1:length(peakIndex)) {
    pbapply::setpb(pb, i)
    peak.i <- peakIndex[i]
    peakScale.i <- peakScale[i]
    
    scales.i <- seq(peakScale.i - 4, peakScale.i + 16, 
                    0.5)
    # 
    # if (peakScale.i + 4 < maxScale) {
    #   scales.i <- seq(peakScale.i - 4, peakScale.i + 4, 
    #                   0.5)
    # }
    # else {
    #   scales.i <- seq(maxScale - 10, maxScale, 0.5)
    # }
    winSize.i <- 16 * max(scales.i)
    if (peak.i - winSize.i/2 < 1) {
      winSize.i <- (peak.i - 1) * 2
      scales.i <- scales.i[scales.i < peak.i/16]
      start.i <- 1
    }
    else {
      start.i <- peak.i - winSize.i/2
    }
    if (peak.i + winSize.i/2 > length(ms)) {
      winSize.i <- (length(ms) - peak.i) * 2
      scales.i <- scales.i[scales.i < (length(ms) - peak.i)/16]
      end.i <- length(ms)
    }
    else {
      end.i <- peak.i + winSize.i/2
    }
    if (length(scales.i) <= 1) {
      peakScale.new <- c(peakScale.new, peakScale[i])
      peakValue.new <- c(peakValue.new, peakValue[i])
      peakCenterIndex.new <- c(peakCenterIndex.new, peakCenterIndex[i])
      unProcessedInd <- c(unProcessedInd, i)
      next
    }
    ms.i <- ms[start.i:end.i]
    wCoefs.i <- cwt(ms.i, scales = scales.i, wavelet = "mexh")
    localMax.i <- getLocalMaximumCWT_J(wCoefs.i, ...)
    if(!any(localMax.i > 0)) next
    keep.i = c(keep.i, i)
    colnames(localMax.i) <- colnames(wCoefs.i)
    ridgeList.i <- getRidge_J(localMax.i, gapTh = 3, skip = NULL, 
                            ...)
    ridgeName.i <- names(ridgeList.i)
    ridgeInfo.i <- matrix(as.numeric(unlist(strsplit(ridgeName.i, 
                                                     "_"))), nrow = 2)
    ridgeLevel.i <- ridgeInfo.i[1, ]
    newPeak.i <- ridgeInfo.i[2, ]
    selInd.i <- which.min(abs(newPeak.i - winSize.i/2))
    newRidgeLine.i <- ridgeList.i[[selInd.i]]
    newRidgeValue.i <- wCoefs.i[cbind(newRidgeLine.i, (1:length(newRidgeLine.i)) + 
                                        ridgeLevel.i[selInd.i] - 1)]
    peakScaleInd.i <- which.max(newRidgeValue.i)
    newPeakValue.i <- max(newRidgeValue.i)
    peakScale.new <- c(peakScale.new, scales.i[peakScaleInd.i])
    peakValue.new <- c(peakValue.new, newPeakValue.i)
    peakCenterIndex.new <- c(peakCenterIndex.new, newRidgeLine.i[peakScaleInd.i] + 
                               start.i - 1)
  }
  if (!is.null(peakSNR)) {
    peakSNR.new <- peakSNR * peakValue.new/peakValue
  }
  else peakSNR.new <- NULL
  names(peakScale.new) <- names(peakValue.new) <- names(peakCenterIndex.new) <- peakName[keep.i]
  unProcessedPeak <- peakName[unProcessedInd]
  return(list(peakIndex = peakIndex, peakValue = peakValue.new, 
              peakCenterIndex = peakCenterIndex.new, peakSNR = peakSNR.new, 
              peakScale = peakScale.new, unProcessedPeak = unProcessedPeak))
}

smth.gaussian_J <- function (x = stop("Numeric Vector 'x' is Required"), window = getOption("smoother.window"), 
                             alpha = getOption("smoother.gaussianwindow.alpha"), ..., 
                             tails = getOption("smoother.tails")) 
{
  if (!is.numeric(x) | !is.numeric(alpha)) {
    stop("argument 'x' and 'alpha' must be numeric", call. = FALSE)
  }
  windowLength = smoother:::.determineWindowLength(x, window)
  makeWindow = function(w, a) {
    hw = abs(w/2)
    e = exp(1)
    a = abs(a)
    ret = sapply(c(0:(w - 1)), function(x) {
      n = x - as.integer(hw)
      k = -0.5 * (a * n/hw)^2
      e^k
    })
    ret
  }
  w = makeWindow(windowLength, alpha[1])
  sizeW = length(w)
  sizeD = length(x)
  w = smoother:::.normalize(w)
  hkwL = as.integer(sizeW/2)
  hkwR = sizeW - hkwL
  pb <- pbapply::startpb(min = 1, max = sizeD)
  ret = sapply(c(1:sizeD), function(i) {
    pbapply::setpb(pb, i)
    ix.d = c((i - hkwL):(i + hkwR - 1))
    ix.w = which(ix.d %in% 1:sizeD)
    ix.d = ix.d[ix.w]
    W.nm = smoother:::.ifthenelse(length(ix.w) != sizeW, smoother:::.normalize(w[ix.w]), 
                       w)
    D.nm = x[ix.d]
    as.numeric(D.nm %*% W.nm)
  })
  if (!tails) {
    ret[c(1:hkwL, (sizeD - hkwR + 1):sizeD)] = NA
  }
  ret
}

