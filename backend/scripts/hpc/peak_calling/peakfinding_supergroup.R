library(xcms)
library(data.table)
library(pbapply)

load(file = "calibrants.RData")


groupval.joanna <- function(peakmat,
                            groupmat,
                            groupidx,
                            sampnames,
                            method = c("medret", "maxint"), 
                            value = "index", intensity = "into"){
  if (length(groupidx) < 
      1) {
    stop("xcmsSet is not been group()ed.")
  }
  groupindex <- groupidx
  sampnum <- seq(length = length(sampnames))
  retcol <- match("rt", colnames(peakmat))
  intcol <- match(intensity, colnames(peakmat))
  sampcol <- match("sample", colnames(peakmat))
  values <- matrix(nrow = length(groupindex), ncol = length(sampnum))
  if (method == "medret") {
    for (i in seq(along = groupindex)) {
      gidx <- groupindex[[i]][order(abs(peakmat[groupindex[[i]], 
                                                retcol] - median(peakmat[groupindex[[i]], 
                                                                         retcol])))]
      values[i, ] <- gidx[match(sampnum, peakmat[gidx, 
                                                 sampcol])]
    }
  }
  else {
    for (i in seq(along = groupindex)) {
      gidx <- groupindex[[i]][order(peakmat[groupindex[[i]], 
                                            intcol], decreasing = TRUE)]
      values[i, ] <- gidx[match(sampnum, peakmat[gidx, 
                                                 sampcol])]
    }
  }
  if (value != "index") {
    values <- peakmat[values, value]
    dim(values) <- c(length(groupindex), length(sampnum))
  }
  colnames(values) <- sampnames
  values_full <- cbind(groupmat[, c("mzmin",
                                    "mzmed",
                                    "mzmax", "npeaks")], values)
  values_full
}


folder = "~/Downloads/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil/RES-2017-10-31_DSM DBS Brazil part 1/MZXML"
myfiles = list.files(folder, full.names = T)

# ------------------------------

# ATTEMPT 1: PEAK CALL EACH SCAN, GROUP X 2 # WITH CALIBRATION FROM INTERNAL STANDARDS

getScanset <- function(file, mode, calibrants){
  print(calibrants)
  print(mode)
  xraw_full = xcmsRaw(file, profstep=0.1)
  print(xraw_full@polarity)
  scans = which(xraw_full@polarity == mode)
  print(scans)
  #scans = 20:30
  # --- PEAK CALL EACH SCAN ---
  scansets <- lapply(scans, FUN=function(scan){
    print(paste("scan", scan))
    res <- NULL
    try({res <- xcmsSet(method="MSW",
                        scales=c(1, seq(2, 30, 2), seq(32, 64, 4)),
                        snames = paste0("scan", scan),
                        files=file,
                        scanrange = c(scan,scan),
                        snthresh=3,
                        nearbyPeak=T)
    try({res <- calibrate(res,
                          calibrants = calibrants,
                          mzppm = 3,
                          method="edgeshift")
        })
    
    })
    print(res)
    # --- return ---
    res
  })
  sampset = NULL
  for(xset in scansets){
    sampset <- if(typeof(sampset) != "S4") c(xset) else c(sampset, xset)
  }
  if(typeof(sampset) != "S4") return(NULL)
  # --- GROUP CALIBRATED SCANS ---
  sampset <- group(sampset,
                   method="mzClust",
                   mzppm=2)
  # --- return ---
  sampset
}



# == UP TO HERE SHOULD BE ARRAY JOB ===

# 'xcms_end' file

xr <- xcmsRaw("/Users/jwolthuis/Downloads/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil/RES-2017-10-31_DSM DBS Brazil part 1/MZXML/RES_20171031_009.mzXML")
rawMat(xr)

plot(xr@tic)

xset_pos_1 <- getScanset("/Users/jwolthuis/Downloads/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil/RES-2017-10-31_DSM DBS Brazil part 1/MZXML/RES_20171031_001.mzXML", 
                         "positive", 
                         calibrants = calibrants$positive[isoprevalence > 99.99]$fullmz)

xset_pos_1

xset_pos_2 <- getScanset("/Users/jwolthuis/Downloads/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil/RES-2017-10-31_DSM DBS Brazil part 1/MZXML/RES_20171031_002.mzXML", 
                         "positive", 
                         calibrants = calibrants$positive[isoprevalence > 99.99]$fullmz)


getSampInfo <- function(sampleXsets){
  sampnames = names(sampleXsets)
  lst <- lapply(1:length(sampleXsets), FUN=function(i){
    xset <- sampleXsets[[i]]
    mat_pscan <- groupval.joanna(xset@peaks,
                                 xset@groups,
                                 xset@groupidx,
                                 sampnames(xset),
                                 method = "maxint",
                                 value = "index",
                                 intensity = "into")
    mat_psamp <- rbindlist(lapply(1:nrow(mat_pscan), FUN=function(i){
      row = mat_pscan[i,]
      scancols <- which(grepl(names(row), pattern = "scan"))
      scansum <- sum(row[scancols], na.rm = T)
      # --- new row ---
      newrow <- data.table(
        mzmin = row[1],
        mzmed = row[2],
        mzmax = row[3],
        npeaks = row[4],
        into = scansum
      )
      # -- return ---
      newrow
    }))
    dt <- as.data.table(mat_psamp)
    dt$sample <-  c(sampnames[[i]])
    dt
    })
    # --- return ---
  names(lst) <- sampnames
  lst
}

samps <- list(xset_pos_1, xset_pos_2)
sampnames = c("TEST1", "TEST2")


# === GROUP SAMPELS TOGETHER TO GET SUPERGROUPS ===

getPeakTable <- function(samples){
  sampints <- getSampInfo(samples)
  # -------------------
  allgroups <- rbindlist(sampints)
  allgroups <- as.data.table(allgroups[order(allgroups[,1]),])
  supergrouped <- mzClustGeneric(as.matrix(allgroups[,-"sample"]), mzppm = 2)
  peakintensities <- pblapply(1:nrow(supergrouped$mat), FUN=function(i){
    # --- get mz boundaries of supergroup ---
    superg_row <- supergrouped$mat[i,]
    superg_mem <- supergrouped$idx[i]
    samps <- unique(allgroups$sample)
    # PER SAMPLE, get info for this supergroup
    sampcols <- lapply(samps, FUN=function(samp){
      # go through members of supergroup
      for(subg in superg_mem){
        # get mz boundaries for this supergroup member group
        subg_rows <- allgroups[subg,]
        # check if current sample has a value here
        if(samp %in% subg_rows$sample){
          int_total = 0
          # go through groups present in this sample and get intensity
          for(i in 1:nrow(subg_rows)){
            subg_row = subg_rows[i,]
            # find intensity in saved group tabels per sample
            if(subg_row$sample == samp){
              sampint <- sampints[[samp]]
              int <- sampints[[samp]][mzmed == subg_row$mzmed, "into"]
              int_total = int_total + int 
            }
          }
          # if not, return na
        }else{ int_total <- NA }
      }
      # return col
      dt <- data.table(A = int_total)
      colnames(dt) <- as.character(samp)
      dt
    })
    row <- do.call(cbind, sampcols)
    row <- cbind(mzmed = superg_row[1],
                 mzmin = superg_row[2], 
                 mzmax = superg_row[3], 
                 row)
    # --- return --- 
    row
    # --------------
  })
  peaktable <- rbindlist(peakintensities)
  # --- return :-) ---
  peaktable
}

peaktable <- getPeakTable(samps, sampnames)

# impuuuute? https://www.r-bloggers.com/imputing-missing-data-with-r-mice-package/

getFilledPeakTable <- function(peaktable){
  sampidx = which(!grepl(colnames(peaktable), pattern = "mz"))
  sampcols = colnames(peaktable)[sampidx]
  # ------------------------------------
  for(samp in sampcols){
    print(samp)
    # --- find missing mzval ranges ---
    nas = which(is.na(peaktable[,..samp]))
    x <- unlist(peaktable[!nas, ..samp])
    filler <- mean(x) + rnorm(length(nas)) * sd(x)
    peaktable[nas, samp] <- filler
  }
  peaktable
}

filledPeakTable <- getFilledPeakTable(peaktable)

