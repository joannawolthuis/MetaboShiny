library(xcms)
library(data.table)
library(pbapply)

# --- get internal standards ---

library(RSQLite)
library(DBI)

conn <- dbConnect(RSQLite::SQLite(), "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/internal.full.db")

IS_table <- dbGetQuery(conn, "SELECT b.compoundname, e.adduct, e.fullmz, e.isoprevalence, e.foundinmode FROM base b
                       JOIN extended e ON b.baseformula = e.baseformula
                       WHERE b.compoundname like '%(IS)%'")

unique(IS_table[,c("adduct", "foundinmode")])
IS_dt <- as.data.table(IS_table)
calib_pos <- IS_dt[foundinmode == "positive"]
calib_neg <- IS_dt[foundinmode == "negative"]
calibrants = list(positive = calib_pos,
                  negative = calib_neg)
save(file = "calibrants.RData", calibrants)

any(is.na(calib_pos$fullmz))
# --------------------------------

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

examplefiles <- list.files(system.file("fticr", package = "msdata"),
                           recursive = TRUE, full.names = TRUE)

scanrange_pos = 20:30

# ATTEMPT 1: PEAK CALL EACH SCAN, GROUP X 2 # WITH CALIBRATION FROM INTERNAL STANDARDS

xsets_pos_1 <- lapply(scanrange_pos, FUN=function(scan){
  try(xcmsSet(method="MSW",
              snames = paste0("scan", scan),
              files=myfiles[1],
              scanrange = c(scan,scan),
              snthresh=3,
              nearbyPeak=T))
})

xset_pos_1 <- c(xsets_pos[[1]])

for(xset in xsets_pos_1[2:length(xsets_pos_1)]){
  print(xset)
  xset_pos_1 <- c(xset_pos_1, xset)
}

xset_pos_1

xsets_pos_2 <- lapply(scanrange_pos, FUN=function(scan){
  try(xcmsSet(method="MSW",
              snames = paste0("scan", scan),
              files=myfiles[2],
              scanrange = c(scan,scan),
              snthresh=3,
              nearbyPeak=T))
})

xset_pos_2 <- c(xsets_pos_2[[1]])

for(xset in xsets_pos_2[2:length(xsets_pos_2)]){
  print(xset)
  xset_pos_2 <- c(xset_pos_2, xset)
}

xset_pos_2

# --- get internal standards ---

library(RSQLite)
library(DBI)

conn <- dbConnect(RSQLite::SQLite(), "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/internal.full.db")
adducts = c("M+H", "M-H")
IS_table <- dbGetQuery(conn, "SELECT b.compoundname, e.adduct, e.fullmz, e.isoprevalence FROM base b
                              JOIN extended e ON b.baseformula = e.baseformula
                              WHERE b.compoundname like '%(IS)%'")

IS_dt <- as.data.table(IS_table)
picked_is <- IS_dt[adduct %in% adducts ]#& isoprevalence > 99.999]
picked_is_pos <- picked_is[adduct == "M+H"]

# ------------------------------

xset_pos1_calib <- calibrate(xset_pos_1, 
                            calibrants = picked_is_pos$fullmz, 
                            mzppm = 2, 
                            method="edgeshift", 
                            plotres=T)
grouped1_scans <- group(xset_pos1_calib, 
                       method="mzClust",mzppm=2)

grouped1_scans@groups

grouped1_mat_pscan <- groupval.joanna(grouped1_scans@peaks,
                                      grouped1_scans@groups,
                                      grouped1_scans@groupidx,
                                      sampnames(grouped1_scans),
                                      method = "maxint",
                                      value = "index",
                                      intensity = "into")

xset_pos2_calib <- calibrate(xset_pos_2, 
                             calibrants = picked_is_pos$fullmz, 
                             mzppm = 2, 
                             method="edgeshift", 
                             plotres=T)
grouped2_scans <- group(xset_pos2_calib, 
                        method="mzClust",mzppm=2)

grouped2_scans@groups

grouped2_mat_pscan <- groupval.joanna(grouped2_scans@peaks,
                                grouped2_scans@groups,
                                grouped2_scans@groupidx,
                                sampnames(grouped2_scans),
                                method = "maxint",
                                value = "index",
                                intensity = "into")


grouped1_mat_psamp <- rbindlist(lapply(1:nrow(grouped1_mat_pscan), FUN=function(i){
  row = grouped1_mat_pscan[i,]
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

grouped1_mat_psamp

grouped2_mat_psamp <- rbindlist(lapply(1:nrow(grouped2_mat_pscan), FUN=function(i){
  row = grouped2_mat_pscan[i,]
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

grouped2_mat_psamp

sampints <- list(grouped1_mat_psamp,
                 grouped2_mat_psamp)

# === GROUP SAMPELS TOGETHER TO GET SUPERGROUPS ===

allgroups <- rbind(cbind(grouped1_scans@groups, 
                         sample=c(1)),
                   cbind(grouped2_scans@groups,
                         sample=c(2)))

allgroups <- as.data.table(allgroups[order(allgroups[,1]),])

supergrouped <- mzClustGeneric(allgroups)

supergrouped$mat

peakintensities <- pblapply(1:nrow(supergrouped$mat), FUN=function(i){
  # --- get mz boundaries of supergroup ---
  sg_row <- supergrouped$mat[i,]
  sg_mem <- supergrouped$idx[i]
  samps <- unique(allgroups$sample)
  # PER SAMPLE, get info for this supergroup
  sampcols <- lapply(samps, FUN=function(samp){
    # go through members of supergroup
    for(g1 in sg_mem){
      # get mz boundaries for this supergroup member group
      g_rows <- allgroups[g1,]
      # check if current sample has a value here
      if(samp %in% g_rows$sample){
        int_total = 0
        # go through groups present in this sample and get intensity
        for(i in 1:nrow(g_rows)){
          row = g_rows[i,]
          # find intensity in saved group tabels per sample
          if(row$sample == samp){
            int <- as.data.table(sampints[samp])[mzmed == row$mzmed, "into"]
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
  row <- cbind(mzmed = sg_row[1],
               mzmin = sg_row[2], 
               mzmax = sg_row[3], 
               row)
  # --- return --- 
  row
  # --------------
})

sg_intvals <- rbindlist(peakintensities)
sg_intvals

# ============================================

# xcmsRaw > xcmsSet???

xr = xcmsRaw(myfiles[1],profstep = 0.01)
plot(xr@tic)
avgspec <- getSpec(xr)
xr_peak <- findPeaks.MSW(xr, snthresh=3, 
                    scales=seq(1,22,3), 
                    nearbyPeak=T, 
                    peakScaleRange=5, 
                    amp.Th=0.01, 
                    ridgeLength=24, 
                    tuneIn=FALSE, 
                    verbose.columns = FALSE)

any(avgspec[,2])


# ============================================

# ATTEMPT 2: AVERAGE OR SUM SPECTRA TO ONE AND USE LOW LEVEL FUNCS

######### function to sum the spectra of two xcmsRaw objects together and return it as a new xcmsRaw object 
getSumSpec = function (object, ...) 
{
  sel <- profRange(object, ...)
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
  # --------------
  points <- cbind(mz = uniquemz, intensity = rowSums(intmat,na.rm = T))
  print(head(points))
  points[order(points[,1]),]
}

# ALL IN ONE FILE 

sumPeakXraw = function(file, profstep=0.01, scales=c(1, seq(2, 30, 2), seq(32, 64, 4)), snr=3, mode="positive"){
  xraw_full = xcmsRaw(file, profstep=profstep)
  scans = which(xraw_full@polarity == mode)
  print(scans)
  xraw = xcmsRaw(file,profstep = profstep, scanrange = c(min(scans), 
                                                         max(scans)))
  print(xraw)
  pos_profile <- xraw@env$profile
  
  p <- pos_profile[,1]
  for(i in 2:ncol(pos_profile)){
    vec <- pos_profile[,i]
    p = p + vec
  }
  pmat <- matrix(p)
  
  scansum <- as.data.frame(getSumSpec(xraw))
  
  print("here")
  xraw@env$mz <- scansum$mz
  xraw@env$intensity <- scansum$intensity
  xraw@env$profile <- pmat
  xraw@scanindex <- as.integer(0)
  xraw@scantime <- xraw@scantime[1]
  xraw@tic <- sum(scansum$intensity)
  xraw@polarity <- as.factor(mode)
  xraw@acquisitionNum <- as.integer(1)
  xraw@scanrange <- 1

  xraw  
}

summed_xraws <- lapply(myfiles[1:10], sumPeakXraw)
specpeaks <- lapply(summed_xraws, FUN=function(xraw){
  peaks <- findPeaks.MSW(xraw, snthresh=3,
                         scales=c(1, seq(2, 30, 2), seq(32, 64, 4)),
                         nearbyPeak=T)
  peaks
  })


specpeaks[2]
summed_xraws[1]
                    
names(spectra) <- gsub(basename(myfiles[1:10]), pattern="\\.mzXML", replacement = "")
getSumPeaks(myfiles[1], mode = "positive")

# === WRITE TO CDF... ===

write.cdf(xr_sum,file=file.path(folder, "test.cdf"))
fp = normalizePath(file.path(folder, "test.cdf"))

reloadset <- xcmsSet(files=fp,
                     method="MSW")
reloadset
xs = xcmsSet(myfiles[1],
             scanrange=c(10,150), 
             method="MSW")
xs
write.cdf(xr_ex, tempfile())

# --- test...
# ============================================

# ATTEMPT 3: FIND HIGHEST INTENSITY SCAN IN RANGE AND GO FROM THERE

xr <- xcmsRaw(myfiles[1],profstep = 0.01)

bestpos <- which.max(xr@tic[1:150])
bestneg <- 150 + which.max(xr@tic[151:300])

xset_pos <- xcmsSet(method="MSW",
        files=myfiles[1],
        scanrange = c(bestpos, bestpos),
        snthresh=3,
        nearbyPeak=T)

xset_neg <- xcmsSet(method="MSW",
                    files=myfiles[1],
                    scanrange = c(bestneg, bestneg),
                    snthresh=3,
                    nearbyPeak=T)

# get seperate xcmssets

# === PEAK CALL EACH SCAN, GROUP X 2 ===
cl = makeCluster(3, "FORK")
stopCluster(cl)

peakGroupScans <- function(file,
                           limit=NULL, 
                           cl=FALSE,
                           profstep=0.01, 
                           scales=c(1, seq(2, 30, 2), seq(32, 64, 4)), 
                           snr=3, 
                           mode="positive", 
                           group_ppm=2){
  
  xraw_full = xcmsRaw(file, profstep=profstep)
  scans = which(xraw_full@polarity == mode)[if(limit) c(1:limit)]

  peaks_p_scan <- pblapply(scans, cl=cl, FUN=function(scan){
    xraw = xcmsRaw(file,
                   profstep = profstep, 
                   scanrange = c(scan,scan))
    xraw_called = findPeaks.MSW(xraw, 
                                snthresh=snr, 
                                scales = scales, #seq(1, 64, 3),#seq(1,22,3), 
                                nearbyPeak=T)
    # --- return --- 
    xraw_called
  })
  
  rows <- lapply(seq_along(peaks_p_scan), FUN=function(i){
    peaks <- as.data.table(peaks_p_scan[[i]]@.Data)
    peak_subset <- peaks[,c("mz", "into")]
    print(peak_subset)
    peak_subset
  })
  
  print(head(rbindlist(rows)))
  uniquepeaks <- as.matrix(aggregate(into ~ mz, 
                                     data=rbindlist(rows), 
                                     FUN=sum))
  print(head(uniquepeaks))
  # GROUPS AND INTENSITIES FOR THIS SINGLE SAMPLE
  peakgroups <- mzClustGeneric(uniquepeaks, 
                               mzppm = group_ppm)
  
  # capture group members for supergroup creation later
  memberlist = lapply(seq_along(peakgroups$idx), 
                      FUN=function(i){
    idx <- peakgroups$idx[[i]]
    membermzs <- c(uniquepeaks[idx, 1])
    # -------
    membermzs
  })
  
  names(memberlist) <- c(peakgroups$mat[,1])
  
  # get this back to a list of mzs and intensities to group again later...
  rows <- lapply(seq_along(peakgroups$idx), 
                 FUN=function(i){
    groupmz <- peakgroups$mat[i,"mzmed"]
    idx <- peakgroups$idx[[i]]
    intsum <- sum(uniquepeaks[idx, "into"])
    # --- return ---
    data.table(mz = groupmz,
               into = intsum)
  })
  file_group_ints <- rbindlist(rows)
  
  # === RETURN ===
  list(mat = file_group_ints, 
       members = memberlist)
}

mode = "positive"
cl=FALSE
peaktables <- lapply(myfiles[1:3], 
                     FUN=function(file){
  peakGroupScans(file,
                 cl=cl, 
                 mode = mode, 
                 limit=5)
})

names(peaktables) <- gsub(basename(myfiles[1:3]), 
                          pattern="\\.mzXML", 
                          replacement = "")

peakoverview <- as.matrix(rbindlist(peaktables))

peakgroups <- mzClustGeneric(peakoverview)

rows <- lapply(seq_along(peakgroups$idx), FUN=function(i){
  groupmz <- peakgroups$mat[i,"mzmed"]
  grouprange <- peakgroups$mat[i,c(1:3)]
  print(grouprange)
  idx <- peakgroups$idx[[i]]
  mzmembers <- peakoverview[idx]
  sampcols <- lapply(seq_along(peaktables), FUN=function(i){
    peaktable = as.data.table(peaktables[[i]])
    sampname = names(peaktables)[i]
    int <- peaktable[mz %in% mzmembers, "into"]
    if(nrow(int) == 0) int <- data.table(NA)
    colnames(int) <- sampname
    # --- return ---
    int
  })
  # --- return ---
  sampints <- cbind(mzmed = grouprange[1], 
                    mzmin = grouprange[2], 
                    mzmax = grouprange[3], do.call(cbind, sampcols))
  print(sampints)
  sampints
})

mean(x) + rnorm(length(missing(x))) * sd(x)

data <- as.data.frame(rbindlist(rows))

# FILL PEAKS

for(i in 4:ncol(data)){
  sampname = colnames(data)[i]
  print(sampname)
  missing = which(is.na(data[,i]))
  # from ORIGINAL LIST - get values in range and integrate...
  fillvals <- sapply(missing, FUN=function(i){
    mzrange = c(data[i, c("mzmin", "mzmax")])
    print(mzrange)
    xraw = xcmsRaw(myfiles[i], profstep=profstep)
    print(xraw@env$mz)
    idx <- which(xraw@env$mz >= mzrange[1] && xraw@env$mz <= mzrange[2])
    print(idx)
  })
}

head(data)

# dct <- list(
#   TEST1 = myfiles[1],
#   TEST2 = myfiles[2]
# )
# getFilledPeakTable <- function(peaktable, filedict){
#   sampidx = which(!grepl(colnames(peaktable), pattern = "mz"))
#   sampcols = colnames(peaktable)[sampidx]
#   # ------------------------------------
#   for(samp in sampcols){
#     file <- filedict[[samp]]
#     xraw <- xcmsRaw(file, profstep = 0.001, scanrange = c(min(scanrange_pos), max(scanrange_pos)))
#     fullspec <- getSpec(xraw)
#     # --- find missing mzval ranges ---
#     na_idx = which(is.na(peaktable[,..samp]))
#     filled_ints <- sapply(na_idx, FUN=function(idx){
#       mzmin = as.numeric(peaktable[idx,"mzmin"])
#       mzmax = as.numeric(peaktable[idx,"mzmax"])
#       # --- find ints in range ---
#       origidx <- which(fullspec[,1] >= (mzmin - mzmin * 2e-6) & fullspec[,1] <= (mzmax +  mzmax * 2e-6))
#       int <- sum(fullspec[origidx, 2],na.rm = T)
#       print(int)
#     })
#     print(na_idx)
#   }
# } # this doesnt work for now...