# ---- source scripts ----

outdir = normalizePath("~/Documents/umc/data/Data/BrSpIt/MZXML")
scriptdir = normalizePath("~/Google Drive/MetaboShiny/miniapps/MassChecker/scripts")
resdir = file.path(outdir, "results")

for(f in list.files(file.path(scriptdir), full.names = T)) source(f)
for(f in list.files(file.path(scriptdir, "AddOnFunctions"), full.names = T)) source(f)

# ---- get internal standards ----

int.db = "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/internal.full.db"

library(RSQLite)
library(data.table)

conn <- RSQLite::dbConnect(RSQLite::SQLite(), int.db)

query = "SELECT DISTINCT b.compoundname, b.baseformula, fullmz, adduct, foundinmode FROM extended e JOIN base b ON e.baseformula = b.baseformula WHERE b.compoundname LIKE '%(IS)%' AND e.isoprevalence > 99"

standards <- as.data.table(RSQLite::dbGetQuery(conn, query))

standards_pos_base <- standards[foundinmode == "positive"]
standards_neg_base <- standards[foundinmode == "negative"]

standards_pos_base <- standards_pos[adduct == "M+H"]
standards_neg_base <- standards_pos[adduct == "M-H"]

# ---- load in files ----

avg_files_pos <- list.files(file.path(resdir, "averaged"), pattern = "\\_pos\\.RData", full.names = T)
avg_files_neg <- list.files(file.path(resdir, "averaged"), pattern = "\\_neg\\.RData", full.names = T)

curr_pos <- avg_files_pos[[1]]
curr_neg <- avg_files_neg[[1]]

### SELECT CURRENT MODE FOR TESTING ###

current <- curr_pos
ppm = 1

# --- peak call these files ---

load(current)

vec <- averaged[[1]]@intensity
names(vec) <- averaged[[1]]@mass

zeros = which(vec == 0)

for(i in 1:length(zeros)){
  curr = vec[i]
  before = vec[i-1]
  after = vec[i + 1]
  if(before > 0 && after > 0) vec[i] = (before + after)/2
}

averaged[[1]]@intensity <- vec

### WAVELET ###

df <- matrix(data = averaged[[1]]@intensity,
             ncol=1)
rownames(df) <- averaged[[1]]@mass


plot(rownames(df)[rng],df[rng,])
lines(rownames(df)[rng],df[rng,])

scales <- seq(1, 64, 2)
sig_noise=6e-7
wCoefs <- MassSpecWavelet::cwt(df, scales=scales, wavelet='mexh')
localMax <- MassSpecWavelet::getLocalMaximumCWT(wCoefs,amp.Th = 500)

res <- MassSpecWavelet::plotLocalMax(localMax, 
                                     wCoefs, 
                                     range = c(min(rng), max(rng)))
abline(h = 10)


# --- find relationship between scale and mz ---

known_peaks <- pastecs::turnpoints(df[rng,])
npeaks <- length(which(known_peaks$peaks))

range = 10
for(peak in known_peaks){
  plot(rownames(df)[peak-range:peak+range],df[peak-range:peak+range,])
  lines(rownames(df)[peak-range:peak+range],df[peak-range:peak+range,])
  Sys.sleep(2)
}
# ----------------------------------------------

ridgeList <- getRidge(localMax,minWinSize = 5)
MassSpecWavelet::plotRidgeList(ridgeList, wCoefs, range = c(min(rng), max(rng)))

majorPeakInfo <- identifyMajorPeaks(df, 
                                    ridgeList, 
                                    wCoefs, 
                                    SNR.Th = sig_noise, 
                                    peakScaleRange = 1,
                                    nearbyPeak = T)

peaks_wav = MALDIquant::createMassPeaks(mass = as.numeric(rownames(df)[majorPeakInfo$peakCenterIndex]),
                                        intensity = majorPeakInfo$peakValue,
                                        snr = majorPeakInfo$peakSNR)

csv <- fread("/Users/jwolthuis/Documents/umc/data/Data/BrSpIt/MZXML/results/specpks_grouped_mdq/grouped_pos.csv")
csv[is.na(csv)] <- 0
mycol <- grepl(colnames(csv), pattern = sampname)
peaks_wav = MALDIquant::createMassPeaks(mass = as.numeric(csv$mzmed),
                                        intensity = csv[,..mycol][[1]])

match_wav <- pbapply::pblapply(1:nrow(standards_pos_base), FUN=function(i){
  std <- standards_pos_base[i,]
  peak <- std$fullmz
  matches <- which(peaks_wav@mass %between% c(peak - peak/(ppm*1e6), peak + peak/(ppm*1e6)))
  if(length(matches) > 0){
    paste_row <- data.table(match_mz = peaks_wav@mass[matches], match_int = peaks_wav@intensity[matches])
    res <- cbind(std, paste_row)
  }else{
    paste_row <- data.table(match_mz = NA, match_int = NA)
    res <- cbind(std, paste_row)
  }
  res
})

res_wav <- rbindlist(match_wav[!is.null(match_wav)])
res_wav

### MALDIQUANT ###

peaks_mald <- MALDIquant::detectPeaks(averaged[[1]], 
                                      snr=sig.noise)

match_mald <- pbapply::pblapply(1:nrow(standards_pos_base), FUN=function(i){
  std <- standards_pos_base[i,]
  peak <- std$fullmz
  matches <- which(peaks_mald@mass %between% c(peak - peak/(ppm*1e6), 
                                               peak + peak/(ppm*1e6)))
  if(length(matches) > 0){
    paste_row <- data.table(match_mz = peaks_mald@mass[matches], 
                            match_int = peaks_mald@intensity[matches])
    res <- cbind(std, paste_row)
  }else{
    paste_row <- data.table(match_mz = NA, 
                            match_int = NA)
    res <- cbind(std, paste_row)
  }
  res
})

res_mald <- rbindlist(match_mald[!is.null(match_mald)])
res_mald

### GAUSSIAN ###

values = list("mean"=NULL,"area"=NULL,"nr"=NULL,"min"=NULL,"max"=NULL,"qual"=NULL,"spikes"=0)
int.factor=1*10^5
factor=FALSE
resol=140000
sampname=names(averaged)
plot=FALSE
thresh=2000
scale = 1
scanmode = "pos"
peaks = searchMZRange(df[,1],values,int.factor,scale,resol,outdir,sampname,scanmode,FALSE,100,100,thresh)  
validPeaks = which(peaks$qual > 0)
mzs <- as.numeric(peaks$mean[validPeaks])
ints <- as.numeric(peaks$area[validPeaks])

peaks_gauss = MALDIquant::createMassPeaks(mass = as.numeric(peaks$mean[validPeaks]),
                                          intensity = peaks$area[validPeaks])

match_gauss <- pbapply::pblapply(1:nrow(standards_pos_base), FUN=function(i){
  std <- standards_pos_base[i,]
  peak <- std$fullmz
  matches <- which(peaks_gauss@mass %between% c(peak - peak/(ppm*1e6), peak + peak/(ppm*1e6)))
  if(length(matches) > 0){
    paste_row <- data.table(match_mz = peaks_gauss@mass[matches], 
                            match_int = peaks_gauss@intensity[matches])
    res <- cbind(std, paste_row)
  }else{
    paste_row <- data.table(match_mz = NA, 
                            match_int = NA)
    res <- cbind(std, paste_row)
  }
  res
})

res_gauss <- rbindlist(match_gauss[!is.null(match_gauss)])
res_gauss

# --- PROcess method ---
library(PROcess)

getPeaksJ <- function (bseoffM, SoN = 2, span = 81, sm.span = 11, 
                       zerothrsh = 2, area.w = 0.003, ratio = 0.2) 
{
  mzs <- as.numeric(rownames(bseoffM))
  n <- dim(bseoffM)[2]
  Spec <- colnames(bseoffM)
  cnts <- rep(0, n)
  for (j in 1:n) {

    bseoff <- cbind(mzs, bseoffM[, j])
    pks <- PROcess::isPeak(bseoff, SoN = SoN, span = span, sm.span = sm.span, 
                  zerothrsh = zerothrsh, area.w = area.w, ratio = ratio)
    print(pks)
    cnts[j] <- sum(pks$peak)
    is.peak <- pks$peak
    if (j > 1) {
      Peak. <- c(Peak., cumsum(is.peak[is.peak]))
      Intensity <- c(Intensity, pks$smooth[is.peak])
      Substance.Mass <- c(Substance.Mass, pks$mz[is.peak])
    }
    else {
      Peak. <- cumsum(is.peak[is.peak])
      Intensity <- pks$smooth[is.peak]
      Substance.Mass <- pks$mz[is.peak]
    }
  }
  pks
}

peaks_proc <- getPeaksJ(df, 0, 81, 11, 2, 0.003, 0.0000001)
is.peak <- peaks_proc$peak

peaks_proc = MALDIquant::createMassPeaks(mass = as.numeric(peaks_proc$mz[is.peak]),
                                          intensity = peaks_proc$smooth[is.peak])
match_proc <- pbapply::pblapply(1:nrow(standards_pos_base), FUN=function(i){
  std <- standards_pos_base[i,]
  peak <- std$fullmz
  matches <- which(peaks_proc@mass %between% c(peak - peak/(ppm*1e6), peak + peak/(ppm*1e6)))
  if(length(matches) > 0){
    paste_row <- data.table(match_mz = peaks_proc@mass[matches], match_int = peaks_proc@intensity[matches])
    res <- cbind(std, paste_row)
  }else{
    paste_row <- data.table(match_mz = NA, match_int = NA)
    res <- cbind(std, paste_row)
  }
  res
})

res_proc <- rbindlist(match_proc[!is.null(match_proc)])

# === SIMPL METHOD ===

vec <- averaged[[1]]@intensity
names(vec) <- averaged[[1]]@mass

zeros = which(vec == 0)

for(i in 1:length(zeros)){
  curr = vec[i]
  before = vec[i-1]
  after = vec[i + 1]
  if(before > 0 && after > 0) vec[i] = (before + after)/2
}

averaged[[1]]@intensity <- vec

peaks <- pastecs::turnpoints(vec)

# -----------

dips <- which(peaks$pits)

peakinfo <- pbapply::pblapply(1:length(dips), FUN=function(i){
  print(i)
  if(i == 1){
    segm <- peaks$points[1:dips[i]]
    names(segm) = names(vec)[1:dips[i]]
    
  }else if(i < length(dips)){
    segm <- peaks$points[dips[i]:dips[i+1]]
    names(segm) = names(vec)[dips[i]:dips[i+1]]
    
  }else{
    NULL
  }
  segm
})

pblapply(1:length(peakinfo), FUN=function(i){
  peak = peakinfo[[i]]
  #plot(peak)
  peak_hist = unlist(sapply(1:length(peak), FUN = function(j){
    counts = peak[j]
    mz = as.numeric(names(peak))[j]
    print(mz)
    # --------------------
    rep(x = mz, times = counts)
  }))
  peaks <- pastecs
  #print(peak_hist)
  fit <- MASS::fitdistr(peak_hist, 
                        densfun = "normal")
  # -------------
  fit
})

# -----------
ncomp = length(which(peaks$peaks))
mu <- as.numeric(names(vec)[which(peaks$peaks)])
fwhm = as.numeric(sapply(mu, getFwhm, resol = resol))
sigma = 2.35482004503 * fwhm

peak_info <- data.table(mz = mu,
                        sigma = sigma,
                        height = peaks$points[which(peaks$peaks)])

peak_info$area <- pbapply::pbsapply(1:nrow(peak_info), FUN=function(i){
  row = peak_info[i,]
  A = row$height*row$sigma/0.3989 
  A
})

peaks_tp = MALDIquant::createMassPeaks(mass = peak_info$mz,
                                       intensity = peak_info$area)

match_tp <- pbapply::pblapply(1:nrow(standards_pos_base), FUN=function(i){
  std <- standards_pos_base[i,]
  peak <- std$fullmz
  matches <- which(peaks_tp@mass %between% c(peak - peak/(ppm*1e6), peak + peak/(ppm*1e6)))
  if(length(matches) > 0){
    paste_row <- data.table(match_mz = peaks_tp@mass[matches], match_int = peaks_tp@intensity[matches])
    res <- cbind(std, paste_row)
  }else{
    paste_row <- data.table(match_mz = NA, match_int = NA)
    res <- cbind(std, paste_row)
  }
  res
})

res_tp <- rbindlist(match_tp[!is.null(match_tp)])

# ----------

# === MIXTOOLS ===

library(mixtools)

# ------ my own peak finder !!! ------


vec <- averaged[[1]]@intensity
names(vec) <- averaged[[1]]@mass

zeros = which(vec == 0)

for(i in 1:length(zeros)){
  curr = vec[i]
  before = vec[i-1]
  after = vec[i + 1]
  if(before > 0 && after > 0) vec[i] = (before + after)/2
}
#smoothed <- MALDIquant::smoothIntensity(averaged, halfWindowSize=5)
#vec <- smoothed[[1]]@intensity
#names(vec) <- smoothed[[1]]@mass


rl <- rle(vec == 0)
i1 <- rl$lengths>1 & rl$values
lst <- split(vec, rep(cumsum(c(TRUE, i1[-length(i1)])), rl$lengths)) 
split_peaks = lapply(lst, function(x) x[which(x!=0)])
#split_peaks = split_peaks[-sapply(split_peaks, is.null)] 
split_peaks = split_peaks[2: length(split_peaks)]

resol = 140000

peak_info <- pbapply::pblapply(split_peaks, FUN=function(peak){
  if(length(peak) > 0){
    tp <- pastecs::turnpoints(peak)
    ncomp = length(which(tp$peaks))
    mu <- as.numeric(names(peak)[tp$tppos])
  }else{
    ncomp = 0
    mu=0
  } 
  fwhm = as.numeric(sapply(mu, getFwhm, resol = resol))
  print(fwhm)
  sigma = 2.35482004503 * fwhm
  # ------
  list(num = ncomp, mu = mu, sigma = sigma)
})


peak_num <- sapply(peak_info, function(x) x$num)
peak_mu <- lapply(peak_info, function(x) x$mu)
peak_sigma <- lapply(peak_info, function(x) x$sigma)


# --------------------------

keep <- which(peak_num > 1)
multipeaks <- split_peaks[keep]
ncomp_list <- peak_num[keep]
mu_list <- peak_mu[keep]
sigma_list <- peak_sigma[keep]

cl = parallel::makeCluster(3, "FORK")

peaks_mix <- pbapply::pblapply(1:length(multipeaks), cl=0, FUN=function(i){
  ncomp <- as.integer(ncomp_list[[i]])
  mu <- mu_list[[i]]
  sigma <- sigma_list[[i]]
  #print(mu)
  peak <- multipeaks[[i]]

  dat = unlist(sapply(1:length(peak), FUN = function(j){
    counts = peak[j]
    mz = as.numeric(names(peak[j]))
    # --------------------
    rep(x = mz, times = counts / 10)
  }))
  #print(ncomp)
  #print(dat)
  # if(ncomp == 2){
    mixture <- mixtools::normalmixEM(dat,
                                     k=ncomp,
                                     fast = TRUE,
                                     mu = mu,
                                     sigma = sigma,
                                     maxit = 1e4
                                     #mean.constr = mu,
                                     #sd.constr = sigma
                                     )
    #mu = mixture$mu
    #lambda = mixture$lambda
  # }else{
    # plot(names(peak),peak)
    # lines(names(peak),peak)
    # abline(v = mu)
    # mixture = nor1mix::norMixMLE(dat, ncomp)
    # lambda = mixture[(2 * length(mu) + 1):(2 * length(mu) + length(mu))]
    # mu = mixture[1:length(mu)]
  #}
  
    # --- return ---
    
    mixture
  #mixtools::plot.mixEM(mixture,whichplots = 2)
  # peak_areas <- lambda * length(dat)
  # names(peak_areas) <- mu
  # print(peak_areas)
  # # --- return ---
  # peak_areas
})


# --- compare all  ---

perc_wav = 100 - (length(which(is.na(res_wav$match_mz))) / nrow(res_wav) * 100.00)
perc_mald = 100 - (length(which(is.na(res_mald$match_mz))) / nrow(res_mald) * 100.00)
perc_gauss = 100 - (length(which(is.na(res_gauss$match_mz))) / nrow(res_gauss) * 100.00)
perc_proc = 100 - (length(which(is.na(res_proc$match_mz))) / nrow(res_proc) * 100.00)
perc_tp = 100 - (length(which(is.na(res_tp$match_mz))) / nrow(res_tp) * 100.00)


{
  cat(paste0("Wavelets method: ", round(perc_wav,1), "% identified, ", round(range(res_wav$match_int, na.rm = T)[2] - range(res_wav$match_int, na.rm = T)[1],2), " intensity range \n"))
  #boxplot(res_wav$match_int)
  cat(paste0("Centroid method: ", round(perc_mald,1), "% identified, ", round(range(res_mald$match_int, na.rm = T)[2] - range(res_mald$match_int, na.rm = T)[1],2), " intensity range \n"))
  #boxplot(res_mald$match_int)
  cat(paste0("Gaussian method: ", round(perc_gauss,1), "% identified, ", round(range(res_gauss$match_int, na.rm = T)[2] - range(res_gauss$match_int, na.rm = T)[1],2), " intensity range \n"))
  #boxplot(res_gauss$match_int)
  #cat(paste0("PROCess method:  ", round(perc_proc,1), "% identified,  ", round(range(res_proc$match_int, na.rm = T)[2] - range(res_proc$match_int, na.rm = T)[1],2), " intensity range \n"))
  #cat(paste0("Turnpoints method:  ", round(perc_tp,1), "% identified,  ", round(range(res_tp$match_int, na.rm = T)[2] - range(res_tp$match_mix, na.rm = T)[1],2), " intensity range \n"))
}
