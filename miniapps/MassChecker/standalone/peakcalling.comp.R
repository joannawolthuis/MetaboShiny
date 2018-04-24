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

standards_pos <- standards[foundinmode == "positive"]
standards_neg <- standards[foundinmode == "negative"]

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

rng = c(500:600)

plot(rownames(df)[rng],df[rng,], xlab = "m/z", ylab="intensity",cex.lab=1.5)
lines(rownames(df)[rng],df[rng,])

scales <- seq(1, 10, 0.1)
sig_noise=6e-7
wCoefs <- MassSpecWavelet::cwt(df, scales=scales, wavelet='mexh')
localMax <- MassSpecWavelet::getLocalMaximumCWT(wCoefs,amp.Th = 500)

res <- MassSpecWavelet::plotLocalMax(localMax, 
                                     wCoefs, 
                                     range = c(min(rng), max(rng)))
abline(h = 10)


# --- find relationship between scale and mz ---

# known_peaks <- pastecs::turnpoints(df[rng,])
# npeaks <- length(which(known_peaks$peaks))
# 
# range = 10
# for(peak in known_peaks){
#   plot(rownames(df)[peak-range:peak+range],df[peak-range:peak+range,])
#   lines(rownames(df)[peak-range:peak+range],df[peak-range:peak+range,])
#   Sys.sleep(2)
# }
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

sampname = gsub(x = basename(current), pattern = "_(pos|neg)\\.RData", replacement = "")
csv <- fread("/Users/jwolthuis/Documents/umc/data/Data/BrSpIt/MZXML/results/specpks_grouped_wavelet/grouped_pos.csv")
csv[,(1:ncol(csv)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]

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


### GAUSS_ADJUSTED ###

vec <- averaged[[1]]@intensity
names(vec) <- averaged[[1]]@mass
rl <- rle(vec == 0)
i1 <- rl$lengths>2 & rl$values
lst <- split(vec, rep(cumsum(c(TRUE, i1[-length(i1)])), rl$lengths)) 
split_peaks = lapply(lst, function(x) x[which(x!=0)])

peak_num <- pbapply::pblapply(split_peaks, FUN=function(peak){
  if(length(peak) > 0){
    tp <- pastecs::turnpoints(peak)
    ncomp = length(which(tp$peaks))
  }else{
    ncomp = 0
  } 
  # ------
  ncomp
})

# ----------

tp_peaks <- pbapply::pblapply(split_peaks, cl=0, FUN=function(peak){
  nested = FALSE
  if(length(peak) > 0){
    tp_all <- pastecs::turnpoints(peak)
    xvals = peak
    mu = peak[tp_all$peak]
    #ncomp = length(which(tp$peaks))
    #peak_mz = as.numeric(names(peak)[which(tp$peaks)])
    if(length(which(tp_all$peaks)) > 1){
      print("multiple")
      nested = TRUE
      # split
      pits = which(tp_all$pits)
      peaks = peak[tp_all$peaks]
      split_split_peaks = lapply(1:length(pits), function(i){
        pit = pits[i]
        if(i == 1){
          xvals = peak[1:pit]
        }else if(i == length(pits)){
          xvals = peak[pit:pits[length(pits)]]
        }else{
          xvals = peak[pits[i-1]:pit]
        }
        mu = intersect(names(xvals), names(peaks))
        list(x = xvals, m = mu)
      })
      return(list(orig_xvals = peak, split_xvals = split_split_peaks))
      #print(tp_all)
    }else{
      print("single")
      return(list(orig_xvals = xvals, split_xvals = FALSE, m = mu))
    }
  }else{
    return(list())
  }
  # ------
  #list(peak, ncomp, peaks, dips)
})

# ---------------

int.factor=1*10^5 # Number of x used to calc area under Gaussian (is not analytic) 

peak_rows <- pbapply::pblapply(tp_peaks, cl=cl, function(res){
  #sig = NULL
  #try({
  if(length(res$split_xvals) > 1){
    rows <- lapply(res$split_xvals, function(peak){
      try({
        if(length(peak[['m']]) == 0){
          NULL
        }else{
          vals = fitG_2(names(res$orig_xvals),
                        res$orig_xvals,
                        mu = peak$m,
                        2,
                        TRUE,
                        lower = c(min(as.numeric(names(peak$x))), .Machine$double.xmin),
                        upper = c(max(as.numeric(names(peak$x))), .Machine$double.xmax)
          )
          fwhm = getFwhm(vals$par[1], 140000)
          sig = 0.42466090 * fwhm
          area = getArea(mu = vals$par[1],
                         resol = 140000,
                         scale = vals$par[2],
                         int.factor = int.factor,
                         sigma = sig)
          qual = getFitQuality(x = as.numeric(names(peak)),
                               y = as.numeric(peak$x),
                               muFirst = peak$m,
                               muLast = vals$par[1],
                               resol = 140000,
                               scale = vals$par[2],
                               sigma = sig,
                               sumFit = vals$value)
          data.table(mz = vals$par[1],
                     intensity = area,
                     qual = qual$fq_new)
        }
      })
    })
  }else{
    if(length(res$m) == 0){
      rows = list()
    }else{
      try({
        vals = fitG_2(names(res$orig_xvals),
                      res$orig_xvals,
                      mu = res$m,
                      2,
                      TRUE,
                      lower = c(min(as.numeric(names(res$orig_xvals))), .Machine$double.xmin),
                      upper = c(max(as.numeric(names(res$orig_xvals))), .Machine$double.xmax)
        )
        fwhm = getFwhm(vals$par[1], 140000)
        sig = 0.42466090 * fwhm
        area = getArea(mu = vals$par[1],
                       resol = 140000,
                       scale = vals$par[2],
                       int.factor = int.factor,
                       sigma = sig)
        qual = getFitQuality(x = as.numeric(names(res$orig_xvals)),
                             y = as.numeric(res$orig_xvals),
                             muFirst = res$m,
                             muLast = vals$par[1],
                             resol = 140000,
                             scale = vals$par[2],
                             sigma = sig,
                             sumFit = vals$value)
        tbl <- data.table(mz = vals$par[1],
                          intensity = area,
                          qual = qual$fq_new)
        rows = list(tbl) 
      })
    }
  }
  if(length(rows) > 0){
    tbl <- rbindlist(rows[sapply(rows, 
                                 is.data.table)])
    tbl    
  }
  
  #})
})

tbl = rbindlist(peak_rows[sapply(peak_rows, is.data.table)])

sampname = gsub(x = basename(current), pattern = "_(pos|neg)\\.RData", replacement = "")
csv <- fread("/Users/jwolthuis/Documents/umc/data/Data/BrSpIt/MZXML/results/specpks_grouped_mdq_gauss/grouped_pos.csv")
csv[,(1:ncol(csv)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]

mycol <- grepl(colnames(csv), pattern = sampname)
peaks_gauss = MALDIquant::createMassPeaks(mass = as.numeric(csv$mzmed),
                                        intensity = csv[,..mycol][[1]])


peaks_gauss = MALDIquant::createMassPeaks(mass = tbl$mz,
                                          intensity = tbl$intensity,
                                          metaData = list(qual=tbl$qual))

match_gauss_2 <- pbapply::pblapply(1:nrow(standards_pos_base), FUN=function(i){
  std <- standards_pos_base[i,]
  peak <- std$fullmz
  matches <- which(peaks_gauss@mass %between% c(peak - peak/(ppm*1e6), 
                                               peak + peak/(ppm*1e6)))
  print(matches)
  if(length(matches) > 0){
    paste_row <- data.table(match_mz = paste(peaks_gauss@mass[matches], collapse = ";"), 
                            match_int = paste(peaks_gauss@intensity[matches],collapse=";"))
    res <- cbind(std, paste_row)
  }else{
    paste_row <- data.table(match_mz = NA, 
                            match_int = NA)
    res <- cbind(std, paste_row)
  }
  res
})

res_gauss_2 <- rbindlist(match_gauss_2[!is.null(match_gauss_2)])
res_gauss_2

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

# --- proFIAT ---

sname = gsub(basename(current), pattern = "_(pos|neg)\\.RData", replacement= "")

sn <- data.table::fread(file.path(outdir,
                                  "sampleNames.txt"))
files <- sn[Sample_Name == sname]

BiocInstaller::biocLite("proFIA")

library(proFIA) # DESIGNED FOR CENTROIDED DATA... MEMORY CANNOT HANDLE PROFILE DATA :-()

if(require(plasFIA)){
  path<-system.file(package="plasFIA","mzML")
  
  #Defining parameters for Orbitrap fusion.
  ppm<-2
  ppmgroup<-1
  paral<-FALSE
  fracGroup<-0.2
  k<-2
  maxo<-FALSE
  
  plasSet<-analyzeAcquisitionFIA(path,ppm=ppm,fracGroup=fracGroup,ppmgroup=ppmgroup,k=k,parallel=paral)
  
}
#set <- analyzeAcquisitionFIA(file.path(outdir, "testpeaks"), ppm = ppm, noiseEstimation = FALSE, SNT=0)
plasSet_2<-analyzeAcquisitionFIA(file.path(outdir, "testpeaks"),
                                 ppm = ppm,
                                 fracGroup = fracGroup,
                                 ppmgroup = ppmgroup,
                                 k = k,
                                 parallel = paral)

ppm = 2

files <- list.files(file.path("/Users/jwolthuis/Documents/umc/data/Data/BrSpIt/MZXML", "testpeaks"),full.names = TRUE)
f = files[1]
mzML


set <- proFIAset(file.path(outdir, "testpeaks"), ppm, parallel = FALSE, BPPARAM = NULL,
          noiseEstimation = TRUE)


# --- compare all  ---

perc_wav = 100 - (length(which(is.na(res_wav$match_mz))) / nrow(res_wav) * 100.00)
perc_mald = 100 - (length(which(is.na(res_mald$match_mz))) / nrow(res_mald) * 100.00)
perc_gauss = 100 - (length(which(is.na(res_gauss$match_mz))) / nrow(res_gauss) * 100.00)
perc_gauss_2 = 100 - (length(which(is.na(res_gauss_2$match_mz))) / nrow(res_gauss_2) * 100.00)

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
