# ============================ LIBRARY ========================

library(MALDIquant)
library(MALDIquantForeign)
library(MassSpecWavelet)
library(pbapply)
library(parallel)
library(data.table)
library(xcms)
library(gsubfn)
# library(snow)
# library(doParallel)
# library(doSNOW)
# hpc
# outdir = "/hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM/BrazilAndSpain/MZXML"
#scriptdir = "/hpc/cog_bioinf/ridder/users/jwolthuis/Pipelines/changed/Direct-Infusion-Pipeline_JOANNA"

# ======================= SET BASE PARAMETERS ======================

thresh = list(pos=2000,
              neg=2000)
snr = 3
resol = 140000
outdir = normalizePath("~/Documents/umc/data/Data/BrSpIt/MZXML/")
outdir = normalizePath("~/Documents/umc/data/Data/Project 2017-029 analysis feces (pig)_DI-MS untargeted - Myrthe Gilbert.WUR/RES_2017-12-07 piglets feces (Gilbert WUR)")

scriptdir = normalizePath("~/Google Drive/MetaboShiny/miniapps/MassChecker/scripts")
resdir = file.path(outdir, "results")
dimsThresh = 100
trim = 0.1
cores = 3
ppm = 2
pos_scans = 10:145
neg_scans = 165:300

for(f in list.files(file.path(scriptdir), full.names = T)) source(f)
for(f in list.files(file.path(scriptdir, "AddOnFunctions"), full.names = T)) source(f)

# ======================= CREATE DIRECTORIES =========================

dir.create(resdir,
           showWarnings = F)

# ========================== CLUSTER ==============================

cl = parallel::makeCluster(cores, "FORK")

# ========================== FWHM ============================



# ================= ORIG PIPELINE 'DIMS' STEP ===================

dir.create(file.path(resdir, "pklist"),
           showWarnings = F)

files = list.files(outdir,
                   pattern = "\\.raw",
                   full.names = T)

cl = parallel::makeCluster(3, "FORK")

parallel::stopCluster(cl)

pbapply::pbsapply(1:length(files), cl=cl, FUN=function(i){
  rawfile = files[[i]]
  filename = file.path(resdir,
                       "pklist",
                       paste0(gsub(basename(rawfile),
                                   pattern="\\.raw",
                                   replacement=""),
                              ".RData"))
  print(filename)
  if(file.exists(filename)) return(NULL)
  try({
    pklist <- dims_go(rawfile,
                      resdir,
                      thresh,
                      trim,
                      resol,
                      scriptdir)
    # --- return ---
    save(x=pklist, 
         file=filename)
    cat("Saved!")
  })
})

dir.create(file.path(resdir, "pklist_v2"),
           showWarnings = F)

pbapply::pbsapply(1:length(files), cl=cl, FUN=function(i){
  rawfile = files[[i]]
  filename = file.path(resdir,
                       "pklist_v2",
                       paste0(gsub(basename(rawfile),
                                   pattern="\\.raw",
                                   replacement=""),
                              ".RData"))
  print(filename)
  if(file.exists(filename)) return(NULL)
  cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/scanpolarity' -raw '$rawfile'")
  polarity <- system(cmd,intern = T)
  
  # Get time values for positive and negative scans
  posInd <- which(polarity == "positive")
  negInd <- which(polarity == "negative")
  
  # get mzvals
  specs_pos = lapply(posInd, FUN=function(scan){
    cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/printspectrum' -raw '$rawfile' -sn $scan")
    spec <- system(cmd,
                   intern = T)
    temp.list <- strsplit(spec, " ")
    #temp.list <- lapply(temp.list, FUN=function(x) if(length(x) == 2) x else NULL)
    mzvals = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[1]])))
    intensities = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[2]])))
    # --- mass spectrum obj ---
    createMassSpectrum(mzvals, 
                       intensities, 
                       metaData = list(sample=gsub(basename(rawfile), 
                                                   pattern="\\.raw", 
                                                   replacement="")))
  })
  specs_neg = lapply(negInd, FUN=function(scan){
    cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/printspectrum' -raw '$rawfile' -sn $scan")
    spec <- system(cmd,intern = T)
    temp.list <- strsplit(spec, " ")
    #temp.list <- lapply(temp.list, FUN=function(x) if(length(x) == 2) x else NULL)
    mzvals = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[1]])))
    intensities = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[2]])))
    # --- mass spectrum obj ---
    createMassSpectrum(mzvals, 
                       intensities, 
                       metaData = list(sample=gsub(basename(rawfile), 
                                                   pattern="\\.raw", 
                                                   replacement="")))
    
  })
  align_pos <- alignSpectra(specs_pos, tolerance = 2E-6)
  align_neg <- alignSpectra(specs_neg, tolerance = 2E-6)
  
  pos <- averageMassSpectra(align_pos,
                            method = "sum")
  neg <- averageMassSpectra(align_neg,
                            method = "sum")
  pklist = list(pos = pos, 
                neg = neg)
  # --- return ---
  save(x=pklist, 
       file=filename)
  cat("Saved!")
})

# ======================= LOAD SAMPLE INFO ===========================

sn <- data.table::fread(file.path(outdir,
                                  "sampleNames.txt"))

sn$batch <- as.numeric(as.factor(gsub(sn$File_Name,
                                      pattern = "_\\d\\d\\d$",
                                      replacement="")))

all_fns <- paste0(file.path(resdir, "pklist"),
                  sn$File_Name,
                  ".RData")

# ===================== LOAD IN SAMPLES ======================

repl.pattern = split(sn,f = sn$Sample_Name)

dir.create(file.path(resdir, "spectra"),
           showWarnings = F)

spectra = pbapply::pblapply(repl.pattern, cl=cl, FUN=function(repl.set){
  repl.spec <- lapply(1:nrow(repl.set), FUN=function(i){
    f = file.path(resdir, "pklist",
                  paste0(repl.set[i, "File_Name"], ".RData"))
    load(f)
    # ----- create mass spectra -----
    pos <- pklist$pos
    neg <- pklist$neg
    list(pos = pos,
         neg = neg)
  })
  repl.spec
})

spec = lapply(spectra, function(x) x$pos)

fn_pos = file.path(resdir, "spectra", "spectra_positive.RData")
save(x=spec, file=fn_pos)

spec = lapply(spectra, function(x) x$neg)

fn_neg = file.path(resdir, "spectra", "spectra_negative.RData")
save(x=spec, file=fn_neg)

spec <- NULL

# ====================== SMOOTH =====================

dir.create(file.path(resdir, "smoothed"),
           showWarnings = F)

files = list.files(file.path(resdir, "pklist"),
                   full.names = T)


pbsapply(files, cl=cl, FUN=function(f){
  load(f)
  for(mode in c("pos","neg")){
    spec = pklist[[mode]]
    smoothed = smoothIntensity(spec, method="SavitzkyGolay", halfWindowSize=10)
    # --- save --- 
    fn = file.path(resdir, "smoothed", paste0(gsub(basename(f), 
                                                   pattern="\\.RData", 
                                                   replacement=""),"_", mode, ".RData"))
    save(x=smoothed, file=fn)
  }
})

# pbsapply(modes, cl=cl, FUN=function(mode){
#   load(file.path(resdir, "spectra", "spectra_", mode,".RData"))
#   smoothed = smoothIntensity(spec, method="SavitzkyGolay", halfWindowSize=10)
#   # --- save --- 
#   fn = file.path(resdir, "smoothed", paste0("smooth_", mode, ".RData"))
#   save(x=smoothed, file=fn)
# })

# =================== BASELINE ===================

dir.create(file.path(resdir, "baseline"),
           showWarnings = F)

# pbsapply(modes, cl=cl, FUN=function(mode){
#   load(file.path(resdir, "smoothed", paste0("smooth_", mode, ".RData")))
#   nobase = removeBaseline(smoothed, method="SNIP", iterations=100)
#   # --- save --- 
#   fn = file.path(resdir, "baseline", paste0("nobase_", mode, ".RData"))
#   save(x=nobase, file=fn)
# })

pbsapply(modes, cl=cl, FUN=function(mode){
  load(file.path(resdir, "smoothed", paste0("smooth_", mode, ".RData")))
  nobase = removeBaseline(smoothed, method="SNIP", iterations=100)
  # --- save --- 
  fn = file.path(resdir, "baseline", paste0("nobase_", mode, ".RData"))
  save(x=nobase, file=fn)
})



# ================= CALIBRATE =================

dir.create(file.path(resdir, "calibrated"),
           showWarnings = F)

pbsapply(modes, cl=cl, FUN=function(mode){
  load(file.path(resdir, "baseline", paste0("nobase_", mode, ".RData")))
  calib = calibrateIntensity(nobase, method="TIC")
  # --- save --- 
  fn = file.path(resdir, "calibrated", paste0("calib_", mode, ".RData"))
  save(x=calib, file=fn)
})

# ================== ALIGN =====================

dir.create(file.path(resdir, "aligned"),
           showWarnings = F)

pbsapply(modes, cl=cl, FUN=function(mode){
  load(file.path(resdir, "calibrated", paste0("calib", mode, ".RData")))
  aligned = alignSpectra(calib, halfWindowSize=20, SNR=3, tolerance=2E-6, warpingMethod="lowess")
  # --- save --- 
  fn = file.path(resdir, "aligned", paste0("align_", mode, ".RData"))
  save(x=aligned, file=fn)
})

# ================= AVERAGE ===================

dir.create(file.path(resdir, "averaged"),
           showWarnings = F)


cl = parallel::makeCluster(3, "FORK")

parallel::stopCluster(cl)

repl.pattern = split(sn,f = sn$Sample_Name)

averaged <- pbapply::pblapply(repl.pattern, cl=cl, FUN=function(repl.set){
  repl.spec <- lapply(1:nrow(repl.set), FUN=function(i){
    f = file.path(resdir, "pklist", paste0(repl.set[i, "File_Name"], ".RData"))
    load(f)
    # ----- create mass spectra -----
    pos <- pklist$pos
    neg <- pklist$neg
    list(pos = pos, 
         neg = neg)
  })
  print(repl.spec)
  # ---------------
  fn_pos = file.path(resdir, "averaged",
                     paste0(unique(repl.set$Sample_Name), "_pos.RData"))
  fn_neg = file.path(resdir, "averaged",
                     paste0(unique(repl.set$Sample_Name), "_neg.RData"))
  # average samples
  averaged <- averageMassSpectra(lapply(repl.spec, function(x) x$pos),
                                 labels=repl.set$Sample_Name,
                                 method="mean")
  save(x=averaged, file=fn_pos)
  averaged <- NULL # RESET
  averaged <- averageMassSpectra(lapply(repl.spec, function(x) x$neg),
                                 labels=repl.set$Sample_Name,
                                 method="mean")
  
  save(x=averaged, file=fn_neg)
})

# pbsapply(modes, cl=cl, FUN=function(mode){
#   load(file.path(resdir, "aligned", paste0("align_", mode, ".RData")))
#   averaged = averageMassSpectra(aligned, labels = sn$Sample_Name, method="mean")
#   # --- save --- 
#   fn = file.path(resdir, "averaged", paste0("average_", mode, ".RData"))
#   save(x=averaged, file=fn)
# })

# ======================== PEAK CALLING ===========================

dir.create(file.path(resdir, "peaks"),
           showWarnings = F)

files = list.files(file.path(resdir, "averaged"),
                   full.names = T)

cl = parallel::makeCluster(3, "FORK")

parallel::stopCluster(cl)


peaks <- pblapply(files, cl=cl, FUN=function(f){
  # METADATA???
  load(f)
  print(f)
  fn = file.path(resdir, "peaks", basename(f))
  if(file.exists(fn)) return(NULL)
  sampname = gsub(basename(f), pattern = "_(pos|neg)\\.RData", replacement = "")
  print(fn)
  # --- smoothed ---
  smoothed <- MALDIquant::smoothIntensity(averaged[[1]])
  # ----------------
  # df <- matrix(data = averaged[[1]]@intensity,
  #              ncol=1)
  # rownames(df) <- averaged[[1]]@mass
  df <- matrix(data = smoothed@intensity,
               ncol=1)
  rownames(df) <- smoothed@mass
  try({
    msw_peaks <- MassSpecWavelet::peakDetectionCWT(df,
                                                   SNR.Th = snr,
                                                   nearbyPeak = T,
                                                   peakThr = 2000)
    
    peakInfo <- msw_peaks
    majorPeakInfo = peakInfo$majorPeakInfo
    peakIndex <- majorPeakInfo$potentialPeakIndex
    betterPeakInfo <- tuneInPeakInfo(df, majorPeakInfo)
    # ----------------
    peaks = createMassPeaks(mass = as.numeric(rownames(df)[betterPeakInfo$peakCenterIndex]),
                            intensity = betterPeakInfo$peakValue,
                            metaData = list(sample = sampname))
    print(paste("(Smooth) Found", length(betterPeakInfo$peakCenterIndex), "peaks!"))
    
    # ----------------
    save(x=peaks, file=fn)
  })
  # ----- return -----    
})

dir.create(file.path(resdir, "peaks_nosmooth"),
           showWarnings = F)

peaks <- pblapply(files, cl=cl, FUN=function(f){
  # METADATA???
  load(f)
  print(f)
  fn = file.path(resdir, "peaks_nosmooth", basename(f))
  if(file.exists(fn)) return(NULL)
  sampname = gsub(basename(f), pattern = "_(pos|neg)\\.RData", replacement = "")
  print(fn)
  # ----------------
  df <- matrix(data = averaged[[1]]@intensity,
               ncol=1)
  rownames(df) <- averaged[[1]]@mass

  try({
    msw_peaks <- MassSpecWavelet::peakDetectionCWT(df,
                                                   SNR.Th = snr,
                                                   nearbyPeak = T)
    
    peakInfo <- msw_peaks
    majorPeakInfo = peakInfo$majorPeakInfo
    peakIndex <- majorPeakInfo$potentialPeakIndex
    betterPeakInfo <- tuneInPeakInfo(df, majorPeakInfo)
    # ----------------
    peaks = createMassPeaks(mass = as.numeric(rownames(df)[betterPeakInfo$peakCenterIndex]),
                            intensity = betterPeakInfo$peakValue,
                            metaData = list(sample = sampname))
    # ----------------
    save(x=peaks, file=fn)
  })
  # ----- return -----    
})

# ------ my own peak finder !!! ------

f <- files[[1]]
load(f)
vec <- averaged[[1]]@intensity
names(vec) <- averaged[[1]]@mass
rl <- rle(vec == 0)
i1 <- rl$lengths>1 & rl$values
lst <- split(vec, rep(cumsum(c(TRUE, i1[-length(i1)])), rl$lengths)) 
split_peaks = lapply(lst, function(x) x[which(x!=0)])

library(mclust)

find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}
# --- nvmfor now ---

library(mclust)

cl = makeCluster(3, "FORK")

res <- pbapply::pblapply(1:length(split_peaks),cl=0, FUN=function(i){
  #print(i/length(split_peaks) * 100.0)
  #print(i)
  peak <- split_peaks[[i]]
  
  if(length(peak)>1){
    #plot(names(peak),peak)
    
    peak_hist = unlist(sapply(1:length(peak), FUN = function(j){
      counts = peak[j]
      mz = as.numeric(names(peak[j]))
      # --------------------
      rep(x=mz, times=counts)
    }))
    
    #hist(peak_hist, prob=T)

    tp<-pastecs::turnpoints(peak)
    ncomp = length(which(tp$peaks))
    
    #points(names(peak)[tp$peaks],peak[tp$peaks],col="red")
    
    if(ncomp > 1){
      try({
        mix <- Mclust(peak_hist, 
                      G = ncomp,verbose = FALSE)   
        
        fit = list(lambda = mix$parameters$pro,
                   mu =  mix$parameters$mean,
                   sigma = if(length(mix$parameters$variance$sigmasq) == 1){
                     rep(sqrt(mix$parameters$variance$sigmasq), length(mix$parameters$mean))
                   }else{
                     sqrt(mix$parameters$variance$sigmasq)
                   })
      })
    }else{
      fit = list(lambda = 1, mu = mean(peak_hist), sigma = sd(peak_hist))
    }

    # for(i in 1:length(fit$mu)){
    #    curve(dnorm(x, mean = fit$mu[[i]],sd = fit$sigma[[i]]), col = i, add = TRUE)
    # }
    
    peakmz = unique(fit$mu)
    
    N <- length(peak_hist)
    draw <- sample(length(fit$lambda), N, replace=TRUE, prob=fit$lambda)

    peakinfo <- lapply(1:length(peakmz), FUN=function(i){
      mz = peakmz[i]
      int = length(which(draw == i))
      # ---------
      data.table(mzmed = mz, 
                 int = int)
    })
    if(length(peakinfo) > 1){
      peakinfo <- rbindlist(peakinfo)
    }
    #print(peakinfo)
    peakinfo
  }
  #Sys.sleep(2)
})


BIC <- mclustBIC(peak)
Mclust(peak, x = BIC)


hist(peak_hist, probability = T)
for(i in 1:length(mix$lambda)){
  curve(dnorm(x, mean = mix$mu[[i]],sd = mix$sigma[[i]]), col = i, add = TRUE)
}

npoints = length(peak_hist)

probs <- mix$lambda
m <- mix$mu
s <- mix$sigma

N <- length(peak_hist)
grp <- sample(length(probs), N, replace=TRUE, prob=probs)
x <- rnorm(N, m[grp], s[grp])
hist(x)
plot(density(x))
hist(peak_hist)

sampleMixture=function(prop,means,dev){
  # Generate a uniformly distributed random number between 0 and 1
  # in order to choose between the two component distributions
  distTest=runif(1)
  if(distTest<prop[1]){
    # Then sample from the first component of the mixture
    sample=rnorm(1,mean=means[1],sd=dev[1])
    names(sample) = "1"
  }else{
    # Sample from the second component of the mixture
    sample=rnorm(1,mean=means[2],sd=dev[2])
    names(sample) = "2"
  }
  return(sample)
}

# Generate a single sample
sampleMixture(probs, m, s)

Mclust(peak_hist)
# Generate 100 samples and plot resulting distribution
samples=replicate(length(peak_hist),sampleMixture(probs, m, s))
plot(density(samples))
names(samples)
length(which(names(samples) == "1"))

# library(matconv)
# hMaps <- makeFuncMaps(
#   pathDict = system.file("extdata", "HiebelerDict.txt", package = "matconv"))
# source(system.file("extdata", "defDataConv.R", package = "matconv"))
# 
# matfiles <- list.files("/Users/jwolthuis/Downloads/S2_File/m_files/",pattern = "\\.m$", full.names = T)
# for(f in matfiles){
#   code <- unlist(readLines(f))
#   out <- mat2r(code, funcConverters = hMaps, dataConverters = dataConvs)
#   newf <- gsub(pattern = "\\.m$", replacement = ".R", x = f)
#   write(x =  out$rCode, file = newf)
# }


pbsapply(split_peaks, function(peak){
  if(length(peak) > 0){
    peak_hist = unlist(sapply(1:length(peak), FUN = function(i){
      counts = peak[i]
      mz = as.numeric(names(peak[i]))
      times = counts/binsize
      #print(mz)
      multiplied = rep(x=mz, times=counts)
      multiplied
    }))
    peakcount = length(find_peaks(peak, m=2))
    print(paste("Found...", peakcount, "peak(s)!"))
    hist(peak_hist, probability = T)
    try({
      mix <- mixtools::normalmixEM(peak_hist)
      print(mix)
      #plot(names(peak), peak)
      #mixtools::plot.mixEM(mix,whichplots = 2)
      for(i in 1:length(mix$lambda)){
        curve(dnorm(x, mean = mix$mu[[i]],sd = mix$sigma[[i]]), col = i, add = TRUE)
      }
    })
    # }else{
    #   fit <- fitdistr(peak_hist, "normal")
    #   para <- fit$estimate
    #   hist(peak_hist, probability = T)
    #  
    # }
    Sys.sleep(1)
  }
})



# ------------------------------------

dir.create(file.path(resdir, "peaks_maldiquant"),
           showWarnings = F)

peaks <- pblapply(files, cl=cl, FUN=function(f){
  # METADATA???
  load(f)
  print(f)
  fn = file.path(resdir, "peaks_maldiquant", basename(f))
  if(file.exists(fn)) return(NULL)
  sampname = gsub(basename(f), pattern = "_(pos|neg)\\.RData", replacement = "")
  print(fn)
  # ----------------
  df <- matrix(data = averaged[[1]]@intensity,
               ncol=1)
  rownames(df) <- averaged[[1]]@mass
  
  peaks <- MALDIquant::detectPeaks(averaged,snr=10)[[1]]

  # ----------------
  
  save(x=peaks, file=fn)
  
  # ----- return -----    
})


# ============================= CREATE PEAK TABLES ==================================

dir.create(file.path(resdir, "specpks_all"),
           showWarnings = F)


for(mode in c("pos", "neg")){
  pk_files <- list.files(file.path(resdir, "peaks_nosmooth"),
                         pattern= paste0("_",mode,"\\.RData"),
                         full.names = T)
  # Gather sub peaktables
  peaktables <- pblapply(pk_files, cl=cl, FUN=function(f){
    load(f)
    tbl <- data.table(mz = peaks@mass,
                      intensity = peaks@intensity)
    names(tbl) = c("mz", gsub(basename(f), 
                              pattern="_(pos|neg)\\.RData", 
                              replacement=""))
    # --- return ---
    tbl
  })
  
  outpgrlist = rbindlist(peaktables, use.names = T,fill = T)

  # set first table to merge w/ later
  #outpgrlist = peaktables[[1]]
  
  # pbsapply(peaktables, FUN=function(pktable){
  #   outpgrlist <<- merge(outpgrlist, 
  #                        pktable, 
  #                        by="mz",
  #                        all=T)
  # })
  
  names(outpgrlist) = c("mzmed", gsub(basename(names(outpgrlist)[-1]), 
                                      pattern="\\.RData", 
                                      replacement=""))
  print(paste0("--- ", mode, " ---"))
  print(dim(outpgrlist))
  
  save(x=outpgrlist, 
       file=file.path(resdir, "specpks_all", paste0("peaks_",
                                                    mode,
                                                    ".RData")))
}

blocky <- function(table, dim=10) table[1:dim, 1:dim]

# ================================ GROUP PEAKS =====================================

modes = c("pos",
          "neg")

sn <- fread(file.path(resdir, "sampleNames.txt"))
sn$batch <- as.numeric(as.factor(gsub(sn$File_Name,
                                      pattern = "_\\d\\d\\d$",
                                      replacement="")))
sn_no_qc <- sn[!(Sample_Name %like% "QC")]
sn_qc <- sn[(Sample_Name %like% "QC")]
sn_adj <- as.data.table(unique(rbind(sn_no_qc, sn_qc)))

# maldiquant method

dir.create(file.path(resdir, "specpks_grouped_mdq"),
           showWarnings = F)

library(pbapply)
library(data.table)
library(MALDIquant)
library(gsubfn)

for(mode in modes){

  # --- load file ---
  fn =  file.path(resdir, "specpks_grouped_mdq", paste0("grouped_",mode,".RData"))
  
  pk_files <- list.files(file.path(resdir, "peaks"),
                         pattern= paste0("_",mode,"\\.RData"),
                         full.names = T)
  
  # Gather sub peaktables
  peaktables <- pblapply(pk_files, cl=0, FUN=function(f){
    load(f)
    peaks@metaData$sampname = gsub(basename(f), 
                                   pattern="_(pos|neg)\\.RData", 
                                   replacement="")
    peaks
  })
  
  binned <- MALDIquant::binPeaks(l = peaktables, tolerance = 3E-6, method = "relaxed") 
  names(binned) <- sapply(peaktables, function(x) x@metaData$sampname)
  
  peaktables <- NULL; gc()
  
  # Gather sub peaktables
  peaktables_binned <- pblapply(1:length(binned), cl=0, FUN=function(i){
    peaks = binned[[i]]
    tbl <- data.table(mz = peaks@mass,
                      intensity = peaks@intensity)
    colnames(tbl)
    sampname = names(binned)[i]
    colnames(tbl) = c("mzmed", sampname)
    # --- return ---
    tbl
  })
  
  binned <- NULL; gc()
  
  outlist = data.table::rbindlist(peaktables_binned, use.names = T, fill = TRUE)

  outpgrlist <- rowsum(outlist, group = outlist$mzmed, na.rm = T, reorder = T)
  outpgrlist$mzmed <- rownames(outpgrlist)
  
  new_colnames = pbsapply(colnames(outpgrlist[,-1]), FUN=function(samp){
    rowsi = sn_adj[Sample_Name == samp]
    batch = unique(rowsi$batch)
    newName = gsubfn::fn$paste("*$batch*$samp")
    newName
  })
  
  colnames(outpgrlist)[2:ncol(outpgrlist)] <- new_colnames
  
  outlist <- NULL; gc()

  outpgrlist[outpgrlist == 0] <- NA
  
  outpgrlist_filt <- outpgrlist[-which(rowMeans(is.na(outpgrlist[,-1])) > 0.6), ]
  
  print(mode)
  print(dim(outpgrlist))
  
  # remove all that only has
  #f =  file.path(resdir, "specpks_grouped_mdq", paste0("grouped_",mode,".RData"))
  #save(outpgrlist, file=f)
  data.table::fwrite(x = outpgrlist_filt,
                     file = file.path(resdir, "specpks_grouped_mdq", paste0("grouped_",
                                                               mode,
                                                               ".csv")))
  
  outpgrlist <- NULL; gc()
}

# -----------------

dir.create(file.path(resdir, "specpks_grouped"),
           showWarnings = F)

for(mode in modes){
  f = file.path(resdir, "specpks_all", paste0("peaks_",
                                              mode,
                                              ".RData"))
  load(f)
  # --- load file ---
  fn =  file.path(resdir, "specpks_grouped", paste0("grouped_",mode,".RData"))
  peakmat <- data.table(mzmed = outpgrlist$mzmed)
  peakmat <- apply(peakmat, MARGIN = 2, as.numeric)
  cl = parallel::makeCluster(3, "FORK")
  groups <- mzClustGeneric(as.matrix(peakmat), 
                           mzppm=ppm,
                           shinyprog=FALSE)
  outlist <- outpgrlist
  groupNames <- colnames(outlist[,-1])
  #parallel::stopCluster(cl)
  # --- cl prep ---
  outrows <- pblapply(1:length(groups$idx), cl=cl, FUN=function(i){
    mzidx <- groups$idx[[i]]
    # FIND PEAKS IN THIS GROUP
    members <- outlist[mzidx, ]
    # --- get intensities ---
    mzmed = as.numeric(mean(members$mzmed))
    mzmin = as.numeric(min(members$mzmed))
    mzmax = as.numeric(max(members$mzmed))
    #fq.worst.pgrp = as.numeric(max(members$fq))
    #fq.best.pgrp = as.numeric(min(members$fq))
    ints.allsamps = rep(0, length(groupNames))
    names(ints.allsamps) = groupNames # same order as sample list!!!
    # # Check for each sample if multiple peaks exists, if so take the sum!
    nrsamples <- nrow(members)
    ints.allsamps <- colSums(outlist[mzidx, -1],
                             na.rm = T)
    #print(ints.allsamps)
    # --- make dt ---
    outpgrlist.a = cbind('mzmed' = mzmed,
                         #"fq.best"=fq.best.pgrp, 
                         #"fq.worst"=fq.worst.pgrp, 
                         nrsamples, 
                         'mzmin' = mzmin, 
                         'mzmax' = mzmax
    )
    outpgrlist.b <- data.frame(as.list(ints.allsamps))
    outpgrlist <- cbind(outpgrlist.a, outpgrlist.b)
    #print(dim(outpgrlist))
    # -----------------------
    as.data.table(outpgrlist)
  })
  outpgrlist <- rbindlist(outrows)
  f =  file.path(resdir, "specpks_grouped", paste0("grouped_",mode,".RData"))
  save(outpgrlist, file=f)
}

# ============================ CORRECT BATCH EFFECT ================================



# ================================ FILL MISSING ====================================


for(f in list.files(file.path(scriptdir, "AddOnFunctions"), full.names = T)) source(f)

dir.create(file.path(resdir, "specpks_filled"),
           showWarnings = F)

cl <<- makeSOCKcluster(3,
                       outfile=file.path(resdir, 
                                         "error_log.txt"))
registerDoSNOW(cl)
stopCluster(cl)
mode="pos"
thresh=2000
snr = 3
resol = 140000

for(mode in c("pos", "neg")){
  f =  file.path(resdir, "grouped", paste0("grouped_",mode,".RData"))
  load(f)
  colnames(outpgrlist) <- gsub(colnames(outpgrlist), 
                               pattern="\\.", 
                               replacement="-")
  colnames(outpgrlist) <- gsub(colnames(outpgrlist), 
                               pattern="^X", 
                               replacement="")
  

  keepcols = pblapply(1:(nrow(outpgrlist)), FUN=function(i){
    vals = as.numeric(outpgrlist[i,-c(1:4)])
    vals[vals < thresh] <- NA 
    missing = which(is.na(vals))
    percmiss = length(missing) / length(vals) * 100
    if(percmiss < 100) outpgrlist[i,]
    #hist(vals, breaks = 100)
    #print(median(vals,na.rm = T))
    #Sys.sleep(3)
  })
  
  outpgrlist = rbindlist(keepcols)
  
  scanmode = mode
  # --- get all unique mz vals ---
  mzvals = unique(outpgrlist$mzmed)
  print(paste("This experiment found", length(mzvals), "peaks in", paste0(mode, " mode.")))
  
  # ---------
  cl <- parallel::makeCluster(cores, "FORK")
  #doSNOW::registerDoSNOW(cl)
  # ---------------
  gaussians <- pbsapply(1:length(mzvals), cl=cl, FUN=function(i){
                         area = generateGaussian(mzvals[i],
                                                 thresh,
                                                 resol,
                                                 FALSE,
                                                 scanmode,
                                                 int.factor = 1*10^5, 1, 1)$area
                         area
                       })
  # --- prep ---
  ref_table <<- data.table(mz = mzvals,
                          int = gaussians)
  print(head(ref_table))
  setkey(ref_table, mz)
  # --- loop 2 ---
  outrows <- pblapply(1:nrow(outpgrlist), cl=cl, FUN=function(i, ref_table, outpgrlist){
                       row = outpgrlist[i,]
                       rowstart = row[,c(1:4)]
                       gauss <- ref_table[row[,1],]
                       area = unlist(gauss$int)
                       # loop and fill
                       filled_row <- sapply(row[,7:length(row)], FUN=function(samp){
                         if(is.na(samp)){
                           rand_gauss <- rnorm(n=1, 
                                               mean=area, 
                                               sd=0.25*area)
                           rand_gauss
                         } else{
                           samp}
                       })
                       res <- cbind(rowstart, 
                                    t(filled_row), 
                                    avg.int=mean(filled_row))
                       # # --- return ---
                       as.data.table(res)
                     }, ref_table = ref_table, outpgrlist = outpgrlist)
  # stop cluster
  stopCluster(cl)
  # combine
  outpgrlist <- rbindlist(outrows)
  # write
  data.table::fwrite(x = outpgrlist,
                     file = file.path(resdir, "filled", paste0("filled_gauss_",
                                                                       mode,
                                                                       ".csv")))
}


# --- fill alternative ---

outdir = normalizePath("~/Documents/umc/data/Data/BrazilAndSpain/MZXML")
scriptdir = normalizePath("~/Google Drive/MetaboShiny/miniapps/MassChecker/scripts")
resdir = file.path(outdir, "results")

sn <- fread(file.path(resdir, "sampleNames.txt"))
sn$batch <- as.numeric(as.factor(gsub(sn$File_Name,
                                      pattern = "_\\d\\d\\d$",
                                      replacement="")))
sn_no_qc <- sn[!(Sample_Name %like% "QC")]
sn_qc <- sn[(Sample_Name %like% "QC")]
sn_adj <- as.data.table(unique(rbind(sn_no_qc, sn_qc)))
perc_keep = 99
fill_mode = "rf"

for(mode in modes){
  f =  file.path(resdir, "grouped", paste0("grouped_",mode,".RData"))
  load(f)
  
  outpgrlist[outpgrlist == 0] <- NA
  
  if(perc_keep < 100){
    keepcols = pblapply(1:(nrow(outpgrlist)), FUN=function(i){
      vals = as.numeric(outpgrlist[i,-c(1:4)])
      missing = which(is.na(vals))
      percmiss = length(missing) / length(vals) * 100
      if(percmiss < perc_keep) outpgrlist[i,]
    })
    filt_outpgrlist = rbindlist(keepcols)
  }else{
    filt_outpgrlist = outpgrlist
  }
  
  print("Left after filtering:")
  print(nrow(filt_outpgrlist))
  
  outpgrlist.1 <- filt_outpgrlist[,c(1:4)]
  
  outpgrlist.2 <- filt_outpgrlist[,-c(1:4)]

  # -------------------------------
  
  # filled = pblapply(5:(ncol(filt_outpgrlist)), FUN=function(i){
  #   sampCol = filt_outpgrlist[,..i][[1]]
  #   vals = as.numeric(sampCol)
  #   missing = which(is.na(vals))
  #   if(length(missing) > 0){
  #     fillval = .5 * min(vals,na.rm = T)
  #     sampCol[missing] <- c(fillval)
  #   }
  #   sampCol
  # })
  # 
  # outpgrlist.2 = as.data.frame(do.call(cbind, filled))

  colnames(outpgrlist.2) = colnames(filt_outpgrlist[,-c(1:4)])
   
  colnames(outpgrlist.2) <- gsub(colnames(outpgrlist.2), 
                                 pattern="\\.", 
                                 replacement="-")
  colnames(outpgrlist.2) <- gsub(colnames(outpgrlist.2), 
                                 pattern="^X", 
                                 replacement="")


  outpgrlist.fill <- switch(fill_mode,
                            gauss = {
                              mzvals = unique(outpgrlist.1$mzmed)
                              print(paste("This experiment found", length(mzvals), "peaks in", paste0(mode, " mode.")))
                              
                              # ---------
                              cl <- parallel::makeCluster(cores, "FORK")
                              #doSNOW::registerDoSNOW(cl)
                              # ---------------
                              gaussians <- pbsapply(1:length(mzvals), cl=cl, FUN=function(i){
                                area = generateGaussian(mzvals[i],
                                                        thresh,
                                                        resol,
                                                        FALSE,
                                                        scanmode,
                                                        int.factor = 1*10^5, 1, 1)$area
                                area
                              })
                              # --- prep ---
                              ref_table <<- data.table(mz = mzvals,
                                                       int = gaussians)
                              print(head(ref_table))
                              setkey(ref_table, mz)
                              # --- loop 2 ---
                              outrows <- pblapply(1:nrow(outpgrlist.2), cl=cl, FUN=function(i, ref_table, outpgrlist){
                                row = outpgrlist.2[i,]
                                mz = mzvals[i]
                                gauss <- ref_table[mz,]
                                area = unlist(gauss$int)
                                # loop and fill
                                filled_row <- sapply(row[,7:length(row)], FUN=function(samp){
                                  if(is.na(samp)){
                                    rand_gauss <- rnorm(n=1, 
                                                        mean=area, 
                                                        sd=0.25*area)
                                    rand_gauss
                                  } else{
                                    samp}
                                })
                                filled_row
                              }, ref_table = ref_table, outpgrlist = outpgrlist)
                              # stop cluster
                              stopCluster(cl)
                              # combine
                              outpgrlist <- rbindlist(outrows)
                            },
                            rf = {
                              library(missForest)
                              missForest(outpgrlist.2,
                                         verbose = TRUE,
                                         maxit=30)$ximp
                            },
                            cart = {
                              library(mice)
                              predictorMatrix <- matrix(0, nrow = ncol(outpgrlist.2), ncol = ncol(outpgrlist.2), dimnames = list(names(outpgrlist.2),
                                                                                                                                 names(outpgrlist.2)))
                              predRows = pblapply(1:nrow(predictorMatrix), FUN=function(i){
                                row = predictorMatrix[i,]
                                samp = rownames(predictorMatrix)[i]
                                sampBatch = unique(sn_adj[Sample_Name == samp, "batch"])
                                batchBuddies = sn_adj[batch == sampBatch, "Sample_Name"]
                                chRows = which(colnames(predictorMatrix) %in% batchBuddies$Sample_Name)
                                row[chRows] <- c(1)
                                # -----------------
                                res = as.data.frame(t(row), row.names = samp)
                                res
                              })

                              predictorMatrix_adj = as.matrix(rbindlist(predRows), ncol(predictorMatrix),
                                                              ncol = ncol(predictorMatrix),
                                                              dimnames = list(names(predictorMatrix),
                                                                              names(predictorMatrix)))
                              diag(predictorMatrix_adj) <- 0
                              # -------------------------
                              myImp = mice(outpgrlist.2,
                                            pred=predictorMatrix_adj,
                                            method="pmm")
                              # -------------------------
                              complete(myImp)
                            }, mi={
                              mdf <- missing_data.frame(outpgrlist.2) # warnings about missingness patterns
                              mdf <- change(mdf, y=colnames(mdf), what = "type", to = "positive-continuous")
                              res <- show(mdf)
                              imputations <- mi(mdf)
                              outpgrlist.fill <- complete(imputations)
                              outpgrlist.fill <- outpgrlist.fill[,-grep(pattern = "missing", x = colnames(outpgrlist.fill))]
                              # ------------------------
                              outpgrlist.fill
                            })

  colnames(outpgrlist.fill) <- colnames(outpgrlist.2)
  print(complete.cases(outpgrlist.fill))
  
  new_colnames = pbsapply(colnames(outpgrlist.fill), FUN=function(samp){
    rowsi = sn_adj[Sample_Name == samp]
    batch = unique(rowsi$batch)
    newName = gsubfn::fn$paste("*$batch*$samp")
    newName
  })
  
  colnames(outpgrlist.fill) = new_colnames
  
  # --------------------------------------------------
  
  outpgrlist = cbind(outpgrlist.1, outpgrlist.fill)
  
  # --------------------------------------------------------

  data.table::fwrite(x = outpgrlist,
                     file = file.path(resdir, "filled", paste0("filled_",
                                                               fill_mode,
                                                               "_",     
                                                               mode,
                                                               ".csv")))
}

# ================ MIXNORM BATCH CORRECTION ===========

library(metabomxtr)

?mixnorm
# ================ COMBAT BATCH CORRECTION =============

outdir = normalizePath("~/Documents/umc/data/Data/BrazilAndSpain/MZXML")
scriptdir = normalizePath("~/Google Drive/MetaboShiny/miniapps/MassChecker/scripts")
resdir = file.path(outdir, "results")
resdir= outdir

for(mode in c("pos", "neg")){
  sn <- fread(file.path(resdir, "sampleNames.txt"))
  f = file.path(resdir, "filled", paste0("filled_", mode, ".RData"))
  load(f)
  outpgrlist = as.data.table(outpgrlist)
  colnames(outpgrlist) <- gsub(colnames(outpgrlist), 
                               pattern="\\.", 
                               replacement="-")
  colnames(outpgrlist) <- gsub(colnames(outpgrlist), 
                               pattern="^X", 
                               replacement="")
  # == zero ones under threshhold ===
  outpgrlist.prefix <- outpgrlist[,c("mzmed", "nrsamples", "mzmin", "mzmax")]
  outpgrlist[,c("mzmed", "nrsamples", "mzmin", "mzmax") := NULL]
  
  colnames(outpgrlist)
  sn_no_qc <- sn[!(Sample_Name %like% "QC")]
  sn_qc <- sn[(Sample_Name %like% "QC")]
  sn_no_qc$outcome <- as.factor(gsub(sn_no_qc$Sample_Name, 
                                                pattern="-.*$", 
                                                replacement=""))
  sn_qc$outcome <- c("QC")
  sn_adj <- as.data.table(unique(rbind(sn_no_qc, sn_qc)))
  sn_adj <- sn_adj[Sample_Name %in% colnames(outpgrlist)]
  csv_pheno <-unique(data.frame(sample = sn_adj$Sample_Name,
                     outcome = as.numeric(as.factor(sn_adj$outcome)),
                     batch = sn_adj$batch))
  nrow(csv_pheno)
  ncol(outpgrlist)
  colnames(outpgrlist)
  mod = model.matrix(~as.factor(outcome), data=csv_pheno)
  mod0 = model.matrix(~1,data=csv_pheno)
  # parametric adjustment
  
  combat_edata = ComBat(dat=outpgrlist, batch=csv_pheno$batch)#, mod=mod, par.prior=TRUE, prior.plots=FALSE)
  
  outpgrlist.norm <- cbind(outpgrlist.prefix, combat_edata)
  
  data.table::fwrite(x = outpgrlist.norm,
                     file = file.path(resdir, "specpks_filled", paste0("batchcorr_",
                                                                       mode,
                                                                       ".csv")))
  }
  
  
