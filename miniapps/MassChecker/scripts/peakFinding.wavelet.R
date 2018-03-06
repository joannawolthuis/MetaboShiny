peakFinding.wavelet <- function(file, scripts, outdir, thresh, resol, scanmode){
# file="./results/average_pklist/S_P53.1_pos.RData"
# file="./results/average_pklist/P53.1_pos.RData"
# scripts="./scripts"
# outdir="./results"

  
  # testing
  outdir = "~/Documents/umc/data/Data/Project 2017-026 DSM feed-10 (IMASDE-Madrid) - Saskia v Mil/RES-2017-11-15_DSM DBS Madrid/results"
  file = "~/Documents/umc/data/Data/Project 2017-026 DSM feed-10 (IMASDE-Madrid) - Saskia v Mil/RES-2017-11-15_DSM DBS Madrid/results/average_pklist/1A_pos.RData"
  thresh=2000
  resol=140000
  scanmode="positive"
  
  load(paste(outdir, "breaks.fwhm.RData", sep="/"))
  load(file)
  
  sampname = strsplit(file, "/")[[1]]
  sampname = sampname[length(sampname)]
  sampname = strsplit(sampname, "_")[[1]]
  sampname = paste(sampname[1:(length(sampname)-1)], collapse = "_")

  options(digits=16)
  int.factor=1*10^5 # Number of x used to calc area under Gaussian (is not analytic) 
  scale=2 # Initial value used to estimate scaling parameter
  width=1024
  height=768
  
  # --- wavelet magic here ---
  res <- peakDetectionCWT(sum_pos,
                          #scales = seq(0.05, 0.10, 0.01),
                          SNR.Th = 3,
                          peakThr = thresh)
  res$majorPeakInfo
  # --------------------------
  # SEE FILE 'peakfinding_mixtools.R' for some more attempts
  
  outlist.persample=cbind("samplenr"=values$nr, "mzmed.pkt"=values$mean, "fq"=values$qual, "mzmin.pkt"=values$min, "mzmax.pkt"=values$max, "height.pkt"=values$area)

  index=which(outlist.persample[,"height.pkt"]==0)
  if (length(index)>0){
    outlist.persample=outlist.persample[-index,]  
  }
  
  save(outlist.persample, file=paste(outdir, "specpks", paste(sampname, "_", scanmode, ".RData", sep=""), sep="/"))
  
  message(paste("There were", values$spikes, "spikes!"))
}

# message("Start")
# cat("==> reading arguments:\n", sep = "")
# 
# cmd_args = commandArgs(trailingOnly = TRUE)
# 
# for (arg in cmd_args) cat("  ", arg, "\n", sep="")
# 
# run(cmd_args[1], cmd_args[2], cmd_args[3], as.numeric(cmd_args[4]), as.numeric(cmd_args[5]), cmd_args[6])
# 
# message("Ready")

peakFinding.wavelet_2 <- function(file, 
                                  scripts, 
                                  outdir, 
                                  thresh, 
                                  sig_noise=6e-7){
    #file <- "/Users/jwolthuis/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/averaged/1-10D_pos.RData"
    # METADATA???
    load(file)

    fn = file.path(outdir, "peaks", basename(file))

    if(file.exists(fn)) return(NULL)
    
    sampname = gsub(basename(file), pattern = "_(pos|neg)\\.RData", replacement = "")
    # ----------------
    df <- matrix(data = averaged[[1]]@intensity,
                 ncol=1)
    rownames(df) <- averaged[[1]]@mass
    
    
    try({
      
      # this has good explanations plz read properly later
      # http://matlab.izmiran.ru/help/toolbox/wavelet/
      
      scales <- seq(1, 64, 3)
      wCoefs <- cwt(df, scales=scales, wavelet='mexh')
      localMax <- getLocalMaximumCWT(wCoefs,amp.Th = 1.5)

      # plotLocalMax(localMax,range = rng)
      # image(rng[[1]]:rng[[2]], 1:5, wCoefs[rng[[1]]:rng[[2]],1:5], 
      #       col=terrain.colors(256), 
      #       axes=TRUE, 
      #       xlab='m/z index', 
      #       ylab='CWT coefficient scale', 
      #       main='CWT coefficients')
      
      ridgeList <- getRidge(localMax,minWinSize = 5)
      
      # plotRidgeList(ridgeList, range = rng)
      
      majorPeakInfo <- identifyMajorPeaks(df, 
                                          ridgeList, 
                                          wCoefs, 
                                          SNR.Th = sig_noise, 
                                          peakScaleRange = 1,
                                          nearbyPeak = T)
      
      # seq(10, 1, by=-2)
      ## Plot the identified peaks
      # peakIndex <- majorPeakInfo$peakIndex
      # plotPeak(df, peakIndex, main=paste('Identified peaks with SNR >', 0))
      # jumps <- seq(upto,nrow(df),by = 100)
      # jumps <- seq(nrow(df), upto, by = -100)
      # 
      # peaki <- as.numeric(gsub(majorPeakInfo$peakCenterIndex, pattern = ".*_", replacement = ""))
      # 
      # for(i in 1:length(jumps)){
      #   rng = c(jumps[i], jumps[i+1])
      #   print(rng)
      #   plot(df[rng[1]:rng[2],])
      #   lines(df[rng[1]:rng[2],])
      #   # --- how many peaks are here? ---
      #   peakc = length(which(peaki %between% c(min(rng), max(rng))))
      #   print(paste("Found", peakc, "peaks in this window."))
      #   # --------------------------------
      #   Sys.sleep(2)
      #   upto = rng[2]
      # }
      
      # rng <- c(6650,
      #          6700)
      
      # plot(df[rng[[1]]:rng[[2]],])
      # lines(df[rng[[1]]:rng[[2]],])
      
      # abline(h = 2000, col="pink")
       
      # msw_peaks <- MassSpecWavelet::peakDetectionCWT(df,
      #                                                SNR.Th = 1.5E-3,
      #                                                nearbyPeak = T,
      #                                                peakScaleRange = 1,
      #                                                amp.Th = 1.5,
      #                                                peakThr = 2000,tuneIn = T)
      #peakInfo <- msw_peaks
      #majorPeakInfo = peakInfo$majorPeakInfo
      
      # peakIndex <- majorPeakInfo$allPeakIndex
      # betterPeakInfo <- MassSpecWavelet::tuneInPeakInfo(df, 
      #                                                   majorPeakInfo, 
      #                                                   peakIndex = peakIndex)
      # ----------------
      peaks = MALDIquant::createMassPeaks(mass = as.numeric(rownames(df)[majorPeakInfo$peakCenterIndex]),
                              intensity = majorPeakInfo$peakValue,
                              snr = majorPeakInfo$peakSNR,
                              metaData = list(sample = sampname))
      # ----------------
      save(x=peaks, file=fn)
    })
    # ----- return -----    
}