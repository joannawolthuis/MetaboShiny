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
  head(sum_pos)
  
  data(exampleMS)
  head(exampleMS)
  head(sum_pos)
  
  res <- peakDetectionCWT(sum_pos,
                          scales = seq(0.05, 0.10, 0.01),
                          SNR.Th = 10,
                          peakThr = thresh)
  res$majorPeakInfo
  # --------------------------
  
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
