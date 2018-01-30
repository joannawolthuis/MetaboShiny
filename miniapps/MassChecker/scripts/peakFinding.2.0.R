peakFinding.2.0 <- function(file, scripts, outdir, thresh, resol, scanmode){
# file="./results/average_pklist/S_P53.1_pos.RData"
# file="./results/average_pklist/P53.1_pos.RData"
# scripts="./scripts"
# outdir="./results"
# thresh=2000
# resol=140000
# scanmode="positive"
  
  dir.create(paste(outdir, "Gaussian_fit", sep="/"),showWarnings = F)
  
  load(paste(outdir, "breaks.fwhm.RData", sep="/"))
  load(file)
  
  sampname = strsplit(file, "/")[[1]]
  sampname = sampname[length(sampname)]
  sampname = strsplit(sampname, "_")[[1]]
  sampname = paste(sampname[1:(length(sampname)-1)], collapse = "_")

#  install.packages("Cairo", lib ="/hpc/local/CentOS6/dbg_mz/R-3.1.2/library")
#  library( "yourLibrary", lib.loc = "/hpc/local/osversion/group/path" )
#  library("Cairo")

  options(digits=16)
  int.factor=1*10^5 # Number of x used to calc area under Gaussian (is not analytic) 
  scale=2 # Initial value used to estimate scaling parameter
  width=1024
  height=768
  
  if (scanmode=="negative"){
    # pklist$neg
    # pklist$breaksFwhm
    findPeaks.Gauss.HPC(sum_neg, breaks.fwhm, int.factor, scale, resol, outdir, sampname, scanmode, FALSE, thresh, width, height)
  } else {
    findPeaks.Gauss.HPC(sum_pos, breaks.fwhm, int.factor, scale, resol, outdir, sampname, scanmode, FALSE, thresh, width, height)
  }
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
