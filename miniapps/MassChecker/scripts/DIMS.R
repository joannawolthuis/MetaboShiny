DIMS <- function(xmlfile, scripts, outdir, trim, dimsThresh, resol){
  # dimsThresh=100
  # nrepl=3
  # xmlfile="./data/BSP20150716_58.mzXML"
  # scripts="./scripts"
  # outdir="./results"
  # trim=0.1
  # resol=140000
  
  dir.create(outdir,showWarnings = F)
  dir.create(paste(outdir, "pklist", sep="/"),showWarnings = F)
  dir.create(paste(outdir, "QC", sep="/"),showWarnings = F)
  
  sampname <- strsplit(xmlfile, "/")[[1]]
  sampname <- sampname[length(sampname)]
  sampname <- strsplit(sampname, ".mzXML")[[1]][1]
  
  
  options(digits=16)
  
  ### process one sample at a time and find peaks FOR BOTH SCAN MODES! #
  int.factor=1*10^5 # Number of x used to calc area under Gaussian (is not analytic) 
  scale=2 # Initial value used to estimate scaling parameter
  width=1024
  height=768
  
  if (!file.exists(paste(paste(outdir, "pklist", sep="/"),"/",sampname, ".RData", sep=""))){
    # Aggregate with dims scipt
    pklist = dims(xmlfile, outdir, dimsThresh, trim, resol)
    #save(pklist, file=paste(paste(outdir, "pklist", sep="/"),"/", sampname, ".RData", sep=""))
    save(pklist, file=paste(paste(outdir, "pklist", sep="/"),"/", sampname, ".RData", sep=""))
    
  }
}

