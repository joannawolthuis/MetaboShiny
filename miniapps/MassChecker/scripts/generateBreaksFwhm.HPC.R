generateBreaksFwhm <- function(file,outdir,trim,resol,nrepl){
  # file="./data/BSP20161128_31.mzXML"
  
  library("xcms")
  withProgress(message = "Initial setup (Step 1)",{
    
    setProgress(0.2, detail="Creating directory...")
    
    trimLeft=NULL
    trimRight=NULL
    breaks.fwhm=NULL
    breaks.fwhm.avg=NULL
    bins=NULL
    posRes=NULL
    negRes=NULL
    
    dir.create(outdir, showWarnings = F)
    x = NULL
    try({x = xcmsRaw(file)}, silent = TRUE)
    if (is.null(x)){
      return(NULL)
    }
    
    setProgress(0.4, detail="Finding ranges...")
    
    trimLeft = round(x@scantime[length(x@scantime)*trim])  
    trimRight = round(x@scantime[length(x@scantime)*(1-trim)])
    message(paste("trimLeft", trimLeft, sep=" "))
    message(paste("trimRight", trimRight, sep=" "))
    
    # Mass range m/z
    lowMZ = x@mzrange[1]
    highMZ = x@mzrange[2]
    message(paste("lowMZ", lowMZ, sep=" "))
    message(paste("highMZ", highMZ, sep=" "))
    
    setProgress(0.6, detail = "Generating breaks...")
    
    # breaks.fwhm <- seq(from=lowMZ, to=highMZ, by=deltaMZ)
    # breaks has fixed distance between min and max of a bin.
    # better if this distance were a function of fwhm=f(mz)
    #segment <- seq(from=lowMZ, to=highMZ, length.out=1001)
    nsegment = 2*(highMZ-lowMZ)
    segment = seq(from=lowMZ, to=highMZ, length.out=nsegment+1)
    breaks.fwhm=NULL
    breaks.fwhm.avg=NULL
    # for (i in 1:2) {
    for (i in 1:nsegment) {
      startsegm <- segment[i]
      endsegm <- segment[i+1]
      resol.mz <- resol*(1/sqrt(2)^(log2(startsegm/200)))
      fwhmsegm <- startsegm/resol.mz
      breaks.fwhm <- c(breaks.fwhm, seq(from=(startsegm+fwhmsegm),to=endsegm, by=0.2*fwhmsegm))
      #breaks.fwhm <- c(breaks.fwhm, seq(from=(startsegm), to=endsegm, by=0.2*fwhmsegm))
      
      # average the m/z instead of start value
      range = seq(from=(startsegm+fwhmsegm),to=endsegm, by=0.2*fwhmsegm)
      deltaMZ = range[2]-range[1] 
      breaks.fwhm.avg <- c(breaks.fwhm.avg, range + 0.5 * deltaMZ)
      
    }
    
    setProgress(0.8, detail="Saving breaks...")
    
    save(breaks.fwhm,breaks.fwhm.avg,trimLeft,trimRight,file=paste(outdir, "breaks.fwhm.RData", sep="/"))
    
    tbl=read.table(file=paste(outdir,"..","sampleNames.txt" ,sep = "/"), header = TRUE, sep="\t")  
    sampleNames=as.vector(unlist(tbl$File_Name))
    groupNames=unique(as.vector(unlist(tbl$Sample_Name)))
    # groupNames=unique(unlist(lapply(strsplit(as.vector(unlist(tbl$Sample_Name)), ".", fixed = TRUE), function(x) x[1])))
    # save(tbl, file=paste(outdir, "sampleNames.RData", sep="/"))
    
    setProgress(0.8, detail="Saving sample info...")
    
    ###############################################################################################
    ####################### experimental design ###################################################
    ###############################################################################################
    # nsampgrps = number of individual biological samples
    # nrepl = number of technical replicates per sample
    nsampgrps = length(sampleNames)/nrepl
    repl.pattern = NULL
    if (nrepl == 3){
      for (x in 1:nsampgrps) { repl.pattern <- c(repl.pattern, list(c(sampleNames[x*nrepl-2],sampleNames[x*nrepl-1],sampleNames[x*nrepl])))}
    } else if (nrepl == 5){
      for (x in 1:nsampgrps) { repl.pattern <- c(repl.pattern, list(c(sampleNames[x*nrepl-4],sampleNames[x*nrepl-3],sampleNames[x*nrepl-2],sampleNames[x*nrepl-1],sampleNames[x*nrepl])))}
    }
    
    save(nsampgrps, repl.pattern, groupNames, sampleNames, file=paste(outdir, "init.RData", sep="/")) 
    ###############################################################################################
    ###############################################################################################
  })
}


# message("Start")
# cat("==> reading arguments:\n", sep = "")
# 
# cmd_args = commandArgs(trailingOnly = TRUE)
# 
# for (arg in cmd_args) cat("  ", arg, "\n", sep="")
# 
# run(cmd_args[1], cmd_args[2], as.numeric(cmd_args[3]), as.numeric(cmd_args[4]), as.numeric(cmd_args[5]))
# 
# message("Ready")