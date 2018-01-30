collectSamples <- function(resultDir, scanmode){
  # resultDir="./results"
  # scanmode="negative"
  # scanmode="positive"
  
  ppm=2
  
  # Check if all jobs terminated correct!
  notRun = NULL
  
  load(paste(resultDir, "init.RData", sep="/"))
  
  object.files = list.files(paste(resultDir, "specpks", sep="/"), full.names=TRUE, pattern="*.RData")
  
  for (i in 1:length(groupNames)) {
    if (!(paste(resultDir, "specpks", paste(paste(groupNames[i],scanmode,sep="_"), ".RData", sep=""), sep="/") %in% object.files)) {
      notRun = c(notRun, paste(resultDir, "specpks", paste(paste(groupNames[i],scanmode,sep="_"), ".RData", sep=""), sep="/"))
    }  
  }
  
  #if (is.null(notRun)){
  message("Collecting samples!")
  
  # negative
  filepath =  paste(resultDir, "specpks", sep="/")
  files = list.files(filepath,recursive=TRUE, full.names=TRUE, pattern=paste("*_",scanmode,".RData",sep=""))
  
  outlist=NULL
  for (i in 1:length(files)) {
    #message(files[i])
    load(files[i])
    
    if (is.null(outlist.persample) || (dim(outlist.persample)[1]==0)){
      tmp=strsplit(files[i], "/")[[1]]
      fname = tmp[length(tmp)]
      #fname = strsplit(files[i], "/")[[1]][8]
      fname = strsplit(fname, ".RData")[[1]][1]
      fname = substr(fname, 13, nchar(fname))
      
      if (i == 1) { outlist <- c(fname, rep("-1",5)) } else { outlist <- rbind(outlist, c(fname, rep("-1",5)))}
    } else {
      if (i == 1) { outlist <- outlist.persample } else { outlist <- rbind(outlist, outlist.persample)}
    }
  }
  
  # remove negative values
  index=which(outlist[,"height.pkt"]<=0)
  if (length(index)>0) outlist = outlist[-index,]
  index=which(outlist[,"mzmed.pkt"]<=0)
  if (length(index)>0) outlist = outlist[-index,]
  
  outdir=paste(resultDir, "specpks_all", sep="/")
  dir.create(outdir)
  save(outlist, file=paste(outdir, paste(scanmode, "RData", sep="."), sep="/"))
  
}
# 
# message("Start")
# cat("==> reading arguments:\n", sep = "")
# 
# cmd_args = commandArgs(trailingOnly = TRUE)
# 
# for (arg in cmd_args) cat("  ", arg, "\n", sep="")
# 
# run(cmd_args[1], cmd_args[2])
# #run("./results")
# 
# message("Ready")