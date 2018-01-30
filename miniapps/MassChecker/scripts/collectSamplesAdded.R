collectSamplesAdded <- function(resultDir, scanmode){
# resultDir="./results"
# scanmode="negative"

  object.files = list.files(paste(resultDir, "adductSums", sep="/"), full.names=TRUE, pattern=paste(scanmode, "_", sep=""))
  
  outlist.tot=NULL
  for (i in 1:length(object.files)) {
    load(object.files[i])
    outlist.tot = rbind(outlist.tot, adductsum)
  }
  
  save(outlist.tot, file=paste(resultDir, "/adductSums_", scanmode, ".RData", sep=""))

}

# message("Start")
# cat("==> reading arguments:\n", sep = "")
# 
# cmd_args = commandArgs(trailingOnly = TRUE)
# 
# for (arg in cmd_args) cat("  ", arg, "\n", sep="")
# 
# run(cmd_args[1], cmd_args[2])
# 
# message("Ready")