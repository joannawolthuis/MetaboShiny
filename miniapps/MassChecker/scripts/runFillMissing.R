fillMissing <- function(file,type, scanmode, resol, outdir, cl=NULL, thresh, scripts) {
# file="./results/grouping_rest/negative_1.RData"
# file="./results/grouping_hmdb/1_negative.RData"
# scanmode= "negative"
# scripts="./scripts"
# resol=140000
# thresh=2000
# outdir="./results"
  print(paste0("filling peaks in file:", file))
  message(file.path(scripts, "AddOnFunctions/replaceZeros.R", sep="/"))
  #source(paste(scripts, "AddOnFunctions/replaceZeros.R", sep="/"))
  replaceZeros(file,type,scanmode,resol,outdir,cl,thresh,scripts)
}

# message("Start")
# cat("==> reading arguments:\n", sep = "")
# 
# cmd_args = commandArgs(trailingOnly = TRUE)
# 
# for (arg in cmd_args) cat("  ", arg, "\n", sep="")
# 
# run(cmd_args[1], cmd_args[2], as.numeric(cmd_args[3]), cmd_args[4], as.numeric(cmd_args[5]), cmd_args[6])
# 
# message("Ready")


