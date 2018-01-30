peakGrouping.2.0 <- function(fileIn, scriptDir, outdir, resol, scanmode){
  # rdata="./results/specpks_all/positive_outlist_i_min_1.197.RData"
  # scriptDir="./scripts"
  # outdir="./results"
  # inpdir="./data"
  # thresh=2000
  # resol=140000
  # scanmode="positive"

  # dir.create(outdir,showWarnings = F)
  # dir.create(paste(outdir, "peak_grouping", sep="/"),showWarnings = F)
  
  options(digits=16)
  
  source(paste(scriptDir, "AddOnFunctions/sourceDir.R", sep="/"))
  sourceDir(paste(scriptDir, "AddOnFunctions", sep="/"))
  
  message(paste("File to group:", fileIn))
  print(scanmode) 
  # load(paste(outdir, "repl.pattern.RData", sep="/"))
  # if (scanmode=="negative") {
  #   sampleNames=groupNames.neg  
  # } else {
  #   sampleNames=groupNames.pos  
  # }
  # 
  # peak.grouping.Gauss.HPC(outdir, fileIn, scanmode, resol, sampleNames)
  # groupingAndIdent(outdir, fileIn, scanmode)
  groupingOnHMDB(outdir, fileIn, scanmode)

}

# message("Start")
# cat("==> reading arguments:\n", sep = "")
# 
# cmd_args = commandArgs(trailingOnly = TRUE)
# 
# for (arg in cmd_args) cat("  ", arg, "\n", sep="")
# 
# run(cmd_args[1], cmd_args[2], cmd_args[3], as.numeric(cmd_args[4]), cmd_args[5])
# 
# message("Ready")
