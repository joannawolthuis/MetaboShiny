sumAdducts <- function(file, scanmode, outdir, adducts) {
# file="./results/hmdb_part/negative_hmdb.1.RData"
# scanmode="negative"
# outdir="./results"
# adducts="1,2"
  
  load(paste(outdir, "repl.pattern.RData", sep="/"))
  
  if (scanmode=="negative"){
    groupNames=groupNames.neg
  } else {
    groupNames=groupNames.pos
  }

  adducts=as.vector(unlist(strsplit(adducts, ",",fixed = TRUE)))
  
  load(file)
  load(paste(outdir, "/outlist2identify_", scanmode, ".RData", sep=""))
  batch = strsplit(file, ".",fixed = TRUE)[[1]][3]
  
outlist.tot=unique(outlist.tot)
  
  source(paste(outdir, "../scripts/AddOnFunctions/sumAdducts.R", sep="/"))
  sumAdducts(outlist.tot, outlist_part, groupNames, adducts, batch, scanmode, outdir)
}
# 
# message("Start")
# cat("==> reading arguments:\n", sep = "")
# 
# cmd_args = commandArgs(trailingOnly = TRUE)
# 
# for (arg in cmd_args) cat("  ", arg, "\n", sep="")
# 
# run(cmd_args[1], cmd_args[2], cmd_args[3], cmd_args[4])
# 
# message("Ready")


