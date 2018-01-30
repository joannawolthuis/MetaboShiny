collectSamplesFilled <- function(resultDir, scanmode, scripts){
# resultDir="./results"
# scanmode="positive"
# scripts="./scripts"  
  
  object.files = list.files(paste(resultDir, "samplePeaksFilled", sep="/"), full.names=TRUE, pattern=scanmode)

  outlist.tot=NULL
  for (i in 1:length(object.files)) {
    load(object.files[i])
    outlist.tot = rbind(outlist.tot, final.outlist.idpat3)
  }
  
  source(paste(scripts, "AddOnFunctions/sourceDir.R", sep="/"))
  sourceDir(paste(scripts, "AddOnFunctions", sep="/"))

  # remove duplicates
  outlist.tot = mergeDuplicatedRows(outlist.tot)
  
  # sort on mass
  outlist.tot = outlist.tot[order(outlist.tot[,"mzmed.pgrp"]),]
  
  save(outlist.tot, file=paste(resultDir, "/outlist2identify_", scanmode, ".RData", sep=""))
  
  # cut HMDB ##########################################################################################################################################

  outdir=paste(resultDir, "hmdb_part_adductSums", sep="/")
  dir.create(outdir, showWarnings = FALSE)
  load(paste(resultDir, "breaks.fwhm.RData", sep="/"))
  
  #if (scanmode=="negative"){
  #  label = "MNeg"
  #  HMDB_add_iso=HMDB_add_iso.Neg
  #} else {
  #  label = "Mpos"
  #  HMDB_add_iso=HMDB_add_iso.Pos
  #}
    
  label = "mz"
    
  if (scanmode=="negative"){
    load(paste(resultDir, "../db/negative.RData", sep="/"))
    HMDB_add_iso=references
  } else {
    load(paste(resultDir, "../db/positive.RData", sep="/"))
    HMDB_add_iso=references
  }
  
  # filter mass range meassured!!!
  HMDB_add_iso = HMDB_add_iso[which(HMDB_add_iso[,label]>=breaks.fwhm[1] & HMDB_add_iso[,label]<=breaks.fwhm[length(breaks.fwhm)]),]
  
  outlist=HMDB_add_iso
  
  #outlist=outlist[grep("HMDB",rownames(outlist),fixed=TRUE),]
  #outlist=outlist[-grep("_",rownames(outlist),fixed=TRUE),]
  
  n=dim(outlist)[1]
  sub=300
  end=0
  check=0
  
  if (n>=sub & (floor(n/sub))>=2){
    for (i in 1:floor(n/sub)){
      start=-(sub-1)+i*sub
      end=i*sub
      
      # message(paste("start:",start))
      # message(paste("end:",end))
      
      outlist_part = outlist[c(start:end),]
      save(outlist_part, file=paste(outdir, paste(scanmode, paste("hmdb",i,"RData", sep="."), sep="_"), sep="/"))
    }
  } 
  
  start = end + 1
  end = n
  
  # message(paste("start:",start))
  # message(paste("end:",end))
  
  outlist_part = outlist[c(start:end),]
  save(outlist_part, file=paste(outdir, paste(scanmode, paste("hmdb",i+1,"RData", sep="."), sep="_"), sep="/"))
  
}

# message("Start")
# cat("==> reading arguments:\n", sep = "")
# 
# cmd_args = commandArgs(trailingOnly = TRUE)
# 
# for (arg in cmd_args) cat("  ", arg, "\n", sep="")
# 
# run(cmd_args[1], cmd_args[2], cmd_args[3])
# 
# message("Ready")