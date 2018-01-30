groupingOnHMDB <- function(outdir, fileIn, scanmode, ppm=2) {
  # fileIn="./results/hmdb_part/negative_hmdb.97.RData"
  # scanmode="negative"
  # outdir="./results"
  # ppm=2
  print(scanmode)
  print(fileIn) 
  options(digits=16)
  load(fileIn) # outlist_part(HMDB_add_iso)
  HMDB_add_iso = outlist_part
  batch = strsplit(fileIn, ".",fixed = TRUE)[[1]][3]
  peaks <- paste(scanmode, "RData", sep=".")
  print(peaks)
  load(paste(outdir, "specpks_all", peaks, sep="/"))
  outlist.copy = outlist.tot
  rm(outlist.tot)
  
  load(paste(outdir, "repl.pattern.RData", sep="/"))
  load(paste(outdir, "breaks.fwhm.RData", sep="/"))
  
  outpgrlist = NULL
  
  if (scanmode=="negative"){
    label = "mz"
    # HMDB_add_iso=HMDB_add_iso.Neg
    repl.pattern=repl.pattern.neg
    groupNames=groupNames.neg
    
  } else {
    label = "mz"
    # HMDB_add_iso=HMDB_add_iso.Pos
    repl.pattern=repl.pattern.pos
    groupNames=groupNames.pos
  }
  
  # n=1
  
  # mass=98.9718582290911 
  # mass.all = as.numeric(HMDB_add_iso[,label])
  # index = which((mass.all > (mass - mtol)) & (mass.all < (mass + mtol)))
  
  # First group on HMDB masses
  while (dim(HMDB_add_iso)[1] > 0) {  
    index = 1
    # message(HMDB_add_iso[index,"CompoundName"])
    mass = as.numeric(HMDB_add_iso[,label][index])
    print(mass)
    mtol = (mass*ppm)/10^6
    mzmed = as.numeric(outlist.copy[,"mzmed.pkt"])
    selp = which((mzmed > (mass - mtol)) & (mzmed < (mass + mtol)))
    tmplist = outlist.copy[selp,,drop=FALSE]
    nrsamples = length(selp)

    if (nrsamples > 0) {
      # message(paste("Bingo ",n))
      print("Found one or more matches!")
      
      mzmed.pgrp = mean(as.numeric(outlist.copy[selp, "mzmed.pkt"]))
      # mzmed.pgrp = median(as.numeric(outlist.copy[selp, "mzmed.pkt"]))
      mzmin.pgrp = mass - mtol
      mzmax.pgrp = mass + mtol
      
      fq.worst.pgrp = as.numeric(max(outlist.copy[selp, "fq"]))
      fq.best.pgrp = as.numeric(min(outlist.copy[selp, "fq"]))
      ints.allsamps = rep(0, length(groupNames))
      names(ints.allsamps) = groupNames # same order as sample list!!!
      
      # # Check for each sample if multiple peaks exists, if so take the sum!
      labels=unique(tmplist[,"samplenr"])
      ints.allsamps[labels] = as.vector(unlist(lapply(labels, function(x) {sum(as.numeric(tmplist[which(tmplist[ , "samplenr"]==x), "height.pkt"]))})))
      # ints.allsamps[labels] = as.vector(unlist(lapply(labels, function(x) {as.numeric(tmplist[which(tmplist[ , "samplenr"]==x), "height.pkt"])})))
 
      outpgrlist = rbind(outpgrlist, 
                         cbind(mzmed.pgrp, 
                               "fq.best"=fq.best.pgrp, 
                               "fq.worst"=fq.worst.pgrp, 
                               nrsamples, 
                               mzmin.pgrp, 
                               mzmax.pgrp,          
                               t(as.matrix(ints.allsamps)
                               ))
                         )
      outlist.copy[selp, "height.pkt"] = -1
    }
    
    HMDB_add_iso = HMDB_add_iso[-index,]
    
    # n=n+1
  }
  
  print(head(outpgrlist))
  
  dir.create(paste(outdir, "grouping_hmdb", sep="/"), showWarnings = FALSE)
  save(outpgrlist, file=paste(paste(outdir, "grouping_hmdb", sep="/"), paste(paste(batch, scanmode, sep="_"), "RData", sep="."), sep="/"))
  
  dir.create(paste(outdir, "grouping_hmdb_done", sep="/"), showWarnings = FALSE)
  save(outlist.copy, file=paste(paste(outdir, "grouping_hmdb_done", sep="/"), paste(paste(batch, scanmode, sep="_"), "RData", sep="."), sep="/"))
  
  # message(Sys.time())
}  
