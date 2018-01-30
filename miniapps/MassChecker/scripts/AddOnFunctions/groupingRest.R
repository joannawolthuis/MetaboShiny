groupingRest <- function(outdir, fileIn, cl=F, scanmode, ppm=2){
  # fileIn="./results/specpks_all_rest/negative_outlist_i_min_1.1.RData"
  # scanmode="negative"
  # outdir="./results"
  # ppm=2
  
  options(digits=16)
  load(fileIn)
  outlist.copy = as.data.frame(unident.peaks)
  colnames(outlist.copy)[which(colnames(outlist.copy) == "intensity")] <- "height.pkt"
  batch = strsplit(fileIn, ".",fixed = TRUE)[[1]][3]
  load(paste(outdir, "repl.pattern.RData", sep="/"))
  # load(paste(outdir, "breaks.fwhm.RData", sep="/"))
  
  outpgrlist = NULL
  
  if (scanmode=="negative"){
    # repl.pattern=repl.pattern.neg
    groupNames=groupNames.neg
  } else {
    # repl.pattern=repl.pattern.pos
    groupNames=groupNames.pos
  }
  
  # Then group on highest peaks
  range = ppm*1e-06
  startcol=7
  
  # while (max(as.numeric(outlist.copy[ , "height.pkt"])) > 0 ) {
  while (dim(outlist.copy)[1] > 0) {
    print(paste("Peaks left:", dim(outlist.copy)[1]))
    sel = which(as.numeric(outlist.copy[ , "height.pkt"]) == max(as.numeric(outlist.copy[ , "height.pkt"])))[1]
    
    # 3ppm range around max
    mzref = as.numeric(outlist.copy[sel, "mzmed.pkt"])
    pkmin = -(range*mzref - mzref)
    pkmax = 2*mzref-pkmin
    
    selp = as.numeric(outlist.copy[ , "mzmed.pkt"]) > pkmin & as.numeric(outlist.copy[ , "mzmed.pkt"]) < pkmax
    tmplist = outlist.copy[selp,,drop=FALSE]
    
    nrsamples = length(unique(tmplist[,"samplenr"]))
    if (nrsamples > 0) {
      
      mzmed.pgrp = mean(as.numeric(outlist.copy[selp, "mzmed.pkt"]))
      mzmin.pgrp = -(range*mzmed.pgrp - mzmed.pgrp)
      mzmax.pgrp = 2*mzmed.pgrp - mzmin.pgrp
      
      selp = as.numeric(outlist.copy[ , "mzmed.pkt"]) > mzmin.pgrp & as.numeric(outlist.copy[ , "mzmed.pkt"]) < mzmax.pgrp
      tmplist = outlist.copy[selp,,drop=FALSE]
      
      # remove used peaks!!!
      tmp = as.vector(which(tmplist[,"height.pkt"]==-1))
      if (length(tmp)>0) tmplist=tmplist[-tmp,,drop=FALSE]
      
      nrsamples = length(unique(tmplist[,"samplenr"]))
      
      fq.worst.pgrp = as.numeric(max(outlist.copy[selp, "fq"]))
      fq.best.pgrp = as.numeric(min(outlist.copy[selp, "fq"]))
      ints.allsamps = rep(0, length(groupNames))
      names(ints.allsamps) = groupNames # same order as sample list!!!  
      
      # Check for each sample if multiple peaks exists, if so take the sum! 
      labels=unique(tmplist[,"samplenr"])
      ints.allsamps[labels] = as.vector(unlist(lapply(labels, function(x) {sum(as.numeric(tmplist[which(tmplist[ , "samplenr"]==x), "height.pkt"]))})))
      
      outpgrlist = rbind(outpgrlist, c(mzmed.pgrp, fq.best.pgrp, fq.worst.pgrp, nrsamples, mzmin.pgrp, mzmax.pgrp, ints.allsamps))#,NA,NA,NA,NA))
    }
    # outlist.copy[selp, "height.pkt"] = -1
    outlist.copy = outlist.copy[-which(selp==TRUE),,drop=FALSE]
  }

  
  outpgrlist = as.data.frame(outpgrlist)  # ignore warnings of duplicate row names
  colnames(outpgrlist)[1:6] = c("mzmed.pgrp", "fq.best", "fq.worst", "nrsamples", "mzmin.pgrp", "mzmax.pgrp")
  #colnames(outpgrlist)[(length(groupNames)+7):ncol(outpgrlist)] = c("assi_HMDB", "iso_HMDB", "HMDB_code", "theormz_HMDB")
  
  #save(outpgrlist_part, file=paste(outdir, paste(scanmode, "_", mzstart, "_", mzend, ".RData", sep=""), sep="/"))
  # save(final.outlist.filt, file=paste(outdir, "peak_grouping", paste(scanmode, "_",batch,".RData", sep=""), sep="/"))
  save(outpgrlist, file=paste(outdir, "grouping_rest", paste(scanmode,".RData", sep=""), sep="/"))

}
