peak.grouping.Gauss.HPC <- function(outdir, fileIn, scanmode, resol, groupNames) {
# fileIn="./results/specpks_all/positive_outlist_i_min_1.197.RData"

  # ppm/2
  #range = 1.5e-06
  range = 2e-06  
  startcol=7

  # outlist.copy <- read.table(file=fileIn, sep="\t", header=TRUE)
  load(fileIn)
  outlist.copy = outlist_i_min_1
  batch = strsplit(fileIn, ".",fixed = TRUE)[[1]][3]

  outpgrlist = NULL
  
  while (max(as.numeric(outlist.copy[ , "height.pkt"])) > 0 ) {
    
    sel = which(as.numeric(outlist.copy[ , "height.pkt"]) == max(as.numeric(outlist.copy[ , "height.pkt"])))[1]
    
    # 3ppm range around max
    mzref = as.numeric(outlist.copy[sel, "mzmed.pkt"])
    pkmin = -(range*mzref - mzref)
    pkmax = 2*mzref-pkmin
    
    selp = as.numeric(outlist.copy[ , "mzmed.pkt"]) > pkmin & as.numeric(outlist.copy[ , "mzmed.pkt"]) < pkmax
    tmplist = outlist.copy[selp, ]
    
    nrsamples = sum(selp)
    if (nrsamples > 1) {
      # 3ppm range around mean
      mzmed.pgrp = mean(as.numeric(outlist.copy[selp, "mzmed.pkt"]))
      mzmin.pgrp = -(range*mzmed.pgrp - mzmed.pgrp)
      mzmax.pgrp = 2*mzmed.pgrp - mzmin.pgrp
      
      selp = as.numeric(outlist.copy[ , "mzmed.pkt"]) > mzmin.pgrp & as.numeric(outlist.copy[ , "mzmed.pkt"]) < mzmax.pgrp
      tmplist = outlist.copy[selp, ]
      
      fq.worst.pgrp = as.numeric(max(outlist.copy[selp, "fq"]))
      fq.best.pgrp = as.numeric(min(outlist.copy[selp, "fq"]))
      ints.allsamps = rep(0, length(groupNames))
      names(ints.allsamps) = groupNames # same order as sample list!!!  
      
      #       if (length(unique(as.vector(tmplist[,"samplenr"]))) != length(as.vector(tmplist[,"samplenr"]))) {
      #         message(paste("bingo", sel))
      #         break
      #       }
      
      # Check for each sample if multiple peaks exists, if so take the sum! 
      labels=unique(tmplist[,"samplenr"])
      ints.allsamps[labels] = as.vector(unlist(lapply(labels, function(x) {sum(as.numeric(tmplist[which(tmplist[ , "samplenr"]==x), "height.pkt"]))})))
      
      outpgrlist = rbind(outpgrlist, c(mzmed.pgrp, fq.best.pgrp, fq.worst.pgrp, nrsamples, mzmin.pgrp, mzmax.pgrp, ints.allsamps))
    }
    outlist.copy[selp, "height.pkt"] = -1
  }

  outpgrlist = as.data.frame(outpgrlist)  # ignore warnings of duplicate row names
  colnames(outpgrlist)[1:6] = c("mzmed.pgrp", "fq.best", "fq.worst", "nrsamples", "mzmin.pgrp", "mzmax.pgrp")
  
#   # filtering ##################################################################################################################
#   final.outlist=outpgrlist[,c("mzmed.pgrp", "fq.best", "fq.worst","nrsamples","mzmin.pgrp","mzmax.pgrp", sampleNames)]
#   # NB: in centroided mode, data files contains many "-1.000" values, from peak finding. Set these to zero.
#   final.outlist[final.outlist == -1] = 0
#   
#   # keep only peaks which occur in 3 out of 3 technical replicates in at least one sample in peak group list
# #   peakFiltering(repl.pattern, final.outlist, nsampgrps, outdir, scanmode, startcol=7)
# #   peakFiltering <- function(repl.pattern, final.outlist, nsampgrps, resultDir, scanmode, startcol){
#   nsamp = length(repl.pattern)
#   nsampgrps = length(repl.pattern[[1]])
#   
#   keep <- rep(0, nrow(final.outlist))
#   for (p in 1:nrow(final.outlist)) {  
#     for (g in 1:nsampgrps) {  # g <- 31
#       if (keep[p] == 0 & sum(final.outlist[p, repl.pattern[[g]]] > 0) == length(repl.pattern[[g]]) ) { keep[p] <- 1 }
#     }
#   }
#   
#   tmp <- cbind(final.outlist, keep)
#   final.outlist.filt <- tmp[keep == 1, ]
# 
#   # omit keep column
#   final.outlist.filt <- final.outlist.filt[ , -ncol(final.outlist.filt)]

  #save(outpgrlist_part, file=paste(outdir, paste(scanmode, "_", mzstart, "_", mzend, ".RData", sep=""), sep="/"))
  # save(final.outlist.filt, file=paste(outdir, "peak_grouping", paste(scanmode, "_",batch,".RData", sep=""), sep="/"))
  save(outpgrlist, file=paste(outdir, "peak_grouping", paste(scanmode, "_",batch,".RData", sep=""), sep="/"))
  
}