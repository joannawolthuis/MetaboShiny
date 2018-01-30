groupingRestBlur <- function(outdir, fileIn, cores=1, scanmode, ppm=2){
  # fileIn="./results/specpks_all_rest/negative_outlist_i_min_1.1.RData"
  # scanmode="negative"
  # outdir="./results"
  # ppm=2

  library(data.table)
  library(gsubfn)
  library(pbapply)
  library(xcms)
  library(meanShiftR)
  
  # fileIn = file.path(outdir, "specpks_all", "positive.RData")
  # scanmode = "positive"
  
  options(digits=16)
  load(fileIn)
  load(paste(outdir, "repl.pattern.RData", sep="/"))
  
  if (scanmode=="negative"){
    # repl.pattern=repl.pattern.neg
    groupNames=groupNames.neg
  } else {
    # repl.pattern=repl.pattern.pos
    groupNames=groupNames.pos
  }
  
  load(paste(outdir, "breaks.fwhm.RData", sep="/"))
  
  outpgrlist = NULL
  
  # Then group on highest peaks
  range = ppm*1e-06
  startcol=7
  set.seed( 2)
  # === MZCLUST ===
  
  outlist <- as.data.table(outlist)
  
  for (i in c(2:5)) {
    print(i)
    set(outlist, NULL, i, as.numeric(outlist[[i]]))
  }
  
  peakmat <- outlist[,c("mzmed.pkt", "height.pkt")]
  
  peakmat <- apply(peakmat, MARGIN = 2, as.numeric)
  mat <- t(peakmat)

  mzvals <- unique(outlist$mzmed.pkt)
  
  message("Grouping through mean shift...")
  
  d = 1:(length(mzvals))
  vec = split(d, 
              ceiling(seq_along(d)/2000))
  
  # --- parallel prep ---
  cl <<- makeSOCKcluster(cores,
                         outfile="~/mclog.txt")
  registerDoSNOW(cl)
  cfun <<- function(a, b) NULL
# ----------------------
  withProgress(message = "Peak grouping",{
    i = 1
    progress <<- function(i) setProgress(value = (i/length(vec)) / 2,
                                         detail = paste("Grouping through mean shift..."))
    opts <- list(progress=progress)
    # -------------------
    groups <- foreach(i=1:length(vec), .options.snow=opts,
                   .packages = "meanShiftR",
                   .verbose = T) %dopar% {
                     if(i < length(vec)){
                       range = c(min(vec[[i]]), max(vec[[i]]) + 1)
                     }else{
                       range = c(min(vec[[i]]), max(vec[[i]]))
                     }
                     a = meanShift(
                       mzvals[range[[1]]:range[[2]]],
                       mzvals[range[[1]]:range[[2]]]
                     )
                     a$assignment <- cbind(group = a$assignment, 
                                           mzidx = c(range[[1]]:range[[2]]))
                     # --- return ---
                     a$assignment
                   }
    # --- prep for next loop ---
    clustering = groups[[1]]
    print(clustering)
    colnames(clustering) = c("group", "mzidx")
    for(i in 2:length(groups)){
      try({
        group <- groups[[i]]
        colnames(group) <- c("group", "mzidx")
        overlap <- merge(clustering, 
                         group, 
                         by='mzidx')
        shift = overlap[,2] - overlap[,3]
        group[,'group'] = group[,'group'] + shift
        clustering = rbind(clustering[1:nrow(clustering)-1,], group)
      })
    }
    group_count = max(clustering[,"group"])
    clustering = as.data.table(clustering)
    # --- get intensities --
    outlist = as.data.table(outlist)
    # --------------------------
    progress <<- function(i) setProgress(value = 0.5 + ((i /group_count) / 2),
                                         detail = paste("Collecting values per group..."))
    opts <- list(progress=progress)
    # --- loop 2 ---
    outrows <- foreach(i=1:group_count, .options.snow=opts, .export = c("outlist", 
                                                                    "clustering", 
                                                                    "mzvals",
                                                                    "groupNames"),
                   .packages = "data.table",
                   .verbose = T) %dopar% {
                     mzidx <- clustering[group == i,]$mzidx
                     # FIND PEAKS IN THIS GROUP
                     mzs <- mzvals[mzidx]
                     members <- outlist[mzmed.pkt %in% mzs,]
                     # --- get intensities ---
                     mzmed = as.numeric(mean(members$mzmed.pkt))
                     mzmin = as.numeric(min(members$mzmin.pkt))
                     mzmax = as.numeric(max(members$mzmax.pkt))
                     fq.worst.pgrp = as.numeric(max(members$fq))
                     fq.best.pgrp = as.numeric(min(members$fq))
                     ints.allsamps = rep(0, length(groupNames))
                     names(ints.allsamps) = groupNames # same order as sample list!!!
                     # # Check for each sample if multiple peaks exists, if so take the sum!
                     labels=unique(members$samplenr)
                     nrsamples <- length(labels)
                     ints.allsamps[labels] = as.vector(unlist(lapply(labels, 
                                                                     function(x){
                                                                       sum(as.numeric(unlist(members[samplenr==x, 'height.pkt'])))
                                                                     })))
                     # --- make dt ---
                     outpgrlist = cbind('mzmed.pgrp' = mzmed,
                                        "fq.best"=fq.best.pgrp, 
                                        "fq.worst"=fq.worst.pgrp, 
                                        nrsamples, 
                                        'mzmin.pgrp' = mzmin, 
                                        'mzmax.pgrp' = mzmax,          
                                        t(as.matrix(ints.allsamps)
                                        ))
                     print(outpgrlist)
                     # -----------------------
                     as.data.table(outpgrlist)
                   }
  })
  stopCluster(cl)
  # --- finish ---
  outpgrlist <- rbindlist(outrows)
  colnames(outpgrlist)[1:6] = c("mzmed.pgrp", "fq.best", "fq.worst", "nrsamples", "mzmin.pgrp", "mzmax.pgrp")
  # --- save result ---
  dir.create(paste(outdir, "grouping_rest", sep="/"), showWarnings = FALSE)
  save(outpgrlist, file=paste(outdir, "grouping_rest", paste(scanmode,".RData", sep=""), sep="/"))
}
