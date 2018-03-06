groupingRestClust_2 <- function(outdir, 
                              cores=1,
                              ppm=2){
  
  library(data.table)
  library(gsubfn)
  library(DBI)
  library(RSQLite)
  library(pbapply)
  library(xcms)
  
  options(digits=16)

  # --- cluster ---
  cl <<- makeSOCKcluster(cores,
                         outfile="~/mclog.txt")
  registerDoSNOW(cl)
  
  j = 0
  
  for(mode in c("pos", "neg")){

    fn = file.path(outdir, "collected", paste0("peaks_",
                                                mode,
                                                ".RData"))
    load(fn)

    peakmat <- data.table(mzmed = outpgrlist$mzmed)
    peakmat <- apply(peakmat, MARGIN = 2, as.numeric)
    groups <- mzClustGeneric(as.matrix(peakmat), 
                             mzppm=ppm,
                             shinyprog=FALSE)
    outlist <- outpgrlist
    groupNames <- colnames(outlist[,-1])
    # --- cl prep ---
    progress <<- function(i) setProgress(value = j + (i / length(groups$idx)) / 2)
    opts <- list(progress=progress)
    # ---------
    outrows <- foreach(i=1:length(groups$idx), .options.snow=opts, .export = c("outdir",
                                                                          "scriptdir",
                                                                          "resol",
                                                                          "thresh_pos",
                                                                          "thresh_neg"),
            .verbose = FALSE,
            .packages = c("MassSpecWavelet", "data.table")) %dopar% {
              mzidx <- groups$idx[[i]]
              # FIND PEAKS IN THIS GROUP
              members <- outlist[mzidx, ]
              # --- get intensities ---
              mzmed = as.numeric(mean(members$mzmed))
              mzmin = as.numeric(min(members$mzmed))
              mzmax = as.numeric(max(members$mzmed))
              
              ints.allsamps = rep(0, length(groupNames))
              names(ints.allsamps) = groupNames # same order as sample list!!!
              # # Check for each sample if multiple peaks exists, if so take the sum!
              nrsamples <- nrow(members)
              ints.allsamps <- colSums(outlist[mzidx, -1],
                                       na.rm = T)
              #print(ints.allsamps)
              # --- make dt ---
              outpgrlist.a = cbind('mzmed' = mzmed, 
                                   nrsamples, 
                                   'mzmin' = mzmin, 
                                   'mzmax' = mzmax
              )
              outpgrlist.b <- data.frame(as.list(ints.allsamps))
              outpgrlist <- cbind(outpgrlist.a, outpgrlist.b)

              # -----------------------
              
              as.data.table(outpgrlist)
            }
    
    outpgrlist <- rbindlist(outrows)
    
    print(dim(outpgrlist))
    
    f =  file.path(outdir, "grouped", paste0("grouped_",mode,".RData"))
    
    save(outpgrlist, file=f)
    
    j = j + 0.5
  }
  stopCluster(cl)
}
