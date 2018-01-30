replaceZeros <- function(file,type,scanmode,resol,outdir,cl,thresh,scriptDir){
# file="./results/grouping_rest/negative_1.RData"
# scanmode= "negative"
# scriptDir="./scripts"
# resol=140000
# thresh=2000
# outdir="./results"
  library(pbapply)
  library(data.table)
  
  control_label="C"
  
  #source(paste(scriptDir, "AddOnFunctions/sourceDir.R", sep="/"))
  #sourceDir(paste(scriptDir, "AddOnFunctions", sep="/"))
  
  dir.create(paste(outdir, "samplePeaksFilled", sep="/"), showWarnings = FALSE)

  # int.factor=1*10^5 # Number of x used to calc area under Gaussian (is not analytic) 
  # scale=2 # Initial value used to estimate scaling parameter
  # width=1024
  # height=768
  
  message(paste("file", file))
  message(paste("scanmode", scanmode))
  message(paste("resol", resol))
  message(paste("outdir", outdir))
  message(paste("thresh", thresh))
  message(paste("scriptDir", scriptDir))
  
  load(paste(outdir, "repl.pattern.RData", sep="/"))
  
  if (scanmode=="negative"){
    # repl.pattern=repl.pattern.neg
    groupNames=groupNames.neg
  } else {
    # repl.pattern=repl.pattern.pos
    groupNames=groupNames.pos
  }

  name = as.vector(unlist(strsplit(file, "/", fixed=TRUE)))
  name = name[length(name)]
  name = paste0(scanmode,"_",type,".RData")
  message(paste("File name: ", name))
  
  # load samplePeaks
  # load  outpgrlist
  load(file)  
  outpgrlist <- as.data.frame(outpgrlist)
  
  # #################################################################################

    outcols <- pblapply(groupNames, 
                      cl=cl, 
                      FUN=function(gname){
    message(gname)
    samplePeaks=outpgrlist[,gname]
    index=which(samplePeaks<=0)
    for (j in 1:length(index)){
      area = generateGaussian(outpgrlist[index[j],"mzmed.pgrp"],
                              thresh,
                              resol,
                              FALSE,
                              scanmode,
                              int.factor=1*10^5,1,1)$area
      print(area)
      # area = area/2
      outpgrlist[index[j],gname] = rnorm(n=1, 
                                         mean=area, 
                                         sd=0.25*area)
    }
    # --- return ---
    as.data.frame(outpgrlist[,gname])
  })
  patcols <- do.call(cbind, outcols)
  outpgrlist <- cbind(outpgrlist[,1:6],
                      patcols)
  ################################################################################
  
  # Add average column
  outpgrlist = cbind(outpgrlist, "avg.int"=apply(outpgrlist[, 7:(ncol(outpgrlist)-4)], 1, mean))

###############################################
  
  message(paste("File saved: ", paste(outdir, "/samplePeaksFilled/", name, sep="")))
  save(outpgrlist, file=paste(outdir, "/samplePeaksFilled/", name, sep=""))
}



replaceZeros_lookup <- function(file,type,scanmode,resol,outdir,cores=1,thresh,scriptDir){
  # file="./results/grouping_rest/negative_1.RData"
  # scanmode= "negative"
  # scriptDir="./scripts"
  # resol=140000
  # thresh=2000
  # outdir="./results"
  library(pbapply)
  library(data.table)
  
  control_label="C"
  
  #source(paste(scriptDir, "AddOnFunctions/sourceDir.R", sep="/"))
  #sourceDir(paste(scriptDir, "AddOnFunctions", sep="/"))
  
  dir.create(paste(outdir, 
                   "samplePeaksFilled", 
                   sep="/"), 
             showWarnings = FALSE)
  
  # int.factor=1*10^5 # Number of x used to calc area under Gaussian (is not analytic) 
  # scale=2 # Initial value used to estimate scaling parameter
  # width=1024
  # height=768
  
  message(paste("file", file))
  message(paste("scanmode", scanmode))
  message(paste("resol", resol))
  message(paste("outdir", outdir))
  message(paste("thresh", thresh))
  message(paste("scriptDir", scriptDir))
  
  load(paste(outdir, "repl.pattern.RData", sep="/"))
  
  if (scanmode=="negative"){
    # repl.pattern=repl.pattern.neg
    groupNames=groupNames.neg
  } else {
    # repl.pattern=repl.pattern.pos
    groupNames=groupNames.pos
  }
  
  name = as.vector(unlist(strsplit(file, "/", fixed=TRUE)))
  name = name[length(name)]
  name = paste0(scanmode,"_",type,".RData")
  message(paste("File name: ", name))
  
  file = file.path(outdir, "grouping_rest", paste0(scanmode, ".RData"))
  # load samplePeaks
  # load  outpgrlist
  load(file)  
  
  outpgrlist <- as.data.table(outpgrlist)
  
  # #################################################################################
  
  mzvals = unique(outpgrlist$mzmed.pgrp)
  
  cl <<- makeSOCKcluster(cores,
                         outfile="~/mclog.txt")
  registerDoSNOW(cl)
  cfun <<- function(a, b) NULL
  # ---------------
  withProgress(message = "Filling missing values...",{
    progress <<- function(i) setProgress(value = (i / length(mzvals))/2,
                                         detail = "Generating reference values...")
    opts <- list(progress=progress)
    cfun <<- function(a, b) c(a, b)
    gaussians <- foreach(i=1:length(mzvals), .options.snow=opts, .export = c("mzvals", 
                                                                      "thresh", 
                                                                      "resol", 
                                                                      "scanmode", 
                                                                      "searchMZRange",
                                                                      "fitGaussianInit",
                                                                      "generateGaussian",
                                                                      "fitGaussian",
                                                                      "getFwhm",
                                                                      "getSD",
                                                                      "optimizeGauss",
                                                                      "getArea",
                                                                      "fit2G_2",
                                                                      "fit4G_2",
                                                                      "fit4peaks",
                                                                      "fitG_2",
                                                                      "fit3G_2",
                                                                      "fit1Peak",
                                                                      "fit2peaks",
                                                                      "fit3peaks",
                                                                      "getFitQuality",
                                                                      "checkOverlap",
                                                                      "isWithinXppm",
                                                                      "sumCurves"),
                   .packages = "xcms",
                   .verbose = T,
                   .inorder = T, # !!! important !!!
                   .combine="c") %dopar% {
                     area = generateGaussian(mzvals[i],
                                             thresh,
                                             resol,
                                             FALSE,
                                             scanmode,
                                             int.factor=1*10^5,1,1)$area
                     area
                   }
    # --- prep ---
    ref_table <- data.table(mz = mzvals,
                            int = gaussians)
    print(head(ref_table))
    setkey(ref_table, mz)
    outpgrlist <- outpgrlist[,-"avg.int"]
    # --- loop 2 ---
    progress <<- function(i) setProgress(value = 0.5 + (i / nrow(outpgrlist))/2,
                                         detail = "Filling values...")
    opts <- list(progress=progress)
    outrows <- foreach(i=1:nrow(outpgrlist), .options.snow=opts, .export = c("outpgrlist", 
                                                                          "ref_table"),
                       .packages="data.table",
                   .verbose = T,
                   .inorder = T) %dopar% {
                     row = outpgrlist[i,]
                     rowstart = row[,c(1:6)]
                     gauss <- ref_table[row[,1],]
                     area = unlist(gauss$int)
                     # loop and fill
                     filled_row <- sapply(row[,7:length(row)], FUN=function(samp){
                       if(samp <= 0){
                         rand_gauss <- rnorm(n=1, 
                                             mean=area, 
                                             sd=0.25*area)
                         rand_gauss
                       } else{
                         samp}
                     })
                     res <- cbind(rowstart, 
                                  t(filled_row), 
                                  avg.int=mean(filled_row))
                     # # --- return ---
                     as.data.table(res)
                   }
  })
  # stop cluster
  stopCluster(cl)
  # combine
  outpgrlist <- rbindlist(outrows)
  
  ###############################################
  
  message(paste("File saved: ", paste(outdir, "/samplePeaksFilled/", name, sep="")))
  save(outpgrlist, file=paste(outdir, "/samplePeaksFilled/", name, sep=""))
}

