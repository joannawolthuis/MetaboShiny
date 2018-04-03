run_step_1 <- function(){
  # -----------------------
 

run_step_2 <- function(){
  # make cluster obj
  
  dir.create(file.path(outdir, "spectra"),showWarnings = F)
  
  # --- cluster ---
  cl <<- makeSOCKcluster(cores,
                         outfile="~/mclog.txt")
  registerDoSNOW(cl)
  # -----------------
  filedir <<- path()
  files = list.files(filedir, 
                     pattern = "\\.raw",
                     full.names = T)
  cfun <<- function(a, b) NULL
  progress <<- function(i) setProgress(value = i / length(files),
                                       detail = paste("Done:", basename(files[i])))
  opts <- list(progress=progress)
  # ---------------
  withProgress(message = "DIMS data extraction",{
    i = 1
    foreach(i=1:length(files), .options.snow=opts, .export = c("dims_raw", 
                                                               "outdir", 
                                                               "scriptdir",
                                                               "dimsThresh",
                                                               "trim",
                                                               "resol"),
            .packages = c("MALDIquant","data.table","gsubfn"),
            .verbose = T,
            .combine=cfun) %dopar% {
              dims_raw(files[[i]],
                       outdir,
                       thresh=dimsThresh,
                       trim=trim,
                       resol=resol, 
                       scriptdir, 
                       mzmin=70, 
                       mzmax=600)
            }
  })
  stopCluster(cl)
  updateCollapse(session, "pipeline", 
                 style = list("DIMS data extraction" = "success"))
}

run_step_3 <- function(){
  
  dir.create(file.path(outdir, "averaged"),showWarnings = F)
  
  averageTechReplicates_jw(outdir, 
                           cores)
}

run_step_4 <- function(){
  
  dir.create(file.path(outdir, "peaks"),showWarnings = F)
  
  # --- cluster ---
  cl <<- makeSOCKcluster(cores,
                         outfile="~/mclog.txt")
  print(cl)
  registerDoSNOW(cl)
  # ---------------
  # PEAKFINDING
  files = list.files(file.path(outdir, 
                               "averaged"),
                     full.names = T)
  #,pattern= if(scanmode == "positive") "pos" else "neg") 
  # -----------------
  cfun <<- function(a, b) NULL
  progress <<- function(i) setProgress(value = i / length(files),
                                       detail = paste("Done:", basename(files[i])))
  opts <- list(progress=progress)
  # ---------------
             # --------------------------
             foreach(i=1:length(files), .options.snow=opts, .export = c("outdir",
                                                                        "scriptdir",
                                                                        "resol",
                                                                        "thresh_pos",
                                                                        "thresh_neg",
                                                                        "peakFinding.wavelet_2"),
                     .verbose = T,
                     .packages = c("MassSpecWavelet", "data.table"),
                     .combine=cfun) %dopar% {
                       scanmode = if(grepl(files[i],pattern = "_pos.RData")) "pos" else "neg"
                       peakFinding.wavelet_2(file = files[i], 
                                             scripts = scriptdir, 
                                             outdir = outdir, 
                                             thresh = if(scanmode == "pos") thresh_pos else{thresh_neg},
                                             sig_noise = 3)
                     }
}

run_step_5 <- function(){
  
  dir.create(file.path(outdir, "collected"),showWarnings = F)
  
  withProgress(message = "Collecting samples", {
    i = 0
    for(scanmode in modes){
      collectSamples_2(outdir,
                       scanmode)
      i = i + 0.5
      setProgress(value = i, detail = toupper(scanmode))
    }
  })
  # ---------------
  updateCollapse(session, "pipeline", 
                 style = list("Collate samples" = "success"))
}

run_step_6 <- function(){
  
  dir.create(file.path(outdir, "grouped"),showWarnings = F)
  
  withProgress(message = "Grouping peaks", {
    
    # -------------------
    switch(input$peak_grouping,
           hclust = {
             groupingRestClust_2(outdir = outdir,
                                 cores=cores,
                                 ppm = ppm)
             # -------------
           },
           meanshift = {
             groupingRestBlur(outdir = outdir,
                              cores = input$cores,
                              ppm = input$ppm,
                              mode = mode)
           })
  })
  updateCollapse(session, "pipeline", 
                 style = list("Peak grouping" = "success"))
}

# run_step_7 <- function(){
#   
#   dir.create(paste(outdir, 
#                    "filled", 
#                    sep="/"), 
#              showWarnings = FALSE)
#   
#   print(paste("Filling missing vals in file", basename(f)))
#   # replaceZeros_lookup(resol = resol,
#   #                     outdir = outdir,
#   #                     cores=cores,
#   #                     thresh_pos = thresh_pos,
#   #                     thresh_neg = thresh_neg,
#   #                     scriptDir = scriptdir)
#   replaceZeros_halfmax(resol, outdir, cores, scriptdir)
#   updateCollapse(session, "pipeline", 
#                  style = list("Fill missing values" = "success"))
#   
# }
# # }