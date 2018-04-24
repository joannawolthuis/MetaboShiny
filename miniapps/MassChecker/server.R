shinyServer(function(input, output, session) {
  
  # dir
  shinyDirChoose(input, 
                 'dir', 
                 roots = c(home = '~'), 
                 filetypes = c('', 'txt'))
  dir <- reactive(input$dir)
  output$dir <- renderPrint(dir())
  
  observe({
    shinyFileChoose(input, id = 'sample_sheet_loc', roots = c(home = '~'), filetypes = c('csv','xls', 'xlt', 'xlsx'))
    sample_sheet <- reactive(input$sample_sheet_loc)
    if(!is.null(sample_sheet())){
      sheetpath <<- as.character(parseFilePaths( c(home = '~'), sample_sheet())$datapath)
      sample_tab <<- as.data.table(read.xlsx(file = sheetpath, sheetIndex = 1))[,1:2]
      print(head(sample_tab))
      colnames(sample_tab) <<- c("File.Name", "Card.ID")
    }
    # -----------------------------------------------------------------
  })
  
  observe({
    shinyFileSave(input, id='report_save_loc', roots = c(home = '~'), filetype=c('csv'))
    report_loc <- reactive(input$report_save_loc)
    adj_loc <- report_loc()
    adj_loc$type <- "csv"
    if(!is.null(report_loc()) & exists("final_report") ){ 
      # -------------------
      savepath <- as.character(parseSavePath(roots =  c(home = '~'), selection = adj_loc )$datapath)
      print(savepath)
      write.csv(x = final_report, file = savepath)
    }
    # -----------------------------------------------------------------
  })
  # check
  
  observe({
    shinyFileSave(input, id='ticreport_save_loc', roots = c(home = '~'), filetype=c('csv'))
    tic_report_loc <- reactive(input$ticreport_save_loc)
    adj_loc <- tic_report_loc()
    adj_loc$type <- "csv"
    if(!is.null(tic_report_loc()) & 
       exists("ticSelection") ){ 
      # -------------------
      savepath <- as.character(parseSavePath(roots =  c(home = '~'), selection = adj_loc )$datapath)
      print(savepath)
      write.csv(x = ticSelection, file = savepath)
    }
    # -----------------------------------------------------------------
  })
  
  observeEvent(input$convert, {
    "
    Requires the following:\
    WINE - 32 BITS INSTALL (please adjust your path below to the correct executable)
    - PROPER WINEPREFIX (see below for example)
    MSFileReader - Thermo program, requires registration. Installs through wine.
    - May need to 'install' the DLL's manually through wine.
    READW - Nice converter that requires MSFileReader. Find a one-file executable please!
    --- This may be really tricky to get wo work, I feel lucky it works currently lol :I ---
    "
    minFiles = 0
    progr = 0
    files = list.files(path(),pattern = "raw",full.names = T)
    done = list.files(path(),pattern = "mzXML",full.names = T)
    maxFiles <- length(files) - length(done)
    READW <- file.path(getwd(),"exec/ReAdW.2016010.msfilereader.exe")
    WINE <- "
    "
    # ----------------------
    # files = c("/Users/jwolthuis/Documents/umc/data/Data/BrSpIt/MZXML/RES_20171103_208.raw",
    #           "/Users/jwolthuis/Documents/umc/data/Data/BrSpIt/MZXML/RES_20171103_209.raw",
    #           "/Users/jwolthuis/Documents/umc/data/Data/BrSpIt/MZXML/RES_20171103_210.raw")
    withProgress(min = 0, max = maxFiles, {
      for(INFILE in files){
        OUTFILE = gsub(INFILE, pattern = "\\.raw$|\\.RAW$", replacement = "\\.mzXML")
        if(file.exists(OUTFILE)) next
        system(fn$paste(
          "env WINEPREFIX=/Users/jwolthuis/.wine32 $WINE '$READW' '$INFILE' '$OUTFILE'"
        ))
        progr <- progr + 1
        setProgress(message = paste(basename(INFILE), "converted", sep=": "),value = progr)
      }
    })
  })
  
  
  observeEvent(input$check, {
    minFiles = 0
    progr = 0
    files = list.files(path(),pattern = "mzXML",full.names = T)
    maxFiles <- length(files)
    # ----------------------
    withProgress(min = 0, max = maxFiles, {
      for(f in files){
        result <- CheckFile(f)
        progr <- progr + 1
        setProgress(message = paste(basename(f), result, sep=": "),value = progr)
      }
    })
  })
  
  
  observeEvent(input$plot, {
    req(projname)
    # -----------
    ticDir <<- file.path(getwd(), "www", projname)
    if(!dir.exists(ticDir)) dir.create(ticDir)
    minFiles = 0
    progr = 0
    xmlfiles <<- list.files(path(),pattern = "mzXML",full.names = T)
    donefiles = list.files(ticDir, pattern="png")
    maxFiles <- length(xmlfiles) - length(donefiles)
    if(maxFiles != 0){
      withProgress(min = 0, max = maxFiles, {
        for(f in xmlfiles){
          imgName <- paste(ticDir, "/", gsub(basename(f), pattern = "\\.mzXML", replacement = "\\.png"), sep="")
          if(file.exists(imgName)){print("already plotted");next;}
          # ----------
          rawF <- NULL
          attempt <- 1
          while( is.null(rawF) && attempt <= 3 ) {
            attempt <- attempt + 1
            try(
              rawF <- xcmsRaw(f, profstep=0.1)
            )
          }
          # -----------
          png(filename=imgName, 320, 240)
          plotTIC(rawF, ident=FALSE, msident=FALSE)
          dev.off()
          progr <- progr + 1
          setProgress(message = paste("plotting",basename(f)),value = progr)
        }
      })
    }
    # ----------------------
    tics <<- list.files(ticDir, pattern = "png")
    #files <<- list.files(path(),pattern = "mzXML|raw")
    output$tics <- renderUI({
      req(projname)
      req(ticDir)
      # ---------------------------
      lapply(1:length(tics), FUN=function(x){
        tic <- tics[x]
        sardine(fadeImageButton(paste0("tic_", x), img.path = file.path(projname, tic),value = T))
      })
    })
    ticTable <<- data.table(filenames = gsub(tics, pattern = "\\.png", replacement = ".mzXML"),
                            judgement = c("OK"))
    output$ticTable <- DT::renderDataTable({
      datatable(ticTable,
                autoHideNavigation = T,
                options = list(lengthMenu = c(10, 30, 50), pageLength = 30,scrollX=TRUE, scrollY=TRUE)
      ) %>% formatStyle(
        'judgement',
        target = 'row',
        backgroundColor = styleEqual(c('OK', "BAD"), c('lightgreen', "lightred"))
      )
    })
  })
  
  
  observeEvent(input$submit, {
    nTics <- length(tics)
    judgement <<- sapply(1:nTics, FUN=function(i){
      ticJudge <- !input[[paste0("tic_", i)]]
      # --- return ---
      if(ticJudge == TRUE){"OK"} else("BAD")
    })
    ticSelection <<- data.table(filenames = gsub(tics, pattern = "\\.png", replacement = ".mzXML"),
                                judgement = judgement)
    output$ticTable <- DT::renderDataTable({
      datatable(ticSelection,
                autoHideNavigation = T,
                options = list(lengthMenu = c(10, 30, 50), pageLength = 30,scrollX=TRUE, scrollY=TRUE)
      ) %>% formatStyle(
        'judgement',
        target = 'row',
        backgroundColor = styleEqual(c('OK', 'BAD'), c('lightgreen', 'red'))
      )
    })
  })
  
  
  observeEvent(input$make_report, {
    # create final report.
    # ------------------
    rawfiles <- data.table(RAW=list.files(path(),pattern = "raw"))
    rawfiles <- data.table(filenames=gsub(pattern = "\\.raw", 
                                          replacement = "",
                                          x = list.files("/Users/jwolthuis/PROCESSING_HPC/",
                                                         pattern = "raw")))
    ticSelection_subbed <- ticSelection
    ticSelection_subbed$filenames <- gsub(ticSelection_subbed$filenames,
                                          pattern = "\\.mzXML", 
                                          replacement = "")
    
    setkey(rawfiles,filenames)
    setkey(ticSelection_subbed,filenames)
    
    res_1 <- merge(rawfiles, ticSelection_subbed, all=TRUE)
    res_2 <- merge(res_1, sample_tab, all=TRUE, by.x = "filenames", by.y='File.Name')    
    # -----------
    library(dplyr)
    set.seed(1)
    res_3 <- as.data.table(
      res_2 %>% 
        mutate(OK=(judgement=="OK")) %>% 
        group_by(Card.ID) %>%
        summarise(OK=sum(OK))
    )
    # ===========
    warning_samples <- res_3[OK < 2]
    warning_samples$result <- c("BAD")
    ok_samples <- res_3[OK >= 2]
    ok_samples$result <- c("OK")
    # -----------
    final_report <<- rbind(ok_samples, warning_samples)
    
    output$final_report <- DT::renderDataTable({
      datatable(final_report,
                autoHideNavigation = T,
                options = list(lengthMenu = c(10, 30, 50), pageLength = 30,scrollX=TRUE, scrollY=TRUE)
      ) %>% formatStyle(
        'result',
        target = 'row',
        backgroundColor = styleEqual(c('OK', 'BAD'), c('lightgreen', 'red'))
      )
    })
  })
  
  observeEvent(input$save_report,{
    req(input$report_title)
    # ---------------------
    info <- final_report
    textplot( info, valign="top",cex = 1)
    title(input$report_title)
  })
  
  # path
  path <- reactive({
    home <- normalizePath("~")
    file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
  })
  
  # files
  observe({
    files <<- list.files(path(),pattern = "mzXML|raw")
    projname <<- gsub(basename(files[1]),pattern = ".*\\/|_\\d*.raw|_\\d*.mzXML", replacement = "")
    outdir <<- file.path(path(), "results")
  })
  
  output$files <- renderPrint(list.files(path(),pattern = "mzXML|raw"))
  
  observeEvent(input$create_sample_names, {
    withProgress({
      inj_path <- input$inj_loc$datapath
      
      #inj_path <- "~/Documents/umc/data/Data/Project 2017-026 DSM feed-10 (IMASDE-Madrid) - Saskia v Mil/Injection list IMASDE (Madrid, Spain).xlsx"
      
      mzfiles <- list.files(path(),pattern="\\.mzXML$")
      
      # mzfiles <- list.files("~/Documents/umc/data/Data/Project 2017-026 DSM feed-10 (IMASDE-Madrid) - Saskia v Mil/RES-2017-11-15_DSM DBS Madrid/",
      #                       pattern="\\.mzXML")
      
      injfile <- xlsx::read.xlsx(file = inj_path,
                                 sheetIndex = 1,
                                 colIndex = c(1:3),
                                 rowIndex = 1:length(mzfiles)/3 + 1)
      df <- as.data.frame(injfile)
      i = 1
      df <- data.table::rbindlist(lapply(1:nrow(df), FUN=function(i){
        row <- df[i,]
        if(row$sample == "QC"){
          row$sample <- paste0("QC", i)
          i <<- i + 1
          }
        # --- return ---
        as.data.frame(row)
      }))
      
      filepfx <- gsub(mzfiles[[1]], 
                      pattern = "_\\d*\\.mzXML", 
                      replacement = "" )
      df.expanded <- df[rep(1:nrow(df),each=input$nrepl),] # replicados
      samplecount = nrow(df)
      filenames <- paste(filepfx, c(1:nrow(df.expanded)), sep = "_")
      filenames <- gsub(filenames,
                        pattern="(_)([0-9])$",
                        replacement="\\100\\2")
      filenames <- gsub(filenames,
                        pattern="(_)([0-9][0-9])$",
                        replacement="\\10\\2")
      df.w.filenames <- df.expanded
      df.w.filenames$filename <- filenames
      dt <- as.data.table(df.w.filenames)
      sampleNames <- dt[,c(3,2)]
      colnames(sampleNames) <- c("File_Name", "Sample_Name")
      output$sample_names <- renderDataTable(datatable(sampleNames))
      fwrite(sampleNames,file = file.path(path(), "sampleNames.txt"),sep = "\t")
    })
  })
  
  observeEvent(input$do_step_1, {
    run_step_1()
  })
  observeEvent(input$do_step_2, {
    run_step_2()
  })
  observeEvent(input$do_step_3, {
    run_step_3()
  })
  observeEvent(input$do_step_4, {
    run_step_4()
  })
  observeEvent(input$do_step_5, {
    run_step_5()
  })
  observeEvent(input$do_step_6, {
    run_step_6()
  })
  observeEvent(input$do_step_7, {
    run_step_7()
  })
  
  observeEvent(input$start_pipeline, {
    ##############
    run_step_1()
    run_step_2()
    run_step_3()
    run_step_4()
    run_step_5()
    run_step_6()
    run_step_7()
  })
  
  observe({
    # --------------
    cores <<- input$cores
    nrepl <<- input$nrepl
    trim <<- input$trim
    dimsThresh <<- input$dimsThresh
    resol <<- input$resol
    thresh_pos <<- input$thresh_pos
    thresh_neg <<- input$thresh_neg
    ppm <<- input$ppm
    # ---------------
  })
  
  run_step_1 <- function(){
    # -----------------------
    if(!dir.exists(outdir)) dir.create(outdir)
    # -----------------------
    sn <- data.table::fread(file.path(outdir, "../sampleNames.txt"))
    sn$batch <- as.factor(gsub(sn$File_Name,
                               pattern = "_\\d\\d\\d$",
                               replacement=""))
    data.table::fwrite(sn, file.path(outdir, "sampleNames.txt"))
    # ----------------------------------
    withProgress(message = "Creating signal bins",{
      
      generateBreaksFwhm_jw(mzmin = 70,
                            mzmax = 600,
                            cores = cores,
                            resol=resol,
                            outdir = outdir)
    })
    # ----------------------------------
    updateCollapse(session, "pipeline", 
                   style = list("Initial setup" = "success"))
  }
  
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
    withProgress(message = paste0("Peak finding"),{
      #length(files)
      switch(input$peak_calling,
             gaussian = { # NOOO
               print("deactivated")
               # foreach(i=1:length(files), .options.snow=opts, .export = c("outdir", 
               #                                                            "scriptdir",
               #                                                            "resol",
               #                                                            "thresh_pos", 
               #                                                            "thresh_neg", 
               #                                                            "findPeaks.Gauss.HPC",
               #                                                            "peakFinding.2.0",
               #                                                            "searchMZRange",
               #                                                            "fitGaussianInit",
               #                                                            "generateGaussian",
               #                                                            "fitGaussian",
               #                                                            "getFwhm",
               #                                                            "getSD",
               #                                                            "optimizeGauss",
               #                                                            "getArea",
               #                                                            "fit2G_2",
               #                                                            "fit4G_2",
               #                                                            "fit4peaks",
               #                                                            "fitG_2",
               #                                                            "fit3G_2",
               #                                                            "fit1Peak",
               #                                                            "fit2peaks",
               #                                                            "fit3peaks",
               #                                                            "getFitQuality",
               #                                                            "checkOverlap",
               #                                                            "isWithinXppm",
               #                                                            "sumCurves"),
               #         .verbose = T,
               #         .combine=cfun) %dopar% {
               #           peakFinding.2.0(file = files[i],
               #                           scripts = scriptdir,
               #                           outdir = outdir,
               #                           scanmode = scanmode,
               #                           thresh = if(scanmode == "positive") thresh_pos else{thresh_neg},
               #                           resol = resol)
               #         } 
             }, 
             wavelet = {
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
             },
             centroid = {
               foreach(i=1:length(files), .options.snow=opts, .export = c("outdir",
                                                                          "scriptdir",
                                                                          "resol",
                                                                          "thresh_pos",
                                                                          "thresh_neg",
                                                                          "peakFinding.wavelet_2"),
                       .verbose = T,
                       .packages = c("MALDIquant", "data.table"),
                       .combine=cfun) %dopar% {
                         # --- check if existing ---
                         fn = file.path(outdir, "peaks", basename(files[i]))
                         if(file.exists(fn)) return(NULL)
                         # -------------------------
                         scanmode = if(grepl(files[i],pattern = "_pos.RData")) "pos" else "neg"
                         load(files[i])
                         peaks <- MALDIquant::detectPeaks(averaged, halfWindowSize = 2)
                         peaks_under_thresh <- which(peaks[[1]]@intensity < switch(scanmode, pos = thresh_pos, neg = thresh_neg))
                         peaks[[1]]@snr <- peaks[[1]]@snr[-peaks_under_thresh]
                         peaks[[1]]@mass <- peaks[[1]]@mass[-peaks_under_thresh]
                         peaks[[1]]@intensity <- peaks[[1]]@intensity[-peaks_under_thresh]
                         # --- save ---
                         save(x=peaks, file=fn)
                         # ------------
                       }
             })
    })
    updateCollapse(session, "pipeline", 
                   style = list("Peak finding" = "success"))
    stopCluster(cl)
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
  
  run_step_7 <- function(){
    
    dir.create(paste(outdir, 
                     "filled", 
                     sep="/"), 
               showWarnings = FALSE)
    
    print(paste("Filling missing vals in file", basename(f)))
      # replaceZeros_lookup(resol = resol,
      #                     outdir = outdir,
      #                     cores=cores,
      #                     thresh_pos = thresh_pos,
      #                     thresh_neg = thresh_neg,
      #                     scriptDir = scriptdir)
      replaceZeros_halfmax(resol, outdir, cores, scriptdir)
    updateCollapse(session, "pipeline", 
                   style = list("Fill missing values" = "success"))
    
  }
  
  # run_step_8 <- function(){
  #   # -- FINAL ---
  #   load(file.path(outdir, "repl.pattern.RData"))
  #   i = 0
  #   withProgress(message = "Collecting samples (II)...",{
  #     for(scanmode in modes){
  #       i <<- i + 0.5
  #       setProgress(i, detail=paste("Collecting samples:", scanmode))
  #       # --------------------
  #       f=file.path(outdir, 
  #                   "samplePeaksFilled", 
  #                   if(scanmode == "positive") "positive_rest.RData" else "negative_rest.RData")
  #       load(f)
  #       colnames(outpgrlist) <- c("mzmed", "fq.best", "fq.worst", "npeaks", "mzmin", "mzmax",
  #                                 switch(scanmode, 
  #                                        positive=groupNames.pos, 
  #                                        negative=groupNames.neg), 
  #                                 "avg.int")
  #       fwrite(x = outpgrlist,
  #              file = file.path(outdir, 
  #                               paste0("outlist_",scanmode,".csv")),
  #              sep = "\t")
  #       
  #     }
  #   })
  #   
  #   updateCollapse(session, "pipeline", 
  #                  style = list("Collect samples (II)" = "success"))
  # }
  
  session$onSessionEnded(function() {
    if(any(!is.na(cl))) parallel::stopCluster(cl)
    R.utils::gcDLLs() # flush dlls
    #save(mSet$dataSet, mSet$analSet, file = file.path(options$work_dir, paste0(options$proj_name, ".RData")))
  })
  
# close server obj (must go last!!)
}) 
