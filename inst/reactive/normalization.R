# this triggers when the user wants to normalize their data to proceed to statistics
shiny::observeEvent(input$initialize, {
  
  shiny::withProgress({
    
    shiny::setProgress(session=session, value= .1)
    
    success = F
    try({
      # read in original CSV file
      metshiCSV <- data.table::fread(lcl$paths$csv_loc,
                                     data.table = TRUE,
                                     header = T)
      
      # create empty mSet with 'stat' as default mode
      shiny::setProgress(session=session, value= .2)
      
      mz.meta <- getColDistribution(metshiCSV)
      exp.vars = mz.meta$meta
      mz.vars = mz.meta$mz
      
      # load batch variable chosen by user
      batches <- input$batch_var
      
      # locate qc containing rows in csv
      qc.rows <- which(grepl("QC", metshiCSV$sample))
      
      # for the non-qc samples, check experimental variables. Which have at least 2 different factors, but as little as possible?
      condition <- MetaboShiny::getDefaultCondition(metshiCSV, 
                                                    excl.rows = qc.rows, 
                                                    exp.vars = exp.vars, 
                                                    excl.cond = c("batch",
                                                                  "injection",
                                                                  "sample",
                                                                  "sampling_date"), 
                                                    min.lev = 2)
      
      # =========================================================================
      
      # if nothing is selected for batch, give empty
      if(is.null(batches)) batches <- ""
      
      # only turn on batch correction if user says so
      batch_corr <- if(length(batches) == 1 & 
                       batches[1] == "") FALSE else TRUE
      
      # if 'batch' is selected, 'injection' is often also present
      # TODO: i can imagine this doesn't work for all  users, please disable this...
      if("batch" %in% batches & "injection" %in% colnames(metshiCSV)){
        batches = c(batches, 
                    "injection")
      }
      
      shiny::setProgress(session=session, value= .3)
      
      # get ppm
      params = gsub(lcl$paths$csv_loc, 
                    pattern="\\.csv",
                    replacement="_params.csv")
      
      ppm <- data.table::fread(params)$ppm
      detPPM <- data.table::fread(params)$ppmpermz
      
      # re-make csv with the corrected data
      metshiCSV <- cbind(metshiCSV[,..exp.vars, with=FALSE], # if 'label' is in excel file remove it, it will clash with the metaboanalystR 'label'
                         "label" = metshiCSV[, ..condition][[1]], # set label as the initial variable of interest
                         metshiCSV[,..mz.vars, with=FALSE])
      
      # remove outliers by making a boxplot and going from there
      if(input$remove_outliers){
        metshiCSV <- MetaboShiny::removeOutliers(metshiCSV, 
                                                 exp.vars)
      }
      
      # if QC present, only keep QCs that share batches with the other samples (may occur when subsetting data/only loading part of the samples)
      if(any(grepl("QC", metshiCSV$sample)) & batch_corr){
        metshiCSV <- MetaboShiny::removeUnusedQC(metshiCSV,
                                                 metshiCSV[,..exp.vars, with=FALSE])
      }
      
      shiny::setProgress(session=session, value= .4)
      
      mz.meta <- getColDistribution(metshiCSV)
      exp.vars = mz.meta$meta
      mz.vars = mz.meta$mz
      
      covar_table <- as.data.frame(metshiCSV[,..exp.vars, with=FALSE])
      
      metshiCSV <- MetaboShiny::asMetaboAnalyst(metshiCSV, 
                                                exp.vars)
      
      # define location to write processed csv to
      csv_loc_final <- tempfile()
      
      # write new csv to new location
      data.table::fwrite(metshiCSV, 
                         file = csv_loc_final)
      
      metshiCSV <- NULL

      gc()
      
      shiny::setProgress(session=session, value= .5)
      
      # = = = = = = = = = = = = = = =
      
      mSet <- NULL
      mSet <- MetaboAnalystR::InitDataObjects(data.type = "pktable",
                                              anal.type = "stat",
                                              paired = FALSE)
      
      
      anal.type <<- "stat"
      mSet$dataSet$paired = F
      
      # load new csv into empty mSet!
      mSet <- MetaboAnalystR::Read.TextData(mSet,
                                            filePath = csv_loc_final,
                                            "rowu",
                                            lbl.type = "disc")  # rows contain samples
      
      
      mSet$metshiParams <- list(
        filt_type = input$filt_type,
        miss_type = input$miss_type,
        norm_type = input$norm_type,
        trans_type = input$trans_type,
        scale_type = input$scale_type,
        max.allow = input$maxMz,
        ref_var = input$ref_var,
        batch_var = input$batch_var,
        batch_use_qcs = input$batch_use_qcs,
        prematched = F,
        rf_norm_parallelize = input$rf_norm_parallel,
        rf_norm_ntree = input$rf_norm_ntree,
        miss_perc = input$miss_perc_2,
        orig.count = length(grep("qc",tolower(rownames(mSet$dataSet$orig)),invert = T))
      )
      
      mSet$dataSet$covars <- data.table::as.data.table(covar_table)
      mSet$dataSet$missing <- is.na(mSet$dataSet$orig)
      mSet$dataSet$start <- mSet$dataSet$orig
      
      print(dim(mSet$dataSet$start))
      mSet <- metshiProcess(mSet, session=NULL, init=T)
      
      # save the used adducts to mSet
      mSet$ppm <- ppm
      
      mSet$storage <- list()
      mSet$settings <- list(subset = list(),
                            exp.var = condition,
                            exp.fac = condition,
                            cls.name = condition,
                            time.var = c(),
                            exp.type =  if(mSet$dataSet$cls.num == 2) "1fb" else "1fm",
                            paired = F,
                            filt.type = input$filt_type,
                            orig.count = nrow(mSet$dataSet$norm))

      if(typeof(mSet) != "double"){
        success = T
      }
    })
    
    if(success){
      
      mSet <<- mSet
      
      qs::qsave(mSet, file = file.path(lcl$paths$proj_dir, 
                                       paste0(lcl$proj_name,"_ORIG.metshi")))
      
      fn <- paste0(tools::file_path_sans_ext(lcl$paths$csv_loc), ".metshi")
      save(mSet, file = fn)
      
      uimanager$refresh <- c("general","statspicker")
      plotmanager$make <- "general"
    }else{
      MetaboShiny::metshiAlert("Normalization failed!")
    }
  })
})