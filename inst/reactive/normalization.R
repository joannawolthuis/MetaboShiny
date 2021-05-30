# triggers when probnorm or compnorm is selected
# let user pick a reference condition
ref.selector <- reactive({
  # -------------
  if(input$norm_type == "ProbNorm" | input$norm_type == "CompNorm"){
    shiny::fluidRow(
      shiny::hr(),
      selectInput('ref_var',
                  'What is your reference condition?',
                  choices = c("")),
      actionButton("check_csv",
                   "Get options",
                   icon=shiny::icon("search")),
      shiny::hr()
    )
  }
})

# triggers when check_csv is clicked - get factors usable for normalization
shiny::observeEvent(input$check_csv, {
  req(lcl$paths$csv_loc)
  switch(input$norm_type,
         ProbNorm=shiny::updateSelectInput(session, "ref_var",
                                           choices = get_ref_vars(fac = "label") # please add options for different times later, not difficult
         ),
         CompNorm=shiny::updateSelectInput(session, "ref_var",
                                           choices = get_ref_cpds() # please add options for different times later, not difficult
         ))
})

# render the created UI
output$ref_select <- shiny::renderUI({ref.selector()})


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
      
      if("label" %in% colnames(metshiCSV)) colnames(metshiCSV)[colnames(metshiCSV) == "label"] <- "orig_label"
      
      keep.samps = !duplicated(metshiCSV$sample)
      metshiCSV = metshiCSV[keep.samps,]
      
      # create empty mSet with 'stat' as default mode
      shiny::setProgress(session=session, value= .2)
      
      mz.meta <- getColDistribution(metshiCSV)
      exp.vars = mz.meta$meta
      mz.vars = mz.meta$mz
      
      # load batch variable chosen by user
      batches <- input$batch_var
      
      # locate qc containing rows in csv
      qc.rows <- which(grepl("QC", 
                             metshiCSV$sample))
      
      # for the non-qc samples, check experimental variables. Which have at least 2 different factors, but as little as possible?
      condition <- getDefaultCondition(metshiCSV, 
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
      
      if(exists("mSet")){
        remove(mSet)
      }
      
      mSet <- MetaboAnalystR::InitDataObjects(data.type = "pktable",
                                              anal.type = "stat",
                                              paired = FALSE)
      
      anal.type <<- "stat"
      
      mSet$dataSet$paired <-  mSet$dataSet$ispaired <- mSet$settings$ispaired <- F
      
      # load new csv into empty mSet!
      mSet <- MetaboAnalystR::Read.TextData(mSet,
                                            filePath = csv_loc_final,
                                            "rowu",
                                            lbl.type = "disc")  # rows contain samples
      
      mSet$dataSet$orig <- qs::qread("data_orig.qs")
      
      mSet$metshiParams <- list(
        filt_type = input$filt_type,
        miss_type = input$miss_type,
        norm_type = input$norm_type,
        trans_type = input$trans_type,
        scale_type = input$scale_type,
        max.allow = input$maxMz,
        ref_var = input$ref_mz,
        batch_var = input$batch_var,
        batch_method_a = input$batch_method_a,
        batch_method_b = input$batch_method_b,
        prematched = F,
        rf_norm_parallelize = input$rf_norm_parallel,
        rf_norm_ntree = input$rf_norm_ntree,
        rf_norm_method = if(input$rf_norm_method) "ranger" else "rf",
        miss_perc = input$miss_perc_2,
        orig.count = length(grep("qc",tolower(rownames(mSet$dataSet$orig)),invert = T)),
        pca_corr = input$pca_corr,
        keep_pcs = input$keep_pcs,
        renorm = input$redo_upon_change
      )
      
      mSet$dataSet$covars <- data.table::as.data.table(covar_table)
      mSet$dataSet$missing <- is.na(mSet$dataSet$orig)
      mSet$dataSet$start <- mSet$dataSet$orig
      
      mSet <- metshiProcess(mSet, session=NULL, init=T, cl=session_cl)
      
      # save the used adducts to mSet
      mSet$ppm <- ppm
      
      mSet$storage <- list()
      mSet$settings <- list(subset = list(),
                            exp.var = condition,
                            exp.fac = condition,
                            cls.name = condition,
                            time.var = c(),
                            exp.type =  if(mSet$dataSet$cls.num == 2) "1fb" else "1fm",
                            ispaired = F,
                            filt.type = input$filt_type,
                            orig.count = nrow(mSet$dataSet$norm))

      if(typeof(mSet) != "double"){
        success = T
        qs::qsave(mSet, file = file.path(lcl$paths$proj_dir, 
                                         paste0(lcl$proj_name,"_ORIG.metshi")))
        
        fn <- paste0(tools::file_path_sans_ext(lcl$paths$csv_loc), ".metshi")
        mSet$dataSet$missing <- mSet$dataSet$start <- NULL
        mSet <<- mSet
        filemanager$do <- "save"
        uimanager$refresh <- c("general","statspicker","ml")
        plotmanager$make <- "general"
      }else{
        MetaboShiny::metshiAlert("Normalization failed!")
      }
    })
  })
})