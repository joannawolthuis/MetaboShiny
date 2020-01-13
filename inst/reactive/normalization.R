# this triggers when the user wants to normalize their data to proceed to statistics
shiny::observeEvent(input$initialize, {
  
  shiny::withProgress({
    
    shiny::setProgress(session=session, value= .1)
    
    # read in original CSV file
    csv_orig <- data.table::fread(lcl$paths$csv_loc,
                                  data.table = TRUE,
                                  header = T)
    
    # create empty mSet with 'stat' as default mode
    mSet <- MetaboAnalystR::InitDataObjects("pktable",
                            "stat",
                            FALSE)
    
    # set default time series mode'
    mSet$dataSet$paired = F
    mSet$dataSet$subset = c()
    
    mz.meta <- MetaboShiny::getColDistribution(csv_orig)
    exp.vars = mz.meta$meta
    mz.vars = mz.meta$mz
    
    csv_orig <- MetaboShiny::cleanCSV(csv_orig, 
                                      exp.vars = exp.vars, mz.vars = mz.vars, 
                                      miss.meta = "unknown", miss.mz = NA)
    
    # load batch variable chosen by user
    batches <- input$batch_var
    
    # locate qc containing rows in csv
    qc.rows <- which(grepl("QC", csv_orig$sample))
    
    # for the non-qc samples, check experimental variables. Which have at least 2 different factors, but as little as possible?
    condition <- MetaboShiny::getDefaultCondition(csv_orig, 
                                                  excl.rows = qc.rows, 
                                                  exp.vars = exp.vars, 
                                                  excl.cond = c("batch", "injection", "sample"), 
                                                  min.lev = 2)

    # =================================|
    
    # if nothing is selected for batch, give empty
    if(is.null(batches)) batches <- ""
    
    # only turn on batch correction if user says so
    batch_corr <- if(length(batches) == 1 & batches[1] == "") FALSE else TRUE
    
    # if 'batch' is selected, 'injection' is often also present
    # TODO: i can imagine this doesn't work for all  users, please disable this...
    if("batch" %in% batches & "injection" %in% colnames(csv_orig)){
      batches = c(batches, "injection")
    }
    
    mz.meta <- MetaboShiny::getColDistribution(csv_orig)
    exp.vars = mz.meta$meta
    mz.vars = mz.meta$mz
    
    # get the part of csv with only the experimental variables
    first_part <- csv_orig[,..exp.vars, with=FALSE]
    # re-make csv with the corrected data
    csv <- cbind(first_part, # if 'label' is in excel file remove it, it will clash with the metaboanalystR 'label'
                 "label" = first_part[,..condition][[1]], # set label as the initial variable of interest
                 csv_orig[,-..exp.vars,with=FALSE])
    
    # remove outliers by making a boxplot and going from there
    if(input$remove_outliers){
      csv <- MetaboShiny::removeOutliers(csv, exp.vars)
    }
    
    # also remove them in the table with covariates
    covar_table <- first_part #[sample %in% keep_samps,]
    
    # if the experimental condition is batch, make sure QC samples are not removed at the end for analysis
    # TODO: this is broken with the new system, move this to the variable switching segment of code
    batchview = if(condition == "batch") TRUE else FALSE
    
    # if QC present, only keep QCs that share batches with the other samples (may occur when subsetting data/only loading part of the samples)
    if(any(grepl("QC", csv$sample))){
      csv <- MetaboShiny::removeUnusedQC(csv, covar_table)
    }
    
    mz.meta <- MetaboShiny::getColDistribution(csv)
    exp.vars = mz.meta$meta
    mz.vars = mz.meta$mz
    
    csv <- MetaboShiny::asMetaboAnalyst(csv, exp.vars)
    
    # define location to write processed csv to
    csv_loc_final <- gsub(pattern = "\\.csv", replacement = "_no_out.csv", x = lcl$paths$csv_loc)
    # remove file if it already exists
    if(file.exists(csv_loc_final)) file.remove(csv_loc_final)
    # write new csv to new location
    data.table::fwrite(csv, file = csv_loc_final)
    
    # rename row names of covariant table to the sample names
    rownames(covar_table) <- covar_table$sample
    
    # load new csv into empty mSet!
    mSet <- MetaboAnalystR::Read.TextData(mSet,
                          filePath = csv_loc_final,
                          "rowu")  # rows contain samples
    
    # add covars to the mSet for later switching and machine learning
    mSet$dataSet$covars <- covar_table
    
    # sanity check data
    mSet <- MetaboAnalystR::SanityCheckData(mSet)
    
    MetaboShiny::mzLeftPostFilt(mSet, input$perc_limit)
    
    # remove metabolites with more than user defined perc missing
    mSet <- MetaboAnalystR::RemoveMissingPercent(mSet,
                                 percent = input$perc_limit/100)
    
    # remove samples with now no one...
    rmv = MetaboShiny::tooEmptySamps(mSet, max.missing.per.samp = input$perc_limit)
    
    # missing value imputation
    if(req(input$miss_type ) != "none"){
      if(req(input$miss_type ) == "rowmin"){ # use sample minimum
        mSet$dataSet$proc <- MetaboShiny::replRowMin(mSet)
      }
      else if(req(input$miss_type ) == "pmm"){ # use predictive mean matching
        # TODO: re-enable, it's very slow
        base <- mSet$dataSet$preproc
        imp <- mice::mice(base, printFlag = TRUE)
        
      }else if(req(input$miss_type ) == "rf"){ # random forest
        mSet$dataSet$proc <- MetaboShiny::replRF(mSet, 
                                                 parallelMode = input$rf_norm_parallelize, 
                                                 ntree = input$rf_norm_ntree,
                                                 cl = session_cl)
        rownames(mSet$dataSet$proc) <- rownames(mSet$dataSet$preproc)
        # - - - - - - - - - - - -
      }else{
        # use built in imputation methods, knn means etc.
        mSet <- MetaboAnalystR::ImputeVar(mSet,
                                          method = input$miss_type
        )
      }
    }
    
    shiny::setProgress(session=session, value= .2)
    
    # if normalizing by a factor, do the below
    if(req(input$norm_type) == "SpecNorm"){
      norm.vec <- mSet$dataSet$covars[match(rownames(mSet$dataSet$preproc),
                                            mSet$dataSet$covars$sample),][[input$samp_var]]
      norm.vec <- scale(x = norm.vec, center = 1) # normalize scaling factor
      
    }else{
      norm.vec <- rep(1, length(mSet$dataSet$cls)) # empty
    }
    
    # write these to mset!!!
    mSet$metshiParams <- list(
      perc_limit = input$perc_limit,
      filt_type = input$filt_type,
      miss_type = input$miss_type,
      norm_type = input$norm_type,
      trans_type = input$trans_type,
      scale_type = input$scale_type,
      ref_var = input$ref_var,
      prematched = F
    )
    
    mSet <- MetaboAnalystR::PreparePrenormData(mSet)
    
    # normalize dataset with user settings(result: mSet$dataSet$norm)
    mSet <- MetaboAnalystR::Normalization(mSet,
                                          rowNorm = input$norm_type,
                                          transNorm = input$trans_type,
                                          scaleNorm = input$scale_type,
                                          ref = input$ref_var)
    
    shiny::setProgress(session=session, value= .4)
    
    # get sample names
    smps <- rownames(mSet$dataSet$norm)
    
    # get which rows are QC samples
    qc_rows <- which(grepl(pattern = "QC", x = smps))
    
    # if at least one row has a QC in it, batch correct
    has.qc <- length(qc_rows) > 0
    
    # lowercase all the covars table column names
    colnames(mSet$dataSet$covars) <- tolower(colnames(mSet$dataSet$covars))
    
    # === check if it does wrong here... ===
    
    if(batch_corr){
      
      if("batch" %in% req(input$batch_var ) & has.qc){
        # save to mSet
        mSet$dataSet$norm <- MetaboShiny::batchCorrQC(mSet, qc_rows)
      }
      
      # remove QC samples if user doesn't use batch as condition
      if(!batchview & has.qc){
        mSet <- MetaboShiny::hideQC(mSet)
      }
      
      # check which batch values are left after initial correction
      left_batch_vars <- grep(input$batch_var,
                              pattern =  ifelse(has.qc, "batch|injection|sample", "injection|sample"),
                              value = T,
                              invert = T)
      
      if(length(left_batch_vars) > 2){
        NULL  # no option for more than 2 other batch variables yet
      } else if(length(left_batch_vars) == 0){
        NULL # if none left, continue after this
      } else{
        
        csv_edata <- MetaboShiny::combatCSV(mSet)
        
        if(length(left_batch_vars) == 1){
          # create a model table
          csv_pheno <- data.frame(sample = 1:nrow(csv),
                                  batch1 = mSet$dataSet$covars[match(smp, 
                                                                     mSet$dataSet$covars$sample),
                                                               left_batch_vars[1], with=FALSE][[1]]
                                  #,batch2 = c(0),
                                  #outcome = as.factor(mSet$dataSet$cls)
                                  )
          # batch correct with comBat
          batch_normalized = t(sva::ComBat(dat = csv_edata,
                                           batch = csv_pheno$batch1
                                           # mod=mod.pheno,
                                           # par.prior=TRUE
          ))
          # fix row names
          rownames(batch_normalized) <- rownames(mSet$dataSet$norm)
        }else{
          # create a model table
          csv_pheno <- data.frame(sample = 1:nrow(csv),
                                  batch1 = mSet$dataSet$covars[smp, match(mSet$dataSet$covars$sample), left_batch_vars[1], with=FALSE][[1]],
                                  batch2 = mSet$dataSet$covars[smp, match(mSet$dataSet$covars$sample), left_batch_vars[2], with=FALSE][[1]],
                                  outcome = as.factor(exp_lbl))
          # batch correct with limma and two batches
          batch_normalized = t(limma::removeBatchEffect(x = csv_edata,
                                                        #design = mod.pheno,
                                                        batch = csv_pheno$batch1,
                                                        batch2 = csv_pheno$batch2))
          rownames(batch_normalized) <- rownames(mSet$dataSet$norm)
        }
        # save normalized table to mSet
        mSet$dataSet$norm <- as.data.frame(batch_normalized)
      }}else{
      # if qcs presnt and user doesn't want to analyse qc samples
      if(!batchview & has.qc){
        # remove QC rows and associated data from mSet
        mSet <- MetaboShiny::hideQC(mSet)
      }
    }
    
    mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
    
    shiny::setProgress(session=session, value= .5)
    
    # make sure covars order is consistent with mset$..$norm order
    mSet$dataSet$covars <- mSet$dataSet$covars[match(rownames(mSet$dataSet$norm),
                                                     mSet$dataSet$covars$sample),]
    
    # set name of variable of interest
    mSet$dataSet$cls.name <- condition
    
    shiny::setProgress(session=session, value= .6)
    
    # generate summary plots and render them in UI
    
    varNormPlots <- MetaboShiny::ggplotNormSummary(mSet = mSet,
                                      plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                      font = lcl$aes$font,
                                      cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
    
    output$var1 <- shiny::renderPlot(varNormPlots$tl)
    output$var2 <- shiny::renderPlot(varNormPlots$bl)
    output$var3 <- shiny::renderPlot(varNormPlots$tr)
    output$var4 <- shiny::renderPlot(varNormPlots$br)
    
    sampNormPlots <- MetaboShiny::ggplotSampleNormSummary(mSet,
                                              plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                              font = lcl$aes$font,
                                              cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
    output$samp1 <- shiny::renderPlot(sampNormPlots$tl)
    output$samp2 <- shiny::renderPlot(sampNormPlots$bl)
    output$samp3 <- shiny::renderPlot(sampNormPlots$tr)
    output$samp4 <- shiny::renderPlot(sampNormPlots$br)
    shiny::setProgress(session=session, value= .8)
    
    # save the used adducts to mSet
    shiny::setProgress(session=session, value= .9)
    
    # get ppm
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), lcl$paths$patdb)
    
    mSet$ppm <- sprintf("%.1f",RSQLite::dbGetQuery(conn, "select ppm from params"))
    
    RSQLite::dbDisconnect(conn)
    
    mSet$paired <- F
    mSet$dataSet$subset <- list()
    mSet$dataSet$exp.var <- condition
    mSet$dataSet$time.var <- c()
    if(mSet$dataSet$cls.num == 2){
      mSet$dataSet$exp.type <- "1fb"
    }else{
      mSet$dataSet$exp.type <- "1fm"
    }
    
    mSet$storage <- list(orig = list(data = mSet$dataSet,
                                     analysis = mSet$analSet))

    
    if(req(input$filt_type ) != "none"){
      shiny::showNotification("Filtering dataset...")
      # TODO; add option to only keep columns that are also in QC ('qcfilter'?)
      mSet <- MetaboAnalystR::FilterVariable(mSet,
                                             filter = input$filt_type,
                                             qcFilter = "F",
                                             rsd = 25)
      keep.mz <- colnames(mSet$dataSet$filt)
      mSet <- MetaboShiny::filt.mSet(mSet, keep.mz)
    }
    
    mSet <<- mSet
    
    datamanager$reload <- c("general","statspicker")
    statsmanager$calculate <- "pca"
    statsmanager$reload <- "pca"
    
  })
})