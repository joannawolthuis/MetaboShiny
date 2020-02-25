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
      
      colnames(metshiCSV) <- gsub(colnames(metshiCSV), pattern="\\|", replacement="/")
      
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
                                                    excl.cond = c("batch", "injection", "sample"), 
                                                    min.lev = 2)
      
      # =================================|
      

      # if nothing is selected for batch, give empty
      if(is.null(batches)) batches <- ""
      
      # only turn on batch correction if user says so
      batch_corr <- if(length(batches) == 1 & batches[1] == "") FALSE else TRUE
      
      # if 'batch' is selected, 'injection' is often also present
      # TODO: i can imagine this doesn't work for all  users, please disable this...
      if("batch" %in% batches & "injection" %in% colnames(metshiCSV)){
        batches = c(batches, "injection")
      }
      
      shiny::setProgress(session=session, value= .3)
      
      # get ppm
      params = gsub(lcl$paths$csv_loc, pattern="\\.csv", replacement="_params.csv")
      ppm <- fread(params)$ppm
      detPPM <- fread(params)$ppmpermz

      # re-make csv with the corrected data
      metshiCSV <- cbind(metshiCSV[,..exp.vars, with=FALSE][, -c("label")], # if 'label' is in excel file remove it, it will clash with the metaboanalystR 'label'
                         "label" = metshiCSV[, ..condition][[1]], # set label as the initial variable of interest
                         metshiCSV[,..mz.vars, with=FALSE])
      
      # remove outliers by making a boxplot and going from there
      if(input$remove_outliers){
        metshiCSV <- MetaboShiny::removeOutliers(metshiCSV, 
                                                 exp.vars)
      }
      
      # if QC present, only keep QCs that share batches with the other samples (may occur when subsetting data/only loading part of the samples)
      if(any(grepl("QC", metshiCSV$sample))){
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
      
      mSet <- MetaboAnalystR::InitDataObjects("pktable",
                                              "stat",
                                              FALSE)
      
      # load new csv into empty mSet!
      mSet <- MetaboAnalystR::Read.TextData(mSet,
                                            filePath = csv_loc_final,
                                            "rowu")  # rows contain samples
      
      # set default time series mode'
      mSet$dataSet$paired = F
      mSet$dataSet$subset = c()
      
      # add covars to the mSet for later switching and machine learning
      mSet$dataSet$covars <- data.table::as.data.table(covar_table)
      #rownames(mSet$dataSet$covars) <- mSet$dataSet$covars$sample
      
      # sanity check data
      mSet <- MetaboAnalystR::SanityCheckData(mSet)
      
      mSet$dataSet$orig <- NULL
      gc()
      
      #print(paste0("Left after missing value check: ", MetaboShiny::mzLeftPostFilt(mSet, input$perc_limit))
      
      shiny::setProgress(session=session, value= .6)
      
      # remove metabolites with more than user defined perc missing
      #mSet <- MetaboAnalystR::RemoveMissingPercent(mSet,
      #                                             percent = 1)#input$perc_limit/100)
      
      # remove samples with now no one...
      #rmv = MetaboShiny::tooEmptySamps(mSet, 
      #                                 max.missing.per.samp = input$perc_limit)
      
      # missing value imputation
      if(req(input$miss_type ) != "none"){
        if(req(input$miss_type ) == "rowmin"){ # use sample minimum
          mSet$dataSet$proc <- replRowMin(mSet)
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
      
      if(req(input$filt_type ) != "none"){
        shiny::showNotification("Filtering dataset...")
        # TODO; add option to only keep columns that are also in QC ('qcfilter'?)
        mSet <- MetaboAnalystR::FilterVariable(mSet,
                                               filter = input$filt_type,
                                               qcFilter = "F",
                                               rsd = 25)
        # keep.mz <- colnames(mSet$dataSet$filt)
        # mSet <- MetaboShiny::filt.mSet(mSet, keep.mz)
      } 
      
      mSet$dataSet$preproc <- NULL
      gc()
      
      shiny::setProgress(session=session, value= .7)
      
      # if normalizing by a factor, do the below
      if(req(input$norm_type) == "SpecNorm"){
        norm.vec <<- mSet$dataSet$covars[match(mSet$dataSet$covars$sample,
                                               rownames(mSet$dataSet$preproc)
                                               ),][[input$samp_var]]
        norm.vec <<- scale(x = norm.vec, center = 1)[,1] # normalize scaling factor
      }else{
        norm.vec <<- rep(1, length(mSet$dataSet$cls)) # empty
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
      
     # mSet$dataSet$proc <- NULL
     # gc()
      
      # normalize dataset with user settings(result: mSet$dataSet$norm)
      mSet <- MetaboAnalystR::Normalization(mSet,
                                            rowNorm = input$norm_type,
                                            transNorm = input$trans_type,
                                            scaleNorm = input$scale_type,
                                            ref = input$ref_var)
      
      
      mSet$dataSet$proc <- NULL
      #mSet$dataSet$prenorm <- NULL
      gc()
      
      shiny::setProgress(session=session, value= .8)
      
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
          smps <- rownames(mSet$dataSet$norm)
          # get which rows are QC samples
          qc_rows <- which(grepl(pattern = "QC", x = smps))
          # get batch for each sample
          batch.idx = as.numeric(as.factor(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$sample),"batch"][[1]]))
          if(length(batch.idx) == 0) return(mSet$dataSet$norm)
          # get injection order for samples
          seq.idx = as.numeric(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$sample),"injection"][[1]])
          # go through all the metabolite columns
          dtNorm <- as.data.table(mSet$dataSet$norm)
          pb <- pbapply::startpb(0, max = ncol(dtNorm))
          i = 0
          dtNorm[,(1:ncol(dtNorm)) := lapply(.SD,function(x){ 
            i <<- i + 1
            pbapply::setpb(pb, i)
            BatchCorrMetabolomics::doBC(Xvec = as.numeric(x),
                                        ref.idx = as.numeric(qc_rows),
                                        batch.idx = batch.idx,
                                        seq.idx = seq.idx,
                                        result = "correctedX",
                                        minBsamp = 1)
            
            }), .SDcols = 1:ncol(dtNorm)]
        }
        
        # remove QC samples if user doesn't use batch as condition
        #if(!batchview & has.qc){
        #}
        
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
            csv_pheno <- data.frame(sample = 1:nrow(mSet$dataSet$covars),
                                    batch1 = mSet$dataSet$covars[match(smps,mSet$dataSet$covars$sample),
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
            csv_pheno <- data.frame(sample = 1:nrow(mSet$dataSet$covars),
                                    batch1 = mSet$dataSet$covars[match(smps,mSet$dataSet$covars$sample), left_batch_vars[1], with=FALSE][[1]],
                                    batch2 = mSet$dataSet$covars[match(smps,mSet$dataSet$covars$sample), left_batch_vars[2], with=FALSE][[1]],
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
      }}
      # else{
      #     # if qcs presnt and user doesn't want to analyse qc samples
      #     if(!batchview & has.qc){
      #       # remove QC rows and associated data from mSet
      #       mSet <- MetaboShiny::hideQC(mSet)
      #     }
      #   }
      
      shiny::setProgress(session=session, value= .9)
      
      mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
      

      # make sure covars order is consistent with mset$..$norm order
      rematch = match(rownames(mSet$dataSet$norm),
                      mSet$dataSet$covars$sample)
      mSet$dataSet$covars <- mSet$dataSet$covars[rematch,]
      
      # set name of variable of interest
      mSet$dataSet$cls.name <- condition
      mSet$dataSet$exp.fac <- condition
      mzs = colnames(mSet$dataSet$norm)
      
      # if(detPPM == "yes"){
      #   colnames(mSet$dataSet$norm) <- gsub(colnames(mSet$dataSet$norm), 
      #                                       pattern = "/.*$", 
      #                                       replacement = "")
      #   colnames(mSet$dataSet$prenorm) <- gsub(colnames(mSet$dataSet$prenorm), 
      #                                          pattern = "/.*$", 
      #                                          replacement = "")
      #   splitMZ = stringr::str_split(mzs, pattern = "/")
      #   mSet$unique_ppm <- cbind(mzs, data.table::rbindlist(lapply(stringr::str_split(mzs, pattern = "/"), as.list)))
      #   colnames(mSet$unique_ppm) <- c("mzname", "mz", "ppm")
      # }else{
      #   mSet$unique_ppm <- data.table::data.table(mzname = colnames(mSet$dataSet$norm),
      #                                             mz = colnames(mSet$dataSet$norm),
      #                                             ppm = c(ppm))
      # }
  
      # save the used adducts to mSet
      mSet$ppm <- ppm
      mSet$paired <- F
      mSet$dataSet$subset <- list()
      mSet$dataSet$exp.var <- condition
      mSet$dataSet$time.var <- c()
      mSet$storage <- list()
      
      if(mSet$dataSet$cls.num == 2){
        mSet$dataSet$exp.type <- "1fb"
      }else{
        mSet$dataSet$exp.type <- "1fm"
      }
      
      mSet$settings <- list(subset = mSet$dataSet$subset,
                            exp.var = mSet$dataSet$exp.var,
                            exp.fac = mSet$dataSet$exp.fac,
                            time.var = mSet$dataSet$time.var,
                            exp.type = mSet$dataSet$exp.type,
                            paired = mSet$paired,
                            orig.count = nrow(mSet$dataSet$prenorm))
      
      if(typeof(mSet) != "double"){
        success = T
      }
    })
    
    if(success){
      saveRDS(list(data = mSet$dataSet,
                   analysis = mSet$analSet,
                   settings = mSet$settings), 
              file = file.path(lcl$paths$proj_dir, 
                               paste0(lcl$proj_name,"_ORIG.metshi")))
     
      # mSet <- readRDS(file.path(lcl$paths$proj_dir, 
      #                           paste0(lcl$proj_name,"_ORIG.metshi")))
      #names(mSet) <- c("dataSet", "analSet", "settings")
      if(has.qc){
        mSet <- MetaboShiny::hideQC(mSet)
      }
      mSet <<- mSet
      
      fn <- paste0(tools::file_path_sans_ext(lcl$paths$csv_loc), ".metshi")
      save(mSet, file = fn)
      
      datamanager$reload <- c("general","statspicker")
      statsmanager$calculate <- "pca"
      statsmanager$reload <- "pca"  
    }else{
      MetaboShiny::metshiAlert("Normalization failed!")
    }
  })
})