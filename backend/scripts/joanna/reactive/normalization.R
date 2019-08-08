# this triggers when the user wants to normalize their data to proceed to statistics
observeEvent(input$initialize, {

  withProgress({

    shiny::setProgress(session=session, value= .1)

    require(MetaboAnalystR)
    
    # read in original CSV file
    csv_orig <- fread(lcl$paths$csv_loc,
                      data.table = TRUE,
                      header = T)

    # create empty mSet with 'stat' as default mode
    mSet <- InitDataObjects("pktable",
                            "stat",
                            FALSE)
    # set default time series mode'
    mSet$timeseries <- FALSE

    # convert all 0's to NA so metaboanalystR will recognize them
    csv_orig[,(1:ncol(csv_orig)) := lapply(.SD,function(x){ ifelse(x == 0, NA, x)})]

    # remove whitespace
    csv_orig$sample <- gsub(csv_orig$sample, pattern=" |\\(|\\)|\\+", replacement="")

    # find experimental variables by converting to numeric
    as.numi <- as.numeric(colnames(csv_orig)[1:100])
    exp.vars <- which(is.na(as.numi))

    # load batch variable chosen by user
    batches <- input$batch_var

    # locate qc containing rows in csv
    qc.rows <- which(grepl("QC", csv_orig$sample))

    # for the non-qc samples, check experimental variables. Which have at least 2 different factors, but as little as possible?
    unique.levels <- apply(csv_orig[!qc.rows,..exp.vars, with=F], MARGIN=2, function(col){
      lvls <- levels(as.factor(col))
      # - - - - - -
      length(lvls)
    })

    # use this as the default selected experimental variable (user can change later)
    which.default <- unique.levels[which(unique.levels == min(unique.levels[which(unique.levels> 1)]))][1]
    condition = names(which.default)

    # =================================|

    # if nothing is selected for batch, give empty
    if(is.null(batches)) batches <- ""

    # only turn on batch correction if user says so
    batch_corr <- if(length(batches) == 1 & batches[1] == "") FALSE else TRUE

    # if 'batch' is selected, 'injection' is often also present
    # TODO: i can imagine this doesn't work for all  users, please disable this...
    if("batch" %in% batches){
      batches = c(batches, "injection")
    }

    # get the part of csv with only the experimental variables
    first_part <- csv_orig[,..exp.vars, with=FALSE]

    # set NULL or missing levels to "unknown"
    first_part[first_part == "" | is.null(first_part)] <- "unknown"

    # re-make csv with the corrected data
    csv <- cbind(first_part[,-c("label")], # if 'label' is in excel file remove it, it will clash with the metaboanalystR 'label'
                 "label" = first_part[,..condition][[1]], # set label as the initial variable of interest
                 csv_orig[,-..exp.vars,with=FALSE])

    if(all(grepl(pattern = "_T\\d", x = first_part$sample))){
      keep.all.samples <- TRUE
      print("Potential for time series - disallowing outlier removal")
    }else{
      keep.all.samples <- FALSE
    }

    # remove outliers by making a boxplot and going from there
    if(input$remove_outliers & !keep.all.samples){
      sums <- rowSums(csv[,-exp.vars,with=FALSE],na.rm = TRUE)
      names(sums) <- csv$sample
      outliers = c(car::Boxplot(as.data.frame(sums)))
      csv <- csv[!(sample %in% outliers),]
    }

    # remove peaks that are missing in all
    csv <- csv[,which(unlist(lapply(csv, function(x)!all(is.na(x))))),with=F]

    # remove samples with really low numbers of peaks
    complete.perc <- rowMeans(!is.na(csv))
    keep_samps <- csv$sample[which(complete.perc > .2)]
    csv <- csv[sample %in% keep_samps,]

    # also remove them in the table with covariates
    covar_table <- first_part[sample %in% keep_samps,]

    # if the experimental condition is batch, make sure QC samples are not removed at the end for analysis
    # TODO: this is broken with the new system, move this to the variable switching segment of code
    batchview = if(condition == "batch") TRUE else FALSE

    # if QC present, only keep QCs that share batches with the other samples (may occur when subsetting data/only loading part of the samples)
    if(any(grepl("QC", csv$sample))){
      samps <- which(!grepl(csv$sample, pattern = "QC"))
      batchnum <- unique(csv[samps, "batch"][[1]])
      keep_samps_post_qc <- covar_table[which(covar_table$batch %in% batchnum),"sample"][[1]]
      covar_table <- covar_table[which(covar_table$batch %in% batchnum),]
      csv <- csv[which(csv$sample %in% keep_samps_post_qc),-"batch"]
    }

    # rename time column or metaboanalyst won't recognize it
    colnames(csv)[which( colnames(csv) == "time")] <- "Time"

    # deduplicate columns
    # TODO: remove the source of the duplicated columns ealier, might already be fixed
    remove = which( duplicated( t(csv[,..exp.vars] )))
    colnms <- colnames(csv)[remove]
    remove.filt <- setdiff(colnms,c("sample","label"))
    csv <- csv[ , -remove.filt, with = FALSE ]

    # find experimental variables
    as.numi <- as.numeric(colnames(csv)[1:100])
    exp.vars <- which(is.na(as.numi))

    # remove all except sample and time in saved csv
    exp_var_names <- colnames(csv)[exp.vars]
    keep_cols <-  c("sample", "label")
    remove <- which(!(exp_var_names %in% keep_cols))

    # define location to write processed csv to
    csv_loc_final <- gsub(pattern = "\\.csv", replacement = "_no_out.csv", x = lcl$paths$csv_loc)

    # remove file if it already exists
    if(file.exists(csv_loc_final)) file.remove(csv_loc_final)

    # write new csv to new location
    fwrite(csv[,-remove,with=F], file = csv_loc_final)

    # rename row names of covariant table to the sample names
    rownames(covar_table) <- covar_table$sample

    # load new csv into empty mSet!
    mSet <- Read.TextData(mSet,
                          filePath = csv_loc_final,
                          "rowu")  # rows contain samples

    # add covars to the mSet for later switching and machine learning
    mSet$dataSet$covars <- covar_table

    # sanity check data
    mSet <- SanityCheckData(mSet)

    # remove metabolites with more than user defined perc missing
    mSet <- RemoveMissingPercent(mSet,
                                 percent = input$perc_limit/100)

    # missing value imputation
    if(req(input$miss_type ) != "none"){
      if(req(input$miss_type ) == "rowmin"){ # use sample minimum
        new.mat <- apply(mSet$dataSet$preproc, 1, function(x) {
          if (sum(is.na(x)) > 0) {
            x[is.na(x)] <- min(x, na.rm = T)/2
          }
          x
        })
        mSet$dataSet$proc <- t(new.mat)
      }
      else if(req(input$miss_type ) == "pmm"){ # use predictive mean matching
        # TODO: re-enable, it's very slow
        require(mice)
        base <- mSet$dataSet$preproc
        imp <- mice::mice(base, printFlag = TRUE)
      }else if(req(input$miss_type ) == "rf"){ # random forest
        samples <- rownames(mSet$dataSet$preproc)

        # convert all to as numeric
        # TODO: remove, should be automatic
        w.missing <- mSet$dataSet$preproc
        w.missing <- apply(w.missing, 2, as.numeric)

        # register other threads as parallel threads
        doParallel::registerDoParallel(session_cl)

        # set amount of tries (defined by missforest package)
        auto.mtry <- floor(sqrt(ncol(mSet$dataSet$preproc)))

        mtry <- ifelse(auto.mtry > 100, 100, auto.mtry)
        print(paste0("mtry=",mtry))
        
        # impute missing values with random forest
        imp <- missForest::missForest(w.missing,
                                      parallelize = "variables", # parallelize over variables, 'forests' is other option
                                      verbose = F,
                                      ntree = 10,
                                      mtry = mtry)

        mSet$dataSet$proc <- imp$ximp
        rownames(mSet$dataSet$proc) <- rownames(mSet$dataSet$preproc)
        # - - - - - - - - - - - -
      }else{
        # use built in imputation methods, knn means etc.
        mSet <- ImputeVar(mSet,
                          method =  #"knn"
                            input$miss_type
        )
      }
    }

    # if user picked a data filter
    if(req(input$filt_type ) != "none"){
      # filter dataset
      # TODO; add option to only keep columns that are also in QC ('qcfilter'?)
      mSet <- FilterVariable(mSet,
                             filter = input$filt_type,
                             qcFilter = "F",
                             rsd = 25)
    }

    shiny::setProgress(session=session, value= .2)

    # if normalizing by a factor, do the below
    if(req(input$norm_type ) == "SpecNorm"){
      norm.vec <- mSet$dataSet$covars[match(mSet$dataSet$covars$sample,
                                            rownames(mSet$dataSet$preproc)),][[input$samp_var]]
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
      ref_var = input$ref_var
    )

    dput(mSet$metshiParams)

    mSet<-PreparePrenormData(mSet)

    # normalize dataset with user settings(result: mSet$dataSet$norm)
    mSet <- Normalization(mSet,
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

    if(batch_corr){

      if("batch" %in% req(input$batch_var ) & has.qc){

        # get batch for each sample
        batch.idx = as.numeric(as.factor(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$sample),"batch"][[1]]))

        # get injection order for samples
        seq.idx = as.numeric(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$sample),"injection"][[1]])

        # go through all the metabolite columns
        corr_cols <- pbapply::pblapply(1:ncol(mSet$dataSet$norm), function(i){
          # fetch non-corrected values
          vec = mSet$dataSet$norm[,i]
          # correct values using QCs and injectiono rder
          corr_vec = BatchCorrMetabolomics::doBC(Xvec = as.numeric(vec),
                                                 ref.idx = as.numeric(qc_rows),
                                                 batch.idx = batch.idx,
                                                 seq.idx = seq.idx,
                                                 result = "correctedX",
                                                 minBsamp = 1) # at least one QC necessary
          corr_vec
        })

        # cbind the corrected columns to re-make table
        qc_corr_matrix <- as.data.frame(do.call(cbind, corr_cols))
        # fix rownames to old rownames
        colnames(qc_corr_matrix) <- colnames(mSet$dataSet$norm)
        rownames(qc_corr_matrix) <- rownames(mSet$dataSet$norm)
        # save to mSet
        mSet$dataSet$norm <- as.data.frame(qc_corr_matrix)

      }

      # remove QC samples if user doesn't use batch as condition
      if(!batchview & has.qc){
        mSet$dataSet$norm <- mSet$dataSet$norm[-qc_rows,]
        mSet$dataSet$cls <- mSet$dataSet$cls[-qc_rows, drop = TRUE]
        mSet$dataSet$covars <- mSet$dataSet$covars[-grep("QC", mSet$dataSet$covars$sample),]
        mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
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
        # get sample names and classes
        smp <- rownames(mSet$dataSet$norm)
        exp_lbl <- mSet$dataSet$cls

        # create csv for comBat
        csv <- as.data.table(cbind(sample = smp,
                                   label = mSet$dataSet$cls,
                                   mSet$dataSet$norm))

        # transpose for combat
        csv_edata <-t(csv[,!c(1,2)])
        colnames(csv_edata) <- csv$sample

        if(length(left_batch_vars) == 1){
          # create a model table
          csv_pheno <- data.frame(sample = 1:nrow(csv),
                                  batch1 = mSet$dataSet$covars[match(smp, mSet$dataSet$covars$sample),left_batch_vars[1], with=FALSE][[1]],
                                  batch2 = c(0),
                                  outcome = as.factor(exp_lbl))
          # batch correct with comBat
          batch_normalized= t(sva::ComBat(dat=csv_edata,
                                          batch=csv_pheno$batch1
                                          #mod=mod.pheno,
                                          #par.prior=TRUE
          ))
          # fix row names
          rownames(batch_normalized) <- rownames(mSet$dataSet$norm)
        }else{
          # create a model table
          csv_pheno <- data.frame(sample = 1:nrow(csv),
                                  batch1 = mSet$dataSet$covars[match(smp, mSet$dataSet$covars$sample), left_batch_vars[1], with=FALSE][[1]],
                                  batch2 = mSet$dataSet$covars[match(smp, mSet$dataSet$covars$sample), left_batch_vars[2], with=FALSE][[1]],
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
      }
    } else{
      # if qcs presnt and user doesn't want to analyse qc samples
      if(!batchview & has.qc){
        # remove QC rows and associated data from mSet
        mSet$dataSet$norm <- mSet$dataSet$norm[-qc_rows,]
        mSet$dataSet$cls <- mSet$dataSet$cls[-qc_rows, drop = TRUE]
        mSet$dataSet$covars <- mSet$dataSet$covars[-grep("QC", mSet$dataSet$covars$sample),]
        mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
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

    varNormPlots <- ggplotNormSummary(mSet = mSet,
                                      colmap = lcl$aes$mycols,
                                      plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                      font = lcl$aes$font)

    output$var1 <- renderPlot(varNormPlots$tl)
    output$var2 <- renderPlot(varNormPlots$bl)
    output$var3 <- renderPlot(varNormPlots$tr)
    output$var4 <- renderPlot(varNormPlots$br)

    sampNormPlots <-  ggplotSampleNormSummary(mSet,
                                              #colmap = lcl$aes$mycols,
                                              plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                              font = lcl$aes$font)
    output$samp1 <- renderPlot(sampNormPlots$tl)
    output$samp2 <- renderPlot(sampNormPlots$bl)
    output$samp3 <- renderPlot(sampNormPlots$tr)
    output$samp4 <- renderPlot(sampNormPlots$br)
    shiny::setProgress(session=session, value= .8)

    # save the used adducts to mSet
    shiny::setProgress(session=session, value= .9)

    mSet$storage <- list()

    # return?
    mSet <<- mSet

    statsmanager$calculate <- "pca"
    statsmanager$reload <- "pca"
    datamanager$reload <- "general"

  })
})
