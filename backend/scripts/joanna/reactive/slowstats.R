# triggers when the 'go' button is pressed on the PLS-DA tab
observeEvent(input$do_plsda, {

  require(e1071)

  # depending on type, do something else
  # TODO: enable sparse and orthogonal PLS-DA
  switch(input$plsda_type,
         normal={
           require(caret)
           withProgress({
             mSet <<- PLSR.Anal(mSet) # perform pls regression
             setProgress(0.3)
             mSet <<- PLSDA.CV(mSet, methodName=if(nrow(mSet$dataSet$norm) < 50) "L" else "T",compNum = 3) # cross validate
             setProgress(0.6)
             mSet <<- PLSDA.Permut(mSet,num = 300, type = "accu") # permute
           })
         },
         sparse ={
           mSet <<- SPLSR.Anal(mSet, comp.num = 3)
         })
  # reload pls-da plots
  datamanager$reload <- "plsda"
})


# triggers if 'go' is pressed in the machine learning tab
observeEvent(input$do_ml, {
  withProgress({

    setProgress(value = 0)

    # get base table to use for process
    curr <- as.data.table(mSet$dataSet$preproc) # the filtered BUT NOT IMPUTED table, ML should be able to deal w/ missing values

    # replace NA's with zero
    curr[,(1:ncol(curr)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]

    # conv to data frame
    curr <- as.data.frame(curr)
    rownames(curr) <- rownames(mSet$dataSet$preproc)

    if(!is.null(input$timecourse_trigger)){
      if(input$timecourse_trigger){
        melted_curr <- melt(as.data.table(curr, keep.rownames = TRUE),id.vars = "rn")
        split_rn = strsplit(melted_curr$rn, split = "_T")
        melted_curr$time <- as.numeric(sapply(split_rn, function(x) x[[2]]))
        melted_curr$rn <- sapply(split_rn, function(x) x[[1]])
        melted_curr$variable <- paste0(melted_curr$variable, "_T", melted_curr$time)
        curr <- reshape::cast(melted_curr[,-"time"])
        curr <- as.data.frame(curr)
        rownames(curr) <- curr$rn
        curr <- curr[,-1]
      }
    }

    # find the qc rows
    is.qc <- grepl("QC|qc", rownames(curr))
    curr <- curr[!is.qc,]

    # reorder according to covars table (will be used soon)
    if(!is.null(input$timecourse_trigger)){
      if(input$timecourse_trigger){
        config <- unique(mSet$dataSet$covars[, c("sample",
                                                 "sex")]) # reorder so both halves match up later
        config$sample <- gsub(x = config$sample, pattern = "_T\\d", replacement = "")
        config <- unique(config)
        order <- match(config$sample,rownames(curr))
        config <- cbind(config[order,], label=mSet$dataSet$cls[order]) # add current experimental condition
        config <- config[,apply(!is.na(config), 2, any), with=FALSE]
      }else{
        order <- match(mSet$dataSet$covars$sample,rownames(curr))
        config <- mSet$dataSet$covars[order, -"label"] # reorder so both halves match up later
        config <- cbind(config, label=mSet$dataSet$cls[order]) # add current experimental condition
        config <- config[,apply(!is.na(config), 2, any), with=FALSE]
      }
    }else{
      order <- match(mSet$dataSet$covars$sample,rownames(curr))
      config <- mSet$dataSet$covars[order, -"label"] # reorder so both halves match up later
      config <- cbind(config, label=mSet$dataSet$cls[order]) # add current experimental condition
      config <- config[,apply(!is.na(config), 2, any), with=FALSE]
    }

    # remove ones w/ every row being different(may be used to identify...)
    covariates <- lapply(1:ncol(config), function(i) as.factor(config[,..i][[1]]))
    names(covariates) <- colnames(config)

    # remove ones with na present
    has.na <- sapply(covariates, function(x) any(is.na(x)))
    has.all.unique <- sapply(covariates, function(x) length(unique(x)) == length(x))

    # rename the variable of interest to 0-1-2 etc.
    char.lbl <- as.character(covariates$label)
    uniques <- unique(char.lbl)
    uniques_new_name <- c(1:length(uniques))
    names(uniques_new_name) = uniques

    remapped.lbl <- uniques_new_name[char.lbl]

    # find which variables are covariant with label, they will be removed
    covariant.with.label <- sapply(covariates, function(x){
      char.x <- as.character(x)
      uniques <- unique(char.x)
      uniques_new_name <- c(1:length(uniques))
      names(uniques_new_name) = uniques
      remapped.x = uniques_new_name[char.x]
      res = if(length(remapped.x) == length(remapped.lbl)){
        all(remapped.x == remapped.lbl)
      }else{ FALSE }
      res
    })

    # now filter out unique, covariate and with missing columns from $covars
    keep_configs <- which(!(names(config) %in% names(config)[unique(c(which(has.na), which(has.all.unique), which(covariant.with.label)))]))
    keep_configs <- c(keep_configs, which(names(config) == "label"))
    #keep_configs <- which(names(config) == "label")

    print("Removing covariates and unique columns. Keeping non-mz variables:")
    print(names(config)[keep_configs])

    config <- config[,..keep_configs,with=F]

    # - - - - - - - - - - - - - - - - - - - - - - -

    # join halves together, user variables and metabolite data
    curr <- cbind(config, curr)
    curr <- as.data.table(curr)

    # remove cols with all NA
    curr <- curr[,colSums(is.na(curr)) < nrow(curr), with=FALSE]
    # remove rows with all NA
    curr <- curr[complete.cases(curr),]

    # how many models will be built? user input
    goes = as.numeric(input$ml_attempts)

    if(is.null(lcl$vectors$ml_train)){
      lcl$vectors$ml_train <<- c("all", "all")
    }
    if(is.null(lcl$vectors$ml_test)){
      lcl$vectors$ml_test <<- c("all", "all")
    }

    if(all(lcl$vectors$ml_test == lcl$vectors$ml_train)){
      if(unique(lcl$vectors$ml_test) == "all"){
        print("nothing selected... continuing in normal non-subset mode")
      }else{
        print("no can do, need to test on something else than train!!!")
        return(NULL)
      }
    }

    # identify which columns are metabolites and which are config/covars
    configCols <- which(!(gsub(x = colnames(curr), pattern = "_T\\d", replacement="") %in% colnames(mSet$dataSet$norm)))
    mzCols <- which(gsub(x = colnames(curr), pattern = "_T\\d", replacement="") %in% colnames(mSet$dataSet$norm))

    # make the covars factors and the metabolites numeric.
    curr[,(configCols):= lapply(.SD, function(x) as.factor(x)), .SDcols = configCols]
    curr[,(mzCols):= lapply(.SD, function(x) as.numeric(x)), .SDcols = mzCols]

    # ============ LOOP HERE ============

    # get results for the amount of attempts chosen
    repeats <- pbapply::pblapply(1:goes,
                                 cl=0,
                                 #cl=session_cl,
                                 function(i,
                                          train_vec = train_vec,
                                          test_vec = test_vec,
                                          configCols = configCols){
      shiny::isolate({

        # get user training percentage
        ml_train_perc <- input$ml_train_perc/100

        if(unique(train_vec)[1] == "all" & unique(test_vec)[1] == "all"){ # BOTH ARE NOT DEFINED
          test_idx = caret::createDataPartition(y = curr$label, p = ml_train_perc, list = FALSE) # partition data in a balanced way (uses labels)
          train_idx = setdiff(1:nrow(curr), test_idx) #use the other rows for testing
          inTrain = train_idx
          inTest = test_idx
        }else if(unique(train_vec)[1] != "all"){ #ONLY TRAIN IS DEFINED
          train_idx <- which(config[,train_vec[1], with=F][[1]] == train_vec[2])
          #train_idx = grep(config$sample, pattern = ml_train_regex) # get training sample ids with regex
          test_idx = setdiff(1:nrow(curr), train_idx) # use the other rows for testing
          reTrain <- caret::createDataPartition(y = config[train_idx, label], p = ml_train_perc) # take a user-defined percentage of the regexed training set
          inTrain <- train_idx[reTrain$Resample1]
          inTest = test_idx
        }else{ # ONLY TEST IS DEFINED
          test_idx = which(config[,test_vec[1], with=F][[1]] == test_vec[2])
          train_idx = setdiff(1:nrow(curr), test_idx) # use the other rows for testing
          reTrain <- caret::createDataPartition(y = config[train_idx, label], p = ml_train_perc) # take a user-defined percentage of the regexed training set
          inTrain <- train_idx[reTrain$Resample1]
          inTest <- test_idx
        }

        # choose predictor "label" (some others are also included but cross validation will be done on this)
        predictor = "label"

        # split training and testing data
        trainY <- curr[inTrain,
                       ..predictor][[1]]
        testY <- curr[inTest,
                      ..predictor][[1]]

        # remove predictive column from training set
        #remove.cols <- c("label") #TODO: make group column removed at start
        training <- curr[inTrain,]#-"label"]
        testing <- curr[inTest,]#-"label"]

        rmv.configCols.tr <- which(sapply(configCols, function(i, training) ifelse(length(levels(as.factor(training[,..i][[1]]))) == 1, TRUE, FALSE), training = training))
        training[, (rmv.configCols.tr) := NULL]

        rmv.configCols.te <- which(sapply(configCols, function(i, testing) ifelse(length(levels(as.factor(testing[,..i][[1]]))) == 1, TRUE, FALSE), testing = testing))
        testing[, (rmv.configCols.te) := NULL]

        # ======= WIP =======

        require(caret)

        # all methods
        caret.mdls <- getModelInfo()
        caret.methods <- names(caret.mdls)
        tune.opts <- lapply(caret.methods, function(mdl) caret.mdls[[mdl]]$parameters)
        names(tune.opts) <- caret.methods

        if(input$ml_folds == "LOOCV"){
          trainCtrl <- trainControl(verboseIter = T,
                                    allowParallel = F,
                                    method="LOOCV") # need something here...

        }else{
          trainCtrl <- trainControl(verboseIter = T,
                                    allowParallel = F,
                                    method=as.character(input$ml_perf_metr),
                                    number=as.numeric(input$ml_folds),
                                    repeats=3) # need something here...
        }

        meth.info <- caret::getModelInfo()[[input$ml_method]]
        params = meth.info$parameters

        #grid.def <- meth.info$grid(training, trainY, len = 1)

        tuneGrid = expand.grid(
          {
            lst = lapply(1:nrow(params), function(i){
              info = params[i,]
              inp.val = input[[paste0("ml_", info$parameter)]]
              # - - check for ranges - -
              if(grepl(inp.val, pattern=":")){
                split = strsplit(inp.val,split = ":")[[1]]
                inp.val <- seq(as.numeric(split[1]),
                               as.numeric(split[2]),
                               as.numeric(split[3]))
              }else if(grepl(inp.val, pattern = ",")){
                split = strsplit(inp.val,split = ",")[[1]]
                inp.val <- split
              }
              # - - - - - - - - - - - - -
              switch(as.character(info$class),
                     numeric = as.numeric(inp.val),
                     character = as.character(inp.val))
            })
            names(lst) = params$parameter
            if(any(sapply(lst,function(x)all(is.na(x))))){
              cat("Missing param, auto-tuning...")
              lst <- list()
              }
            #lst <- lst[sapply(lst,function(x)all(!is.na(x)))]
            lst
          })

        fit <- train(
          label ~ .,
          data = training,
          method = input$ml_method,
          ## Center and scale the predictors for the training
          ## set and all future samples.
          preProc = input$ml_preproc,
          tuneGrid = if(nrow(tuneGrid) > 0) tuneGrid else NULL,
          trControl = trainCtrl
        )

        result.predicted.prob <- predict(fit, testing, type="prob") # Prediction

        data = list(predictions = result.predicted.prob$`Population`, labels = testing$label)

        # train and cross validate model
        # return list with mode, prediction on test data etc.s
         list(type = input$ml_method,
              model = fit,
              prediction = result.predicted.prob[,2],
              labels = testing$label)
      })
    },
    train_vec = lcl$vectors$ml_train,
    test_vec = lcl$vectors$ml_test,
    configCols = configCols
    ) # for session_cl
    
    # check if a storage list for machine learning results already exists
    if(!"ml" %in% names(mSet$analSet)){
      mSet$analSet$ml <<- list() # otherwise make it
    }
    
    # save the summary of all repeats (will be used in plots) TOO MEMORY HEAVY
    pred <- ROCR::prediction(lapply(repeats, function(x) x$prediction), 
                             lapply(repeats, function(x) x$labels))
    perf <- ROCR::performance(pred, "tpr", "fpr")
    perf_auc <- ROCR::performance(pred, "auc")
    
    perf.long <- data.table::rbindlist(lapply(1:length(perf@x.values), function(i){
      xvals <- perf@x.values[[i]]
      yvals <- perf@y.values[[i]]
      aucs <- signif(perf_auc@y.values[[i]][[1]], digits = 2)
      
      res <- data.table::data.table(attempt = c(i),
                                    FPR = xvals,
                                    TPR = yvals,
                                    AUC = aucs)
      res
    }))
    
    mean.auc <- mean(unlist(perf_auc@y.values))
    
    roc_data <- list(m_auc = mean.auc,
                     perf = perf.long)
    
    # roc_data <- list(type = {unique(lapply(repeats, function(x) x$type))},
    #                  models = {lapply(repeats, function(x) x$model)},
    #                  predictions = {lapply(repeats, function(x) x$prediction)},
    #                  labels = {lapply(repeats, function(x) x$labels)})

    bar_data <- rbindlist(lapply(1:length(repeats), function(i){
      x = repeats[[i]]
      tbl = as.data.table(varImp(x$model)$importance, keep.rownames=T)
      tbl$rep = c(i)
      colnames(tbl) = c("mz",
                        "importance",
                        "rep")
      # - - -
      tbl
    }))

    # save results to mset
    mSet$analSet$ml[[input$ml_method]][[input$ml_name]] <<- list("roc" = roc_data,
                                                                 "bar" = bar_data)
    mSet$analSet$ml$last <<- list(name = input$ml_name,
                                  method = input$ml_method)

    # render plots for UIs
    datamanager$reload <- "ml"
  })
})

# mummichog

observeEvent(input$do_mummi, {

  peak_tbl <- if(mSet$dataSet$cls.num == 2){
    if("tt" %in% names(mSet$analSet)){
      continue = T
      data.table(
        `p.value` = mSet$analSet$tt$sig.mat[,"p.value"],
        `m.z` = rownames(mSet$analSet$tt$sig.mat),
        `t.score` = mSet$analSet$tt$sig.mat[,if("V" %in% colnames(mSet$analSet$tt$sig.mat)) "V" else "t.stat"]
      )
    }else{continue=F;
    NULL}
  }else{
    if("aov" %in% names(mSet$analSet)){
      continue = T
      data.table(
        `p.value` = mSet$analSet$aov$sig.mat[,"p.value"],
        `m.z` = rownames(mSet$analSet$aov$sig.mat),
        `t.score` = mSet$analSet$aov$sig.mat[,"F.stat"]
      )
    }else{continue=F;
    NULL}
  }

  if(!continue) NULL

  # seperate in pos and neg peaks..
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), lcl$paths$patdb)
  pospeaks <- DBI::dbGetQuery(conn, "SELECT DISTINCT mzmed FROM mzvals WHERE foundinmode = 'positive'")
  negpeaks <- DBI::dbGetQuery(conn, "SELECT DISTINCT mzmed FROM mzvals WHERE foundinmode = 'negative'")
  peak_tbl_pos <- peak_tbl[`m.z` %in% unlist(pospeaks)]
  peak_tbl_neg <- peak_tbl[`m.z` %in% unlist(negpeaks)]
  DBI::dbDisconnect(conn)

  for(mode in c("positive", "negative")){

    path <- tempfile()
    fwrite(x = peak_tbl, file = path, sep = "\t")
    mummi<-InitDataObjects("mass_all", "mummichog", FALSE)
    mummi<-Read.PeakListData(mSetObj = mummi, filename = path);
    mummi<-UpdateMummichogParameters(mummi, as.character(input$mummi_ppm), mode, input$mummi_sigmin);
    mummi<-SanityCheckMummichogData(mummi)
    mummi<-PerformMummichog(mummi, input$mummi_org, "fisher", "gamma")

    lcl$vectors[[paste0("mummi_", substr(mode, 1, 3))]] <<- list(sig = mummi$mummi.resmat,
                                                                    pw2cpd = {
                                                                      lst = mummi$pathways$cpds
                                                                      names(lst) <- mummi$pathways$name
                                                                      # - - -
                                                                      lst
                                                                    },
                                                                    cpd2mz = mummi$cpd2mz_dict)
    lcl$tables$mummichog <<- mummi$mummi.resmat
    output[[paste0("mummi_", substr(mode, 1, 3), "_tab")]] <- DT::renderDataTable({
      DT::datatable(mummi$mummi.resmat,selection = "single")
    })
  }
  output[[paste0("mummi_detail_tab")]] <- DT::renderDataTable({
    DT::datatable(data.table("no pathway selected"="Please select a pathway!"))
  })
})
