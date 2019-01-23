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
    # prepare matrix
    
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
    
    #curr2=curr
    
    # join halves together, user variables and metabolite data
    curr <- cbind(config, curr)
    curr <- as.data.table(curr)
    
    # remove cols with all NA
    curr <- curr[,colSums(is.na(curr)) < nrow(curr), with=FALSE]
    # remove rows with all NA
    curr <- curr[complete.cases(curr),]
    # identify which columns are metabolites and which are config/covars
    configCols <- which(!(gsub(x = colnames(curr), pattern = "_T\\d", replacement="") %in% colnames(mSet$dataSet$norm)))
    mzCols <- which(gsub(x = colnames(curr), pattern = "_T\\d", replacement="") %in% colnames(mSet$dataSet$norm))
    
    # make the covars factors and the metabolites numeric.
    curr[,(configCols):= lapply(.SD, function(x) as.factor(x)), .SDcols = configCols]
    curr[,(mzCols):= lapply(.SD, function(x) as.numeric(x)), .SDcols = mzCols]
    
    # how many models will be built? user input
    goes = as.numeric(input$ml_attempts)
    
    # ============ LOOP HERE ============
    
    # get results for the amount of attempts chosen
    repeats <- pbapply::pblapply(1:goes, 
                                 cl=session_cl, 
                                 function(i, ...){
      shiny::isolate({
        
        # get regex user input for filtering testing and training set
        ml_train_regex <- "" #input$ml_train_regex
        ml_test_regex <- "" #input$ml_test_regex
        
        # get user training percentage
        ml_train_perc <- input$ml_train_perc/100
        
        if(ml_train_regex == "" & ml_test_regex == ""){ # BOTH ARE NOT DEFINED
          test_idx = caret::createDataPartition(y = curr$label, p = ml_train_perc, list = FALSE) # partition data in a balanced way (uses labels)
          train_idx = setdiff(1:nrow(curr), test_idx) #use the other rows for testing
          inTrain = train_idx
          inTest = test_idx
        }else if(ml_train_regex != ""){ #ONLY TRAIN IS DEFINED
          train_idx = grep(config$sample, pattern = ml_train_regex) # get training sample ids with regex
          test_idx = setdiff(1:nrow(curr), train_idx) # use the other rows for testing
          reTrain <- caret::createDataPartition(y = config[train_idx, label], p = ml_train_perc) # take a user-defined percentage of the regexed training set
          inTrain <- train_idx[reTrain$Resample1]
          inTest = test_idx
        }else{ # ONLY TEST IS DEFINED
          test_idx = grep(config$sample, pattern = ml_test_regex) # get training sample ids with regex
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
        remove.cols <- c("label") #TODO: make group column removed at start
        remove.idx <- which(colnames(curr) %in% remove.cols)
        training <- curr[inTrain,-remove.idx, with=FALSE]
        testing <- curr[inTest,-remove.idx, with=FALSE]
        
        # get covar indices
        predIdx <- which(colnames(curr) %in% colnames(config))
        
        # remove unused levels in the factor part of the table after filtering
        training <- data.matrix(gdata::drop.levels(training))
        testing <- data.matrix(gdata::drop.levels(testing))
        
        #shiny::setProgress(value = i/goes)
        
        # train and cross validate model
        switch(input$ml_method,
               rf = { # random forest
                 model = randomForest::randomForest(x = training, 
                                                    y = trainY, 
                                                    ntree = 500, # amount of trees made TODO: make user choice
                                                    importance=TRUE) # include variable importance in model
                 
                 prediction <- stats::predict(model, 
                                              testing, 
                                              "prob")[,2]
                 # get importance table
                 importance = as.data.table(model$importance, 
                                            keep.rownames = T)
                 rf_tab <- importance[which(MeanDecreaseAccuracy > 0), c("rn", 
                                                                         "MeanDecreaseAccuracy")]
                 rf_tab <- rf_tab[order(MeanDecreaseAccuracy, decreasing = T)] # reorder for convenience
                 rf_tab <- data.frame(MDA = rf_tab$MeanDecreaseAccuracy, row.names = rf_tab$rn) 
                 # return list with model, prediction on test data etc.
                 list(type="rf",
                      feats = as.data.table(rf_tab, 
                                            keep.rownames = T), 
                      model = model,
                      prediction = prediction,
                      labels = testY)
               }, 
               ls = { # lasso
                 # user x cross validation input
                 nfold = switch(input$ml_folds, 
                                "5" = 5,
                                "10" = 10,
                                "20" = 20,
                                "50" = 50,
                                "LOOCV" = length(trainY)) # leave one out CV
                 
                 family = "binomial" #TODO: enable multinomial for multivariate data!!
                 
                 #training <- scale(training, T, T)
                 #testing <- scale(testing, T, T)
                 
                 #  make model (a. internal cross validation until optimized)
                 cv1 <- glmnet::cv.glmnet(training, 
                                          trainY, 
                                          family = family, 
                                          type.measure = "auc", 
                                          alpha = 1, 
                                          keep = TRUE, 
                                          nfolds=nfold)
                 # pick the best model (other options, lambda.1se)
                 cv2 <- data.frame(cvm = cv1$cvm[cv1$lambda == cv1[["lambda.min"]]], lambda = cv1[["lambda.min"]], alpha = 1)
                 
                 # save final model
                 model <- glmnet::glmnet(as.matrix(training), 
                                         trainY, 
                                         family = family, 
                                         lambda = cv2$lambda, 
                                         alpha = cv2$alpha)
                 
                 # test on testing data and save prediction
                 prediction <- stats::predict(model,
                                              type = "response", 
                                              newx = testing, 
                                              s = "lambda.min")#[,2] # add if necessary
                 # return list with mode, prediction on test data etc.s
                 list(type = "ls",
                      model = model,
                      prediction = prediction, 
                      labels = testY)
               }, 
               gls = {
                 NULL
               })
      })
    })#, input, config, curr) # for session_cl
    # check if a storage list for machine learning results already exists
    if(!"ml" %in% names(mSet$analSet)){
      
      mSet$analSet$ml <<- list(ls=list(), 
                               rf=list()) # otherwise make it
      
      
    }
    # save the summary of all repeats (will be used in plots)
    roc_data <- list(type = {unique(lapply(repeats, function(x) x$type))},
                     models = {lapply(repeats, function(x) x$model)},
                     predictions = {lapply(repeats, function(x) x$prediction)},
                     labels = {lapply(repeats, function(x) x$labels)})
    
    bar_data <- switch(input$ml_method,
                       rf = {
                         res <- aggregate(. ~ rn, rbindlist(lapply(repeats, function(x) as.data.table(x$feats, keep.rownames=T))), mean)
                         data <- res[order(res$MDA, decreasing = TRUE),]
                         colnames(data) <- c("mz", "mda")
                         # - - -
                         data
                       },
                       ls = {
                         feat_count <- lapply(repeats, function(x){
                           beta <- x$model$beta
                           feats <- which(beta[,1] > 0)
                           names(feats)
                         })
                         feat_count_tab <- table(unlist(feat_count))
                         feat_count_dt <- data.table::data.table(feat_count_tab)
                         colnames(feat_count_dt) <- c("mz", "count")
                         data <- feat_count_dt[order(feat_count_dt$count, decreasing = T)]
                         # - - -
                         data
                       })
    bar_data$mz <- factor(bar_data$mz, levels=bar_data$mz)
    
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
        `t.score` = mSet$analSet$tt$sig.mat[,"t.stat"]
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
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), global$paths$patdb)
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
    
    global$vectors[[paste0("mummi_", substr(mode, 1, 3))]] <<- list(sig = mummi$mummi.resmat,
                                                                    pw2cpd = {
                                                                      lst = mummi$pathways$cpds
                                                                      names(lst) <- mummi$pathways$name
                                                                      # - - -
                                                                      lst
                                                                    },
                                                                    cpd2mz = mummi$cpd2mz_dict)
    global$tables <<- mummi$mummi.resmat
    output[[paste0("mummi_", substr(mode, 1, 3), "_tab")]] <- DT::renderDataTable({
      DT::datatable(mummi$mummi.resmat,selection = "single")
    })
  }
  output[[paste0("mummi_detail_tab")]] <- DT::renderDataTable({
    DT::datatable(data.table("no pathway selected"="Please select a pathway!")) 
  })
})