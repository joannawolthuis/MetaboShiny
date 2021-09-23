pcaCorr <- function(curr, center, scale, start_end_pcs){
  res <- prcomp(curr, 
                center = center,
                scale = scale)
  pc.use <- as.numeric(start_end_pcs[1]:start_end_pcs[2]) # explains 93% of variance
  trunc <- res$x[,pc.use] %*% t(res$rotation[,pc.use])
  as.data.frame(trunc)
}

getMLperformance = function(ml_res, pos.class, x.metric, y.metric, silent = F){
  
  if(!("Resample" %in% names(ml_res$train.performance))){
    is.loocv = TRUE
  }else{
    is.loocv = FALSE
  }
  
  if(is.loocv){
    if(!silent) print("Cannot estimate train performance.")
    coord.collection = list()
  }else{
    spl.fold.performance = split(ml_res$train.performance,
                                 ml_res$train.performance$Resample)
    coord.collection = lapply(spl.fold.performance, function(l){
      prediction = ROCR::prediction(predictions = l[,pos.class],
                                    labels = l$obs)
      coords = ROCR::performance(prediction,
                                 x.measure = x.metric,
                                 measure = y.metric)
      coords
    })  
  }
  
  # alternative precision recall...
  
  
  
  # -------------------------------
  
  
  prediction = ROCR::prediction(ml_res$prediction[,pos.class], 
                                ml_res$labels)
  
  coords = ROCR::performance(prediction,
                             x.measure = x.metric,
                             measure = y.metric) 
  
  coord.collection$Test = coords
  
  coords.dt = data.table::rbindlist(lapply(1:length(coord.collection), function(i){
    coords = coord.collection[[i]]
    which.test = names(coord.collection)[[i]]
    data.table::data.table(x = coords@x.values[[1]],
                           y = coords@y.values[[1]],
                           cutoff = coords@alpha.values[[1]],
                           `Test set` = c(which.test))
  }))
  
  for(col in c("x", "y")){
    if(Inf %in% coords.dt[[col]]){
      if(!silent) print("WARNING: Inf changed to 1")
      coords.dt[[col]][coords.dt[[col]] == Inf] <- 1
    }
    if(NaN %in% coords.dt[[col]]){
      if(!silent) print("WARNING: NaN changed to 1")
      coords.dt[[col]][is.nan(coords.dt[[col]])] <- 1
    }
  }
  
  list(coords = coords.dt,
       names = list(x = coords@x.name,
                    y = coords@y.name,
                    alpha = "Cutoff"))
}

#' @title Run machine learning
#' @description Large wrapper function to run machine learning
#' @param curr Non-normalized peak table (generally mSet$dataSet$proc or something similar)
#' @param config Configuration table (metadata table)
#' @param train_vec Metadata group to train on
#' @param test_vec Metadata group to test on
#' @param configCols Columns representing metadata
#' @param ml_method caret ML method
#' @param ml_perf_metr Performance measuring method
#' @param ml_folds Cross validation folds
#' @param ml_preproc Preproc table object
#' @param tuneGrid Table of settings to try out to optimize parameters
#' @return List of machine learning model, importance, labels and prediction.
#' @seealso 
#'  \code{\link[caret]{createDataPartition}},\code{\link[caret]{trainControl}},\code{\link[caret]{train}},\code{\link[caret]{varImp}}
#'  \code{\link[stats]{predict}}
#' @rdname runML
#' @export 
#' @importFrom caret createDataPartition trainControl train varImp
#' @importFrom stats predict
runML <- function(training,
                  testing,
                  config,
                  train_vec,
                  test_vec,
                  ml_method,
                  ml_perf_metr,
                  ml_folds,
                  ml_preproc,
                  tuneGrid,
                  fold_variable,
                  folds,
                  maximize=T,
                  cl=0,
                  shuffle = F,
                  n_permute = 10,
                  shuffle_mode = "train",
                  silent = F){
  
  # get user training percentage
  need.rm = c("split")
  training[, (need.rm) := NULL]
  testing[, (need.rm) := NULL]
  
  if(length(cl) > 0){
    doParallel::registerDoParallel(cl)  
  }

  # shuffle if shuffle
  trainOrders = list(1:nrow(training))
  if(shuffle & shuffle_mode == "train"){
    for(i in 1:n_permute){
      shuffled_order = sample(1:nrow(training))
      trainOrders = append(trainOrders, list(shuffled_order))
    }
  }
  
  if(ml_method == "glm (logistic)"){
    ml_method = "glm"
    is_logit = T
  }else{
    is_logit = F
  }
  
  hasProb = !is.null(caret::getModelInfo(paste0("^",ml_method,"$"),regex = T)[[1]]$prob)
  
  iterations = length(trainOrders)
  
  results = lapply(trainOrders,  function(train.order,
                                          train.set = train.set,
                                          test.set = test.set)
  {
    orig.lbl = train.set[['label']]
    train.set[['label']] <- as.factor(train.set[['label']][train.order])
    def_scoring = ifelse(ifelse(is.factor(train.set[["label"]]), 
                                "Accuracy", "RMSE") %in% c("RMSE", "logLoss", "MAE"), 
                         FALSE,
                         TRUE)
    success=F
    
    shuffled = !all(train.set$label == orig.lbl)
    
    ## does not divide classes properly :(
    # if(length(fold_variable) == 0){
    #   fold_variable = train.set[['label']]
    # }else{
    #   fold_variable = paste0(train.set[['label']], fold_variable)
    # }
    # 
    # folds <- if(ml_folds != "LOOCV") caret::groupKFold(fold_variable,
    #                                                    k = min(as.numeric(ml_folds), 
    #                                                            length(unique(fold_variable)))) else NULL
    trainCtrl <- caret::trainControl(verboseIter = T,
                                     allowParallel = T,
                                     method = if(ml_folds == "LOOCV") "LOOCV" else as.character(ml_perf_metr),
                                     number = as.numeric(ml_folds),
                                     trim = TRUE, 
                                     returnData = FALSE,
                                     classProbs = if(is.null(caret::getModelInfo(paste0("^",ml_method,"$"),regex = T)[[1]]$prob)) FALSE else TRUE,
                                     #index = folds,
                                     savePredictions = "final")
    if(ml_method == "glm"){
      fit <- caret::train(
        label ~ .,
        data = train.set,
        method = ml_method,
        ## Center and scale the predictors for the training
        ## set and all future samples.
        preProc = ml_preproc,
        maximize = if(maximize) def_scoring else !def_scoring,
        tuneGrid = if(nrow(tuneGrid) > 0) tuneGrid else NULL,
        #trControl = trainCtrl,
        family = if(is_logit) "binomial" else NULL
      )  
    }else{
      fit <- caret::train(
        label ~ .,
        data = train.set,
        method = ml_method,
        ## Center and scale the predictors for the training
        ## set and all future samples.
        preProc = ml_preproc,
        maximize = if(maximize) def_scoring else !def_scoring,
        importance = if(ml_method == c("ranger")) 'permutation' else TRUE,
        tuneGrid = if(nrow(tuneGrid) > 0) tuneGrid else NULL,
        trControl = trainCtrl
      )
    }
    
    result.predicted.prob <- stats::predict(fit, 
                                            test.set,
                                            type = if(hasProb) "prob" else "raw") # Prediction
    
    l <- list(#model = fit,
      type = ml_method,
      best.model = fit$bestTune,
      train.performance = fit$pred,
      importance = caret::varImp(fit)$importance,
      labels = testing$label,
      distr = list(train = rownames(train.set),
                   test = rownames(test.set)),
      prediction = result.predicted.prob,
      shuffled = shuffled)
    return(l)
  }, train.set = training, test.set = testing)
  
  # train and cross validate model
  # return list with mode, prediction on test data etc.s
  results
}

ml_prep_data <- function(settings, mSet, input, cl){
  olddir = getwd()
  tmpdir = file.path(tempdir()) # needed for 
  if(!dir.exists(tmpdir)) dir.create(tmpdir)
  setwd(tmpdir)
  
  ### PIPELINE ###
  # pick source table
  pickedTbl <- settings$ml_used_table
  
  if(pickedTbl == "pca" & !("pca" %in% names(mSet$analSet))){
    stop("Please run PCA first!")
  }
  
  # covars needed
  keep.config = setdiff(unique(c(settings$ml_include_covars, settings$ml_batch_covars,
                                 settings$ml_train_subset[[1]], settings$ml_test_subset[[1]])),
                        "label")
  if(length(keep.config) > 0){
    config = mSet$dataSet$covars[, ..keep.config,drop=F]
  }else{
    config = data.table::data.table()
  }
  config$label = mSet$dataSet$cls
  
  # train/test split
  ## get indices
  if(!is.null(settings$ml_train_subset) | !is.null(settings$ml_test_subset)){
    # add clause for same train_test
    test_idx = NULL
    train_idx = NULL
    if(!is.null(settings$ml_test_subset)){
      test_idx = which(config[[settings$ml_test_subset[[1]]]] %in% settings$ml_test_subset[[2]])
    }
    if(!is.null(settings$ml_train_subset)){
      train_idx = which(config[[settings$ml_train_subset[[1]]]] %in% settings$ml_train_subset[[2]])
    }
    if(is.null(train_idx)){
      train_idx = setdiff(1:nrow(config), test_idx)  
    }else if(is.null(test_idx)){
      test_idx = setdiff(1:nrow(config), train_idx)
    }
  }else{
    # make joined label of label+batch and split that train/test
    split_label = if(length(settings$ml_batch_covars) > 0){
      #print("splitting tr/te % per batch")
      covars = c("label", settings$ml_batch_covars)
      apply(config[, ..covars],
            MARGIN = 1, 
            function(x) paste0(x, collapse="_"))
    }else{
      #print("unbiased split over pool")
      config$label
    }
    train_idx = caret::createDataPartition(y = split_label, 
                                           p = settings$ml_train_perc/100,
                                           list = FALSE)[,1] # partition data in a balanced way (uses labels)
    
    test_idx = setdiff(1:nrow(config), train_idx)
    #table(split_label[test_idx])
  }
  
  if(mSet$metshiParams$renorm & pickedTbl != "pca"){
    samps_train = rownames(mSet$dataSet$norm)[train_idx]
    samps_test = rownames(mSet$dataSet$norm)[test_idx]
    mSet.settings = mSet$settings
    reset_mSet <- reset.mSet(mSet,
                             fn = file.path(lcl$paths$proj_dir, 
                                            paste0(lcl$proj_name,
                                                   "_ORIG.metshi")))
    # --- GET TRAIN ---
    mSet_train = subset_mSet(reset_mSet, "sample", samps_train)
    mSet_train = change.mSet(mSet_train, 
                             stats_var = mSet.settings$exp.var, 
                             time_var =  mSet.settings$time.var,
                             stats_type = mSet.settings$exp.type)
    
    mSet_train$dataSet$orig <- mSet_train$dataSet$start
    mSet_train$dataSet$start <- mSet_train$dataSet$preproc <- mSet_train$dataSet$proc <- mSet_train$dataSet$prenorm <- NULL
    mSet_train = metshiProcess(mSet_train, init = F, cl = 0)
    # --- GET TEST ---
    mSet_test = subset_mSet(reset_mSet, "sample", samps_test)
    mSet_test = change.mSet(mSet_test, 
                            stats_var = mSet.settings$exp.var, 
                            time_var =  mSet.settings$time.var,
                            stats_type = mSet.settings$exp.type)
    mSet_test$dataSet$orig <- mSet_test$dataSet$start
    mSet_test$dataSet$start <- mSet_test$dataSet$preproc <- mSet_test$dataSet$proc <- mSet_test$dataSet$prenorm <- NULL
    mSet_test = metshiProcess(mSet_test, init = F, cl = 0)
    
    # ------- rejoin and create curr -------
    config_train = mSet_train$dataSet$covars[, ..keep.config,drop=F]
    config_train$label = mSet_train$dataSet$cls
    config_test = mSet_test$dataSet$covars[, ..keep.config,drop=F]
    config_test$label = mSet_test$dataSet$cls
    
    mz.in.both = intersect(colnames(mSet_train$dataSet$norm),
                           colnames(mSet_test$dataSet$norm))
    
    curr = rbind(mSet_train$dataSet$norm[,mz.in.both,drop=F],
                 mSet_test$dataSet$norm[,mz.in.both,drop=F])
    config = rbind(config_train, 
                   config_test)
    
    mSet_test <- mSet_train <- config_test <- config_train <- mz.in.both <- mSet$storage <- NULL
  }else{
    curr = as.data.frame(switch(pickedTbl, 
                                orig = mSet$dataSet$orig,
                                norm = mSet$dataSet$norm,
                                pca = mSet$analSet$pca$x))
  }
  
  test_sampnames = rownames(curr)[test_idx]
  
  # PCA correct
  if(settings$ml_pca_corr & pickedTbl != 'pca'){
    #print("Performing PCA and subtracting PCs...")
    curr <- pcaCorr(curr, 
                    center = if(pickedTbl == "norm") F else T,
                    scale = if(pickedTbl == "norm") F else T, 
                    start_end_pcs = settings$ml_keep_pcs)
  }
  
  #@ split
  training_data = list(curr = curr[train_idx,,drop=F],
                       config = config[train_idx,,drop=F])
  
  testing_data = list(curr = curr[test_idx,,drop=F],
                      config = config[test_idx,,drop=F])
  
  curr = NULL
  
  if(length(settings$ml_batch_covars) == 0){
    settings$ml_batch_balance <- settings$ml_batch_sampling <- F
  }
  
  # resampling
  
  if(settings$ml_sampling != "none"){
    
    # split on factor (either batch or placeholder to create one result)
    if(settings$ml_batch_balance){
      split.fac = training_data$config[, settings$ml_batch_covars, with=F][[1]]
    }else{
      split.fac = rep(1, nrow(training_data$config))
    }
    split.fac = if(settings$ml_batch_balance){
      training_data$config[, settings$ml_batch_covars, with=F][[1]]
    }else{rep(1, nrow(training_data$config))
    }
    
    spl.testing.idx = split(1:nrow(training_data$curr), split.fac)
    balance.overview = table(split.fac)
    biggest.group.overall = max(balance.overview)
    smallest.group.overall = min(balance.overview)
    
    orig.samp.distr = table(training_data$config$label)
    
    training_data$config = training_data$config
    testing_data$config = testing_data$config
    
    size.global = if(settings$ml_sampling != "down") biggest.group.overall else smallest.group.overall
    size.preset = settings$ml_groupsize
    
    # resample
    resampled.data.list = lapply(spl.testing.idx, function(idx){
      
      curr.subset = training_data$curr[idx,,drop=F]
      config.subset = training_data$config[idx,,drop=F]
      config.top.row = config.subset[1,,drop=F]
      
      sampling = settings$ml_sampling
      
      size.local = if(sampling != "down") max(table(config.subset$label)) else min(table(config.subset$label))
      
      # all upsampling except "down"
      group.size = 
        if(settings$ml_batch_balance){
          if(settings$ml_batch_size_sampling){ if(size.preset > 0) size.preset else size.global} else size.local
        }else size.local
      
      curr.group.sizes = table(config.subset$label)
      
      if(sum(curr.group.sizes > group.size) == 1){
        majority.idxs = which(config.subset$label == names(which(curr.group.sizes > group.size)))
        minority.idxs = setdiff(1:nrow(config.subset), majority.idxs)
        keep.samps = c(minority.idxs, sample(majority.idxs, size = group.size))
        config.subset = config.subset[keep.samps,,drop=F]
        curr.subset = curr.subset[keep.samps,,drop=F]
      }
      
      # total group size within this loop?
      
      ## K for the k-fold methods
      K = min(min(table(config.subset$label))-1, 10)
      mz.names = colnames(curr.subset)
      
      if(ncol(curr.subset) > 0){
        colnames(curr.subset) <- paste0("mz",1:ncol(curr.subset))
      }
      
      switch(sampling,
             up = {
               nconfig = ncol(config.subset)
               new_data = upsample.adj(cbind(config.subset, curr.subset), 
                                       as.factor(config.subset$label), 
                                       maxClass = group.size) 
               config.subset = new_data[,1:nconfig,drop=F]
               curr.subset = new_data[,!(colnames(new_data) %in% c("Class", colnames(config.subset))),drop=F]
             },
             adasyn = {
               resampled = smotefamily::ADAS(X = curr.subset, 
                                             target = config.subset$label,
                                             K = K)
               new_data = resampled$data
               curr.subset = new_data[, colnames(new_data) != "class",drop=F]
               config.subset = data.table(label = new_data$class)
             },
             smote = {
               resampled = smotefamily::SMOTE(X = curr.subset, 
                                              target = config.subset$label,
                                              K = K)
               new_data = resampled$data
               curr.subset = new_data[, colnames(new_data) != "class",drop=F]
               config.subset = data.table(label = new_data$class)
             },
             rose = {
               rose.dat = cbind(label = config.subset$label, 
                                curr.subset)
               resampled = ROSE::ROSE(label ~ ., 
                                      data = rose.dat,
                                      N = group.size * length(unique(config.subset$label)))
               new_data = resampled$data
               curr.subset = new_data[,2:(ncol(new_data)),drop=F]
               config.subset = data.table(label = new_data[[1]])
             },
             down = {
               nconfig = ncol(config.subset)
               new_data = downsample.adj(cbind(config.subset, curr.subset), 
                                         as.factor(config.subset$label), 
                                         minClass = group.size) 
               config.subset = new_data[,1:nconfig,drop=F]
               curr.subset = new_data[,!(colnames(new_data) %in% c("Class", colnames(config.subset))),drop=F]
               
             }
      )
      colnames(curr.subset) <- mz.names
      
      if(settings$ml_batch_balance & settings$ml_sampling %in% c("rose", "smote", "adasyn")){
        config.subset[[settings$ml_batch_covars]] <- c(config.top.row[[settings$ml_batch_covars]])
      }
      
      list(curr = curr.subset,
           config = config.subset)
    })
    training_data = list(
      curr = data.table::rbindlist(lapply(resampled.data.list, function(x) x$curr),use.names = T),
      config = data.table::rbindlist(lapply(resampled.data.list, function(x) x$config), use.names = T)
    )
    resampled.data.list <- NULL
  }
  setwd(olddir)
  list(train = training_data, test = testing_data)
}

ml_run <- function(settings, mSet, input, cl){
  res = list()
  #({
  {
    if("for_ml" %in% names(mSet$dataSet)){
      jobi = mSet$dataSet$for_ml$mapper[ml_name == settings$ml_name,]$unique_data_id
      training_data = mSet$dataSet$for_ml$datasets[[jobi]]$train
      testing_data = mSet$dataSet$for_ml$datasets[[jobi]]$test
    }else{
      tr_te = ml_prep_data(settings = settings, 
                           mSet = mSet,
                           input = input, cl=0)
      training_data = tr_te$train
      testing_data = tr_te$test 
    }
    
    if(settings$ml_specific_mzs != "no"){
      if(pickedTbl != "pca"){
        msg = "Using user-specified m/z set."
        if(settings$ml_specific_mzs != "none"){
          if(length(settings$ml_mzs) > 0){
            training_data$curr <- training_data$curr[, ..settings$ml_mzs]  
          }else{
            mzs = getAllHits(mSet = mSet,
                             expname = settings$ml_specific_mzs,
                             randomize = settings$ml_mzs_rand)
            mzs = mzs$m.z[1:settings$ml_mzs_topn]
            mzs = gsub("^X", "", mzs)
            mzs = gsub("\\.$", "-", mzs)
            training_data$curr <- as.data.frame(training_data$curr[, ..mzs])
          }
        }else{
          curr <- data.table::data.table() # only use configs
        }
      }
    }
    
    # replace training data with the new stuff
    training_data$config$split <- "train"
    testing_data$config$split <- "test"
    
    # remove columns that should not be in prediction
    keep.config = unique(c("label", "split", settings$ml_include_covars))
    
    fold_variable = if(length(settings$ml_batch_covars) > 0){
      apply(training_data$config[, settings$ml_batch_covars,with=F],
            MARGIN = 1, 
            function(x) paste0(x, collapse="_"))
    }else{
      c()
    }
    
    testing_data$config = testing_data$config[, ..keep.config]
    training_data$config = training_data$config[, ..keep.config]
    
    # merge back into one
    training = cbind(training_data$config,
                     training_data$curr)
    testing = cbind(testing_data$config[,colnames(training_data$config),with=F],
                    testing_data$curr)
    
    training_data <- testing_data <- NULL
    
    # CARET SETTINGS
    caret.mdls <- caret::getModelInfo()
    caret.methods <- names(caret.mdls)
    tune.opts <- lapply(caret.methods, function(mdl) caret.mdls[[mdl]]$parameters)
    names(tune.opts) <- caret.methods
    
    meth.info <- caret.mdls[[settings$ml_method]]
    params = meth.info$parameters
    
    tuneGrid = if(length(params) == 0){
      data.frame()
    }else{
      expand.grid(
        {
          lst = lapply(1:nrow(params), function(i){
            info = params[i,]
            inp.val = settings[[paste0("ml_", info$parameter)]]
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
            #cat("Missing param, auto-tuning...")
            lst <- list()
          }
          lst
        })  
    }
    
    # make sure levels of predicted class aren't numeric
    levels(training$label) <- paste0("class",  Hmisc::capitalize(levels(training$label)))
    levels(training$label) <- ordered(levels(training$label))
    levels(testing$label) <- paste0("class",  Hmisc::capitalize(levels(testing$label)))
    levels(testing$label) <- ordered(levels(testing$label))
    
    # correct mzs in case model cannot handle numeric column names
    colnames(training) <- make.names(colnames(training))
    colnames(testing) <- make.names(colnames(testing))
    
    # run ML
    ml_res = runML(training = training,
                   testing = testing,
                   ml_method = settings$ml_method,
                   ml_perf_metr = settings$ml_perf_metr,
                   ml_folds = settings$ml_folds,
                   ml_preproc = settings$ml_preproc,
                   tuneGrid = tuneGrid,
                   fold_variable = fold_variable,
                   maximize = T,
                   shuffle = settings$ml_label_shuffle,
                   n_permute = settings$ml_n_shufflings,
                   shuffle_mode = if(settings$ml_shuffle_mode) "train" else "test",
                   cl = cl)
    
    res = list(res = ml_res, params = settings)
  }
  res
}
