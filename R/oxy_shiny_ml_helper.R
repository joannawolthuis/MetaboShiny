pcaCorr <- function(curr, center, scale, start_end_pcs){
  res <- prcomp(curr, 
                center = center,
                scale = scale)
  pc.use <- as.numeric(start_end_pcs[1]:start_end_pcs[2]) # explains 93% of variance
  trunc <- res$x[,pc.use] %*% t(res$rotation[,pc.use])
  as.data.frame(trunc)
}

getMLperformance = function(ml_res, pos.class, 
                            x.metric, y.metric,
                            silent = F){
  
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
                  folds,
                  maximize=T,
                  cl=0,
                  shuffle = F,
                  n_permute = 10,
                  shuffle_mode = "train",
                  silent = F,
                  tmpdir,
                  use_slurm=F){
  
  # get user training percentage
  need.rm = c("split")
  training[, (need.rm) := NULL]
  testing <- testing[, colnames(testing) != "split"]
  
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
  
  # save train and test
  train_fn = normalizePath(file.path(tmpdir, basename(tempfile())))
  test_fn = normalizePath(file.path(tmpdir, basename(tempfile())))
  qs::qsave(training, train_fn)
  qs::qsave(testing, test_fn)
  
  iterations = length(trainOrders)
  params = data.frame(
    train_fn = rep(train_fn, iterations),
    test_fn = rep(test_fn, iterations),
    ml_perf_metr = rep(ml_perf_metr, iterations),
    ml_folds = rep(ml_folds, iterations),
    ml_method = rep(ml_method, iterations),
    #ml_preproc = rep(ml_preproc, iterations),
    maximize = rep(maximize, iterations),
    trainOrder = I(trainOrders))
  
  has_slurm = Sys.getenv("SLURM_CPUS_ON_NODE") != ""
  
  folds <<- folds
  tuneGrid <<- tuneGrid
  
  if(has_slurm & use_slurm){
    batch_job = rslurm::slurm_apply(ml_single_run,
                                    params = params, 
                                    nodes = iterations,
                                    cpus_per_node = 1,
                                    pkgs = c("MetaboShiny",
                                             "caret",
                                             "data.table"),
                                    global_objects = c("folds", "tuneGrid"),
                                    slurm_options = list(time = "00:30:00"))
    completed = F
    
    print("Waiting on cluster to finish jobs...")
    
    while(!completed){
      Sys.sleep(5)
      completed = slurm_job_complete(batch_job)
    }
    
    # my ver has a progress bar
    print("Cluster batch job complete! Collecting results...")
    results <- get_slurm_out_jw(batch_job,outtype = "raw")
    #rslurm::cleanup_files(batch_job) #cleanup files
    
  }else{
    print("parallel mode")
    if(length(cl) > 0){
      doParallel::registerDoParallel(cl)  
    }
    results = lapply(trainOrders, 
                     ml_single_run,
                     train_fn = train_fn, 
                     test_fn = test_fn,
                     ml_perf_metr = ml_perf_metr,
                     ml_folds = ml_folds,
                     ml_method = ml_method,
                     ml_preproc = ml_preproc,
                     maximize = maximize)
  }
  # train and cross validate model
  # return list with mode, prediction on test data etc.s
  results
}

ml_single_run <- function(trainOrder,
                          train_fn,
                          test_fn,
                          ml_perf_metr, 
                          ml_folds,
                          ml_method,
                          ml_preproc=NULL,
                          maximize){
  training = qs::qread(train_fn)
  testing = qs::qread(test_fn)
  
  training$label <- as.factor(training$label)
  testing$label <- as.factor(testing$label)
  
  levels(training$label) = make.names(levels(training$label))
  levels(testing$label) = make.names(levels(testing$label))
  
  hasProb = !is.null(caret::getModelInfo(paste0("^",ml_method,"$"),regex = T)[[1]]$prob)
  
  orig.lbl = training[['label']]
  training[['label']] <- as.factor(training[['label']][trainOrder])
  def_scoring = ifelse(ifelse(is.factor(training[["label"]]), 
                              "Accuracy", "RMSE") %in% c("RMSE", "logLoss", "MAE"), 
                       FALSE,
                       TRUE)
  success=F
  
  shuffled = !all(training$label == orig.lbl)
  ml_folds <- if(!is.null(folds)) length(folds) else ml_folds
  
  trainCtrl <- caret::trainControl(verboseIter = T,
                                   allowParallel = T,
                                   method = if(ml_folds == "LOOCV") "LOOCV" else as.character(ml_perf_metr),
                                   number = as.numeric(ml_folds),
                                   trim = TRUE, 
                                   returnData = FALSE,
                                   classProbs = if(is.null(caret::getModelInfo(paste0("^",ml_method,"$"),regex = T)[[1]]$prob)) FALSE else TRUE,
                                   index = folds,
                                   returnResamp = "all",
                                   savePredictions = "final")
  if(ml_method == "glm"){
    fit <- caret::train(
      label ~ .,
      data = training,
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
      data = training,
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
                                          testing,
                                          type = if(hasProb) "prob" else "raw") # Prediction
  
  train.performance = fit$pred

  caret.mdls <- caret::getModelInfo()
  has.importance = "varImp" %in% names(caret.mdls[[ml_method]])
  
  l <- list(
    type = ml_method,
    best.model = fit$bestTune,
    train.performance = train.performance,
    importance = if(has.importance) caret::varImp(fit)$importance else data.table::data.table(unavailable = "method has no variable importance!"),
    labels = testing$label,
    in_test = rownames(testing),
    prediction = result.predicted.prob,
    shuffled = shuffled)
  return(l)
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
  
  if(length(settings$ml_batch_covars) == 0) settings$ml_batch_covars <- c("")
  if(settings$ml_batch_covars %in% c("", " ")) settings$ml_batch_covars <- c()
  
  # covars needed
  keep.config = setdiff(unique(c(settings$ml_include_covars, settings$ml_batch_covars,
                                 settings$ml_train_subset[[1]], settings$ml_test_subset[[1]])),
                        c("label", "", " "))
 
  sample_names = mSet$dataSet$covars$sample
  
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
      covars = c("label", setdiff(settings$ml_batch_covars, "individual"))
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
                       config = config[train_idx,,drop=F],
                       samples = sample_names[train_idx])
  
  testing_data = list(curr = curr[test_idx,,drop=F],
                      config = as.data.frame(config[test_idx,,drop=F]),
                      samples = sample_names[test_idx])
  
  rownames(testing_data$curr) <- rownames(testing_data$config) <- mSet$dataSet$covars$sample[test_idx]
  
  curr = NULL
  
  if(length(settings$ml_batch_covars) == 0){
    settings$ml_batch_balance <- settings$ml_batch_sampling <- F
  }
  
  # resampling
  
  if(settings$ml_sampling != "none"){
    print(settings$ml_sampling)
    # split on factor (either batch or placeholder to create one result)
    if(settings$ml_batch_balance){
      split.fac = training_data$config[, setdiff(settings$ml_batch_covars, "individual"), with=F][[1]]
    }else{
      split.fac = rep(1, nrow(training_data$config))
    }
    split.fac = if(settings$ml_batch_balance){
      training_data$config[, setdiff(settings$ml_batch_covars, "individual"), with=F][[1]]
    }else{rep(1, nrow(training_data$config))
    }
    
    spl.testing.idx = split(1:nrow(training_data$curr), split.fac)
    balance.overview = table(split.fac)
    biggest.group.overall = max(balance.overview)
    smallest.group.overall = min(balance.overview)
    
    orig.samp.distr = table(training_data$config$label)
    
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
      
      # equal groups but still upsampling wanted...
      if(length(unique(curr.group.sizes))==1 & settings$ml_sampling == "rose"){
        if(group.size > unique(curr.group.sizes)){
          # remove first sample
          rmv.smp = 1 #sample(1:nrow(curr.subset),size = 1)
          curr.subset = curr.subset[-rmv.smp,,drop=F]
          config.subset = config.subset[-rmv.smp,, drop=F]   
        }
      }
      
      ## K for the k-fold methods
      K = min(min(table(config.subset$label))-1, 10)
      mz.names = colnames(curr.subset)
      
      if(ncol(curr.subset) > 0){
        colnames(curr.subset) <- paste0("mz",1:ncol(curr.subset))
      }
      
      train.rownames = training_data$config$sample
      
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
        config.subset[[setdiff(settings$ml_batch_covars,"individual")]] <- c(config.top.row[[setdiff(settings$ml_batch_covars,
                                                                                                     "individual")]])
      }
      
      list(curr = curr.subset,
           config = config.subset)
    })
    training_data = list(
      curr = data.table::rbindlist(lapply(resampled.data.list, function(x) x$curr),use.names = T),
      config = data.table::rbindlist(lapply(resampled.data.list, function(x) x$config), use.names = T)
      #samples = training.rownames
      )
    resampled.data.list <- NULL
  }
  setwd(olddir)
  list(train = training_data, test = testing_data)
}

ml_run <- function(settings, mSet, input, cl, tmpdir, use_slurm = F){
  res = list()
  #({
  {
    
    if(length(settings$ml_batch_covars) == 1 & 
       settings$ml_batch_covars[1] == ""){
      settings$ml_batch_covars <- c()
      settings$ml_batch_balance = F
    }
    
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
    
    settings$in_train = training_data$samples
    
    if(settings$ml_specific_mzs != "no"){
      if(settings$ml_used_table != "pca"){
        msg = "Using user-specified m/z set."
        if(settings$ml_specific_mzs != "none"){
          if(length(settings$ml_mzs) > 0){
            training_data$curr <- training_data$curr[, ..settings$ml_mzs]  
          }else{
            mzs = getAllHits(mSet = mSet,
                             expname = settings$ml_specific_mzs,
                             randomize = settings$ml_mzs_rand)
            #print(mzs)
            #print(settings$ml_mzs_topn)
            
            mzs = mzs$m.z[1:settings$ml_mzs_topn]
            #print(mzs[1:10])
            
            mzs = gsub("^X", "", mzs)
            #print(mzs[1:10])
            
            mzs = gsub("\\.$", "-", mzs)
            #print(mzs[1:10])
            training_data$curr <- as.data.frame(training_data$curr[, mzs])
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
    
    # divvy folds based on batch (doensn't work for now)    
    fold_variable = if(length(settings$ml_batch_covars) > 0){
      as.factor(apply(training_data$config[, unique(c("label", 
                                                      settings$ml_batch_covars)),with=F],
                      MARGIN = 1,
                      function(x) paste0(x, collapse="_")))
    }else{
      training_data$config$label
    }
    
    folds <- if(settings$ml_folds == "LOOCV"){
      NULL
    } else {
       if(!settings$ml_covar_fold_seperate & length(settings$ml_batch_covars) == 0){
        # caret::createFolds(fold_variable, 
        #                    k = as.numeric(settings$ml_folds),
        #                    list = TRUE,returnTrain = T)
        NULL
      }else if(length(settings$ml_batch_covars) > 0){
        nfold = min(as.numeric(settings$ml_folds),
                    length(table(fold_variable)))
        caret::groupKFold(fold_variable,
                          k = nfold)
      }
    }
    
    print(folds)
    
    if(!is.null(folds)){
      fold_translation <- lapply(1:length(folds), function(i){
        idxs = folds[[i]]
        fold_variable[idxs]
      })  
      names(fold_translation) <- names(folds)
    }else{
      fold_translation = NULL
    }
    
    settings$samps_per_fold = fold_translation
    
    training_data$config = training_data$config[, ..keep.config]
    testing_data$config = testing_data$config[, keep.config]
    
    # merge back into one
    training = cbind(training_data$config,
                     training_data$curr)
    testing = cbind(testing_data$config[, colnames(training_data$config)],
                    testing_data$curr)
    
    if(!is.null(settings$ml_mtry)){
      if(!is.na(settings$ml_mtry)){
        if(settings$ml_mtry == "sqrt"){
          settings$ml_mtry = as.character(ceiling(sqrt(ncol(training_data$curr))))
        } 
      }
    }
    
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
                   folds = folds,
                   maximize = T,
                   shuffle = settings$ml_label_shuffle,
                   n_permute = settings$ml_n_shufflings,
                   shuffle_mode = if(settings$ml_shuffle_mode) "train" else "test",
                   cl = cl,
                   tmpdir=tmpdir,
                   use_slurm = use_slurm)
    
    res = list(res = ml_res, params = settings)
  }
  res
}

get_slurm_out_jw <- function (slr_job, outtype = "raw", wait = TRUE, ncores = NULL) 
{
  if (!(class(slr_job) == "slurm_job")) {
    stop("slr_job must be a slurm_job")
  }
  outtypes <- c("table", "raw")
  if (!(outtype %in% outtypes)) {
    stop(paste("outtype should be one of:", paste(outtypes, 
                                                  collapse = ", ")))
  }
  if (!(is.null(ncores) || (is.numeric(ncores) && length(ncores) == 
                            1))) {
    stop("ncores must be an integer number of cores")
  }
  if (wait) {
    rslurm:::wait_for_job(slr_job)
  }
  res_files <- paste0("results_", 0:(slr_job$nodes - 1), ".RDS")
  tmpdir <- paste0("_rslurm_", slr_job$jobname)
  missing_files <- setdiff(res_files, dir(path = tmpdir))
  if (length(missing_files) > 0) {
    missing_list <- paste(missing_files, collapse = ", ")
    warning(paste("The following files are missing:", missing_list))
  }
  res_files <- file.path(tmpdir, setdiff(res_files, missing_files))
  if (length(res_files) == 0) 
    return(NA)
  if (is.null(ncores)) {
    slurm_out <- pbapply::pblapply(res_files, readRDS)
  }
  else {
    cl = parallel::makeCluster(ncores)
    slurm_out <- pbapply::pblapply(res_files, readRDS, cl = cl)
    parallel::stopCluster(cl)
  }
  slurm_out <- do.call(c, slurm_out)
  if (outtype == "table") {
    slurm_out <- as.data.frame(do.call(rbind, slurm_out))
  }
  slurm_out
}

slurm_job_complete <- function (slr_job) 
{
  if (!(class(slr_job) == "slurm_job")) 
    stop("input must be a slurm_job")
  squeue_out <- suppressWarnings(system(paste("squeue -n", 
                                              slr_job$jobname), intern = TRUE))
  queue <- read.table(text = squeue_out, header = TRUE)
  completed <- nrow(queue) == 0
  return(completed)
}

slurm_apply_metshi <- function (f, params, ..., jobname = NA, nodes = 2, cpus_per_node = 2, 
                                processes_per_node = cpus_per_node, preschedule_cores = TRUE, 
                                global_objects = NULL, add_objects = NULL, pkgs = rev(.packages()), 
                                libPaths = NULL, rscript_path = NULL, r_template = NULL, 
                                sh_template = NULL, slurm_options = list(), submit = TRUE){
  print("altered MetShi slurm_apply")
  #print("a")
  if (!is.function(f)) {
    stop("first argument to slurm_apply should be a function")
  }
  if (!is.data.frame(params)) {
    stop("second argument to slurm_apply should be a data.frame")
  }
  if (is.null(names(params)) || (!is.primitive(f) && !"..." %in% 
                                 names(formals(f)) && any(!names(params) %in% names(formals(f))))) {
    stop("column names of params must match arguments of f")
  }
  if (!is.numeric(nodes) || length(nodes) != 1) {
    stop("nodes should be a single number")
  }
  if (!is.numeric(cpus_per_node) || length(cpus_per_node) != 
      1) {
    stop("cpus_per_node should be a single number")
  }
  if (!missing("add_objects")) {
    warning("Argument add_objects is deprecated; use global_objects instead.", 
            .call = FALSE)
    global_objects <- add_objects
  }
  if (is.null(r_template)) {
    r_template <- system.file("templates/slurm_run_R.txt", 
                              package = "rslurm")
  }
  if (is.null(sh_template)) {
    sh_template <- system.file("templates/submit_sh.txt", 
                               package = "rslurm")
  }
  #print("b")
  jobname <- rslurm:::make_jobname(jobname)
  tmpdir <- paste0("_rslurm_", jobname)
  dir.create(tmpdir, showWarnings = F)
  #print(normalizePath(tmpdir))
  more_args <- list(...)
  #print("c")
  paramfile = file.path(tmpdir, "params.RDS")
  saveRDS(params, file = paramfile)
  #print("c2")
  argfile = file.path(tmpdir, "more_args.RDS")
  saveRDS(more_args, file = argfile)
  #print("c3")
  funcfile = file.path(tmpdir, "f.RDS")
  #print(funcfile)
  saveRDS(f, file = funcfile)
  #print("d")
  if (!is.null(global_objects)) {
    save(list = global_objects, file = file.path(tmpdir, 
                                                 "add_objects.RData"), envir = environment(f))
  }
  if (nrow(params) < cpus_per_node * nodes) {
    nchunk <- cpus_per_node
  }
  else {
    nchunk <- ceiling(nrow(params)/nodes)
  }
  #print("e")
  nodes <- ceiling(nrow(params)/nchunk)
  template_r <- readLines(r_template)
  script_r <- whisker::whisker.render(template_r, list(pkgs = pkgs, 
                                                       add_obj = !is.null(global_objects), nchunk = nchunk, 
                                                       cpus_per_node = cpus_per_node, processes_per_node = processes_per_node, 
                                                       preschedule_cores = preschedule_cores, libPaths = libPaths))
  writeLines(script_r, file.path(tmpdir, "slurm_run.R"))
  template_sh <- readLines(sh_template)
  slurm_options <- rslurm:::format_option_list(slurm_options)
  #print("f")
  if (is.null(rscript_path)) {
    rscript_path <- file.path(R.home("bin"), "Rscript")
  }
  script_sh <- whisker::whisker.render(template_sh, list(max_node = nodes - 
                                                           1, cpus_per_node = cpus_per_node, jobname = jobname, 
                                                         flags = slurm_options$flags, options = slurm_options$options, 
                                                         rscript = rscript_path))
  writeLines(script_sh, file.path(tmpdir, "submit.sh"))
  #print("g")
  if (submit && system("squeue", ignore.stdout = TRUE)) {
    submit <- FALSE
    cat("Cannot submit; no Slurm workload manager found\n")
  }
  if (submit) {
    jobid <- rslurm:::submit_slurm_job(tmpdir)
  }
  else {
    jobid <- NA
    cat(paste("Submission scripts output in directory", 
              tmpdir, "\n"))
  }
  #print("h")
  return(rslurm:::slurm_job(jobname, jobid, nodes))
}

ml_slurm <- function(i, 
                     mloc, 
                     sloc,
                     tmpdir){
  small_mSet <- qs::qread(mloc)
  queue <- qs::qread(sloc)
  settings = queue[[i]]
  res = list()
  try({
    maxcores = if(Sys.getenv("SLURM_CPUS_ON_NODE") != ""){
      as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
    }
    #if(file.exists(logfile)) file.remove(logfile)
    ml_session_cl <- parallel::makeCluster(maxcores, 
                                           outfile="")#,setup_strategy = "sequential") # leave 1 core for general use and 1 core for shiny session
    # send specific functions/packages to other threads
    parallel::clusterEvalQ(ml_session_cl, {
      library(data.table)
      library(iterators)
      library(MetaboShiny)
      library(MetaDBparse)
    }) 
    res = ml_run(settings = settings, 
                              mSet = small_mSet,
                              cl = ml_session_cl,
                              tmpdir = tmpdir, 
                              use_slurm = T)
  })
  res
}

#' @title Get performance for multi-comparison ML model
#' @description ROC curves can be a bit tricky for multivariate models. This evaluates each possible pair of categories to generate individual and average AUC.
#' @param model ML model
#' @return FPR,TPR,average AUC,AUC for a given pair, and the name of the comparison
#' @seealso 
#'  \code{\link[pROC]{multiclass.roc}},\code{\link[pROC]{auc}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname getMultiMLperformance
#' @export 
#' @importFrom pROC multiclass.roc auc
#' @importFrom data.table rbindlist
getMultiMLperformance <- function(x, type="roc"){
  try({
    mroc = pROC::multiclass.roc(x$labels,
                                x$prediction)
  },silent = F)
  # try({
  #   mroc = pROC::multiclass.roc(x$labels, factor(x$prediction,
  #                                                ordered = T))
  # }, silent=F)
  
  data.table::rbindlist(lapply(mroc$rocs, function(roc.pair){
    try({
      dt = data.table(FPR = sapply(roc.pair$specificities, function(x) 1-x),
                      TPR = roc.pair$sensitivities,
                      AUC_AVG = as.numeric(mroc$auc),
                      AUC_PAIR = as.numeric(pROC::auc(roc.pair)),
                      comparison = paste0(roc.pair$levels,collapse=" vs. "))
    },silent=T)
    try({
      dt = data.table(FPR = sapply(roc.pair[[1]]$specificities, function(x) 1-x),
                      TPR = roc.pair[[1]]$sensitivities,
                      AUC_AVG = as.numeric(mroc$auc),
                      AUC_PAIR = as.numeric(pROC::auc(roc.pair[[1]])),
                      comparison = paste0(roc.pair[[1]]$levels,collapse=" vs. "))  
    },silent=T)
    dt
  })) 
}
