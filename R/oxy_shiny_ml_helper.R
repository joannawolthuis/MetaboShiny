pcaCorr <- function(curr, center, scale, start_end_pcs){
  res <- prcomp(curr, 
                center = center,
                scale = scale)
  pc.use <- as.numeric(start_end_pcs[1]:start_end_pcs[2]) # explains 93% of variance
  trunc <- res$x[,pc.use] %*% t(res$rotation[,pc.use])
  as.data.frame(trunc)
}

getMLperformance = function(ml_res, pos.class, x.metric, y.metric){
  
  if(!("Resample" %in% names(ml_res$train.performance))){
    is.loocv = TRUE
  }else{
    is.loocv = FALSE
  }
  
  if(is.loocv){
    print("Cannot estimate ROC for LOOCV 'folds'.")
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
                  silent = F){
  
  # get user training percentage
  need.rm = c("split")
  training[, (need.rm) := NULL]
  testing[, (need.rm) := NULL]
  
  # shuffle if shuffle
  trainOrders = list(1:nrow(training))
  if(shuffle & shuffle_mode == "train"){
    for(i in 1:n_permute){
      trainOrders = append(trainOrders, list(sample(1:nrow(training))))
    }
  }
  
  hasProb = !is.null(caret::getModelInfo(paste0("^",ml_method,"$"),regex = T)[[1]]$prob)
  
  iterations = length(trainOrders)
  
  results = lapply(trainOrders,  function(train.order,
                                          train.set = train.set,
                                          test.set = test.set)
  {
    orig.lbl = train.set[['label']]
    train.set[['label']] <- train.set[['label']][train.order] 
    def_scoring = ifelse(ifelse(is.factor(train.set[["label"]]), 
                                "Accuracy", "RMSE") %in% c("RMSE", "logLoss", "MAE"), 
                         FALSE,
                         TRUE)
    success=F
    
    shuffled = !all(train.set$label == orig.lbl)
    
    trainCtrl <- caret::trainControl(verboseIter = T,
                                     allowParallel = T,
                                     method = if(ml_folds == "LOOCV") "LOOCV" else as.character(ml_perf_metr),
                                     number = as.numeric(ml_folds),
                                     #repeats = 3,
                                     trim = TRUE, 
                                     returnData = FALSE,
                                     classProbs = if(is.null(caret::getModelInfo(paste0("^",ml_method,"$"),regex = T)[[1]]$prob)) FALSE else TRUE,
                                     #index = if(!shuffled) folds else NULL,
                                     savePredictions = "final")
    
    if(silent){
      msg = capture.output(
        fit <- caret::train(
          label ~ .,
          data = train.set,
          method = ml_method,
          ## Center and scale the predictors for the training
          ## set and all future samples.
          preProc = ml_preproc,
          maximize = if(maximize) def_scoring else !def_scoring,
          importance = if(ml_method %in% c("ranger")) 'permutation' else TRUE,
          tuneGrid = if(nrow(tuneGrid) > 0) tuneGrid else NULL,
          trControl = trainCtrl
        )
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
        importance = if(ml_method %in% c("ranger")) 'permutation' else TRUE,
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