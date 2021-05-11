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
    spl.fold.performance = split(ml_res$train.performance, ml_res$train.performance$Resample)
    coord.collection = lapply(spl.fold.performance, function(l){
      
      prediction = l[[pos.class]]
      labels = l[["obs"]]
      prediction = ROCR::prediction(prediction,
                                    labels)
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
runML <- function(curr,
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
                  n_permute = 10){
  
  # get user training percentage
  if(unique(train_vec)[1] != "all"){ #ONLY TRAIN IS DEFINED
    train_idx <- which(curr[,train_vec[1], with=F][[1]] == train_vec[2])
    test_idx = setdiff(1:nrow(curr), train_idx) # use the other rows for testing
    inTrain <- train_idx
    inTest = test_idx
  }else{ # ONLY TEST IS DEFINED
    test_idx = which(curr[,test_vec[1], with=F][[1]] == test_vec[2])
    train_idx = setdiff(1:nrow(curr), test_idx) # use the other rows for testing
    inTrain = train_idx
    inTest <- test_idx
  }
  
  trainSamps = rownames(curr)[inTrain]
  testSamps = rownames(curr)[inTest]
  
  need.rm = unique(c(train_vec[1],test_vec[1],"split"))
  curr <- curr[,-..need.rm]
  
  # choose predictor "label" (some others are also included but cross validation will be done on this)
  predictor = "label"
  
  # split training and testing data
  trainY <- curr[inTrain,
                 ..predictor][[1]]
  testY <- curr[inTest,
                ..predictor][[1]]
  
  training <- curr[inTrain,]
  testing <- curr[inTest,]
  
  hasProb = !is.null(caret::getModelInfo(paste0("^",ml_method,"$"),regex = T)[[1]]$prob)
  
  if(cl != 0){
    doParallel::registerDoParallel(cl)
  }
  
  trainCtrl <- caret::trainControl(verboseIter = T,
                                   allowParallel = T,
                                   method = if(ml_folds == "LOOCV") "LOOCV" else as.character(ml_perf_metr),
                                   number = as.numeric(ml_folds),
                                   repeats = 3,
                                   trim = TRUE, 
                                   returnData = FALSE,
                                   classProbs = if(!hasProb) FALSE else TRUE,
                                   index = folds,
                                   savePredictions = "all")
  
  # shuffle if shuffle
  trainOrders = list(1:nrow(training))
  if(shuffle){
    for(i in 1:n_permute){
      trainOrders = append(trainOrders, list(sample(1:nrow(training))))
    }
  }
  # ------------------
  reps = pbapply::pblapply(trainOrders, function(ordr){
    training[['label']] = training[['label']][ordr]
    def_scoring = ifelse(ifelse(is.factor(training[["label"]]), "Accuracy", "RMSE") %in% c("RMSE", "logLoss", "MAE"), FALSE, TRUE)
    success=F
    fit <- caret::train(
      label ~ .,
      data = training,
      method = ml_method,
      ## Center and scale the predictors for the training
      ## set and all future samples.
      preProc = ml_preproc,
      maximize = if(maximize) def_scoring else !def_scoring,
      importance = if(ml_method %in% c("ranger")) 'permutation' else TRUE,
      tuneGrid = if(nrow(tuneGrid) > 0) tuneGrid else NULL,
      trControl = trainCtrl
    )
    result.predicted.prob <- stats::predict(fit, 
                                            testing,
                                            type = if(hasProb) "prob" else "raw") # Prediction
    list(type = ml_method,
         train.performance = fit$pred,
         importance = caret::varImp(fit)$importance,
         prediction = result.predicted.prob,
         labels = testing$label,
         distr = list(train = trainSamps, test = testSamps),
         shuffled = !all(ordr == 1:nrow(training)))
  })
  # train and cross validate model
  # return list with mode, prediction on test data etc.s
  reps
}