library(caret)
library(data.table)
library(xlsx)
library(glmnet)
library(pbapply)
library(parallel)

# ============== OTHER ================

brazil <- fread(input="/Users/jwolthuis/Analysis/SP/Brazil1_filled.csv",header = T)
spain <- fread(input="/Users/jwolthuis/Analysis/SP/Spain1_filled.csv",header = T)
brasp <-  fread(input="/Users/jwolthuis/Analysis/SP/BrazilAndSpain_filled.csv",header = T)

brazil <- fread(input="/Users/jwolthuis/Analysis/SP/Brazil_W.csv",header = T)
spain <- fread(input="/Users/jwolthuis/Analysis/SP/Spain_W.csv",header = T)
brasp <-  fread(input="/Users/jwolthuis/Analysis/SP/BrazilAndSpain_W.csv",header = T)
#ned <- fread(input="/Users/jwolthuis/Analysis/SP/Ned2.csv",header = T)

spain[,(1:ncol(spain)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]
brazil[,(1:ncol(brazil)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]
#ned[is.na(ned)] <- 0
brasp[,(1:ncol(brasp)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]

# ==================================

set.seed(849)

cl <- makeCluster(3, "FORK")

run_simulation <- function(curr, amin=0, amax=1, predictor = "Label"){
  
  nvars <- length(curr[,1:(which(colnames(curr) == "Injection"))])
  
  # ----------- remove qcs ----------
  
  curr <- curr[which(!grepl(curr$Sample, pattern = "QC"))]
  
  # ---------------------------------
  
  print(curr[,..predictor])
  
  curr[,predictor] <- as.factor(c(curr[,predictor,
                                       with=FALSE])[[1]])
  #curr <- data.table(curr)
  
  inTrain <- createDataPartition(y = c(curr[,predictor, 
                                            with=FALSE])[[1]],
                                 ## the outcome data are needed
                                 p = .6,
                                 ## The percentage of data in the training set
                                 list = FALSE)
  trainY <- curr[inTrain, 
                 predictor, 
                 with=FALSE][[1]]
  
  family <- if(length(levels(as.factor(trainY))) > 2){"multinomial"} else("binomial")
  
  training <- curr[ inTrain,]
  testing <- curr[-inTrain, ]
  
  trainX <- apply(training[, -c(1:nvars), with=FALSE], 2, as.numeric)
  testX <- apply(testing[, -c(1:nvars), with=FALSE], 2, as.numeric)
  
  testY <- testing[,predictor, 
                   with=FALSE][[1]]
  
  # Using glmnet to directly perform CV
  
  a <- seq(amin, amax, 0.1)
  
  search <- rbindlist(pbapply::pblapply(a,  cl=cl, FUN=function(i){
    cv <- cv.glmnet(trainX, trainY, family = family, nfold = 5, type.measure = "auc", alpha = i)
    data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
  }))
  
  my_alpha <- which(search$cvm == min(search$cvm))
  cv_final <- search[my_alpha, ]
  md <- glmnet(as.matrix(trainX), trainY, family = family, lambda = cv2$lambda.1se, alpha = cv2$alpha)
  
  
  print(prediction)
  
  if(length(levels(as.factor(trainY))) == 2){

    pred <- ROCR::prediction(prediction, 
                             testY)
    # calculate probabilities for TPR/FPR for predictions
    perf <- ROCR::performance(pred,"tpr","fpr")
    auc <- ROCR::performance(pred,"auc")@y.values[[1]] # shows calculated AUC for model
    plot(perf@x.values[[1]], perf@y.values[[1]])
    co <- md_final$beta[,1]
    
  }else{
    perf = NULL
    auc <- c(pROC::multiclass.roc(response = rep(testY,3), 
                         predictor = c(prediction), 
                         levels = c(1,2,3))$auc)
    co <- md_final$beta[[1]][,1]
    
  }

  # ---------------------
  
  inds <- which(co > 0)
  print(co[inds])
  variables <- names(co)[inds]
  variables <- variables[!(variables %in% '(Intercept)')]
  variables <- variables[order(variables, decreasing = TRUE)]
  
  # --------------
  list(auc = auc, alpha=a[my_alpha], roc = perf, cpds = variables, model = md_final, testX = as.matrix(testX), test_y = testY)
}

run_simulation_all_alpha <- function(curr, amin=0, amax=1, predictor = "Label", l = "min"){
  
  nvars <- length(curr[,1:(which(colnames(curr) == "Injection"))])
  
  # first 4 are identifiers, leave alone 
  

  for (i in (4:nvars)) {
    print(i)
    curr[,i] <- as.numeric(as.factor(curr[,..i][[1]]))-1
  }
  
  # ----------- remove qcs ----------
  
  curr <- curr[which(!grepl(curr$Sample, pattern = "QC"))]
  
  # ---------------------------------
  
  curr[,predictor] <- as.factor(c(curr[,..predictor])[[1]])

  inTrain <- createDataPartition(y = c(curr[,..predictor])[[1]],
                                 ## the outcome data are needed
                                 p = .6,
                                 ## The percentage of data in the training set
                                 list = FALSE)
  trainY <- curr[inTrain, 
                 ..predictor][[1]]
  
  family <- if(length(levels(as.factor(trainY))) > 2){"multinomial"} else("binomial")
  print(family)
  
  training <- curr[ inTrain,]
  testing <- curr[-inTrain, ]
  
  predIdx <- which(colnames(curr) == predictor | colnames(curr) %in% c("Group","Sample","$","X","Injection","Faecal_consistency_score"))
  
  trainX <- apply(training[, -c(1:4,predIdx), with=FALSE], 2, as.numeric)
  testX <- apply(testing[, -c(1:4,predIdx), with=FALSE], 2, as.numeric)

  print(colnames(trainX)[1:20])
  print(colnames(testX)[1:20])
  
  testY <- testing[,predictor, 
                   with=FALSE][[1]]
  
  # Using glmnet to directly perform CV
  
  a <- seq(amin, amax, 0.2)
  
  models <- pbapply::pblapply(a,  cl=cl, FUN=function(i){
    
    print(i)
    lambda = paste0('lambda.', l)
    
    cv1 <- cv.glmnet(trainX, trainY, family = family, nfold = 5, type.measure = "auc", alpha = i)
    cv2 <- data.frame(cvm = cv1$cvm[cv1$lambda == cv1[[lambda]]], lambda = cv1[[lambda]], alpha = i)
    # --------

    md <- glmnet(as.matrix(trainX), trainY, family = family, lambda = cv2$lambda, alpha = cv2$alpha)
    
    # --- roc ---
    
    prediction <- predict(md,
                          type = "response", 
                          newx = testX, 
                          s = lambda)
    
    # roc = switch(family,
    #              binomial = {
    #                roc <- ROCR::prediction(prediction,
    #                                        testY)  
    #              }, multinomial = {
    #                
    #              })
    # --------

    list(alpha = i, lambda = cv2[[lambda]], model = md, prediction = prediction, labels = testY)
  })
  models
}
# https://www.r-bloggers.com/variable-selection-with-elastic-net/ <<-- used this tutorial

res_brasp_stl <- run_simulation_all_alpha(brasp, 0, 1, "Stool_condition")
res_bra_stl <- run_simulation_all_alpha(brazil, 0, 1, "Stool_condition")
res_sp_stl <- run_simulation_all_alpha(spain, 0, 1, "Stool_condition")

res_brasp_bth <- run_simulation_all_alpha(brasp, 0, 1, "Batch")
res_bra_bth <- run_simulation_all_alpha(brazil, 0, 1, "Batch")
res_sp_bth <- run_simulation_all_alpha(spain, 0, 1, "Batch")

res_bra_frm <- run_simulation_all_alpha(brazil, 0, 1, "Farm")
res_sp_frm <- run_simulation_all_alpha(spain, 0, 1, "Farm")

res_brasp_cnt <- run_simulation_all_alpha(brasp, 0, 1, "Country")



plot.many <- function(res.obj, which_alpha = 1){
  predictions <- do.call("cbind", lapply(res.obj, function(x) x$prediction))
  colnames(predictions) <- sapply(res.obj, function(x) x$alpha)
  testY = res.obj[[1]]$labels
  
  if(length(unique(testY)) > 2){
    
    return("not supported yet")
    # https://stats.stackexchange.com/questions/112383/roc-for-more-than-2-outcome-categories
    #
    # for (type.id in length(unique(testY))) {
    #   type = as.factor(iris.train$Species == lvls[type.id])
    #   
    #   nbmodel = NaiveBayes(type ~ ., data=iris.train[, -5])
    #   nbprediction = predict(nbmodel, iris.test[,-5], type='raw')
    #   
    #   score = nbprediction$posterior[, 'TRUE']
    #   actual.class = iris.test$Species == lvls[type.id]
    #   
    #   pred = prediction(score, actual.class)
    #   nbperf = performance(pred, "tpr", "fpr")
    #   
    #   roc.x = unlist(nbperf@x.values)
    #   roc.y = unlist(nbperf@y.values)
    #   lines(roc.y ~ roc.x, col=type.id+1, lwd=2)
    #   
    #   nbauc = performance(pred, "auc")
    #   nbauc = unlist(slot(nbauc, "y.values"))
    #   aucs[type.id] = nbauc
    # }
  }else{
    data <- data.frame(D = as.numeric(as.factor(testY))-1,
                       D.str = testY)
    data <- cbind(data, predictions)
    roc_coord <- plotROC::melt_roc(data, "D", m = 3:ncol(data))
  }
  names(roc_coord)[which(names(roc_coord) == "name")] <- "alpha"
  
  roc_coord <- roc_coord[roc_coord$alpha %in% which_alpha,]
  # plot
  ggplot2::ggplot(roc_coord, 
                  ggplot2::aes(d = D, m = M, color = alpha)) + 
    plotROC::geom_roc(labelsize=0,show.legend = TRUE) + 
    plotROC::style_roc() + 
    theme(axis.text=element_text(size=19),
          axis.title=element_text(size=19,face="bold"),
          legend.title=element_text(size=19),
          legend.text=element_text(size=19))
}

plot.many(res_brasp_stl, which_alpha = c(0, 1))
plot.many(res_bra_stl)
plot.many(res_sp_stl)

#plot.many(res_brasp_bth)
plot.many(res_bra_bth)
plot.many(res_sp_bth)
plot.many(res_brasp_cnt, which_alpha = c(0,1))


cpds <- res_brasp_w_stl[[4]]$model$beta
used <- which(cpds > 0)
mz_ab <- rownames(cpds)[used]

cpds <- res_bra_stl[[3]]$model$beta
used <- which(cpds > 0)
mz_a <- rownames(cpds)[used]

cpds <- res_sp_stl[[2]]$model$beta
used <- which(cpds > 0)
mz_b <- rownames(cpds)[used]

intersect(mz_ab, mz_a)
# --- look up joined ones in db ---

db_search_list = paste0("/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/", c("kegg", "chebi", "hmdb", "internal", "noise", "metacyc", "wikipathways"), ".full.db")

cpds <- mz_ab

searchid = "mz"

cpd_info <- pbapply::pblapply(cpds, FUN=function(cpd){
  print(cpd)
  res <- R.utils::withTimeout(multimatch(cpd, db_search_list, searchid), timeout = 3, onTimeout = "warning")
  res
})

for(i in 1:length(cpd_info)){
  print(i)
  print(cpds[[i]])
  print(cpd_info[[i]])
  readline(prompt="Press [enter] to continue")
}

