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
#ned <- fread(input="/Users/jwolthuis/Analysis/SP/Ned2.csv",header = T)

brazil[is.na(brazil)] <- 0
spain[is.na(spain)] <- 0
#ned[is.na(ned)] <- 0
brasp[is.na(brasp)] <- 0

# ==================================

set.seed(849)

cl <- makeCluster(3, "FORK")

test_situation <- function(curr, amin=0.1, amax=0.9, predictor = "Label"){
  
  # ---------------------------------
  
  curr[,predictor] <- as.factor(c(curr[,predictor,
                                       with=FALSE])[[1]])
  curr <- data.table(curr)
  
  inTrain <- createDataPartition(y = c(curr[, predictor, 
                                            with=FALSE])[[1]],
                                 ## the outcome data are needed
                                 p = .5,
                                 ## The percentage of data in the training set
                                 list = FALSE)
  trainY <- curr[inTrain, 
                 predictor, 
                 with=FALSE][[1]]
  
  family <- if(length(levels(as.factor(trainY))) > 2){"multinomial"} else("binomial")
  
  training <- curr[ inTrain,]
  testing <- curr[-inTrain, ]
  
  remove <-  c("Sample", "Label")
  
  trainX <- apply(training[, -remove, with=FALSE], 2, as.numeric)
  testX <- apply(testing[, -remove, with=FALSE], 2, as.numeric)
  testY <- testing[,predictor, 
                   with=FALSE][[1]]
  
  # Using glmnet to directly perform CV
  
  a <- seq(amin, amax, 0.05)
  
  search <- rbindlist(pbapply::pblapply(a,  cl=cl, FUN=function(i){
    cv <- cv.glmnet(as.matrix(trainX), trainY, family = family, nfold = 10, type.measure = "auc", alpha = i)
    data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
  }))
  
  my_alpha <- which(search$cvm == min(search$cvm))
  cv_final <- search[my_alpha, ]
  md_final <- glmnet(as.matrix(trainX), trainY, family = family, lambda = cv_final$lambda.1se, alpha = cv_final$alpha)
  
  # --- roc ---
  
  prediction <- predict(md_final,
                        type = "response", 
                        newx = testX, 
                        s = 'lambda.1se')
  
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
  variables <- names(co)[inds]
  variables <- variables[!(variables %in% '(Intercept)')]
  
  # --------------
  list(auc = auc, alpha=a[my_alpha], roc = perf, cpds = variables, model = md_final, testX = as.matrix(testX), test_y = testY)
}

# https://www.r-bloggers.com/variable-selection-with-elastic-net/ <<-- used this tutorial

res_brasp_batch <- test_situation(brasp, 0.1, 0.9, "Batch")
res_brasp_label <- test_situation(brasp, 0.1, 0.9, "Label")

res_brasp_label$cpds
res_brasp_batch$cpds


# --- look up joined ones in db ---

db_search_list = paste0("/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/", c("kegg", "chebi", "hmdb", "internal", "noise", "metacyc", "wikipathways"), ".full.db")

cpds <- res_brasp_label$cpds
res_brasp$alpha

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
