library(caret)
library(data.table)
library(xlsx)
library(glmnet)

# --------------

brazil <- fread(input="/Users/jwolthuis/Analysis/SP/Brazil1_noBatch.csv",header = T)
spain <- fread(input="/Users/jwolthuis/Analysis/SP/Spain1_noBatch.csv",header = T)

ned <- fread(input="/Users/jwolthuis/Analysis/SP/Ned2_noBatch.csv",header = T)

curr = ned

# ==================================

curr$Label <- as.factor(curr$Label)

inTrain <- createDataPartition(y = curr$Label,
                               ## the outcome data are needed
                                  p = .6,
                                ## The percentage of data in the training set
                                 list = FALSE)
training <- curr[ inTrain,]
testing <- curr[-inTrain,]

# train the model
set.seed(849)

control <- trainControl(method="LOOCV", 
                        number=10, 
                        returnResamp="all",
                        classProbs=TRUE, 
                        summaryFunction=twoClassSummary)

control <- trainControl(method = "repeatedcv", repeats = 5) #for ridge

# train the model
model <- train(x = training[,-c(1:2)], 
               y = training$Label,
               method = "glmnet",
               trControl = control)

pred <- predict(model, testing)

comp = data.table(known = testing$Label, prediction = pred)
pred_corr = length(which(comp$known == comp$prediction)) / nrow(comp) * 100.0

# estimate variable importance
importance <- varImp(model, scale=FALSE)

# summarize importance
print(importance)

# plot importance
plot(importance)

top20_brazil = importance$importance


# ============ PLAYING W/ LASSO AND ELASTICNET ============


cv_lasso <- cv.glmnet(as.matrix(trainX), 
                      trainY, 
                      family = "binomial", 
                      nfold = 10, 
                      type.measure = "deviance", 
                      alpha = 1)

md_lasso <- glmnet(as.matrix(trainX), 
                   trainY, 
                   family = "binomial", 
                   lambda = cv_lasso$lambda.1se, 
                   alpha = 1)
coef(md_lasso)

cv_ridge <- cv.glmnet(as.matrix(trainX), 
                      trainY, 
                      family = "binomial", 
                      nfold = 10, 
                      type.measure = "deviance", 
                      alpha = 0)

md_ridge <- glmnet(as.matrix(trainX), 
                   trainY, 
                   family = "binomial", 
                   lambda = cv_ridge$lambda.1se, 
                   alpha = 0)

# ============== OTHER ================

library(glmnet)
library(caret)
library(pbapply)

brazil <- fread(input="/Users/jwolthuis/Analysis/BR/Brazil1_noBatch.csv",header = T)
spain <- fread(input="/Users/jwolthuis/Analysis/SP/Spain1_noBatch.csv",header = T)

brazil[is.na(brazil)] <- 0
spain[is.na(spain)] <- 0


datasets = list(brazil, spain)

# ==================================

set.seed(849)

test_situation <- function(curr){
  
  #curr[is.na(curr)] <- 0
  
  # ---------------------------------
  curr$Label <- as.factor(curr$Label)
  
  inTrain <- createDataPartition(y = curr$Label,
                                 ## the outcome data are needed
                                 p = .6,
                                 ## The percentage of data in the training set
                                 list = FALSE)
  training <- curr[ inTrain,]
  testing <- curr[-inTrain,]
  
  trainX <- training[, -c(1,2)]
  testX <- testing[, -c(1,2)]
  trainY <- training$Label
  
  # Using glmnet to directly perform CV
  
  a <- seq(0.1, 0.9, 0.05)
  
  search <- rbindlist(pblapply(a, FUN=function(i){
    cv <- cv.glmnet(as.matrix(trainX), trainY, family = "binomial", nfold = 10, type.measure = "deviance", alpha = i)
    data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
  }))
  
  print("here")
  
  cv_final <- search[search$cvm == min(search$cvm), ]
  md_final <- glmnet(as.matrix(trainX), trainY, family = "binomial", lambda = cv_final$lambda.1se, alpha = cv_final$alpha)
  
  #md_final <- glmnet(as.matrix(trainX), trainY, family = "binomial", lambda = cv_final$lambda.1se, alpha = cv_final$alpha)
  
  pred <- predict(md_final, s='lambda.1se', newx = as.matrix(testX), type="class")
  comp <- data.table(known = testing$Label, prediction = pred)
  pred_corr <- length(which(comp$known == comp$prediction)) / nrow(comp) * 100.0
  
  co <- coef(md_final)
  inds <- which(co[,1]>0)
  variables <- row.names(co)[inds]
  variables <- variables[!(variables %in% '(Intercept)')]
  
  # --------------
  list(perc = pred_corr, cpds = variables, model = md_final, testX = as.matrix(testX), test_y = testing$Label)
}

# https://www.r-bloggers.com/variable-selection-with-elastic-net/ <<-- used this tutorial

attempts = 1

best_perc = list(b2s = 0, s2b = 0)
best_model = NULL

for(attempt in 1:attempts){
  
  res_bra = test_situation(brazil)
  
  # lasso.model <- cv.glmnet(x = res_bra$testX, y = res_bra$test_y, 
  #                          family = 'binomial', type.measure = 'auc',nfolds = 3)
  # 
  # # Apply model to testing dataset
  # res_bra$lasso.prob <- predict(lasso.model,type="response", 
  #                            newx = res_bra$testX, s = 'lambda.1se')
  # pred <- prediction(res_bra$lasso.prob, res_bra$test_y)
  # 
  # # calculate probabilities for TPR/FPR for predictions
  # perf <- performance(pred,"tpr","fpr")
  # performance(pred,"deviance") # shows calculated AUC for model
  # plot(perf,colorize=FALSE, col="black") # plot ROC curve
  # lines(c(0,1),c(0,1),col = "gray", lty = 4 )
  
  res_spa = test_situation(spain)
  
  brazil2spain = {
    pred = predict(res_bra$model, s='lambda.1se', newx = res_spa$testX, type="class")
    comp = data.table(known = res_spa$test_y, prediction = pred)
    pred_corr = length(which(comp$known == comp$prediction)) / nrow(comp) * 100.0
    # ------
    pred_corr
  }

  spain2brazil = {
    pred = predict(res_spa$model, s='lambda.1se', newx = res_bra$testX, type="class")
    comp = data.table(known = res_bra$test_y, prediction = pred)
    pred_corr = length(which(comp$known == comp$prediction)) / nrow(comp) * 100.0
    # ------
    pred_corr
  }  
  
  new_perc <- list(b2s = brazil2spain, s2b = spain2brazil)
  print(new_perc)
  
  if(brazil2spain > best_perc$b2s || spain2brazil > best_perc$s2b){
    print("yea")
    best_perc <<- new_perc
    best_model = list(br=res_bra, sp=res_spa)
  }
}

pred_bra_self = predict(best_model$br$model, s='lambda.1se', newx = best_model$br$testX, type="class")
comp_bra_self = data.table(known = best_model$br$test_y, prediction = pred_bra_self)
pred_corr_bra_self = length(which(comp_bra_self$known == comp_bra_self$prediction)) / nrow(comp_bra_self) * 100.0
# ------
pred_corr_bra_self

pred_spa_self = predict(best_model$sp$model, s='lambda.1se', newx = best_model$sp$testX, type="class")
comp_spa_self = data.table(known = best_model$sp$test_y, prediction = pred_spa_self)
pred_corr_spa_self = length(which(comp_spa_self$known == comp_spa_self$prediction)) / nrow(comp_spa_self) * 100.0
# ------
pred_corr_spa_self

in_both = intersect(best_model$br$cpds, best_model$sp$cpds)
in_both

# --- look up joined ones in db ---

db_search_list = paste0("/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/", c("kegg", "chebi", "hmdb", "internal", "noise", "smpdb", "wikipathways"), ".full.db")
in_both <- variables
results_in_both <- pblapply(in_both, FUN=function(cpd){
  match_list <- lapply(db_search_list, FUN=function(match.table){
    get_matches(cpd, match.table, searchid="baseformula")
  })
  match_table <<- unique(as.data.table(rbindlist(match_list))[Name != ""])
})

results_in_both[[1]]

for(i in 1:length(results_in_both)){
  print(i)
  print(in_both[[i]])
  print(results_in_both[[i]])
  readline(prompt="Press [enter] to continue")
}

# --- deconstruct other two ---

cpds_brazil = best_model$br$cpds
results_brazil <- pblapply(cpds_brazil, FUN=function(cpd){
  match_list <- lapply(db_search_list, FUN=function(match.table){
    get_matches(cpd, match.table, searchid="baseformula")
  })
  match_table <<- unique(as.data.table(rbindlist(match_list))[Name != ""])
})

for(i in 1:length(cpds_brazil)){
  print(i)
  print(cpds_brazil[[i]])
  print(results_brazil[[i]])
  readline(prompt="Press [enter] to continue")
}

cpds_spain = best_model$sp$cpds
results_spain <- pblapply(cpds_spain, FUN=function(cpd){
  match_list <- lapply(db_search_list, FUN=function(match.table){
    get_matches(cpd, match.table, searchid="baseformula")
  })
  match_table <<- unique(as.data.table(rbindlist(match_list))[Name != ""])
})

for(i in 1:length(cpds_spain)){
  print(i)
  print(cpds_spain[[i]])
  print(results_spain[[i]])
  readline(prompt="Press [enter] to continue")
}

save(best_perc, best_model, file="lasso_feb272018.RData")
