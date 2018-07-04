# GOAL 1: SUBSAMPLING

caret::downSample()

library(caret)

data(oil)
table(oilType)
downSample(fattyAcids, oilType)

upSample(fattyAcids, oilType)

# - - - - - - - - - - - -
input <- list(ml_attempts = 50, 
              ml_saturation = TRUE, 
              ml_train_perc = 80, 
              ml_train_regex = "", ## TRAINING SET
              ml_test_regex = "", ## TESTING SET
              ml_name = "test",
              ml_method = "ls",
              ml_top_x = 10,
              ml_sampling = "up")

curr <- as.data.table(mSet$dataSet$preproc)

curr[,(1:ncol(curr)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]

config <- mSet$dataSet$batches[match(rownames(mSet$dataSet$preproc),mSet$dataSet$batches$Sample),]
config <- config[!is.na(config$Sample),]
config <- cbind(config, Label=mSet$dataSet$cls)

keep_curr <- match(mSet$dataSet$batches$Sample,rownames(mSet$dataSet$preproc))

curr <- cbind(config, curr[keep_curr])

curr <- curr[which(!grepl(curr$Sample,
                          pattern = "QC"))]

configCols <- which(!(colnames(curr) %in% colnames(mSet$dataSet$norm)))
mzCols <- which(colnames(curr) %in% colnames(mSet$dataSet$norm))

curr[,(configCols):= lapply(.SD, function(x) as.factor(x)), .SDcols = configCols]
curr[,(mzCols):= lapply(.SD, function(x) as.numeric(x)), .SDcols = mzCols]

print(input$ml_attempts)
print(input$ml_saturation)

if(input$ml_saturation){
  goes = 9
}else{
  goes = as.numeric(input$ml_attempts)
}

print(goes)

repeats <- pbapply::pblapply(1:goes, function(i){
  
  # - train / test - 
  
  ml_train_regex <<- input$ml_train_regex
  ml_test_regex <<- input$ml_test_regex
  
  old_idxes <- 1:nrow(config)
  config_adj <- cbind(config, idx=old_idxes)
  
  if(input$ml_saturation){
    print("Sat mode")
    ml_train_perc <- 0.1 * i
    #config_filled <- config_adj[caret::createDataPartition(y = config_adj[,Label],p = use_perc)$Resample1,]
    #ml_train_perc <- (input$ml_train_perc/100)*(.05*i)
  }else{
    ml_train_perc <- input$ml_train_perc/100
    classes <- as.factor(grepl(config_adj$Sample, pattern = input$ml_train_regex))
    downSampled <- downSample(config_adj, classes)
    upSampled <- upSample(config_adj, classes)
    
    # UPSAMPLE W/ SMOTE
    #ncols<-ncol(dm)
    #stats <- table(curr$Label)
    #minority <- which(stats == min(stats))
    #majority <- which(stats == max(stats))
    #stats[[minority]]/stats[[majority]] * 100
    
    curr2 <- DMwR::SMOTE(Label ~ . , 
                         curr,
                         k=5, 
                         perc.over = 120, 
                         perc.under = 100)
    
    table(curr$Label)
    table(curr2$Label)
    
    # - - - - - - - - -
    config_filled <- as.data.table(switch(input$ml_sampling, up=upSampled, down=downSampled))
  }
  
  print(i)
  print(ml_train_perc)
 
  if(input$ml_train_regex != ""){
    print("a")
    train_idx <- grep(config_filled$Sample, pattern = input$ml_train_regex)
    reTrain <- caret::createDataPartition(y = config_filled[train_idx, Label],p = ml_train_perc)
    inTrain <- train_idx[reTrain$Resample1]
    inTrain <- config_filled[inTrain, "idx"][[1]]
    
    if(input$ml_test_regex != ""){
      print("b")
      test_idx <- grep(config_filled$Sample, pattern = input$ml_test_regex)
      reTest <- caret::createDataPartition(y = config_filled[test_idx, Label],p = input$ml_train_perc/100)
      inTest <- test_idx[reTest$Resample1]
      inTest <- config_filled[inTest, "idx"][[1]]
    }else{
      inTest <- setdiff(1:nrow(curr), train_idx)
    }
  }else{
    if(input$ml_test_regex != ""){

      test_idx <- grep(config_adj$Sample, pattern = input$ml_test_regex)
      train_idx <- setdiff(1:nrow(config_adj), test_idx)

      reTrain <- caret::createDataPartition(y = config[train_idx, Label],p = ml_train_perc)
      inTrain <- train_idx[reTrain$Resample1]
      inTest <- config_adj[test_idx, "idx"][[1]]

    }else{
      print("c")
      inTrain <- caret::createDataPartition(y = config_filled$Label,
                                            ## the outcome data are needed
                                            p = ml_train_perc,#input$ml_train_perc/100,
                                            ## The percentage of data in the training set
                                            list = FALSE)
      
      inTest <- setdiff(1:nrow(config_filled), inTrain)
      
      inTrain <- config_filled[inTrain, "idx"][[1]]
      
      inTest <- config_filled[inTest, "idx"][[1]]
      
    }
    }
  
  print(inTrain)
  print(inTest)
  
  print(length(inTrain))
  print(length(inTest))
  # - - divide - -
  
  predictor = "Label"
  
  trainY <- curr[inTrain, 
                 ..predictor][[1]]
  testY <- curr[inTest,
                ..predictor][[1]]
  
  training <- curr[inTrain,-c("Sample", "Label")]
  testing <- curr[inTest,-c("Sample", "Label")]
  
  predIdx <- which(colnames(curr) %in% colnames(config))
  
  training <- data.matrix(gdata::drop.levels(training))
  testing <- data.matrix(gdata::drop.levels(testing))
  
  switch(input$ml_method,
         rf = {
           model = randomForest::randomForest(x = training, 
                                              y = trainY, 
                                              ntree = 500,
                                              importance=TRUE)
           
           prediction <- stats::predict(model, 
                                        testing, 
                                        "prob")[,2]
           
           importance = as.data.table(model$importance, keep.rownames = T)
           rf_tab <- importance[which(MeanDecreaseAccuracy > 0), c("rn", "MeanDecreaseAccuracy")]
           rf_tab <- rf_tab[order(MeanDecreaseAccuracy, decreasing = T)]
           rf_tab <- data.frame(MDA = rf_tab$MeanDecreaseAccuracy, row.names = rf_tab$rn) 
           list(type="rf",
                feats = as.data.table(rf_tab, keep.rownames = T), 
                model = model,
                prediction = prediction,
                labels = testY)
         }, 
         ls = {
           nfold = length(trainY)
           family = "binomial"
           
           cv1 <- glmnet::cv.glmnet(training, trainY, family = family, nfold = nfold, type.measure = "auc", alpha = 1, keep = TRUE)
           cv2 <- data.frame(cvm = cv1$cvm[cv1$lambda == cv1[["lambda.min"]]], lambda = cv1[["lambda.min"]], alpha = 1)
           
           model <- glmnet::glmnet(as.matrix(training), trainY, family = family, lambda = cv2$lambda, alpha = cv2$alpha)
           
           prediction <- stats::predict(model,
                                        type = "response", 
                                        newx = testing, 
                                        s = "lambda.min")#[,2] # add if necessary
           list(type = "ls",
                model = model,
                prediction = prediction, 
                labels = testY)
         }, 
         gls = {
           NULL
         })
})

ml_name <- input$ml_name

xvals <- list(type = {unique(lapply(repeats, function(x) x$type))},
               models = {lapply(repeats, function(x) x$model)},
               predictions = {lapply(repeats, function(x) x$prediction)},
               labels = {lapply(repeats, function(x) x$labels)})

plotly::ggplotly(ggPlotROC(xvals, input$ml_attempts, cf))
plotly::ggplotly(ggPlotBar(repeats, input$ml_attempts, cf, input$ml_top_x))


# GOAL 2: SATURATION CURVE

# set training and test set to static


# set seed for reproducibility
set.seed(7)

# randomize mtcars
curr2 <- curr[sample(nrow(curr)),]

configCols <- which(!(colnames(curr2) %in% colnames(mSet$dataSet$norm)))
mzCols <- which(colnames(curr2) %in% colnames(mSet$dataSet$norm))

curr2[,(configCols):= lapply(.SD, function(x) as.factor(x)), .SDcols = configCols]
curr2[,(mzCols):= lapply(.SD, function(x) as.numeric(x)), .SDcols = mzCols]

# split iris data into training and test sets
mtcarsIndex <- createDataPartition(config$Label, p = .8, list = F)

mtcarsTrain <- curr2[mtcarsIndex,]
mtcarsTest <- curr2[-mtcarsIndex,]

training <- data.matrix(gdata::drop.levels(mtcarsTrain))
testing <- data.matrix(gdata::drop.levels(mtcarsTest))


#create empty data frame
learnCurve <- data.frame(m = integer(40),
                         trainRMSE = integer(40),
                         cvRMSE = integer(40))

# test data response feature
trainY <- as.factor(training[,"Label"])

testY <- as.factor(testing[,"Label"])

mse = function(x,y) { mean((x-y)^2)}

# loop over training examples
for (i in 10:40) {
  
  learnCurve$m[i] <- i
  
  print(i)
  
  family = "binomial"
  
  rmv <- which(colnames(training) %in% c("Sample", "Label"))
  
  print(training[1:i,-rmv][,1:10])
  cv1 <- glmnet::cv.glmnet(training[1:i,-rmv], 
                           trainY[1:i], 
                           family = family, 
                           nfold = i, 
                           type.measure = "auc", 
                           alpha = 1, 
                           keep = TRUE)
  
  lambda.id <- which(cv1$lambda == cv1$lambda.min)
  mse.1 = mse(cv1$fit[,lambda.id], as.numeric(trainY[1:i])-1) 
  cat("MSE (cross-validation): ", mse.1, "\n")
  rmse.1 <- sqrt(mse.1)
  
  learnCurve$trainRMSE[i] <- rmse.1
  
  cv2 <- data.frame(cvm = cv1$cvm[cv1$lambda == cv1[["lambda.min"]]], lambda = cv1[["lambda.min"]], alpha = 1)
  
  model <- glmnet::glmnet(training[1:i,-rmv], trainY[1:i], family = family, lambda = cv2$lambda, alpha = cv2$alpha)
  
  prediction <- predict(model,
                        type = "response",
                        newx = testing[,-rmv],
                        s = "lambda.min")#[,2] # add if necessary
  
  error <- as.numeric(testY) - 1
  mse.2 = mse(prediction, as.numeric(testY)-1) 
  cat("MSE (testing): ", mse.2, "\n")
  rmse.2 <- sqrt(mse.2)
  
  learnCurve$cvRMSE[i] <- rmse.2
}

data <- melt(learnCurve, id.vars = "m", measure.vars = c("cvRMSE", "trainRMSE"))

ggplot2::ggplot(data = data) + 
  ggplot2::geom_line(aes(x=m, y=value, color=variable, group=variable), cex=1.7) + 
  plot.theme(base_size = 10) +
  theme(#legend.position="none",
        axis.text=element_text(size=10),
        axis.title=element_text(size=13,face="bold"),
        legend.title = element_text(size=15, face="bold"),
        legend.text = element_text(size=13)
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank()
        )+
  geom_hline(yintercept=0,
              colour="gray",
              alpha = 0.5) +
  labs(x="Samples in training set",y="RMSE")+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
  coord_cartesian(xlim = c(11.5,38.5))


# plot learning curves of training set size vs. error measure
# for training set and test set


# ---