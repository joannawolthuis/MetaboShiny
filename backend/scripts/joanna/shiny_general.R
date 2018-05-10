#' @export
get_ref_vars <- function(fac="Label"){
  req(csv_loc)
  # csv_loc = "backend/appdata/euronutrition/Test.csv"
  csv <- fread(csv_loc, sep="\t", header = T)[,1:5]
  # --- return ---
  c(unique(csv[,..fac]))[[fac]]
}

get_ref_cpds <- function(){
  req(csv_loc)
  # csv_loc = "backend/appdata/euronutrition/Test.csv"
  csv <- fread(csv_loc, sep="\t", header = T)[1,]
  keep.colnames <- colnames(csv)[!colnames(csv) %in% c("Sample", "Label", "Time")]
  # --- return ---
  c(keep.colnames)
}

getProfile <- function(varName, title=varName, mode="stat"){
  # ---------------
  print(varName)
  print(colnames(mSet$dataSet$norm))
  # ---------------
  varInx <- colnames(mSet$dataSet$norm) == varName;
  var <- as.data.table(mSet$dataSet$norm, 
                       keep.rownames = T)[,varInx, with=FALSE];
  samp.names <- rownames(mSet$dataSet$norm)
  # ---------------
  if(mode == "time"){
    time.fac <<- mSet$dataSet$time.fac;
    translator <- data.table(
      index = 1:length(samp.names),
      Sample = gsub(x = samp.names, pattern = "_T\\d$", replacement=""),
      Group = mSet$dataSet$facB,
      Time = time.fac,
      Abundance = mSet$dataSet$norm[,varInx]
    )
  }else if(mode == "stat"){
    translator <- data.table(
      index = 1:length(samp.names),
      Sample = gsub(x = samp.names, pattern = "T\\d$", replacement=""),
      Group = mSet$dataSet$cls,
      Abundance = mSet$dataSet$norm[,varInx]
    )
  }
  # ---------------
  return(translator)
}

kegg.charge <- function(atomlist){
  charges  <-regmatches(
    atomlist,
    regexpr(atomlist, pattern = "#[+-]|#\\d*[+-]",perl = T)
  ) 
  formal_charge = 0
  for(ch in charges[!is.na(charges)]){
    ch.base <- gsub(ch, pattern = "#", replacement = "")
    ch.base <- if(ch.base == "-" | ch.base == "+") paste0(ch.base, 1) else(ch.base)
    ch.base <- gsub(ch.base, pattern = "\\+", replacement = "")
    ch.base <- as.numeric(sub("([0-9.]+)-$", "-\\1", ch.base))
    # -------
    formal_charge = formal_charge + ch.base
  }
  formal_charge
}

glmnet_all_alpha <- function(curr, nvars, cl = 0, a = c(0, 0.2, 0.4, 0.6, 0.8, 1), predictor = "Label", l = "min", perc_train = 0.6){
  
  # first 4 are identifiers, leave alone 

  curr_orig <- curr
  
  curr <- curr_orig
  
  # mSetObj <- mSet
  
  # - - - - - - - - - - -
  
  for (i in (1:nvars)) {
    print(i)
    curr[,i] <- as.numeric(as.factor(curr[,..i][[1]]))-1
  }
  
  # ----------- remove qcs ----------
  
  curr <- curr[which(!grepl(curr$Sample, pattern = "QC"))]
  
  # ---------------------------------
  
  curr[,predictor] <- as.factor(c(curr[,..predictor])[[1]])
  
  print(perc_train)
  
  inTrain <- caret::createDataPartition(y = c(curr[,..predictor])[[1]],
                                 ## the outcome data are needed
                                 p = perc_train,
                                 ## The percentage of data in the training set
                                 list = FALSE)
  trainY <- curr[inTrain, 
                 ..predictor][[1]]
  
  family <- if(length(levels(as.factor(trainY))) > 2){"multinomial"} else("binomial")
  print(family)
  
  training <- curr[ inTrain,]
  testing <- curr[-inTrain, ]
  
  print(dim(testing))
  print(dim(training))
  
  predIdx <- which(colnames(curr) == predictor | colnames(curr) == "Sample")
  
  trainX <- apply(training[, -c(predIdx), with=FALSE], 2, as.numeric)
  testX <- apply(testing[, -c(predIdx), with=FALSE], 2, as.numeric)
  
  testY <- testing[,predictor, 
                   with=FALSE][[1]]

  models <- pbapply::pblapply(a,  cl=cl, FUN=function(i){
    
    print(i)
    
    lambda = paste0('lambda.', l)
    
    nfold = length(trainY)
    
    cv1 <- glmnet::cv.glmnet(trainX, trainY, family = family, nfold = nfold, type.measure = "auc", alpha = i, keep = TRUE)
    
    used_vars_per_lambda <- lapply(cv1$lambda, function(l){
      c<-coef(cv1, s = l,exact=TRUE)
      inds<-which(c[,1]>0)
      c(names(inds)) 
    })
    
    var_stability <- table(unlist(used_vars_per_lambda)) / 99 * 100
    var_stability <- var_stability[2:length(var_stability)]
    
    # - - - - - - - - - - - - - - - - - - - - - - - - -
    
    cv2 <- data.frame(cvm = cv1$cvm[cv1$lambda == cv1[[lambda]]], lambda = cv1[[lambda]], alpha = i)
    
    # --------
    
    md <- glmnet::glmnet(as.matrix(trainX), trainY, family = family, lambda = cv2$lambda, alpha = cv2$alpha)
    
    # --- roc ---
    
    prediction <- stats::predict(md,
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
    
    list(alpha = i, lambda = cv2[[lambda]], model = md, prediction = prediction, labels = testY, int_cv_feat = var_stability)
  })
  models
}

