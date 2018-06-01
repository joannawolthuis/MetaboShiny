### GOAL A: FILTER VARIABLES MORE EFFICIENTLY

# BASED ON VARIANCE

m <- mSet$dataSet$norm

#BiocInstaller::biocLite("genefilter")

m[is.na(m)] <- 0

filtered <- t(genefilter::varFilter(t(m)))
rownames(filtered) <- rownames(m)

curr <- as.data.table(filtered)

config <- mSet$dataSet$batches[match(rownames(filtered),mSet$dataSet$batches$Sample),]
config <- config[!is.na(config$Sample),]
keep_curr <- match(mSet$dataSet$batches$Sample,rownames(filtered))

config <- cbind(config, Label=mSet$dataSet$cls)

curr <- cbind(config, curr[keep_curr])

nvars = ncol(config) + 1
cl = 0
predictor = "Label"

# BASED ON PERCENTAGE MISSING
# COMBINE?

### GOAL B: RE-RUN LASSO AND RANDOM FOREST

# RANDOM FOREST

input = list(trees = 500, rf_train_perc = 60)


# -- RF: PER COUNTRY --

br <- grep(config$Sample, pattern = "^BR")
sp <- grep(config$Sample, pattern = "^BR",invert = T)

# change which country here!!
inTrain <- sp

# -- RF: NORMAL (PERC PARTITION) --

inTrain <- caret::createDataPartition(y = config$Label,
                                      ## the outcome data are needed
                                      p = 60/100,
                                      ## The percentage of data in the training set
                                      list = FALSE)
# === RF: RUN ===

trainY <- curr[inTrain, 
               ..predictor][[1]]

training <- curr[ inTrain,]
testing <- curr[-inTrain, ]

trainX <- apply(training[, -c(1:ncol(config)), with=FALSE], 2, as.numeric)
testX <- apply(testing[, -c(1:ncol(config)), with=FALSE], 2, as.numeric)

testY <- testing[,predictor, 
                 with=FALSE][[1]]

model = randomForest::randomForest(trainX, trainY, ntree=input$trees,importance=TRUE)
#result <- randomForest::rfcv(trainX, trainY, cv.fold=10, recursive=TRUE)    
#with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

importance = as.data.table(model$importance, keep.rownames = T)
rf_tab <- importance[which(MeanDecreaseAccuracy > 0), c("rn", "MeanDecreaseAccuracy")]
rf_tab <- rf_tab[order(MeanDecreaseAccuracy, decreasing = T)]
rf_tab <<- data.frame(MDA = rf_tab$MeanDecreaseAccuracy, row.names = rf_tab$rn) 

prediction <- stats::predict(model,
                             testX, "prob")[,2]

data <- data.frame(D = as.numeric(as.factor(testY))-1,
                   D.str = testY)
data <- cbind(data, prediction)

roc_coord <- data.frame(D = rep(data[, "D"], length(3)), M = data[, 3], name = rep(names(data)[3], each = nrow(data)), stringsAsFactors = FALSE)
#roc_coord <- plotROC::melt_roc(data, "D", m = 3:ncol(data))

ggplot2::ggplot(roc_coord, 
                ggplot2::aes(d = D, m = M, color = name)) + 
  plotROC::geom_roc(labelsize=0,show.legend = TRUE) + 
  plotROC::style_roc() + 
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=19,face="bold"),
        legend.title=element_text(size=19),
        legend.text=element_text(size=19))

rf_general <- importance

rfg_sorted <- rf_general[order(rf_general$MeanDecreaseAccuracy,decreasing = T)]

ggplotSummary(cpd = rfg_sorted$rn[[6]])

### LASSO

input = list(lasnet_alpha = 1, lasnet_trainfrac = 60)

# SET PARAMS

a =  input$lasnet_alpha
perc_train = input$lasnet_trainfrac/100
predictor = "Label"
l = "min"


repeats <- pbapply::pblapply(1:30, cl=cl, function(i){
  
  # --- LASSO START ---
  
  for (i in (1:nvars)) {
    print(i)
    curr[,i] <- as.numeric(as.factor(curr[,..i][[1]]))-1
  }
  
  # --- remove qcs ---
  
  curr <- curr[which(!grepl(curr$Sample, pattern = "QC"))]
  
  curr[,predictor] <- as.factor(c(curr[,..predictor])[[1]])
  
  # --- SET COUNTRY HERE ---
  
  inTrain <- sp
  
  # ------------------------
  
  inTrain <- caret::createDataPartition(y = config$Label,
                                        ## the outcome data are needed
                                        p = input$lasnet_trainfrac/100,
                                        ## The percentage of data in the training set
                                        list = FALSE)
  # ------------------------
  
  trainY <- curr[inTrain, 
                 ..predictor][[1]]
  
  family <- if(length(levels(as.factor(trainY))) > 2){"multinomial"} else("binomial")
  print(family)
  
  training <- curr[ inTrain,]
  testing <- curr[-inTrain, ]
  
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
    
    var_stability <- table(unlist(used_vars_per_lambda)) / length(used_vars_per_lambda) * 100
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
    
    list(alpha = i, lambda = cv2[[lambda]], model = md, prediction = prediction, labels = testY, int_cv_feat = var_stability)
  })
  
  plot.many(models, which_alpha = alphas)
  
  lasnet_tables <- lapply(models, function(x){
    tab = x$model$beta
    keep = which(tab[,1] > 0)
    tab_new = data.frame("beta" = tab[keep,1],
                         "absbeta" = abs(tab[keep,1]),
                         row.names = rownames(tab)[keep])
    colnames(tab_new) <- c("beta", "abs_beta")
    tab_final <- tab_new[order(tab_new[,1],decreasing = TRUE),]
    # - - -
    tab_final
  })
  
  names(lasnet_tables) <- alphas
  
  lasnet_tables
  
  })

feat_count <- sapply(repeats, function(x) rownames(x[[1]]))
feat_count_tab <- table(unlist(feat_count))
feat_count_dt <- data.table::data.table(feat_count_tab)
colnames(feat_count_dt) <- c("mz", "count")

feat_count_dt_sorted <- feat_count_dt[order(feat_count_dt$count, decreasing = T)]
feat_count_dt_sorted$mz <- factor(feat_count_dt_sorted$mz, levels=feat_count_dt_sorted$mz)


ggplotSummary(cpd = feat_count_dt_sorted$mz[[1]])


### GOAL C: T-SNE ON FILTERED DATASET

### GOAL D: SPARSE GROUP LASSO

library(grpregOverlap)

curr <- as.data.table(mSet$dataSet$preproc)

br <- grep(config$Sample, pattern = "^BR")
sp <- grep(config$Sample, pattern = "^BR",invert = T)

curr[,(1:ncol(curr)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]

config <- mSet$dataSet$batches[match(rownames(mSet$dataSet$preproc),mSet$dataSet$batches$Sample),]
config <- config[!is.na(config$Sample),]
keep_curr <- match(mSet$dataSet$batches$Sample,rownames(mSet$dataSet$preproc))

config <- cbind(config, Label=mSet$dataSet$cls)

curr <- cbind(config, curr[keep_curr])

test <- get_lasso_groups(patdb)
test <- get_all_matches(RSQLite::dbConnect(RSQLite::SQLite(), patdb),
                        which_dbs = list.files("/Users/jwolthuis/Google Drive/MetaboShiny/backend/db",
                                               full.names = T,
                                               pattern = "\\.full\\.db$"),
                        group_by = "mz")
