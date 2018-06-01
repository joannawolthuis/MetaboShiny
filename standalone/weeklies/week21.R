# --- FUNC DEFINITION ---

do_lasso <- function(){
  models <- glmnet_all_alpha(curr = curr,
                             nvars = ncol(config) + 1, 
                             cl = 0,
                             a = alphas,
                             perc_train = input$lasnet_trainfrac/100)
  
  
  # - - store results in mset - - -
  
  plot.many(models, which_alpha = alphas)
  
  lasnet_tables <<- lapply(models, function(x){
    table = x$model$beta
    keep = which(table[,1] > 0)
    table = data.frame("beta" = table[keep,1],
                       "absbeta" = abs(table[keep,1]),
                       row.names = rownames(table)[keep])
    colnames(table) <- c("beta", "abs_beta")
    # - - -
    table
  })
  
  lasnet_stab_tables <<- lapply(models, function(x){
    table = data.frame("perc_chosen" = c(x$int_cv_feat))
    rownames(table) <- names(x$int_cv_feat)
    colnames(table) <- "perc_chosen"
    # - - -
    table
  })
  
  names(lasnet_tables) <<- alphas
  names(lasnet_stab_tables) <<- alphas
  
  list(features=lasnet_tables, model = models)
}

# --- 30 LASSO RUNS TO FIND COMMON CHOSEN FEATURES ---

curr <- as.data.table(mSet$dataSet$preproc)
curr[,(1:ncol(curr)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]

config <- mSet$dataSet$batches[match(rownames(mSet$dataSet$preproc),mSet$dataSet$batches$Sample),]
config <- config[!is.na(config$Sample),]
keep_curr <- match(mSet$dataSet$batches$Sample,rownames(mSet$dataSet$preproc))

config <- cbind(config, Label=mSet$dataSet$cls)

curr <- cbind(config, curr[keep_curr])

input = list(lasnet_alpha = 1, lasnet_trainfrac = 60)

alphas <- input$lasnet_alpha 

maxruns = 50

cl = parallel::makeCluster(3, "FORK")
repeats <- pbapply::pblapply(1:maxruns, cl=cl, function(i){
do_lasso()
})

feat_count <- sapply(repeats, function(x) rownames(x$features[[1]]))
feat_count_tab <- table(unlist(feat_count))
feat_count_dt <- data.table::data.table(feat_count_tab)
colnames(feat_count_dt) <- c("mz", "count")

data <- feat_count_dt[order(feat_count_dt$count, decreasing = T)]
data$mz <- round(as.numeric(data$mz), digits = 2)
data$mz <- factor(data$mz, levels=data$mz)

LASSO_30 <- data$mz[1:20]

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- ggplot(data[1:20,], aes(mz,count))

p + geom_bar(stat = "identity", aes(fill = count)) +
  geom_hline(aes(yintercept=maxruns)) + 
  scale_fill_gradientn(colours = rainbow(10))


library(ROCR)
data(ROCR.xval)
# plot ROC curves for several cross-validation runs (dotted
# in grey), overlaid by the vertical average curve and boxplots
# showing the vertical spread around the average.

pred <- prediction(ROCR.xval$predictions, ROCR.xval$labels)
perf <- performance(pred,"tpr","fpr")
plot(perf,col="grey82",lty=3)
plot(perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE)

xvals <- list()
lapply(repeats, function(x) x$model[[1]]$prediction)
xvals$predictions <- lapply(repeats, function(x) x$model[[1]]$prediction)
xvals$labels <- lapply(repeats, function(x) x$model[[1]]$label)

pred <- prediction(xvals$predictions, xvals$labels)
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="grey82",lty=3, cex.lab=1.3)
plot(perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE)


# --- LASSO: TRAIN ON SPAIN, TEST ON BRASIL ---

curr <- as.data.table(mSet$dataSet$preproc)

br <- grep(config$Sample, pattern = "^BR")
sp <- grep(config$Sample, pattern = "^BR",invert = T)

curr[,(1:ncol(curr)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]

config <- mSet$dataSet$batches[match(rownames(mSet$dataSet$preproc),mSet$dataSet$batches$Sample),]
config <- config[!is.na(config$Sample),]
keep_curr <- match(mSet$dataSet$batches$Sample,rownames(mSet$dataSet$preproc))

config <- cbind(config, Label=mSet$dataSet$cls)

curr <- cbind(config, curr[keep_curr])

# - - - - - - - - - - - - - - - - - - - - - - - - - 

input = list(lasnet_alpha = 1, lasnet_trainfrac = 60)

# SET PARAMS
curr = curr
nvars = ncol(config) + 1
cl = 0
a =  input$lasnet_alpha
perc_train = input$lasnet_trainfrac/100
predictor = "Label"
l = "min"

repeats_spbr <- pbapply::pblapply(1:50, cl=cl, function(i, country=sp){
  # - - LASSO START - -
  
  for (i in (1:nvars)) {
    print(i)
    curr[,i] <- as.numeric(as.factor(curr[,..i][[1]]))-1
  }
  
  # ----------- remove qcs ----------
  
  curr <- curr[which(!grepl(curr$Sample, pattern = "QC"))]
  
  curr[,predictor] <- as.factor(c(curr[,..predictor])[[1]])
  
  reTrain <- caret::createDataPartition(y = config[country,Label],p = 0.6)
  inTrain <- country[reTrain$Resample1]
  
  print(inTrain)
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
    tab = data.frame("beta" = tab[keep,1],
                     "absbeta" = abs(tab[keep,1]),
                     row.names = rownames(tab)[keep])
    colnames(tab) <- c("beta", "abs_beta")
    tab <- tab[sort(tab[,2],decreasing = TRUE)]
    # - - -
    tab
  })
  
  names(lasnet_tables) <- alphas
  
  # - - - -
  list(table = lasnet_tables, model = models)
})

repeats_brsp <- pbapply::pblapply(1:50, cl=cl, function(i, country=br){
  # - - LASSO START - -
  
  for (i in (1:nvars)) {
    print(i)
    curr[,i] <- as.numeric(as.factor(curr[,..i][[1]]))-1
  }
  
  # ----------- remove qcs ----------
  
  curr <- curr[which(!grepl(curr$Sample, pattern = "QC"))]
  
  curr[,predictor] <- as.factor(c(curr[,..predictor])[[1]])
  
  reTrain <- caret::createDataPartition(y = config[country,Label],p = 0.6)
  inTrain <- country[reTrain$Resample1]
  print(inTrain)
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
    tab = data.frame("beta" = tab[keep,1],
                     "absbeta" = abs(tab[keep,1]),
                     row.names = rownames(tab)[keep])
    colnames(tab) <- c("beta", "abs_beta")
    tab <- tab[sort(tab[,2],decreasing = TRUE)]
    # - - -
    tab
  })
  
  names(lasnet_tables) <- alphas
  
  # - - - -
  list(table = lasnet_tables, model = models)
})

xvals <- list()
xvals$predictions <- lapply(repeats_brsp, function(x) x$model[[1]]$prediction)
xvals$labels <- lapply(repeats_brsp, function(x) x$model[[1]]$label)

pred <- prediction(xvals$predictions, xvals$labels)
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="grey82",lty=3, cex.lab=1.3)
plot(perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE)

# - - - - -
feat_count <- lapply(repeats_spbr, function(x){
  beta <- x$model[[1]]$model$beta
  feats <- which(beta[,1] > 0)
  names(feats)
})

feat_count_tab <- table(unlist(feat_count))
feat_count_dt <- data.table::data.table(feat_count_tab)
colnames(feat_count_dt) <- c("mz", "count")

data <- feat_count_dt[order(feat_count_dt$count, decreasing = T)]
data$mz <- round(as.numeric(data$mz), digits = 2)
data$mz <- factor(data$mz, levels=data$mz)

LASSO_12 <- data$mz[1:10]

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- ggplot(data[1:20,], aes(mz,count))

p + geom_bar(stat = "identity", aes(fill = count)) +
  geom_hline(aes(yintercept=maxruns)) + 
  scale_fill_gradientn(colours = rainbow(10))

# -- venn diagram ---

LS_ALL <- unlist(lapply(repeats, function(x) rownames(x$features[[1]])))
LS_12 <- unlist(lapply(repeats_spbr, function(x) rownames(x$table[[1]])))
LS_21 <- unlist(lapply(repeats_brsp, function(x) rownames(x$table[[1]])))

considered_all <- Reduce(intersect, list(LS_ALL,LS_12,LS_21))

ggplotSummary(cpd = considered_all[[6]])

matchies <- lapply(considered_all, multimatch,
           dbs =list.files("/Users/jwolthuis/Google Drive/MetaboShiny/backend/db",
                            full.names = T,
                            pattern = "\\.full\\.db$"), 
           searchid="mz")

DT::datatable(matchies[[6]][,-c("structure", "Description")])
VennDiagram::venn.diagram(
  x = list(LASSO_ALL, LASSO_12, LASSO_21),
  category.names = c("LS_ALL", "LS_12", "LS_21"),
  filename = '#18_venn_diagramm.png',
  output = TRUE ,
  imagetype="png" ,
  height = 980 ,
  width = 980 ,
  resolution = 300,
  # compression = "lzw",
  lwd = 0.5,
  lty = 'blank',
  fill = c('yellow', 'red', 'blue'),
  cex = 0.6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  # cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 135, 140),
  # cat.dist = c(0.055, 0.055, 0.085, 0.055),
  cat.fontfamily = "sans"
  # rotation = 1
)
# - - - save results - - -

LASSO_BRSP <- rownames(lasnet_tables[[1]])

LASSO_SPBR <- rownames(lasnet_tables[[1]])


# --- LASSO: TRAIN ON BRASIL, TEST ON SPAIN (just change above inTrain to 'sp')---
...

# === RF: 1 ===

curr <- as.data.table(mSet$dataSet$preproc)

curr[,(1:ncol(curr)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]

config <- mSet$dataSet$batches[match(rownames(mSet$dataSet$preproc),mSet$dataSet$batches$Sample),]
config <- config[!is.na(config$Sample),]
keep_curr <- match(mSet$dataSet$batches$Sample,rownames(mSet$dataSet$preproc))

config <- cbind(config, Label=mSet$dataSet$cls)

curr <- cbind(config, curr[keep_curr])

input = list(trees = 5000, rf_train_perc = 60)

nvars = ncol(config) + 1
cl = 0

predictor = "Label"

# -- RF: PER COUNTRY --

br <- grep(config$Sample, pattern = "^BR")
sp <- grep(config$Sample, pattern = "^BR",invert = T)

# change which country here!!
inTrain <- br

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

rf_brsp <- importance
rf_spbr <- importance

rfg_sorted <- rf_general[order(rf_general$MeanDecreaseAccuracy,decreasing = T)]
rfbrsp_sorted <- rf_brsp[order(rf_brsp$MeanDecreaseAccuracy,decreasing = T)]
rfspbr_sorted <- rf_spbr[order(rf_spbr$MeanDecreaseAccuracy,decreasing = T)]

RF_1 <- rfg_sorted[1:20,]$rn
RF_BS <- rfbrsp_sorted[1:20,]$rn
RF_SB <- rfspbr_sorted[1:20,]$rn
FOLD_CHANGE = rownames(mSet$analSet$fc$sig.mat)
T_TEST = rownames(mSet$analSet$tt$sig.mat[order(mSet$analSet$tt$sig.mat[,"p.value"], decreasing = FALSE),][1:20,])
PLSDA = plsda_tab[1:20,]$rn


VennDiagram::venn.diagram(
  x = list(FOLD_CHANGE, T_TEST , LASSO_30, LASSO_BRSP, LASSO_SPBR),
  category.names = c("FC" , "TT" , "LA", "L_BS", "L_SB"),
  filename = 'lasso_venn_diagramm.png',
  output = TRUE ,
  imagetype="png" ,
  height = 980 ,
  width = 980 ,
  resolution = 300,
  lwd = 0.5,
  lty = 'blank',
  fill = c('yellow', 'pink', 'purple', 'blue', 'green'),
  cex = 0.6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.fontfamily = "sans"
)


VennDiagram::venn.diagram(
  x = list(FOLD_CHANGE, T_TEST , LASSO_30, RF_1, PLSDA),
  category.names = c("FC" , "TT" , "LASSO", "RF_ALL", "PLSDA"),
  filename = 'las_rf_venn_diagramm.png',
  output = TRUE ,
  imagetype="png" ,
  height = 980 ,
  width = 980 ,
  resolution = 300,
  # compression = "lzw",
  lwd = 0.5,
  lty = 'blank',
  fill = c('yellow', 'pink', 'purple', 'blue', 'green'),
  cex = 0.6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  # cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 135, 140),
  # cat.dist = c(0.055, 0.055, 0.085, 0.055),
  cat.fontfamily = "sans"
  # rotation = 1
)


VennDiagram::venn.diagram(
  x = list(RF_1, RF_BS, RF_SB),
  category.names = c("ALL", "RF_21", "RF_12"),
  filename = '#16_venn_diagramm.png',
  output = TRUE ,
  imagetype="png" ,
  height = 980 ,
  width = 980 ,
  resolution = 300,
  # compression = "lzw",
  lwd = 0.5,
  lty = 'blank',
  fill = c('yellow', 'red', 'blue'),
  cex = 0.6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  # cat.default.pos = "outer",
  cat.pos = c(0, 0, 0),
  # cat.dist = c(0.055, 0.055, 0.085, 0.055),
  cat.fontfamily = "sans"
  # rotation = 1
)

VennDiagram::venn.diagram(
  x = list(RF_1, RF_BS, RF_SB, T_TEST, FOLD_CHANGE),
  category.names = c("RF_ALL", "RF_BS", "RF_SB", "TT", "FC"),
  filename = '#11_venn_diagramm.png',
  output = TRUE ,
  imagetype="png" ,
  height = 980 ,
  width = 980 ,
  resolution = 300,
  # compression = "lzw",
  lwd = 0.5,
  lty = 'blank',
  fill = c('yellow', 'red', 'purple', 'blue', 'green'),
  cex = 0.6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  # cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 135, 140),
  # cat.dist = c(0.055, 0.055, 0.085, 0.055),
  cat.fontfamily = "sans"
  # rotation = 1
)
