#' @export
set_shiny_plot_height = function(session, output_width_name, isSquare){
  
  width = session$clientData[[output_width_name]] 
  
  # Do something with the width
  width/if(isSquare) 1 else 2
}


calcHeatMap <- function(mSet, signif.only, 
                        source.anal, 
                        top.hits, 
                        cols){
  sigvals = NULL
  
  #similarly to venn diagram
  flattened <- getTopHits(mSet, 
                          source.anal,
                          top.hits)
  
  sigvals = flattened[[1]]
  
  # change top hits used in heatmap depending on time series / bivariate / multivariate mode
  # reordering of hits according to most significant at the top
  if(!is.null(sigvals)){
    # reorder matrix used
    
    x <- mSet$dataSet$norm[,sigvals]
    final_matrix <- t(x) # transpose so samples are in columns
    
    # check if the sample order is correct - mSet$..$ norm needs to match the matrix
    sample_order <- match(colnames(final_matrix), rownames(mSet$dataSet$norm))
    
    if(mSet$dataSet$exp.type %in% c("2f", "t1f", "t")){
      
      # create convenient table with the ncessary info
      translator <- data.table::data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],GroupA=mSet$dataSet$facA[sample_order], GroupB=mSet$dataSet$facB[sample_order])
      hmap.lvls <- c(levels(mSet$dataSet$facA), levels(mSet$dataSet$facB))
      
      # reorder first by time, then by sample
      split.translator <- split(translator, by = c("GroupB"))
      split.translator.ordered <- lapply(split.translator, function(tbl) tbl[order(tbl$GroupA)])
      translator <- data.table::rbindlist(split.translator.ordered)
      
      # ensure correct sample order
      final_matrix <- final_matrix[,match(translator$Sample, colnames(final_matrix))]
      
      # disable automatic ordering of samples through clustering
      my_order=F
      
    }else{
      # no complicated reordering necessary
      translator <- data.table::data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],
                                           Group=mSet$dataSet$cls[sample_order])
      hmap.lvls <- levels(mSet$dataSet$cls)
      my_order = T # enable sorting through dendrogram
    }
    
    # create name - to - color mapping vector for the plotting functions
    color.mapper <- {
      classes <- hmap.lvls
      cols <- sapply(1:length(classes), function(i) cols[i]) # use user-defined colours
      names(cols) <- classes
      # - - -
      cols
    }
    
    # write to mSet
    mSet$analSet$heatmap <- list(
      matrix = final_matrix,
      translator = translator,
      colors = color.mapper,
      my_order = my_order)
    
    # - - - - -
    
    mSet
  }
}

runML <- function(curr,
                  config,
                  train_vec,
                  test_vec,
                  configCols,
                  ml_method,
                  ml_perf_metr,
                  ml_folds,
                  ml_preproc,
                  tuneGrid,
                  ml_train_perc,
                  sampling = "none"){
  
  # get user training percentage
  ml_train_perc <- ml_train_perc/100
  
  if(unique(train_vec)[1] == "all" & unique(test_vec)[1] == "all"){ # BOTH ARE NOT DEFINED
    test_idx = caret::createDataPartition(y = curr$label, p = ml_train_perc, list = FALSE) # partition data in a balanced way (uses labels)
    train_idx = setdiff(1:nrow(curr), test_idx) #use the other rows for testing
    inTrain = train_idx
    inTest = test_idx
  }else if(unique(train_vec)[1] != "all"){ #ONLY TRAIN IS DEFINED
    train_idx <- which(config[,train_vec[1], with=F][[1]] == train_vec[2])
    test_idx = setdiff(1:nrow(curr), train_idx) # use the other rows for testing
    reTrain <- caret::createDataPartition(y = curr[train_idx, label], p = ml_train_perc) # take a user-defined percentage of the regexed training set
    inTrain <- train_idx[reTrain$Resample1]
    inTest = test_idx
  }else{ # ONLY TEST IS DEFINED
    test_idx = which(config[,test_vec[1], with=F][[1]] == test_vec[2])
    train_idx = setdiff(1:nrow(curr), test_idx) # use the other rows for testing
    reTrain <- caret::createDataPartition(y = curr[train_idx, label], p = ml_train_perc) # take a user-defined percentage of the regexed training set
    inTrain <- train_idx[reTrain$Resample1]
    inTest <- test_idx
  }
  
  if("split" %in% colnames(curr)){
    curr <- curr[,-"split"]
  }
  
  # choose predictor "label" (some others are also included but cross validation will be done on this)
  predictor = "label"
  
  # split training and testing data
  trainY <- curr[inTrain,
                 ..predictor][[1]]
  testY <- curr[inTest,
                ..predictor][[1]]
  
  training <- curr[inTrain,]
  testing <- curr[inTest,]
  
  require(caret)
  
  if(ml_folds == "LOOCV"){
    trainCtrl <- caret::trainControl(verboseIter = T,
                                     allowParallel = F,
                                     method="LOOCV",
                                     trim=TRUE, 
                                     returnData = FALSE,
                                     sampling = sampling) # need something here...
    
  }else{
    trainCtrl <- caret::trainControl(verboseIter = T,
                                     allowParallel = F,
                                     method = as.character(ml_perf_metr),
                                     number = as.numeric(ml_folds),
                                     repeats=3,
                                     trim=TRUE, 
                                     returnData = FALSE,
                                     sampling = sampling) # need something here...
  }
  
  fit <- caret::train(
    label ~ .,
    data = training,
    method = ml_method,
    ## Center and scale the predictors for the training
    ## set and all future samples.
    preProc = ml_preproc,
    tuneGrid = if(nrow(tuneGrid) > 0) tuneGrid else NULL,
    trControl = trainCtrl
  )
  
  result.predicted.prob <- stats::predict(fit, testing, type="prob") # Prediction
  
  # train and cross validate model
  # return list with mode, prediction on test data etc.s
  list(type = ml_method,
       # model = fit,
       importance = caret::varImp(fit)$importance,
       prediction = result.predicted.prob[,2],
       labels = testing$label)
  
}

getMultiMLperformance <- function(model){
  roc = pROC::multiclass.roc(model$labels, model$prediction)
  data.table::rbindlist(lapply(roc$rocs, function(roc.pair){
    data.table(FPR = sapply(roc.pair$specificities, function(x) 1-x),
               TPR = roc.pair$sensitivities,
               AUC_AVG = as.numeric(roc$auc),
               AUC_PAIR = as.numeric(pROC::auc(roc.pair)),
               comparison = paste0(roc.pair$levels,collapse=" vs. "))
  })) 
}

getColDistribution <- function(csv){
  suppressWarnings({
    gsubbed = gsub(x = colnames(csv),
                   pattern="[+|\\-].*$",
                   replacement="")
    as.numi <- as.numeric(gsubbed)
    exp.vars <- which(is.na(as.numi))
    mz.vars <- which(!is.na(as.numi))  
  })
  list(meta = exp.vars, mz = mz.vars)
}

cleanCSV <- function(csv, regex=" |\\(|\\)|\\+",exp.vars, mz.vars, miss.meta, miss.mz){
  # remove whitespace
  csv$sample <- gsub(csv$sample, pattern=regex, replacement="")
  
  # convert all 0's to NA so metaboanalystR will recognize them
  csv[,(exp.vars) := lapply(.SD,function(x){ ifelse(x == "" | is.na(x) | x == "Unknown", miss.meta, x)}), .SDcols = exp.vars]
  csv[,(mz.vars) := lapply(.SD,function(x){ ifelse(x == 0, miss.mz, x)}), .SDcols = mz.vars]
  csv <- csv[,which(unlist(lapply(csv, function(x)!all(is.na(x))))),with=F]
  colnames(csv)[which(colnames(csv) == "time")] <- "Time"
  csv$sample <- gsub("[^[:alnum:]./_-]", "", csv$sample)
  suppressWarnings({
    csv <- csv[,-"label"]
  })
  # - - - 
  csv
}

getDefaultCondition <- function(csv, exp.vars, excl.rows, excl.cond, min.lev){
  unique.levels <- apply(csv[!excl.rows, ..exp.vars, with=F], MARGIN=2, function(col){
    lvls <- levels(as.factor(col))
    # - - - - - -
    length(lvls)
  })
  unique.levels <- unique.levels[!(names(unique.levels) %in% excl.cond)]
  # use this as the default selected experimental variable (user can change later)
  which.default <- unique.levels[which(unique.levels == min(unique.levels[which(unique.levels > min.lev)]))][1]
  condition = names(which.default)
  condition
}

removeOutliers <- function(csv, exp.vals){
  sums <- rowSums(csv[,-exp.vars,with=FALSE],na.rm = TRUE)
  names(sums) <- csv$sample
  outliers = c(car::Boxplot(as.data.frame(sums)))
  csv[!(sample %in% outliers),]  
}

removeUnusedQC <- function(csv, covar_table){
  samps <- which(!grepl(covar_table$sample, pattern = "QC"))
  batchnum <- unique(covar_table[samps, "batch"][[1]])
  keep_samps_post_qc <- covar_table[which(covar_table$batch %in% batchnum),"sample"][[1]]
  covar_table <- covar_table[which(covar_table$batch %in% batchnum),]
  csv[which(csv$sample %in% keep_samps_post_qc),]  
}

asMetaboAnalyst <- function(csv, exp.vars){
  # remove all except sample and time in saved csv
  exp_var_names <- colnames(csv)[exp.vars]
  keep_cols <-  c("sample", "label")
  remove <- which(!(exp_var_names %in% keep_cols))
  csv[,-remove,with=F]
}

mzLeftPostFilt <- function(mSet, perc_limit){
  int.mat <- mSet$dataSet$preproc
  minConc <- mSet$dataSet$minConc
  missvals = apply(is.na(int.mat), 2, sum)/nrow(int.mat)
  good.inx <- missvals < perc_limit/100
  if(length(which(good.inx))==0){
    MetaboShiny::metshiAlert(paste("No m/z left after filtering, please make your missing value correction more lenient... Recommended minumum to retain at least 1 m/z value:", paste0(min(missvals)*100, "%")))
    return(NULL)
  }  
}

tooEmptySamps <- function(mSet, max.missing.per.samp){
  w.missing <- mSet$dataSet$preproc
  miss.per.samp = rowSums(is.na(w.missing))
  miss.per.samp.perc = sapply(miss.per.samp, function(x)( x / ncol(w.missing) ) * 100)
  which(miss.per.samp.perc >= max.missing.per.samp)
}

replRowMin <- function(mSet){
  samples <- rownames(mSet$dataSet$preproc)
  w.missing <- data.table::as.data.table(mSet$dataSet$preproc)
  pb <- pbapply::startpb(0, ncol(w.missing))
  i=0
  w.missing = w.missing[,as.data.table(apply(.SD, 2, function(x) {i<<-i+1; pbapply::setpb(pb, i); x[is.na(x)] <- mean(x, na.rm=T)/2; x}))]
  w.missing <- as.data.frame(w.missing)
  rownames(w.missing) <- samples
  return(w.missing)
  }

replRF <- function(mSet, parallelMode, ntree, cl){
  samples <- rownames(mSet$dataSet$preproc)
  # convert all to as numeric
  # TODO: remove, should be automatic
  w.missing <- mSet$dataSet$preproc
  w.missing <- apply(w.missing, 2, as.numeric)
  
  # register other threads as parallel threads
  doParallel::registerDoParallel(cl)
  
  # set amount of tries (defined by missforest package)
  auto.mtry <- floor(sqrt(ncol(mSet$dataSet$preproc)))
  
  mtry <- ifelse(auto.mtry > 100, 
                 100, 
                 auto.mtry)
  
  # impute missing values with random forest
  imp <- missForest::missForest(w.missing,
                                parallelize = parallelMode, # parallelize over variables, 'forests' is other option
                                verbose = F,
                                ntree = ntree,
                                mtry = mtry)
  imp$ximp
}

batchCorrQC <- function(mSet){
  smps <- rownames(mSet$dataSet$norm)
  # get which rows are QC samples
  qc_rows <- which(grepl(pattern = "QC", x = smps))
  # get batch for each sample
  batch.idx = as.numeric(as.factor(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$sample),"batch"][[1]]))
  if(length(batch.idx) == 0) return(mSet$dataSet$norm)
  # get injection order for samples
  seq.idx = as.numeric(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$sample),"injection"][[1]])
  # go through all the metabolite columns
  corr_cols <- pbapply::pblapply(1:ncol(mSet$dataSet$norm), function(i){
    # fetch non-corrected values
    vec = mSet$dataSet$norm[,i]
    # correct values using QCs and injectiono rder
    corr_vec = BatchCorrMetabolomics::doBC(Xvec = as.numeric(vec),
                                           ref.idx = as.numeric(qc_rows),
                                           batch.idx = batch.idx,
                                           seq.idx = seq.idx,
                                           result = "correctedX",
                                           minBsamp = 1) # at least one QC necessary
    corr_vec
  })
  
  # cbind the corrected columns to re-make table
  qc_corr_matrix <- as.data.frame(do.call(cbind, corr_cols))
  # fix rownames to old rownames
  colnames(qc_corr_matrix) <- colnames(mSet$dataSet$norm)
  rownames(qc_corr_matrix) <- rownames(mSet$dataSet$norm)
  as.data.frame(qc_corr_matrix)
}

hideQC <- function(mSet){
  smps <- rownames(mSet$dataSet$norm)
  # get which rows are QC samples
  qc_rows <- which(grepl(pattern = "QC", x = smps))
  mSet$dataSet$norm <- mSet$dataSet$norm[-qc_rows,]
  mSet$dataSet$cls <- mSet$dataSet$cls[-qc_rows, drop = TRUE]
  mSet$dataSet$covars <- mSet$dataSet$covars[-grep("QC", mSet$dataSet$covars$sample),]
  mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
  mSet
}

combatCSV <- function(mSet){
  # get sample names and classes
  smp <- rownames(mSet$dataSet$norm)
  exp_lbl <- mSet$dataSet$cls
  
  # create csv for comBat
  csv <- data.table::as.data.table(cbind(sample = smp,
                                         label = mSet$dataSet$cls,
                                         mSet$dataSet$norm))
  
  # transpose for combat
  csv_edata <-t(csv[,!c(1,2)])
  colnames(csv_edata) <- csv$sample
  csv_edata
}

metshiTable <- function(content, options=NULL){
  opts = list(deferRender = TRUE, 
              scrollY = 200,
              searching = TRUE,
              scrollCollapse = FALSE,
              rownames = FALSE,
              scroller = TRUE,
              scrollX = T,
              dom = 'Bftip',
              buttons = 
                list(list(
                  extend = 'collection',
                  buttons = c('csv', 'excel', 'copy'),
                  text = "<i class='fa fa-save'></i>"
                )))
  if(!is.null(options)){
    opts <- append(opts, options)      
  }
  mz_rownames = stringr::str_match(rownames(content), "(\\d+\\.\\d+)")[,2]
  if(!is.na(mz_rownames[1])){
    rownames(content) <- paste0(rownames(content), sapply(rownames(content), function(mz) if(grepl("\\-", mz)) "" else "+"))
  }
  DT::datatable(content,
                selection = 'single',
                autoHideNavigation = T,
                class = 'compact', height = "500px",
                extensions = c("FixedColumns", "Scroller", "Buttons"), 
                options = opts
  )
}

getTopHits <- function(mSet, expnames, top){
  
  experiments <- stringr::str_match(expnames, 
                                    pattern = "\\(.*\\)")[,1]
  
  experiments <- unique(gsub(experiments, pattern = "\\(\\s*(.+)\\s*\\)", replacement="\\1"))
  
  table_list <- lapply(experiments, function(experiment){
    
    analysis = mSet$storage[[experiment]]$analysis
    
    rgx_exp <- gsub(experiment, pattern = "\\(", replacement = "\\\\(")
    rgx_exp <- gsub(rgx_exp, pattern = "\\)", replacement = "\\\\)")
    rgx_exp <- gsub(rgx_exp, pattern = "\\-", replacement = "\\\\-")
    rgx_exp <- gsub(rgx_exp, pattern = "\\+", replacement = "\\\\+")
    
    categories = grep(unlist(expnames),
                      pattern = paste0("\\(",rgx_exp, "\\)"), value = T)
    
    categories = gsub(categories, pattern = " \\(\\s*(.+)\\s*\\)", replacement = "")
    
    # go through the to include analyses
    
    tables <- lapply(categories, function(name){
      
      base_name <- search_name <- gsub(name, pattern = " -.*$| ", replacement="")
      
      if(base_name %in% gbl$constants$ml.models){
        search_name <- "ml"
      }
      
      # fetch involved mz values
      tbls <- switch(search_name,
                     ml = {
                       which.ml <- gsub(name, pattern = "^.*- | ", replacement="")
                       mzvals = analysis$ml[[base_name]][[which.ml]]$bar[order(analysis$ml[[base_name]][[which.ml]]$bar$importance,
                                                                               decreasing = T),]$mz
                       mzvals <- type.convert(gsub(mzvals, pattern = "'|`", replacement=""))
                       res <- list(mzvals)
                       names(res) <- paste0(which.ml, " (", base_name, ")")
                       # - - -
                       res
                     },
                     aov = {
                       res = list(rownames(analysis$aov$sig.mat[order(analysis$aov$sig.mat[,2],
                                                                                 decreasing = F),]))
                       names(res) = base_name
                       res
                     },
                     aov2 = {
                       res = list(rownames(analysis$aov2$sig.mat[order(analysis$aov2$sig.mat[,"Interaction(adj.p)"],
                                                                                  decreasing = F),]))
                       names(res) = base_name
                       res
                     },
                     asca = {
                       res = list(rownames(analysis$asca$sig.list$Model.ab[order(analysis$asca$sig.list$Model.ab[,1],
                                                                                            decreasing = T),]))
                       names(res) = base_name
                       res
                     },
                     MB = {
                       res = list(rownames(analysis$MB$stats)[order(analysis$MB$stats[,1],
                                                                                decreasing = T)])
                       names(res) = base_name
                       res
                     },
                     tt = {
                       res = list(rownames(analysis$tt$sig.mat[order(analysis$tt$sig.mat[,2],
                                                                                decreasing = F),]))
                       names(res) = base_name
                       res
                     },
                     fc = {
                       res = list(rownames(analysis$fc$sig.mat[order(abs(analysis$fc$sig.mat[,2]),
                                                                                decreasing = F),]))
                       names(res) = base_name
                       res
                     },
                     volcano = {
                       res = list(rownames(analysis$volcano$sig.mat))
                       names(res) = base_name
                       res
                     },
                     plsda = {
                       which.plsda <- gsub(name, pattern = "^.*- | ", replacement="")
                       
                       compounds_pc <- data.table::as.data.table(analysis$plsda$vip.mat,keep.rownames = T)
                       colnames(compounds_pc) <- c("rn", paste0("PC", 1:(ncol(compounds_pc)-1)))
                       ordered_pc <- setorderv(compounds_pc, which.plsda, -1)
                       
                       res <- list(ordered_pc$rn)
                       names(res) <- paste0(which.plsda, " (PLS-DA)")
                       # - - -
                       res
                     },
                     pca = {
                       which.pca <- gsub(name, pattern = "^.*- | ", replacement="")
                       
                       compounds_pc <- data.table::as.data.table(analysis$pca$rotation,keep.rownames = T)
                       ordered_pc <- setorderv(compounds_pc, which.pca, -1)
                       res <- list(ordered_pc$rn)
                       names(res) <- paste0(which.pca, " (PCA)")
                       # - - -
                       res
                     },
                     volc = {
                       res <- list(rownames(analysis$volcano$sig.mat))
                       names(res) = base_name
                       res
                     },
                     {MetaboShiny::metshiAlert("Not currently supported...")
                       return(NULL)})
      
      if(is.null(tbls)) return(NULL)
      
      # user specified top hits only
      tbls_top <- lapply(tbls, function(tbl){
        if(length(tbl) < top){
          tbl
        }else{
          tbl[1:top]
        }
      })
      names(tbls_top) <- paste0(experiment, ": ", names(tbls_top))
      tbls_top
    })
    
    # unnest the nested lists
    flattened <- flattenlist(tables)
    
    # remove NAs
    flattened <- lapply(flattened, function(x) x[!is.na(x)])
    
    #rename and remove regex-y names
    names(flattened) <- gsub(x = names(flattened), pattern = "(.*\\.)(.*$)", replacement = "\\2")
    # return
    flattened
  })
  
  flattened <- flattenlist(table_list)
  names(flattened) <- gsub(x = names(flattened), pattern = "(.*\\.)(.*$)", replacement = "\\2")
  flattened <- lapply(flattened, function(x) x[!is.na(x)])
  return(flattened)
}
