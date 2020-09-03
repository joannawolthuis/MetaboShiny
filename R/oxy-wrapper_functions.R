globalVariables(c("-log(p)", "-log10P", "..change_var", "..col.fac", "..count..", "..exp.vars", "..keep.cols", "..matches", "..predictor", "..rmcols", "..scaled..", "..shape.fac", "..stats_var", "..time_var", "..x", "Abundance", "Color", "FPR", "Group", "GroupA", "GroupB", "Individual", "Label", "Metric", "PC", "Peak", "Sample", "Shape", "TPR", "Text", "Value", "abstract", "acc", "adduct", "aes", "attempt", "color", "comparison", "coord_flip", "correlation", "count", "debug_browse_content", "debug_enrich", "debug_input", "debug_lcl", "debug_mSet", "debug_matches", "debug_pieinfo", "debug_result_filters", "debug_selection", "exp.vars", "expand_limits", "extremity", "facet_grid", "facet_wrap", "freq", "fullformula", "gbl", "geom_point", "ggplot", "ggtitle", "group", "importance.mean", "individual", "label", "labs", "lcl", "log2FC", "log2fc", "m/z", "p-value", "pathway", "pc", "position_dodge", "position_jitterdodge", "samples", "scale_size_area", "scale_x_discrete", "scale_y_discrete", "searchRev", "shape", "value", "variable", "x", "xlab", "y"))

#' @title Collect info needed for heatmap
#' @description Given an mSet and an analysis of interest, generates the matrix etc. needed for heatmap creation later on.
#' @param mSet mSet object
#' @param signif.only Only include significant hits?
#' @param source.anal Source analysis (tt, aov, etc.)
#' @param top.hits Top n hits to include in heatmap
#' @param cols Color vector
#' @return List with matrix and other required information for heatmap creation.
#' @seealso 
#'  \code{\link[data.table]{rbindlist}}
#' @rdname calcHeatMap
#' @export 
#' @importFrom data.table data.table rbindlist
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
#' @param ml_train_perc Percentage in training
#' @param sampling Up- or downsampling settings, Default: 'none'
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

#' @title Get performance for multi-comparison ML model
#' @description ROC curves can be a bit tricky for multivariate models. This evaluates each possible pair of categories to generate individual and average AUC.
#' @param model ML model
#' @return FPR,TPR,average AUC,AUC for a given pair, and the name of the comparison
#' @seealso 
#'  \code{\link[pROC]{multiclass.roc}},\code{\link[pROC]{auc}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname getMultiMLperformance
#' @export 
#' @importFrom pROC multiclass.roc auc
#' @importFrom data.table rbindlist
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

#' @title Get metadata/mz column distribution
#' @description Which columns are metadata, which are m/z values?
#' @param csv CSV to evaluate
#' @return List with meta = which columns (in numbers) are metadata, mz = which are mz.
#' @rdname getColDistribution
#' @export 
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

#' @title Clean csv
#' @description Remove whitespace, replace zeros with user preferred missing value filler, cleans sample names with regex.
#' @param csv CSV to clean
#' @param regex Regex to apply to sample names., Default: ' |\(|\)|\+'
#' @param exp.vars Column numbers that are metadata
#' @param mz.vars Column numbers that are m/z
#' @param miss.meta What to fill missing values with in metadata
#' @param miss.mz What to fill missing values with in m/z values
#' @return Data table
#' @rdname cleanCSV
#' @export 
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

#' @title Get default experimental condition
#' @description Searches for experimental condition that has at least 'min.lev' categories, and then the least amount of categories within that.
#' @param csv CSV to search in
#' @param exp.vars Which columns are metadata?
#' @param excl.rows Rows to exclude from evaluation
#' @param excl.cond Conditions to exclude from evaluation
#' @param min.lev Minimum kevels a condition needs to have to be considered. 
#' @return Metadata column name that is chosen as initial experimental variable.
#' @rdname getDefaultCondition
#' @export 
getDefaultCondition <- function(csv, exp.vars, excl.rows, excl.cond, min.lev){
  unique.levels <- apply(csv[!excl.rows, ..exp.vars, with=F], MARGIN=2, function(col){
    lvls <- levels(as.factor(col))
    # - - - - - -
    length(lvls)
  })
  unique.levels <- unique.levels[!(names(unique.levels) %in% excl.cond)]
  min.grpsize = apply(csv[,names(unique.levels), with=F],2,function(x)min(table(x)))
  # use this as the default selected experimental variable (user can change later)
  which.default <- unique.levels[which(unique.levels == min(unique.levels[which(unique.levels > min.lev)]) &
                                       min.grpsize >= 3)][1]
  condition = names(which.default)
  condition
}

#' @title Remove outliers
#' @description Uses boxplot function to exclude outlier samples from csv.
#' @param csv Data table
#' @param exp.vals Which columns are metadata?
#' @return Data table with outliers removed
#' @seealso 
#'  \code{\link[car]{Boxplot}}
#' @rdname removeOutliers
#' @export 
#' @importFrom car Boxplot
removeOutliers <- function(csv, exp.vals){
  sums <- rowSums(csv[,-exp.vars,with=FALSE],na.rm = TRUE)
  names(sums) <- csv$sample
  outliers = c(car::Boxplot(as.data.frame(sums)))
  csv[!(sample %in% outliers),]  
}

#' @title Remove unused QC samples
#' @description Removes QC samples that don't have any matching samples in their batch, thus not useful for batch correction.
#' @param csv Data table
#' @param covar_table Metadata table (which has the batch data)
#' @return Data table without unused QC samples
#' @rdname removeUnusedQC
#' @export 
removeUnusedQC <- function(csv, covar_table){
  samps <- which(!grepl(covar_table$sample, pattern = "QC"))
  batchnum <- unique(covar_table[samps, "batch"][[1]])
  keep_samps_post_qc <- covar_table[which(covar_table$batch %in% batchnum),"sample"][[1]]
  covar_table <- covar_table[which(covar_table$batch %in% batchnum),]
  csv[which(csv$sample %in% keep_samps_post_qc),]  
}

#' @title Reformat peaktable to be in MetaboAnalyst format
#' @description Takes the MetShi peaktable and reformats it to be accepted by MetaboAnalystR as input.
#' @param csv Data table
#' @param exp.vars Which columns are metadata?
#' @return Data frame that can be imported
#' @rdname asMetaboAnalyst
#' @export 
asMetaboAnalyst <- function(csv, exp.vars){
  # remove all except sample and time in saved csv
  exp_var_names <- colnames(csv)[exp.vars]
  keep_cols <-  c("sample", "label")
  remove <- which(!(exp_var_names %in% keep_cols))
  csv[,-remove,with=F]
}

#' @title Check if any m/z are left after missing value based filtering
#' @description After excluding m/z values with more than 'perc_limit' samples missing, are any m/z values left for analysis?
#' @param mSet mSet object
#' @param perc_limit Max allowed missing samples for a m/z value in percent. 
#' @rdname mzLeftPostFilt
#' @export 
mzLeftPostFilt <- function(mSet, perc_limit){
  int.mat <- mSet$dataSet$preproc
  minConc <- mSet$dataSet$minConc
  missvals = apply(is.na(int.mat), 2, sum)/nrow(int.mat)
  good.inx <- missvals < perc_limit/100
  if(length(which(good.inx))==0){
    metshiAlert(paste("No m/z left after filtering, please make your missing value correction more lenient... Recommended minumum to retain at least 1 m/z value:", paste0(min(missvals)*100, "%")))
    return(NULL)
  }  
}

#' @title Which samples have too many m/z missing?
#' @description Checks which samples have more than 'max.missing.per.samp' percent of m/z values missing. This can cause problems for KNN means missing value imputation (if >80% is missing).
#' @param mSet mSet object
#' @param max.missing.per.samp Max percentage of m/z allowed to be missing for a given sample.
#' @return Indices of which samples should be removed according to this threshold.
#' @rdname tooEmptySamps
#' @export 
tooEmptySamps <- function(mSet, max.missing.per.samp){
  w.missing <- mSet$dataSet$preproc
  miss.per.samp = rowSums(is.na(w.missing))
  miss.per.samp.perc = sapply(miss.per.samp, function(x)( x / ncol(w.missing) ) * 100)
  which(miss.per.samp.perc >= max.missing.per.samp)
}

#' @title Fill missing values with half minimum of this sample
#' @description Performs missing value filling with the sample half minimum pre-normalization.
#' @param mSet mSet object
#' @return mSet object
#' @seealso 
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[pbapply]{pboptions}}
#' @rdname replRowMin
#' @export 
#' @importFrom data.table as.data.table
#' @importFrom pbapply startpb setpb
replRowMin <- function(mSet){
  samples <- rownames(mSet$dataSet$preproc)
  mz <- colnames(mSet$dataSet$preproc)
  w.missing <- data.table::as.data.table(mSet$dataSet$preproc)
  pb <- pbapply::startpb(0, nrow(w.missing))
  i=0
  w.missing = w.missing[,data.table::as.data.table(apply(.SD, 
                                                         1, 
                                                         function(x) {i<<-i+1; 
                                                         pbapply::setpb(pb, i); 
                                                         x[is.na(x)] <- min(x, na.rm=T)/2; 
                                                         return(x)
                                                         }))]
  w.missing <- t(as.data.frame(w.missing))
  rownames(w.missing) <- samples
  colnames(w.missing) <- mz
  return(w.missing)
  }

#' @title Impute missing values with Random Forest
#' @description Uses the missForest package to impute missing values. The most time-consuming but accurate imputation function.
#' @param mSet mSet object
#' @param parallelMode Parallelize 'variables' or 'forests'?
#' @param ntree How many trees per m/z value?
#' @param cl parallel::makeCluster object for multithreading
#' @return mSet object
#' @seealso 
#'  \code{\link[doParallel]{registerDoParallel}}
#'  \code{\link[missForest]{missForest}}
#' @rdname replRF
#' @export 
#' @importFrom doParallel registerDoParallel
#' @importFrom missForest missForest
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

#' @title Batch correction using QC samples
#' @description Using QC samples, attempt to correct existing batch effect.
#' @param mSet mSet object
#' @return mSet object
#' @seealso 
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[BatchCorrMetabolomics]{doBC}}
#' @rdname batchCorrQC
#' @export 
#' @importFrom pbapply pblapply
#' @importFrom BatchCorrMetabolomics doBC
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

#' @title Remove QC values from dataset
#' @description Post-batch correction, it's likely that users don't need the QC samples anymore. This removes them from the mSet, including their metadata rows.
#' @param mSet mSet object
#' @return mSet object
#' @rdname hideQC
#' @export 
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

#' @title Prepare COMBAT table
#' @description COMBAT batch correction requires a certain table format. This adjusts the mSet tables to fit that standard.
#' @param mSet mSet object
#' @return CSV ready for use in COMBAT.
#' @seealso 
#'  \code{\link[data.table]{as.data.table}}
#' @rdname combatCSV
#' @export 
#' @importFrom data.table as.data.table
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

#' @title Generate interactive table
#' @description Generates interactive table in MetaboShiny. Adjust this function to change all table display functionality!
#' @param content Table content (generally a data table)
#' @param options Table options (in DT format), Default: NULL
#' @param rownames Use row names?, Default: T
#' @return DT data table object for use in shiny.
#' @seealso 
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[DT]{datatable}}
#' @rdname metshiTable
#' @export 
#' @importFrom stringr str_match
#' @importFrom DT datatable
metshiTable <- function(content, options=NULL, rownames= T){
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
  mz_rownames = stringr::str_match(rownames(content),
                                   "(\\d+\\.\\d+)")[,2]
  if(!is.na(mz_rownames[1])){
    rownames(content) <- paste0(rownames(content), sapply(rownames(content), function(mz) if(grepl("\\-", mz)) "" else "+"))
  }
  DT::datatable(content,
                selection = 'single',
                autoHideNavigation = T,
                class = 'compact', height = "500px",
                extensions = c("FixedColumns", "Scroller", "Buttons"), 
                options = opts,
                rownames = rownames
  )
}

#' @title Get top hits of an experiment
#' @description Goes through mSet storage to find experiments of interest, then takes the top x m/z values and returns them. Used mostly for Venn Diagram creation.
#' @param mSet mSet object
#' @param expnames Experiment names
#' @param top Number of top m/z values to return to user
#' @return List object with top hits per experiment
#' @seealso 
#'  \code{\link[stringr]{str_match}}
#'  \code{\link[data.table]{as.data.table}}
#' @rdname getTopHits
#' @export 
#' @importFrom stringr str_match
#' @importFrom data.table as.data.table
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
                     {metshiAlert("Not currently supported...")
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

getPlots <- function(do, mSet, input, gbl, lcl, venn_yes, my_selection){
  toWrap <- switch(do,
                   general = {
           # make sidebar
           # make pca, plsda, ml(make plotmanager do that)
           # update select input bars with current variable and covariables defined in excel
           if(!is.null(mSet)){
             varNormPlots <- MetaboShiny::ggplotNormSummary(mSet = mSet,
                                                            cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             
             sampNormPlots <- MetaboShiny::ggplotSampleNormSummary(mSet,
                                                                   cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             
             list(var1=varNormPlots$tl, var2=varNormPlots$bl,
                  var3=varNormPlots$tr, var4=varNormPlots$br,
                  samp1=sampNormPlots$tl, samp2=sampNormPlots$bl, 
                  samp3=sampNormPlots$tr, samp4=sampNormPlots$br)
             
           }},
         venn = {
           # get user input for how many top values to use for venn
           top = input$venn_tophits
           if(nrow(venn_yes$now) > 4 | nrow(venn_yes$now) <= 1){
             MetaboShiny::metshiAlert("Can only take more than 2 and less than five analyses!")
             list()
           }else{
             p <- MetaboShiny::ggPlotVenn(mSet = mSet,
                                          venn_yes = as.list(venn_yes),
                                          top = input$venn_tophits,
                                          cols = lcl$aes$mycols,
                                          cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             lcl$vectors$venn_lists <- p$info
             list(venn_plot = p$plot)
           }
         },
         enrich = {
           p = MetaboShiny::ggPlotMummi(mSet$analSet$enrich, 
                                        if(input$mummi_enr_method) "mummichog" else "gsea",
                                        cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
           list(enrich_plot = p)
         },
         summary = {
           p = MetaboShiny::ggplotSummary(mSet, my_selection$mz, 
                                          shape.fac = input$shape_var, 
                                          cols = lcl$aes$mycols, 
                                          cf=gbl$functions$color.functions[[lcl$aes$spectrum]],
                                          styles = input$ggplot_sum_style,
                                          add_stats = input$ggplot_sum_stats, 
                                          color.fac = input$col_var, 
                                          text.fac = input$txt_var)
           
           list(summary_plot = p)
         },
         pattern = {
           p = MetaboShiny::ggPlotPattern(mSet,
                                          n = input$pattern_topn,
                                          cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
           
           list(pattern_plot = p)
         },
         aov = { # render manhattan-like plot for UI
           p = MetaboShiny::ggPlotAOV(mSet,
                                      cf = gbl$functions$color.functions[[lcl$aes$spectrum]], 20)
           
           list(aov_plot = p)
         },
         volc = {
           # render volcano plot with user defined colours
           p = MetaboShiny::ggPlotVolc(mSet,
                                       cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                       20)
           
           list(volc_plot = p)
         },
         tsne = {
           mode <- if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
             "ipca" # interactive PCA (old name, i like tpca more :P )
           }else{
             "normal" # normal pca
           }
           
           if("tsne" %in% names(mSet$analSet)){
             if(input$tsne_2d3d | !input$ggplotly){ # check if switch button is in 2d or 3d mode
               # render 2d plot
               p = MetaboShiny::plotPCA.2d(mSet, 
                                           cols = lcl$aes$mycols,
                                           pcx = 1,
                                           pcy = 2, 
                                           type = "tsne",
                                           mode = mode,
                                           shape.fac = input$shape_var,
                                           col.fac = input$col_var
               )
               
             }else{
               # render 3d plot
               p = MetaboShiny::plotPCA.3d(mSet, 
                                           lcl$aes$mycols,
                                           pcx = 1,
                                           pcy = 2,
                                           pcz = 3,
                                           type = "tsne",
                                           mode = mode,
                                           shape.fac = input$shape_var,
                                           font = lcl$aes$font,
                                           col.fac = input$col_var,
                                           cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             }
           }
           list(tsne_plot = p)
         },
         pca = {
           if("pca" %in% names(mSet$analSet)){
             scree = MetaboShiny::ggPlotScree(mSet,
                                              cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             # chekc which mode we're in
             mode <- if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
               "ipca" # interactive PCA (old name, i like tpca more :P )
             }else{
               "normal" # normal pca
             }
             
             if(input$pca_2d3d | !input$ggplotly){ # check if switch button is in 2d or 3d mode
               # render 2d plot
               pca = MetaboShiny::plotPCA.2d(mSet, 
                                             cols = lcl$aes$mycols,
                                             pcx = input$pca_x,
                                             pcy = input$pca_y, 
                                             mode = mode,
                                             type = "pca",
                                             col.fac = input$col_var,
                                             shape.fac = input$shape_var)
               loadings = MetaboShiny::plotPCAloadings.2d(mSet,pcx = input$pca_x,
                                                          pcy = input$pca_y, 
                                                          type = "pca",
                                                          cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             }else{
               # render 3d plot
               pca = MetaboShiny::plotPCA.3d(mSet, 
                                             pcx = input$pca_x,
                                             pcy = input$pca_y,
                                             pcz = input$pca_z, 
                                             type = "pca",
                                             col.fac = input$col_var,
                                             mode = mode,
                                             cols = lcl$aes$mycols,
                                             shape.fac = input$shape_var,
                                             font = lcl$aes$font,
                                             cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
               loadings = MetaboShiny::plotPCAloadings.3d(mSet,
                                                          pcx = input$pca_x,
                                                          pcy = input$pca_y,
                                                          pcz = input$pca_z, 
                                                          type = "pca",
                                                          font = lcl$aes$font,
                                                          cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             }
             list(plot_pca = pca, plot_pca_loadings = loadings, pca_scree = scree)
           }else{
             NULL
           }
         },
         plsda = {
           
           if("plsda" %in% names(mSet$analSet)){ # if plsda has been performed...
             
             mode <- if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
               "ipca" # interactive PCA (old name, i like tpca more :P )
             }else{
               "normal" # normal pca
             }
             
             # render cross validation plot
             cv <- MetaboShiny::ggPlotClass(mSet, cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             # render permutation plot
             perm <- MetaboShiny::ggPlotPerm(mSet,cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             # see PCA - render 2d or 3d plots, just with plsda as mode instead
             if(input$plsda_2d3d | !input$ggplotly){
               # 2d
               plsda <- MetaboShiny::plotPCA.2d(mSet, lcl$aes$mycols,
                                                pcx = input$plsda_x,
                                                pcy = input$plsda_y, 
                                                type = "plsda",
                                                mode = mode,
                                                col.fac = input$col_var,
                                                shape.fac = input$shape_var,
                                                cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
               
               loadings <- MetaboShiny::plotPCAloadings.2d(mSet,pcx = input$plsda_x,
                                                           pcy = input$plsda_y, 
                                                           type = "plsda",
                                                           cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             }else{
               # 3d
               plsda <- MetaboShiny::plotPCA.3d(mSet, lcl$aes$mycols,
                                                pcx = input$plsda_x,
                                                pcy = input$plsda_y,
                                                pcz = input$plsda_z, 
                                                type = "plsda",
                                                mode = mode,
                                                col.fac = input$col_var,
                                                shape.fac = input$shape_var,
                                                font = lcl$aes$font, 
                                                cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
               
               loadings <- MetaboShiny::plotPCAloadings.3d(mSet,
                                                           pcx = input$plsda_x,
                                                           pcy = input$plsda_y,
                                                           pcz = input$plsda_z, 
                                                           type = "plsda",
                                                           font = lcl$aes$font,
                                                           cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             }
             list(plot_plsda = plsda, 
                  plot_plsda_loadings = loadings,
                  plsda_cv_plot = cv,
                  plsda_perm_plot = perm)
           }else{NULL}
         },
         ml = {
           if("ml" %in% names(mSet$analSet)){
             ml_roc = MetaboShiny::ggPlotROC(data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$roc,
                                             attempts = input$ml_attempts,
                                             cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             
             barplot_data <- MetaboShiny::ggPlotBar(data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$bar,
                                                     attempts = input$ml_attempts,
                                                     cf =gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                     topn = input$ml_top_x,
                                                     ml_name = mSet$analSet$ml$last$name,
                                                     ml_type = mSet$analSet$ml$last$method)
             
             ml_barplot <- barplot_data$plot
             lcl$tables$ml_bar <- barplot_data$mzdata
           }
           list(ml_roc = ml_roc, 
                ml_bar = ml_barplot)
         },
         multigroup = {
           p = MetaboShiny::ggplotSummary(mSet, my_selection$mz, 
                                          shape.fac = input$shape_var, 
                                          cols = lcl$aes$mycols, 
                                          cf=gbl$functions$color.functions[[lcl$aes$spectrum]], 
                                          mode = "multi",
                                          styles = input$ggplot_sum_style,
                                          add_stats = input$ggplot_sum_stats, color.fac = input$col_var, text.fac = input$txt_var)
           list(summary_plot = p)
         },
         meba = {
           p = MetaboShiny::ggplotMeba(mSet, my_selection$mz,
                                       draw.average = T,
                                       cols = lcl$aes$mycols,
                                       cf=gbl$functions$color.functions[[lcl$aes$spectrum]])
           list(summary_plot = p)
         },
         tt = {
           # render manhattan-like plot for UI
           p = MetaboShiny::ggPlotTT(mSet,
                                     cf = gbl$functions$color.functions[[lcl$aes$spectrum]], 
                                     20)
           
           list(tt_plot = p)
         },
         fc = {
           # render manhattan-like plot for UI
           p <- MetaboShiny::ggPlotFC(mSet,
                                      gbl$functions$color.functions[[lcl$aes$spectrum]], 20)
           
           list(fc_plot = p)
         },
         network = {
           
           pval = input$network_sign
           
           mat = mSet$analSet$network$rcorr
           matp = mat$P
           matr = mat$r
           matr[matp < pval] <- 0
           igr = igraph::graph.adjacency(adjmatrix = matr,
                                         weighted = T,
                                         diag = F)
           netw = visNetwork::visIgraph(igr)
           
           cf = gbl$functions$color.functions[[lcl$aes$spectrum]]
           cols = cf(nrow(netw$x$nodes))
           
           netw$x$nodes$color <- sapply(netw$x$nodes$id, function(id){
             cols[which(mSet$analSet$network$order == id)]
           })
           netw$x$nodes$title = netw$x$nodes$label
           
           netw$x$edges$value <- netw$x$edges$weight
           netw$x$edges$title <- paste0("r=",netw$x$edges$weight)
           
           netw$x$edges$color <- c("black")
           
           p = netw %>% 
             visNetwork::visOptions(highlightNearest = TRUE,
                                    collapse = T,
                                    autoResize = T,
                                    nodesIdSelection = TRUE) %>% 
             visNetwork::visNodes(borderWidth = 2,font = list(size=22,
                                                              face=lcl$aes$font$family)) %>%
             visNetwork::visEdges(scaling = list(min = 0.5,
                                                 max = 6)) %>%
             visNetwork::visInteraction(#navigationButtons = TRUE,
               keyboard = T)
           if(!input$network_auto){
             if(input$network_style == "hierarchical"){
               p = p %>% visNetwork::visHierarchicalLayout()  
             }else{
               p = p %>% visNetwork::visIgraphLayout(layout = input$network_style)
             }
           }
           # - - - - - - -
           matr[matr==1 | lower.tri(matr)] <- NA
           
           p2 = heatmaply::heatmaply(matr,
                                     Colv = T,
                                     Rowv = T,
                                     branches_lwd = 0.3,
                                     margins = c(0, 0, 0, 0),
                                     col = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                     column_text_angle = 90,
                                     ylab = "m/z\n",
                                     showticklabels = if(ncol(matr) <= 95) c(F,T) else c(F,F),
                                     symm=T,
                                     symbreaks=T,
                                     dendrogram="none"
           )
           lcl$vectors$network_heatmap <- p2$x$layout$yaxis$ticktext
           
           # - - - - - - - 
           list(network = p, 
                network_heatmap = p2)
         },
         heatmap = {
           
           breaks = seq(min(mSet$dataSet$norm), 
                        max(mSet$dataSet$norm), 
                        length = 256/2)
           
           p = {
             
             if(!is.null(mSet$analSet$heatmap$matrix)){
               # create heatmap object
               hmap <- suppressWarnings({
                 if(input$heatlimits){
                   heatmaply::heatmaply(mSet$analSet$heatmap$matrix[1:if(input$heatmap_topn < nrow(mSet$analSet$heatmap$matrix)) input$heatmap_topn else nrow(mSet$analSet$heatmap$matrix),],
                                        Colv = mSet$analSet$heatmap$my_order,
                                        Rowv = T,
                                        branches_lwd = 0.3,
                                        margins = c(60, 0, NA, 50),
                                        col = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                        col_side_colors = mSet$analSet$heatmap$translator[,!1],
                                        col_side_palette = mSet$analSet$heatmap$colors,
                                        subplot_widths = c(.9,.1),
                                        subplot_heights = if(mSet$analSet$heatmap$my_order) c(.1, .05, .85) else c(.05,.95),
                                        column_text_angle = 90,
                                        xlab = "Sample",
                                        ylab = "m/z",
                                        showticklabels = c(T,F),
                                        limits = c(min(mSet$dataSet$norm), max(mSet$dataSet$norm)),
                                        #symm=F,symkey=F,
                                        symbreaks=T
                                        #label_names = c("m/z", "sample", "intensity") #breaks side colours...
                   )
                 }else{
                   heatmaply::heatmaply(mSet$analSet$heatmap$matrix[1:if(input$heatmap_topn < nrow(mSet$analSet$heatmap$matrix)) input$heatmap_topn else nrow(mSet$analSet$heatmap$matrix),],
                                        Colv = mSet$analSet$heatmap$my_order,
                                        Rowv = T,
                                        branches_lwd = 0.3,
                                        margins = c(60, 0, NA, 50),
                                        colors = gbl$functions$color.functions[[lcl$aes$spectrum]](256),
                                        col_side_colors = mSet$analSet$heatmap$translator[,!1],
                                        col_side_palette = mSet$analSet$heatmap$colors,
                                        subplot_widths = c(.9,.1),
                                        subplot_heights = if(mSet$analSet$heatmap$my_order) c(.1, .05, .85) else c(.05,.95),
                                        column_text_angle = 90,
                                        xlab = "Sample",
                                        ylab = "m/z",
                                        #showticklabels = c(T,F),
                                        #symm=F,symkey=F,
                                        symbreaks=T
                                        #label_names = c("m/z", "sample", "intensity") #breaks side colours...
                   )
                 }
               })
               hmap$x$layout$annotations[[1]]$text <- ""
               # save the order of mzs for later clicking functionality
               lcl$vectors$heatmap <- hmap$x$layout[[if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")) "yaxis2" else "yaxis3"]]$ticktext 
               # return
               hmap
             }else{
               data = data.frame(text = "No significant hits available!\nPlease try alternative source statistics below.")
               ggplot2::ggplot(data) + ggplot2::geom_text(ggplot2::aes(label = text), x = 0.5, y = 0.5, size = 10) +
                 ggplot2::theme(text = ggplot2::element_text(family = lcl$aes$font$family)) + ggplot2::theme_bw()
             }
           }
           list(heatmap = p)
         }, 
         power = {
           p = {
             if("power" %in% names(mSet$analSet)){
               MetaboShiny::ggPlotPower(mSet, 
                                        max_samples = max(mSet$analSet$power[[1]]$Jpred),
                                        comparisons = names(mSet$analSet$power),
                                        cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             }else{
               NULL
             }
           }
           list(power_plot = p)
         },
         wordcloud = {
           if(nrow(lcl$tables$wordcloud_filt) > 0){
             topWords = if(input$wordcloud_topWords > nrow(lcl$tables$wordcloud_filt)) nrow(lcl$tables$wordcloud_filt) else input$wordcloud_topWords
             output$wordcloud <-  wordcloud2::renderWordcloud2({ 
               wordcloud2::wordcloud2(data.table::as.data.table(lcl$tables$wordcloud_filt)[order(n, decreasing = T)][1:topWords,], color = "random-light", size=.7, shape = "circle")
             })
             p = {
               wcdata = data.table::as.data.table(lcl$tables$wordcloud_filt)[order(n, decreasing = T)][1:topWords,]
               colnames(wcdata)[2] <- "freq"
               MetaboShiny::ggPlotWordBar(wcdata = wcdata,
                                          cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
             }
           }
           list(wordbar = p)
         }
  )
  
  finalPlots <- mapply(function(myplot, plotName){
    
    if(grepl("aov|tt|fc|pattern|asca|volc", plotName)){
      whichAnal <- stringr::str_match(plotName, "aov|tt|fc|pattern|asca|volc")[,1]
      if(is.null(mSet$analSet[[whichAnal]]$sig.mat)){
        data = data.frame(text = "No significant hits!")
        myplot = ggplot2::ggplot(data) + ggplot2::geom_text(ggplot2::aes(label = text), x = 0.5, y = 0.5, size = 10)
      }
    }
    
    isSquare <- grepl("pca|plsda|tsne|roc|heatmap|var|samp|network", plotName) & !grepl("scree|cv|perm|venn", plotName)
    
    # === WRAPPER ===
    
    canBe3D <- grepl("pca|plsda|tsne", plotName) & !grepl("scree|perm|cv", plotName)
    if(canBe3D){
      whichAnal <- stringr::str_match(plotName, "pca|plsda|tsne")[,1]
      is3D <- !input[[paste0(whichAnal, "_2d3d")]]
    }else{
      is3D <- plotName %in% c("heatmap", "network", "network_heatmap")
    }
    
    if(!is3D){
      myplot <- myplot + 
        gbl$functions$plot.themes[[lcl$aes$theme]](base_size = 15) + 
        ggplot2::theme(legend.position="none",
                       axis.line = ggplot2::element_line(colour = 'black', size = .5),
                       plot.title = ggplot2::element_text(hjust = 0.5,
                                                          vjust = 0.1,
                                                          size=lcl$aes$font$title.size*1.2),
                       text = ggplot2::element_text(family = lcl$aes$font$family))
      if(grepl("venn", plotName)){
        myplot <- myplot + 
          ggplot2::theme_void() +
          ggplot2::theme(panel.grid = ggplot2::element_blank(),
                         legend.position="none",
                         text = ggplot2::element_text(family = lcl$aes$font$family))
        
      }
      
      if(length(myplot$data) > 0){
        data = myplot$data
      }else{
        data = myplot[["layers"]][[1]][["data"]]
      }
      
      if(any(mSet$report$mzStarred$star) & any(grepl("mz|m/z", names(data)))){
        #if(grepl("tt|fc|volc|pca_load|plsda_load|ml_bar", plotName)){
        myX = myplot[["layers"]][[1]][["mapping"]][["x"]][[2]]
        myY = myplot[["layers"]][[1]][["mapping"]][["y"]][[2]]
        myText = myplot[["layers"]][[1]][["mapping"]][["text"]][[2]]
        flip = grepl("tt|fc|aov|var|samp|pattern", plotName)
        
        myplot <- myplot + 
          ggplot2::geom_text(ggplot2::aes(x = data[[myX]],
                        y = data[[myY]],
                        label = sapply(data[[myText]], function(mz){
                          ifelse(mz %in% mSet$report$mzStarred[star == TRUE]$mz,
                                 as.character(mz),'')
                        }),size=2),
                    nudge_y = {if(flip) 0 else .04} * max(data[[myY]]),
                    nudge_x = {if(flip) .04 else 0} * max(data[[myX]]),
                    #vjust="inward",hjust="inward"
          )  
        #}
      }
    }
    
     try({
        if(plotName %in% c("heatmap", "network_heatmap", "network")){
          data = data.frame(text = "Currently only available in 'plotly' mode!\nPlease switch in the sidebar.")
          myplot <- ggplot2::ggplot(data) + ggplot2::geom_text(ggplot2::aes(label = text), x = 0.5, y = 0.5, size = 10) +
            ggplot2::theme(text = ggplot2::element_text(family = lcl$aes$font$family)) + ggplot2::theme_bw()
        }
      }, silent = F)

    finalPlot = list(myplot)
    finalPlot
  }, toWrap, names(toWrap))
  
  list(lcl = lcl, plots = finalPlots)
}
# # old remotes field
# Remotes:
#   github::UMCUGenetics/MetaDBparse,
# github::xia-lab/MetaboAnalystR,
# github::rwehrens/BatchCorrMetabolomics,
# github::yixuan/showtext,
# github::joannawolthuis/ggVennDiagram,
# bioc::xcms, 
# bioc::CAMERA, 
# bioc::MSnbase,
# bioc::impute, 
# bioc::pcaMethods, 
# bioc::siggenes, 
# bioc::globaltest, 
# bioc::GlobalAncova, 
# bioc::Rgraphviz, 
# bioc::KEGGgraph, 
# bioc::preprocessCore, 
# bioc::genefilter, 
# bioc::SSPA, 
# bioc::sva, 
# bioc::ChemmineR, 
# bioc::KEGGREST,
# bioc::fgsea,
# bioc::Rdisop,
# bioc::Biostrings,
# bioc::GOSemSim,
# bioc::fmcsR