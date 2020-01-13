#' @export
make.metshi.csv <-
  function(patdb,
           csv = gsub(patdb, 
                      pattern = "\\.db", 
                      replacement = "\\.csv"),
           max.vals = -1,
           group_adducts = F,
           which_dbs = c(),
           which_adducts = c(),
           groupfac = "mz"
  ){
    
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), normalizePath(patdb))
    
    shiny::showNotification("Checking for mismatches between peak tables and metadata...")
    
    fn_meta <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT sample FROM individual_data")[,1]
    fn_int <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT filename FROM mzintensities")[,1]
    
    MetaboShiny::metshiAlert(paste0(
      "Missing from metadata:", paste0(setdiff(fn_int,fn_meta), collapse=", ")),
      "Missing from peaklist:", paste0(setdiff(fn_meta,fn_int), collapse=", "))
    
    if(DBI::dbExistsTable(conn, "setup")){
      query <- strwrap(gsubfn::fn$paste("select distinct d.*, s.*
                                         from mzintensities i
                                         join individual_data d
                                         on i.filename = d.sample
                                         join setup s on d.[group] = s.[group]"),
                       width=10000,
                       simplify=TRUE)  
    }else{
      query <- strwrap(gsubfn::fn$paste("select distinct d.*
                                        from mzintensities i
                                        join individual_data d
                                        on i.filename = d.sample"),
                       width=10000,
                       simplify=TRUE)
    }
    
    RSQLite::dbExecute(conn, "PRAGMA journal_mode=WAL;")
    RSQLite::dbExecute(conn, "CREATE INDEX IF NOT EXISTS filenames ON mzintensities(filename)")
    all_mz = RSQLite::dbGetQuery(conn, "SELECT DISTINCT mzmed FROM mzvals")[,1]
    
    RSQLite::dbDisconnect(conn)
    
    if(file.exists(csv)) file.remove(csv)
    
    # write rows to csv
    pbapply::pblapply(fn_meta, 
                      #cl = session_cl, 
                      function(filename){
      # connect
      conn <- RSQLite::dbConnect(RSQLite::SQLite(), normalizePath(patdb))
      
      # adjust query
      query_add = gsubfn::fn$paste(" WHERE i.filename = '$filename'")
      
      z.meta = data.table::as.data.table(RSQLite::dbGetQuery(conn, paste0(query, query_add)))
      duplicate <- which(duplicated(gsub(colnames(z.meta), pattern = "\\.+\\d+$", replacement="")))
      z.meta <- z.meta[1,-duplicate,with=F]
      
      if(nrow(z.meta)==0) return(NA)
      
      colnames(z.meta) <- tolower(colnames(z.meta))
      z.int = data.table::as.data.table(RSQLite::dbGetQuery(conn, 
                                          paste0("SELECT DISTINCT 
                                                i.mzmed as identifier,
                                                i.intensity
                                                FROM mzintensities i", query_add)))
      if(nrow(z.int)==0) return(NA)
      
      missing_mz <- setdiff(all_mz, z.int$identifier)
      
      # cast to wide
      cast.dt <- data.table::dcast.data.table(z.int,
                                              formula = ... ~ identifier,
                                              fun.aggregate = sum,
                                              value.var = "intensity")
      
      complete = as.numeric(cast.dt[1,])
      names(complete) = colnames(cast.dt)
      
      missing = rep(NA, length(missing_mz))
      names(missing) <- missing_mz
      
      complete.row = c(complete[-1], missing)
      reordered <- order(as.numeric(names(complete.row)))
      
      RSQLite::dbDisconnect(conn)
      
      # write
      data.table::fwrite(c(z.meta, "."=NA, complete.row[reordered]), 
             file = csv,
             append = T)
    })
    
    # - - measure file size - -
    disk_size = file.info(csv)$size
    size <- utils:::format.object_size(disk_size, "Mb")
    cat(paste("... Resulting file is approximately"),size,"...")
    
    }

calcHeatMap <- function(mSet, signif.only, 
                        source.tbl, 
                        top.hits, cols){
  sigvals = NULL
  
  try({
    switch(source.tbl,
           tt={
             decreasing <- F
             if(signif.only){
               tbl = as.data.frame(mSet$analSet$tt$sig.mat)
               sigvals = tbl$p.value
             }else{
               tbl = as.data.frame(mSet$analSet$tt$p.value)
               sigvals = tbl[,1]
             }
           },
           fc={
             decreasing <- T
             if(signif.only){
               tbl <- as.data.frame(mSet$analSet$fc$sig.mat)
               sigvals = abs(tbl$`log2(FC)`)
             }else{
               tbl = as.data.frame(abs(mSet$analSet$fc$fc.log))
               sigvals = tbl[,1]
             }},
           aov = {
             decreasing=F
             tbl=as.data.frame(mSet$analSet$aov$sig.mat)
             sigvals = tbl$p.value
           },
           aov2 = {
             decreasing=F
             tbl = as.data.frame(mSet$analSet$aov2$sig.mat)
             sigvals = tbl$`Interaction(adj.p)`
           },
           asca = {
             decreasing=T
             tbl = as.data.frame(mSet$analSet$asca$sig.list$Model.ab)
             sigvals = tbl$Leverage
           },
           meba = {
             decreasing=T
             tbl = as.data.frame(mSet$analSet$MB$stats)
             sigvals = tbl$`Hotelling-T2`
           })
  })
  
  # change top hits used in heatmap depending on time series / bivariate / multivariate mode
  # reordering of hits according to most significant at the top
  if(!is.null(sigvals)){
    #check top x used (slider bar in UI), if more than total matches use total matches
    topn = if(length(sigvals) < top.hits) length(sigvals) else top.hits
    mzorder <- order(sigvals, decreasing = decreasing)
    mzsel <- rownames(tbl)[mzorder]#[1:topn]
    
    # reorder matrix used
    x <- mSet$dataSet$norm[,mzsel]
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
                  train_vec,
                  test_vec,
                  configCols,
                  ml_method,
                  ml_perf_metr,
                  ml_folds,
                  ml_preproc,
                  tuneGrid,
                  ml_train_perc){
  
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
    reTrain <- caret::createDataPartition(y = config[train_idx, label], p = ml_train_perc) # take a user-defined percentage of the regexed training set
    inTrain <- train_idx[reTrain$Resample1]
    inTest = test_idx
  }else{ # ONLY TEST IS DEFINED
    test_idx = which(config[,test_vec[1], with=F][[1]] == test_vec[2])
    train_idx = setdiff(1:nrow(curr), test_idx) # use the other rows for testing
    reTrain <- caret::createDataPartition(y = config[train_idx, label], p = ml_train_perc) # take a user-defined percentage of the regexed training set
    inTrain <- train_idx[reTrain$Resample1]
    inTest <- test_idx
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
                                     returnData = FALSE) # need something here...
    
  }else{
    trainCtrl <- caret::trainControl(verboseIter = T,
                                     allowParallel = F,
                                     method=as.character(ml_perf_metr),
                                     number=as.numeric(ml_folds),
                                     repeats=3,
                                     trim=TRUE, 
                                     returnData = FALSE) # need something here...
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
               AUC = as.numeric(roc$auc),
               comparison = paste0(roc.pair$levels,collapse=" vs. "))
  })) 
}

getColDistribution <- function(csv){
  suppressWarnings({
    as.numi <- as.numeric(colnames(csv))
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
  unique.levels <- apply(csv[!excl.rows,..exp.vars, with=F], MARGIN=2, function(col){
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
  samps <- which(!grepl(csv$sample, pattern = "QC"))
  batchnum <- unique(csv[samps, "batch"][[1]])
  keep_samps_post_qc <- covar_table[which(covar_table$batch %in% batchnum),"sample"][[1]]
  covar_table <- covar_table[which(covar_table$batch %in% batchnum),]
  csv[which(csv$sample %in% keep_samps_post_qc),-"batch"]  
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
  w.missing <- mSet$dataSet$preproc
  w.missing <- apply(w.missing, 2, as.numeric)
  new.mat <- apply(w.missing, 1, function(x) {
    if(all(is.na(x))){
      x = c(0)
    }else{
      if (sum(is.na(x)) > 0) {
        x[is.na(x)] <- c(min(x[!is.na(x)], na.rm = T)/2)
      }  
    }
    x
  })
  t(new.mat)
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

batchCorrQC <- function(mSet, qc_rows){
  smps <- rownames(mSet$dataSet$norm)
  # get which rows are QC samples
  qc_rows <- which(grepl(pattern = "QC", x = smps))
  # get batch for each sample
  batch.idx = as.numeric(as.factor(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$sample),"batch"][[1]]))
  
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
