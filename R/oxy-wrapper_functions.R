globalVariables(c("-log(p)", "-log10P", "..change_var", "..col.fac", "..count..", "..exp.vars", 
                  "..keep.cols", "..matches", "..predictor", "..rmcols", "..scaled..", "..shape.fac", 
                  "..stats_var", "..time_var", "..x", "Abundance", "Color", "FPR", "Group", "GroupA", 
                  "GroupB", "Individual", "Label", "Metric", "PC", "Peak", "Sample", "Shape", "TPR", 
                  "Text", "Value", "abstract", "acc", "adduct", "aes", "attempt", "color", "comparison", 
                  "coord_flip", "correlation", "count", "debug_browse_content", "debug_enrich", "debug_input", 
                  "debug_lcl", "debug_mSet", "debug_matches", "debug_pieinfo", "debug_result_filters", "debug_selection", 
                  "exp.vars", "expand_limits", "extremity", "facet_grid", "facet_wrap", "freq", "fullformula", "gbl", 
                  "geom_point", "ggplot", "ggtitle", "group", "importance.mean", "individual", "label", "labs", "lcl", 
                  "log2FC", "log2fc", "m/z", "p-value", "pathway", "pc", "position_dodge", "position_jitterdodge", 
                  "samples", "scale_size_area", "scale_x_discrete", "scale_y_discrete", "searchRev", "shape",
                  "value", "variable", "x", "xlab", "y"))

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
                        cols,
                        which.data){
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
    
    inTbl = switch(which.data, 
                   original = mSet$dataSet$prog,
                   "pre-batch correction" = mSet$dataSet$prebatch,
                   normalized = mSet$dataSet$norm)
    
    x <- inTbl[,as.character(sigvals)]
    final_matrix <- t(x) # transpose so samples are in columns
    
    # check if the sample order is correct - mSet$..$ norm needs to match the matrix
    sample_order <- match(colnames(final_matrix), rownames(mSet$dataSet$norm))
    
    if(mSet$settings$exp.type %in% c("2f", "t1f", "t")){
      
      # create convenient table with the necessary info
      translator <- data.table::data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],
                                           GroupA=mSet$dataSet$facA[sample_order], 
                                           GroupB=mSet$dataSet$facB[sample_order])
      colnames(translator) <- c("Sample", mSet$dataSet$facA.lbl, mSet$dataSet$facB.lbl)
      
      hmap.lvls <- c(levels(mSet$dataSet$facA), levels(mSet$dataSet$facB))
      
      # reorder first by time, then by sample
      split.translator <- split(translator, by = c(mSet$dataSet$facB.lbl))
      split.translator.ordered <- lapply(split.translator, function(tbl) tbl[order(tbl[[mSet$dataSet$facA.lbl]])])
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
getMultiMLperformance <- function(x, type="roc"){
  
  try({
    mroc = pROC::multiclass.roc(x$labels,
                                x$prediction)
  },silent = F)
  # try({
  #   mroc = pROC::multiclass.roc(x$labels, factor(x$prediction,
  #                                                ordered = T))
  # }, silent=F)
  
  data.table::rbindlist(lapply(mroc$rocs, function(roc.pair){
    try({
      dt = data.table(FPR = sapply(roc.pair$specificities, function(x) 1-x),
                      TPR = roc.pair$sensitivities,
                      AUC_AVG = as.numeric(mroc$auc),
                      AUC_PAIR = as.numeric(pROC::auc(roc.pair)),
                      comparison = paste0(roc.pair$levels,collapse=" vs. "))
    },silent=T)
    try({
      dt = data.table(FPR = sapply(roc.pair[[1]]$specificities, function(x) 1-x),
                      TPR = roc.pair[[1]]$sensitivities,
                      AUC_AVG = as.numeric(mroc$auc),
                      AUC_PAIR = as.numeric(pROC::auc(roc.pair[[1]])),
                      comparison = paste0(roc.pair[[1]]$levels,collapse=" vs. "))  
    },silent=T)
    dt
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
  min.grpsize = apply(csv[,names(unique.levels), with=F], 2,function(x) min(table(x)))
  # use this as the default selected experimental variable (user can change later)
  good.grpsize = min.grpsize >= 3
  unique.levels = unique.levels[good.grpsize]
  which.default <- unique.levels[which(unique.levels == min(unique.levels[which(unique.levels >= min.lev)]))][1]
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
  int.mat <- qs::qread("preproc.qs")
  mSet$dataSet$proc <- pbapply::pbapply(int.mat, 1, function(x) {
    if (sum(is.na(x)) > 0) {
      x[is.na(x)] <- min(x, na.rm = T)/2
    }
    x
  })
  mSet$dataSet$proc <- t(mSet$dataSet$proc)
  return(mSet)
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
replRF <- function(mSet, parallelMode, ntree, cl, rf.method){
  print(rf.method)
  w.missing <- qs::qread("preproc.qs")
  samples <- rownames(w.missing)
  # convert all to as numeric
  # TODO: remove, should be automatic
  w.missing <- apply(w.missing, 2, as.numeric)
  
  # register other threads as parallel threads
  doParallel::registerDoParallel(cl)
  
  imp <- switch(rf.method,
         ranger = {
           print(w.missing[1:10,1:10])
           w.missing <<- w.missing
           w.missing.df <- as.data.frame(w.missing)
           colnames(w.missing.df) <- paste0("mz", 1:ncol(w.missing.df))
           imp = missRanger.joanna(data = w.missing.df ,formula =  . ~ ., verbose = 0, num.threads = length(cl))
           colnames(imp) <- colnames(w.missing)
         },
         rf = {
           auto.mtry <- floor(sqrt(ncol(w.missing)))
         
         mtry <- ifelse(auto.mtry > 100, 
                        100, 
                        auto.mtry)
         
         # impute missing values with random forest
         imp <- missForest::missForest(w.missing,
                                       parallelize = parallelMode, # parallelize over variables, 'forests' is other option
                                       verbose = F,
                                       ntree = ntree,
                                       mtry = mtry)
         imp$ximp})
  return(imp)
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
  qc_rows <- grep(pattern = "QC", x = smps)
  # save QCs (may need for dimred)
  mSet$dataSet$qc_norm <- mSet$dataSet$norm[qc_rows,]
  mSet$dataSet$qc_cls <- mSet$dataSet$cls[qc_rows, drop = TRUE]
  mSet$dataSet$qc_covars <- mSet$dataSet$covars[grep("QC",sample),]
  
  # hide from stats
  mSet$dataSet$norm <- mSet$dataSet$norm[-qc_rows,]
  mSet$dataSet$cls <- mSet$dataSet$cls[-qc_rows, drop = TRUE]
  mSet$dataSet$covars <- mSet$dataSet$covars[grep("QC",sample,invert = T),]
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
combatCSV <- function(mSet, tbl="norm"){
  # get sample names and classes
  smp <- rownames(mSet$dataSet[[tbl]])
  exp_lbl <- mSet$dataSet$cls
  csv = mSet$dataSet[[tbl]]
  csv_edata <- t(csv)
  colnames(csv_edata) <- rownames(mSet$dataSet[[tbl]])
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
                class = 'compact', height = "500px",
                extensions = c("FixedColumns", "Scroller", "Buttons"), 
                options = opts,
                rownames = rownames,
                escape = FALSE
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
getTopHits <- function(mSet, expnames, top, thresholds=c(), filter_mode="top"){

  experiments <- stringr::str_match(expnames, 
                                    pattern = "(all m\\/z)|\\(.*\\)")[,1]
  
  experiments <- unique(gsub(experiments, pattern = "\\(\\s*(.+)\\s*\\)", replacement="\\1"))

  exp_table = data.frame(name = expnames, threshold = c(if(length(thresholds)>0) thresholds else 0))
  
  table_list <- lapply(experiments, function(experiment){
    if(experiment == "all m/z"){
      flattened = list(colnames(mSet$dataSet$norm))
      names(flattened) = c("all m/z")
      flattened
    }else if(experiment == "random"){
      flattened = list(sample(colnames(mSet$dataSet$norm)))
      names(flattened) = c("random")
      flattened
    }else{
      analysis = mSet$storage[[experiment]]$analSet
      
      rgx_exp <- gsub(experiment, pattern = "\\(", replacement = "\\\\(")
      rgx_exp <- gsub(rgx_exp, pattern = "\\)", replacement = "\\\\)")
      rgx_exp <- gsub(rgx_exp, pattern = "\\-", replacement = "\\\\-")
      rgx_exp <- gsub(rgx_exp, pattern = "\\+", replacement = "\\\\+")
      
      categories = grep(unlist(expnames),
                        pattern = paste0("\\(",rgx_exp, "\\)"), value = T)
      # go through the to include analyses
      
      tables <- lapply(categories, function(name){
        
        name_orig = name
        filter = exp_table[exp_table$name == name_orig, "threshold"]
        sign = stringr::str_extract(filter,pattern = ">|<|=")
        thresh = as.numeric(gsub(sign, "", filter))
        
        name = gsub(name, pattern = " \\(\\s*(.+)\\s*\\)", replacement = "")
        
        base_name <- search_name <- gsub(name, pattern = " -.*$| ", replacement="")
        
        if(base_name %in% gbl$constants$ml.models){
          search_name <- "ml"
        }
        
        # fetch involved mz values
        tbls <- switch(search_name,
                       ml = {
                         which.ml <- gsub(name, pattern = "^.*- ", replacement="")
                         
                         data = analysis$ml[[base_name]][[which.ml]]
                         if(!is.null(data$res$prediction)){
                           data$res$shuffled = FALSE
                           data$res = list(data$res)
                         }
                         data.dt = data.table::as.data.table(data$res[[which(unlist(sapply(data$res, function(x) !x$shuffled)))]]$importance, keep.rownames=T)
                         
                         colnames(data.dt)[1:2] <- c("m/z","importance")
                         
                         data.ordered <- data.dt[order(importance, decreasing = T),]
                         
                         data.ordered$`m/z` <- gsub("^X","",data.ordered$`m/z`)
                         
                         res = list(data.frame(`m/z`=data.ordered$`m/z`,
                                               value=data.ordered$"importance"))
                         names(res) <- paste0(which.ml, " (", base_name, ")")
                         res
                       },
                       venn = {
                         res = list(data.frame(`m/z`=analysis$venn$mzs,
                                               value=c(0)))
                         names(res) = base_name
                         res
                       },
                       corr = {
                         values = analysis$corr$cor.mat[order(abs(analysis$corr$cor.mat[,1]),
                                                                         decreasing = F),]
                         res = list(data.frame(`m/z` = rownames(values),
                                               value = values[,2])
                         )
                         names(res) = base_name
                         res
                       },
                       aov = {
                         values = analysis$aov$sig.mat[order(analysis$aov$sig.mat[,2],
                                                             decreasing = F),]
                         res = list(data.frame(`m/z` = rownames(values),
                                               value = values[,2])
                         )
                         names(res) = base_name
                         res
                       },
                       aov2 = {
                         values = analysis$aov2$sig.mat[order(analysis$aov2$sig.mat[,"Interaction(adj.p)"],
                                                              decreasing = F),]
                         res = list(data.frame(`m/z` = rownames(values),
                                               value = values[,2])
                         )
                         names(res) = base_name
                         res
                       },
                       asca = {
                         values = analysis$asca$sig.list$Model.ab[order(analysis$asca$sig.list$Model.ab[,1],
                                                                        decreasing = T),]
                         res = list(data.frame(`m/z` = rownames(values),
                                               value = values[,2])
                         )
                         names(res) = base_name
                         res
                       },
                       MB = {
                         values = analysis$MB$stats[order(analysis$MB$stats[,1],
                                                            decreasing = F),]
                         res = list(data.frame(`m/z` = rownames(values),
                                               value = values[,2])
                         )
                         names(res) = base_name
                         res
                       },
                       tt = {
                         values = analysis$tt$sig.mat[order(analysis$tt$sig.mat[,2],
                                                          decreasing = F),]
                         res = list(data.frame(`m/z` = rownames(values),
                                               value = values[,2])
                                    )
                         names(res) = base_name
                         res
                       },
                       fc = {
                         values = analysis$fc$sig.mat[order(analysis$fc$sig.mat[,2],
                                                            decreasing = F),]
                         res = list(data.frame(`m/z` = rownames(values),
                                               value = values[,2])
                         )
                         names(res) = base_name
                         res
                       },
                       combi = {
                         # special opt for volcano... TODO: somehow make this work for the others
                         is.volc = all(colnames(analysis$combi$sig.mat) %in% c("rn",  "log2(FC)", "-log10(p)"))
                         res = if(is.volc){
                           print("special volcano plot mode")
                           abs.fc = abs(analysis$combi$sig.mat$`log2(FC)`)
                           comb.vals = abs.fc * analysis$combi$sig.mat$`-log10(p)`
                           res = analysis$combi$sig.mat[order(comb.vals,
                                                                 decreasing = T),]
                           res[,2:3] <- NULL
                           res$value <- c(0)
                           res = list(res)
                         }else{
                           res = list(data.frame("m/z"=analysis$combi$sig.mat$rn,
                                                 value=c(0)))
                         }
                         
                         names(res) = base_name
                         #print(res)
                         res
                       },
                       plsda = {
                         which.plsda <- gsub(name, pattern = "^.*- | ", replacement="")
                         
                         compounds_pc <- data.table::as.data.table(analysis$plsda$vip.mat,keep.rownames = T)
                         colnames(compounds_pc) <- c("rn", paste0("Component ", 1:(ncol(compounds_pc)-1)))
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
                       featsel = {
                         decision = analysis$featsel[[1]]$finalDecision
                         res = list(names(decision[decision != "Rejected"]))
                         names(res) <- "featsel"
                         res
                       },
                       return(NULL))
        
        if(is.null(tbls)) return(NULL)
        
        # user specified top hits only
        tbls_top <- lapply(tbls, function(tbl){
          filt_tbl = switch(filter_mode,
                 top = if(nrow(tbl) < top){
                   tbl[,1]
                 }else{
                   tbl[1:top, 1]
                 },
                 threshold = tbl[switch(sign,
                                        ">" = {tbl$value > thresh},
                                        "=" = {tbl$value == thresh},
                                        "<" = {tbl$value < thresh}),1])
          
        })
        keep_tbls = unlist(lapply(tbls_top, function(t) length(t)>0))
        tbls_top=tbls_top[keep_tbls]
        if(length(tbls_top)>0){
          names(tbls_top) <- paste0(experiment, ": ", names(tbls_top))  
          tbls_top
        }else{
          list()
        }
      })
      
      # unnest the nested lists
      flattened <- flattenlist(tables)
      
      # remove NAs
      flattened <- lapply(flattened, function(x) x[!is.na(x)])
      
      #rename and remove regex-y names
      names(flattened) <- gsub(x = names(flattened), pattern = "(.*\\.)(.*$)", replacement = "\\2")
      # return
      flattened      
    }
  })
  flattened <- flattenlist(table_list)
  names(flattened) <- gsub(x = names(flattened), pattern = "(.*\\.)(.*$)", replacement = "\\2")
  flattened <- lapply(flattened, function(x) x[!is.na(x)])
  return(flattened)
}

getAllHits <- function(mSet, expname, randomize = F){
  
  experiment <- stringr::str_match(expname, 
                                    pattern = "(all m\\/z)|\\(.*\\)")[,1]
  
  experiment <- unique(gsub(experiment, pattern = "\\(\\s*(.+)\\s*\\)", replacement="\\1"))
  
  exp_table = data.frame(name = expname, threshold = c(0))
  
  if(experiment == "all m/z"){
    flattened = list(colnames(mSet$dataSet$norm))
    names(flattened) = c("all m/z")
    flattened
  }else{
    analysis = mSet$storage[[experiment]]$analSet
    
    rgx_exp <- gsub(experiment, pattern = "\\(", replacement = "\\\\(")
    rgx_exp <- gsub(rgx_exp, pattern = "\\)", replacement = "\\\\)")
    rgx_exp <- gsub(rgx_exp, pattern = "\\-", replacement = "\\\\-")
    rgx_exp <- gsub(rgx_exp, pattern = "\\+", replacement = "\\\\+")
    
    name_orig = grep(expname,
                     pattern = paste0("\\(",rgx_exp, "\\)"), value = T)
    # go through the to include analyses
      
    filter = exp_table[exp_table$name == name_orig, "threshold"]
    sign = stringr::str_extract(filter,pattern = ">|<|=")
    thresh = as.numeric(gsub(sign, "", filter))
    
    name = gsub(name_orig, pattern = " \\(\\s*(.+)\\s*\\)", replacement = "")
    
    base_name <- search_name <- gsub(name, pattern = " -.*$| ", replacement="")
    
    if(base_name %in% gbl$constants$ml.models){
      search_name <- "ml"
    }
    
    # fetch involved mz values
    tbl <- switch(search_name,
                   venn = {
                     venn.mzs = analysis$venn$mzs
                     not.in.venn = setdiff(colnames(mSet$dataSet$norm), venn.mzs)
                     res = data.frame(`m/z` = c(venn.mzs, not.in.venn),
                                      value = c(rep(0, length(venn.mzs)),
                                                rep(1, length(not.in.venn)))
                                      #,statistic = c(0)
                                      )
                     
                     # -------------------------
                     res
                   },
                   tt = {
                     res = data.frame(`m/z` = names(analysis$tt$p.value),
                                      value = analysis$tt$p.value,
                                      statistic = analysis$tt$t.score
                     )
                     res$significant = sapply(res$m.z, function(mz) mz %in% rownames(analysis$tt$sig.mat))
                     res
                   },
                   fc = {
                     res = data.frame(`m/z` = names(analysis$fc$fc.all),
                                      value = analysis$fc$fc.all)
                     res$significant = sapply(res$m.z, function(mz) mz %in% rownames(analysis$fc$sig.mat))
                     names(res) = base_name
                     res
                   },
                   combi = {
                     # --- only volc for now ---
                     res = data.frame(`m/z` = names(analysis$combi$all.vals$x),
                                      value = sapply(names(analysis$combi$all.vals$x), function(x) if(x %in% analysis$combi$sig.mat$rn) 0 else 1),
                                      statistic = analysis$combi$all.vals$x * analysis$combi$all.vals$y)
                     res$significant = sapply(res$m.z, function(mz) mz %in% analysis$combi$sig.mat$rn)
                     # -------------------------
                     res
                   },
                  featsel = {
                    decision = analysis$featsel[[1]]$finalDecision
                    res = data.frame(`m/z` = names(decision),
                                     value = sapply(names(decision), function(x) ifelse(as.character(decision[x]) == "Rejected", 1, 0)),
                                     statistic = sapply(names(decision), function(x) switch(as.character(decision[x]), 
                                                                                            Rejected = 0, 
                                                                                            Tentative = 1, 
                                                                                            Confirmed = 2))
                                     )
                    res
                  },
                   ml = {
                     ml_name = gsub(paste0(base_name, " - "), "", name)
                     mdls = analysis$ml[[base_name]][[ml_name]]$res
                     final_nonshuffle = mdls[[which(sapply(mdls, function(mdl) !mdl$shuffled))]]
                     res = data.frame(`m/z` = gsub("\\.$", "-", gsub("^X", "", rownames(final_nonshuffle$importance))),
                                      value =  c(0),
                                      stastistic = final_nonshuffle$importance[,1])
                     res
                   },
                   {metshiAlert("Not currently supported...")
                     return(NULL)})

    try({
      tbl = tbl[order(abs(tbl$statistic),decreasing = T),]
    }, silent = F)
    
    if(nrow(tbl)>0){
      if(randomize){
        tbl[sample(1:nrow(tbl)),]
      }else{
        tbl  
      }
    }else{
      data.table::data.table()
    }
  }
}

getPlots <- function(do, mSet, input, gbl, lcl, venn_yes, my_selection){
                  toWrap <- switch(do,
                   general = {
                     # make sidebar
                     # make pca, plsda, ml(make plotmanager do that)
                     # update select input bars with current variable and covariables defined in excel
                     if(!is.null(mSet)){
                       varNormPlots <- ggplotNormSummary(mSet = mSet,
                                                         cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       
                       sampNormPlots <- ggplotSampleNormSummary(mSet,
                                                                cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       
                       list(var1=varNormPlots$tl, var2=varNormPlots$bl,
                            var3=varNormPlots$tr, var4=varNormPlots$br,
                            samp1=sampNormPlots$tl, samp2=sampNormPlots$bl, 
                            samp3=sampNormPlots$tr, samp4=sampNormPlots$br)
                       
                     }},
                    venn = {
                     # get user input for how many top values to use for venn
                     top = input$venn_tophits
                     if(nrow(venn_yes$now) > 7 | nrow(venn_yes$now) <= 1){
                       metshiAlert("Can only take more than 2 and less than seven analyses!")
                       list()
                     }else{
                       p <- ggPlotVenn(mSet = mSet,
                                       venn_yes = as.list(venn_yes),
                                       filter_mode = input$venn_filter_mode,
                                       top = input$venn_tophits,
                                       cols = lcl$aes$mycols,
                                       plot_mode = ifelse(input$venn_plot_mode, "upset", "venn"),
                                       cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       lcl$vectors$venn_lists <- p$info
                       list(venn_plot = p$plot)
                     }
                   },
                   enrich = {
                     p = ggPlotMummi(mSet, 
                                     cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                     plot_mode = if(input$enrich_plot_mode) "volclike" else "gsea",
                                     show_nonsig = T)
                     list(enrich_plot = p)
                   },
                   summary = {
                     p = ggplotSummary(mSet, my_selection$mz, 
                                       shape.fac = input$shape_var, 
                                       cols = lcl$aes$mycols, 
                                       cf=gbl$functions$color.functions[[lcl$aes$spectrum]],
                                       styles = input$ggplot_sum_style,
                                       add_stats = input$ggplot_sum_stats, 
                                       color.fac = input$col_var, 
                                       fill.fac = input$fill_var,
                                       text.fac = input$txt_var)
                     
                     list(summary_plot = p)
                   },
                   corr = {
                     p = ggPlotPattern(mSet,
                                       n = input$corr_topn,
                                       cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                     
                     list(corr_plot = p)
                   },
                   asca = {
                     p = ggPlotASCA(mSet,
                                    cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                     
                     list(asca_plot = p)
                   },
                   aov = { # render manhattan-like plot for UI
                     p = ggPlotAOV(mSet,
                                   cf = gbl$functions$color.functions[[lcl$aes$spectrum]], 20,
                                   topn=input$aov_topn)
                     
                     list(aov_plot = p)
                   },
                   volcano = {
                     # render volcano plot with user defined colours
                     p = ggPlotVolc(mSet,
                                    cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                    20)
                     
                     list(volcano_plot = p)
                   },
                   ica = {
                     mode <- if(mSet$settings$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
                       "ipca" # interactive PCA (old name, i like tpca more :P )
                     }else{
                       "normal" # normal pca
                     }
                     
                     if("ica" %in% names(mSet$analSet)){
                       if(input$ica_2d3d | !input$ggplotly){ # check if switch button is in 2d or 3d mode
                         # render 2d plot
                         p = plotPCA.2d(mSet, 
                                        cols = lcl$aes$mycols,
                                        pcx = input$ica_x,
                                        pcy = input$ica_y,
                                        type = "ica",
                                        mode = mode,
                                        shape.fac = input$shape_var,
                                        col.fac = input$col_var,
                                        fill.fac = input$fill_var, 
                                        ellipse = input$ica_ellipse)
                         
                       }else{
                         # render 3d plot
                         p = plotPCA.3d(mSet, 
                                        lcl$aes$mycols,
                                        pcx = input$ica_x,
                                        pcy = input$ica_y,
                                        pcz = input$ica_z,
                                        type = "ica",
                                        mode = mode,
                                        shape.fac = input$shape_var,
                                        font = lcl$aes$font,
                                        col.fac = input$col_var,
                                        fill.fac = input$fill_var, 
                                        cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                        ellipse = input$ica_ellipse)
                       }
                     }
                     list(ica_plot = p)
                   },
                   tsne = {
                     mode <- if(mSet$settings$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
                       "ipca" # interactive PCA (old name, i like tpca more :P )
                     }else{
                       "normal" # normal pca
                     }
                     
                     if("tsne" %in% names(mSet$analSet)){
                       if(input$tsne_2d3d | !input$ggplotly){ # check if switch button is in 2d or 3d mode
                         # render 2d plot
                         p = plotPCA.2d(mSet, 
                                        cols = lcl$aes$mycols,
                                        pcx = input$tsne_x,
                                        pcy = input$tsne_y,
                                        type = "tsne",
                                        mode = mode,
                                        shape.fac = input$shape_var,
                                        col.fac = input$col_var,
                                        fill.fac = input$fill_var, 
                                        ellipse = input$tsne_ellipse)
                         
                       }else{
                         # render 3d plot
                         p = plotPCA.3d(mSet, 
                                        lcl$aes$mycols,
                                        pcx = input$tsne_x,
                                        pcy = input$tsne_y,
                                        pcz = input$tsne_z,
                                        type = "tsne",
                                        mode = mode,
                                        shape.fac = input$shape_var,
                                        font = lcl$aes$font,
                                        col.fac = input$col_var,
                                        fill.fac = input$fill_var, 
                                        cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                        ellipse = input$tsne_ellipse)
                       }
                     }
                     list(tsne_plot = p)
                   },
                   umap = {
                     mode <- if(mSet$settings$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
                       "ipca" # interactive PCA (old name, i like tpca more :P )
                     }else{
                       "normal" # normal pca
                     }
                     
                     if("umap" %in% names(mSet$analSet)){
                       if(input$umap_2d3d | !input$ggplotly){ # check if switch button is in 2d or 3d mode
                         # render 2d plot
                         p = plotPCA.2d(mSet, 
                                        cols = lcl$aes$mycols,
                                        pcx = input$umap_x,
                                        pcy = input$umap_y,
                                        type = "umap",
                                        mode = mode,
                                        shape.fac = input$shape_var,
                                        col.fac = input$col_var,
                                        fill.fac = input$fill_var, 
                                        ellipse = input$umap_ellipse)
                         
                       }else{
                         # render 3d plot
                         p = plotPCA.3d(mSet, 
                                        lcl$aes$mycols,
                                        pcx = input$umap_x,
                                        pcy = input$umap_y,
                                        pcz = input$umap_z,
                                        type = "umap",
                                        mode = mode,
                                        shape.fac = input$shape_var,
                                        font = lcl$aes$font,
                                        col.fac = input$col_var,
                                        fill.fac = input$fill_var, 
                                        cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                        ellipse = input$umap_ellipse)
                       }
                     }
                     list(umap_plot = p)
                   },
                   
                   pca = {
                     if("pca" %in% names(mSet$analSet)){
                       scree = ggPlotScree(mSet,
                                           cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       # chekc which mode we're in
                       mode <- if(mSet$settings$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
                         "ipca" # interactive PCA (old name, i like tpca more :P )
                       }else{
                         "normal" # normal pca
                       }
                       
                       if(input$pca_2d3d | !input$ggplotly){ # check if switch button is in 2d or 3d mode
                         # render 2d plot
                         pca = plotPCA.2d(mSet, 
                                          cols = lcl$aes$mycols,
                                          pcx = input$pca_x,
                                          pcy = input$pca_y, 
                                          mode = mode,
                                          type = "pca",
                                          col.fac = input$col_var,
                                          fill.fac = input$fill_var, 
                                          shape.fac = input$shape_var, 
                                          ellipse = input$pca_ellipse)
                         loadings = plotPCAloadings.2d(mSet,pcx = input$pca_x,
                                                       pcy = input$pca_y, 
                                                       type = "pca",
                                                       cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       }else{
                         # render 3d plot
                         pca = plotPCA.3d(mSet, 
                                          pcx = input$pca_x,
                                          pcy = input$pca_y,
                                          pcz = input$pca_z, 
                                          type = "pca",
                                          col.fac = input$col_var,
                                          mode = mode,
                                          cols = lcl$aes$mycols,
                                          fill.fac = input$fill_var, 
                                          shape.fac = input$shape_var,
                                          font = lcl$aes$font,
                                          cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                          ellipse = input$pca_ellipse)
                         loadings = plotPCAloadings.3d(mSet,
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
                       
                       mode <- if(mSet$settings$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
                         "ipca" # interactive PCA (old name, i like tpca more :P )
                       }else{
                         "normal" # normal pca
                       }
                       
                       # render cross validation plot
                       cv <- ggPlotClass(mSet, cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       # render permutation plot
                       perm <- ggPlotPerm(mSet,cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       # see PCA - render 2d or 3d plots, just with plsda as mode instead
                       if(input$plsda_2d3d | !input$ggplotly){
                         # 2d
                         plsda <- plotPCA.2d(mSet, lcl$aes$mycols,
                                             pcx = input$plsda_x,
                                             pcy = input$plsda_y, 
                                             type = "plsda",
                                             mode = mode,
                                             col.fac = input$col_var,
                                             shape.fac = input$shape_var,
                                             fill.fac = input$fill_var, 
                                             cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                             ellipse = input$plsda_ellipse)
                         
                         loadings <- plotPCAloadings.2d(mSet,pcx = input$plsda_x,
                                                        pcy = input$plsda_y, 
                                                        type = "plsda",
                                                        cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       }else{
                         # 3d
                         plsda <- plotPCA.3d(mSet, lcl$aes$mycols,
                                             pcx = input$plsda_x,
                                             pcy = input$plsda_y,
                                             pcz = input$plsda_z, 
                                             type = "plsda",
                                             mode = mode,
                                             col.fac = input$col_var,
                                             shape.fac = input$shape_var,
                                             fill.fac = input$fill_var, 
                                             font = lcl$aes$font, 
                                             cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                             ellipse = input$plsda_ellipse)
                         
                         loadings <- plotPCAloadings.3d(mSet,
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
                     if("ml" %in% names(mSet$analSet) & 
                        !(input$ml_plot_posclass %in% c("placeholder", "")) & 
                          input$ml_plot_x != "" &
                          input$ml_plot_y != ""){

                       if(length(mSet$analSet$ml) > 0){
                         
                         data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]

                         if(!is.null(data$res$prediction)){
                           data$res$shuffled = FALSE
                           data$res = list(data$res)
                         }
                         
                         # PLOT #
                         ml_performance_rows = lapply(1:length(data$res), function(i){
                           res = data$res[[i]]
                           ml_performance = getMLperformance(ml_res = res, 
                                                             pos.class = input$ml_plot_posclass,
                                                             x.metric = input$ml_plot_x,
                                                             y.metric = input$ml_plot_y)
                           ml_performance$coords$shuffled = c(res$shuffled)
                           ml_performance$coords$run = i
                           plot_coords = ml_performance$coords[`Test set`=="Test"]
                           ml_performance
                         })
                         
                         coords = data.table::rbindlist(lapply(ml_performance_rows,
                                                               function(x) x$coords))
                         
                         ml_performance = list(coords = coords,
                                               names = ml_performance_rows[[1]]$names)
                         
                         ml_roc = ggPlotCurves(ml_performance,
                                               cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                         
                         no_shuffle_imp = data$res[[which(unlist(sapply(data$res, function(x) !x$shuffle)))]]$importance
                         
                         barplot_data <- ggPlotBar(data = no_shuffle_imp,
                                                   cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                   topn = input$ml_topn,
                                                   ml_name = data$params$ml_name,
                                                   ml_type = data$params$ml_method)
                         
                         ml_barplot <- barplot_data$plot
                         lcl$tables$ml_bar <- barplot_data$mzdata
                         
                         list(ml_roc = ml_roc, 
                              ml_bar = ml_barplot)
                         
                       }} else list()
                   },
                   ml_mistake = {
                     if("ml" %in% names(mSet$analSet)){
                       if(length(mSet$analSet$ml) > 0){
                         data =  mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]
                         mistake_plot <- if(input$ml_mistake_var != "") ggPlotMLMistakes(labels = data$roc$labels,
                                                                                         predictions = data$roc$predictions,
                                                                                         test_sampnames = data$roc$inTest,
                                                                                         covars = mSet$dataSet$covars,
                                                                                         metadata_focus = input$ml_mistake_var,
                                                                                         cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                                                         smooth_line = input$ml_mistake_smooth,
                                                                                         show_reps = input$ml_mistake_per_rep#,show_reps = T
                         ) else NULL
                         list(ml_mistake = mistake_plot)
                       }
                      }else list()
                   },
                   multigroup = {
                     p = ggplotSummary(mSet, my_selection$mz, 
                                       shape.fac = input$shape_var, 
                                       cols = lcl$aes$mycols, 
                                       cf=gbl$functions$color.functions[[lcl$aes$spectrum]], 
                                       mode = "multi",
                                       styles = input$ggplot_sum_style,
                                       add_stats = input$ggplot_sum_stats, 
                                       color.fac = input$col_var, 
                                       text.fac = input$txt_var,
                                       fill.fac = input$fill_var)
                     list(summary_plot = p)
                   },
                   meba = {
                     p = ggPlotMeba(mSet, 
                                    #my_selection$mz,
                                    #draw.average = T,
                                    #cols = lcl$aes$mycols,
                                    cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                    topn=input$meba_topn)
                     list(meba_plot = p)
                   },
                   tt = {
                     # render manhattan-like plot for UI
                     p = ggPlotTT(mSet,
                                  cf = gbl$functions$color.functions[[lcl$aes$spectrum]], 
                                  20,topn=input$tt_topn)
                     
                     list(tt_plot = p)
                   },
                   fc = {
                     # render manhattan-like plot for UI
                     p <- ggPlotFC(mSet,
                                   gbl$functions$color.functions[[lcl$aes$spectrum]], 20, 
                                   topn=input$fc_topn)
                     
                     list(fc_plot = p)
                   },
                   diffcorr = {
                     
                     matr = mSet$analSet$diffcorr[,c("Gene1", "Gene2", "zScoreDiff")]
                     
                     matr=matr[1:100,]
                     colnames(matr) <- c("From", "To", "Weight")
                     mygraph <- igraph::graph.data.frame(matr)
                     
                     adjmatr = igraph::get.adjacency(mygraph, 
                                                     sparse = TRUE, 
                                                     attr='Weight')
                     
                     igr = igraph::graph.adjacency(adjmatrix = adjmatr,
                                                   weighted = T,
                                                   diag = F)
                     
                     netw = visNetwork::visIgraph(igr)
                     
                     cf = gbl$functions$color.functions[[lcl$aes$spectrum]]
                     
                     degr = igraph::degree(igr)
                     cols = cf(max(degr))
                     
                     netw$x$nodes$color <- sapply(netw$x$nodes$id, function(id){
                       cols[degr[names(degr) == id]]
                     })
                     
                     netw$x$nodes$title = netw$x$nodes$label
                     
                     netw$x$edges$color <- sapply(netw$x$edges$weight, function(w){
                       if(w < 0) "red" else "blue"
                     })
                     
                     netw$x$edges$title <- paste0("z-score=",netw$x$edges$weight)
                     
                     netw$x$edges$weight <- abs(netw$x$edges$weight)
                     netw$x$edges$value <- netw$x$edges$weight
                     
                     p = netw %>% 
                       visNetwork::visOptions(highlightNearest = TRUE,
                                              collapse = T,
                                              autoResize = T,
                                              nodesIdSelection = TRUE) %>% 
                       visNetwork::visNodes(borderWidth = 2,font = list(size=22,
                                                                        face=lcl$aes$font$family)) %>%
                       visNetwork::visEdges(scaling = list(min = 3,
                                                           max = 12)) %>%
                       visNetwork::visInteraction(#navigationButtons = TRUE,
                         keyboard = T)
                     if(!input$network_auto){
                       if(input$network_style == "hierarchical"){
                         p = p %>% visNetwork::visHierarchicalLayout()  
                       }else{
                         p = p %>% visNetwork::visIgraphLayout(layout = input$network_style)
                       }
                     }
                     # --- heatmap ---
                     hmap_matr = as.data.frame(as.matrix(adjmatr))
                     hmap_matr[hmap_matr==0 | lower.tri(hmap_matr)] <- NA
                     
                     p2 = heatmaply::heatmaply(hmap_matr,
                                               Colv = T,
                                               Rowv = T,
                                               branches_lwd = 0.3,
                                               margins = c(0, 0, 0, 0),
                                               col = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                               column_text_angle = 90,
                                               ylab = "m/z\n",
                                               showticklabels = if(ncol(hmap_matr) <= 100) c(T,T) else c(F,F),
                                               symm = T,
                                               symbreaks = T,
                                               dendrogram="none"
                     )
                     lcl$vectors$diffcorr_heatmap <- p2$x$layout$yaxis$ticktext
                     list(diffcorr = p, 
                          diffcorr_heatmap = p2)
                   },
                   network = {
                     
                     pval = input$network_sign
                     
                     mat = mSet$analSet$network$rcorr
                     matp = mat$P
                     matr = mat$r
                     matr[matp < pval] <- 0
                     
                     #palette = colorRampPalette(c("green", "white", "red")) (20)
                     #heatmap(x = matr, col = palette, symm = TRUE)
                     # #
                     # justMZ = matr[,"240.04920"]
                     # #
                     
                     igr = igraph::graph.adjacency(adjmatrix = matr,
                                                   weighted = T,
                                                   diag = F)
                     netw = visNetwork::visIgraph(igr)
                     
                     cf = gbl$functions$color.functions[[lcl$aes$spectrum]]
                     degr = igraph::degree(igr)
                     cols = cf(max(degr))
                     
                     netw$x$nodes$color <- sapply(netw$x$nodes$id, function(id){
                       cols[degr[names(degr) == id]]
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
                     
                     mat = mSet$analSet$heatmap$matrix[1:if(input$heatmap_topn < nrow(mSet$analSet$heatmap$matrix)) input$heatmap_topn else nrow(mSet$analSet$heatmap$matrix),]
                     
                     sideLabels = if(mSet$settings$exp.type %in% c("2f", "t1f")){
                       as.data.frame(mSet$analSet$heatmap$translator[,!1])
                     }else{
                       varOrder = match(mSet$dataSet$covars$sample, colnames(mat))
                       as.data.frame(mSet$dataSet$covars[, input$fill_var, with=F])
                     }
                     
                     sidePalette = if(mSet$settings$exp.type %in% c("2f", "t1f")){
                       mSet$analSet$heatmap$colors
                     }else{
                       hmap.lvls = unlist(unique(sideLabels))
                       cols = if(length(hmap.lvls) < length(lcl$aes$mycols)) lcl$aes$mycols else gbl$functions$color.functions[[lcl$aes$spectrum]](length(hmap.lvls))
                       color.mapper <- {
                         classes <- hmap.lvls
                         cols <- sapply(1:length(classes), function(i) cols[i]) # use user-defined colours
                         names(cols) <- classes
                         # - - -
                         cols
                       }
                       color.mapper
                     }
                     
                     p1 = {
                       
                       if(!is.null(mSet$analSet$heatmap$matrix)){
                         if(input$heatmap_topn > 2000){
                           data = data.frame(text = "Huge heatmap!\nDeactivating interactivity to avoid crashing.\nPlease apply a filter.")
                           hmap_int <- ggplot2::ggplot(data) + ggplot2::geom_text(ggplot2::aes(label = text), x = 0.5, y = 0.5, size = 10) +
                             ggplot2::theme(text = ggplot2::element_text(family = lcl$aes$font$family)) + ggplot2::theme_bw()
                          hmap_int
                          }else{
                           hmap_int <- suppressWarnings({
                             if(input$heatlimits){
                               heatmaply::heatmaply(mat,
                                                    Colv = mSet$analSet$heatmap$my_order,
                                                    Rowv = T,
                                                    branches_lwd = 0.3,
                                                    margins = c(60, 0, NA, 50),
                                                    col = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                    col_side_colors = sideLabels,
                                                    col_side_palette = sidePalette,
                                                    subplot_widths = c(.9,.1),
                                                    subplot_heights = if(mSet$analSet$heatmap$my_order) c(.1, .05, .85) else c(.05,.95),
                                                    column_text_angle = 90,
                                                    xlab = "Sample",
                                                    ylab = "m/z",
                                                    showticklabels = c(if(ncol(mat) <= 31) T else F, if(nrow(mat) <= 31) T else F),
                                                    limits = c(min(mSet$dataSet$norm), max(mSet$dataSet$norm)),
                                                    symbreaks = T
                               )
                             }else{
                               heatmaply::heatmaply(mat,
                                                    Colv = mSet$analSet$heatmap$my_order,
                                                    Rowv = T,
                                                    branches_lwd = 0.3,
                                                    margins = c(60, 0, NA, 50),
                                                    colors = gbl$functions$color.functions[[lcl$aes$spectrum]](256),
                                                    col_side_colors = sideLabels,
                                                    col_side_palette = sidePalette,
                                                    subplot_widths = c(.9,.1),
                                                    subplot_heights = if(mSet$analSet$heatmap$my_order) c(.1, .05, .85) else c(.05,.95),
                                                    column_text_angle = 90,
                                                    xlab = "Sample",
                                                    ylab = "m/z",
                                                    showticklabels = c(if(ncol(mat) <= 31) T else F, if(nrow(mat) <= 31) T else F),
                                                    symbreaks=T
                               )
                             }
                           })
                           # create heatmap object
                           hmap_int$x$layout$annotations[[1]]$text <- ""
                           # save the order of mzs for later clicking functionality
                           lcl$vectors$heatmap <- hmap_int$x$layout[[if(mSet$settings$exp.type %in% c("2f", "t", "t1f")) "yaxis2" else "yaxis3"]]$ticktext 
                           
                           }
                         # return
                         hmap_int
                       }else{
                         data = data.frame(text = "No significant hits available!\nPlease try alternative source statistics below.")
                         ggplot2::ggplot(data) + ggplot2::geom_text(ggplot2::aes(label = text), x = 0.5, y = 0.5, size = 10) +
                           ggplot2::theme(text = ggplot2::element_text(family = lcl$aes$font$family)) + ggplot2::theme_bw()
                       }
                     }
                     
                     p2 = 
                       function(){
                         if(!is.null(mSet$analSet$heatmap$matrix)){
                           suppressWarnings({
                             colSide = sideLabels
                             colnames(colSide) <- ""
                             pal = sidePalette
                             samples = as.character(colSide[[1]])
                             sampCols = as.character(sapply(samples, function(x) pal[names(pal) == x]))
                             hmap_stat = heatmap3::heatmap3(mat,
                                                            Colv = mSet$analSet$heatmap$my_order,
                                                            Rowv = T,
                                                            col = gbl$functions$color.functions[[lcl$aes$spectrum]](256),
                                                            ColSideColors = sampCols,
                                                            xlab = "Sample",
                                                            ColSideLabs = colnames(sideLabels),
                                                            ylab = "m/z",
                                                            balanceColor = if(input$heatlimits) T else F,
                                                            labRow = if(nrow(mat) <= 31) rownames(mat) else rep("", nrow(mat)),
                                                            labCol = if(ncol(mat) <= 31) colnames(mat) else rep("", ncol(mat)),
                                                            scale = "none",
                                                            legendfun = function() plot(0, xaxt = "n", bty = "n", yaxt = "n", type = "n", xlab = "", 
                                                                                        ylab = "")
                             )
                           })
                           # create heatmap object
                           # save the order of mzs for later clicking functionality
                           #lcl$vectors$heatmap <- c()
                           # return
                           hmap_stat
                         }else{
                           data = data.frame(text = "No significant hits available!\nPlease try alternative source statistics below.")
                           ggplot2::ggplot(data) + ggplot2::geom_text(ggplot2::aes(label = text), x = 0.5, y = 0.5, size = 10) +
                             ggplot2::theme(text = ggplot2::element_text(family = lcl$aes$font$family)) + ggplot2::theme_bw()
                         }
                       }
                     
                     list(heatmap_interactive = p1,
                          heatmap_static = p2)
                   }, 
                   power = {
                     p = {
                       if("power" %in% names(mSet$analSet)){
                         ggPlotPower(mSet, 
                                     max_samples = max(mSet$analSet$power[[1]]$Jpred),
                                     comparisons = names(mSet$analSet$power),
                                     cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       }else{
                         NULL
                       }
                     }
                     list(power_plot = p)
                   },
                   combi = {
                     p = {
                       if("combi" %in% names(mSet$analSet)){
                         ggPlotCombi(mSet,
                                     cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       }else{
                         NULL
                      }
                     }
                     list(combi_plot = p)
                   },
                   wordcloud = {
                     if(nrow(lcl$tables$wordcloud_filt) > 0){
                       topWords = if(input$wordcloud_topWords > nrow(lcl$tables$wordcloud_filt)) nrow(lcl$tables$wordcloud_filt) else input$wordcloud_topWords
                       wordcloud = wordcloud2::wordcloud2(data.table::as.data.table(lcl$tables$wordcloud_filt)[order(n, decreasing = T)][1:topWords,], color = "random-light", size=.7, shape = "circle")
                       wordbar = {
                         wcdata = data.table::as.data.table(lcl$tables$wordcloud_filt)[order(n, decreasing = T)][1:topWords,]
                         colnames(wcdata)[2] <- "freq"
                         ggPlotWordBar(wcdata = wcdata,
                                       cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       }
                     }
                     list(
                       wordcloud = wordcloud, 
                       wordbar = wordbar)
                   }
  )
  
  finalPlots <- mapply(function(myplot, plotName){
    
    targets = "aov|tt|fc|corr|asca|volcano|meba"
    
    if(grepl(targets, plotName)){
      
      whichAnal <- stringr::str_match(plotName, targets)[,1]
      
      if(whichAnal == "aov"){
        whichAnal = if(mSet$settings$exp.type %in% c("t", "2f", "t1f")) "aov2" else "aov"
      }
      if(whichAnal == "meba") whichAnal <- "MB"
      
      if(is.null(mSet$analSet[[whichAnal]][[if(whichAnal == "corr") "cor.mat" else if(whichAnal == "asca") "sig.list" else if(whichAnal == "MB") "stats" else "sig.mat"]])){
        data = data.frame(text = "No significant hits!")
        myplot = ggplot2::ggplot(data) + ggplot2::geom_text(ggplot2::aes(label = text), x = 0.5, y = 0.5, size = 10)
      }
    }
    
    isSquare <- grepl("pca|plsda|tsne|roc|heatmap|var|samp|network|umap|ica", plotName) & !grepl("scree|cv|perm|venn", plotName)
    
    # === WRAPPER ===
    
    canBe3D <- grepl("pca|plsda|tsne|umap|ica", plotName) & !grepl("scree|perm|cv", plotName)
    if(canBe3D){
      whichAnal <- stringr::str_match(plotName, "pca|plsda|tsne|umap|ica")[,1]
      is3D <- !input[[paste0(whichAnal, "_2d3d")]]
    }else{
      is3D <- plotName %in% c("network")
    }
    
    if(!is3D & !grepl("heatmap", plotName)){
      myplot <- myplot + guides(fill = guide_legend(ncol = 1),
                                shape = guide_legend(ncol = 1),
                                color = guide_legend(ncol = 1))
   
      myplot <- myplot + 
        gbl$functions$plot.themes[[lcl$aes$theme]](base_size = lcl$aes$font$plot.font.size) + 
        ggplot2::theme(legend.position = if(input$legend) "right" else "none",
                       legend.key.size = unit(.5,"line"),
                       legend.title = element_text(size=15),
                       legend.text = element_text(size=12),
                       axis.line = ggplot2::element_line(colour = 'black',
                                                         size = .5),
                       plot.title = ggplot2::element_text(hjust = 0.5,
                                                          vjust = 0.1,
                                                          size=lcl$aes$font$plot.font.size * 1.2),
                       text = ggplot2::element_text(family = lcl$aes$font$family))
      
      if(grepl("venn", plotName) & !input$venn_plot_mode){
        myplot <- myplot +
          ggplot2::theme_void(base_size = lcl$aes$font$plot.font.size) +
          ggplot2::theme(panel.grid = ggplot2::element_blank(),
                         legend.position="none",
                         text = ggplot2::element_text(family = lcl$aes$font$family))
        
      }
      
      
      if(grepl("ml_bar", plotName)){
        myplot <- myplot + 
          ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                         axis.ticks.x=ggplot2::element_blank(),
                         text = ggplot2::element_text(family = lcl$aes$font$family))
      }
      
      if(length(myplot$data) > 0){
        data = myplot$data
      }else{
        data = myplot[["layers"]][[1]][["data"]]
      }
      
      if(any(mSet$report$mzStarred$star) & any(grepl("mz|m/z", names(data))) & !grepl("ml", plotName)){
        #if(grepl("tt|fc|volc|pca_load|plsda_load|ml_bar", plotName)){
        myX = rlang::quo_get_expr(myplot[["layers"]][[1]][["mapping"]][['x']])
        myY = rlang::quo_get_expr(myplot[["layers"]][[1]][["mapping"]][['y']])
        myText = rlang::quo_get_expr(myplot[["layers"]][[1]][["mapping"]][['text']])
        myCol=rlang::quo_get_expr(myplot[["layers"]][[1]][["mapping"]][['colour']])
        flip = grepl("tt|fc|aov|var|samp|corr", plotName)
        matchMe = match(data[[myText]], mSet$report$mzStarred[star == TRUE]$mz)
        isMatch = which(!is.na(matchMe))
        xVals = data[[myX]][isMatch]
        yVals = data[[myY]][isMatch]
        labelVals = data[[myText]][isMatch]
        colVals = data[[myCol]][isMatch]
        
        data = data.frame(x = xVals,
                          y = yVals,
                          text = labelVals,
                          symb = c("\u2605"),
                          col = colVals)
        
        if(is.numeric(xVals[1]) & is.numeric(yVals[1])){
          myplot <- 
            myplot + 
            ggplot2::geom_text(data = data,
                               aes(x = x,
                                   y = y,
                                   label = symb),
                               color = "black",
                               show.legend = F,
                               size=7) + 
            ggplot2::geom_text(data = data,
                               aes(x = x,
                                   y = y,
                                   color = col,
                                   label = symb),
                               show.legend = F,
                               size=5) +
            # ggrepel::geom_text_repel(data = data,
            #                          aes(x = x,# + .04*min(x),
            #                              y = y,#+ .04*min(y),
            #                              label = text),point.padding = 1,
            #                          size=4) 
          ggplot2::geom_text(data = data,
                             aes(x = x, #+ 0.04 * max(x),
                                 y = y,#+ 0.04 * max(x),
                                 label = text),
                             position = 	
                               position_jitter(),
                             size=4)
        }
      }
    }
    finalPlot = list(myplot)
    finalPlot
  }, toWrap, names(toWrap))
  res = list(lcl = lcl, plots = finalPlots)
  res
}

metshiProcess <- function(mSet, session, init=F, cl=0){
  
  sums = colSums(mSet$dataSet$missing)
  missing.per.mz.perc = sums/nrow(mSet$dataSet$missing)*100
  good.inx <- missing.per.mz.perc < mSet$metshiParams$miss_perc
  mSet$dataSet$orig <- as.data.frame(mSet$dataSet$orig[, good.inx, drop = FALSE])
  
  qs::qsave(mSet$dataSet$orig, "data_orig.qs")
  
  if(!init) mSet$dataSet$missing <- NULL
  
  if(mSet$metshiParams$filt_type != "none" & (ncol(mSet$dataSet$orig) > mSet$metshiParams$max.allow)){
    
    # TODO; add option to only keep columns that are also in QC ('qcfilter'?)
    keep.mz <- colnames(FilterVariableMetshi(mSet,
                                             filter = mSet$metshiParams$filt_type,
                                             qcFilter = "F", #TODO: mSet$metshiParams$useQCs
                                             rsd = 25,
                                             max.allow = mSet$metshiParams$max.allow
    )$dataSet$filt)  
    if(mSet$metshiParams$norm_type == "ProbNorm"){
      keep.mz = unique(c(keep.mz, mSet$metshiParams$ref_var))
    }
    mSet$dataSet$orig <- mSet$dataSet$orig[,keep.mz]
    mSet$dataSet$filt <- NULL
  }
  
  qs::qsave(mSet$dataSet$orig, "data_orig.qs")
  
  # sanity check data
  mSet <- MetaboAnalystR::SanityCheckData(mSet)
  
  #shiny::setProgress(session=session, value= .6)
  
  # missing value imputation
  if(req(mSet$metshiParams$miss_type) != "none"){
    if(req(mSet$metshiParams$miss_type) == "rowmin"){ # use sample minimum
      mSet <- replRowMin(mSet)
    }
    else if(req(mSet$metshiParams$miss_type ) == "pmm"){ # use predictive mean matching
      # TODO: re-enable, it's very slow
      base <- mSet$dataSet$orig
      imp <- mice::mice(base, printFlag = TRUE)
      
    }else if(req(mSet$metshiParams$miss_type ) == "rf"){ # random forest
      mSet$dataSet$proc <- replRF(mSet, 
                                               parallelMode = mSet$metshiParams$rf_norm_parallelize, 
                                               ntree = mSet$metshiParams$rf_norm_ntree,
                                               cl = cl,
                                               rf.method = mSet$metshiParams$rf_norm_method)
      w.missing <- qs::qread("preproc.qs")
      rownames(mSet$dataSet$proc) <- rownames(w.missing)
      # - - - - - - - - - - - -
    }else{
      # use built in imputation methods, knn means etc.
      mSet <- MetaboAnalystR::ImputeMissingVar(mSet,
                                               method = mSet$metshiParams$miss_type
      )
    }
  }
  
  # if normalizing by a factor, do the below
  if(req(mSet$metshiParams$norm_type) == "SpecNorm"){
    rematch = match(
      rownames(mSet$dataSet$preproc),
      mSet$dataSet$covars$sample
    )
    mSet$dataSet$covars <- mSet$dataSet$covars[rematch,]
    
    norm.vec <<- mSet$dataSet$covars[[mSet$metshiParams$samp_var]]
    norm.vec <<- scale(x = norm.vec, center = 1)[,1] # normalize scaling factor
  }else{
    norm.vec <<- rep(1, length(mSet$dataSet$cls)) # empty
  }
  
  mSet <- MetaboAnalystR::PreparePrenormData(mSet)

  if(mSet$metshiParams$norm_type == "QcNorm"){
    data <- qs::qread("prenorm.qs")  
    rematch = match(
      rownames(data),
      mSet$dataSet$covars$sample
    )
    mSet$dataSet$covars <- mSet$dataSet$covars[rematch,]
    
    batches = mSet$dataSet$covars$batch
    normalized_blocks = pbapply::pblapply(unique(batches),function(lvl){
      rows = data[which(batches == lvl),]
      is_qc = grep("^qc", tolower(rownames(rows)))
      if(length(is_qc) == nrow(rows)){
        rows
      }else{
        qcs = rows[is_qc,]
        avg_qc_sample = colMeans(qcs)
        non_qcs = rows[-is_qc,]
        qc_norm_rows = lapply(1:nrow(non_qcs), function(i){
          x = non_qcs[i,]
          as.list(x/median(as.numeric(x/avg_qc_sample), na.rm = T))
          #as.list(non_qcs[i,]/avg_qc_sample)
        })
        res = as.data.frame(data.table::rbindlist(qc_norm_rows, use.names = T))
        rownames(res) = rownames(non_qcs)
        rbind(qcs, res)  
      }
    })
    qc_norm_table = do.call("rbind", normalized_blocks)
    mSet$dataSet$norm <- qc_norm_table
  }else{
    data <- qs::qread("prenorm.qs")  
    # normalize dataset with user settings(result: mSet$dataSet$norm)
    mSet <- MetaboAnalystR::Normalization(mSet,
                                          rowNorm = mSet$metshiParams$norm_type,
                                          transNorm = mSet$metshiParams$trans_type,
                                          scaleNorm = mSet$metshiParams$scale_type,
                                          ref = mSet$metshiParams$ref_var) 

  }
  
  mSet$dataSet$prenorm <- NULL
  
  #shiny::setProgress(session=session, value= .8)
  
  # get sample names
  smps <- rownames(mSet$dataSet$norm)
  # get which rows are QC samples
  qc_rows <- which(grepl(pattern = "QC", x = smps))
  #print(qc_rows)
  # if at least one row has a QC in it, batch correct
  has.qc <- length(qc_rows) > 0
  
  rematch = match(
    rownames(mSet$dataSet$norm),
    mSet$dataSet$covars$sample
  )
  
  mSet$dataSet$covars <- mSet$dataSet$covars[rematch,]
  
  # lowercase all the covars table column names
  colnames(mSet$dataSet$covars) <- tolower(colnames(mSet$dataSet$covars))
  
  mSet$dataSet$prebatch <- mSet$dataSet$norm
  
  left_batch_vars = mSet$metshiParams$batch_var
  
  batch_method_a = mSet$metshiParams$batch_method_a
  batch_method_b = mSet$metshiParams$batch_method_b
  
  # IN CASE SUBSETTING ELIMINATES BATCH EFFECT (ONLY ONE BATCH SELECTED)
  # keep var 1 ?
  keep.batch.1 = length(unique(unlist(mSet$dataSet$covars[, left_batch_vars[1], with=F]))) > 1
  # keep var 2 ? 
  keep.batch.2 = if(length(left_batch_vars) > 1) length(unique(unlist(mSet$dataSet$covars[, left_batch_vars[2], with=F]))) > 1 else F
  if(!keep.batch.1 & !keep.batch.2){
    left_batch_vars = c()
  }else if(!keep.batch.1 & keep.batch.2){
    left_batch_vars = left_batch_vars[2]
    batch_method_a = batch_method_b 
  }else if(keep.batch.1 & !keep.batch.2){
    left_batch_vars = left_batch_vars[1]
  }
  
  if(length(left_batch_vars)>0){
    
    # APPLY THE FIRST METHOD ONLY FOR BATCH + INJECTION
    
    if(batch_method_a == "limma" & 
       batch_method_b == "limma" & 
       length(left_batch_vars) == 2){
      # create a model table
      csv_pheno <- data.frame(sample = 1:nrow(mSet$dataSet$covars),
                              batch1 = mSet$dataSet$covars[, left_batch_vars[1], with=FALSE][[1]],
                              batch2 = mSet$dataSet$covars[, left_batch_vars[2], with=FALSE][[1]]
                              #,outcome = as.factor(exp_lbl)
      )
      
      csv_edata <- combatCSV(mSet, tbl = "norm")
      
      # batch correct with limma and two batches
      batch_normalized = t(limma::removeBatchEffect(x = csv_edata,
                                                    batch = csv_pheno$batch1
                                                    ,batch2 = csv_pheno$batch2))
      rownames(batch_normalized) <- rownames(mSet$dataSet$norm)
      mSet$dataSet$norm <- as.data.frame(batch_normalized)
    }else{
      if("batch" %in% left_batch_vars){# & mSet$metshiParams$batch_use_qcs){# & has.qc){
        
        # get batch for each sample
        batch.idx = as.numeric(as.factor(mSet$dataSet$covars$batch))
        
        if(length(batch.idx) == 0) return(mSet$dataSet$norm)
        # get injection order for samples
        hasRT = any(grepl(pattern = "RT", colnames(mSet$dataSet$proc)))
        
        if(hasRT & batch_method_a == "batchCorr"){
          metshiAlert("Only available for LC-MS data! Defaulting to WaveICA.")
          batch_method_a <- "waveica"
        }
        
        mSet$dataSet$norm <- batchCorr_mSet(mSet, batch_method_a, batch_var = left_batch_vars, cl=cl)
        
        left_batch_vars <- grep(left_batch_vars,
                                pattern = "batch|injection|sample",
                                value = T,
                                invert = T)
      }
      
      # check which batch values are left after initial correction
      if(length(left_batch_vars) == 0){
        NULL # if none left, continue after this
      } else{
        mSet$dataSet$norm <- batchCorr_mSet(mSet, batch_method_b, batch_var = left_batch_vars, cl=cl) 
      }}
  }
  
  if(mSet$metshiParams$pca_corr){
    print("Performing PCA and subtracting PCs...")
    res <- prcomp(mSet$dataSet$norm, 
                  center = F,
                  scale = F)
    pc.use <- as.numeric(mSet$metshiParams$keep_pcs[1]:mSet$metshiParams$keep_pcs[2]) # explains 93% of variance
    trunc <- res$x[,pc.use] %*% t(res$rotation[,pc.use])
    mSet$dataSet$norm <- as.data.frame(trunc)
  }
  
  #shiny::setProgress(session=session, value= .9)
  
  mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
  
  # make sure covars order is consistent with mset$..$norm order
  rematch = match(
    rownames(mSet$dataSet$norm),
    mSet$dataSet$covars$sample
  )
  
  mSet$dataSet$covars <- mSet$dataSet$covars[rematch,]
  
  mSet$report <- list(mzStarred = data.table::data.table(mz = colnames(mSet$dataSet$norm),
                                                         star = c(FALSE)))  
  data.table::setkey(mSet$report$mzStarred, mz)
  
  if(has.qc & !init){
    mSet <- hideQC(mSet)
  }
  
  rematch = match(
    rownames(mSet$dataSet$norm),
    mSet$dataSet$covars$sample
  )
  
  mSet$dataSet$covars <- mSet$dataSet$covars[rematch,]
  
  if(!init){
    mSet$dataSet$missing <- mSet$dataSet$start <- NULL 
  }
  
  mSet$analSet <- list(type = "stat")
  mSet
  # TODO: WHERE IS THE REORDERING FAILING?
}

render.kegg.node.jw <- function (plot.data, cols.ts, img, same.layer = TRUE, type = c("gene", 
                                                                                      "compound")[1], text.col = "black", cex = 0.25) 
{
  width = ncol(img)
  height = nrow(img)
  nn = nrow(plot.data)
  pwids = plot.data$width
  if (!all(pwids == max(pwids))) {
    message("Info: ", "some node width is different from others, and hence adjusted!")
    wc = table(pwids)
    pwids = plot.data$width = as.numeric(names(wc)[which.max(wc)])
  }
  if (type == "gene") {
    if (same.layer != T) {
      rect.out = pathview:::sliced.shapes(plot.data$x + 0.5, height - 
                                            plot.data$y, plot.data$width/2 - 0.5, plot.data$height/2 - 
                                            0.25, cols = cols.ts, draw.border = F, shape = "rectangle")
      text(plot.data$x + 0.5, height - plot.data$y, labels = as.character(plot.data$labels), 
           cex = cex, col = text.col)
      return(invisible(1))
    }
    else {
      img2 = img
      pidx = cbind(ceiling(plot.data$x - plot.data$width/2) + 
                     1, floor(plot.data$x + plot.data$width/2) + 
                     1, ceiling(plot.data$y - plot.data$height/2) + 
                     1, floor(plot.data$y + plot.data$height/2) + 
                     1)
      cols.ts = cbind(cols.ts)
      ns = ncol(cols.ts)
      brk.x = sapply(plot.data$width/2, function(wi) seq(-wi, 
                                                         wi, length = ns + 1))
      for (k in 1:ns) {
        col.rgb = col2rgb(cols.ts[, k])/255
        pxr = t(apply(pidx[, 1:2], 1, function(x) x[1]:x[2])) - 
          plot.data$x - 1
        sel = pxr >= ceiling(brk.x[k, ]) & pxr <= floor(brk.x[k + 
                                                                1, ])
        for (i in 1:nn) {
          sel.px = (pidx[i, 1]:pidx[i, 2])[sel[i, ]]
          node.rgb = img[pidx[i, 3]:pidx[i, 4], sel.px, 
                         1:3]
          node.rgb.sum = apply(node.rgb, c(1, 2), sum)
          blk.ind = which(node.rgb.sum == 0 | node.rgb.sum == 
                            1, arr.ind = T)
          node.rgb = array(col.rgb[, i], dim(node.rgb)[3:1])
          node.rgb = aperm(node.rgb, 3:1)
          for (j in 1:3) node.rgb[cbind(blk.ind, j)] = 0
          img2[pidx[i, 3]:pidx[i, 4], sel.px, 1:3] = node.rgb
        }
      }
      return(img2)
    }
  }
  else if (type == "compound") {
    plot.data <<- plot.data
    if (same.layer != T) {
      nc.cols = ncol(cbind(cols.ts))
      if (nc.cols > 2) {
        na.cols = rep("#FFFFFF", nrow(plot.data))
        cir.out = pathview:::sliced.shapes(plot.data$x, height - 
                                             plot.data$y, plot.data$width[1], plot.data$width[1], 
                                           cols = na.cols, draw.border = F, shape = "ellipse", 
                                           lwd = 0.2)
      }
      cir.out = pathview:::sliced.shapes(plot.data$x, height - plot.data$y, 
                                         plot.data$width[1], plot.data$width[1], cols = cols.ts, 
                                         shape = "ellipse", blwd = 0.2)
      return(invisible(1))
    }
    else {
      blk = c(0, 0, 0)
      img2 = img
      
      w = ncol(img)
      h = nrow(img)
      
      cidx = rep(1:w, each = h)
      ridx = rep(1:h, w)
      
      plot.data <<- plot.data
      
      pidx = lapply(1:nn, function(i) {
        ii = which((cidx - plot.data$x[i])^2 + (ridx - 
                                                  plot.data$y[i])^2 < (plot.data$width[i])^2)
        imat = cbind(cbind(ridx, cidx)[rep(ii, each = 3), 
        ], 1:3)
        imat[, 1:2] = imat[, 1:2] + 1
        ib = which(abs((cidx - plot.data$x[i])^2 + (ridx - 
                                                      plot.data$y[i])^2 - (plot.data$width[i])^2) <= 
                     8)
        ibmat = cbind(cbind(ridx, cidx)[rep(ib, each = 3), 
        ], 1:3)
        ibmat[, 1:2] = ibmat[, 1:2] + 1
        return(list(fill = imat, border = ibmat))
      })
      cols.ts = cbind(cols.ts)
      
      rows_corr_list = lapply(1:nrow(cols.ts), function(i){
        row = cols.ts[i,]
        has.hit = !(all(row == "#FFFFFF"))
        if(has.hit){
          difvals = unique(setdiff(row, "#FFFFFF"))
          uniq.hits = length(difvals)
          rep(difvals, each = floor(length(row)/uniq.hits))#, length.out = length(row)) 
        }else{
          row
        }
      })
      cols.ts = do.call("rbind", rows_corr_list)
      rows_corr_list <<- rows_corr_list
      
      ns = ncol(cols.ts)
      brk.x = sapply(plot.data$width, function(wi) seq(-wi, 
                                                       wi, length = ns + 1))
      for (i in 1:nn) {
        pxr = pidx[[i]]$fill[, 2] - 1 - plot.data$x[i]
        cols.ts.adj = cols.ts[i, ]
        # if(!all(cols.ts.adj == "#FFFFFF")){
        #   cols.ts.adj = cols.ts.adj[cols.ts[i, ] != "#FFFFFF"]
        # }
        # ns = length(cols.ts.adj)
        # print(ns)
        col.rgb = col2rgb(cols.ts.adj)/255
        for (k in 1:ns) {
          sel = pxr >= brk.x[k, i] & pxr <= brk.x[k + 
                                                    1, i]
          img2[pidx[[i]]$fill[sel, ]] = col.rgb[, k]
        }
        img2[pidx[[i]]$border] = blk
      }
      print("uwu")
      return(img2)
    }
  }
  else stop("unrecognized node type!")
}

ml_loop_wrapper <- function(mSet_loc, gbl, jobs, ml_session_cl=0){
  
  parallel::clusterExport(ml_session_cl, c("ml_run",
                                           "gbl", 
                                           "mSet_loc"),
                          envir = environment())
  
  parallel::clusterEvalQ(ml_session_cl,{
    small_mSet <- qs::qread(mSet_loc)
  })
  
  pbapply::pblapply(jobs, 
                    cl = if(length(jobs) > 1) ml_session_cl else 0, 
                    function(settings, ml_cl){
                      res = list()
                      try({
                        res = ml_run(settings = settings, 
                                     mSet = small_mSet,
                                     input = input,
                                     cl = ml_cl)  
                      })
                      res
                    },
                    ml_cl = if(length(ml_queue$jobs) > 1) 0 else ml_session_cl)
}

assignInNamespace(x = "render.kegg.node", value = render.kegg.node.jw, ns = "pathview")

