# create listener
statsmanager <- shiny::reactiveValues()

shiny::observe({
  
  if(is.null(statsmanager$calculate)){
    
    NULL # if not reloading anything, nevermind
    
  }else{

    if(!is.null(mSet)){
      
      mSet.old <- mSet
      success = F
      
      try({
        switch(statsmanager$calculate,
               vennrich = {
                 if("storage" %not in% names(mSet)){
                   mSet$storage <- list()
                 }
                 # TODO: use this in venn diagram creation
                 mSet <- MetaboShiny::store.mSet(mSet)
               },
               pattern = {
                 # pearson kendall spearman
                 lvls = levels(mSet$dataSet$cls)
                 pat = input$pattern_seq_order
                 pat_order = match(lvls,pat)
                 pattern = paste0(pat_order-1, collapse="-")
                 mSet <- MetaboAnalystR::Match.Pattern(mSet, input$pattern_corr, pattern)
                 names(mSet$analSet)[names(mSet$analSet) == "corr"] <- "pattern" 
               },
               pca = {
                 shiny::withProgress({
                   mSet <- MetaboAnalystR::PCA.Anal(mSet) # perform PCA analysis
                 })
                 success = T
               },
               meba = {
                 shiny::withProgress({
                   mSet <- MetaboAnalystR::performMB(mSet, 10) # perform MEBA analysis
                 })
               },
               asca = {
                 # perform asca analysis
                 shiny::withProgress({
                   mSet <- MetaboAnalystR::Perform.ASCA(mSet, 1, 1, 2, 2)
                   mSet <- MetaboAnalystR::CalculateImpVarCutoff(mSet, 0.05, 0.9)
                 })
               },
               enrich = {
                 shiny::withProgress({
                   
                   enr_mSet <- MetaboAnalystR::InitDataObjects("mass_all", 
                                                               "mummichog",
                                                               FALSE)
                   
                   MetaboAnalystR::SetPeakFormat("rmp")
                   
                   enr_mSet <- MetaboAnalystR::UpdateInstrumentParameters(enr_mSet, 
                                                                          mSet$ppm, 
                                                                          "mixed");
                   
                   #similarly to venn diagram
                   flattened <- getTopHits(mSet, 
                                           input$mummi_anal,
                                           input$mummi_topn)
                   
                   setProgress(0.1)
                   
                   myFile <- tempfile(fileext = ".csv")
                   tbl = data.table::data.table("m.z" = as.numeric(gsub(flattened[[1]], pattern="+|-", replacement="")),
                                                mode = sapply(flattened[[1]], function(mz){
                                                  if(grepl(pattern="-",x=mz)) "negative" else "positive"
                                                }))
                   
                   tbl[, "p.value"] = c(0)
                   tbl[, "t.score"] = c(0)
                   
                   enr_mSet$dataSet$mummi.orig <- cbind(tbl$p.value,
                                                        tbl$m.z,
                                                        tbl$t.score)
                   colnames(enr_mSet$dataSet$mummi.orig) = c("p.value", 
                                                             "m.z", 
                                                             "t.score")
                   enr_mSet$dataSet$pos_inx <- tbl$mode == "positive"
                   enr_mSet$dataSet$mumType <- "list"
                   
                   mumDataContainsPval <<- 0#if(grepl(input$mummi_anal, pattern = "tt|aov|aov2")) 1 else 0
                   
                   shiny::setProgress(0.2)
                   
                   enr_mSet <- MetaboAnalystR::SanityCheckMummichogData(enr_mSet)
                   enr_mSet <- MetaboAnalystR::Setup.AdductData(enr_mSet, input$mummi_adducts);
                   
                   shiny::setProgress(0.3)
                   
                   require(enviPat)
                   
                   elecMass = 0.000548579909
                   mummi_adducts <- adducts[Name %in% input$mummi_adducts]
                   add_db_custom_rows <- lapply(1:nrow(mummi_adducts), function(i){
                     row = mummi_adducts[i,]
                     if(is.na(row$AddAt)) row$AddAt <- ""
                     if(is.na(row$AddEx)) row$AddEx <- ""
                     if(is.na(row$RemAt)) row$RemAt <- ""
                     if(is.na(row$RemEx)) row$RemEx <- ""
                     addForm <- mergeform(row$AddAt, row$AddEx)
                     remForm <- mergeform(row$RemAt, row$RemEx)
                     addMass <- enviPat::check_chemform(isotopes, addForm)$monoisotopic_mass
                     remMass <- enviPat::check_chemform(isotopes, remForm)$monoisotopic_mass
                     if(addMass < 0) addMass <- 0
                     if(remMass < 0) remMass <- 0
                     netChange = (addMass - remMass + (row$Charge * -elecMass)) / abs(row$Charge)
                     mzSpec <- if(row$xM > 1) paste0("(",row$xM,"*mw)") else if(abs(row$Charge)>1) paste0("mw/",abs(row$Charge)) else "mw"
                     addOn <- if(netChange < 0){
                       paste("-", abs(netChange))
                     }else{
                       paste("+", abs(netChange))
                     }
                     data.table::data.table(Ion_Name = row$Name, Ion_Mass = paste(mzSpec, addOn))
                   })
                   add_db_cust <- data.table::rbindlist(add_db_custom_rows)
                   
                   shiny::setProgress(0.4)
                   
                   # === PerformAdductMapping ===
                   myAdds <- enr_mSet$dataSet$adduct.list
                   hit.inx <- match(tolower(myAdds), tolower(add_db_cust$Ion_Name))
                   match.values <- add_db_cust[hit.inx, ]
                   sel.add <- nrow(match.values)
                   if (sel.add > 0) {
                     shiny::showNotification(paste("A total of ", sel.add, 
                                                   " adducts were successfully selected!", sep = ""))
                   }else{
                     shiny::showNotification("No adducts were selected!")
                   }
                   
                   enr_mSet$add.map <- match.values
                   
                   shiny::setProgress(0.5)
                   
                   # === GETLIB ===
                   lib = input$mummi_org
                   libVersion <- "current"#input$mummi_db_ver
                   
                   filenm <- paste(lib, ".rds", sep = "")
                   biocyc <- grepl("biocyc", lib)
                   
                   if (!is.null(enr_mSet$curr.cust)) {
                     if (biocyc) {
                       user.curr <- enr_mSet$curr.map$BioCyc
                     }else {
                       user.curr <- enr_mSet$curr.map$KEGG
                     }
                     currency <<- user.curr
                     if (length(currency) > 0) {
                       shiny::showNotification("Currency metabolites were successfully uploaded!")
                     }else {
                       shiny::showNotification("Errors in currency metabolites uploading!")
                     }
                   }
                   if (libVersion == "old" && end.with(lib, "kegg")) {
                     mum.url <- paste("https://www.metaboanalyst.ca/resources/libs/mummichog/kegg_2018/", 
                                      filenm, sep = "")
                   }else {
                     mum.url <- paste("https://www.metaboanalyst.ca/resources/libs/mummichog/", 
                                      filenm, sep = "")
                   }
                   download.file(mum.url, destfile = filenm, method = "libcurl", 
                                 mode = "wb")
                   mummichog.lib <- readRDS(filenm)
                   
                   shiny::setProgress(0.6)
                   
                   prematchies <- MetaboShiny::get_prematches(if(biocyc) "metacyc" else "kegg", 
                                                              "source", 
                                                              lcl$paths$patdb)
                   
                   iden.vs.add <- unique(prematchies[,c("identifier", "adduct")])
                   iden.vs.add$i <- match(iden.vs.add$identifier,mummichog.lib$cpd.lib$id)
                   iden.vs.add <- data.table::as.data.table(iden.vs.add[complete.cases(iden.vs.add),])
                   iden.vs.add <- iden.vs.add[adduct %in% mummi_adducts$Name]
                   
                   shiny::setProgress(0.7)
                   
                   # ==== NEW_ADDUCT_MZLIST ===
                   
                   use.rules <<- input$mummi_rules
                   
                   new_adduct_mzlist <- function (mSetObj = NA, mw){
                     
                     mode <- mSetObj$dataSet$mode
                     ion.name <- mSetObj$add.map$Ion_Name
                     ion.mass <- mSetObj$add.map$Ion_Mass
                     mw_modified <- NULL
                     neg.ions <- mummi_adducts[Ion_mode == "negative"]$Name
                     ion.name.neg <- intersect(ion.name, neg.ions)
                     ion.mass.neg <- ion.mass[which(ion.name %in% neg.ions)]
                     ion.name.pos <- setdiff(ion.name, neg.ions)
                     ion.mass.pos <- ion.mass[which(ion.name %in% ion.name.pos)]
                     mass.list.neg <- as.list(ion.mass.neg)
                     mass.user.neg <- lapply(mass.list.neg, function(x) eval(parse(text = paste(gsub("PROTON", 
                                                                                                     1.007825, x)))))
                     mw_modified.neg <- do.call(cbind, mass.user.neg)
                     colnames(mw_modified.neg) <- ion.name.neg
                     mass.list.pos <- as.list(ion.mass.pos)
                     mass.user.pos <- lapply(mass.list.pos, function(x) eval(parse(text = paste(gsub("PROTON", 
                                                                                                     1.007825, x)))))
                     mw_modified.pos <- do.call(cbind, mass.user.pos)
                     colnames(mw_modified.pos) <- ion.name.pos
                     mw_modified <- list(mw_modified.neg, mw_modified.pos)
                     
                     if(use.rules){
                       print("rules")
                       mw_modified <- lapply(mw_modified, function(mw_adds){
                         for(i in 1:nrow(mw_adds)){
                           ok.adducts <- iden.vs.add[i,]$adduct
                           mw_adds[i, !(colnames(mw_adds) %in% ok.adducts)] <- 0
                         }
                         mw_adds
                       }) 
                     }
                     names(mw_modified) <- c("neg", "pos")
                     return(mw_modified)}
                   
                   shiny::setProgress(0.8)
                   
                   #unlockBinding("new_adduct_mzlist", as.environment("package:MetaboAnalystR"))
                   assignInNamespace("new_adduct_mzlist", new_adduct_mzlist, ns="MetaboAnalystR", 
                                     envir=as.environment("package:MetaboAnalystR"))
                   
                   enr_mSet <- MetaboAnalystR::SetPeakEnrichMethod(enr_mSet, if(input$mummi_enr_method) "mummichog" else "gsea")
                   
                   shiny::setProgress(0.9)
                   
                   enr_mSet <- MetaboAnalystR::SetMummichogPval(enr_mSet, 0.0)#input$mummi_pval)
                   
                   enr_mSet <- MetaboAnalystR::PerformPSEA(enr_mSet, 
                                                           lib, 
                                                           libVersion, 
                                                           100)#input$mummi_perm)
                   mSet$analSet$enrich <- enr_mSet
                   
                   # === PLOT ====
                   
                   output$enrich_plot <- plotly::renderPlotly({
                     # --- ggplot ---
                     MetaboShiny::ggPlotMummi(mSet$analSet$enrich, 
                                              if(input$mummi_enr_method) "mummichog" else "gsea",
                                              cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                              plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                              plotlyfy=TRUE,font = lcl$aes$font)
                   })
                   # render results table
                   enrich$overview <- mSet$analSet$enrich$mummi.resmat
                 })
               },
               ml = {
                 try({
                   shiny::setProgress(value = 0)
                   
                   # get base table to use for process
                   curr <- data.table::as.data.table(mSet$dataSet$proc)
                   
                   # replace NA's with zero
                   curr <- curr[,(1:ncol(curr)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]
                   
                   # conv to data frame
                   curr <- as.data.frame(curr)
                   rownames(curr) <- rownames(mSet$dataSet$proc)
                   
                   # find the qc rows and remove them
                   is.qc <- grepl("QC|qc", rownames(curr))
                   if(sum(is.qc) > 0){
                     curr <- curr[!is.qc,]
                   }
                   # reorder according to covars table (will be used soon)
                   
                   order <- match(rownames(curr), mSet$dataSet$covars$sample)
                   if("label" %in% colnames(mSet$dataSet$covars)){
                     config <- mSet$dataSet$covars[order, -"label"]
                   }else{
                     config <- mSet$dataSet$covars[order, ]
                   }
                   
                   config <- config[, input$ml_include_covars,with=F]# reorder so both halves match up later
                   
                   if(mSet$dataSet$exp.type %in% c("2f", "t1f")){
                     # just set to facA for now..
                     if(nrow(config)==0){
                       config <- data.frame(label=mSet$dataSet$facA)
                     }else{
                       config <- cbind(config, label=mSet$dataSet$facA) # add current experimental condition
                     }
                   }else{
                     if(nrow(config)==0){
                       config <- data.frame(label=mSet$dataSet$cls)
                     }else{
                       config <- cbind(config, 
                                       label=mSet$dataSet$cls) # add current experimental condition
                     }
                   }
                   
                   if(!input$ml_random_split){
                     try({
                       shiny::showNotification("Using same train/test split for all repeats...")
                     })
                     # make split for all repeats
                     train_idx = caret::createDataPartition(y = config$label, p = input$ml_train_perc/100, list = FALSE) # partition data in a balanced way (uses labels)
                     # add column to config for this split
                     config$split <- c("train")
                     config$split[train_idx] <- "test"
                     # set test_vec and train_vec c(split, train), c(split, test)
                     lcl$vectors$ml_train <<- c("split", "train")
                     lcl$vectors$ml_test <<- c("split", "test")
                   }else{
                     if(is.null(lcl$vectors$ml_train)){
                       lcl$vectors$ml_train <<- c("all", "all")
                     }
                     if(is.null(lcl$vectors$ml_test)){
                       lcl$vectors$ml_test <<- c("all", "all")
                     }
                     if(all(lcl$vectors$ml_test == lcl$vectors$ml_train)){
                       if(unique(lcl$vectors$ml_test) == "all"){
                         shiny::showNotification("No subset selected... continuing in normal non-subset mode")
                       }else{
                         MetaboShiny::metshiAlert("Cannot test on the training set!")
                         return(NULL)
                       }
                     }
                   }
                   
                   config <- data.table::as.data.table(config)
                   config <- config[,apply(!is.na(config), 2, any), with=FALSE]
                   
                   predictor = config$label
                   predict_idx <- which(colnames(config)== "label")
                   exact_matches <- which(unlist(lapply(config, function(col) all(col == predictor))))
                   remove = setdiff(exact_matches, predict_idx)
                   
                   # remove ones w/ every row being different(may be used to identify...)
                   #covariates <- lapply(1:ncol(config), function(i) as.factor(config[,..i][[1]]))
                   #names(covariates) <- colnames(config)
                   
                   # # remove ones with na present
                   has.na <- apply(config, MARGIN=2, FUN=function(x) any(is.na(x) | tolower(x) == "unknown"))
                   has.all.unique <- apply(config, MARGIN=2, FUN=function(x) length(unique(x)) == length(x))
                   remove = colnames(config)[which(has.na | has.all.unique)]
                   
                   #keep_configs <- which(names(config) == "label")
                   remove <- unique(c(remove, "sample",  
                                      "individual", 
                                      colnames(config)[caret::nearZeroVar(config)]))
                   
                   keep_configs <- which(!(colnames(config) %in% remove))
                   
                   try({
                     shiny::showNotification(paste0("Keeping non-mz variables after NA/unique filtering: ",
                                                    paste0(names(config)[keep_configs],collapse = ", ")))
                   })
                   
                   config <- config[,..keep_configs,with=F]
                   
                   # rename the variable of interest to 0-1-2 etc.
                   char.lbl <- as.character(config$label)
                   uniques <- unique(char.lbl)
                   uniques_new_name <- c(1:length(uniques))
                   names(uniques_new_name) = uniques
                   
                   remapped.lbl <- uniques_new_name[char.lbl]
                   
                   # join halves together, user variables and metabolite data
                   curr <- cbind(config, curr)
                   curr <- data.table::as.data.table(curr)
                   
                   # how many models will be built? user input
                   goes = as.numeric(input$ml_attempts)
                   
                   # identify which columns are metabolites and which are config/covars
                   configCols <- which(!(gsub(x = colnames(curr), pattern = "_T\\d", replacement="") %in% colnames(mSet$dataSet$norm)))
                   mzCols <- which(gsub(x = colnames(curr), pattern = "_T\\d", replacement="") %in% colnames(mSet$dataSet$norm))
                   
                   # make the covars factors and the metabolites numeric.
                   nums <- which(unlist(lapply(curr, is.numeric))) 
                   configCols <- setdiff(configCols, nums)
                   
                   curr[,(configCols):= lapply(.SD, function(x) as.factor(x)), .SDcols = configCols]
                   curr[,(mzCols):= lapply(.SD, function(x) as.numeric(x)), .SDcols = mzCols]
                   
                   # = tuning = 
                   require(caret)
                   
                   # all methods
                   caret.mdls <- caret::getModelInfo()
                   caret.methods <- names(caret.mdls)
                   tune.opts <- lapply(caret.methods, function(mdl) caret.mdls[[mdl]]$parameters)
                   names(tune.opts) <- caret.methods
                   
                   meth.info <- caret.mdls[[input$ml_method]]
                   params = meth.info$parameters
                   #grid.def <- meth.info$grid(training, trainY, len = 1)
                   
                   tuneGrid = expand.grid(
                     {
                       lst = lapply(1:nrow(params), function(i){
                         info = params[i,]
                         inp.val = input[[paste0("ml_", info$parameter)]]
                         # - - check for ranges - -
                         if(grepl(inp.val, pattern=":")){
                           split = strsplit(inp.val,split = ":")[[1]]
                           inp.val <- seq(as.numeric(split[1]),
                                          as.numeric(split[2]),
                                          as.numeric(split[3]))
                         }else if(grepl(inp.val, pattern = ",")){
                           split = strsplit(inp.val,split = ",")[[1]]
                           inp.val <- split
                         }
                         # - - - - - - - - - - - - -
                         switch(as.character(info$class),
                                numeric = as.numeric(inp.val),
                                character = as.character(inp.val))
                       })
                       names(lst) = params$parameter
                       if(any(sapply(lst,function(x)all(is.na(x))))){
                         cat("Missing param, auto-tuning...")
                         lst <- list()
                       }
                       #lst <- lst[sapply(lst,function(x)all(!is.na(x)))]
                       lst
                     })
                   
                   # ============ LOOP HERE ============
                   
                   # get results for the amount of attempts chosen
                   
                   shiny::withProgress(message = "Running...", {
                     repeats <- pbapply::pblapply(1:goes,
                                                  cl = session_cl,
                                                  function(i, 
                                                           train_vec,
                                                           test_vec,
                                                           config,
                                                           configCols,
                                                           ml_method,
                                                           ml_perf_metr,
                                                           ml_folds,
                                                           ml_preproc,
                                                           tuneGrid,
                                                           ml_train_perc,
                                                           sampling){
                                                    MetaboShiny::runML(curr,
                                                                       train_vec = train_vec,
                                                                       test_vec = test_vec,
                                                                       config = config,
                                                                       configCols = configCols,
                                                                       ml_method = ml_method,
                                                                       ml_perf_metr = ml_perf_metr,
                                                                       ml_folds = ml_folds,
                                                                       ml_preproc = ml_preproc,
                                                                       tuneGrid = tuneGrid,
                                                                       ml_train_perc = ml_train_perc,
                                                                       sampling = sampling)
                                                  },
                                                  train_vec = lcl$vectors$ml_train,
                                                  test_vec = lcl$vectors$ml_test,
                                                  config = config,
                                                  configCols = configCols,
                                                  ml_method = input$ml_method,
                                                  ml_perf_metr = input$ml_perf_metr,
                                                  ml_folds = input$ml_folds,
                                                  ml_preproc = input$ml_preproc,
                                                  tuneGrid = tuneGrid,
                                                  ml_train_perc = input$ml_train_perc,
                                                  sampling = if(input$ml_sampling == "none") NULL else input$ml_sampling
                     )
                   })
                   
                   # check if a storage list for machine learning results already exists
                   if(!"ml" %in% names(mSet$analSet)){
                     mSet$analSet$ml <- list() # otherwise make it
                   }
                   
                   mz.imp <- lapply(repeats, function(x) x$importance)
                   # aucs
                   if(length(levels(mSet$dataSet$cls)) > 2){
                     perf <- lapply(1:length(repeats), function(i){
                       x = repeats[[i]]
                       res = MetaboShiny::getMultiMLperformance(x)
                       res$attempt = c(i)
                       res
                     })
                     perf.long <- data.table::rbindlist(perf)
                     mean.auc <- mean(perf.long$AUC_AVG)
                   }else{
                     # save the summary of all repeats (will be used in plots) TOO MEMORY HEAVY
                     pred <- ROCR::prediction(lapply(repeats, function(x) x$prediction), 
                                              lapply(repeats, function(x) x$labels))
                     perf <- ROCR::performance(pred, "tpr", "fpr")
                     perf_auc <- ROCR::performance(pred, "auc")
                     perf.long <- data.table::rbindlist(lapply(1:length(perf@x.values), function(i){
                       xvals <- perf@x.values[[i]]
                       yvals <- perf@y.values[[i]]
                       aucs <- signif(perf_auc@y.values[[i]][[1]], digits = 2)
                       res <- data.table::data.table(attempt = c(i),
                                                     FPR = xvals,
                                                     TPR = yvals,
                                                     AUC_AVG = aucs,
                                                     AUC_PAIR = aucs)
                       res
                     }))
                     perf.long$comparison <- paste0(levels(mSet$dataSet$cls),collapse=" vs. ")
                     mean.auc <- mean(unlist(perf_auc@y.values))
                   }
                   
                   roc_data <- list(m_auc = mean.auc,
                                    perf = perf.long,
                                    imp = mz.imp)
                   
                   bar_data <- data.table::rbindlist(lapply(1:length(repeats), function(i){
                     x = repeats[[i]]
                     tbl = data.table::as.data.table(x$importance, keep.rownames=T)
                     tbl$rep = c(i)
                     colnames(tbl) = c("mz",
                                       "importance",
                                       "rep")
                     # - - - - - - -
                     tbl
                   }))
                   
                   if(input$ml_method == "glmnet"){
                     bar_filt = bar_data[importance > 0]
                     all_mz <- table(bar_filt$mz)
                     tbl = data.table::as.data.table(t(all_mz))[,2:3]
                     colnames(tbl) = c("mz", "importance")
                     tbl$dummy <- c(NA)
                     bar_data <- tbl
                   }
                   # save results to mset
                   if(input$ml_method %not in% mSet$analSet$ml){
                     mSet$analSet$ml[[input$ml_method]] <- list()
                   }
                   
                   mSet$analSet$ml[[input$ml_method]][[input$ml_name]] <- list("roc" = roc_data,
                                                                               "bar" = bar_data)
                   mSet$analSet$ml$last <- list(name = input$ml_name,
                                                method = input$ml_method)
                   
                   #mSet_ml <<- mSet
                   })
               },
               heatmap = {
                 # reset
                 mSet <- MetaboShiny::calcHeatMap(mSet, 
                                                  signif.only = input$heatsign,
                                                  source.anal = input$heattable,
                                                  top.hits = input$heatmap_topn,
                                                  cols = lcl$aes$mycols)
                 output$heatmap_now <- shiny::renderText(input$heattable)
               },
               tt = {
                 withProgress({
                   mSet <- MetaboAnalystR::Ttests.Anal(mSet,
                                                       nonpar = input$tt_nonpar,
                                                       threshp = 0.1, # TODO: make the threshold user defined...
                                                       paired = mSet$dataSet$paired,
                                                       equal.var = input$tt_eqvar
                   )
                 })
               },
               fc = {
                 withProgress({
                   if(mSet$dataSet$paired){
                     mSet <- MetaboAnalystR::FC.Anal.paired(mSet,
                                                             1.5, # TODO: make this threshold user defined
                                                             1)  
                   }else{
                     mSet <- MetaboAnalystR::FC.Anal.unpaired(mSet,
                                                               1.5, # TODO: make this threshold user defined
                                                               1) 
                   }
                   if(!is.null(mSet$analSet$fc$sig.mat)){
                     rownames(mSet$analSet$fc$sig.mat) <- gsub(rownames(mSet$analSet$fc$sig.mat), 
                                                               pattern = "^X",
                                                               replacement = "")
                   }else{
                     mSet$analSet$fc$sig.mat <- data.frame()
                   }
                 })
               },
               aov = {
                 aovtype = if(mSet$dataSet$exp.type %in% c("t", "2f", "t1f")) "aov2" else "aov"
                 redo = aovtype %not in% names(mSet$analSet)
                 if(redo){ # if done, don't redo
                   shiny::withProgress({
                     mSet <- switch(mSet$dataSet$exp.type,
                                    "1fm"=MetaboAnalystR::ANOVA.Anal(mSet, thresh=0.1,post.hoc = "fdr",nonpar = F),
                                    "2f"=MetaboAnalystR::ANOVA2.Anal(mSet, 0.1, "fdr", "", 1, 1),
                                    "t"=MetaboAnalystR::ANOVA2.Anal(mSet, 0.1, "fdr", "time0", 1, 1),
                                    "t1f"=MetaboAnalystR::ANOVA2.Anal(mSet, 0.1, "fdr", "time", 1, 1))
                   })
                 }
               },
               volc = {
                 shiny::withProgress({
                   mSet <-  MetaboAnalystR::Volcano.Anal(mSet,
                                                         paired = mSet$dataSet$paired, 
                                                         1.5, 0,
                                                         0.75, F, 0.1,
                                                         TRUE, "raw") # TODO: make thresholds user-defined
                 })
               },
               tsne = {
                 shiny::withProgress({
                   coords = tsne::tsne(mSet$dataSet$norm, k = 3,
                                       initial_dims = input$tsne_dims,
                                       perplexity = input$tsne_perplex,
                                       max_iter = input$tsne_maxiter)
                   colnames(coords) <- paste("t-SNE dimension", 1:3)
                   mSet$analSet$tsne <- list(x = coords)
                 })
               },
               plsda = {
                 library(e1071)
                 library(pls)
                 # depending on type, do something else
                 # TODO: enable sparse and orthogonal PLS-DA
                 switch(input$plsda_type,
                        normal={
                          require(caret)
                          withProgress({
                            mSet <- MetaboAnalystR::PLSR.Anal(mSet) # perform pls regression
                            setProgress(0.3)
                            mSet <- MetaboAnalystR::PLSDA.CV(mSet, methodName=if(nrow(mSet$dataSet$norm) < 50) "L" else "T",compNum = 3) # cross validate
                            setProgress(0.6)
                            mSet <- MetaboAnalystR::PLSDA.Permut(mSet,num = 300, type = "accu") # permute
                          })
                        },
                        sparse ={
                          mSet <- MetaboAnalystR::SPLSR.Anal(mSet, comp.num = 3)
                        })
               },
               power = {
                 shiny::withProgress({
                   pwr.analyses = lapply(input$power_comps, function(combi){
                     mSet.temp <- MetaboAnalystR::InitPowerAnal(mSet, combi)
                     mSet.temp <- MetaboAnalystR::PerformPowerProfiling(mSet.temp, 
                                                                        fdr.lvl = input$power_fdr, 
                                                                        smplSize = input$power_nsamp)  
                     shiny::incProgress(amount = 1 / length(input$power_comps))
                     mSet.temp$analSet$power
                   })   
                 }, max = length(input$power_comps))
                 names(pwr.analyses) <- input$power_comps
                 mSet$analSet$power <- pwr.analyses 
               },
               wordcloud = {
                 if(input$wordcloud_manual){
                   shiny::withProgress({ # progress bar
                     abstracts = MetaboShiny::getAbstracts(searchTerms = input$wordcloud_searchTerm,
                                                           mindate = input$wordcloud_dateRange[1],
                                                           maxdate = input$wordcloud_dateRange[2],
                                                           retmax = input$wordcloud_absFreq)
                     lcl$tables$abstracts <- abstracts
                     shiny::setProgress(0.5)
                     lcl$tables$wordcloud_orig <- MetaboShiny::getWordFrequency(abstracts$abstract)
                     lcl$tables$wordcloud_filt <- lcl$tables$wordcloud_orig
                   }, message = "Searching...", max = 1)
                 }else{
                   if(nrow(shown_matches$forward_full) > 0){
                     lcl$tables$wordcloud_orig <- MetaboShiny::getWordFrequency(shown_matches$forward_full$description)
                     lcl$tables$wordcloud_filt <- lcl$tables$wordcloud_orig
                   }
                 }
               }) 
        if(typeof(mSet) != "double"){
          success <- T
        }
      })
      
      if(success){
        mSet <<- mSet
        lcl$hasChanged <<- TRUE
        shinyjs::show(selector = paste0("div.panel[value=collapse_", statsmanager$calculate, "_plots]"))
        shinyjs::show(selector = paste0("div.panel[value=collapse_", statsmanager$calculate, "_tables]"))
        shinyBS::updateCollapse(session, paste0("collapse_",input$statistics),open = paste0("collapse_", 
                                                                                            statsmanager$calculate, 
                                                                                            c("_tables","_plots")))
      }else{
        MetaboShiny::metshiAlert("Analysis failed!")
        mSet <<- mSet.old
      }
      lcl <<- lcl
    }
    # - - - -
    statsmanager$calculate <- NULL # set reloading to 'off'
  }
})
