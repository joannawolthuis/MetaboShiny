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
                 mSet <- MetaboShiny::store.mSet(mSet, name = mSet$settings$cls.name, proj.folder = lcl$paths$work_dir)
               },
               corr = {
                 # pearson kendall spearman
                 lvls = levels(mSet$dataSet$cls)
                 pat = input$corr_seq_order
                 pat_order = match(lvls,pat)
                 pattern = paste0(pat_order-1, collapse="-")
                 mSet <- MetaboAnalystR::Match.Pattern(mSet, input$corr_corr, pattern)
                 # == filter ===
                 dt = mSet$analSet$corr$cor.mat
                 pthresh = input$corr_p_thresh #0.1
                 corrthresh = input$corr_r_thresh #0.1
                 for(col in colnames(dt)){
                   mSet$analSet$corr[[col]] <- dt[,col]
                 }
                 keepMe = abs(mSet$analSet$corr$cor.mat[,'correlation']) >= corrthresh & 
                   mSet$analSet$corr$cor.mat[,'p-value'] <= pthresh
                 mSet$analSet$corr$cor.mat <- mSet$analSet$corr$cor.mat[keepMe,]
               },
               diffcorr = {
                 library(DGCA, quietly = TRUE)
                 #data(darmanis)
                 data(design_mat)
                 #ddcor_res = ddcorAll(inputMat = darmanis, design = design_mat,
                 #                     compare = c("oligodendrocyte", "neuron"),
                 #                     adjust = "none", nPerm = 0, nPairs = 100)
                 inMat = t(mSet$dataSet$norm)
                 vars = unique(unlist(mSet$dataSet$covars[,mSet$settings$exp.var,with=F]))
                 dMat = matrix(ncol = length(vars), 
                               nrow=nrow(mSet$dataSet$covars), 
                               data = rep(0, length(vars)),
                               dimnames = list(mSet$dataSet$covars$sample, vars))
                 
                 for(var in vars){
                   dMat[,var] <- as.numeric(unlist(mSet$dataSet$covars[,mSet$settings$exp.var,with=F] == var))
                 }
                 
                 mSet$analSet$diffcorr = DGCA::ddcorAll(inputMat = inMat, design = dMat,
                                                        compare = colnames(dMat),
                                                        adjust = "fdr", nPerm = 0)#10)# nPairs = 100)
               },
               pca = {
                 shiny::withProgress({
                   if(input$pca_source != "normalized"){
                     mSet_orig = mSet
                     mSet$dataSet$norm <- switch(input$pca_source,
                                                 "pre-batch correction" = mSet$dataSet$prebatch,
                                                 original = mSet$dataSet$proc)
                     pcaRes <- MetaboAnalystR::PCA.Anal(mSet)$analSet$pca # perform PCA analysis
                     mSet = mSet_orig
                     mSet$analSet$pca <- pcaRes  
                   }else{
                     mSet <- MetaboAnalystR::PCA.Anal(mSet) # perform PCA analysis
                   }
                 })
                 success = T
               },
               ica = {
                 shiny::withProgress({
                   # ica package
                   inTbl = switch(input$ica_source, 
                                  original = mSet$dataSet$proc,
                                  "pre-batch correction" = mSet$dataSet$prebatch,
                                  normalized = mSet$dataSet$norm)
                   nc = input$ica_ncomp
                   icaRes = switch(input$ica_method,
                                   fast = ica::icafast(inTbl, nc = input$ica_ncomp, center = F,maxit = input$ica_maxiter),
                                   imax = ica::icaimax(inTbl, nc = input$ica_ncomp, center = F,maxit = input$ica_maxiter),
                                   jade = ica::icajade(inTbl, nc = input$ica_ncomp, center = F,maxit = input$ica_maxiter))
                   mSet$analSet$ica <- icaRes
                 })
                 success = T
               },
               umap = {
                 shiny::withProgress({
                   # umap package
                   inTbl = switch(input$umap_source, 
                                  original = mSet$dataSet$proc,
                                  "pre-batch correction" = mSet$dataSet$prebatch,
                                  normalized = mSet$dataSet$norm)
                   umapRes = umap::umap(d = inTbl,
                                        method = "naive",
                                        n_components = input$umap_ncomp,
                                        n_neighbors = input$umap_neighbors)
                   mSet$analSet$umap <- umapRes
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
               network = {
                 
                 if(input$network_sel){
                   useHits = colnames(mSet$dataSet$norm)
                 }else{
                   # aov
                   flattened <- getTopHits(mSet, 
                                           input$network_table,
                                           input$network_topn)
                   
                   useHits = flattened[[1]]  
                 }
                 
                 
                 # ---
                 #TODO: gaussian graphical model
                 # ---
                 rcorrMat = Hmisc::rcorr(x = as.matrix(mSet$dataSet$norm[, useHits]),
                                         #y = as.matrix(mSet$dataSet$norm[, useHits]),
                                         type = input$network_corr)
                 mSet$analSet$network <- list(rcorr = rcorrMat,
                                              order = useHits)
                 output$network_now <- shiny::renderText(input$network_table)
               },
               enrich = {
                 shiny::withProgress({
                   
                   #similarly to venn diagram
                   flattened <- getTopHits(mSet, 
                                           input$mummi_anal,
                                           input$mummi_topn)
                   
                   hasP = grepl("tt|aov|asca",input$mummi_anal)
                   
                   setProgress(0.1)
                   
                   myFile <- tempfile(fileext = ".csv")
                   tbl = data.table::data.table("m.z" = as.numeric(gsub(flattened[[1]], pattern="(\\+|\\-|RT).*$", replacement="")),
                                                mode = sapply(flattened[[1]], function(mz){
                                                  if(grepl(pattern="-",x=mz)) "negative" else "positive"
                                                }))
                   tbl <- tbl[complete.cases(tbl)]
                   
                   hasT = grepl("tt", input$mummi_anal)
                   
                   anal = gsub(" \\(.*$|", "", input$mummi_anal)
                   subset = gsub("\\(|\\)|.*\\(", "", input$mummi_anal)
                   
                   tbl[, "p.value"] = if(hasP) mSet$storage[[subset]]$analysis[[anal]]$sig.mat[match(flattened[[1]], 
                                                                                                     rownames(mSet$storage[[subset]]$analysis[[anal]]$sig.mat)),
                                                                                               if(anal == "aov2") "Interaction(adj.p)" else "p.value"] else c(0)
                   tbl[, "t.score"] = if(hasT) mSet$storage[[subset]]$analysis[[anal]]$sig.mat[match(flattened[[1]], 
                                                                                                     rownames(mSet$storage[[subset]]$analysis[[anal]]$sig.mat)),
                                                                                               "t.stat"] else c(0)
                   
                   if(hasP) if(all(is.na(tbl$p.value))) tbl$p.value <- c(0)
                   if(hasT) if(all(is.na(tbl$t.score))) tbl$t.score <- c(0)
                   
                   tmpfile <- tempfile()
                   
                   fwrite(tbl, file=tmpfile)
                   
                   enr_mSet <- MetaboAnalystR::InitDataObjects("mass_all",
                                                               "mummichog",
                                                               FALSE)
                   MetaboAnalystR::SetPeakFormat("mpt")
                   enr_mSet <- MetaboAnalystR::UpdateInstrumentParameters(enr_mSet,
                                                                          mSet$ppm,
                                                                          "mixed",
                                                                          "yes",
                                                                          0.02);
                   
                   enr_mSet <- MetaboAnalystR::Read.PeakListData(enr_mSet, tmpfile);
                   
                   shiny::setProgress(0.2)
                   
                   enr_mSet <- MetaboAnalystR::SanityCheckMummichogData(enr_mSet)
                   enr_mSet <- MetaboAnalystR::Setup.AdductData(enr_mSet, input$mummi_adducts);
                   
                   shiny::setProgress(0.3)
                   
                   #===
                   
                   enr_mSet<-MetaboAnalystR::SetPeakEnrichMethod(enr_mSet, if(input$mummi_enr_method | !hasT) "mum" else "gsea", "v2")
                   enr_mSet<-MetaboAnalystR::SetMummichogPval(enr_mSet, if(hasP) as.numeric(gsub(",",".",input$mummi_pval)) else 1)
                   
                   #===
                   
                   elecMass = 0.000548579909
                   mummi_adducts <- adducts[Name %in% input$mummi_adducts]
                   add_db_custom_rows <- lapply(1:nrow(mummi_adducts), function(i){
                     row = mummi_adducts[i,]
                     if(is.na(row$AddAt)) row$AddAt <- ""
                     if(is.na(row$AddEx)) row$AddEx <- ""
                     if(is.na(row$RemAt)) row$RemAt <- ""
                     if(is.na(row$RemEx)) row$RemEx <- ""
                     addForm <- enviPat::mergeform(row$AddAt, row$AddEx)
                     remForm <- enviPat::mergeform(row$RemAt, row$RemEx)
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
                   
                   shiny::setProgress(0.6)
                   
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
                     if(!is.null(mw_modified.neg)) colnames(mw_modified.neg) <- ion.name.neg
                     mass.list.pos <- as.list(ion.mass.pos)
                     mass.user.pos <- lapply(mass.list.pos, function(x) eval(parse(text = paste(gsub("PROTON", 
                                                                                                     1.007825, x)))))
                     mw_modified.pos <- do.call(cbind, mass.user.pos)
                     if(!is.null(mw_modified.pos)) colnames(mw_modified.pos) <- ion.name.pos
                     mw_modified <- list(mw_modified.neg, mw_modified.pos)
                     
                     if(use.rules){
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
                   
                   assignInNamespace("new_adduct_mzlist", new_adduct_mzlist, ns="MetaboAnalystR", 
                                     envir=as.environment("package:MetaboAnalystR"))
                   
                   shiny::setProgress(0.9)
                   enr_mSet$dataSet$N <- 20
                   
                   enr_mSet <- MetaboAnalystR::PerformPSEA(mSetObj = enr_mSet, 
                                                           lib = input$mummi_org,
                                                           libVersion = "current",
                                                           permNum = 100) 
                   
                   filenm <- if(input$mummi_enr_method | !hasT) "mummichog_matched_compound_all.csv" else "mummichog_fgsea_pathway_enrichment.csv"
                   enr_mSet$dataSet$mumResTable <- data.table::fread(filenm)
                   
                   if(!input$mummi_enr_method){
                     tbl.rows <- lapply(1:length(enr_mSet$path.hits), function(i){
                       l = enr_mSet$path.hits[[i]]
                       row = enr_mSet$dataSet$mumResTable[i,]
                       row$Cpd.Hits <- paste0(l, collapse=";")
                       row
                     })
                     tbl <- data.table::rbindlist(tbl.rows)
                     tbl <- tidyr::separate_rows(tbl,
                                                 "Cpd.Hits",
                                                 sep = ";")
                     tbl <- as.data.frame(tbl)
                     tbl$Matched.Compound <- tbl$Cpd.Hits
                     
                     add.tbl = data.frame(Matched.Compound = names(unlist(enr_mSet$cpd_form_dict)),
                                          adduct = unlist(enr_mSet$cpd_form_dict))
                     mz.tbl = data.frame(Query.Mass = names(unlist(enr_mSet$mz2cpd_dict)),
                                         Matched.Compound = unlist(enr_mSet$mz2cpd_dict))
                     mzorig.tbl = data.frame(Matched.Compound = names(unlist(enr_mSet$cpd2mz_dict)),
                                             Orig.Mass = unlist(enr_mSet$cpd2mz_dict))
                     mergy = merge(tbl, add.tbl)
                     mergy = merge(mergy, mz.tbl)
                     mergy = merge(mergy, mzorig.tbl)
                     mergy$Mass.Diff = abs(as.numeric(mergy$Query.Mass) - as.numeric(mergy$Orig.Mass))
                     mergy = unique(mergy[,c("Query.Mass", "Matched.Compound","adduct","Mass.Diff")])
                     enr_mSet$dataSet$mumResTable <- mergy
                   }
                   
                   mSet$analSet$enrich <- list(mummi.resmat = enr_mSet$mummi.resmat,
                                               mummi.gsea.resmat = enr_mSet$mummi.gsea.resmat,
                                               mumResTable = enr_mSet$dataSet$mumResTable,
                                               path.nms = enr_mSet$path.nms,
                                               path.hits = enr_mSet$path.hits)
                   enr_mSet <- NULL
                 })
               },
               ml = {
                 try({
                   
                   {
                     assign("ml_queue", shiny::isolate(shiny::reactiveValuesToList(ml_queue)), envir = .GlobalEnv)
                     assign("input", shiny::isolate(shiny::reactiveValuesToList(input)), envir = .GlobalEnv)
                    
                     if(session_cl != 0){
                       parallel::clusterExport(session_cl, c("input", "ml_queue", "gbl", "lcl"))
                     }
                     
                     # make subsetted mset for ML so it's not as huge in memory
                     small_mSet = mSet
                     small_mSet$dataSet = small_mSet$dataSet[c("cls", "orig.cls", "orig", "norm", "covars")]
 
                     ml_run <- function(settings, mSet, input, cl){
                       res = list()
                       #({
                       {
                         tmpdir = file.path(tempdir(), settings$ml_name) # needed for 
                         if(!dir.exists(tmpdir)) dir.create(tmpdir)
                         setwd(tmpdir)
                         
                         # PIPELINE
                         # pick source table
                         pickedTbl <- settings$ml_used_table
                         
                         if(pickedTbl == "pca" & !("pca" %in% names(mSet$analSet))){
                           stop("Please run PCA first!")
                         }
                         
                         # covars needed
                         keep.config = setdiff(c(settings$ml_include_covars, settings$ml_batch_covars,
                                                 settings$ml_train_subset[[1]], settings$ml_test_subset[[1]]),
                                               "label")
                         config = mSet$dataSet$covars[,..keep.config]
                         config$label = mSet$dataSet$cls
                         
                         # train/test split
                         ## get indices
                         if(!is.null(settings$ml_train_subset) | !is.null(settings$ml_test_subset)){
                           # add clause for same train_test
                           test_idx = NULL
                           train_idx = NULL
                           if(!is.null(settings$ml_test_subset)){
                             test_idx = which(config[[settings$ml_test_subset[[1]]]] %in% settings$ml_test_subset[[2]])
                           }
                           if(!is.null(settings$ml_train_subset)){
                             train_idx = which(config[[settings$ml_train_subset[[1]]]] %in% settings$ml_train_subset[[2]])
                           }
                           if(is.null(train_idx)){
                             train_idx = setdiff(1:nrow(config), test_idx)  
                           }else if(is.null(test_idx)){
                             test_idx = setdiff(1:nrow(config), train_idx)
                           }
                         }else{
                           # make joined label of label+batch and split that train/test
                           split_label = if(length(settings$ml_batch_covars) > 0){
                             #print("splitting tr/te % per batch")
                             covars = c("label", settings$ml_batch_covars)
                             apply(config[, ..covars],
                                   MARGIN = 1, 
                                   function(x) paste0(x, collapse="_"))
                           }else{
                             #print("unbiased split over pool")
                             config$label
                           }
                           train_idx = caret::createDataPartition(y = split_label, 
                                                                  p = settings$ml_train_perc/100,
                                                                  list = FALSE)[,1] # partition data in a balanced way (uses labels)
                           
                           test_idx = setdiff(1:nrow(config), train_idx)
                           #table(split_label[test_idx])
                         }
                         
                         # subset to specific m/z values used
                         mzs = c()
                         if(pickedTbl != "pca"){
                           if(settings$ml_specific_mzs != "no"){
                             msg = "Using user-specified m/z set."
                             if(!is.null(settings$ml_mzs)){
                               curr <- curr[,settings$ml_mzs, with=F]
                             }else{
                               mzs = getTopHits(mSet = mSet,
                                                expnames = settings$ml_specific_mzs, 
                                                top = settings$ml_mzs_topn, 
                                                filter_mode = "top")[[1]]
                               mzs = gsub("^X", "", mzs)
                               mzs = gsub("\\.$", "-", mzs)
                             }
                           }   
                         }
                         
                         if(mSet$metshiParams$renorm & pickedTbl != "pca"){
                           samps_train = rownames(mSet$dataSet$norm)[train_idx]
                           samps_test = rownames(mSet$dataSet$norm)[test_idx]
                           mSet.settings = mSet$settings
                           reset_mSet <- reset.mSet(mSet,
                                                    fn = file.path(lcl$paths$proj_dir, 
                                                                   paste0(lcl$proj_name,
                                                                          "_ORIG.metshi")))
                           # --- GET TRAIN ---
                           mSet_train = subset_mSet(reset_mSet, "sample", samps_train)
                           mSet_train = change.mSet(mSet_train, 
                                                    stats_var = mSet.settings$exp.var, 
                                                    time_var =  mSet.settings$time.var,
                                                    stats_type = mSet.settings$exp.type)
                           if(length(mzs) > 0){
                             mSet_train = subset_mSet_mz(mSet_train, mzs)
                           }
                           mSet_train$dataSet$orig <- mSet_train$dataSet$start
                           mSet_train$dataSet$start <- mSet_train$dataSet$preproc <- mSet_train$dataSet$proc <- mSet_train$dataSet$prenorm <- NULL
                           mSet_train = metshiProcess(mSet_train, init = F, cl = 0)
                           # --- GET TEST ---
                           mSet_test = subset_mSet(reset_mSet, "sample", samps_test)
                           if(length(mzs) > 0){
                             mSet_test = subset_mSet_mz(mSet_test, mzs)
                           }
                           mSet_test = change.mSet(mSet_test, 
                                                   stats_var = mSet.settings$exp.var, 
                                                   time_var =  mSet.settings$time.var,
                                                   stats_type = mSet.settings$exp.type)
                           mSet_test$dataSet$orig <- mSet_test$dataSet$start
                           mSet_test$dataSet$start <- mSet_test$dataSet$preproc <- mSet_test$dataSet$proc <- mSet_test$dataSet$prenorm <- NULL
                           mSet_test = metshiProcess(mSet_test, init = F, cl = 0)
                           # ------- rejoin and create curr -------
                           config_train = mSet_train$dataSet$covars[,..keep.config]
                           config_train$label = mSet_train$dataSet$cls
                           config_test = mSet_test$dataSet$covars[,..keep.config]
                           config_test$label = mSet_test$dataSet$cls
                           
                           mz.in.both = intersect(colnames(mSet_train$dataSet$norm),
                                                  colnames(mSet_test$dataSet$norm))
                           
                           curr = rbind(mSet_train$dataSet$norm[,mz.in.both],
                                        mSet_test$dataSet$norm[,mz.in.both])
                           config = rbind(config_train, 
                                          config_test)
                           
                           mSet_test <- mSet_train <- config_test <- config_train <- mz.in.both <- mSet$storage <- NULL
                         }else{
                           curr = as.data.frame(switch(pickedTbl, 
                                                       orig = mSet$dataSet$orig,
                                                       norm = mSet$dataSet$norm,
                                                       pca = mSet$analSet$pca$x))
                           if(length(mzs) > 0){
                             curr <- curr[, mzs]
                           }
                         }

                         mSet$storage <- NULL
                         mSet = NULL

                         test_sampnames = rownames(curr)[test_idx]
                         
                         # PCA correct
                         if(settings$ml_pca_corr & pickedTbl != 'pca'){
                           #print("Performing PCA and subtracting PCs...")
                           curr <- pcaCorr(curr, 
                                           center = if(pickedTbl == "norm") F else T,
                                           scale = if(pickedTbl == "norm") F else T, 
                                           start_end_pcs = settings$ml_keep_pcs)
                         }
                         
                         
                         #@ split
                         training_data = list(curr = curr[train_idx,],
                                              config = config[train_idx,])
                         
                         testing_data = list(curr = curr[test_idx,],
                                             config = config[test_idx,])
                         
                         curr = NULL

                         if(length(settings$ml_batch_covars) == 0){
                           settings$ml_batch_balance <- settings$ml_batch_sampling <- F
                         }
                         
                         # resampling
                         if(settings$ml_sampling != "none"){
                           
                           # split on factor (either batch or placeholder to create one result)
                           if(settings$ml_batch_balance){
                             split.fac = training_data$config[, settings$ml_batch_covars, with=F][[1]]
                           }else{
                             split.fac = rep(1, nrow(training_data$config))
                           }
                           split.fac = if(settings$ml_batch_balance){
                             training_data$config[, settings$ml_batch_covars, with=F][[1]]
                           }else{rep(1, nrow(training_data$config))
                           }
                           
                           spl.testing.idx = split(1:nrow(training_data$curr), split.fac)
                           balance.overview = sapply(spl.testing.idx, 
                                                     function(idx) table(training_data$config[idx, settings$ml_batch_covars, with=F]))
                           biggest.group.overall = max(balance.overview)
                           smallest.group.overall = min(balance.overview)
                           
                           orig.samp.distr = table(training_data$config$label)
                           
                           training_data$config = training_data$config
                           testing_data$config = testing_data$config
                           
                           size.global = if(settings$ml_sampling != "down") biggest.group.overall else smallest.group.overall
                           
                           # resample
                           resampled.data.list = lapply(spl.testing.idx, function(idx){
                             
                             curr.subset = training_data$curr[idx,]
                             config.subset = training_data$config[idx,]
                             config.top.row = config.subset[1,]
                             
                             # total group size within this loop?
                             size.local = if(settings$ml_sampling != "down") max(table(config.subset$label)) else min(table(config.subset$label))
                             
                             # all upsampling except "down"
                             group.size = 
                               if(settings$ml_batch_balance){
                                 if(settings$ml_batch_size_sampling) size.global else size.local
                               }else size.local
                             
                             # resample
                             ## K for the k-fold methods
                             K = min(min(table(config.subset$label))-1, 10)
                             mz.names = colnames(curr.subset)
                             colnames(curr.subset) <- paste0("mz",1:ncol(curr.subset))
                             
                             switch(settings$ml_sampling,
                                    up = {
                                      nconfig = ncol(config.subset)
                                      new_data = upsample.adj(cbind(config.subset, curr.subset), 
                                                              as.factor(config.subset$label), 
                                                              maxClass = group.size) 
                                      config.subset = new_data[,1:nconfig]
                                      curr.subset = new_data[,!(colnames(new_data) %in% c("Class", colnames(config.subset)))]
                                    },
                                    adasyn = {
                                      resampled = smotefamily::ADAS(X = curr.subset, 
                                                                    target = config.subset$label,
                                                                    K = K)
                                      new_data = resampled$data
                                      curr.subset = new_data[, colnames(new_data) != "class"]
                                      config.subset = data.table(label = new_data$class)
                                    },
                                    smote = {
                                      resampled = smotefamily::SMOTE(X = curr.subset, 
                                                                     target = config.subset$label,
                                                                     K = K)
                                      new_data = resampled$data
                                      curr.subset = new_data[, colnames(new_data) != "class"]
                                      config.subset = data.table(label = new_data$class)
                                    },
                                    rose = {
                                      resampled = ROSE::ROSE(label ~ ., 
                                                             data = cbind(label = config.subset$label, curr.subset),
                                                             N = group.size)
                                      new_data = resampled$data
                                      curr.subset = new_data[,2:(ncol(new_data))]
                                      config.subset = data.table(label = new_data[[1]])
                                    },
                                    down = {
                                      nconfig = ncol(config.subset)
                                      new_data = downsample.adj(cbind(config.subset, curr.subset), 
                                                                as.factor(config.subset$label), 
                                                                minClass = group.size) 
                                      config.subset = new_data[,1:nconfig]
                                      curr.subset = new_data[,!(colnames(new_data) %in% c("Class", colnames(config.subset)))]
                                      
                                    }
                             )
                             colnames(curr.subset) <- mz.names
                             
                             if(settings$ml_batch_balance & settings$ml_sampling %in% c("rose", "smote", "adasyn")){
                               config.subset[[settings$ml_batch_covars]] <- c(config.top.row[[settings$ml_batch_covars]])
                             }
                             
                             list(curr = curr.subset,
                                  config = config.subset)
                           })
                           training_data = list(
                             curr = data.table::rbindlist(lapply(resampled.data.list, function(x) x$curr),use.names = T),
                             config = data.table::rbindlist(lapply(resampled.data.list, function(x) x$config), use.names = T)
                           )
                           resampled.data.list <- NULL
                         }
                         
                         # divvy folds for cross validation
                         folds <- if(length(settings$ml_batch_covars) > 0 &
                                     settings$ml_folds != "LOOCV" &
                                     (settings$ml_sampling %in% c("up","down","none") | settings$ml_batch_balance)){
                           #print("Creating CV folds based on batch factor(s) and label.")
                           covars = c("label", settings$ml_batch_covars)
                           split_label = apply(training_data$config[, ..covars],
                                               MARGIN = 1, 
                                               function(x) paste0(x, collapse="_"))
                           caret::groupKFold(split_label, k = min(as.numeric(settings$ml_folds), 
                                                                  length(unique(split_label))))
                         }else NULL
                         
                         # for(fold in folds){
                         #   print("---")
                         #   print(table(training_data$config$country[fold]))
                         #   print(table(training_data$config$label[fold]))
                         # }
                         
                         # replace training data with the new stuff
                         training_data$config$split <- "train"
                         testing_data$config$split <- "test"
                         
                         # remove columns that should not be in prediction
                         keep.config = unique(c("label", "split", settings$ml_include_covars))
                         testing_data$config = testing_data$config[, ..keep.config]
                         training_data$config = training_data$config[, ..keep.config]
                         
                         # merge back into one
                         training = cbind(training_data$config,
                                          training_data$curr)
                         testing = cbind(testing_data$config[,colnames(training_data$config),with=F],
                                         testing_data$curr)
                         
                         training_data <- testing_data <- NULL
                         
                         # CARET SETTINGS
                         caret.mdls <- caret::getModelInfo()
                         caret.methods <- names(caret.mdls)
                         tune.opts <- lapply(caret.methods, function(mdl) caret.mdls[[mdl]]$parameters)
                         names(tune.opts) <- caret.methods
                         
                         meth.info <- caret.mdls[[settings$ml_method]]
                         params = meth.info$parameters
                         
                         tuneGrid = expand.grid(
                           {
                             lst = lapply(1:nrow(params), function(i){
                               info = params[i,]
                               inp.val = settings[[paste0("ml_", info$parameter)]]
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
                               #cat("Missing param, auto-tuning...")
                               lst <- list()
                             }
                             #lst <- lst[sapply(lst,function(x)all(!is.na(x)))]
                             lst
                           })
                         
                         # make sure levels of predicted class aren't numeric
                         levels(training$label) <- paste0("class",  Hmisc::capitalize(levels(training$label)))
                         levels(training$label) <- ordered(levels(training$label))
                         levels(testing$label) <- paste0("class",  Hmisc::capitalize(levels(testing$label)))
                         levels(testing$label) <- ordered(levels(testing$label))
                         
                         # correct mzs in case model cannot handle numeric column names
                         colnames(training) <- make.names(colnames(training))
                         colnames(testing) <- make.names(colnames(testing))
                         
                         # run ML
                         ml_res = runML(training = training,
                                        testing = testing,
                                        ml_method = settings$ml_method,
                                        ml_perf_metr = settings$ml_perf_metr,
                                        ml_folds = settings$ml_folds,
                                        ml_preproc = settings$ml_preproc,
                                        tuneGrid = tuneGrid,
                                        folds = folds,
                                        maximize = T,
                                        shuffle = settings$ml_label_shuffle,
                                        n_permute = settings$ml_n_shufflings,
                                        shuffle_mode = if(settings$ml_shuffle_mode) "train" else "test",
                                        cl = cl)
                         
                         res = list(res = ml_res, params = settings)
                       }
                       #})
                       res
                     }
                    
                     # # set static train/test
                     # joint_lbl = paste0(small_mSet$dataSet$covars$country,"_",small_mSet$dataSet$covars$new_group)
                     # train = caret::createDataPartition(joint_lbl, p = 0.8)[[1]]
                     # train_samps = small_mSet$dataSet$covars$sample[train]
                     # basejob = ml_queue$jobs[[1]]
                     # jobs = list()
                     # for( i in 1:26){
                     #   for(j in 1:30){
                     #     job = basejob
                     #     job$ml_mzs_topn = i
                     #     job$ml_name = gsub("1$", paste(i, paste0("#", j)), job$ml_name)
                     #     jobs[[job$ml_name]] = job
                     #   }
                     # }
                     # ml_queue$jobs = jobs
                     
                     basejob = ml_queue$jobs[[1]]
                     jobs = list()
                     for(i in 1:110){
                       for(j in 1:10){
                         job = basejob
                         job$ml_mzs_topn = i
                         job$ml_name = gsub("1$", paste(i, paste0("#", j)), job$ml_name)
                         jobs[[job$ml_name]] = job
                       }
                     }
                     
                     ml_queue$jobs = jobs
                     
                     parallel::clusterExport(session_cl, c("ml_run", "small_mSet", "gbl"))
                     
                     ml_queue_res <- pbapply::pblapply(ml_queue$jobs, cl=session_cl, function(settings, ml_cl){
                       res = list()
                       try({
                         res = ml_run(settings = settings, 
                                      mSet = small_mSet,
                                      input = input,
                                      cl = ml_cl)  
                       })
                       res
                     }, ml_cl = session_cl)
                   }
                   
                   print("Done!")
                   
                   closeAllConnections()
                   
                   shiny::showNotification("Gathering results...")
                   
                   ml_queue_res = ml_queue_res[unlist(sapply(ml_queue_res, function(l) length(l) > 0))]
                   
                   for(res in ml_queue_res){
                     mSet$analSet$ml[[res$params$ml_method]][[res$params$ml_name]] <- res
                     settings <- res$params
                   }
                   
                   ########
                   mSet$analSet$ml$last <- list(name = settings$ml_name,
                                                method = settings$ml_method)
                   
                   shiny::showNotification("Done! Restoring parallel setup from before ML...")
                   logfile = file.path(lcl$paths$work_dir, "metshiLog.txt")
                   if(file.exists(logfile)) file.remove(logfile)
                   session_cl <<- parallel::makeCluster(input$ncores,outfile=logfile)#,setup_strategy = "sequential") # leave 1 core for general use and 1 core for shiny session
                   # send specific functions/packages to other threads
                   parallel::clusterEvalQ(session_cl, {
                     library(data.table)
                     library(iterators)
                     library(MetaboShiny)
                     library(MetaDBparse)
                   })
                   lcl$vectors$ml_train <- lcl$vectors$ml_train <<- NULL
                   print("ML done and saved.")
                 })
               },
               heatmap = {
                 # reset
                 mSet <- calcHeatMap(mSet, 
                                     signif.only = input$heatsign,
                                     source.anal = input$heattable,
                                     top.hits = input$heatmap_topn,
                                     cols = lcl$aes$mycols,
                                     which.data = input$heatmap_source)
                 output$heatmap_now <- shiny::renderText(input$heattable)
               },
               tt = {
                 withProgress({
                   mSet <- MetaboAnalystR::Ttests.Anal(mSet,
                                                       nonpar = input$tt_nonpar,
                                                       threshp = input$tt_p_thresh,
                                                       paired = mSet$dataSet$ispaired,
                                                       equal.var = input$tt_eqvar
                   )
                 })
               },
               fc = {
                 withProgress({
                   mSet$dataSet$combined.method = T
                   if(mSet$dataSet$ispaired){
                     mSet <- MetaboAnalystR::FC.Anal.paired(mSet,
                                                            as.numeric(input$fc_thresh),
                                                            1)  
                   }else{
                     mSet <- MetaboAnalystR::FC.Anal.unpaired(mSet,
                                                              as.numeric(input$fc_thresh), # TODO: make this threshold user defined
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
                 aovtype = if(mSet$settings$exp.type %in% c("t", "2f", "t1f")) "aov2" else "aov"
                 redo = aovtype %not in% names(mSet$analSet)
                 if(redo){ # if done, don't redo
                   shiny::withProgress({
                     mSet <- switch(mSet$settings$exp.type,
                                    "1fm"=MetaboAnalystR::ANOVA.Anal(mSet, thresh=0.1,post.hoc = "fdr",nonpar = F),
                                    "2f"=MetaboAnalystR::ANOVA2.Anal(mSet, 0.1, "fdr", "", 1, 1),
                                    "t"=MetaboAnalystR::ANOVA2.Anal(mSet, 0.1, "fdr", "time0", 1, 1),
                                    "t1f"=MetaboAnalystR::ANOVA2.Anal(mSet, 0.1, "fdr", "time", 1, 1))
                   })
                 }
               },
               combi = {
                 
                 anal1 = input$combi_anal1 #"corr"
                 anal2 = input$combi_anal2 #"aov"
                 anal1_col = input$combi_anal1_var #"correlation"
                 anal2_col = input$combi_anal2_var #"-log10(p)"
                 anal1_res = mSet$analSet[[anal1]]
                 anal2_res = mSet$analSet[[anal2]]
                 anal1_res_table = data.table::as.data.table(anal1_res[grepl("\\.mat", names(anal1_res))][[1]],keep.rownames=T)
                 anal2_res_table = data.table::as.data.table(anal2_res[grepl("\\.mat", names(anal2_res))][[1]],keep.rownames=T)
                 
                 translator=list(
                   "t.stat" = "t.score",
                   "p.value" = "p.value",
                   "-log10(p)" = "p.log",
                   "Fold Change" = "fc.all",
                   "log2(FC)" = "fc.log",
                   "correlation" = "correlation"
                 )
                 
                 if(anal1_col %in% names(translator)){
                   anal1_all_res = anal1_res[[translator[[anal1_col]]]]
                 }else{
                   anal1_all_res = anal1_res_table[[anal1_col]]
                   names(anal1_all_res) <- anal1_res_table$rn
                 }
                 
                 if(anal2_col %in% names(translator)){
                   anal2_all_res = anal2_res[[translator[[anal2_col]]]]
                 }else{
                   anal2_all_res = anal2_res_table[[anal2_col]]
                   names(anal2_all_res) <- anal2_res_table$rn
                 }
                 
                 anal1_trans = "none"
                 anal2_trans = "none"
                 
                 # if(input$combi_anal1_trans != "none"){
                 #   try({
                 #     transFun= switch(input$combi_anal1_trans,
                 #                      "log10"=log10,
                 #                      "-log10"=function(x) -log10(x),
                 #                      "abs"=abs) 
                 #     anal1_res_table[[anal1_col]] <- transFun(anal1_res_table[[anal1_col]])
                 #     anal1_trans = input$combi_anal1_trans
                 #   })
                 # }
                 # if(input$combi_anal2_trans != "none"){
                 #   try({
                 #     transFun= switch(input$combi_anal2_trans,
                 #                      "log10"=log10,
                 #                      "-log10"=function(x) -log10(x),
                 #                      "abs"=abs) 
                 #     anal2_res_table[[anal2_col]] <- transFun(anal2_res_table[[anal2_col]])
                 #     anal2_trans = input$combi_anal2_trans
                 #   })
                 # }
                 
                 mzInBoth = intersect(rownames(anal1_res_table),rownames(anal2_res_table))
                 if(length(mzInBoth) > 0){
                   combined_anal = merge(
                     anal1_res_table,
                     anal2_res_table,
                     by = "rn"
                   )
                   keep.cols = c(1, which(colnames(combined_anal) %in% c(anal1_col,anal2_col)))
                   dt <- as.data.frame(combined_anal)[,keep.cols]
                   mSet$analSet$combi <- list(sig.mat = dt, 
                                              all.vals = list(x=anal1_all_res, y=anal2_all_res),
                                              trans = list(x=anal1_trans, y=anal2_trans),
                                              source = list(x=anal1, y=anal2))
                 }
                 
               },
               volcano = {
                 shiny::withProgress({
                   mSet <- MetaboAnalystR::Volcano.Anal(mSetObj = mSet,
                                                        paired = mSet$dataSet$paired, 
                                                        fcthresh = 1.1, cmpType = 0,
                                                        percent.thresh = 0.75,
                                                        nonpar = F, threshp = 0.1,
                                                        equal.var = TRUE,
                                                        pval.type = "fdr") # TODO: make thresholds user-defined
                 })
               },
               tsne = {
                 shiny::withProgress({
                   inTbl = switch(input$tsne_source, 
                                  original = mSet$dataSet$prog,
                                  "pre-batch correction" = mSet$dataSet$prebatch,
                                  normalized = mSet$dataSet$norm)
                   coords = tsne::tsne(inTbl, k = 3,
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
        save_info$has_changed <- TRUE
        shinyjs::show(selector = paste0("div.panel[value=collapse_", statsmanager$calculate, "_plots]"))
        shinyjs::show(selector = paste0("div.panel[value=collapse_", statsmanager$calculate, "_tables]"))
        shinyBS::updateCollapse(session, paste0("collapse_",input$statistics),open = paste0("collapse_", 
                                                                                            statsmanager$calculate, 
                                                                                            c("_tables","_plots")))
        if(lcl$beep){
          beepr::beep(sound = lcl$aes$which_beep)
          Sys.sleep(0.6)
          beepr::beep(sound = lcl$aes$which_beep)
          Sys.sleep(0.6)
          beepr::beep(sound = lcl$aes$which_beep)
        }
        
      }else{
        MetaboShiny::metshiAlert("Analysis failed!")
        #shiny::showNotification(msg.vec)
        mSet <<- mSet.old
      }
      lcl <<- lcl
    }
    # - - - -
    statsmanager$calculate <- NULL # set reloading to 'off'
  }
})