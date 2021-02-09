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
               corr = {
                 # pearson kendall spearman
                 lvls = levels(mSet$dataSet$cls)
                 pat = input$corr_seq_order
                 pat_order = match(lvls,pat)
                 pattern = paste0(pat_order-1, collapse="-")
                 mSet <- MetaboAnalystR::Match.Pattern(mSet, input$corr_corr, pattern)
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
                   shiny::setProgress(value = 0)
                   
                   pickedTbl <- if(input$ml_run_on_norm) "norm" else "orig"
                   
                   # get base table to use for process
                   curr <- data.table::as.data.table(mSet$dataSet[[pickedTbl]])
                   
                   if(input$ml_specific_mzs != "no"){
                     shiny::showNotification("Using user-specified m/z set.")
                     if(!is.null(input$ml_mzs)){
                       curr <- curr[,input$ml_mzs, with=F]
                     }else{
                       mzs = getTopHits(mSet, 
                                        input$ml_specific_mzs, 
                                        input$ml_mzs_topn)[[1]]
                       mzs = gsub("^X", "", mzs)
                       mzs = gsub("\\.$", "-", mzs)
                       curr <- curr[,..mzs]
                     }
                   }
                   
                   # replace NA's with zero
                   for (j in seq_len(ncol(curr))){
                     set(curr,which(is.na(curr[[j]])),j,0)
                   }
                   
                   # conv to data frame
                   curr <- as.data.frame(curr)
                   rownames(curr) <- rownames(mSet$dataSet[[pickedTbl]])
                   
                   # find the qc rows and remove them
                   is.qc <- grepl("QC|qc", rownames(curr))
                   if(sum(is.qc) > 0){
                     curr <- curr[!is.qc,]
                   }
                   
                   order <- match(rownames(curr), mSet$dataSet$covars$sample)
                   
                   if("label" %in% colnames(mSet$dataSet$covars)){
                     config <- mSet$dataSet$covars[order, -"label"]
                   }else{
                     config <- mSet$dataSet$covars[order, ]
                   }
                   
                   if(!is.null(lcl$vectors$ml_train)){
                     if(unique(lcl$vectors$ml_train) %in% c("split","all")){
                       lcl$vectors$ml_train <- NULL
                     }
                   }
                   if(!is.null(lcl$vectors$ml_test)){
                     if(unique(lcl$vectors$ml_test)%in% c("split","all")){
                       lcl$vectors$ml_test <- NULL
                     }
                   }
                   
                   batch_sampling = input$ml_batch_sampling
                   batches = input$ml_batch_covars
                   
                   if(input$ml_samp_distr != " "){
                     spl.name = stringr::str_split(input$ml_samp_distr, " - ")[[1]]
                     ml.method = spl.name[1]
                     ml.name = spl.name[2]
                     ml.anal = mSet$analSet$ml[[ml.method]][[ml.name]]
                     if("distr" %in% names(ml.anal)){
                       batch_sampling <- "none"
                       distrs = unique(ml.anal$distr)
                       if(length(distrs) > 1){
                         shiny::showNotification("You selected a model with multiple repeats with random split every repeat!
                                               Will only use the first.")
                       }
                       ml.dist = distrs[[1]]
                       vec = rep("train", nrow(mSet$dataSet$covars))
                       vec[as.numeric(ml.dist$test)] <- "test"
                       config$split <- vec
                     }else{
                       shiny::showNotification("Not available in this model (pre-update). Defaulting to selected settings!")
                     }
                   }
                   
                   needed_for_subset <- c(lcl$vectors$ml_train[1],
                                          lcl$vectors$ml_test[1])
                   
                   label = unlist(config[,mSet$settings$exp.var,with=F])
                   
                   if(input$ml_label_shuffle){
                     print("Shuffling labels...")
                     label = sample(label)
                   }
                   
                   config <- config[, unique(c(input$ml_include_covars,
                                               if(length(batches)>0) input$ml_batch_covars else c(),
                                               needed_for_subset)),with=F]# reorder so both halves match up later  
                   
                   
                   if(mSet$settings$exp.type %in% c("2f", "t1f")){
                     cls <- if(input$ml_run_on_norm) "" else "orig."
                     
                     # just set to facA for now..
                     if(nrow(config)==0){
                       config <- data.frame(label=mSet$dataSet[[paste0(cls,"facA")]])
                     }else{
                       config <- cbind(config, label=mSet$dataSet[[paste0(cls,"facA")]]) # add current experimental condition
                     }
                   }else{
                     cls <- "cls" #if(input$ml_run_on_norm) "cls" else "orig.cls"
                     if(nrow(config)==0){
                       config <- data.frame(label=label)
                     }else{
                       config <- cbind(config, 
                                       label=label) # add current experimental condition
                     }
                   }
                   
                   # TRAIN/TEST SPLIT
                   t = as.data.table(cbind(config, curr))
                   
                   if(!is.null(lcl$vectors$ml_train) | !is.null(lcl$vectors$ml_test)){
                     # add clause for same train_test
                     test_idx = NULL
                     train_idx = NULL
                     if(!is.null(lcl$vectors$ml_test)){
                       test_idx = which(t[[lcl$vectors$ml_test[1]]] == lcl$vectors$ml_test[2])
                     }
                     if(!is.null(lcl$vectors$ml_train)){
                       train_idx = which(t[[lcl$vectors$ml_train[1]]] == lcl$vectors$ml_train[2])
                     }
                     if(is.null(train_idx)){
                       train_idx = setdiff(1:nrow(t), test_idx)  
                     }else if(is.null(test_idx)){
                       test_idx = setdiff(1:nrow(t), train_idx)
                     }
                   }else{
                     split_label = if(length(batches) == 1){
                       print("splitting tr/te % per batch")
                       paste0(t$label,"AND",t[,..batches][[1]])
                     }else{
                       print("unbiased split over pool")
                       t$label
                     }
                     train_idx = caret::createDataPartition(y = split_label, 
                                                            p = input$ml_train_perc/100,
                                                            list = FALSE)[,1] # partition data in a balanced way (uses labels)
                     test_idx = setdiff(1:nrow(t), train_idx)
                   }
                   
                   test_sampnames = rownames(curr)[test_idx]
                   
                   if(length(batches) > 0){
                     batch_sizes = table(t[,..batches][[1]])
                     size_per_group = if(batch_sampling == "down") min(batch_sizes) else max(batch_sizes)
                     n_batches = length(unique(t[,..batches][[1]]))
                     #size_per_group = round(size_per_group / length(unique(config$label)), digits = 0)
                   }
                   
                   training_data = t[train_idx,]
                   testing_data = t[-train_idx,]
                   
                   configCols = 1:ncol(config)
                   mzCols = setdiff(1:ncol(training_data), configCols)
                   
                   
                   ## ONLY APPLY TO TRAINING
                   ## IF BATCHES ARE PRESENT 
                   if(length(batches) > 0 & batch_sampling != "none"){
                     
                     batch_is_label = all(t[[batches]] == t$label)
                     
                     #if(batch_sampling != "none"){
                     print("resampling classes per batch group")
                     # RESAMPLE BASED ON BATCH VARIABLE
                     mzs=colnames(training_data)[mzCols]
                     colnames(training_data)[mzCols]=paste0("mz",1:length(mzCols))
                     labelCol = which(colnames(training_data) == "label")
                     
                     spl.t = if(!batch_is_label) split(training_data, training_data[,..batches]) else list(training_data)
                     
                     if(input$ml_batch_size_sampling){
                       #### RESAMPLE TO EQUALIZE BATCH CATEGORY SIZE (ALL COUNTRIES IN TRAIN SAME SIZE EACH)
                       library(plyr)
                       balanced.spl.t <- lapply(spl.t, function(l){
                         r = switch(batch_sampling,
                                    up = upsample.adj(l, as.factor(l$label), maxClass = size_per_group),
                                    rose = {
                                      resampled = ROSE::ROSE(label ~ ., 
                                                             data = {
                                                               dat = l[, c(labelCol, mzCols), with=F]
                                                               cols = colnames(dat)[2:ncol(dat)]
                                                               dat[, (cols) := lapply(.SD, as.numeric), 
                                                                   .SDcols = cols]
                                                               dat
                                                             }, 
                                                             N = size_per_group * n_batches) # TODO MAKE THIS LEVELS OF CONFIG
                                      resampled.data = resampled$data
                                      colnames(resampled.data) <- c("label", mzs)
                                      keep.config = sapply(configCols, function(i){
                                        length(unique(l[[i]] )) == 1
                                      })
                                      resampled.config = data.table::as.data.table(lapply(which(keep.config), 
                                                                                          function(i){
                                                                                            c(rep(unique(l[[i]]), 
                                                                                                  nrow(resampled.data)))
                                                                                          }))
                                      colnames(resampled.config) <- colnames(l)[which(keep.config)]
                                      joined.data = cbind(resampled.config, resampled.data)
                                      joined.data
                                    },
                                    down = downsample.adj(l, as.factor(l$label), minClass = size_per_group))
                        as.data.table(r)   
                       })
                       r = rbindlist(balanced.spl.t)
                       
                       if("Class" %in% names(r)){
                         r$Class <- NULL  
                       }
                       
                     }else if(length(batches)==1){
                       
                       ### BALANCE LABELS WITHIN BATCHES, BUT DON'T MAKE THEM ALL THE SAME SIZE
                       balanced.spl.t <- lapply(spl.t, function(l){
                         size_per_group = if(batch_sampling == "down") min(table(l$label)) else max(table(l$label))
                         r = switch(batch_sampling,
                                    up = upsample.adj(l, as.factor(l$label), maxClass = size_per_group),
                                    rose = {
                                      resampled = ROSE::ROSE(label ~ ., 
                                                             data = {
                                                               dat = l[, c(labelCol, mzCols), with=F]
                                                               cols = colnames(dat)[2:ncol(dat)]
                                                               dat[, (cols) := lapply(.SD, as.numeric), 
                                                                   .SDcols = cols]
                                                               dat
                                                             }, 
                                                             N = size_per_group)
                                      resampled.data = resampled$data
                                      colnames(resampled.data) <- c("label", mzs)
                                      keep.config = sapply(configCols, function(i){
                                        length(unique(l[[i]] )) == 1
                                      })
                                      resampled.config = data.table::as.data.table(lapply(which(keep.config), 
                                                                                          function(i){
                                                                                            c(rep(unique(l[[i]]), 
                                                                                                  nrow(resampled.data)))
                                                                                          }))
                                      colnames(resampled.config) <- colnames(l)[which(keep.config)]
                                      joined.data = cbind(resampled.config, resampled.data)
                                      joined.data
                                    },
                                    down = downsample.adj(l, as.factor(l$label), minClass = size_per_group))
                         if("Class" %in% names(r)){
                           r$Class <- NULL  
                         }
                         as.data.table(r) 
                       })
                       r = rbindlist(balanced.spl.t)
                     }else{
                       r = as.data.table(cbind(config, curr))
                     }
                     configCols = colnames(config)
                     configCols = intersect(configCols, colnames(r))
                     config_adj = r[, ..configCols]
                     training_data_adj = r[, -..configCols]
                     colnames(training_data_adj) <- mzs
                     training_cleaned = cbind(config_adj, training_data_adj)
                     
                     testing_data_adj = cbind(split = "test", testing_data)
                     training_data_cleaned_adj = cbind(split = "train", training_cleaned)
                     
                     configCols = c(configCols, "split")
                     curr = rbind(testing_data_adj, 
                                  training_data_cleaned_adj)
                     config = curr[,..configCols]
                     curr = curr[,-configCols,with=F]
                     
                   }else{
                     if(input$ml_sampling != "none"){
                       print('resampling based on class...')
                       l = training_data
                       mzs = colnames(l)[mzCols]
                       colnames(l)[mzCols]=paste0("mz",1:length(mzCols))
                       labelCol = which(colnames(training_data) == "label")
                       size_per_group = if(input$ml_sampling == "down") min(table(l[[labelCol]])) else max(table(l[[labelCol]]))
                       r = switch(input$ml_sampling,
                                  up = upsample.adj(l, as.factor(l$label), maxClass = size_per_group),
                                  rose = {
                                    resampled = ROSE::ROSE(label ~ ., 
                                                           data = {
                                                             dat = l[, c(labelCol, mzCols), with=F]
                                                             cols = colnames(dat)[2:ncol(dat)]
                                                             dat[, (cols) := lapply(.SD, as.numeric), 
                                                                 .SDcols = cols]
                                                             dat
                                                           }, 
                                                           N = size_per_group * 2)
                                    resampled.data = resampled$data
                                    colnames(resampled.data) <- c("label", mzs)
                                    keep.config = sapply(configCols, function(i){
                                      length(unique(l[[i]] )) == 1
                                    })
                                    resampled.config = data.table::as.data.table(lapply(which(keep.config), 
                                                                                        function(i){
                                                                                          c(rep(unique(l[[i]]), 
                                                                                                nrow(resampled.data)))
                                                                                        }))
                                    colnames(resampled.config) <- colnames(l)[which(keep.config)]
                                    joined.data = cbind(resampled.config, resampled.data)
                                    joined.data
                                  },
                                  down = downsample.adj(l, as.factor(l$label), minClass = size_per_group))
                       configCols = colnames(config)
                       configCols = intersect(configCols, colnames(r))
                       r = as.data.table(r)
                       if("Class" %in% names(r)){
                         r$Class <- NULL  
                       }
                       config_adj = r[, ..configCols]
                       training_data_adj = r[, -..configCols]
                       colnames(training_data_adj) <- mzs
                       training_data = cbind(config_adj, training_data_adj)
                     }
                     
                     # ========
                     train.config = training_data[,..configCols]
                     train.mzs = training_data[,..mzCols]
                     test.config = testing_data[,..configCols]
                     test.mzs = testing_data[,..mzCols]
                     test.config$split = "test"
                     train.config$split = "train"
                     # add the right columns and define curr
                     curr = rbind(test.mzs, 
                                  train.mzs)
                     config = rbind(test.config, 
                                    train.config)
                     removeCols = setdiff(needed_for_subset, input$ml_include_covars)
                     
                     folds = NULL
                     if(length(batches) == 1){
                       if(input$ml_folds != "LOOCV"){
                         if(as.numeric(input$ml_folds) <= 5){
                           group_fac = config[split == "train",][[batches]]
                           ml_folds = min(c(as.numeric(input$ml_folds), length(unique(group_fac))))
                           folds = caret::groupKFold(group_fac, k = ml_folds)
                         }
                       }  
                     }
                     
                     if(length(removeCols) > 0){
                       curr = curr[,-removeCols,with=F]
                       config = config[,-removeCols, with=F]
                     }
                   }
                   
                   lcl$vectors$ml_train = c("split","train")
                   lcl$vectors$ml_test = c("split","test")
                   
                   config <- droplevels(config)
                   config <- data.table::as.data.table(config)
                   config <- config[,apply(!is.na(config), 2, any), with=FALSE]
                   
                   curr <- as.data.table(curr)
                   
                   predictor <- config$label
                   predict_idx <- which(colnames(config) == "label")
                   exact_matches <- which(unlist(lapply(config, function(col) all(as.numeric(as.factor(col)) == as.numeric(as.factor(predictor))))))
                   remove <- setdiff(exact_matches, predict_idx)
                   
                   # # remove ones with na present
                   has.na <- apply(config, MARGIN=2, FUN=function(x) any(is.na(x) | tolower(x) == "unknown"))
                   has.all.unique <- apply(config, MARGIN=2, FUN=function(x) length(unique(x)) == length(x))
                   remove = colnames(config)[which(has.na | has.all.unique)]
                   
                   outersect = function(a,b) setdiff(union(a,b), intersect(a,b))
                   
                   if(length(batches) > 0 ){
                     alsoRemove = input$ml_batch_covars[sapply(input$ml_batch_covars, function(b)
                       !(b %in% input$ml_include_covars))]
                     
                     remove = c(remove, alsoRemove)
                   }
                   
                   #keep_configs <- which(names(config) == "label")
                   remove <- unique(c(remove, 
                                      "sample",  
                                      "individual", 
                                      colnames(config)[caret::nearZeroVar(config)]))
                   
                   remove = setdiff(remove, "label")
                   remove = setdiff(remove, "split")
                   
                   keep_configs <- which(!(colnames(config) %in% remove))
                   
                   try({
                     msg = paste0("Keeping non-mz variables after NA/unique filtering: ",
                                  paste0(names(config)[keep_configs],
                                         collapse = ", "))
                     print(msg)
                     shiny::showNotification(msg)
                   })
                   
                   # remove duplicated config columns
                   
                   config = config[,!duplicated(colnames(config)),with=F]
                   
                   # ================================
                   
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
                   
                   levels(curr$label) <- paste0("class",  Hmisc::capitalize(levels(curr$label)))
                   levels(curr$label) <- ordered(levels(curr$label))
                   
                   # ============ LOOP HERE ============
                   
                   #colnames(curr) <- paste0("X", gsub("\\-","min",colnames(curr)))
                   colnames(curr) <- make.names(colnames(curr))
                   
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
                                                           sampling,
                                                           batch_sampling,
                                                           batches,
                                                           folds){
                                                    runML(curr,
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
                                                          sampling = sampling,
                                                          batch_sampling = batch_sampling,
                                                          batches = batches,
                                                          folds = folds)
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
                                                  sampling = NULL,#if(input$ml_sampling == "none") NULL else input$ml_sampling,
                                                  batch_sampling = if(input$ml_batch_sampling == "none") NULL else input$ml_batch_sampling,
                                                  batches = input$ml_batch_covars,
                                                  folds = folds
                     )
                   })
                   
                   # check if a storage list for machine learning results already exists
                   if(!("ml" %in% names(mSet$analSet))){
                     mSet$analSet$ml <- list() # otherwise make it
                   }
                   
                   # ---------
                   
                   labels = lapply(repeats, function(x) x$labels)
                   predictions = lapply(repeats, function(x) x$prediction)
                   
                   perf.per.curve = lapply(c("roc", "precrec"), function(curve_type){
                     if(length(levels(mSet$dataSet$cls)) > 2){
                       #curve_type = "roc" # don't have precrec here implemented yet
                       perf <- lapply(1:length(labels), function(i){
                         prediction = predictions[[i]]
                         classes = labels[[i]]
                         rmv.classes = setdiff(levels(classes), unique(as.character(classes)))
                         if(length(rmv.classes) > 0){
                           prediction = prediction[,!(colnames(prediction) %in% rmv.classes)]
                           classes = droplevels(classes)
                         }
                         #x = list(prediction=predictions[[i]],
                         #         labels=as.factor(labels[[i]]))
                         colnames(prediction) <- paste0(colnames(prediction), "_pred_m1")
                         samprows = lapply(1:length(classes),function(i){
                           default = rep(0, length(levels(classes)))
                           default = data.table::as.data.table(t(data.table::data.table(default)))
                           colnames(default) = paste0(levels(classes), "_true")
                           lbl = as.character(classes[i])
                           match_lbl = paste0(lbl,"_true") 
                           default[, match_lbl] <- 1
                           default
                         })
                         true_cols = data.table::rbindlist(samprows)
                         my_data = cbind(true_cols, prediction)
                         multi_performance <- if(curve_type == "roc") multiROC::multi_roc(my_data) else multiROC::multi_pr(my_data)
                         plottable = if(curve_type == "roc"){
                           data.table::as.data.table(multiROC::plot_roc_data(multi_performance))
                         }else{
                           data.table::as.data.table(multiROC::plot_pr_data(multi_performance)) 
                         }
                         res = plottable[Method == "m1", -"Method"]
                         #res = getMultiMLperformance(x)
                         res$attempt = c(i)
                         res
                       })
                       
                       perf.long <- data.table::rbindlist(perf)
                       perf.long <- perf.long[!(Group %in% c("Micro", "Macro")),]
                       
                       #TPR is sensitivity
                       #FPR = 1 - specificity
                       if(curve_type == "roc"){
                         TPR = perf.long$Sensitivity
                         FPR = 1 - perf.long$Specificity
                         perf.long[[1]] <- FPR
                         perf.long[[2]] <- TPR
                       }
                       colnames(perf.long) <- c("x", "y", 
                                                "comparison", 
                                                "AUC_PAIR", 
                                                "attempt")
                       mean.auc <- mean(unique(perf.long$AUC_PAIR))
                       perf.long$AUC_AVG <- c(mean.auc)
                     }else{
                       # save the summary of all repeats (will be used in plots) TOO MEMORY HEAVY
                       pred <- ROCR::prediction(predictions = lapply(predictions, function(x) x[[2]]), 
                                                labels = lapply(labels, function(x) x))
                       y = switch(curve_type, roc = "tpr", precrec = "prec")
                       x = switch(curve_type, roc = "fpr", precrec = "rec")
                       
                       perf <- ROCR::performance(prediction.obj = pred, 
                                                 measure = y, 
                                                 x.measure = x)
                       perf_auc <- ROCR::performance(pred, switch(curve_type, 
                                                                  roc = "auc", 
                                                                  precrec = "aucpr"))
                       perf.long <- data.table::rbindlist(lapply(1:length(perf@x.values), function(i){
                         xvals <- perf@x.values[[i]]
                         yvals <- perf@y.values[[i]]
                         aucs <- signif(perf_auc@y.values[[i]][[1]], digits = 2)
                         res = data.table::data.table(attempt = c(i),
                                                      x = xvals,
                                                      y = yvals,
                                                      AUC_AVG = aucs,
                                                      AUC_PAIR = aucs)
                         res
                       }))
                       perf.long$comparison <- paste0(levels(mSet$dataSet$cls),collapse=" vs. ")
                       mean.auc <- mean(unlist(perf_auc@y.values))
                       perf.long$AUC_AVG <- c(mean.auc)
                     }
                     perf.long$metric = curve_type
                     perf.long
                   })
                   
                   perf.long <- data.table::rbindlist(perf.per.curve)
                   
                   # ---------
                   
                   roc_data <- list(perf = perf.long,
                                    labels = lapply(repeats, function(x) x$labels),
                                    predictions = lapply(repeats, function(x) x$prediction),
                                    imp = lapply(repeats, function(x) x$importance),
                                    inTest = test_sampnames)
                   
                   bar_data <- data.table::rbindlist(lapply(1:length(repeats), function(i){
                     x = repeats[[i]]
                     tbl = data.table::as.data.table(x$importance, keep.rownames=T)[,1:2]
                     tbl$rep = c(i)
                     colnames(tbl) = c("mz",
                                       "importance",
                                       "rep")
                     tbl$mz <- gsub("^X","",tbl$mz)
                     tbl$mz <- gsub("\\.$","-",tbl$mz)
                     
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
                   if(input$ml_method %not in% names(mSet$analSet$ml)){
                     mSet$analSet$ml[[input$ml_method]] <- list()
                   }
                   
                   mSet$analSet$ml[[input$ml_method]][[input$ml_name]] <- list("roc" = roc_data,
                                                                               "bar" = bar_data,
                                                                               "distr" = lapply(repeats, function(x) x$distr))
                   mSet$analSet$ml$last <- list(name = input$ml_name,
                                                method = input$ml_method)
                   
                   lcl$vectors$ml_train <- lcl$vectors$ml_train <<- NULL
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
                 
                 anal1_trans = "none"
                 if(input$combi_anal1_trans != "none"){
                   try({
                     transFun= switch(input$combi_anal1_trans,
                                      "log10"=log10,
                                      "-log10"=function(x) -log10(x),
                                      "abs"=abs) 
                     anal1_res_table[[anal1_col]] <- transFun(anal1_res_table[[anal1_col]])
                     anal1_trans = input$combi_anal1_trans
                   })
                 }
                 
                 anal2_trans = "none"
                 if(input$combi_anal2_trans != "none"){
                   try({
                     transFun= switch(input$combi_anal2_trans,
                                      "log10"=log10,
                                      "-log10"=function(x) -log10(x),
                                      "abs"=abs) 
                     anal2_res_table[[anal2_col]] <- transFun(anal2_res_table[[anal2_col]])
                     anal2_trans = input$combi_anal2_trans
                   })
                 }
                 
                 mzInBoth = intersect(rownames(anal1_res_table),rownames(anal2_res_table))
                 if(length(mzInBoth) > 0){
                   combined_anal = merge(
                     anal1_res_table,
                     anal2_res_table,
                     by = "rn"
                   )
                   keep.cols = c(1, which(colnames(combined_anal) %in% c(anal1_col,anal2_col)))
                   dt <- as.data.frame(combined_anal)[,keep.cols]
                   mSet$analSet$combi <<- list(sig.mat = dt, 
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