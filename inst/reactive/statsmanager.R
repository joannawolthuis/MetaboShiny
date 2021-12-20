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
                 mSet <- MetaboShiny::store.mSet(mSet, name = mSet$settings$cls.name, proj.folder = lcl$paths$proj_dir)
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
                   flattened <- list(getAllHits(mSet,
                                                input$mummi_anal))
                   if(input$mummi_enr_method & grepl("combi", input$mummi_anal)){
                     head(flattened)
                     flattened[[1]]$value <- 1
                     flattened[[1]]$significant <- F
                     flattened[[1]] <- flattened[[1]][order(abs(flattened[[1]]$significant),decreasing = T),]
                     flattened[[1]]$value[1:input$mummi_topn] <- 0
                     flattened[[1]]$significant[1:input$mummi_topn] <- T
                   }
                   
                   hasP = T#grepl("tt|aov|asca|combi|venn",input$mummi_anal)
                   setProgress(0.1)
                   
                   myFile <- tempfile(fileext = ".csv")
                   tbl = data.table::data.table("m.z" = as.numeric(gsub(flattened[[1]]$m.z, pattern="(\\+|\\-|RT).*$", replacement="")),
                                                mode = sapply(flattened[[1]]$m.z, function(mz){
                                                  if(grepl(pattern="-",x=mz)) "negative" else "positive"
                                                }))
                   
                   
                   hasT = ncol(flattened[[1]]) >= 3
                   
                   anal = gsub(" \\(.*$|", "", input$mummi_anal)
                   subset = gsub("\\(|\\)|.*\\(", "", input$mummi_anal)
                   
                   tbl[, "p.value"] = if(hasP) flattened[[1]][,2] else c(0)
                   tbl[, "t.score"] = if(hasT) flattened[[1]][,3] else c(NA)
                   
                   if(hasP) if(all(is.na(tbl$p.value))) tbl$p.value <- c(0)
                   if(hasT) if(all(is.na(tbl$t.score))) tbl$t.score <- c(0)
                   
                   tmpfile <- tempfile()
                   
                   print("Preview of input table:")
                   print(head(tbl))
                   
                   fwrite(if(hasT) tbl else tbl[,1:3], file=tmpfile)
                   
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
                   
                   enr_mSet <- MetaboAnalystR::SetPeakEnrichMethod(enr_mSet, if(input$mummi_enr_method | !hasT) "mum" else "gsea", "v2")
                   enr_mSet <- MetaboAnalystR::SetMummichogPval(enr_mSet, if(hasP) as.numeric(gsub(",",".",input$mummi_pval)) else 0.9)
                   
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
                     row$Charge = as.numeric(row$Charge)
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
                   
                   # ==== BUILD CUSTOM DB HERE ====
                   
                   map_id = input$mummi_org
                   
                   print(map_id)
                   
                   lib_name <- paste0(map_id, "_kegg")
                   file_name <- paste0(lib_name, ".qs")
                   
                   if(!file.exists(file_name)){
                     mummichog.lib <- build.enrich.KEGG(map_id)
                     qs::qsave(mummichog.lib, file = file_name)
                   }
                   
                   mummichog.lib = qs::qread(file_name)
                   
                   mummi.adducts <- new_adduct_mzlist(enr_mSet,
                                                      mw = mummichog.lib$cpd.lib$mw)
                   
                   mummi.adducts <- list(dpj_positive = mummi.adducts$pos,
                                         positive = mummi.adducts$pos,
                                         negative = mummi.adducts$neg)
                   
                   mummichog.lib$cpd.lib$adducts <- mummi.adducts
                   cpd.tree <- list()
                   ms_modes <- c("dpj_positive", "positive", "negative")
                   
                   for (ms_mode in ms_modes) {
                     l2 <- list()
                     l2[[49]] <- ""
                     l2[[2001]] <- ""
                     mz.mat <- mummichog.lib$cpd.lib$adducts[[ms_mode]]
                     floor.mzs <- floor(mz.mat)
                     for (i in 1:nrow(floor.mzs)) {
                       neighbourhood <- floor.mzs[i, ]
                       for (n in neighbourhood) {
                         if ((n > 50) & (n < 2000)) {
                           l2[[n]] <- append(l2[[n]], i)
                         }
                       }
                     }
                     cpd.tree[[ms_mode]] <- lapply(l2, unique)
                   }
                   
                   org = map_id
                   mummichog.lib$cpd.tree <- cpd.tree
                   qs::qsave(mummichog.lib, file = file_name)
                   
                   print(paste0(org, " mummichog library updated with chosen adducts!"))
                   
                   # ==============================
                   
                   enr_mSet <- MetaboAnalystR::PerformPSEA(mSetObj = enr_mSet, 
                                                           lib = lib_name,
                                                           libVersion = "current",
                                                           permNum = 100)                      
                   
                   filenm <- "mummichog_matched_compound_all.csv"
                   enr_mSet$dataSet$mumResTable <- data.table::fread(filenm)
                   
                   flattened[[1]]$significant = flattened[[1]]$value < if(hasP) as.numeric(gsub(",",".",input$mummi_pval)) else 0.9
                   
                   # names to compounds
                   cpd2name <- data.table::data.table(Matched.Compound = mummichog.lib$cpd.lib$id,
                                                      Compound.Name = mummichog.lib$cpd.lib$name)
                   enr_mSet$dataSet$mumResTable = merge(cpd2name, enr_mSet$dataSet$mumResTable)
                   
                   mSet$analSet$enrich <- list(mummi.resmat = enr_mSet$mummi.resmat,
                                               mummi.gsea.resmat = enr_mSet$mummi.gsea.resmat,
                                               mumResTable = enr_mSet$dataSet$mumResTable,
                                               mummi.input = enr_mSet$dataSet$mummi.proc,
                                               path.nms = enr_mSet$path.nms,
                                               path.hits = enr_mSet$path.hits,
                                               path.all = enr_mSet$pathways,
                                               path.lib = enr_mSet$lib.organism,
                                               cpd.value = enr_mSet$cpd_exp_dict,
                                               orig.input = flattened[[1]],
                                               enr.method = if(input$mummi_enr_method | !hasT) "mum" else "gsea")
                   enr_mSet <- NULL
                 })
               },
               featsel = {
                 print("Feature selection start...")
                 curr = cbind(label = mSet$dataSet$cls, 
                              mSet$dataSet$norm)
                 boruta_res = Boruta::Boruta(x = curr[,2:ncol(curr)],
                                             y = as.factor(curr[[1]]))
                 mSet$analSet$featsel <- list(boruta_res)
               },
               ml = {
                   {
                     assign("ml_queue", shiny::isolate(shiny::reactiveValuesToList(ml_queue)), envir = .GlobalEnv)
                     assign("input", shiny::isolate(shiny::reactiveValuesToList(input)), envir = .GlobalEnv)
                     
                     # make subsetted mset for ML so it's not as huge in memory
                     small_mSet <- list(metshiParams=mSet$metshiParams)
                     small_mSet$dataSet <- mSet$dataSet[c("cls", "orig.cls", 
                                                          "orig", "norm", 
                                                          "covars")]
                     
                     ml_queue$jobs <- ml_queue$jobs[!(names(ml_queue$jobs) %in% mSet$analSet$ml$rf)]
                     
                     #ml_queue$jobs = ml_queue$jobs[grep("sa", names(ml_queue$jobs))]
                     resource_saver_mode = F # only make the datasets once, for example with resampling
                     if(resource_saver_mode){
                       vars.require.dataset.change <- c("ml_test_subset",
                                                        "ml_train_subset",
                                                        "ml_sampling",
                                                        "ml_batch_balance",
                                                        "ml_batch_covars",
                                                        "ml_batch_size_sampling",
                                                        "ml_groupsize",
                                                        "ml_include_covars",
                                                        "ml_used_table",
                                                        "ml_pca_corr",
                                                        "ml_train_perc",
                                                        "ml_keep_pcs",
                                                        "ml_pca_corr"
                       )
                       rows.changes.dataset <- lapply(ml_queue$jobs, function(job){
                         job[vars.require.dataset.change]
                       })
                       unique.dataset.change.jobs = unique(rows.changes.dataset)
                       print(paste("preparing", length(unique.dataset.change.jobs), "dataset(s) for jobs"))
                       train.test.unique <- pbapply::pblapply(unique.dataset.change.jobs, 
                                                              cl=0, # TODO: make parallel
                                                              function(job){
                                                                tr_te = ml_prep_data(settings = job, 
                                                                                     mSet = small_mSet,
                                                                                     input = input, cl=0)
                                                              })
                       
                       mapper = data.table::rbindlist(lapply(1:length(rows.changes.dataset), function(i){
                         job = rows.changes.dataset[[i]]
                         jobi = which(sapply(unique.dataset.change.jobs, 
                                             function(uniq.job) identical(job, uniq.job)))
                         data.table::data.table(ml_name = ml_queue$jobs[[i]]$ml_name, unique_data_id=jobi)
                       }))
                       small_mSet$dataSet$for_ml <- list(datasets = train.test.unique, 
                                                         mapper = mapper)
                     }
                     
                     uses.specific.mzs <- any(sapply(ml_queue$jobs, function(settings) settings$ml_specific_mzs != "no"))
                     if(uses.specific.mzs){
                       keep.analyses <- gsub(" \\(.*$", "", sapply(ml_queue$jobs, function(settings) settings$ml_specific_mzs))
                       keep.analyses <- unique(keep.analyses[keep.analyses != "no"])
                       #if("pca" %in% names(small_mSet$analSet)) keep.analyses <- unique(c("pca", keep.analyses))
                       small_mSet$analSet <- mSet$analSet[keep.analyses]
                       small_mSet$storage <- lapply(mSet$storage, function(store){
                         store$analSet <- store$analSet[keep.analyses]
                         store
                       })
                     }else{
                       small_mSet$analSet <- NULL
                       small_mSet$storage <- list()
                     }
                     
                     has_slurm = Sys.getenv("SLURM_CPUS_ON_NODE") != ""
                     use_slurm = T
                     
                     net_cores = input$ncores# - 1
                     if(net_cores > 0 & (!has_slurm | !use_slurm)){
                       try({
                         parallel::stopCluster(session_cl)
                         parallel::stopCluster(ml_session_cl)
                       })
                       logfile <- tempfile()
                       print("Log file at:")
                       print(logfile)
                       #if(file.exists(logfile)) file.remove(logfile)
                       ml_session_cl <- parallel::makeCluster(net_cores, outfile=logfile)#,setup_strategy = "sequential") # leave 1 core for general use and 1 core for shiny session
                       # send specific functions/packages to other threads
                       parallel::clusterEvalQ(ml_session_cl, {
                         library(data.table)
                         library(iterators)
                         library(MetaboShiny)
                         library(MetaDBparse)
                       })  
                     }else{
                       ml_session_cl <- 0
                     }
                     
                     mSet_loc <- tempfile()
                     qs::qsave(small_mSet, mSet_loc)

                     if(has_slurm & use_slurm){
                       
                       job_time = "24:00:00"
                       
                       print("Attempting to submit jobs through slurm!")
                       
                       settings_loc <- tempfile()
                       first_job_parallel_count = if(ml_queue$jobs[[1]]$ml_label_shuffle){
                         ml_queue$jobs[[1]]$ml_n_shufflings + 1
                       }else{
                         1
                       }
                       qs::qsave(ml_queue$jobs, file = settings_loc)
                       
                       pars = data.frame(i = 1:length(ml_queue$jobs),
                                         mloc = c(mSet_loc),
                                         sloc = c(settings_loc),
                                         tmpdir = c(dirname(tempfile())))
                       
                       print(head(pars))
                       success = F
                       
                       try({
                         fn <- normalizePath(paste0(tools::file_path_sans_ext(lcl$paths$csv_loc), 
                                                    ".metshi"))
                         qs::qsave(mSet, file = fn)
                
                         mem_gb = "50G"
                         pars_filt = pars#[9981:10000,]
                         jobname=paste0("METSHI_ML_",
                                        lcl$proj_name,
                                        "_",
                                        gsub("file","",
                                             basename(tempfile())))
                         
                         maxjobs = 500
                         nodecount = floor(maxjobs / first_job_parallel_count)
                         nodecount =  min(nodecount, nrow(pars_filt))
                         
                         batch_job <- slurm_apply_metshi(ml_slurm, 
                                                         pars_filt,#[5000,], 
                                                         cpus_per_node = 1,
                                                         jobname = jobname,
                                                         nodes = nodecount,#[1],
                                                         global_objects = "gbl",
                                                         slurm_options = list(time = job_time,
                                                                              mem = mem_gb))
                         
                         completed = F
                         
                         print("Waiting on cluster to finish jobs...")
                          
                         jobs_ntot = length(ml_queue$jobs)
                         
                         pb = pbapply::startpb(max = jobs_ntot)
                         while(!completed){
                           Sys.sleep(5)
                           running_jobs = list.files(paste0("_rslurm_", 
                                                            batch_job$jobname),
                                                     pattern = "rslurm")
                           pbapply::setpb(pb, value = length(running_jobs))
                           completed = slurm_job_complete(batch_job)
                         }

                         #rslurm::print_job_status(batch_job)
                         # my ver has a progress bar
                         print("Cluster batch job complete! Collecting results...")
                         ml_queue_res <- get_slurm_out_jw(batch_job,
                                                          outtype = "raw")
                         
                         names(ml_queue_res) <- sapply(ml_queue_res, function(x) x$params$ml_name)
                         rslurm::cleanup_files(batch_job) #cleanup files
                         
                         success = T
                       })
                       if(success){
                         print("reloading mSet...")
                         # load in mSet back
                         mSet <- qs::qread(fn)
                       }else{
                         stop("ml failed")
                       }
                     }else{
                       try({
                         parallel::clusterExport(ml_session_cl, c("ml_run",
                                                                  "gbl", 
                                                                  "mSet_loc"),
                                                 envir = environment())
                         
                         read_in = parallel::clusterEvalQ(ml_session_cl,{
                           small_mSet <- qs::qread(mSet_loc)
                         })  
                       })
                       
                       ml_queue_res <- pbapply::pblapply(ml_queue$jobs, 
                                         cl = if(length(ml_queue$jobs) > 1) ml_session_cl else 0, 
                                         function(settings, ml_cl){
                                           data = list()
                                           try({
                                             data = ml_run(settings = settings, 
                                                          mSet = small_mSet,
                                                          input = input,
                                                          cl = ml_cl,
                                                          tmpdir = dirname(tempfile()), 
                                                          use_slurm=F)  
                                           })
                                           data
                                         },
                                         ml_cl = if(length(ml_queue$jobs) > 1) 0 else ml_session_cl)
                     }
                   }
                   
                   print("Done!")
                   
                   shiny::showNotification("Gathering results...")
                   
                   ml_queue_res = ml_queue_res[unlist(sapply(ml_queue_res, function(l) length(l) > 0))]
                   
                   if(!("ml" %in% names(mSet$analSet))){
                     mSet$analSet$ml <- list()
                   }
                   
                   for(res in ml_queue_res){
                     mSet$analSet$ml[[res$params$ml_method]][[res$params$ml_name]] <- res
                     settings <- res$params
                   }
                   
                   if(length(mSet$analSet$ml[[res$params$ml_method]][[res$params$ml_name]]) > 0){
                     mSet$analSet$ml$last <- list(name = settings$ml_name,
                                                  method = settings$ml_method)  
                     shiny::showNotification("Done!")
                   }else{
                     shiny::showNotification("Failed...")
                   }
                   try({
                     parallel::stopCluster(ml_session_cl)
                   })
                   lcl$vectors$ml_train <- lcl$vectors$ml_train <<- NULL
                   print("ML done and saved.")
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
                   mSet <- Ttests.Anal.JW(mSet,
                                       nonpar = input$tt_nonpar,
                                       threshp = input$tt_p_thresh,
                                       paired = mSet$dataSet$ispaired,
                                       equal.var = input$tt_eqvar,
                                       multicorr_method = input$tt_multi_test)
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
                 
                 if(anal1 != anal2){
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
                 }else{
                   mSet$analSet$combi <- list(sig.mat = anal1_res_table[, c("rn",
                                                                            anal1_col,
                                                                            anal2_col), with=F], 
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