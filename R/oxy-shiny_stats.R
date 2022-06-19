metshiCorr <- function(mSet, input){
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
  return(mSet)
}

metshiDiffCorr <- function(mSet, input){
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
  return(mSet)
}

metshiPCA <- function(mSet, input){
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
}

metshiICA <- function(mSet, input){
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
  return(mSet)
}

metshiUMAP <- function(mSet, input){
  inTbl = switch(input$umap_source, 
                 original = mSet$dataSet$proc,
                 "pre-batch correction" = mSet$dataSet$prebatch,
                 normalized = mSet$dataSet$norm)
  umapRes = umap::umap(d = inTbl,
                       method = "naive",
                       n_components = input$umap_ncomp,
                       n_neighbors = input$umap_neighbors)
  mSet$analSet$umap <- umapRes
  return(mSet)
}

metshiGetEnrichInputTable <- function(mSet, input){
  flattened <- list(getAllHits(mSet,
                               input$mummi_anal))
  
  if(input$mummi_enr_method){
    flattened[[1]]$value <- 1
    flattened[[1]]$significant <- F
    flattened[[1]] <- flattened[[1]][order(abs(flattened[[1]]$statistic),
                                           decreasing = T),]
    flattened[[1]]$value[1:input$mummi_topn] <- 0
    flattened[[1]]$significant[1:input$mummi_topn] <- T
    if(input$enrich_use_corr){
      rcorrMat = Hmisc::rcorr(x = as.matrix(mSet$dataSet$norm),
                              type = input$enrich_corr_method)
      of_interest = flattened[[1]][flattened[[1]]$value == 0,][["m.z"]]
      corr_mz = unique(unlist(sapply(of_interest, function(mz){
        names(which(rcorrMat$r[,mz] >= input$enrich_min_corr))
      })))
      flattened[[1]][flattened[[1]]$m.z %in% corr_mz,]$value <- 0
    }
  }
  
  hasP = T#grepl("tt|aov|asca|combi|venn",input$mummi_anal)

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
  # ----------
  return(tbl)
}

doEnrich <- function(input, tempfile, ppm){
  
  tbl <- data.table::fread(tempfile)
  dataSet <- list()
  dataSet$type <- "mass_all"
  dataSet$design.type <- "regular"
  dataSet$cls.type <- "disc"
  dataSet$format <- "rowu"
  dataSet$paired <- FALSE
  analSet <- list()
  analSet$type <- "mummichog"
  enr_mSet <- list()
  enr_mSet$dataSet <- dataSet
  enr_mSet$analSet <- analSet
  enr_mSet$imgSet <- list()
  enr_mSet$msgSet <- list()
  enr_mSet$msgSet$msg.vec <- vector(mode = "character")
  enr_mSet$cmdSet <- vector(mode = "character")
  MetaboAnalystR:::.init.global.vars("mummichog")
  shiny::showNotification("MetaboAnalyst R objects initialized ...")
  
  hasP <- all(tbl$`p.value` == 0)
  hasT <- all(is.na(tbl$`t.score`))
  
  # -----------------
  
  MetaboAnalystR::SetPeakFormat("mpt")
  
  enr_mSet <- MetaboAnalystR::UpdateInstrumentParameters(enr_mSet,
                                                         ppm,
                                                         "mixed",
                                                         "yes",
                                                         0.02);
  
  enr_mSet <- MetaboAnalystR::Read.PeakListData(enr_mSet, tmpfile);
  
  enr_mSet <- MetaboAnalystR::SanityCheckMummichogData(enr_mSet)
  enr_mSet <- MetaboAnalystR::Setup.AdductData(enr_mSet, input$mummi_adducts);
  
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
  
  # === PerformAdductMapping ===
  myAdds <- enr_mSet$dataSet$adduct.list
  hit.inx <- match(tolower(myAdds), 
                   tolower(add_db_cust$Ion_Name))
  match.values <- add_db_cust[hit.inx, ]
  sel.add <- nrow(match.values)
  if (sel.add > 0) {
    shiny::showNotification(paste("A total of ", sel.add, 
                                  " adducts were successfully selected!", sep = ""))
  }else{
    shiny::showNotification("No adducts were selected!")
  }
  
  enr_mSet$add.map <- match.values
  
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
    return(mw_modified)
  }
  
  assignInNamespace("new_adduct_mzlist", new_adduct_mzlist, ns="MetaboAnalystR", 
                    envir=as.environment("package:MetaboAnalystR"))
  
  enr_mSet$dataSet$N <- 20
  
  # ==== BUILD CUSTOM DB HERE ====
  
  map_id = input$mummi_org

  lib_name <- paste0(map_id, "_kegg")
  file_name <- paste0(lib_name, ".qs")
  
  if(!file.exists(file_name)){
    mummichog.lib <- build.enrich.KEGG(map_id)
    qs::qsave(mummichog.lib, file = file_name)
  }  
  
  
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
  return(enr_mSet)
}

metshiNetwork <- function(mSet, input){
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
  return(mSet)
}