metshiCliffD <- function(mSet, input){
  matrices = split(mSet$dataSet$norm, 
                   mSet$dataSet$cls)
  scores = MetaboShiny::cliffDelta(x = matrices[[1]], 
                                   y = matrices[[2]])
  res = data.frame(cliffD = scores)
  rownames(res) <- colnames(mSet$dataSet$norm)
  names(scores) <- colnames(mSet$dataSet$norm)
  mSet$analSet$cliffd <- list(cliffd = scores, 
                              sig.mat = res)
  return(mSet)
}

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
    mSet_pca = mSet
    mSet_pca$dataSet$norm <- switch(input$pca_source,
                                "pre-batch correction" = mSet$dataSet$prebatch,
                                original = mSet$dataSet$proc)
    pcaRes <- MetaboAnalystR::PCA.Anal(mSet_pca)$analSet$pca # perform PCA analysis
    mSet$analSet$pca <- pcaRes  
  }else{
    mSet <- MetaboAnalystR::PCA.Anal(mSet) # perform PCA analysis
  }
  return(mSet)
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

metshiMultirank <- function(input, mSet, lcl, multirank_yes){
  analyses = unlist(multirank_yes$now$name)
  print(analyses)
  results = data.table::rbindlist(lapply(analyses, function(an){
    hits = MetaboShiny::getAllHits(mSet, an)
    hits$ranking = 1:nrow(hits)
    hits$analysis=an
    return(hits)
  }))
  results.avg <- results[,lapply(.SD, mean, na.rm=TRUE),by=`m.z`,.SDcols=c("ranking")]
  results.avg <- results.avg[order(ranking, decreasing = F)]
  results.avg$ranking <- rank(results.avg$ranking,ties.method = "average")
  results.avg$group = "mean"
  results$group = results$analysis
  results.merged = rbind(results.avg,
                         results[,c("m.z", "ranking", "group")])
  mSet$analSet$multirank = list(result_table = results.merged)
  mSet
}

doEnrich <- function(input, tempfile, ppm, lcl){
  
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
  try({
    shiny::showNotification("MetaboAnalyst R objects initialized ...")
  })
  
  hasP <- all(tbl$`p.value` == 0)
  hasT <- all(is.na(tbl$`t.score`))
  
  # -----------------
  
  MetaboAnalystR::SetPeakFormat("mpt")
  
  enr_mSet <- MetaboAnalystR::UpdateInstrumentParameters(enr_mSet,
                                                         ppm,
                                                         "mixed",
                                                         "yes",
                                                         0.02);
  
  enr_mSet <- MetaboAnalystR::Read.PeakListData(enr_mSet, tempfile);
  
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
  try({
    if (sel.add > 0) {
      shiny::showNotification(paste("A total of ", sel.add, 
                                    " adducts were successfully selected!", sep = ""))
    }else{
      shiny::showNotification("No adducts were selected!")
    }  
  })
  
  
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
    mw_modified <- list(mw_modified.neg, 
                        mw_modified.pos)
    
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
  
  cpds = mummichog.lib$cpd.lib$id
  mummi.adducts <- list(dpj_positive = mummi.adducts$pos,
                        positive = mummi.adducts$pos,
                        negative = mummi.adducts$neg)
  
  if(use.rules){
    kegg_user_db = MetaDBparse::showAllBase(lcl$paths$db_dir, "kegg")
    keep = data.table::as.data.table(kegg_user_db[kegg_user_db$identifier %in% cpds,])
    setkey(keep, "identifier")
    
    adduct_rows = pbapply::pblapply(cpds, 
                                   function(cpd){
      row = keep[identifier == cpd]
      if(nrow(row) > 0){
        ext_info = MetaDBparse::searchRev(row$structure,
                                          ext.dbname = "extended",
                                          outfolder = lcl$paths$db_dir)
        if(nrow(ext_info) > 0){
          use_info = unique(ext_info[ext_info$isoprevalence == 100, c("adduct", "fullmz")])
          use_info$identifier = cpd  
          reshape2::dcast(use_info, identifier ~ adduct, value.var ="fullmz")  
        }else{
          data.table::data.table()
        }
      }else{
        data.table::data.table()
      }
      })
    adducts_from_db = data.table::rbindlist(adduct_rows, fill = T)
    remaining_cpds = adducts_from_db$identifier
    adducts_from_db$identifier=NULL
    adducts_from_db[is.na(adducts_from_db)] <- -9999
    # split in positive/negative
    pos_cols = colnames(mummi.adducts$pos)
    neg_cols = colnames(mummi.adducts$neg)
    pos_adducts = adducts_from_db[,..pos_cols]
    neg_adducts = adducts_from_db[,..neg_cols]
    
    mummi.adducts <- list(dpj_positive = pos_adducts,
                          positive = pos_adducts,
                          negative = neg_adducts)
    
    mummichog.lib$cpd.lib$id <- remaining_cpds
    mummichog.lib$cpd.lib$name
    print("in progress...")

  }
  
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

#' Calculate Cliff's delta values (for continuous data)
#'
#' @rdname cliffDelta
#'
#' @description Calculate effect sizes using Cliff's delta values for pairwise continuous variable
#' comparisons (e.g. comparisons that would be done by a wilcox test or t-test)
#'
#' @param x A numeric vector for the 1st group of observations. Alternatively, a matrix, where
#' comparisons will be performed by row.
#' @param y Same as x but for the 2nd group of observations.
#'
#' @return A numeric vector of Cliff's delta values
#' @export
#'
cliffDelta <- function (x, ...) {
  UseMethod("cliffDelta", x)
}

#' @rdname cliffDelta
#' @method cliffDelta default
#' @export
cliffDelta.default <- function(x, y, use.r.implementation=USE_CLIFF_DELTA_R){
  
  if(!is.numeric(x) | !is.numeric(y)){ stop('x and y must be numeric matrices') }
  
  if(!use.r.implementation){
    cliffd(x,y)
  } else {
    signs <- sign(outer(x, y, FUN="-"))
    sum(signs, na.rm=T) / length(signs)
  }
}

#' @rdname cliffDelta
#' @method cliffDelta matrix
#' @export
cliffDelta.matrix  <- function(x, y, use.r.implementation=USE_CLIFF_DELTA_R){
  
  if(ncol(x)!=ncol(y)){ stop('x and y must have the sample number of columns') }
  x <- as.matrix(x); dimnames(x) <- NULL
  y <- as.matrix(y); dimnames(y) <- NULL
  
  if(!is.numeric(x) & !is.logical(x)){ stop('x must be a numeric or logical matrix') }
  if(!is.numeric(y) & !is.logical(y)){ stop('y must be a numeric or logical matrix') }
  
  if(!use.r.implementation){
    ## Cpp implementation
    ## Calculate cliff delta for every col between x and y
    
    unlist(pbapply::pblapply(1L:ncol(x), function(i){
      #i=1
      cliffd(x[,i],y[,i])
    }), use.names=F)
    
  } else {
    ## R implementation (~6-7x slower than C++ implementation)
    n_comparisons <- nrow(x)*nrow(y)
    
    ## Rows of x are compared with cols of y
    ## Therefore need to transpose
    y <- t(y)
    
    ## For every row of x, calculate how many x values are larger than values in the y matrix
    gt_sums <- rowSums(
      pbapply::pbapply(x,1,function(i){
        #i=x[1,]
        rowSums(i > y)
      })
    )
    
    ## ... same for less than
    lt_sums <- rowSums(
      pbapply::pbapply(x,1,function(i){
        rowSums(i < y)
      })
    )
    
    sign_sums <- gt_sums - lt_sums
    sign_sums / n_comparisons
  }
}

#' @rdname cliffDelta
#' @method cliffDelta data.frame
#' @export
cliffDelta.data.frame <- cliffDelta.matrix
