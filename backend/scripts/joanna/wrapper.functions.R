#' @export
get.csv <- function(patdb, 
                    time.series = T, 
                    #exp.condition = "diet",
                    max.vals = -1,
                    group_adducts = T,
                    which_dbs = file.path(options$db_dir, "kegg.full.db"),
                    which_adducts = c("M+H", "M-H", "M"),
                    group_by = "mz"
                    #,var_table = "setup",
                    #batches = NULL
                    ){
  library(data.table)
  # --- announce some stuff ---
  adducts <- paste(which_adducts, collapse = ", ")
  groupfac <- group_by
  max.cols <- if(max.vals == -1) "unlimited" else max.vals + 3
  cat(gsubfn::fn$paste("Creating csv for metabolomics analysis with max $max.cols columns.
                - Grouping by $groupfac
                - Using adducts $adducts"))
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)
  #batches <- batches[-which(batches == "run")]
  #if(exp.condition == "run") var_table <- "run"
  # ----------------
  if(group_adducts){
    z = get_all_matches(#exp.condition, 
                        pat.conn = conn,
                        which_dbs,
                        which_adducts,
                        group_by
                        #,var_table,
                        #batches
                        )

  }else{
    query <- strwrap(gsubfn::fn$paste("select distinct d.*, s.*, b.*,
                              i.mzmed as identifier,
                              i.intensity
                              from mzintensities i
                              join individual_data d
                              on i.filename = d.card_id
                              join setup s on d.[Group] = s.[Group]
                              join batchinfo b on b.sample = d.card_id
                              group by d.card_id, 
                              d.sampling_date, 
                              i.mzmed"),
                     width=10000,
                     simplify=TRUE)
    z = RSQLite::dbGetQuery(conn, query)
  }
  RSQLite::dbDisconnect(conn)
  z.dt <- as.data.table(z)
  nvars = ncol(z.dt) - 3
  cast.dt <- dcast.data.table(z.dt, 
                              formula = ... ~ identifier,
                              fun.aggregate = sum, 
                              value.var = "intensity") # what to do w/ duplicates? 
  # --- cast to right format ---
  small.set <- cast.dt[,1:nvars,]
  # --- name for metaboanalyst ---
  colnames(small.set)[which(colnames(small.set) == "label")] <- "#"
  colnames(small.set)[which(colnames(small.set) == "card_id")] <- "Sample"
  colnames(small.set)[which(colnames(small.set) == "sampling_date")] <- "Time"
  names(small.set) <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(colnames(small.set)), perl=TRUE)
  small.set[,which(unlist(lapply(small.set, function(x)!all(is.na(x))))),with=F]
  small.set <- cbind(small.set, "$" = c(0))
  
  # --- rejoin w/ rest ---
  print(head(small.set))
  
  small.set <- cbind(small.set, cast.dt[,-c(1:nvars), with=FALSE])
  
  # --- make time series if necessary (this factorizes sampling date) ---
  if(time.series){
      small.set$Time <- as.numeric(as.factor(as.Date(small.set$Time)))
      small.set$Sample <- paste(small.set$Sample, as.character(small.set$Time), sep="_T")
  }
  # --- measure file size ---
  size <- object.size(small.set)
  print(small.set[1:5,1:20])
  cat(paste("Resulting file will be approximately "))
  print(size, units = "MB")
  # ------- return -------
  return(small.set)
}

Read.TextData.J <- function (mSetObj = NA, filePath, format = "rowu", lbl.type = "disc") 
{
  mSetObj <- MetaboAnalystR:::.get.mSet(mSetObj)
  mSetObj$dataSet$cls.type <- lbl.type
  mSetObj$dataSet$format <- format
  dat <- MetaboAnalystR:::.readDataTable(filePath)
  print(dat[1:10,1:10])
  if (class(dat) == "try-error") {
    print("here1")
    AddErrMsg(mSetObj, "Data format error. Failed to read in the data!")
    AddErrMsg(mSetObj, "Please check the followings: ")
    AddErrMsg(mSetObj, "Either sample or feature names must in UTF-8 encoding; Latin, Greek letters are not allowed.")
    AddErrMsg(mSetObj, "We recommend using a combination of English letters, underscore, and numbers for naming purpose")
    AddErrMsg(mSetObj, "Make sure sample names and feature (peak, compound) names are unique;")
    AddErrMsg(mSetObj, "Missing values should be blank or NA without quote.")
    return(0)
  }
  if (ncol(dat) == 1) {
    print("here2")
    AddErrMsg(mSetObj, "Error: Make sure the data table is saved as comma separated values (.csv) format!")
    AddErrMsg(mSetObj, "Please also check the followings: ")
    AddErrMsg(mSetObj, "Either sample or feature names must in UTF-8 encoding; Latin, Greek letters are not allowed.")
    AddErrMsg(mSetObj, "We recommend to use a combination of English letters, underscore, and numbers for naming purpose.")
    AddErrMsg(mSetObj, "Make sure sample names and feature (peak, compound) names are unique.")
    AddErrMsg(mSetObj, "Missing values should be blank or NA without quote.")
    return(0)
  }
  msg <- NULL
  if (substring(format, 4, 5) == "ts") {
    if (substring(format, 1, 3) == "row") {
      msg <- c(msg, "Samples are in rows and features in columns")
      smpl.nms <- dat[, 1]
      all.nms <- colnames(dat)
      facA.lbl <- all.nms[2]
      cls.lbl <- facA <- dat[, 2]
      facB.lbl <- all.nms[3]
      facB <- dat[, 3]
      conc <- dat[, -c(1:3)]
      var.nms <- colnames(conc)
    }
    else {
      msg <- c(msg, "Samples are in columns and features in rows.")
      all.nms <- dat[, 1]
      facA.lbl <- all.nms[1]
      cls.lbl <- facA <- dat[1, -1]
      facB.lbl <- all.nms[2]
      facB <- dat[2, -1]
      var.nms <- dat[-c(1:2), 1]
      conc <- t(dat[-c(1:2), -1])
      smpl.nms <- rownames(conc)
    }
    facA <- as.factor(as.character(facA))
    facB <- as.factor(as.character(facB))
    if (mSetObj$dataSet$design.type == "time" | mSetObj$dataSet$design.type == 
        "time0") {
      if (!(tolower(facA.lbl) == "time" | tolower(facB.lbl) == 
            "time")) {
        AddErrMsg(mSetObj, "No time points found in your data")
        AddErrMsg(mSetObj, "The time points group must be labeled as <b>Time</b>")
        print("here3")
        return(0)
      }
    }
  }
  else {
    print(dat[1:10,1:10])
    if (substring(format, 1, 3) == "row") {
      msg <- c(msg, "Samples are in rows and features in columns")
      smpl.nms <- dat[, 1]
      dat[, 1] <- NULL
      if (lbl.type == "qc") {
        rownames(dat) <- smpl.nms
        mSetObj$dataSet$orig <- dat
        mSetObj$dataSet$cmpd <- colnames(dat)
        return(1)
      }
      cls.lbl <- dat[, 1]
      conc <- dat[, -1]
      var.nms <- colnames(conc)
    }
    else {
      msg <- c(msg, "Samples are in columns and features in rows.")
      var.nms <- dat[-1, 1]
      dat[, 1] <- NULL
      smpl.nms <- colnames(dat)
      cls.lbl <- dat[1, ]
      conc <- t(dat[-1, ])
    }
  }
  print(cls.lbl)
  dat <- NULL
  msg <- c(msg, "The uploaded file is in comma separated values (.csv) format.")
  empty.inx <- is.na(smpl.nms) | smpl.nms == ""
  if (sum(empty.inx) > 0) {
    msg <- c(msg, paste(sum(empty.inx), 
                        "empty rows were detected and excluded from your data."))
    smpl.nms <- smpl.nms[!empty.inx]
    cls.lbl <- cls.lbl[!empty.inx]
    conc <- conc[!empty.inx, ]
  }
  empty.inx <- is.na(cls.lbl) | cls.lbl == ""
  if (sum(empty.inx) > 0) {
    if (mSetObj$analSet$type != "roc") {
      msg <- c(msg, paste(sum(empty.inx), 
                          "were detected and excluded from your data."))
      smpl.nms <- smpl.nms[!empty.inx]
      cls.lbl <- cls.lbl[!empty.inx]
      conc <- conc[!empty.inx, ]
    }
    else {
      cls.lbl[is.na(cls.lbl)] <- ""
      msg <- c(msg, paste("<font color=\"orange\">", sum(empty.inx), 
                          "new samples</font> were detected from your data."))
    }
  }
  if (mSetObj$analSet$type == "roc") {
    print("here4")
    if (length(unique(cls.lbl[!empty.inx])) > 2) {
      AddErrMsg(mSetObj, "ROC analysis is only defined for two-group comparisions!")
      return(0)
    }
  }
  empty.inx <- is.na(smpl.nms) | smpl.nms == ""
  if (sum(empty.inx) > 0) {
    msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx), 
                        "empty samples</font> were detected and excluded from your data."))
    smpl.nms <- smpl.nms[!empty.inx]
    cls.lbl <- cls.lbl[!empty.inx]
    conc <- conc[!empty.inx, ]
  }
  if (length(unique(smpl.nms)) != length(smpl.nms)) {
    dup.nm <- paste(smpl.nms[duplicated(smpl.nms)], collapse = " ")
    AddErrMsg(mSetObj, "Duplicate sample names are not allowed!")
    AddErrMsg(mSetObj, dup.nm)
    print("here5")
    return(0)
  }
  empty.inx <- is.na(var.nms) | var.nms == ""
  if (sum(empty.inx) > 0) {
    msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx), 
                        "empty features</font> were detected and excluded from your data."))
    var.nms <- var.nms[!empty.inx]
    conc <- conc[, !empty.inx]
  }
  if (length(unique(var.nms)) != length(var.nms)) {
    dup.nm <- paste(var.nms[duplicated(var.nms)], collapse = " ")
    AddErrMsg(mSetObj, "Duplicate feature names are not allowed!")
    AddErrMsg(mSetObj, dup.nm)
    print("here6")
    return(0)
  }
  if (sum(is.na(iconv(smpl.nms))) > 0) {
    na.inx <- is.na(iconv(smpl.nms))
    nms <- paste(smpl.nms[na.inx], collapse = "; ")
    AddErrMsg(mSetObj, paste("No special letters (i.e. Latin, Greek) are allowed in sample names!", 
                             nms, collapse = " "))
    print("here7")
    return(0)
  }
  if (sum(is.na(iconv(var.nms))) > 0) {
    na.inx <- is.na(iconv(var.nms))
    nms <- paste(var.nms[na.inx], collapse = "; ")
    AddErrMsg(mSetObj, paste("No special letters (i.e. Latin, Greek) are allowed in feature names!", 
                             nms, collapse = " "))
    print("here8")
    return(0)
  }
  smpl.nms <- gsub("[^[:alnum:]./_-]", "", smpl.nms)
  var.nms <- gsub("[^[:alnum:][:space:],'./_-]", "", var.nms)
  cls.lbl <- MetaboAnalystR:::ClearStrings(as.vector(cls.lbl))
  rownames(conc) <- smpl.nms
  colnames(conc) <- var.nms
  if (mSetObj$dataSet$paired) {
    mSetObj$dataSet$orig.cls <- mSetObj$dataSet$pairs <- cls.lbl
  }
  else {
    if (lbl.type == "disc") {
      if (min(table(cls.lbl)) < 2) {
        print(table(cls.lbl))
        print(paste("A total of", length(levels(as.factor(cls.lbl))), 
                                 "groups found with", length(smpl.nms), "samples."))
        print(paste("A total of", length(levels(as.factor(cls.lbl))), 
              "groups found with", length(smpl.nms), "samples."))
        AddErrMsg(mSetObj, "At least three replicates are required in each group!")
        AddErrMsg(mSetObj, "Or maybe you forgot to specify the data format?")
        print(cls.lbl)
        print("here9")
        return(0)
      }
      mSetObj$dataSet$orig.cls <- mSetObj$dataSet$cls <- as.factor(as.character(cls.lbl))
      if (substring(format, 4, 5) == "ts") {
        mSetObj$dataSet$facA <- as.factor(as.character(facA))
        mSetObj$dataSet$orig.facA <- as.factor(as.character(facA))
        mSetObj$dataSet$facA.lbl <- facA.lbl
        mSetObj$dataSet$facB <- as.factor(as.character(facB))
        mSetObj$dataSet$orig.facB <- as.factor(as.character(facB))
        mSetObj$dataSet$facB.lbl <- facB.lbl
      }
    }
    else {
      mSetObj$dataSet$orig.cls <- mSetObj$dataSet$cls <- as.numeric(cls.lbl)
    }
  }
  if (mSetObj$dataSet$type == "conc") {
    mSetObj$dataSet$cmpd <- var.nms
  }
  mSetObj$dataSet$orig <- conc
  mSetObj$msgSet$read.msg <- c(msg, paste("The uploaded data file contains ", 
                                          nrow(conc), " (samples) by ", ncol(conc), " (", tolower(GetVariableLabel(mSetObj)), 
                                          ") data matrix.", sep = ""))
  print(mSetObj$msgSet$read.msg)
  return(MetaboAnalystR:::.set.mSet(mSetObj))
}

Ttests.Anal.J <- function (mSetObj = NA, nonpar = F, threshp = 0.05, paired = FALSE, 
                           equal.var = TRUE, multicorr) 
{
  mSetObj <- MetaboAnalystR:::.get.mSet(mSetObj)
  res <- MetaboAnalystR:::GetTtestRes(mSetObj, paired, equal.var, nonpar)
  t.stat <- res[, 1]
  p.value <- res[, 2]
  print(any(p.value < threshp))
  names(t.stat) <- names(p.value) <- colnames(mSetObj$dataSet$norm)
  p.log <- -log10(p.value)
  fdr.p <- p.adjust(p.value, multicorr)
  inx.imp <- which(fdr.p <= threshp)
  if(length(inx.imp)==0) inx.imp = 0
  #print(inx.imp)
  msg <- NULL
  mSetObj$msgSet$current.msg <- paste(c(msg, "A total of", 
                                        sum(inx.imp), "significant features were found."), collapse = " ")
  sig.num <- sum(inx.imp)
  print(sig.num)
  if (sig.num > 0) {
    sig.t <- t.stat[inx.imp]
    sig.p <- p.value[inx.imp]
    lod <- -log10(sig.p)
    sig.q <- fdr.p[inx.imp]
    sig.mat <- cbind(sig.t, sig.p, lod, sig.q)
    colnames(sig.mat) <- c("t.stat", "p.value", "-log10(p)", 
                           "FDR")
    ord.inx <- order(sig.p)
    sig.mat <- sig.mat[ord.inx, , drop = F]
    sig.mat <- signif(sig.mat, 5)
    if (nonpar) {
      tt.nm = "Wilcoxon Rank Test"
      file.nm <- "wilcox_rank.csv"
      colnames(sig.mat) <- c("V", "p.value", "-log10(p)", 
                             "FDR")
    }
    else {
      tt.nm = "T-Tests"
      file.nm <- "t_test.csv"
      colnames(sig.mat) <- c("t.stat", "p.value", "-log10(p)", 
                             "FDR")
    }
    write.csv(sig.mat, file = file.nm)
    tt <- list(tt.nm = tt.nm, sig.nm = file.nm, sig.num = sig.num, 
               paired = paired, raw.thresh = threshp, p.value = sort(p.value), 
               p.log = p.log, thresh = -log10(threshp), inx.imp = inx.imp, 
               sig.mat = sig.mat)
  }
  else {
    tt <- list(sig.num = sig.num, paired = paired, raw.thresh = threshp, 
               p.value = sort(p.value), p.log = p.log, thresh = -log10(threshp), 
               inx.imp = inx.imp)
  }
  mSetObj$analSet$tt <- tt
  return(MetaboAnalystR:::.set.mSet(mSetObj))
}

doBC_J <- function (Xvec, ref.idx, batch.idx, seq.idx, result = c("correctedX", 
                                                                  "corrections"), method = c("lm", "rlm", "tobit"), correctionFormula = formula("X ~ S * B"), 
                    minBsamp = ifelse(is.null(seq.idx), 2, 4), imputeVal = NULL, 
                    ...) 
{
  result <- match.arg(result)
  method <- match.arg(method)
  if (is.null(imputeVal) & method == "tobit") 
    stop("Tobit regression requires a value for 'imputeVal'")
  batch.idx <- factor(batch.idx)
  nbatches <- nlevels(batch.idx)
  if (is.factor(seq.idx)) 
    seq.idx <- as.numeric(levels(seq.idx))[seq.idx]
  if (is.null(seq.idx) & method != "lm") {
    warning("Using method = 'lm' since seq.idx equals NULL")
  }
  if (is.logical(ref.idx)) 
    ref.idx <- which(ref.idx)
  Xref <- Xvec[ref.idx]
  Bref <- batch.idx[ref.idx]
  Sref <- seq.idx[ref.idx]
  glMean <- mean(Xref, na.rm = TRUE)
  nNonNA <- tapply(Xref, Bref, function(x) sum(!is.na(x)))
  tooFew <- names(nNonNA)[nNonNA < minBsamp]
  if (length(tooFew) > length(levels(Bref)) - 2) 
    return(rep(NA, length(Xvec)))
  if (is.null(seq.idx)) {
    if (!is.null(imputeVal)) 
      Xref[is.na(Xref)] <- imputeVal
    Bmod <- lm(Xref ~ Bref - 1)
    Bcorrections <- (glMean - coef(Bmod))[batch.idx]
    switch(result, correctedX = Xvec + Bcorrections, Bcorrections)
  }
  else {
    if (!is.null(imputeVal)) 
      Xref[is.na(Xref)] <- imputeVal
    fitdf <- data.frame(S = Sref, B = Bref, X = Xref)
    if (length(tooFew) > 0) {
      fitdf <- fitdf[!(fitdf$B %in% tooFew), ]
      fitdf$B <- factor(fitdf$B)
    }
    print("here??")
    Bmods2 <- switch(method, lm = lm(correctionFormula, 
                                     data = fitdf), rlm = MASS::rlm(correctionFormula, 
                                                                    data = fitdf, maxit = 50), tobit = crch::crch(correctionFormula, 
                                                                                                                  data = fitdf, left = imputeVal))
    predictdf <- data.frame(S = seq.idx, B = batch.idx)
    predictdf$B[predictdf$B %in% tooFew] <- NA
    predictions <- rep(NA, length(Xvec))
    predictions[!(predictdf$B %in% tooFew)] <- predict(Bmods2, 
                                                       newdata = predictdf)
    switch(result, correctedX = Xvec + glMean - predictions, 
           glMean - predictions)
  }
}

