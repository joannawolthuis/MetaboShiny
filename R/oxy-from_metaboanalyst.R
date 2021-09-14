#' @export
Ttests.Anal.JW <- function (mSetObj = NA, nonpar = F, threshp = 0.05, paired = FALSE, 
          equal.var = TRUE, all_results = FALSE, multicorr_method="fdr"){
  mSetObj <- MetaboAnalystR:::.get.mSet(mSetObj)
  res <- MetaboAnalystR:::GetTtestRes(mSetObj, paired, equal.var, nonpar)
  t.stat <- res[, 1]
  p.value <- res[, 2]
  names(t.stat) <- names(p.value) <- colnames(mSetObj$dataSet$norm)
  p.log <- -log10(p.value)
  fdr.p <- p.adjust(p.value, multicorr_method)
  if (all_results == TRUE) {
    all.mat <- data.frame(signif(t.stat, 5), signif(p.value, 
                                                    5), signif(p.log, 5), signif(fdr.p, 5))
    if (nonpar) {
      tt.nm = "Wilcoxon Rank Test"
      file.nm <- "wilcox_rank_all.csv"
      colnames(all.mat) <- c("V", "p.value", "-log10(p)", 
                             multicorr_method)
    }
    else {
      tt.nm = "T-Tests"
      file.nm <- "t_test_all.csv"
      colnames(all.mat) <- c("t.stat", "p.value", "-log10(p)", 
                             multicorr_method)
    }
    MetaboAnalystR:::fast.write.csv(all.mat, file = file.nm)
  }
  inx.imp <- fdr.p <= threshp
  sig.num <- sum(inx.imp)
  if (is.na(sig.num)) {
    MetaboAnalystR:::AddMsg(paste("No significant features were found."))
    return(0)
  }
  else {
    MetaboAnalystR:::AddMsg(paste("A total of", sig.num, "significant features were found."))
  }
  if (sig.num > 0) {
    sig.t <- t.stat[inx.imp]
    sig.p <- p.value[inx.imp]
    lod <- -log10(sig.p)
    sig.q <- fdr.p[inx.imp]
    sig.mat <- cbind(sig.t, sig.p, lod, sig.q)
    colnames(sig.mat) <- c("t.stat", "p.value", "-log10(p)", 
                           multicorr_method)
    ord.inx <- order(sig.p)
    sig.mat <- sig.mat[ord.inx, , drop = F]
    sig.mat <- signif(sig.mat, 5)
    if (nonpar) {
      tt.nm = "Wilcoxon Rank Test"
      file.nm <- "wilcox_rank.csv"
      colnames(sig.mat) <- c("V", "p.value", "-log10(p)", 
                             multicorr_method)
    }
    else {
      tt.nm = "T-Tests"
      file.nm <- "t_test.csv"
      colnames(sig.mat) <- c("t.stat", "p.value", "-log10(p)", 
                             multicorr_method)
    }
    MetaboAnalystR:::fast.write.csv(sig.mat, file = file.nm)
    tt <- list(tt.nm = tt.nm, sig.nm = file.nm, sig.num = sig.num, 
               paired = paired, raw.thresh = threshp, t.score = sort(t.stat), 
               p.value = sort(p.value), p.log = p.log, thresh = -log10(threshp), 
               inx.imp = inx.imp, sig.mat = sig.mat)
  }
  else {
    tt <- list(sig.num = sig.num, paired = paired, raw.thresh = threshp, 
               t.score = sort(t.stat), p.value = sort(p.value), 
               p.log = p.log, thresh = -log10(threshp), inx.imp = inx.imp)
  }
  mSetObj$analSet$tt <- tt

  return(MetaboAnalystR:::.set.mSet(mSetObj))
}
