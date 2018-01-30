# remove duplicates, peak groups with mz within a few ppm
# ppmdef should be 2 times bigger than during identification!!!
remove.dupl <- function(peaklist, ppmdef=4, tolint=0.01) {

#   peaklist = outpgrlist
#   ppmdef = 2
#   tolint = 0.01

  int.cols = 7:ncol(peaklist)
  
  p <- 1
  while (p < nrow(peaklist)) {
    mzref <- peaklist[p, "mzmed.pgrp"]
    # print(mzref)
    dist.ppm <- ppmdef * mzref / 1000000 
    sel <- peaklist[ , "mzmed.pgrp"] > (mzref - dist.ppm) & peaklist[ , "mzmed.pgrp"] < (mzref + dist.ppm)
    subset <- peaklist[sel, ]
    if (nrow(subset) > 1) {
      avi <- rep(1, max(int.cols))
      for (c in int.cols) { avi[c] <- max(subset[ ,c])/mean(subset[ ,c]) }
      
      # remove NaN
      avi[which(is.nan(avi))] = 1
      
      if (mean(avi) > (1-tolint) & mean(avi < (1+tolint))) { 
        peaklist <- peaklist[-which(sel), ]
        newline <- subset[1, ]
        newline[ , "mzmed.pgrp"] <- mean(subset[ , "mzmed.pgrp"])
        newline[ , "mzmin.pgrp"] <- min(subset[ , "mzmin.pgrp"])
        newline[ , "mzmax.pgrp"] <- max(subset[ , "mzmax.pgrp"])
        newline[ , int.cols] <- apply(subset[ , int.cols], 2, max)
        # newline[ , "avg.int"] <- mean(as.numeric(newline[ , int.cols]))
        peaklist <- rbind(peaklist, newline)
      }
    }
    p <- p + 1
  }
  peaklist
}
