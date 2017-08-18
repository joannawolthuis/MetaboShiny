#' @export
import.outlist <- function(mode,
                           resDir){
  if(mode == "negative" ){
    load(file.path(resDir, "outlist2identify_negative.RData")) # From HPC cluster
  } else if(mode == "positive") {
    load(file.path(resDir, "outlist2identify_positive.RData")) # From HPC cluster
  }
  # --- return ---
  outlist
}

#' @export
filter.outlist <- function(outlist, sample.characteristic = "FR1", threshold = .8){ # must be in 80% of samples
  # datatable-ify
  outlist <- as.data.table(outlist)
  # get sample cols
  samp.cols <- grep(sample.characteristic, colnames(outlist), value = TRUE)
  # get min presence needed
  min.nsamp <- threshold * length(samp.cols)
  # filter
  keep <- outlist[nrsamples > min.nsamp]
  # --- return ---
  keep
}

#' @export
remove.outliers <- function(outlist.pos, outlist.neg, prefix){
  # conv to data table
  pos <- as.data.table(outlist.pos)
  neg <- as.data.table(outlist.neg)
  # join both temporarily
  outlist <- rbind(pos, neg)
  # find outliers
  sample.cols <- grep(pattern = prefix, x = names(outlist))
  sample.total.intensities <- colSums(outlist[sample.cols])
  outlier.results <- boxplot.stats(sample.total.intensities)$out
  outlier.samples <- grep(pattern = prefix, x = names(outlier.results))
  # --- return ---
  results <- list(pos = pos[,-outlier.samples, with=F],
                  neg = neg[,-outlier.samples, with=F])
  results
}

#' @export
rename.old.to.new <- function(outlist, 
                              old.prefix="BSP", 
                              replicates=3, 
                              new.prefix){
  # get sample names
  samples <- grep("BSP", colnames(outlist), value = T)
  # sort properly
  before.names <- mixedsort(samples, decreasing=TRUE)
  # create a translation key
  unique.sample.amount <- length(before.names) / replicates
  after.names <- paste(new.prefix, rep(1:unique.sample.amount, each=replicates), c(1:replicates), sep = "-")
  # replace names
  outlist.colnames <- mgsub(before.names, after.names, names(outlist)) # YAY IT WORKS
  names(outlist) <- outlist.colnames
  # --- return ---
  outlist
}