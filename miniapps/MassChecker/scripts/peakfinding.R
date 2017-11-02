library(xlsx)
library(xcms)
library(data.table)
library(snow)
library(BiocParallel)
library(parallel)
library(pbapply)
library(gsubfn)

groupval.joanna <- function(peakmat,
                            groupmat,
                            groupidx,
                            sampnames,
                            method = c("medret", "maxint"), 
                            value = "index", intensity = "into"){
  if (length(groupidx) < 
      1) {
    stop("xcmsSet is not been group()ed.")
  }
  groupindex <- groupidx
  sampnum <- seq(length = length(sampnames))
  retcol <- match("rt", colnames(peakmat))
  intcol <- match(intensity, colnames(peakmat))
  sampcol <- match("sample", colnames(peakmat))
  values <- matrix(nrow = length(groupindex), ncol = length(sampnum))
  if (method == "medret") {
    for (i in seq(along = groupindex)) {
      gidx <- groupindex[[i]][order(abs(peakmat[groupindex[[i]], 
                                                retcol] - median(peakmat[groupindex[[i]], 
                                                                         retcol])))]
      values[i, ] <- gidx[match(sampnum, peakmat[gidx, 
                                                 sampcol])]
    }
  }
  else {
    for (i in seq(along = groupindex)) {
      gidx <- groupindex[[i]][order(peakmat[groupindex[[i]], 
                                            intcol], decreasing = TRUE)]
      values[i, ] <- gidx[match(sampnum, peakmat[gidx, 
                                                 sampcol])]
    }
  }
  if (value != "index") {
    values <- peakmat[values, value]
    dim(values) <- c(length(groupindex), length(sampnum))
  }
  colnames(values) <- sampnames
  values_full <- cbind(groupmat[, c("mzmin",
                                    "mzmed",
                                    "mzmax", "npeaks")], values)
  values_full
}

# ====== SCALE EXPERIMENTATION ========

scales = list(
  a = seq(1,64,2),
  b = seq(1,22,3),
  c = seq(1,5,9),
  d = seq(1,7)
)

## ======= SOME EXPERIMENTATION =======
# scales = list(
#   a = seq(1,64,2),
#   b = seq(1,22,3),
#   c = seq(1,5,9),
#   d = seq(1,7)
# )
# a, snthresh = 10, nearbyPeak=T: 12109
# b, snthresh = 10, nearbyPeak=T: 10877
# c, snthresh = 10, nearbyPeak=T: phail
# d, snthresh = 10, nearbyPeak=T: 13637

# a, snthresh = 10, nearbyPeak=F: 10251
# b, snthresh = 10, nearbyPeak=F: 9749
# c, snthresh = 10, nearbyPeak=F: phail
# d, snthresh = 10, nearbyPeak=F: 13450

# a, snthresh = 3, nearbyPeak=T: 22890 <<<<<<< for now go w/ this
# b, snthresh = 3, nearbyPeak=T: 20609
# c, snthresh = 3, nearbyPeak=T: phail
# d, snthresh = 3, nearbyPeak=T: 26071

# a, snthresh = 3, nearbyPeak=F: 13496
# b, snthresh = 3, nearbyPeak=F: 13965
# c, snthresh = 3, nearbyPeak=F: phail
# d, snthresh = 3, nearbyPeak=F: 23949

testfile <- allfiles[1]

for(scale in scales){
  print(scale)
  attempt <- 0
  xset = NULL
  startscans = c(10, 15, 20)
  while( is.null(xset) && attempt < 3 ) {
    attempt <- attempt + 1
    try(
      xset <- xcmsSet(method="MSW",
                    files=testfile,
                    scales=scale,
                    snthresh=3,
                    nearbyPeak=F)
    )
  }
  if(!is.null(xset)) print(nrow(xset@peaks))
}

# =====================================

# hpc : 
nslots <- Sys.getenv( "NSLOTS" )
nslots = parallel::detectCores()

print(fn$paste("Doing MSW peak calling with $nslots slots!"))

allfiles = list.files("/Users/jwolthuis/PROCESSING_HPC", pattern = "\\.mzXML", full.names = T)
#bpparam = SnowParam(workers = nslots, type = "SOCK")

scales = list(
  a = seq(1,64,2),
  b = seq(1,22,3),
  c = seq(1,5,9)
)

scale = seq(1,64,2)

cl = parallel::makeCluster(4, "FORK", outfile="")
parallel::stopCluster(cl)

xcmsItems_pos <- pblapply(seq_along(allfiles), cl=cl, FUN=function(i){
  # --------------------------------------
  f <- allfiles[i]
  attempt <- 0
  # --------------------------------------
  xset = NULL
  startscans = c(10, 15, 20)
  # --------------------------------------
  while( is.null(xset) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
      xset <- xcmsSet(method="MSW",
                      files=f,
                      scales=scale,
                      snthresh=3,
                      nearbyPeak=T)
    )
  }
  if(!is.null(xset)) print(nrow(xset@peaks))
  xset
})


xcmsItems_neg <- pblapply(seq_along(allfiles), cl=cl, FUN=function(i){
  # --------------------------------------
  print(i)
  f <- filenames[i]
  # --------------------------------------
  xset = NULL
  attempt <- 0
  endscans = c(300, 290, 280)
  # --------------------------------------
  while( is.null(xset) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
      xset <- xcmsSet(f,
                      method="MSW",
                      scales = scale,
                      scanrange = c(151:endscans[attempt]),
                      snthresh=3
      )  
    )
  }
  # --- return ---
  if(!is.null(xset)) print(nrow(xset@peaks))
  xset
})

# ----------------------------

xset_pos <- xcmsItems_pos[[1]]
xset_neg <- xcmsItems_neg[[1]]

for(i in 2:length(xcmsItems_pos)){
  print(i)
  # -----
  if(!is.null(xcmsItems_neg[[i]])) xset_neg <- c(xset_neg, xcmsItems_neg[[i]])else{print(xcmsItems_neg[[i]])}
  if(!is.null(xcmsItems_pos[[i]])) xset_pos <- c(xset_pos, xcmsItems_pos[[i]])else{print(xcmsItems_pos[[i]])}
}

# --------------------------------------

print(fn$paste("Grouping..."))

xset_neg@peaks

xset_neg_g <- group(xset_neg,
                    method="mzClust",
                    mzppm=2)
xset_neg_g@groups

xset_pos_g <- group(xset_pos,
                    method="mzClust")
xset_pos_g@groups

## ====================================

xset_pos_g <- group(xset_pos,
                    method="nearest")
xset_pos_g@groups

xset_neg_g <- group(xset_neg,
                    method="nearest")
xset_neg_g@groups


xset_pos_g_clust <- group(xset_pos,
                         method="mzClust")
xset_pos_g_clust@groups

xset_pos_g_density <- group(xset_pos,
                          method="density")
xset_pos_g_density@groups


xset_neg_g <- group(xset_neg, 
                    method="density",
                    mzppm=2)

xset_pos_g@groups

# --------------------------------------

print(fn$paste("Filling missing peaks..."))

xset_pos_gf <-fillPeaks(xset_pos_g, 
                        method="MSW")

xset_neg_gf <-fillPeaks(xset_neg_g, 
                        method="MSW")

# --------------------------------------

# mzCalib <- calibrate(xset2, calibrants,method="linear", mzabs=0.0001, mzppm=5, neighbours=3, plotres=FALSE)

# --------------------------------------

dat_pos <- as.data.table(
  groupval.joanna(xset_pos_gf@peaks,
                  xset_pos_gf@groups,
                  xset_pos_gf@groupidx,
                  sampnames(xset_pos_gf),
                  method = "maxint", 
                  value = "index", 
                  intensity = "into")
)

dat_neg <- as.data.table(
  groupval.joanna(xset_neg_gf@peaks,
                  xset_neg_gf@groups,
                  xset_neg_gf@groupidx,
                  sampnames(xset_neg_gf),
                  method = "maxint", 
                  value = "index", 
                  intensity = "into")
)

outlist_pos <- dat_pos
outlist_neg <- dat_neg

print(fn$paste("Saving results :-)"))

# --------------------------------------

sheetpath = file.path(folder, "File names vs Card IDs.xlt")

sample_tab <- as.data.table(read.xlsx(file = sheetpath, 
                                      sheetIndex = 1))[, 1:2]

colnames(sample_tab) <- c("File.Name", "Card.ID")

qcloc <- grep("QC|Quality",
              x = levels(sample_tab$Card.ID)
              )

levels(sample_tab$Card.ID)[qcloc] <- "QC"

## --------------------------------------------

to_replace_name_idx <- grep(colnames(outlist_pos), pattern = "RES_BSP")
to_replace_names <- colnames(outlist_pos)[to_replace_name_idx]
new_names <- as.character(factor(to_replace_names, levels = sample_tab$File.Name, labels = sample_tab$Card.ID))
colnames(outlist_pos)[to_replace_name_idx] <- new_names
unique_new_names <- ave(new_names, new_names, FUN=function(x) if (length(x)>1) paste0(x[1], "_", seq_along(x)) else x[1])
colnames(outlist_pos)[to_replace_name_idx] <- unique_new_names

## --------------------------------------------

to_replace_name_idx <- grep(colnames(outlist_neg), pattern = "RES_BSP")
to_replace_names <- colnames(outlist_neg)[to_replace_name_idx]
new_names <- as.character(factor(to_replace_names, levels = sample_tab$File.Name, labels = sample_tab$Card.ID))
colnames(outlist_neg)[to_replace_name_idx] <- new_names
unique_new_names <- ave(new_names, new_names, FUN=function(x) if (length(x)>1) paste0(x[1], "_", seq_along(x)) else x[1])
colnames(outlist_neg)[to_replace_name_idx] <- unique_new_names

# --------------------------------------------

# filter...
save(outlist_pos, 
     outlist_neg, 
     xset_neg, 
     xset_pos,
     xset_neg_g,
     xset_pos_g,
     xset_neg_gf,
     xset_pos_gf,
     file="/Users/jwolthuis/Documents/outlists.turkey.RData")

save(outlist_pos, 
     file="/Users/jwolthuis/Documents/outlist.pos.turkey.RData")

save(outlist_neg, 
     file="/Users/jwolthuis/Documents/outlist.neg.turkey.RData")

print("Done! <3")