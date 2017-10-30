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

# =======================================================

folder = "/Users/jwolthuis/PROCESSING_HPC/"

sample_tab <- fread(file = "/Users/jwolthuis/Turkey_OKfiles.csv"
)[,2:3]

ok_filenames <- c(sample_tab[judgement == "OK", filenames])

filenames <- sapply(ok_filenames, FUN=function(fn){
  paste0(folder, fn, sep="")
})

# hpc : 
nslots <- Sys.getenv( "NSLOTS" )
nslots = parallel::detectCores()

print(fn$paste("Doing MSW peak calling with $nslots slots!"))

#bpparam = SnowParam(workers = nslots, type = "SOCK")

cl = parallel::makeCluster(4, "FORK", outfile="")
stopCluster(cl)

xcmsItems_pos <- pblapply(seq_along(filenames), cl=cl, FUN=function(i){
  # --------------------------------------
  print(i)
  f <- filenames[i]
  # # --------------------------------------
  # while( is.null(rawF) && attemptA <= 3 ) {
  # rawF <- NULL
  # attemptA <- 1
  #   attemptA <- attemptA + 1
  #   try(
  #     rawF <- xcmsRaw(f)
  #   )
  # }
  #plot(rawF@tic)
  #attemptB <- 0
  # --------------------------------------
  xset = NULL
  startscans = c(1, 10, 20)
  # --------------------------------------
  while( is.null(xset) && attemptB <= 3 ) {
    attemptB <- attemptB + 1
    try(
      xset <- xcmsSet(f,
                      method="MSW",
                      scanrange = c(startscans[attemptB]:150)
      )  
    )
  }
  # --- return ---
  xset
})


xcmsItems_neg <- pblapply(seq_along(filenames), cl=cl, FUN=function(i){
  # --------------------------------------
  print(i)
  f <- filenames[i]
  # --------------------------------------
  # rawF <- NULL
  # attemptA <- 1
  # while( is.null(rawF) && attemptA <= 3 ) {
  #   attemptA <- attemptA + 1
  #   try(
  #     rawF <- xcmsRaw(f)
  #   )
  # }
  # plot(rawF@tic)
  # --------------------------------------
  xset = NULL
  attemptB <- 0
  endscans = c(310, 305, 300)
  # --------------------------------------
  while( is.null(xset) && attemptB <= 3 ) {
    attemptB <- attemptB + 1
    try(
      xset <- xcmsSet(f,
                      method="MSW",
                      scanrange = c(151:endscans[attemptB])
      )  
    )
  }
  # --- return ---
  xset
})

# --- why are some of these null??? SCANRANGE ---

null_names <- which(sapply(xcmsItems_pos, is.null))
null_filenames <- filenames[null_names]
rawF <- xcmsRaw(null_filenames[1])
xcmsItems_pos_nully <- pblapply(seq_along(null_filenames), cl=FALSE, FUN=function(i){
  # --------------------------------------
  print(i)
  f <- null_filenames[i]
  rawF <- NULL
  attemptA <- 1
  # --------------------------------------
  while( is.null(rawF) && attemptA <= 3 ) {
    attemptA <- attemptA + 1
    try(
      rawF <- xcmsRaw(f)
     )
   }
  # --------------------------------------
  plot(rawF@tic)
  print(rawF@tic)
  attemptB <- 1
  xset = NULL
  startscans = c(1, 10, 20)
  # --------------------------------------
  while( is.null(xset) && attemptB <= 3 ) {
    attemptB <- attemptB + 1
    try(
      xset <- xcmsSet(f,
                      method="MSW",
                      scanrange = c(startscans[attemptB]:145),
                      scales=c(1,4,9)
      )  
    )
  }
  # --- return ---
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

for(xs in xcmsItems_pos_nully){
  xset_pos <- c(xset_pos, xs)
}

# --------------------------------------

print(fn$paste("Grouping..."))

xset_pos_g <- group(xset_pos, 
                    method="mzClust",
                    mzppm=2)

xset_neg_g <- group(xset_neg, 
                    method="mzClust",
                    mzppm=2)

# --------------------------------------

# nslots <- Sys.getenv( "NSLOTS" )
# nslots = parallel::detectCores()
# print(fn$paste("Doing MSW peak filling with $nslots slots!"))
# bpparam = SnowParam(workers = nslots, type = "SOCK")

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