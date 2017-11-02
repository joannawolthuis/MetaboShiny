#$ -S /home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -l h_rt=72:00:00
#$ -l h_vmem=30G
#$ -m beas

args = commandArgs(trailingOnly=TRUE)
folder = args[1]

# --- load libraries ---

library(xcms)
library(snow)
library(BiocParallel)
library(gsubfn)
library(data.table)

# --- set defaults ---

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

scale <- seq(1,64,2)

files = list.files(folder, pattern="\\.RData", full.names = T)

# ----------------------------

load(files[[1]])

xset_pos_tot <- xset_pos
xset_neg_tot <- xset_neg

for(i in 2:length(files)){
    print(i)
    # --- RELOAD ---
    load(files[[i]])
    # --------------
    if(!is.null(xset_neg)) xset_neg_tot <- c(xset_neg_tot, xset_neg) else{print(xset_neg)}
    if(!is.null(xset_pos)) xset_pos_tot <- c(xset_pos_tot, xset_pos) else{print(xset_pos)}
}

# --------------------------------------

print("Grouping...")

xset_neg_g <- group(xset_neg_tot,
method="mzClust",
mzppm=2)

xset_pos_g <- group(xset_pos_tot,
method="mzClust",
mzppm=2)

# --------------------------------------

print("Filling missing peaks...")

xset_pos_gf <-fillPeaks(xset_pos_g,
method="MSW")

xset_neg_gf <-fillPeaks(xset_neg_g,
method="MSW")

# --------------------------------------

save(outlist_pos,
     outlist_neg,
     xset_neg,
     xset_pos,
     xset_neg_g,
     xset_pos_g,
     xset_neg_gf,
     xset_pos_gf,
     file=file.path(folder,"ritems.RData"))

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

print("Saving results :-)")


save(outlist_pos,
     file=file.path(folder,"outlist.pos.RData"))

save(outlist_neg,
     file=file.path(folder, "outlist.neg.RData"))

print("Done! Please do the renaming manually just to be sure <3")
