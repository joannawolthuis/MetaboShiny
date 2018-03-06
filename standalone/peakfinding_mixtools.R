library(mixtools)

test_file = "~/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/average_pklist/BR1-1_pos.RData"
load(test_file)

head(sum_pos)

data(faithful)
attach(faithful)

plot(sum_pos[180:230,])

hist(waiting, main="Time between Old Faithful eruptions",
     xlab="Minutes", ylab="", cex.main=1.5, cex.lab=1.5, cex.axis=1.4)

wait1 <- normalmixEM(waiting, lambda = .5, mu = c(55, 80), sigma = 5)
plot(wait1, which = 2, cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.8,
     main2 = "Time between Old Faithful eruptions", xlab2 = "Minutes")

wait2 <- spEMsymloc(waiting, mu0 = c(55, 80))
plot(wait2, lty = 2, newplot = FALSE, addlegend = FALSE)

test1 <- normalmixEM(sum_pos[180:230,]
                     , lambda = .5, mu = c(55, 80), sigma = 5)

plot(test1, which=2)

test1$ft

# === second package ===
library(Peaks)
library.dynam('Peaks', 'Peaks', lib.loc=NULL) 

test <- SpectrumSearch(sum_pos[100:2000,],
                       background = T)

x = as.numeric(sum_pos)

ll <- with(rle(x == 0), {
  ifelse(x == 0 & (seq_along(x) != cumsum(lengths)[lengths <= 3 & values]), NA, x)
})

nonzero <- split(x, with(rle(is.na(ll)), rep(1:length(lengths), lengths) + ll * 0))

test$pos
plot(test$y)

# --- PREPROCESSING ?? ---

library(MALDIquant)
library(MassSpecWavelet)
library(data.table)
library(pbapply)
library(parallel)

specdir = "~/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/spectra"
dir.create(specdir,
           showWarnings = F)

# --- load sample names ---

sn <- data.table::fread(file.path("~/Documents/umc/data/Data/BrazilAndSpain/MZXML",
                                  "sampleNames.txt"))

sn$batch <- as.numeric(as.factor(gsub(sn$File_Name, 
                                      pattern = "_\\d\\d\\d$", 
                                      replacement="")))

qcs <- sn[grep(Sample_Name, pattern="QC")]
qc_fns <- paste0("~/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/pklist/",
                 qcs$File_Name, ".RData")

# load qcs

spectra_qc <- pbapply::pblapply(qc_fns, FUN=function(f){
  print(f)
  load(f)
  spec <- createMassSpectrum(mass=as.numeric(rownames(pklist$pos)),
                             intensity=as.numeric(pklist$pos),
                             metaData=list(name=gsub(basename(f),
                                                     pattern = "\\mzXML",
                                                     replacement = "")))
  # neg <- createMassSpectrum(mass=as.numeric(rownames(pklist$neg)),
  #                                        intensity=as.numeric(pklist$neg),
  #                                        metaData=list(name=gsub(basename(f),
  #                                                                pattern = "\\mzXML",
  #                                                                replacement = "")))
  pklist <<- NULL # free memory??
  spec
})

# -------------------------

# smoothing
# smoothed <- smoothIntensity(subset, method="SavitzkyGolay",
#                             halfWindowSize=10)

aligned_qc <- alignSpectra(spectra_qc,
                           halfWindowSize=20,
                           SNR=2,
                           tolerance=2E-6,
                           warpingMethod="lowess")

groups_qc <- as.factor(qcs$Sample_Name)

# average

averaged_qc <- averageMassSpectra(aligned_qc, 
                                  labels=groups_qc,
                                  method="mean")

# peaks?
# peaks <- detectPeaks(averaged[[1]], method="MAD",
#                      halfWindowSize=20, SNR=10)
# 
# peaks <- binPeaks(peaks, tolerance=2E-6)
# points(peaks, col="red", pch=4)

# peaks in here


stopCluster(cl)
cl = parallel::makeCluster(3, "FORK")

peaks_qc <- pblapply(averaged_qc, cl=cl, FUN=function(qc){
  df <- matrix(data = qc@intensity, 
               ncol=1)
  
  rownames(df) <- qc@mass
  msw_peaks <- MassSpecWavelet::peakDetectionCWT(df,
                                                 SNR.Th = 3,
                                                 nearbyPeak = T,
                                                 peakThr = 2000)
  
  peakInfo <- msw_peaks
  majorPeakInfo = peakInfo$majorPeakInfo
  peakIndex <- majorPeakInfo$potentialPeakIndex
  # MassSpecWavelet::plotPeak(test_df, 
  #                           peakIndex, 
  #                           main = paste('Identified peaks with SNR >', 3)) 	
  betterPeakInfo <- tuneInPeakInfo(df, majorPeakInfo)
  # ----------------
  createMassPeaks(mass = as.numeric(rownames(df)[betterPeakInfo$peakCenterIndex]), 
                  intensity = betterPeakInfo$peakValue)
})

save(spectra_qc, 
     aligned_qc, 
     averaged_qc, 
     peaks_qc, 
     file=paste0("~/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/spectra/qcs.RData"))

# USE LATER TO CALIBRATE
reference_peaks <- referencePeaks(peaks_qc, 
                                  method="strict", 
                                  minFrequency = 0.9, 
                                  tolerance=2E-6)

# =================================================================


library(MALDIquant)
library(MassSpecWavelet)
library(data.table)
library(pbapply)
library(parallel)

# LOAD

load("~/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/spectra/qcs.RData")

# DO PCA

library(MetaboAnalystR)

# --- create csv from peaks obj ---

qc_peaktables = lapply(names(peaks_qc), FUN=function(sample){
  peaks = peaks_qc[[sample]]
  tbl <- data.table(mz = peaks@mass,
                    intensity = peaks@intensity)
  #names(tbl)[2] <- sample
  tbl
})

peaktable_qc = qc_peaktables[1]

for(i in 2:length(names(peaks_qc))){
  pktable = qc_peaktables[[i]]
  peaktable_qc <<- merge(peaktable_qc, 
                         pktable, 
                         by="mz",
                         all=T)
}

names(peaktable_qc) = c("mzmed", names(peaks_qc))

options(digits=22)

# FILL MISSING?

library(parallel)
library(snow)
library(doParallel)
library(doSNOW)

outpgrlist <- as.data.table(peaktable_qc)

for(f in list.files("~/Google Drive/MetaboShiny/miniapps/MassChecker/scripts/AddOnFunctions", full.names = T)) source(f)

mzvals = unique(outpgrlist$mzmed)

cl <<- makeSOCKcluster(3,
                       outfile="~/mclog.txt")
registerDoSNOW(cl)

# ---------------

thresh = 2000
resol = 140000
scanmode = "positive"
# withProgress(message = "Filling missing values...",{
#   progress <<- function(i) setProgress(value = (i / length(mzvals))/2,
#                                        detail = "Generating reference values...")
# opts <- list(progress=progress)
cfun <<- function(a, b) c(a, b)
gaussians <- foreach(i=1:length(mzvals), .export = c("mzvals", 
                                                     "thresh", 
                                                     "resol", 
                                                     "scanmode", 
                                                     "searchMZRange",
                                                     "fitGaussianInit",
                                                     "generateGaussian",
                                                     "fitGaussian",
                                                     "getFwhm",
                                                     "getSD",
                                                     "optimizeGauss",
                                                     "getArea",
                                                     "fit2G_2",
                                                     "fit4G_2",
                                                     "fit4peaks",
                                                     "fitG_2",
                                                     "fit3G_2",
                                                     "fit1Peak",
                                                     "fit2peaks",
                                                     "fit3peaks",
                                                     "getFitQuality",
                                                     "checkOverlap",
                                                     "isWithinXppm",
                                                     "sumCurves"),
                     .packages = "xcms",
                     .verbose = T,
                     .inorder = T, # !!! important !!!
                     .combine="c") %dopar% {
                       area = generateGaussian(mzvals[i],
                                               thresh,
                                               resol,
                                               FALSE,
                                               scanmode,
                                               int.factor=1*10^5,1,1)$area
                       area
                     }
# --- prep ---
ref_table <- data.table(mz = mzvals,
                        int = gaussians)
print(head(ref_table))
setkey(ref_table, mz)
setkey(outpgrlist, mzmed)

# --- loop 2 ---

# loop and fill
filled_row <- sapply(row[2:length(row)], FUN=function(samp){
  print(samp)
  if(is.na(samp)){
    rand_gauss <- rnorm(n=1, 
                        mean=area, 
                        sd=0.25*area)
    rand_gauss
  } else{
    samp}
})

outrows <- foreach(i=1:nrow(outpgrlist), 
                   .export = c("outpgrlist",
                               "ref_table"),
                   .packages="data.table",
                   .verbose = T,
                   .inorder = T) %dopar% {
                     row = outpgrlist[i,]
                     rowstart = row[,1]
                     gauss <- ref_table[mz == row[,1],]
                     area = unlist(gauss$int)
                     # loop and fill
                     filled_row <- sapply(row[,2:length(row)], FUN=function(samp){
                       print(samp)
                       if(is.na(samp)){
                         rand_gauss <- rnorm(n=1, 
                                             mean=area, 
                                             sd=0.25*area)
                         rand_gauss
                       } else{
                         samp}
                     })
                     res <- cbind(rowstart, 
                                  t(filled_row))
                     # # --- return ---
                     as.data.table(res)
                   }

# stop cluster
stopCluster(cl)

# combine
outpgrlist <- rbindlist(outrows)
colnames(outpgrlist)
csv <- t(outpgrlist)
colnames(csv) = csv[1,]
csv = csv[-1,]

batch_qc = sapply(rownames(csv), FUN=function(x){
  qc_num = as.numeric(gsub(x, pattern="QC", replacement=""))
  if(qc_num > 0 & qc_num < 7){
    1
  }else if(qc_num > 7 && qc_num < 13){
    2
  }else{
    3
  }
})

blocky <- function(table, dim=10) table[1:dim, 1:dim]

metabo_csv = as.data.table(cbind(Sample= rownames(csv), 
                                 Label=batch_qc, 
                                 csv))
blocky(metabo_csv)
fwrite(metabo_csv, file="~/qcs.csv")

# === COMBAT NORMALIZATION ===

csv_pheno <- data.frame(sample = 1:nrow(metabo_csv),
                        #outcome = rep(1, nrow(metabo_csv)),
                        batch = metabo_csv$Label)

# mod = model.matrix(~as.factor(outcome), 
#                    data=csv_pheno)

# parametric adjustment

csv_edata <-as.data.frame(t(metabo_csv[,!c(1,2)]))
csv_edata <- data.matrix(csv_edata, rownames.force = NA)

csv_edata
colnames(csv_edata) <- metabo_csv$Sample

batch_normalized = sva::ComBat(dat=csv_edata, batch=csv_pheno$batch)

metabo_csv_combat = as.data.table(cbind(Sample= colnames(batch_normalized), 
                                        Label=batch_qc, 
                                        t(batch_normalized)))

blocky(metabo_csv_combat)  
fwrite(metabo_csv_combat, file="~/qcs_combat.csv")

# === T-SNE VISUALIZATION ===

library(Rtsne)
library(plotly)
library(ggplot2)

csv = t(csv_edata)
csv_cb = t(batch_normalized)
blocky(csv)
t.sne_before = Rtsne(csv, dims = 3, perplexity=4)
t.sne_after = Rtsne(csv_cb, dims = 3, perplexity=4)

# plot tsne

df <- t.sne_before$Y
df <- t.sne_after$Y

fac.lvls <- batch_qc
chosen.colors <<- rainbow(length(unique(fac.lvls)))
# --- add ellipses ---
classes <- as.factor(as.numeric(fac.lvls))

plots <- plotly::plot_ly()

for(class in levels(classes)){
  row = which(classes == class)
  print(row)
  # ---------------------
  xc=df[row, 1]
  yc=df[row, 2]
  zc=df[row, 3]
  # --- plot ellipse ---
  o <- rgl::ellipse3d(cov(cbind(xc,yc,zc)), 
                      centre=c(mean(xc), 
                               mean(yc), 
                               mean(zc)), 
                      level = 0.95)
  mesh <- c(list(x = o$vb[1, o$ib]/o$vb[4, o$ib], 
                 y = o$vb[2, o$ib]/o$vb[4, o$ib], 
                 z = o$vb[3, o$ib]/o$vb[4, o$ib]))
  plots = plots %>% add_trace(
    x=mesh$x, 
    y=mesh$y, 
    z=mesh$z, 
    type='mesh3d',
    alphahull=0,
    opacity=0.1
  )
}

adj_plot <- plotly_build(plots)
rgbcols <- toRGB(chosen.colors)
c = 1
for(i in seq_along(adj_plot$x$data)){
  print(i)
  item = adj_plot$x$data[[i]]
  if(item$type == "mesh3d"){
    adj_plot$x$data[[i]]$color <- rgbcols[c]
    adj_plot$x$data[[i]]$visible <- TRUE
    c = c + 1
  }
}

library(gsubfn)

# --- return ---
pca_plot <<- adj_plot %>% add_trace(
  hoverinfo = 'text',
  text = rownames(df),
  x = df[,1], 
  y = df[,2], 
  z = df[,3], 
  visible = rep(T, times=length(unique(fac.lvls))),
  type = "scatter3d",
  opacity=1,
  color= c(classes), colors=chosen.colors
) %>%  layout(scene = list(
  xaxis = list(
    title = fn$paste("Dim1")),
  yaxis = list(
    title = fn$paste("Dim2")),
  zaxis = list(
    title = fn$paste("Dim3")))) 
# --- return ---

pca_plot

# =================================================================

# brs <- sn[grep(Sample_Name, pattern="BR")]
# br_fns <- paste0("~/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/pklist/",
#                  brs$File_Name, ".RData")
sn <- data.table::fread(file.path("~/Documents/umc/data/Data/BrazilAndSpain/MZXML",
                                  "sampleNames.txt"))

sn$batch <- as.numeric(as.factor(gsub(sn$File_Name, 
                                      pattern = "_\\d\\d\\d$", 
                                      replacement="")))

all_fns <- paste0("~/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/pklist/",
                  sn$File_Name, ".RData")

cl = parallel::makeCluster(3, "FORK")

# --- LOAD IN ALL SAMPLES ---

repl.pattern = split(sn,f = sn$Sample_Name)

library(MALDIquant)
library(MassSpecWavelet)

averaged <- pbapply::pblapply(repl.pattern, cl=cl, FUN=function(repl.set){
  repl.spec <- lapply(1:nrow(repl.set), FUN=function(i){
    f = paste0("~/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/pklist/",
               repl.set[i, "File_Name"], ".RData")
    load(f)
    # ----- create mass spectra -----
    pos <- createMassSpectrum(mass=as.numeric(rownames(pklist$pos)),
                              intensity=as.numeric(pklist$pos),
                              metaData=list(name=gsub(basename(f),
                                                      pattern = "\\mzXML",
                                                      replacement = "")))
    neg <- createMassSpectrum(mass=as.numeric(rownames(pklist$neg)),
                              intensity=as.numeric(pklist$neg),
                              metaData=list(name=gsub(basename(f),
                                                      pattern = "\\mzXML",
                                                      replacement = "")))
    list(pos = pos, 
         neg = neg)
  })
  print(repl.spec)
  # average samples
  averaged_pos <- averageMassSpectra(lapply(repl.spec, function(x) x$pos),
                                 labels=repl.set$Sample_Name,
                                 method="mean")
  averaged_neg <- averageMassSpectra(lapply(repl.spec, function(x) x$neg),
                                     labels=repl.set$Sample_Name,
                                     method="mean")
  print("here")
  # call peaks
  fn_pos = paste0("~/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/spectra/",
                  unique(repl.set$Sample_Name), "_pos.RData")
  fn_neg = paste0("~/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/spectra/",
         unique(repl.set$Sample_Name), "_neg.RData")
  print(fn_pos)
  save(x=averaged_pos, file=fn_pos)
  save(x=averaged_neg, file=fn_neg)
})

library(pbapply)

parallel::stopCluster(cl)
cl = parallel::makeCluster(3, "FORK")

spec_files = list.files("~/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/spectra/",
                        pattern = paste0(mode, "\\.RData"),
                        full.names = T)

pblapply(spec_files, cl=cl, FUN=function(f){
  # --------------== mode get ==--------------><
  mode = if(grepl("pos",f)) "pos" else "neg"  #|
  # ------------------------------------------><
  sampname = gsub(basename(f), 
                  pattern="_(pos|neg)\\.RData", 
                  replacement = "")
  fn = paste0("~/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/specpks/",
              sampname, "_", mode,".RData")
  if(file.exists(fn)) return()
  # -----------------------------------------------------------
  load(f)
  switch(mode,
         pos = {df <- matrix(data = averaged_pos[[1]]@intensity, 
                             ncol=1)
         rownames(df) <- averaged_pos[[1]]@mass
         },
         neg = {df <- matrix(data = averaged_neg[[1]]@intensity, 
                             ncol=1)
         rownames(df) <- averaged_neg[[1]]@mass
         })
  try({
    msw_peaks <- MassSpecWavelet::peakDetectionCWT(df,
                                                   SNR.Th = 3,
                                                   nearbyPeak = T,
                                                   peakThr = 2000)
    
    peakInfo <- msw_peaks
    majorPeakInfo = peakInfo$majorPeakInfo
    peakIndex <- majorPeakInfo$potentialPeakIndex
    betterPeakInfo <- tuneInPeakInfo(df, majorPeakInfo)
    # ----------------
    peaks = createMassPeaks(mass = as.numeric(rownames(df)[betterPeakInfo$peakCenterIndex]),
                            intensity = betterPeakInfo$peakValue,
                            metaData = list(sample = sampname))
    save(x=peaks, 
         file=fn)
  })
})


