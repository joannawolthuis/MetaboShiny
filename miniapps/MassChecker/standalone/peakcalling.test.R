# PEAK CALLING DIFFERENT METHODS...

library(gsubfn)

thresh = list(pos=2000,
              neg=2000)
snr = 3
resol = 140000
outdir = normalizePath("~/Documents/umc/data/Data/BrazilAndSpain/MZXML")
scriptdir = normalizePath("~/Google Drive/MetaboShiny/miniapps/MassChecker/scripts")
resdir = file.path(outdir, "results")
dimsThresh = 100
trim = 0.1
cores = 3
ppm = 2
pos_scans = 10:145
neg_scans = 165:300

# RAWFILE
files = list.files(outdir, pattern="\\.raw",full.names = T)
rawfile <- files[[1]]

# 0: SOURCE...

spectra <- pblapply(files[1:3], FUN=function(rawfile){
  cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/scantimes' -raw '$rawfile'")
  scantimes <- as.numeric(system(cmd,intern = T))
  
  # Generate a matrix
  cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/scanpolarity' -raw '$rawfile'")
  polarity <- system(cmd,intern = T)
  
  # Get time values for positive and negative scans
  posTimes <- scantimes[polarity == "positive"]
  negTimes <- scantimes[polarity == "negative"]
  
  # Generate an index with which to select values for each mode
  posInd <- which(scantimes %in% posTimes)
  negInd <- which(scantimes %in% negTimes)
  
  # get mzvals
  pos = lapply(posInd, FUN=function(scan){
    cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/printspectrum' -raw '$rawfile' -sn $scan")
    spec <- system(cmd,
                   intern = T)
    temp.list <- strsplit(spec, " ")
    #temp.list <- lapply(temp.list, FUN=function(x) if(length(x) == 2) x else NULL)
    mzvals = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[1]])))
    intensities = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[2]])))
    # --- mass spectrum obj ---
    data.table(mz = mzvals, 
               i = intensities)
  })
  neg = lapply(negInd, FUN=function(scan){
    cmd = gsubfn::fn$paste("'$scriptdir/ms/unthermo/tools/printspectrum' -raw '$rawfile' -sn $scan")
    spec <- system(cmd,intern = T)
    temp.list <- strsplit(spec, " ")
    #temp.list <- lapply(temp.list, FUN=function(x) if(length(x) == 2) x else NULL)
    mzvals = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[1]])))
    intensities = as.numeric(unlist(lapply(temp.list, FUN=function(x) x[[2]])))
    # --- mass spectrum obj ---
    data.table(mz = mzvals, 
               i = intensities)
  })
  
  # Separate each mode into its own matrix
  posY <- rbindlist(pos)
  negY <- rbindlist(neg)
  # --- return ---
  list(pos = posY, neg=negY)
})


# 0.1 - binned according to fwhm

# === GENERATE FWHM BREAK RANGE ===
mzrange = c(70, 600)

lowMZ = mzrange[1]
highMZ = mzrange[2]

nsegment = 2*(highMZ-lowMZ)
segment = seq(from=lowMZ, to=highMZ, length.out=nsegment+1)
breaks.fwhm=NULL
breaks.fwhm.avg=NULL
# for (i in 1:2) {
stopCluster(cl)
cl = parallel::makeCluster(3, "FORK")

breaks.fwhm <- unlist(pbsapply(1:nsegment, cl=cl, FUN=function(i){
  startsegm <- segment[i]
  endsegm <- segment[i+1]
  resol.mz <- resol*(1/sqrt(2)^(log2(startsegm/200)))
  fwhmsegm <- startsegm/resol.mz
  break.fwhm <- seq(from=(startsegm+fwhmsegm),to=endsegm, by=0.2*fwhmsegm)
  # --- return ---
  break.fwhm
}))

breaks.fwhm.avg <- unlist(pbsapply(1:nsegment, cl=cl, FUN=function(i){
  startsegm <- segment[i]
  endsegm <- segment[i+1]
  resol.mz <- resol*(1/sqrt(2)^(log2(startsegm/200)))
  fwhmsegm <- startsegm/resol.mz
  breaks.fwhm <- c(breaks.fwhm, seq(from=(startsegm+fwhmsegm),to=endsegm, by=0.2*fwhmsegm))
  #breaks.fwhm <- c(breaks.fwhm, seq(from=(startsegm), to=endsegm, by=0.2*fwhmsegm))
  
  # average the m/z instead of start value
  range = seq(from=(startsegm+fwhmsegm),to=endsegm, by=0.2*fwhmsegm)
  deltaMZ = range[2]-range[1] 
  break.fwhm.avg <- range + 0.5 * deltaMZ
  # --- return ---
  break.fwhm.avg
}))

# ==============================

# AGGREGATE 

fwhm_binned <- pblapply(spectra, FUN=function(spec){
  pos_spec <- spec$pos
  neg_spec <- spec$neg
  setkey(pos_spec, mz)
  setkey(neg_spec, mz)
  
  posY <<- pos_spec[mz %between% mzrange]
  negY <<- neg_spec[mz %between% mzrange]
  
  posY[, bin:= cut(mz, breaks.fwhm, include.lowest=TRUE, right=TRUE, labels=FALSE)]
  negY[, bin:= cut(mz, breaks.fwhm, include.lowest=TRUE, right=TRUE, labels=FALSE)]
  
  #yp2 <- cut(posY$mz, breaks.fwhm, include.lowest=TRUE, right=TRUE, labels=FALSE)
  #yn <- cut(negY$mz, breaks.fwhm, include.lowest=TRUE, right=TRUE, labels=FALSE)
  
  # Empty the bins
  posBins<- rep(0,length(breaks.fwhm)-1)
  
  negBins<- rep(0,length(breaks.fwhm)-1)
  
  
  # Get the list of intensity values for each bin, and add the
  # intensity values which are in the same bin
  
  if (nrow(posY) > 0) {
    posY[, aggr := sum(i), by=bin]
    ap <- posY[!is.na(bin),c("bin", "aggr")]
    posBins[ap$bin] <- ap$aggr
  }
  if (nrow(negY) > 0) {
    negY[, aggr := sum(i), by=bin]
    an <- negY[!is.na(bin),c("bin", "aggr")]
    negBins[an$bin] <- an$aggr
  }
  
  # Zero any values that are below the threshold
  posBins[posBins < thresh$pos] <- 0
  negBins[negBins < thresh$neg] <- 0
  
  posRes = posBins
  negRes = negBins
  
  posRes = t(posRes)
  negRes = t(negRes)
  
  label=gsub(".+/(.+)\\.raw$", "\\1", rawfile, ignore.case=TRUE)
  
  # Add in file names as row names
  rownames(posRes) = label
  rownames(negRes) = label
  
  # Add 0.5 to the values in breaks.fwhm, and delete the last value
  a <- breaks.fwhm.avg[-length(breaks.fwhm.avg)]  # + 0.5*deltaMZ
  
  # Format as string and show precision of float to 2 digits
  b <- sprintf("%.5f",a)
  
  # Use this as the column names
  colnames(posRes) <- b
  colnames(negRes) <- b
  
  posResT <- t(posRes)
  negResT <- t(negRes)
  
  pos_spec_fwhm = MALDIquant::createMassSpectrum(as.numeric(rownames(posResT)), c(posResT), metaData = list(sample=label))
  neg_spec_fwhm = MALDIquant::createMassSpectrum(as.numeric(rownames(negResT)), c(negResT), metaData = list(sample=label))
  # --- return ---
  list(pos = pos_spec_fwhm, 
       neg = neg_spec_fwhm)
})

# ==============================
# 0.2 - unbinned

un_binned <- pblapply(spectra, FUN=function(spec){
  pos_spec <- MALDIquant::createMassSpectrum(as.numeric(posY$mz), as.numeric(posY$i), metaData = list(sample="TEST"))
  neg_spec <- MALDIquant::createMassSpectrum(as.numeric(negY$mz), as.numeric(negY$i), metaData = list(sample="TEST"))
  list(pos = pos_spec, 
       neg = neg_spec)
})

# JOINING SAMPLES

fwhm_triplet_pos <- lapply(fwhm_binned, function(x) x$pos)
fwhm_triplet_neg <- lapply(fwhm_binned, function(x) x$neg)
  
unb_triplet_pos <- lapply(un_binned, function(x) x$pos)
unb_triplet_neg <- lapply(un_binned, function(x) x$neg)

# ALIGNING

fwhm_aligned_pos <- MALDIquant::alignSpectra(fwhm_triplet_pos)
fwhm_aligned_neg <- MALDIquant::alignSpectra(fwhm_triplet_neg)

unb_aligned_pos <- MALDIquant::alignSpectra(unb_triplet_pos)
unb_aligned_neg <- MALDIquant::alignSpectra(unb_triplet_neg)

# AVERAGING

avg_fwhm_algn_pos <- MALDIquant::averageMassSpectra(fwhm_aligned_pos, method="mean")
avg_fwhm_algn_neg <- MALDIquant::averageMassSpectra(fwhm_aligned_neg, method="mean")

avg_fwhm_pos <- MALDIquant::averageMassSpectra(fwhm_triplet_pos, method="mean")
avg_fwhm_neg <- MALDIquant::averageMassSpectra(fwhm_triplet_neg, method="mean")

# -------------------
avg_unb_algn_pos <- MALDIquant::averageMassSpectra(unb_aligned_pos, method="mean")
avg_unb_algn_neg <- MALDIquant::averageMassSpectra(unb_aligned_neg, method="mean")

avg_unb_pos <- MALDIquant::averageMassSpectra(unb_triplet_pos, method="mean")
avg_unb_neg <- MALDIquant::averageMassSpectra(unb_triplet_neg, method="mean")

# SMOOTHING OR NOT??

smo_fwhm_algn_pos <- MALDIquant::smoothIntensity(avg_fwhm_algn_pos)
smo_fwhm_algn_neg <- MALDIquant::smoothIntensity(avg_fwhm_algn_neg)

smo_unb_algn_pos <- MALDIquant::smoothIntensity(avg_unb_algn_pos)
smo_unb_algn_neg <- MALDIquant::smoothIntensity(avg_unb_algn_neg)

smo_fwhm_pos <- MALDIquant::smoothIntensity(avg_fwhm_pos)
smo_fwhm_neg <- MALDIquant::smoothIntensity(avg_fwhm_neg)

smo_unb_pos <- MALDIquant::smoothIntensity(avg_unb_pos)
smo_unb_neg <- MALDIquant::smoothIntensity(avg_unb_neg)

# 1: Generic MALDIQUANT

peaks_maldi_pos <- pblapply(list(FWHM_NOALIGN_NOSMOOTH=avg_fwhm_pos,
                                 FWHM_NOALIGN_SMOOTH=smo_fwhm_pos,
                                 FWHM_ALIGN_NOSMOOTH=avg_fwhm_algn_pos,
                                 FWHM_ALIGN_SMOOTH=smo_fwhm_algn_pos,
                                 UNBINNED_NOALIGN_NOSMOOTH=avg_unb_pos,
                                 UNBINNED_NOALIGN_SMOOTH=smo_unb_pos,
                                 UNBINNED_ALIGN_NOSMOOTH=avg_unb_algn_pos,
                                 UNBINNED_ALIGN_SMOOTH=smo_unb_algn_pos), FUN=function(spec){
                                   MALDIquant::detectPeaks(spec, snr=10)
                                 })

peaks_maldi_pos

# 2: GAUSSIAN

range = c(500.000, 501)


a <- avg_fwhm_pos@mass[which(avg_fwhm_pos@mass %between% range)]
b <- avg_fwhm_pos@intensity[which(avg_fwhm_pos@mass %between% range)]
spec <- MALDIquant::createMassSpectrum(a, b)

plot(a,b)

jumprange <- seq(70, 600, by = 0.01)

for(i in 10100:length(jumprange)){
  print(i)
  range = c(jumprange[i-1], jumprange[i])
  a <- avg_fwhm_pos@mass[which(avg_fwhm_pos@mass %between% range)]
  b <- avg_fwhm_pos@intensity[which(avg_fwhm_pos@mass %between% range)]
  df <- matrix(data = b,
               ncol=1)
  rownames(df) <- a
  try({
    msw_peaks <- MassSpecWavelet::peakDetectionCWT(df,
                                                   SNR.Th = 3,
                                                   nearbyPeak = FALSE,
                                                   scales = c(1, seq(2, 30, 2), seq(32, 64, 4))
    )
    MassSpecWavelet::plotPeak(df, msw_peaks$majorPeakInfo$allPeakIndex)
  })
  Sys.sleep(1)
}

range <- c(jumprange[44], jumprange[45])
a <- avg_fwhm_pos@mass[which(avg_fwhm_pos@mass %between% range)]
b <- avg_fwhm_pos@intensity[which(avg_fwhm_pos@mass %between% range)]
plot(a,b)
spec <- MALDIquant::createMassSpectrum(a, b)


# ------------------------

a <- avg_unb_pos@mass[which(avg_unb_algn_pos@mass %between% range)]
b <- avg_unb_pos@intensity[which(avg_unb_algn_pos@mass %between% range)]
spec <- MALDIquant::createMassSpectrum(a, b)

plot(a,b)

#spec <- avg_fwhm_pos
base_ions_pos <- c(94.080061918,103.075356012,105.054620568,112.08692377,117.054620568,119.070270632,126.136469038,132.101905118,135.109707592,135.120735362,137.06361441,146.081169674,151.079264474,153.077156028,155.033885124,156.07675301,160.096819738,170.092403074,171.168960534,178.115521254,187.071333272,193.034279048,207.14186473,221.157514794,227.055014492,229.216206588,235.173164858,246.169984678,255.22647541,262.164899298,291.235765114,323.144059172,341.00332532,372.310835254,381.367325986,403.360965626,416.082341002,646.086255054,662.081169674,1030.086255054,1046.081169674)
base_ions_neg <- c(92.065509014,101.060803108,103.040067664,110.072370866,115.040067664,117.055717728,124.121916134,130.087352214,133.095154688,133.106182458,135.049061506,144.06661677,149.06471157,151.062603124,153.01933222,154.062200106,158.082266834,168.07785017,169.15440763,176.10096835,185.056780368,191.019726144,205.127311826,219.14296189,225.040461588,227.201653684,233.158611954,244.155431774,253.211922506,260.150346394,289.22121221,321.129506268,338.988772416,370.29628235,379.352773082,401.346412722,414.067788098,644.07170215,660.06661677,1028.07170215,1044.06661677)


for(ion in base_ions_pos){
  range = c(ion - (4E-6 * ion),ion + (4E-6 * ion))
  a <- avg_fwhm_pos@mass[which(avg_fwhm_pos@mass %between% range)]
  b <- avg_fwhm_pos@intensity[which(avg_fwhm_pos@mass %between% range)]
  plot(a,b)
  df <- matrix(data = b,
               ncol=1)
  rownames(df) <- a
  try({
    msw_peaks <- MassSpecWavelet::peakDetectionCWT(df,
                                                   SNR.Th = 3,
                                                   nearbyPeak = FALSE
                                                   ,scales = c(1:3)
    )
    MassSpecWavelet::plotPeak(df, msw_peaks$majorPeakInfo$allPeakIndex)
  })
  Sys.sleep(1)
}

# === FULL SPECTRUM CHECK FOR IS ===

range = c(70, 600)
a <- avg_fwhm_pos@mass[which(avg_fwhm_pos@mass %between% range)]
b <- avg_fwhm_pos@intensity[which(avg_fwhm_pos@mass %between% range)]
plot(a,b)
df <- matrix(data = b,
             ncol=1)
rownames(df) <- a
msw_peaks <- MassSpecWavelet::peakDetectionCWT(df,
                                               SNR.Th = 3,
                                               nearbyPeak = FALSE)

msw_peaks$majorPeakInfo
MassSpecWavelet::plotPeak(df, msw_peaks$majorPeakInfo$peakIndex)
points(1:length(a),b)

# === FIND INTERNAL STANDARDS IN SPECTRUM ===

betterPeakInfo <- tuneInPeakInfo(df, msw_peaks$majorPeakInfo)
peaks = createMassPeaks(mass = as.numeric(rownames(df)[betterPeakInfo$peakCenterIndex]),
                        intensity = betterPeakInfo$peakValue,
                        snr = betterPeakInfo$peakSNR)

peakmat <- data.table(mzmed = peaks@mass)
peakmat <- apply(peakmat, MARGIN = 2, as.numeric)
groups <- mzClustGeneric(as.matrix(peakmat), 
                         mzppm=ppm,
                         shinyprog=FALSE)

for(ion in base_ions_pos){
  range = c(ion - (2E-6 * ion),ion + (2E-6 * ion))
  print(any(groups$mat[,1] %between% range))
}


# SAME FOR GAUSS...

range = as.numeric(df[,1])
names(range) = rownames(df)
range[34:43]
options(digits=16)
int.factor=1*10^5 # Number of x used to calc area under Gaussian (is not analytic) 
scale=2 # Initial value used to estimate scaling parameter

values = list("mean"=NULL,"area"=NULL,"nr"=NULL,"min"=NULL,"max"=NULL,"qual"=NULL,"spikes"=0)

values = searchMZRange(range,
                       values,
                       int.factor,
                       scale,
                       140000,
                       resdir,
                       "test",
                       "positive",
                       FALSE,
                       600,
                       600,
                       2000)  

peaks = createMassPeaks(mass = as.numeric(values$mean),
                        intensity = values$height.pkt)

outlist.persample=cbind("samplenr"=values$nr, 
                        "mzmed.pkt"=values$mean, 
                        "fq"=values$qual, 
                        "mzmin.pkt"=values$min, 
                        "mzmax.pkt"=values$max, 
                        "height.pkt"=values$area)
index = which(outlist.persample[,"height.pkt"]==0)
if (length(index)>0){
  outlist.persample=outlist.persample[-index,]  
}

save(outlist.persample, file=paste(outdir, "specpks", paste(sampname, "_", scanmode, ".RData", sep=""), sep="/"))

message(paste("There were", values$spikes, "spikes!"))

install.packages("mclust")
library(mclust)
fit2 = mclust::Mclust(iris[,1:4])
range = 139000:140000
dim(df)
data <- data.frame(rownames(df)[range],df[range])
plot(data)
values = searchMZRange(df[range,],
                       values,
                       int.factor,
                       scale,
                       140000,
                       resdir,
                       "test",
                       "pos",
                       FALSE,
                       600,
                       600,
                       2000)  

plot(data)
points(values$mean) # 8 peaks??

df2 <- matrix(data = data[,2],
             ncol=1)
rownames(df2) <- data[,1]
msw_peaks <- MassSpecWavelet::peakDetectionCWT(df2,
                                               SNR.Th = 3,
                                               nearbyPeak = FALSE,
                                               scales=c(1:4))
msw_peaks$majorPeakInfo
MassSpecWavelet::plotPeak(df2, msw_peaks$majorPeakInfo$peakIndex)

# === PEAK CALLING COMPARISON ===

peaks_gauss_pos <- pblapply(list(FWHM_NOALIGN_NOSMOOTH=avg_fwhm_pos,
                                 FWHM_NOALIGN_SMOOTH=smo_fwhm_pos,
                                 FWHM_ALIGN_NOSMOOTH=avg_fwhm_algn_pos,
                                 FWHM_ALIGN_SMOOTH=smo_fwhm_algn_pos,
                                 UNBINNED_NOALIGN_NOSMOOTH=avg_unb_pos,
                                 UNBINNED_NOALIGN_SMOOTH=smo_unb_pos,
                                 UNBINNED_ALIGN_NOSMOOTH=avg_unb_algn_pos,
                                 UNBINNED_ALIGN_SMOOTH=smo_unb_algn_pos), FUN=function(spec){
                                   # return MassPeaks obj
                                   NULL
                                 })
# 3: WAVELET DECONVOLUTION
library(MassSpecWavelet)

peaks_wave_pos <- pblapply(list(FWHM_NOALIGN_NOSMOOTH=avg_fwhm_pos,
                                FWHM_NOALIGN_SMOOTH=smo_fwhm_pos,
                                FWHM_ALIGN_NOSMOOTH=avg_fwhm_algn_pos,
                                FWHM_ALIGN_SMOOTH=smo_fwhm_algn_pos,
                                UNBINNED_NOALIGN_NOSMOOTH=avg_unb_pos,
                                UNBINNED_NOALIGN_SMOOTH=smo_unb_pos,
                                UNBINNED_ALIGN_NOSMOOTH=avg_unb_algn_pos,
                                UNBINNED_ALIGN_SMOOTH=smo_unb_algn_pos), cl=cl,FUN=function(spec){
                                  # return MassPeaks obj
                                  df <- matrix(data = spec@intensity,
                                               ncol=1)
                                  rownames(df) <- spec@mass
                                  msw_peaks <- MassSpecWavelet::peakDetectionCWT(df,
                                                                                 SNR.Th = 10,
                                                                                 nearbyPeak = T,
                                                                                 peakThr = 2000)
                                  
                                  peakInfo <- msw_peaks
                                  majorPeakInfo = peakInfo$majorPeakInfo
                                  peakIndex <- majorPeakInfo$potentialPeakIndex
                                  betterPeakInfo <- tuneInPeakInfo(df, majorPeakInfo)
                                  # ----------------
                                  peaks = createMassPeaks(mass = as.numeric(rownames(df)[betterPeakInfo$peakCenterIndex]),
                                                          intensity = betterPeakInfo$peakValue)
                                  peaks
                                })
peaks_wave_pos

library(bladderbatch)
data(bladderdata)


pheno = pData(bladderEset)
# add fake age variable for numeric
pheno$age = c(1:7, rep(1:10, 5))
write.table(data.frame(cel=rownames(pheno), pheno), row.names=F, quote=F, sep="\t", file="bladder-pheno.txt")

edata = exprs(bladderEset)
write.table(edata, row.names=T, quote=F, sep="\t", file="bladder-expr.txt")
# use dataframe instead of matrix
mod = model.matrix(~as.factor(cancer) + age, data=pheno)
t = Sys.time()
cdata = sva::ComBat(dat=edata, batch=as.factor(pheno$batch), mod=mod)
print(Sys.time() - t)
print(cdata[1:5, 1:5])

# ================ compare wavelet and gaussian peak calling ===============

fp = outdir

"/Users/jwolthuis/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/averaged/1-1B_neg.RData"

load("/Users/jwolthuis/Documents/umc/data/Data/BrazilAndSpain/MZXML/results/averaged/1-1B_pos.RData")

tryRange = 1:5000 # good example region

library(MALDIquant)
# smooth or not???
smoothed = smoothIntensity(averaged)
spec = smoothed[[1]]
df <- matrix(data = spec@intensity,
             ncol=1)
rownames(df) <- spec@mass

# ---------------

plot(averaged$`1-1B`@intensity[tryRange])
spec = averaged$`1-1B`
df <- matrix(data = spec@intensity,
             ncol=1)
rownames(df) <- spec@mass


# wavelet
msw_peaks <- MassSpecWavelet::peakDetectionCWT(df,
                                               SNR.Th = 3,
                                               nearbyPeak = FALSE,
                                               #scales =  c(1:15),
                                               peakThr = 100)
length(msw_peaks$majorPeakInfo$allPeakIndex)
# gauss

head(df)
values = list("mean"=NULL,"area"=NULL,"nr"=NULL,"min"=NULL,"max"=NULL,"qual"=NULL,"spikes"=0)
int.factor=1*10^5
factor=FALSE
resol=140000
sampname=names(averaged)
plot=F
thresh=2000
scale = 1
peaks = searchMZRange(df[tryRange,],values,int.factor,scale,resol,outdir,sampname,scanmode,FALSE,100,100,thresh)  

# plot
validPeaks = which(peaks$qual > 0)
where = peaks$mean[validPeaks]
plot(averaged$`1-1B`@mass[tryRange], averaged$`1-1B`@intensity[tryRange])

MassSpecWavelet::plotPeak(df[1:5000], peakIndex=msw_peaks$majorPeakInfo$allPeakIndex,mz = names(df[tryRange,]))

MassSpecWavelet::plotPeak(df[tryRange], peakIndex=msw_peaks$majorPeakInfo$allPeakIndex,mz = names(df[tryRange,]),method = "I")

points(averaged$`1-1B`@mass[tryRange], averaged$`1-1B`@intensity[tryRange])


abline(v = where, untf = FALSE, col="red")
