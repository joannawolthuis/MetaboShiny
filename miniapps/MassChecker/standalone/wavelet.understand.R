library(MassSpecWavelet)

outdir = "/Users/jwolthuis/Documents/umc/data/Data/BrSpIt/MZXML"
resdir = file.path(outdir, "results")

avg_files_pos <- list.files(file.path(resdir, "averaged"), pattern = "\\_pos\\.RData", full.names = T)
avg_files_neg <- list.files(file.path(resdir, "averaged"), pattern = "\\_neg\\.RData", full.names = T)

curr_pos <- avg_files_pos[[1]]
curr_neg <- avg_files_neg[[1]]

### SELECT CURRENT MODE FOR TESTING ###

current <- curr_pos
ppm = 2

# --- peak call these files ---

load(current)

vec <- averaged[[1]]@intensity
names(vec) <- averaged[[1]]@mass

zeros = which(vec == 0)

for(i in 1:length(zeros)){
  curr = vec[i]
  before = vec[i-1]
  after = vec[i + 1]
  if(before > 0 && after > 0) vec[i] = (before + after)/2
}

averaged[[1]]@intensity <- vec

testF <- list.files(file.path(resdir, "spectra"),full.names = T)[[1]]
load(testF)
# == SMOOTH ==

windowRange <- 2:20

for(i in windowRange){
  print(paste("---", i, "---"))
  pbefore <- length(which(pastecs::turnpoints(averaged[[1]]@intensity)$peaks))
  smooth <- MALDIquant::smoothIntensity(averaged, halfWindowSize = i)
  pafter <- length(which(pastecs::turnpoints(smooth[[1]]@intensity)$peaks))
  print(pafter - pbefore)
}

pbefore <- length(which(pastecs::turnpoints(averaged[[1]]@intensity)$peaks))
smooth <- MALDIquant::smoothIntensity(averaged,halfWindowSize=3)
pafter <- length(which(pastecs::turnpoints(smooth[[1]]@intensity)$peaks))

### WAVELET ###

# --------------------

smooth <- MALDIquant::smoothIntensity(averaged,
                                      halfWindowSize=3)

# --- another exoctic spectrum ---

plot(smooth[[1]]@intensity[150:200])


df <- matrix(data = smooth[[1]]@intensity,
             ncol=1)
rownames(df) <- smooth[[1]]@mass

df <- matrix(data = averaged[[1]]@intensity,
             ncol=1)
rownames(df) <- averaged[[1]]@mass

# --- check spec ---

vec <- MassSpecWavelet:::smoothDWT(df, nLevel = 4, wf = "mexh")


plot(df[80:130])
lines(df[80:130])

plot(vec[80:130])
lines(vec[80:130])

df <- matrix(data = vec, ncol = 1)

data <- rbind(head(df,n = 1000), tail(df, n = 1000))[,1]
data <- df[,1]

tp <- pastecs::turnpoints(data)

# ------------------------------
p = res_scale_tune$par

peaksTP = switch(tp$firstispeak,
                 "TRUE" = tp$tppos[seq(1, length(tp$tppos), 2)],
                 "FALSE" = tp$tppos[seq(2, length(tp$tppos), 2)])
scales = c(1,64)

res_scale_tune <- optim(par = c(0.2,0,2,2,6), fn = function(p){
  
  # fnp = Inf
  # fpp = Inf
  # fpa = Inf
  
  res <- c(fnp, fpp, fpa)
  
  p = c(0.5, 300, 0, 3, 5)
  
  try({
    
    jump = p[1]
    amp = p[2]
    snr = p[3]
    ws = p[4]
    max_scale = p[5]

    wCoefs <- MassSpecWavelet::cwt(data, scales = seq(1, 
                                                      max_scale, 
                                                      jump) , 
                                   wavelet='mexh')
    
    localMax <- getLocalMaximumCWT_J(wCoefs 
                                     ,amp.Th = amp,
                                     minWinSize = ws
                                     ) # need some rounding to deal w/ the decimal scales...
    
    ridgeList <- getRidge_J(localMax)
    
    majorPeakInfo <- identifyMajorPeaks(data, 
                                        ridgeList, 
                                        wCoefs, 
                                        SNR.Th = snr,
                                        peakScaleRange = 1,
                                        nearbyPeak = TRUE)
    
    peaksMSW <- majorPeakInfo$peakCenterIndex
    
    # --- explore ---
    jumps = seq(1, length(data), 100)
    
    for(i in 2:length(jumps)){
      a = jumps[i-1]
      b = jumps[i]
      MassSpecWavelet::plotPeak(data,
                                peaksMSW,range = c(a,b))
      Sys.sleep(1)
    }
    # ----------------
    
    tune <- tuneInPeakInfo_J(df,majorPeakInfo,peakScale = majorPeakInfo$peakScale,maxScale = 32)
    peaksMSW <- tune$peakCenterIndex
    
    tpp = 0
    fnp = 0
    fpp = 0
    
    window = 1
    
    for(peak in peaksTP){
      match = which(peaksMSW %between% c(peak - window, peak + window))
      #print(match)
      if(length(match > 0) == 0) fnp <- fnp + 1 else if(length(match > 0) > 1){fpp = fpp + length(match > 0) - 1; tpp = tpp + 1;}else tpp = tpp + 1 
    }
    
    fps <- peaksMSW[which(data[peaksMSW] == 0)]
    length(fps)
    
    print(paste("TRUE POS:", tpp))
    # --- return ---
    
    MassSpecWavelet::plotPeak(data,
                              peaksMSW,range = c(300,320))
    
    missing = length(peaksTP) - tpp
    
    res <- c(fnp, fpp, fpa)
    
  })
  
  summed_res <- sum(sapply(res, abs)) + (64 - max_scale)

  print(summed_res)
  
  summed_res
  
},method = "SANN", 
gr = function(x) c(sample(x = seq(0.1, 1.0, 0.05), 
                          size = 1), 
                   sample(x = seq(0, 3000, 500), 
                          size = 1),
                   sample(x = seq(0, 10, 1), 
                          size = 1),
                   sample(x = seq(2, 5, 1), 
                          size = 1),
                   sample(x = seq(2,64, 1),
                          size = 1)),
control = list(maxit = 1000))

res_scale_tune
p = res_scale_tune$par
# --- CRS ---

library(nloptr)


res_scale_tune <- nloptr::crs2lm(x0 = c(0.5,1000,0,5), fn = function(p){
  
  fnp = Inf
  fpp = Inf
  fpa = Inf
  #missing = Inf
  
  res <- c(fnp, fpp, fpa)
  
  print(p)
  
  try({
    
    jump = p[1]
    amp = p[2]
    snr = p[3]
    ws = p[4]
    
    wCoefs <- MassSpecWavelet::cwt(data, scales = seq(scales[1], scales[2], jump) , wavelet='mexh')
    
    localMax <- MassSpecWavelet::getLocalMaximumCWT(wCoefs, amp.Th = amp,minWinSize = ws)
    
    ridgeList <- getRidge(localMax,minWinSize = ws)
    
    majorPeakInfo <- identifyMajorPeaks(data, 
                                        ridgeList, 
                                        wCoefs, 
                                        SNR.Th = snr,
                                        peakScaleRange = 1,
                                        nearbyPeak = TRUE)
    peaksMSW <- majorPeakInfo$peakCenterIndex
    
    
    tpp = 0
    fnp = 0
    fpp = 0
    
    window = 1
    
    for(peak in peaksTP){
      match = which(peaksMSW %between% c(peak - window, peak + window))
      #print(match)
      if(length(match > 0) == 0) fnp <- fnp + 1 else if(length(match > 0) > 1){fpp = fpp + length(match > 0) - 1; tpp = tpp + 1;}else tpp = tpp + 1 
    }
    
    print(paste("TRUE POS:", tpp))
    # --- return ---
    
    MassSpecWavelet::plotPeak(data,
                              peaksMSW,range = c(100,200))
    
    missing = length(peaksTP) - tpp
    
    res <- c(fnp, fpp, fpa)
    
  })
  
  summed_res <- sum(sapply(res, abs))
  
  summed_res
  
}, lower = c(0.1, 0, 0, 2),
upper = c(1.0, 3000, 10, 5),
xtol_rel = 1)




peaks_wav = MALDIquant::createMassPeaks(mass = as.numeric(rownames(df)[majorPeakInfo$peakCenterIndex]),
                                        intensity = majorPeakInfo$peakValue,
                                        snr = majorPeakInfo$peakSNR)

