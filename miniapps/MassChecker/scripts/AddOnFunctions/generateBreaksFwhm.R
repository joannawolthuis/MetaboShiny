# options(digits=16)
# resol = 140000
# breaks.fwhm = 1 
# i = 1
# breaks.fwhm.avg = NULL
# 
# Sys.time()
# 
# while (breaks.fwhm[length(breaks.fwhm)]<1000){
# 
#   resol.mz = resol*(1/sqrt(2)^(log2(breaks.fwhm[i]/200)))
#   fwhm.0.1 = (breaks.fwhm[i]/resol.mz)/10
#   breaks.fwhm = c(breaks.fwhm, breaks.fwhm[i] + fwhm.0.1) 
#   breaks.fwhm.avg = c(breaks.fwhm.avg,(breaks.fwhm[i] +breaks.fwhm[i+1])/2)
#   
#   if (i %% 10000 == 0){ 
#     message(paste("i =", i)) 
#     message(paste("breaks.fwhm =", breaks.fwhm[length(breaks.fwhm)]))
#   }
#   
#   i = i + 1
# }
# 
# Sys.time()
# save(list=c("breaks.fwhm", "breaks.fwhm.avg"), file="breaks.RData")
# message("Breaks saved!")
