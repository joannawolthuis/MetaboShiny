getSD <- function(x,y,resol=140000) {   
  
  index = which(y==max(y))
  mean = x[index]
  resol.mz = resol*(1/sqrt(2)^(log2(mean/200)))
  fwhm = mean/resol.mz
  sd = (fwhm/2)*0.85 
  return(sd)
}
