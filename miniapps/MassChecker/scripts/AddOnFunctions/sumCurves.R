sumCurves <- function(mean1, mean2, scale1, scale2, sigma1, sigma2, x2, x, resol, plot) {
#   mean1=mean[i-1]
#   mean2=mean[i]
#   scale1=scale[i-1]
#   scale2=scale[i]
#   sigma1=sigma[i-1]
#   sigma2=sigma[i]

  # message("=============================================================> sum 2 curves!")
  
  sumFit=(scale1*dnorm(x2,mean1,sigma1))+(scale2*dnorm(x2,mean2,sigma2))
  if (plot) lines(x2,sumFit,col="brown")
  
  #mean1plus2 = mean(c(mean1,mean2))
  mean1plus2 = weighted.mean(c(mean1,mean2),c(max(scale1*dnorm(x2,mean1,sigma1)),max(scale2*dnorm(x2,mean2,sigma2))))
  
  if (plot) abline(v = mean1plus2, col="brown")
  fwhm = getFwhm(mean1plus2, resol)
  half_max = max(sumFit)*0.5
  if (plot) lines(c(mean1plus2 - 0.5*fwhm, mean1plus2 + 0.5*fwhm),c(half_max,half_max),col="orange")
  
  # sumFit=(scale1*dnorm(x,mean1,sigma1))+(scale2*dnorm(x,mean2,sigma2))
#  fq=abs(sum(y) - sum(sumFit))/sum(y)
  
#   h2=c(paste("mean =", mean1plus2, sep=" "),
#        paste("fq =", fq, sep=" "))
#   
#   legend("topright", legend=h2)
  
  # I assume that the sum of the distributions if also normal, which is not
  #area = sum(scale1*dnorm(x2,mean1,sigma1))+sum(scale2*dnorm(x2,mean2,sigma2))
  #area = max(scale1*dnorm(x2,mean1,sigma1))+max(scale2*dnorm(x2,mean2,sigma2))
  area = max(sumFit) 
  scale = scale1 + scale2
  sigma = (fwhm/2)*0.85
   
  return(list("mean"= mean1plus2,"area"=area, "scale"=scale, "sigma"=sigma)) # "qual"=fq
  
} 
