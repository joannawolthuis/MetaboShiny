getArea <- function(mu,resol,scale,sigma,int.factor){

#   mu=p3[1]
#   scale=p3[2]
#   sigma=sigma1
  
#   mu=p3[1]
#   scale=p3[2]
#   sigma=sigma1

  # avoid too big vector (cannot allocate vector of size ...)
  if (mu>1000) return(0)
  
  fwhm = getFwhm(mu,resol)
  mzMin = mu - 2*fwhm
  mzMax = mu + 2*fwhm
  mz.range = mzMax - mzMin 
  x_int = seq(mzMin,mzMax,length=mz.range*int.factor)
  
  #plot(x_int,scale*dnorm(x_int,mu,sigma),col="red",type="l")

  return(sum(scale*dnorm(x_int,mu,sigma))/100)  
}
