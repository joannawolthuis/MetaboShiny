getDeltaMZ <- function(mass, breaks.fwhm){
  index = which(breaks.fwhm<mass)
  tmp = index[length(index)]
  return(breaks.fwhm[tmp+1] - breaks.fwhm[tmp]) 
}
