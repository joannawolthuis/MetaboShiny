isWithinXppm <- function(mean, scale, sigma, area, x2, x, ppm=4, resol, plot) {
#   mean=retVal$mean
#   scale=retVal$scale
#   sigma=retVal$sigma
#   area=retVal$area
#   ppm=3

  # sort!!!!!!!!!!!!!!!!
  index = order(mean)
  mean = mean[index]
  scale = scale[index]
  sigma = sigma[index]
  area = area[index]
  
  summed = NULL
  remove = NULL
  
  if (length(mean)>1){
    for (i in 2:length(mean)){
      if ((abs(mean[i-1]-mean[i])/mean[i-1])*10^6 < ppm) {
        
        # avoid double occurance in sum
        if ((i-1) %in% summed) next
        
        retVal = sumCurves(mean[i-1], mean[i], scale[i-1], scale[i], sigma[i-1], sigma[i], x2, x, resol, plot)
        summed = c(summed, i-1, i)
        if (is.nan(retVal$mean)) retVal$mean=0 
        mean[i-1] = retVal$mean
        mean[i] = retVal$mean
        area[i-1] = retVal$area
        area[i] = retVal$area
        scale[i-1] = retVal$scale
        scale[i] = retVal$scale
        sigma[i-1] = retVal$sigma
        sigma[i] = retVal$sigma
        
        remove = c(remove, i)
      }
    }
  }
  
  if (length(remove)!=0){
    mean=mean[-c(remove)]
    area=area[-c(remove)]
    scale=scale[-c(remove)]
    sigma=sigma[-c(remove)]
  }
  
  return(list("mean"=mean, "area"=area, "scale"=scale, "sigma"=sigma, "qual"=NULL))
  
}