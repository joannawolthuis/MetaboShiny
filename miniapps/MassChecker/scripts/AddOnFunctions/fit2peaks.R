fit2peaks <- function(x2,x,y,index,scale,resol,useBounds=FALSE,plot=FALSE,FQ,int.factor){
  
  peak.mean = NULL
  peak.area = NULL
  peak.scale = NULL
  peak.sigma = NULL

  range1=c(index[1]-2,index[1]-1,index[1],index[1]+1,index[1]+2)
  if (range1[1]==0) range1=range1[-1]

  range2=c(index[2]-2,index[2]-1,index[2],index[2]+1,index[2]+2)
    
  if (length(x)<range2[length(range2)]) range2=range2[-length(range2)]
    
  range1=checkOverlap(range1,range2)[[1]]
  range2=checkOverlap(range1,range2)[[2]]
  
  # check for negative or 0
  remove=which(range1<1)
  if (length(remove)>0) range1=range1[-remove]
  remove=which(range2<1)
  if (length(remove)>0) range2=range2[-remove]
  
  # remove NA
  if (length(which(is.na(y[range1])))!=0) range1=range1[-which(is.na(y[range1]))]
  if (length(which(is.na(y[range2])))!=0) range2=range2[-which(is.na(y[range2]))]
  
  mu1 = weighted.mean(x[range1],y[range1])
  sigma1 = getSD(x[range1],y[range1])
  
#   message(paste("fit2peaks mu =>", mu1))
#   message(paste("fit2peaks sigma =>", sigma1))
#   message(paste("fit2peaks scale =>", scale))
  
  fitP = fitG_2(x[range1],y[range1],sigma1,mu1,scale,useBounds)
  p = fitP$par
  
  mu2 = weighted.mean(x[range2],y[range2])
  sigma2 = getSD(x[range2],y[range2])
  fitP = fitG_2(x[range2],y[range2],sigma2,mu2,scale,useBounds)
  p2 = fitP$par
    
  fit2P = fit2G_2(x, y, sigma1, sigma2, p[1], p[2], p2[1], p2[2],useBounds)
  p3 = fit2P$par

  if (is.null(sigma2)) sigma2=sigma1
  
  
  # plot ###################
  sumFit2 = (p3[2]*dnorm(x2,p3[1],sigma1))+(p3[4]*dnorm(x2,p3[3],sigma2))
  sumFit = (p3[2]*dnorm(x,p3[1],sigma1))+(p3[4]*dnorm(x,p3[3],sigma2))
  fq=getFitQuality(x,y,sort(c(p3[1],p3[3]))[1],sort(c(p3[1],p3[3]))[2],resol,sumFit=sumFit)$fq_new

  fwhm = getFwhm(p3[1],resol)
  half_max = max(p3[2]*dnorm(x2,p3[1],sigma1))*0.5
  if (plot & (fq < FQ)) lines(c(p3[1] - 0.5*fwhm, p3[1] + 0.5*fwhm),c(half_max,half_max),col="orange")
  if (plot & (fq < FQ)) lines(x2,p3[4]*dnorm(x2,p3[3],sigma2),col="grey")
  if (plot & (fq < FQ)) abline(v = p3[3], col="grey")

  fwhm = getFwhm(p3[3],resol)
  half_max = max(p3[4]*dnorm(x2,p3[3],sigma2))*0.5
  if (plot & (fq < FQ)) lines(c(p3[3] - 0.5*fwhm, p3[3] + 0.5*fwhm),c(half_max,half_max),col="orange")

  if (plot & (fq < FQ)) lines(x2,sumFit2,col="black")

  if (plot & (fq < FQ)) lines(x2,p3[2]*dnorm(x2,p3[1],sigma1),col="grey")
  if (plot & (fq < FQ)) abline(v = p3[1], col="grey")
  
    
  h2=c(paste("mean =", p3[1], sep=" "),
       paste("mean =", p3[3], sep=" "),
       paste("fq =", fq, sep=" "))
    
  if (plot & (fq < FQ)) legend("topright", legend=h2)
  ##########################

#  lines(x2,p3[4]*dnorm(x2,p3[3],sigma2),col="red")
#   area1 = sum(p3[2]*dnorm(x2,p3[1],sigma1))
#   area2 = sum(p3[4]*dnorm(x2,p3[3],sigma2))
  
#   area1 = max(p3[2]*dnorm(x2,p3[1],sigma1))
#   area2 = max(p3[4]*dnorm(x2,p3[3],sigma2))

  area1 = getArea(p3[1],resol,p3[2],sigma1,int.factor)
  area2 = getArea(p3[3],resol,p3[4],sigma2,int.factor)

  peak.area = c(peak.area, area1)
  peak.area = c(peak.area, area2)

  peak.mean = c(peak.mean, p3[1])
  peak.mean = c(peak.mean, p3[3])
    
  peak.scale = c(peak.scale, p3[2])
  peak.scale = c(peak.scale, p3[4])

  peak.sigma = c(peak.sigma, sigma1)
  peak.sigma = c(peak.sigma, sigma2)
    
  return(list("mean"=peak.mean, "scale"=peak.scale, "sigma"=peak.sigma, "area"=peak.area, "qual"=fq))
}