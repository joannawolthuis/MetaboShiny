fit3peaks <- function(x2,x,y,index,scale,resol,useBounds=FALSE,plot=FALSE,FQ,int.factor){
  
  peak.mean = NULL
  peak.area = NULL
  peak.scale = NULL
  peak.sigma = NULL

  range1=c(index[1]-2,index[1]-1,index[1],index[1]+1,index[1]+2)
  range2=c(index[2]-2,index[2]-1,index[2],index[2]+1,index[2]+2)
  range3=c(index[3]-2,index[3]-1,index[3],index[3]+1,index[3]+2)
  
  remove=which(range1<1)
  if (length(remove)>0) {
    range1=range1[-remove]
  }
  remove=which(range2<1)
  if (length(remove)>0) {
    range2=range2[-remove]
  }

  if (length(x)<range3[length(range3)]) range3=range3[-length(range3)]
  
  range1=checkOverlap(range1,range2)[[1]]
  range2=checkOverlap(range1,range2)[[2]]
  range2=checkOverlap(range2,range3)[[1]]
  range3=checkOverlap(range2,range3)[[2]]
  
  # check for negative or 0
  remove=which(range1<1)
  if (length(remove)>0) range1=range1[-remove]
  remove=which(range2<1)
  if (length(remove)>0) range2=range2[-remove]
  remove=which(range3<1)
  if (length(remove)>0) range3=range3[-remove]
  
  # remove NA
  if (length(which(is.na(y[range1])))!=0) range1=range1[-which(is.na(y[range1]))]
  if (length(which(is.na(y[range2])))!=0) range2=range2[-which(is.na(y[range2]))]
  if (length(which(is.na(y[range3])))!=0) range3=range3[-which(is.na(y[range3]))]
  
  mu1 = weighted.mean(x[range1],y[range1])
  sigma1 = getSD(x[range1],y[range1])
  fitP = fitG_2(x[range1],y[range1],sigma1,mu1,scale,useBounds)
  p = fitP$par
  
  mu2 = weighted.mean(x[range2],y[range2])
  sigma2 = getSD(x[range2],y[range2])
  fitP = fitG_2(x[range2],y[range2],sigma2,mu2,scale,useBounds)
  p2 = fitP$par
  
  mu3 = weighted.mean(x[range3],y[range3])
  sigma3 = getSD(x[range3],y[range3])
  fitP = fitG_2(x[range3],y[range3],sigma3,mu3,scale,useBounds)
  p3 = fitP$par
  
  fit3P = fit3G_2(x, y, sigma1, sigma2, sigma3, p[1], p[2], p2[1], p2[2], p3[1], p3[2], useBounds)
  p4 = fit3P$par
  
  # plot ##############################
  sumFit2 = (p4[2]*dnorm(x2,p4[1],sigma1))+(p4[4]*dnorm(x2,p4[3],sigma2))+(p4[6]*dnorm(x2,p4[5],sigma3))
  sumFit = (p4[2]*dnorm(x,p4[1],sigma1))+(p4[4]*dnorm(x,p4[3],sigma2))+(p4[6]*dnorm(x,p4[5],sigma3))
  fq=getFitQuality(x,y,sort(c(p4[1],p4[3],p4[5]))[1],sort(c(p4[1],p4[3],p4[5]))[3],resol,sumFit=sumFit)$fq_new

  if (plot & (fq < FQ)) lines(x2,p4[2]*dnorm(x2,p4[1],sigma1),col="yellow")
  if (plot & (fq < FQ)) abline(v = p4[1], col="yellow")
  fwhm = getFwhm(p4[1],resol)
  half_max = max(p4[2]*dnorm(x2,p4[1],sigma1))*0.5
  if (plot & (fq < FQ)) lines(c(p4[1] - 0.5*fwhm, p4[1] + 0.5*fwhm),c(half_max,half_max),col="orange")
  
  if (plot & (fq < FQ)) lines(x2,p4[4]*dnorm(x2,p4[3],sigma2),col="yellow")
  if (plot & (fq < FQ)) abline(v = p4[3], col="yellow")
  fwhm = getFwhm(p4[3],resol)
  half_max = max(p4[4]*dnorm(x2,p4[3],sigma2))*0.5
  if (plot & (fq < FQ)) lines(c(p4[3] - 0.5*fwhm, p4[3] + 0.5*fwhm),c(half_max,half_max),col="orange")
  
  if (plot & (fq < FQ)) lines(x2,p4[6]*dnorm(x2,p4[5],sigma3),col="yellow")
  if (plot & (fq < FQ)) abline(v = p4[5], col="yellow")
  fwhm = getFwhm(p4[5],resol)
  half_max = max(p4[6]*dnorm(x2,p4[5],sigma3))*0.5
  if (plot & (fq < FQ)) lines(c(p4[5] - 0.5*fwhm, p4[5] + 0.5*fwhm),c(half_max,half_max),col="orange")
  
  if (plot & (fq < FQ)) lines(x2,sumFit2,col="red")

  h2=c(paste("mean =", p4[1], sep=" "),
       paste("mean =", p4[3], sep=" "),
       paste("mean =", p4[5], sep=" "),
       paste("fq =", fq, sep=" "))
  
  if (plot & (fq < FQ)) legend("topright", legend=h2)
  #########################################
  
#   area1 = sum(p4[2]*dnorm(x2,p4[1],sigma1))
#   area2 = sum(p4[4]*dnorm(x2,p4[3],sigma2))
#   area3 = sum(p4[6]*dnorm(x2,p4[5],sigma3))

#   area1 = max(p4[2]*dnorm(x2,p4[1],sigma1))
#   area2 = max(p4[4]*dnorm(x2,p4[3],sigma2))
#   area3 = max(p4[6]*dnorm(x2,p4[5],sigma3))
  
  area1 = getArea(p4[1],resol,p4[2],sigma1,int.factor)
  area2 = getArea(p4[3],resol,p4[4],sigma2,int.factor)
  area3 = getArea(p4[5],resol,p4[6],sigma3,int.factor)

  peak.area = c(peak.area, area1)
  peak.area = c(peak.area, area2)
  peak.area = c(peak.area, area3)
  
  peak.mean = c(peak.mean, p4[1])
  peak.mean = c(peak.mean, p4[3])
  peak.mean = c(peak.mean, p4[5])
  
  peak.scale = c(peak.scale, p4[2])
  peak.scale = c(peak.scale, p4[4])
  peak.scale = c(peak.scale, p4[6])
  
  peak.sigma = c(peak.sigma, sigma1)
  peak.sigma = c(peak.sigma, sigma2)
  peak.sigma = c(peak.sigma, sigma3)
  
  return(list("mean"=peak.mean, "scale"=peak.scale, "sigma"=peak.sigma, "area"=peak.area, "qual"=fq))

}