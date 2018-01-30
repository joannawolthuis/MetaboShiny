fit4peaks <- function(x2,x,y,index,scale,resol,useBounds=FALSE,plot=FALSE,FQ,int.factor) {
  
  peak.mean = NULL
  peak.area = NULL
  peak.scale = NULL
  peak.sigma = NULL
  
  range1=c(index[1]-2,index[1]-1,index[1],index[1]+1,index[1]+2)
  range2=c(index[2]-2,index[2]-1,index[2],index[2]+1,index[2]+2)
  range3=c(index[3]-2,index[3]-1,index[3],index[3]+1,index[3]+2)
  range4=c(index[4]-2,index[4]-1,index[4],index[4]+1,index[4]+2)
  if (range1[1]==0) range1=range1[-1]
  if (length(x)<range4[length(range4)]) range4=range4[-length(range4)]
  
  range1=checkOverlap(range1,range2)[[1]]
  range2=checkOverlap(range1,range2)[[2]]
  range2=checkOverlap(range2,range3)[[1]]
  range3=checkOverlap(range2,range3)[[2]]
  range3=checkOverlap(range3,range4)[[1]]
  range4=checkOverlap(range3,range4)[[2]]
  
  remove=which(range4>length(x))
  if (length(remove)>0) {
    range4=range4[-remove]
    # message(length(range4))
  }
  
  # check for negative or 0
  remove=which(range1<1)
  if (length(remove)>0) range1=range1[-remove]
  remove=which(range2<1)
  if (length(remove)>0) range2=range2[-remove]
  remove=which(range3<1)
  if (length(remove)>0) range3=range3[-remove]
  remove=which(range4<1)
  if (length(remove)>0) range4=range4[-remove]
  
  # remove NA
  if (length(which(is.na(y[range1])))!=0) range1=range1[-which(is.na(y[range1]))]
  if (length(which(is.na(y[range2])))!=0) range2=range2[-which(is.na(y[range2]))]
  if (length(which(is.na(y[range3])))!=0) range3=range3[-which(is.na(y[range3]))]
  if (length(which(is.na(y[range4])))!=0) range4=range4[-which(is.na(y[range4]))]

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
  
  mu4 = weighted.mean(x[range4],y[range4])
  sigma4 = getSD(x[range4],y[range4])
  fitP = fitG_2(x[range4],y[range4],sigma4,mu4,scale,useBounds)
  p4 = fitP$par
  
  fit4P = fit4G_2(x, y, sigma1, sigma2, sigma3, sigma3, p[1], p[2], p2[1], p2[2], p3[1], p3[2],  p4[1], p4[2], useBounds)
  p5 = fit4P$par
  
  # plot #####################################
  sumFit2 = (p5[2]*dnorm(x2,p5[1],sigma1))+(p5[4]*dnorm(x2,p5[3],sigma2))+(p5[6]*dnorm(x2,p5[5],sigma3))+(p5[8]*dnorm(x2,p5[7],sigma3))
  sumFit = (p5[2]*dnorm(x,p5[1],sigma1))+(p5[4]*dnorm(x,p5[3],sigma2))+(p5[6]*dnorm(x,p5[5],sigma3))+(p5[8]*dnorm(x,p5[7],sigma3))
  fq=getFitQuality(x,y,sort(c(p5[1],p5[3],p5[5],p5[7]))[1],sort(c(p5[1],p5[3],p5[5],p5[7]))[4],resol,sumFit=sumFit)$fq_new
  
  if (plot & (fq < FQ)) lines(x2,p5[2]*dnorm(x2,p5[1],sigma1),col="purple")
  if (plot & (fq < FQ)) abline(v = p5[1], col="purple")
  fwhm = getFwhm(p5[1],resol)
  half_max = max(p5[2]*dnorm(x2,p5[1],sigma1))*0.5
  if (plot & (fq < FQ)) lines(c(p5[1] - 0.5*fwhm, p5[1] + 0.5*fwhm),c(half_max,half_max),col="orange")
  
  if (plot & (fq < FQ)) lines(x2,p5[4]*dnorm(x2,p5[3],sigma2),col="purple")
  if (plot & (fq < FQ)) abline(v = p5[3], col="purple")
  fwhm = getFwhm(p5[3],resol)
  half_max = max(p5[4]*dnorm(x2,p5[3],sigma2))*0.5
  if (plot & (fq < FQ)) lines(c(p5[3] - 0.5*fwhm, p5[3] + 0.5*fwhm),c(half_max,half_max),col="orange")
  
  if (plot & (fq < FQ)) lines(x2,p5[6]*dnorm(x2,p5[5],sigma3),col="purple")
  if (plot & (fq < FQ)) abline(v = p5[5], col="purple")
  fwhm = getFwhm(p5[5],resol)
  half_max = max(p5[6]*dnorm(x2,p5[5],sigma3))*0.5
  if (plot & (fq < FQ)) lines(c(p5[5] - 0.5*fwhm, p5[5] + 0.5*fwhm),c(half_max,half_max),col="orange")
  
  if (plot & (fq < FQ)) lines(x2,p5[8]*dnorm(x2,p5[7],sigma3),col="purple")
  if (plot & (fq < FQ)) abline(v = p5[7], col="purple")
  fwhm = getFwhm(p5[7],resol)
  half_max = max(p5[8]*dnorm(x2,p5[7],sigma3))*0.5
  if (plot & (fq < FQ)) lines(c(p5[7] - 0.5*fwhm, p5[7] + 0.5*fwhm),c(half_max,half_max),col="orange")
  
  if (plot & (fq < FQ)) lines(x2,sumFit2,col="blue")
  
  #fq = abs(sum(y) - sum(sumFit))/sum(y)
  #fq=abs(sum(y) - sum(sumFit))/sum(y)
  #fq=mean(abs(sumFit - y)/sumFit)
  
  
  h2=c(paste("mean =", p5[1], sep=" "),
       paste("mean =", p5[3], sep=" "),
       paste("mean =", p5[5], sep=" "),
       paste("mean =", p5[7], sep=" "),
       paste("fq =", fq, sep=" "))
  
  if (plot & (fq < FQ)) legend("topright", legend=h2)
  #############################################
 
#   area1 = sum(p5[2]*dnorm(x2,p5[1],sigma1))
#   area2 = sum(p5[4]*dnorm(x2,p5[3],sigma2))
#   area3 = sum(p5[6]*dnorm(x2,p5[5],sigma3))
#   area4 = sum(p5[8]*dnorm(x2,p5[7],sigma4))

#   area1 = max(p5[2]*dnorm(x2,p5[1],sigma1))
#   area2 = max(p5[4]*dnorm(x2,p5[3],sigma2))
#   area3 = max(p5[6]*dnorm(x2,p5[5],sigma3))
#   area4 = max(p5[8]*dnorm(x2,p5[7],sigma4))
  
  area1 = getArea(p5[1],resol,p5[2],sigma1,int.factor)
  area2 = getArea(p5[3],resol,p5[4],sigma2,int.factor)
  area3 = getArea(p5[5],resol,p5[6],sigma3,int.factor)
  area4 = getArea(p5[7],resol,p5[8],sigma4,int.factor)
  
  peak.area = c(peak.area, area1)
  peak.area = c(peak.area, area2)
  peak.area = c(peak.area, area3)
  peak.area = c(peak.area, area4)
  
  peak.mean = c(peak.mean, p5[1])
  peak.mean = c(peak.mean, p5[3])
  peak.mean = c(peak.mean, p5[5])
  peak.mean = c(peak.mean, p5[7])
  
  peak.scale = c(peak.scale, p5[2])
  peak.scale = c(peak.scale, p5[4])
  peak.scale = c(peak.scale, p5[6])
  peak.scale = c(peak.scale, p5[8])
  
  peak.sigma = c(peak.sigma, sigma1)
  peak.sigma = c(peak.sigma, sigma2)
  peak.sigma = c(peak.sigma, sigma3)
  peak.sigma = c(peak.sigma, sigma4)
  
  return(list("mean"=peak.mean, "scale"=peak.scale, "sigma"=peak.sigma, "area"=peak.area, "qual"=fq))
  
}