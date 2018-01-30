fit1Peak <- function(x2,x,y,index,scale,resol,plot,FQ,useBounds) {
#FQ=FQ1
  
  if (length(y)<3){
    message("Range to small, no fit possible!")
  } else {
    
    if ((length(y)==4)) {
      mu = weighted.mean(x,y)
      sigma = getSD(x,y)
      fitP = fitG_2(x,y,sigma,mu,scale,useBounds)
    } else {
      
      if ((length(x) - length(index)) < 2) {
        range1=c((length(x)-4):length(x))
      } else if (length(index) < 2) {
        range1=c(1:5)
      } else {
        range1=c(index[1]-2,index[1]-1,index[1],index[1]+1,index[1]+2)
      }
      
      if (range1[1]==0) range1=range1[-1]
      
      # remove NA
      if (length(which(is.na(y[range1])))!=0) range1=range1[-which(is.na(y[range1]))]

      mu = weighted.mean(x[range1],y[range1])
      sigma = getSD(x[range1],y[range1])
      fitP = fitG_2(x,y,sigma,mu,scale,useBounds)
    }  
    
    p2 = fitP$par
    
    #fq_new = abs(sum(y) - sum(p2[2]*dnorm(x,p2[1],sigma)))/sum(y)
    fq_new = getFitQuality(x,y,p2[1],p2[1],resol,p2[2],sigma)$fq_new
    
    if (plot & (fq_new < FQ)) lines(x2,p2[2]*dnorm(x2,p2[1],sigma), col="green")
    
    scale_new = 1.2*scale
    # message(fq_new)
    
    if (fq_new > FQ) { # <=== bad fit?
      # optimize scaling factor
      fq = 0
      scale = 0
      
      if (sum(y)>sum(p2[2]*dnorm(x,p2[1],sigma))){
        while ((round(fq, digits = 3) != round(fq_new, digits = 3)) & (scale_new<10000))  { 
          fq = fq_new
          scale = scale_new
          
          #message(scale)
          fitP = fitG_2(x,y,sigma,mu,scale,useBounds)
          p2 = fitP$par
          
          #fq_new = abs(sum(y) - sum(p2[2]*dnorm(x,p2[1],sigma)))/sum(y)
          fq_new = getFitQuality(x,y,p2[1],p2[1],resol,p2[2],sigma)$fq_new
          scale_new=1.2*scale
          
          if (plot & (fq_new < FQ)) lines(x2,p2[2]*dnorm(x2,p2[1],sigma), col="green")
#           message(paste("fq_new: ", fq_new))
#           message(paste("scale_new: ", scale_new))
          #           (round(fq, digits = 4) != round(fq_new, digits = 4))
        }
      } else {
        while ((round(fq, digits = 3) != round(fq_new, digits = 3)) & (scale_new<10000)) { 
          fq = fq_new
          scale = scale_new
          
          # message(scale)
          fitP = fitG_2(x,y,sigma,mu,scale,useBounds)
          p2 = fitP$par
          
          #fq_new = abs(sum(y) - sum(p2[2]*dnorm(x,p2[1],sigma)))/sum(y)
          fq_new = getFitQuality(x,y,p2[1],p2[1],resol,p2[2],sigma)$fq_new
          scale_new=0.8*scale
          
          if (plot & (fq_new < FQ)) lines(x2,p2[2]*dnorm(x2,p2[1],sigma), col="green")
          # message(paste("fq_new: ", fq_new))
        }
      }
      
      if (fq < fq_new) {
#         message(paste("fq_new: ", fq_new))
#         message(paste("fq: ", fq))
#         message(paste("scale_new: ", scale_new))
#         message(paste("scale: ", scale))
        
        fitP = fitG_2(x,y,sigma,mu,scale,useBounds)
        p2 = fitP$par
        fq_new = fq
        # message(paste("==> fq_new: ", fq_new))
        if (plot & (fq_new < FQ)) lines(x2,p2[2]*dnorm(x2,p2[1],sigma), col="dark green")
        
      }
    }

    if (plot & (fq_new < FQ)) {
      # plot ###################
      #lines(x2,p2[2]*dnorm(x2,p2[1],sigma), col="green")
      fwhm = getFwhm(p2[1],resol)
      half_max = max(p2[2]*dnorm(x2,p2[1],sigma))*0.5
      lines(c(p2[1] - 0.5*fwhm, p2[1] + 0.5*fwhm),c(half_max,half_max),col="orange")
      abline(v = p2[1], col="green")
      h=c(paste("mean =", p2[1], sep=" "),
          paste("fq =", fq_new, sep=" "))
      legend("topright", legend=h)
      ##########################
      
      
#       abline(v = x[6], col="red")
#       fwhm = getFwhm(x[6])
#       abline(v =x[6] + 0.6*fwhm, col="red")
#       abline(v =x[6] - 0.6*fwhm, col="red")
      
    }
  }
  
  return(list("mean"=p2[1], "scale"=p2[2], "sigma"=sigma, "qual"=fq_new))
}  