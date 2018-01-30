fitGaussian <- function(x2,x,y,index,scale,resol,outdir,force,useBounds,plot,scanmode,int.factor,width,height) {
  # force=length(index)
  # useBounds=FALSE
   
  peak.mean = NULL
  peak.area = NULL
  peak.qual = NULL
  peak.min = NULL
  peak.max = NULL
  
  FQ1 = 0.15
  FQ = 0.2
  
  # One local max
  if (force==1){

    retVal = fit1Peak(x2,x,y,index,scale,resol,plot,FQ1,useBounds)

    scale = 2
    
    if (retVal$mean[1]<x[1] | retVal$mean[1]>x[length(x)]) {   # <=== mean outside range

      # do it again
      return(fitGaussian(x2,x,y,index,scale,resol,outdir,force=1,useBounds=TRUE,plot,scanmode,int.factor,width,height))
        
    } else { # <=== mean within range
    
      if (retVal$qual > FQ1){ # <=== bad fit
          
        # Try to fit two curves
          
        # diff(diff(x)) essentially computes the discrete analogue of the second derivative, so should be negative at local maxima.
        # The +1 below takes care of the fact that the result of diff is shorter than the input vector.
        new_index=which(diff(sign(diff(y)))==-2)+1
          
        if (length(new_index)!=2) {
          new_index = round(length(x)/3)
          new_index = c(new_index, 2*new_index)
        }
        
        #length(new_index)  
        return(fitGaussian(x2,x,y,new_index,scale,resol,outdir,force=2,useBounds=FALSE,plot,scanmode,int.factor,width,height))
          
      } else { # <=== good fit
          
        peak.mean = c(peak.mean, retVal$mean)
        #peak.area = c(peak.area, sum(retVal$scale*dnorm(x2,retVal$mean,retVal$sigma)))
        # "Centroid mode"
#         peak.area = c(peak.area, max(retVal$scale*dnorm(x2,retVal$mean,retVal$sigma)))
        peak.area = c(peak.area, getArea(retVal$mean,resol,retVal$scale,retVal$sigma,int.factor))
        peak.qual = retVal$qual
        peak.min = x[1]
        peak.max = x[length(x)]
      }
    }
    
  # Two local max  
  } else if (force==2 & (length(x)>6)) {
    
    # fit two curves
    retVal = fit2peaks(x2,x,y,index,scale,resol,useBounds,plot,FQ,int.factor) # <=== fit 2 curves

    if (retVal$mean[1]<x[1] | retVal$mean[1]>x[length(x)] |   # <=== one of means outside range
        retVal$mean[2]<x[1] | retVal$mean[2]>x[length(x)]) {
      
      # check quality
      if (retVal$qual > FQ) { # <=== bad fit 
        
        # do it again
        return(fitGaussian(x2,x,y,index,scale,resol,outdir,force=2,useBounds=TRUE,plot,scanmode,int.factor,width,height))
        
        # good fit
      } else {
        
        # check which mean is outside range
        # Todo ========> check this!!!!!!!!!!!!!!!
        for (i in 1:length(retVal$mean)){
          if (retVal$mean[i]<x[1] | retVal$mean[i]>x[length(x)] ) {
            peak.mean = c(peak.mean, -i)
            peak.area = c(peak.area, -i)
          } else {
            peak.mean = c(peak.mean, retVal$mean[i])
            peak.area = c(peak.area, retVal$area[i])
          }
        }
        peak.qual = retVal$qual
        peak.min = x[1]
        peak.max = x[length(x)]
      }
    } else { # <=== all means within range
          
      if (retVal$qual > FQ) { # <=== bad fit 
              
        # Try to fit three curves
        new_index=which(diff(sign(diff(y)))==-2)+1
              
        if (length(new_index)!=3) {
          new_index = round(length(x)/4)
          new_index = c(new_index, 2*new_index, 3*new_index)
        }
        
        #length(new_index)      
        return(fitGaussian(x2,x,y,new_index,scale,resol,outdir,force=3,useBounds=FALSE,plot,scanmode,int.factor,width,height))
            
      } else { # <======== good fit
        
        # check if means are within 3 ppm and sum if so  
        tmp = retVal$qual 
          
        nMeanNew = -1
        nMean = length(retVal$mean)
        while (nMean!=nMeanNew){
          nMean = length(retVal$mean)
          retVal = isWithinXppm(retVal$mean, retVal$scale, retVal$sigma, retVal$area, x2, x, ppm=4, resol, plot)
          nMeanNew = length(retVal$mean)
        }
          
        retVal$qual = tmp
        
        h2=NULL

        for (i in 1:length(retVal$mean)){
          h2 = c(h2, paste("mean =", retVal$mean[i], sep=" "))
            
          peak.mean = c(peak.mean, retVal$mean[i])
          peak.area = c(peak.area, retVal$area[i])
        }
        peak.qual = retVal$qual
        peak.min = x[1]
        peak.max = x[length(x)]
        
        h2 = c(h2, paste("fq =", retVal$qual, sep=" "))
        if (plot) legend("topright", legend=h2)
      }  
    }

  # Three local max  
  } else if (force==3 & (length(x)>6)){
    
    retVal = fit3peaks(x2,x,y,index,scale,resol,useBounds,plot,FQ,int.factor)

    # outside range
    if (retVal$mean[1]<x[1] | retVal$mean[1]>x[length(x)] | # <=== one of means outside range
        retVal$mean[2]<x[1] | retVal$mean[2]>x[length(x)] |
        retVal$mean[3]<x[1] | retVal$mean[3]>x[length(x)]) {
      
      # check quality
      if (retVal$qual > FQ) { # <=== bad fit 
        
        # do it again
        return(fitGaussian(x2,x,y,index,scale,resol,outdir,force,useBounds=TRUE,plot,scanmode,int.factor,width,height))
      
      # good fit
      } else {
        
        # check which mean is outside range
        # Todo ========> check this!!!!!!!!!!!!!!!
        for (i in 1:length(retVal$mean)){
          if (retVal$mean[i]<x[1] | retVal$mean[i]>x[length(x)] ) {
            peak.mean = c(peak.mean, -i)
            peak.area = c(peak.area, -i)
          } else {
            peak.mean = c(peak.mean, retVal$mean[i])
            peak.area = c(peak.area, retVal$area[i])
          }
        }
        peak.qual = retVal$qual
        peak.min = x[1]
        peak.max = x[length(x)]
      }
      
    } else { # <=== all means within range
      
      if (retVal$qual > FQ) { # <=== bad fit 
        
        # Try to fit four curves
        new_index=which(diff(sign(diff(y)))==-2)+1
        
        if (length(new_index)!=4) {
          new_index = round(length(x)/5)
          new_index = c(new_index, 2*new_index, 3*new_index, 4*new_index)
        }

        #length(new_index)
        return(fitGaussian(x2,x,y,new_index,scale,resol,outdir,force=4,useBounds=FALSE,plot,scanmode,int.factor,width,height))
      
      } else { # <======== good fit
        
        # check if means are within 3 ppm and sum if so  
        tmp = retVal$qual 
        
        nMeanNew = -1
        nMean = length(retVal$mean)
        while (nMean!=nMeanNew){
          nMean = length(retVal$mean)
          retVal = isWithinXppm(retVal$mean, retVal$scale, retVal$sigma, retVal$area, x2, x, ppm=4, resol, plot)
          nMeanNew = length(retVal$mean)
        }
        
        retVal$qual = tmp

        h2=NULL
#         peak.mean=NULL
#         peak.area=NULL
        
        for (i in 1:length(retVal$mean)){
          h2 = c(h2, paste("mean =", retVal$mean[i], sep=" "))
            
          peak.mean = c(peak.mean, retVal$mean[i])
          peak.area = c(peak.area, retVal$area[i])
        }
        peak.qual = retVal$qual
        peak.min = x[1]
        peak.max = x[length(x)]
        
        h2 = c(h2, paste("fq =", retVal$qual, sep=" "))
        if (plot) legend("topright", legend=h2)
      
        }  
    }

  # Four local max  
  } else if (force==4 & (length(x)>6)){
    
    retVal = fit4peaks(x2,x,y,index,scale,resol,useBounds,plot,FQ,int.factor)

    if (retVal$mean[1]<x[1] | retVal$mean[1]>x[length(x)] |
        retVal$mean[2]<x[1] | retVal$mean[2]>x[length(x)] |
        retVal$mean[3]<x[1] | retVal$mean[3]>x[length(x)] |
        retVal$mean[4]<x[1] | retVal$mean[4]>x[length(x)]) {
      
      # check quality
      if (retVal$qual > FQ) { # <=== bad fit 
        
        # do it again
        return(fitGaussian(x2,x,y,index,scale,resol,outdir,force,useBounds=TRUE,plot,scanmode,int.factor,width,height))

      # good fit
      } else {
        
        # check which mean is outside range
        # Todo ========> check this!!!!!!!!!!!!!!!
        for (i in 1:length(retVal$mean)){
          if (retVal$mean[i]<x[1] | retVal$mean[i]>x[length(x)] ) {
            peak.mean = c(peak.mean, -i)
            peak.area = c(peak.area, -i)
          } else {
            peak.mean = c(peak.mean, retVal$mean[i])
            peak.area = c(peak.area, retVal$area[i])
          }
        }
        peak.qual = retVal$qual
        peak.min = x[1]
        peak.max = x[length(x)]
      }
      
    } else {
      
      if (retVal$qual > FQ) { # <=== bad fit 
        
        # Generate 1 curve
        return(fitGaussian(x2,x,y,index,scale,resol,outdir,force=5,useBounds=FALSE,plot,scanmode,int.factor,width,height))
      
      } else { # <======== good fit

        # check if means are within 3 ppm and sum if so  
        tmp = retVal$qual 
        
        nMeanNew = -1
        nMean = length(retVal$mean)
        while (nMean!=nMeanNew){
          nMean = length(retVal$mean)
          retVal = isWithinXppm(retVal$mean, retVal$scale, retVal$sigma, retVal$area, x2, x, ppm=4, resol, plot)
          nMeanNew = length(retVal$mean)
        }
        
        retVal$qual = tmp
        
        h2=NULL
        
        for (i in 1:length(retVal$mean)){
          h2 = c(h2, paste("mean =", retVal$mean[i], sep=" "))
            
          peak.mean = c(peak.mean, retVal$mean[i])
          peak.area = c(peak.area, retVal$area[i])
        }
        peak.qual = retVal$qual
        peak.min = x[1]
        peak.max = x[length(x)]
        
        h2 = c(h2, paste("fq =", retVal$qual, sep=" "))
        if (plot) legend("topright", legend=h2)
      }  
    }

  # More then four local max  
  } else {
    
    scale=2
    FQ1=0.40
    useBounds=TRUE
    index=which(y==max(y))
    retVal = fit1Peak(x2,x,y,index,scale,resol,plot,FQ1,useBounds)
    
    if (retVal$qual > FQ1){ # <=== bad fit
      
      if (plot) dev.off()    
      
      rval = generateGaussian(x,y,resol,plot,scanmode,int.factor, width, height)
      peak.mean = c(peak.mean, rval$mean)
      peak.area = c(peak.area, rval$area)
      peak.min = rval$min
      peak.max = rval$max
      peak.qual = 0
      
    } else { # <=== good fit
      
      peak.mean = c(peak.mean, retVal$mean)
      peak.area = c(peak.area, getArea(retVal$mean,resol,retVal$scale,retVal$sigma,int.factor))
      peak.qual = retVal$qual
      peak.min = x[1]
      peak.max = x[length(x)]
    }
  }
  
  return(list("mean"=peak.mean, "area"=peak.area, "qual" = peak.qual, "min"=peak.min, "max"=peak.max))
}
