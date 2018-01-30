fitGaussianInit <- function(x,y,int.factor,scale,resol,outdir,sampname,scanmode,plot,width,height,i,start,end) {
# scanmode="negative"
  
  mz.range = x[length(x)] - x[1]
  x2 = seq(x[1],x[length(x)],length=mz.range*int.factor)
  
#   # diff(diff(x)) essentially computes the discrete analogue of the second derivative, so should be negative at local maxima.
#   # The +1 below takes care of the fact that the result of diff is shorter than the input vector.
#   index=which(diff(sign(diff(y)))==-2)+1
  
  # Alway try to fit one curve first
  index = which(y==max(y))

  if (scanmode=="positive"){
    plot_label="pos_fit.png"
  } else {
    plot_label="neg_fit.png"
  }

  # if (plot) {
  #   CairoPNG(filename=paste(outdir,"Gaussian_fit",paste(sampname, x[1], plot_label, sep="_"), sep="/"), width, height)
  #   plot(x,y,xlab="m/z",ylab="I", ylim=c(0,1.5*max(y)),main=paste(i,start,end, sep=" ")) #, ylim=c(0,1e+05)
  # }
  
  retVal = fitGaussian(x2,x,y,index,scale,resol,outdir,force=length(index),useBounds=FALSE,plot,scanmode,int.factor,width,height)

  if (plot) {
    if (length(retVal$mean)==1) {
      
      result = tryCatch(dev.off(), warning = function(w){},
                                   error=function(e){},
                                   finally = {})
      
      # file.rename(paste(outdir,"Gaussian_fit", paste(sampname, x[1], plot_label, sep="_"), sep="/"),
      #             paste(outdir,"Gaussian_fit", paste(sampname, round(retVal$mean,digits = 6), plot_label, sep="_"), sep="/"))
    } else {
      
      
      h2=NULL
      for (i in 1:length(retVal$mean)){
        h2=c(h2, paste("mean =", retVal$mean[i], sep=" "))
      }
      h2=c(h2, paste("fq =", retVal$qual, sep=" "))
      legend("topright", legend=h2)
      
      dev.off()
      
      for (i in 1:length(retVal$mean)){
        if (retVal$mean[i]!=-1){
          file.copy(paste(outdir,"Gaussian_fit", paste(sampname, x[1], plot_label, sep="_"), sep="/"),
                    paste(outdir,"Gaussian_fit", paste(sampname, round(retVal$mean[i],digits = 6), plot_label, sep="_"), sep="/"))
        }
      }
      file.remove(paste(outdir,"Gaussian_fit", paste(sampname, x[1], plot_label, sep="_"), sep="/"))
    }
  }

  return(list("mean"=retVal$mean, "area"=retVal$area, "qual"=retVal$qual, "min"=retVal$min, "max"=retVal$max))
  
}
  