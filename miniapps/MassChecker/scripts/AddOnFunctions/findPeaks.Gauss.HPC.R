### fit Gaussian estimate mean and integrate to obtain intensity
findPeaks.Gauss.HPC <- function(plist, breaks.fwhm, int.factor, scale, resol, outdir, sampname, scanmode, plot, thresh, width, height) {
   # plist=pklist$neg
   # breaks.fwhm=pklist$breaksFwhm
   # label="/specpks/Neg_specpks"
   # scanmode="negative"
   # plot=TRUE

  range = as.vector(plist)
  names(range) = rownames(plist)
  #range[34:43]

  values = list("mean"=NULL,"area"=NULL,"nr"=NULL,"min"=NULL,"max"=NULL,"qual"=NULL,"spikes"=0)
  
  values = searchMZRange(range,values,int.factor,scale,resol,outdir,sampname,scanmode,plot,width,height,thresh)  

  outlist.persample=NULL
  outlist.persample=cbind("samplenr"=values$nr, "mzmed.pkt"=values$mean, "fq"=values$qual, "mzmin.pkt"=values$min, "mzmax.pkt"=values$max, "height.pkt"=values$area)
  #outlist.persample=cbind("samplenr"=sample.nr, "mzmed.pkt"=peak.mean, "fq"=peak.qual, "mzmin.pkt"=peak.min, "mzmax.pkt"=peak.max, "height.pkt"=peak.area)
  index=which(outlist.persample[,"height.pkt"]==0)
  if (length(index)>0){
    outlist.persample=outlist.persample[-index,]  
  }

  save(outlist.persample, file=paste(outdir, "specpks", paste(sampname, "_", scanmode, ".RData", sep=""), sep="/"))

  message(paste("There were", values$spikes, "spikes!"))
  #return(peaklist.all)
}  
