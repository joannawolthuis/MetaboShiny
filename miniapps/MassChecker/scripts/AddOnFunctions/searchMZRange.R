searchMZRange <- function(range,values,int.factor,scale,resol,outdir,sampname,scanmode,plot,width,height,thresh){
# range=sub_range

  end=NULL
  index=as.vector(which(range!=0))
  
  # bad infusion
  if (length(index)==0) return(values)
  
  start=index[1]
  subRangeLength = 15
  
  pb <- pbapply::startpb(1, length(index))

  for (i in 1:length(index)){
    pbapply::setpb(pb, i)
    #print(paste(i / length(index) * 100.00, "%"))
  # for (i in 1:129000){

    if (i<length(index) & (index[i+1] - index[i]) > 1){
      
      end=index[i]
      
      # start=395626
      # end=395640
      # i=25824

      # gaan met de banaan
      # 128836
      # start 1272853
      # mz start 316.83302
      # end 1272870
      # mz end 316.84269
      
      x = as.numeric(names(range)[c(start:end)])
      y = as.vector(range[c(start:end)])
      
#     # Trim zeros
#     x = as.vector(trimZeros(x,y)[[1]])
#     y = as.vector(trimZeros(x,y)[[2]])
      
      if (length(y)!=0) {
        
        # check if intensity above thresh
        if (max(y) < thresh | is.nan(max(y))) {
          start=index[i+1]
          next
        }

        # message("gaan met de banaan")
        # message(i)
        # message(paste("start", start, sep=" "))
        # message(paste("mz start", x[1], sep=" "))
        # message(paste("end", end, sep=" "))
        # message(paste("mz end", x[length(x)], sep=" "))
        
        if (length(y)>subRangeLength) {
          
          y[which(y<min(y)*1.1)] = 0
          sub_range = y
          names(sub_range) = x
          
          # Range to long!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          values = searchMZRange(sub_range,values,int.factor,scale,resol,outdir,sampname,scanmode,plot,width,height,thresh)

        } else if (length(y)>3) {
          # Check only zeros
          if (sum(y)==0) next

          rval = fitGaussianInit(x,y,int.factor,scale,resol,outdir,sampname, scanmode, plot,width,height,i,start,end)
          
          if (rval$qual[1]==1) {

            rval = generateGaussian(x,y,resol,plot,scanmode,int.factor, width, height)

            values$mean = c(values$mean, rval$mean)
            values$area = c(values$area, rval$area)
            values$nr = c(values$nr, sampname)
            values$min = c(values$min, rval$min)
            values$max = c(values$max, rval$max)
            values$qual = c(values$qual, 0)

            values$spikes = values$spikes + 1

          } else {
            for (j in 1:length(rval$mean)){
              values$mean = c(values$mean, rval$mean[j])
              values$area = c(values$area, rval$area[j])
              values$nr = c(values$nr, sampname)
              values$min = c(values$min, rval$min[1])
              values$max = c(values$max, rval$max[1])
              values$qual = c(values$qual, rval$qual[1])
            }
          }

        } else {

          rval = generateGaussian(x,y,resol,plot,scanmode,int.factor, width, height)

          values$mean = c(values$mean, rval$mean)
          values$area = c(values$area, rval$area)
          values$nr = c(values$nr, sampname)
          values$min = c(values$min, rval$min)
          values$max = c(values$max, rval$max)
          values$qual = c(values$qual, 0)

          values$spikes = values$spikes + 1
        }
      }
      start=index[i+1]
    } 
  }
  
  # last little range
  end = index[length(index)]
  x = as.numeric(names(range)[c(start:end)])
  y = as.vector(range[c(start:end)])
  
#   x = as.vector(trimZeros(x,y)[[1]])
#   y = as.vector(trimZeros(x,y)[[2]])

  if (length(y)!=0) {
    
    # check if intensity above thresh
    if (max(y) < thresh | is.nan(max(y))) {
      #start=index[i+1]
      # do nothing!!
    } else {
      
#       message("gaan met de banaan")
#       message(paste("start", start, sep=" "))
#       message(paste("mz start", x[1], sep=" "))
#       message(paste("end", end, sep=" "))
#       message(paste("mz end", x[length(x)], sep=" "))
      
      if (length(y)>subRangeLength) {
        
        y[which(y<min(y)*1.1)] = 0
        sub_range = y
        names(sub_range) = x
        
        # Range to long!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        values = searchMZRange(sub_range,values,int.factor,scale,resol,outdir,sampname,scanmode,plot,width,height,thresh)

      } else if (length(y)>3) {
        
        # Check only zeros
        if (sum(y)==0) next 
        
        rval = fitGaussianInit(x,y,int.factor,scale,resol,outdir,sampname,scanmode,plot,width,height,i,start,end)
        
        if (rval$qual[1]==1) {
          #message("Quality = 1!!!")
          rval = generateGaussian(x,y,resol,plot,scanmode,int.factor, width, height)
          
          values$mean = c(values$mean, rval$mean)
          values$area = c(values$area, rval$area)
          values$nr = c(values$nr, sampname)
          values$min = c(values$min, rval$min)
          values$max = c(values$max, rval$max)
          values$qual = c(values$qual, 0)
          
          values$spikes = values$spikes + 1
          
        } else {
          for (j in 1:length(rval$mean)){
            values$mean = c(values$mean, rval$mean[j])
            values$area = c(values$area, rval$area[j])
            values$nr = c(values$nr, sampname)
            values$min = c(values$min, rval$min[1])
            values$max = c(values$max, rval$max[1])
            values$qual = c(values$qual, rval$qual[1])
          }
        }
      } else {
        rval = generateGaussian(x,y,resol,plot,scanmode,int.factor, width, height)
        
        values$mean = c(values$mean, rval$mean)
        values$area = c(values$area, rval$area)
        values$nr = c(values$nr, sampname)
        values$min = c(values$min, rval$min)
        values$max = c(values$max, rval$max)
        values$qual = c(values$qual, 0)
        
        values$spikes = values$spikes + 1
      }
    }
  }
  
  return(values)
  #return(list("mean"=retVal$mean, "area"=retVal$area, "qual"=retVal$qual, "min"=retVal$min, "max"=retVal$max, "nr"=sample.nr, "spikes"=spikes))
}  
