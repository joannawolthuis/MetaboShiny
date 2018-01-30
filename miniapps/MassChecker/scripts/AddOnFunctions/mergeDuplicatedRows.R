mergeDuplicatedRows <- function(peaklist) {
# peaklist = outlist.tot
# resultDir = "./results"
# scanmode = "positive"

  # peaklist_index=which(peaklist[,"mzmed.pgrp"]=="94.9984524624111")
  # peaklist[peaklist_index,]

  collapse <- function(label,pklst,index){
    # label = "iso_HMDB"
    # pklst = peaklist
    # index = peaklist_index
    tmp2=as.vector(pklst[index,label])
    if (length(which(is.na(tmp2)))>0) tmp2=tmp2[-which(is.na(tmp2))]
    return(paste(tmp2,collapse = ";"))
  } 
  
  options(digits=16)
  
  collect=NULL
  remove=NULL
  
  index = which(duplicated(peaklist[, "mzmed.pgrp"]))
  
  while (length(index) > 0){
    
    peaklist_index = which(peaklist[, "mzmed.pgrp"] == peaklist[index[1], "mzmed.pgrp"])
    # peaklist[peaklist_index,"iso_HMDB",drop=FALSE]
    tmp=peaklist[peaklist_index[1],,drop=FALSE]
    
    tmp[,"assi_HMDB"]=collapse("assi_HMDB",peaklist,peaklist_index)
    tmp[,"iso_HMDB"]=collapse("iso_HMDB",peaklist,peaklist_index)
    tmp[,"HMDB_code"]=collapse("HMDB_code",peaklist,peaklist_index)
    tmp[,"assi_noise"]=collapse("assi_noise",peaklist,peaklist_index)
    if (tmp[,"assi_noise"]==";") tmp[,"assi_noise"]=NA
    tmp[,"theormz_noise"]=collapse("theormz_noise",peaklist,peaklist_index)
    if (tmp[,"theormz_noise"]=="0;0") tmp[,"theormz_noise"]=NA
    
    collect = rbind(collect, tmp)
    remove = c(remove, peaklist_index)

    index=index[-which(peaklist[index, "mzmed.pgrp"] == peaklist[index[1], "mzmed.pgrp"])]
  } 
  
  if (!is.null(remove)) peaklist = peaklist[-remove,]
  peaklist = rbind(peaklist,collect)
  
  return(peaklist)
}
