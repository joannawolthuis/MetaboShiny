remove.dupl.2.1 <- function(peaklist) {
# peaklist = outlist.tot
# resultDir = "./results"  
# scanmode = "positive"  

  # peaklist=peaklist[1:100000,]
  
  collect=NULL
  remove=NULL
  
  index = which(duplicated(peaklist[, "mzmed.pgrp"]))
  
  while (length(index) > 0){
    
    peaklist_index = which(peaklist[, "mzmed.pgrp"] == peaklist[index[1], "mzmed.pgrp"])
    tmp=peaklist[peaklist_index[1],,drop=FALSE]
    if (!is.na(peaklist[peaklist_index[1],"assi_HMDB"])) tmp[,"assi_HMDB"]=paste(peaklist[peaklist_index,"assi_HMDB"],collapse = ";") else tmp[,"assi_HMDB"]=NA
    if (!is.na(peaklist[peaklist_index[1],"iso_HMDB"])) tmp[,"iso_HMDB"]=paste(peaklist[peaklist_index,"iso_HMDB"],collapse = ";") else tmp[,"iso_HMDB"]=NA
    if (!is.na(peaklist[peaklist_index[1],"HMDB_code"])) tmp[,"HMDB_code"]=paste(peaklist[peaklist_index,"HMDB_code"],collapse = ";") else tmp[,"HMDB_code"]=NA
    if (peaklist[peaklist_index[1],"assi_noise"]!="") tmp[,"assi_noise"]=paste(peaklist[peaklist_index,"assi_noise"],collapse = ";") else tmp[,"assi_noise"]=""
    
    collect = rbind(collect, tmp)
    remove = c(remove, peaklist_index)

    index=index[-which(index==index[1])]
  } 
  
  peaklist = peaklist[-remove,]
  peaklist = rbind(peaklist,collect)
  
  return(peaklist)
}
