sumAdducts <- function(peaklist, theor.MZ, grpnames.long, adducts, batch, scanmode, outdir){
# theor.MZ = outlist_part
# grpnames.long = groupNames
# peaklist = outlist.tot
  
  dir.create(paste(outdir, "adductSums", sep="/"), showWarnings = FALSE)
  
  hmdb_codes = rownames(theor.MZ)
  hmdb_names = theor.MZ[,1]
  
  # remove isotopes!!!
  index = grep("HMDB",hmdb_codes,fixed=TRUE)
  hmdb_codes = hmdb_codes[index]
  hmdb_names = hmdb_names[index]
  index = grep("_",hmdb_codes,fixed=TRUE)
  if (length(index)>0) hmdb_codes = hmdb_codes[-index]
  if (length(index)>0) hmdb_names = hmdb_names[-index]
  
  #i=which(hmdb_codes=="HMDB41792")
  
  # negative
  names=NULL
  adductsum=NULL
  names_long=NULL
  
  if (length(hmdb_codes)!=0) {
  
  # assign("last.warning", NULL, envir = baseenv())
  # result = tryCatch(
  #   {

      for(i in 1:length(hmdb_codes)){
      # for(i in 1:10){
        
        #compound="HMDB00045"
        compound=hmdb_codes[i]
        compound_plus=c(compound,paste(compound, adducts, sep = "_"))
        
        # x=peaklist$HMDB_code[1]
        metab=unlist(lapply(peaklist$HMDB_code, function(x) {(length(intersect(unlist(strsplit(as.vector(x),";")),compound_plus))>0)}))
        # peaklist[metab, "assi.hmdb"]
        # which(metab==TRUE)
        
        # if (length(which(metab==TRUE))>0) message("Bingo found something")
        
        total=c()
        
        # peaklist[metab, c("mzmed.pgrp", "HMDB_code", "C34.1")]
        ints=peaklist[metab, c(7:(length(grpnames.long)+6))]
        total=apply(ints, 2, sum)
        
        if (sum(total)!=0){
          # message(i)
          names = c(names, compound)
          adductsum<-rbind(adductsum,total)
          names_long = c(names_long, hmdb_names[i])
        }  
      }
  
    # }
    # , warning=function(w) {
    #   message(paste("CATCHED", w))
    # }
    # , error = function(e) {
    #   message(paste("CATCHED", e))
    # })
  
    if (!is.null(adductsum)){ 
      colnames(adductsum)=grpnames.long
      rownames(adductsum)=names
      adductsum = cbind(adductsum, "HMDB_name"=names_long)
      save(adductsum, file=paste(outdir, "adductSums", paste(scanmode, "_",batch,".RData", sep=""), sep="/"))
    }
  }  
}