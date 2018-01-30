iden.code <- function(peaklist, db, ppm, theor_mass_label) {
  # theor_mass_label = {"MNeg", "Mpos"}
  
  mcol <- peaklist[ , "mzmed.pgrp"]
  theor.mcol <- db[,theor_mass_label]
  assi_HMDB <- iso_HMDB <- HMDB_code <- c(rep("",length(mcol)))
  theormz_HMDB <- c(rep(0,length(mcol)))
  
  for(t in 1:length(mcol)){
    mz<-mcol[t]
    mtol <- mz*ppm/1000000
    selp <- which(theor.mcol > (mz - mtol) & theor.mcol < (mz + mtol))
    assinames <- isonames <- HMDBcodes <- ""
    
    if(length(selp)!=0){
      for(i in 1:length(selp)){

        if(grepl(" iso ", db[selp[i], "CompoundName"])){
          mainpeak <- strsplit(db[selp[i],"CompoundName"]," iso ")[[1]][1]
          
          # Check if peak without isotope occurs, this assumes that peaklist is ordered on mass!!! 
          if(length(grep(mainpeak, assi_HMDB, fixed=TRUE)) > 0){
            isonames <- as.character(paste(isonames,as.character(db[selp[i],"CompoundName"]), sep=";"))
          }
          
        } else {
          assinames <- as.character(paste(assinames,as.character(db[selp[i], "CompoundName"]), sep=";"))
          HMDBcodes <- as.character(paste(HMDBcodes, as.character(rownames(db)[selp[i]]), sep=";"))
        }  
      }
    }  
    
    assi_HMDB[t] <- as.character(assinames)
    iso_HMDB[t] <-  as.character(isonames)
    if ((assinames=="") & (isonames=="")){ 
      theormz_HMDB[t] <- ""
    } else {
      theormz_HMDB[t] <- theor.mcol[selp[1]]
    }
    HMDB_code[t] <- as.character(HMDBcodes)
  }
  
  return(cbind(peaklist, assi_HMDB, iso_HMDB, theormz_HMDB, HMDB_code))
}

