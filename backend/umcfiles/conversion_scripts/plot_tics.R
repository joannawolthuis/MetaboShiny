library(xcms)
library(pbapply)
library(parallel)
library(data.table)
library(xlsx)


cl = makeCluster(4, type="FORK") 
##  NB: peak finding, peak group finding and fill missing done on HPC !
CheckCDFfile<-function(file, type=".cdf"){
  print(paste("Loading File:", file, sep=""))
  xr<-xcmsRaw(file, profstep=0)
  for(i in 1:length(xr@scanindex)){
    scan<-getScan(xr, scan=i)
    if(is.unsorted(scan[,"mz"]) == TRUE){
      print(" corrupt, fixing ")
      newfile<-sub(type, "-Fixed.mzdata", file, ignore.case=TRUE, fixed=TRUE)
      write.mzdata(xr, newfile)
      file.copy(file, sub(type, ".OLD", file, ignore.case=TRUE))
      unlink(file)
      return(1)
    }
    if(i == length(xr@scanindex)){
      print(" is ok ")
      return(0)
    }
  }
}
sapply(xmlfiles, CheckCDFfile)

# check sorting of files - can be corrupt
outdir <- "/Users/jwolthuis/PROCESSING_HPC"   # project dir E:\Metabolomics\Project2016_03_biggen

xmlfiles<-list.files(outdir,recursive=TRUE, pattern="mzXML", ignore.case=TRUE, full.names=TRUE)

rawfiles<-list.files(outdir,recursive=TRUE, pattern="raw", ignore.case=TRUE, full.names=TRUE)
TICdir <- paste(outdir, "/TICs", sep="")
dir.create(TICdir)


plotTICtoFolder <- function(file, TICdir=TICdir){
  print(file)
  # -------------
  rawF <- NULL
  attempt <- 1
  while( is.null(rawF) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
     rawF <- xcmsRaw(file, profstep=0.1)
    )
  }
  # -----------
  png(filename=paste(TICdir, "/", basename(file), ".png", sep=""), 320, 240)
  plotTIC(rawF, ident=FALSE, msident=FALSE)
  dev.off()
}

getMLrows <- function(file){
  print(file)
  # -------------
  rawF <- NULL
  attempt <- 1
  while( is.null(rawF) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
      rawF <- xcmsRaw(file, profstep=0.1)
    )
  }
  judgeN <- 0
  plotTIC(rawF, ident=FALSE, msident=FALSE)
  while(!judgeN %in% c(1, 2, 3))
  {
    judgeN <- as.integer(readline(prompt="Judge tic (OK=1, DOUBT=2, BAD=3): "))
  }
  judgement <- if(judgeN == 1) "OK" else if(judgeN == 2) "DOUBT" else if(judgeN == 3) "BAD"
  MLrow <- as.data.table(t(c(basename(file),judgement, rawF@tic)))
  colnames(MLrow) = c("filename", "judgement", paste0("sec", 1:length(rawF@tic)))
  # -----------------
  return(MLrow)
}

x11()
ML_list <- lapply(xmlfiles, getMLrows)

save(ML_list, file = "TIC_ml.RData")
