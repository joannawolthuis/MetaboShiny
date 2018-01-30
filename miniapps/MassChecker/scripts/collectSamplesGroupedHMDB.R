collectSamplesGroupedHMDB <- function(resultDir, scanmode, ppm=2){
  # filepath =  paste(resultDir, "grouping_hmdb_done", sep="/")
  # files = list.files(filepath,recursive=TRUE, full.names=TRUE, pattern=paste("*_",scanmode,".RData",sep=""))
  # 
  # index = NULL
  # for (i in 1:length(files)) {
  #   #message(files[i])
  #   load(files[i])
  #   index = c(index, which(outlist.copy[,"height.pkt"]==-1))
  # }
  # 
  # load(paste(resultDir, "specpks_all", paste(scanmode, "RData", sep="."), sep="/"))
  # 
  # remove = which(outlist.tot[,"mzmed.pkt"] %in% outlist.copy[index,"mzmed.pkt"])
  # outlist.rest = outlist.tot[-remove,]
  # 
  # outdir=paste(resultDir, "specpks_all_rest", sep="/")
  # dir.create(outdir)
  
  f <- file.path(outdir, "specpks_all_rest", paste0(scanmode, ".RData"))

  # sort on mass
  outlist = outlist.rest[order(as.numeric(outlist.rest[,"mzmed.pkt"])),]

  n=dim(outlist)[1]
  sub=10000
  end=0
  min_1_last=sub
  check=0
  outlist_i_min_1=NULL

  if (n>=sub & (floor(n/sub)-1)>=2){
    for (i in 2:floor(n/sub)-1){
      start=-(sub-1)+i*sub
      end=i*sub

      if (i>1){
        outlist_i = outlist[c(start:end),]

        n_moved = 0

        # Calculate 3ppm and replace border, avoid cut within peakgroup!
        while ((as.numeric(outlist_i[1,"mzmed.pkt"]) - as.numeric(outlist_i_min_1[min_1_last,"mzmed.pkt"]))*1e+06/as.numeric(outlist_i[1,"mzmed.pkt"]) < 3) {
          outlist_i_min_1 = rbind(outlist_i_min_1, outlist_i[1,])
          outlist_i = outlist_i[-1,]
          n_moved = n_moved + 1
        }

        # message(paste("Process", i-1,":", dim(outlist_i_min_1)[1]))
        save(outlist_i_min_1, file=paste(outdir, paste(scanmode, paste("outlist_i_min_1",i-1,"RData", sep="."), sep="_"), sep="/"))
        check=check+dim(outlist_i_min_1)[1]

        outlist_i_min_1 = outlist_i
        min_1_last = dim(outlist_i_min_1)[1]

      } else {
        outlist_i_min_1 = outlist[c(start:end),]
      }
    }
  }

  start = end + 1
  end = n
  outlist_i = outlist[c(start:end),]
  n_moved = 0

  if(!is.null(outlist_i_min_1)){
    # Calculate 4ppm and replace border, avoid cut within peakgroup!
    while ((as.numeric(outlist_i[1,"mzmed.pkt"]) - as.numeric(outlist_i_min_1[min_1_last,"mzmed.pkt"]))*1e+06/as.numeric(outlist_i[1,"mzmed.pkt"]) < 2*ppm) {
      outlist_i_min_1 = rbind(outlist_i_min_1, outlist_i[1,])
      outlist_i = outlist_i[-1,]
      n_moved = n_moved + 1
    }

    # message(paste("Process", i+1-1,":", dim(outlist_i_min_1)[1]))
    save(outlist_i_min_1, file=paste(outdir, paste(scanmode, paste("outlist_i_min_1",i,"RData", sep="."), sep="_"), sep="/"))
    check=check+dim(outlist_i_min_1)[1]
  }

  outlist_i_min_1=outlist_i
  # message(paste("Process", i+2-1,":", dim(outlist_i_min_1)[1]))
  save(outlist_i_min_1, file=paste(outdir, paste(scanmode, paste("outlist_i_min_1",i+1,"RData", sep="."), sep="_"), sep="/"))
  check=check+dim(outlist_i_min_1)[1]
  
  if (check==dim(outlist)[1]){
    message(paste("Check is oke!"))
  } else {
    message(paste("Check is failed!"))
  }

}

# message("Start")
# cat("==> reading arguments:\n", sep = "")
# 
# cmd_args = commandArgs(trailingOnly = TRUE)
# 
# for (arg in cmd_args) cat("  ", arg, "\n", sep="")
# 
# run(cmd_args[1], cmd_args[2])
# #run("./results")
# 
# message("Ready")