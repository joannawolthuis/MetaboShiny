averageTechReplicates <- function(outdir, nrepl, dimsThresh=100, cores=1){
# outdir="./results"
# nrepl=5
# dimsThresh=100
  
  #thresh2remove = 1*10^9 # plasma
  thresh2remove = 5*10^8 # blood spots
  #thresh2remove = 1*10^8 # research (Mia)  

  removeFromRepl.pat <- function(bad_samples, repl.pattern, groupNames, nrepl) {
    # bad_samples=remove_pos
    
    tmp = repl.pattern
    
    removeFromGroup=NULL
    
    for (i in 1:length(tmp)){
      tmp2 = repl.pattern[[i]]
      
      remove=NULL
      
      for (j in 1:length(tmp2)){
        if (tmp2[j] %in% bad_samples){
          message(tmp2[j])
          message(paste("remove",tmp2[j]))
          message(paste("remove i",i))
          message(paste("remove j",j))
          
          remove = c(remove, j)
        }
      }
      
      if (length(remove)==nrepl) removeFromGroup=c(removeFromGroup,i) 
      if (!is.null(remove)) repl.pattern[[i]]=repl.pattern[[i]][-remove]
    }
    
    if (length(removeFromGroup)!=0) groupNames=groupNames[-removeFromGroup]

    return(list("pattern"=repl.pattern, "groupNames"=groupNames))
    
  }
  
  dir.create(paste(outdir, "specpks", sep="/"),showWarnings = F)
  dir.create(paste(outdir, "Gaussian_fit", sep="/"),showWarnings = F)
  dir.create(paste(outdir, "average_pklist", sep="/"),showWarnings = F)
  
  # get repl.pattern
  load(paste(outdir, "init.RData", sep="/"))

  progress <<- function(i) setProgress(value = i / length(repl.pattern),
                                       detail = paste("Averaged replicates for", basename(groupNames[i])))
  opts <- list(progress=progress)
  # --- cluster ---
  cl <<- makeSOCKcluster(cores,
                         outfile="~/mclog.txt")
  registerDoSNOW(cl)
  # ---------------
  withProgress(message = "Average technical replicates",{
    to.remove <- foreach(i=1:length(repl.pattern), .options.snow=opts,
                   .verbose = T) %dopar% {
                     techRepsArray.pos = NULL
                     techRepsArray.neg = NULL
                     
                     remove_neg = NULL
                     remove_pos = NULL
                     
                     tech_reps = as.vector(unlist(repl.pattern[i]))
                     sum_neg = 0
                     sum_pos = 0
                     n_pos = 0
                     n_neg = 0
                     
                     for (j in 1:length(tech_reps)){
                       load(paste(paste(outdir, "pklist/", sep="/"), tech_reps[j], ".RData", sep=""))
                       
                       if (sum(pklist$neg[,1])<thresh2remove){
                         remove_neg=c(remove_neg, tech_reps[j])
                       } else {
                         n_neg=n_neg+1
                         sum_neg=sum_neg+pklist$neg  
                       }
                       
                       techRepsArray.neg = cbind(techRepsArray.neg, pklist$neg)
                       
                       if (sum(pklist$pos[,1])<thresh2remove){
                         remove_pos=c(remove_pos, tech_reps[j])
                       } else {
                         n_pos=n_pos+1
                         sum_pos=sum_pos+pklist$pos
                       }
                       
                       techRepsArray.pos = cbind(techRepsArray.pos, pklist$pos)
                     }
                     
                     # filter within bins on at least signal in more than one tech. rep.!!!
                     if (!is.null(dim(sum_pos))) sum_pos[apply(techRepsArray.pos,1,function(x) length(which(x>dimsThresh))==1),1]=0
                     if (!is.null(dim(sum_neg))) sum_neg[apply(techRepsArray.neg,1,function(x) length(which(x>dimsThresh))==1),1]=0
                     
                     if (n_neg!=0){
                       sum_neg[,1]=sum_neg[,1]/n_neg
                       colnames(sum_neg)=groupNames[i]
                       save(sum_neg, file=paste(paste(outdir, "average_pklist", sep="/"),"/", groupNames[i], "_neg.RData", sep=""))
                     }
                     if (n_pos!=0){
                       sum_pos[,1]=sum_pos[,1]/n_pos
                       colnames(sum_pos)=groupNames[i]
                       save(sum_pos, file=paste(paste(outdir, "average_pklist", sep="/"),"/", groupNames[i], "_pos.RData", sep=""))
                     }
                     
                     # --- return ---
                     list(pos=remove_pos, 
                          neg=remove_neg)
                   }
  })
  stopCluster(cl)
  
  remove_pos <- sapply(to.remove, function(x) x$pos)
  remove_pos <- unlist(remove_pos[!sapply(remove_pos, is.null)])
  remove_neg <- sapply(to.remove, function(x) x$neg)
  remove_neg <- unlist(remove_pos[!sapply(remove_neg, is.null)])

  retVal = removeFromRepl.pat(remove_pos, repl.pattern, groupNames, nrepl)
  repl.pattern.pos = retVal$pattern
  groupNames.pos = retVal$groupNames
  write.table(remove_pos, file=paste(outdir, "miss_infusions_pos.txt", sep="/"), row.names=FALSE, col.names=FALSE ,sep= "\t")

  retVal = removeFromRepl.pat(remove_neg, repl.pattern, groupNames, nrepl)
  repl.pattern.neg = retVal$pattern
  groupNames.neg = retVal$groupNames
  write.table(remove_neg, file=paste(outdir, "miss_infusions_neg.txt", sep="/"), row.names=FALSE, col.names=FALSE ,sep= "\t")
  
  save(repl.pattern.pos,repl.pattern.neg,groupNames.pos,groupNames.neg, file=paste(outdir, "repl.pattern.RData", sep="/"))
}
# 
# message("Start")
# cat("==> reading arguments:\n", sep = "")
# 
# cmd_args = commandArgs(trailingOnly = TRUE)
# 
# for (arg in cmd_args) cat("  ", arg, "\n", sep="")
# 
# run(cmd_args[1], as.numeric(cmd_args[2]))
# 
# message("Ready")
