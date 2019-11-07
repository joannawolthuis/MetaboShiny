prev.mSet <- function(mSet, what, input=input){
  
  mSet <- switch(what,
                 change = {
                   mSet$dataSet$exp.var <- input$stats_var
                   mSet$dataSet$time.var <- input$time_var
                   mSet$dataSet$exp.type <- input$stats_type
                   mSet
                 },
                 subset = {
                   mSet$dataSet$subset[[input$subset_var]] <- input$subset_group
                   mSet
                 },
                 pair = {
                   mSet$dataSet$paired <- TRUE
                   mSet
                 },
                 unpair = {
                   mSet$dataSet$paired <- FALSE
                   mSet
                 },
                 unsubset = {
                   mSet$dataSet$subset <- list()
                   mSet
                 })
  return(mSet)
}

# function to generate names for msets
name.mSet <- function(mSet){
  info_vec = c()
  if(mSet$dataSet$exp.type %in% c("t1f", "t")){
    info_vec = c("timeseries")
  }
  if(mSet$dataSet$paired){
    info_vec = c(info_vec, "paired")
  }
  if(length(mSet$dataSet$subset) > 0){
    subsetgroups = sapply(1:length(mSet$dataSet$subset), function(i){
      if(names(mSet$dataSet$subset)[i] == "sample"){
        NULL
      }else{
        paste0(names(mSet$dataSet$subset)[i], "=", paste0(mSet$dataSet$subset[[i]], collapse="+"))
      }
      })  
  }else{
    subsetgroups = NULL
  }
  
  if(length(info_vec)>0){
    extra_info <- paste0("(", info_vec, ")", collapse = " & ")
  }else{
    extra_info <- ""
  }
  
  if(grepl(mSet$dataSet$exp.type, pattern = "^1f")){
    mSet$dataSet$exp.type <- "1f"
  }
  
  prefix.name <- switch(mSet$dataSet$exp.type,
                        "1f"={
                          change_var <- if(length(mSet$dataSet$exp.var)>1) mSet$dataSet$exp.var[1] else mSet$dataSet$exp.var
                          change_var
                        },
                        "2f"={
                          # facB should be time if it is there...
                          time.check = grepl("time", mSet$dataSet$exp.var)
                          is.time = any(time.check)
                          if(is.time){
                            print("time series potential")
                            idx2 = which(time.check)
                            idx1 = setdiff(c(1,2), idx2)
                          }else{
                            idx1=1
                            idx2=2
                          }
                          paste0(mSet$dataSet$exp.var[c(idx1,idx2)],collapse="+")
                        },
                        "t"={
                          "time"
                        },
                        "t1f"={
                          change_var <- if(length(mSet$dataSet$exp.var)>1) mSet$dataSet$exp.var[1] else mSet$dataSet$exp.var
                          paste0(change_var," in time")
                        })
  
  mset_name = paste0(prefix.name, extra_info, if(!is.null(subsetgroups)) ":" else "", tolower(paste0(subsetgroups, collapse=",")))
  
  mset_name
}

load.mSet <- function(mSet, name){
  mSet$dataSet <- mSet$storage[[name]]$data
  mSet$analSet <- mSet$storage[[name]]$analysis
  return(mSet)
}

store.mSet <- function(mSet, name=mSet$dataSet$cls.name){
  mSet$storage[[name]] <- list(data = mSet$dataSet,
                               analysis = mSet$analSet)
  return(mSet)
}

change.mSet <- function(mSet, stats_type, stats_var=NULL, time_var=NULL){
  
  mSet$dataSet$exp.var <- stats_var
  mSet$dataSet$time.var <- time_var
  
  if(grepl(mSet$dataSet$exp.type, pattern = "^1f")){
    mSet$dataSet$exp.type <- "1f"
  }
  
  mSet <- switch(mSet$dataSet$exp.type,
                   "1f"={
                     change_var <- if(length(stats_var)>1) stats_var[1] else stats_var
                     # change current variable of interest to user pick from covars table
                     mSet$dataSet$cls <- as.factor(mSet$dataSet$covars[,change_var, with=F][[1]])
                     # adjust bivariate/multivariate (2, >2)...
                     mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
                     # - - - 
                     mSet
                   },
                   "2f"={
                     time.check = grepl("time", mSet$dataSet$exp.var)
                     is.time = any(time.check)
                     if(is.time){
                       print("time series potential")
                       idx2 = which(time.check)
                       idx1 = setdiff(c(1,2), idx2)
                     }else{
                       idx1=1
                       idx2=2
                     }
                     mSet$dataSet$facA <- as.factor(mSet$dataSet$covars[,stats_var, with=F][[idx1]])
                     mSet$dataSet$facB <- as.factor(mSet$dataSet$covars[,stats_var, with=F][[idx2]])
                     mSet$dataSet$facA.lbl <- stats_var[idx1]
                     mSet$dataSet$facB.lbl <- stats_var[idx2]
                     mSet$dataSet$exp.type <- "2f"
                     # ONLY ANOVA2
                     mSet$dataSet$paired <- F
                     # - - - 
                     mSet
                   },
                   "t"={
                     mSet <- MetaboAnalystR::SetDesignType(mSet, "time0")
                     mSet$dataSet$sbj <- as.factor(mSet$dataSet$covars$individual)
                     if(!any(duplicated(mSet$dataSet$sbj))){
                       print("Won't work, need multiple of the same sample in the 'individual' metadata column!")
                       return(NULL)
                     }
                     mSet$dataSet$time.fac <- as.factor(mSet$dataSet$covars[,time_var, with=F][[1]])
                     mSet$dataSet$exp.type <- "t"
                     mSet$paired <- F
                     
                     # - - - 
                     mSet
                   },
                   "t1f"={
                     mSet <- MetaboAnalystR::SetDesignType(mSet, "time")
                     if(!any(duplicated(as.factor(mSet$dataSet$covars$individual)))){
                       print("Won't work, need multiple of the same sample in the 'individual' metadata column!")
                       return(NULL)
                     }
                     change_var <- if(length(stats_var)>1) stats_var[1] else stats_var
                     
                     mSet$dataSet$facA <- as.factor(mSet$dataSet$covars[,change_var, with=F][[1]])
                     mSet$dataSet$facB <- as.factor(mSet$dataSet$covars[,time_var, with=F][[1]])
                     mSet$dataSet$facA.lbl <- change_var
                     mSet$dataSet$facB.lbl <- "time"
                     mSet$dataSet$exp.fac <- mSet$dataSet$facA
                     mSet$dataSet$time.fac <- mSet$dataSet$facB
                     mSet$dataSet$exp.type <- "t1f"
                     mSet$paired <- F
                     # - - - 
                     mSet
                   })
  mSet$analSet <- NULL
  return(mSet)
}

subset.mSet <- function(mSet, subset_var, subset_group, name=mSet$dataSet$cls.name){
  mSet$dataSet$cls.name <- name
  if(!is.null(subset_var)){
    keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[subset_var]] %in% subset_group)]
    mSet$dataSet$covars <- mSet$dataSet$covars[sample %in% keep.samples]
    keep.proc <- rownames(mSet$dataSet$proc) %in% keep.samples
    mSet$dataSet$proc <- mSet$dataSet$proc[keep.proc,]
    keep.norm <- rownames(mSet$dataSet$norm) %in% keep.samples
    mSet$dataSet$norm <- mSet$dataSet$norm[keep.norm,]
    keep.preproc <- rownames(mSet$dataSet$preproc) %in% keep.samples
    mSet$dataSet$preproc <- mSet$dataSet$preproc[keep.preproc,]
    mSet$dataSet$cls <- droplevels(mSet$dataSet$cls[keep.norm])
    mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
    if("facA" %in% names(mSet$dataSet)){
      mSet$dataSet$facA <- mSet$dataSet$facA[keep.norm]
      mSet$dataSet$facB <- mSet$dataSet$facB[keep.norm]
      mSet$dataSet$time.fac <- mSet$dataSet$time.fac[keep.norm]
      mSet$dataSet$exp.fac <- mSet$dataSet$exp.fac[keep.norm]
    }
    mSet$dataSet$subset[[subset_var]] <- subset_group
    print(paste0("Samples left: ", length(keep.samples)))
  }
  mSet
}