prev.mSet <- function(mSet, what, input=input){
  
  mSet <- switch(what,
                 change = {
                   mSet$dataSet$exp.var <- input$stats_var
                   mSet$dataSet$time.var <- input$time_var
                   mSet$dataSet$exp.type <- input$stats_type
                   if(input$stats_type %in% c("t", "t1f") | input$paired){
                     mSet$dataSet$paired <- TRUE 
                   }else{
                     mSet$dataSet$paired <- FALSE
                   }
                   mSet
                 },
                 subset = {
                   mSet$dataSet$subset[[input$subset_var]] <- input$subset_group
                   mSet
                 },
                 unsubset = {
                   mSet$dataSet$subset <- list()
                   mSet
                 })
  return(mSet)
}

name.mSet <- function(mSet){
  info_vec = c()
  if(mSet$dataSet$exp.type %in% c("t1f", "t")){
    info_vec = c("timeseries")
  }else if(mSet$dataSet$paired){
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
  
  subsetgroups <- unlist(subsetgroups)
  
  if(length(info_vec)>0){
    extra_info <- paste0("(", paste0(info_vec, collapse = " & "), ")")
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
                          change_var
                        })
  mset_name = paste0(prefix.name, extra_info, 
                     if(!is.null(subsetgroups)) tolower(paste0(":", subsetgroups, collapse=",")) else "")
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
  
  mSet.orig <- mSet
  
  mSet$dataSet$exp.type <- stats_type
  mSet$dataSet$exp.var <- stats_var
  mSet$dataSet$time.var <- time_var
  
  if(grepl(mSet$dataSet$exp.type, pattern = "^1f")){
    mSet$dataSet$exp.type <- "1f"
  }
  
  print(mSet$dataSet$exp.type)
  
  mSet <- switch(mSet$dataSet$exp.type,
                   "1f"={
                     change_var <- if(length(stats_var)>1) stats_var[1] else stats_var
                     mSet$dataSet$exp.lbl <- change_var
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
                       shiny::showNotification("One variable may be time-related. Using for visualisation...")
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
                     mSet$dataSet$exp.lbl <- stats_var
                     # - - - 
                     mSet
                   },
                   "t"={
                     mSet <- MetaboAnalystR::SetDesignType(mSet, "time0")
                     mSet$dataSet$exp.fac <- as.factor(mSet$dataSet$covars$individual)
                     if(!any(duplicated(mSet$dataSet$exp.fac))){
                       MetaboShiny::metshiAlert("This analysis needs multiple of the same sample in the 'individual' metadata column!")
                       return(NULL)
                     }
                     mSet$dataSet$exp.lbl <- "sample"
                     mSet$dataSet$time.fac <- as.factor(mSet$dataSet$covars[,time_var, with=F][[1]])
                     mSet$dataSet$exp.type <- "t"

                     # - - - 
                     mSet
                   },
                   "t1f"={
                     mSet <- MetaboAnalystR::SetDesignType(mSet, "time")
                     if(!any(duplicated(as.factor(mSet$dataSet$covars$individual)))){
                       MetaboShiny::metshiAlert("This analysis needs multiple of the same sample in the 'individual' metadata column!")
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
                     mSet$dataSet$exp.lbl <- change_var
                     # - - - 
                     mSet
                   })
  mSet$analSet <- NULL
  if(is.null(mSet)) mSet <- mSet.orig 
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
      mSet$dataSet$facA <- droplevels(mSet$dataSet$facA[keep.norm])
      mSet$dataSet$facB <- droplevels(mSet$dataSet$facB[keep.norm])
    }
    if("time.fac" %in% names(mSet$dataSet)){
      mSet$dataSet$time.fac <- droplevels(mSet$dataSet$time.fac[keep.norm])
      mSet$dataSet$exp.fac <- droplevels(mSet$dataSet$exp.fac[keep.norm])
    }
    mSet$dataSet$subset[[subset_var]] <- subset_group
  }
  mSet
}

pair.mSet <- function(mSet, name){
  
  stats_var = mSet$dataSet$exp.var  
  
  overview.tbl <- data.table::as.data.table(cbind(sample = mSet$dataSet$covars$sample,
                                                  individual = mSet$dataSet$covars$individual,
                                                  variable = as.character(
                                                    mSet$dataSet$covars[, ..stats_var][[1]])
                                                  ))
  # which samples should we keep?
  # A. which individuals are present multiple times in the dataset?
  indiv.samp <- unique(overview.tbl[, c("individual","sample")])
  indiv.var <-  unique(overview.tbl[, c("individual","variable")])
  
  indiv.present.multiple <- unique(indiv.samp$individual[c(which(duplicated(indiv.samp$individual)))])

  # A2: mandatory for time series 1 factor. which individuals remain in the same category?
  indiv.cat.change = indiv.var$individual[c(which(duplicated(indiv.var$individual)))] 
  indiv.no.cat.change <- setdiff(indiv.present.multiple, indiv.cat.change)

  # B: we want to balance the classes...
  # DOWNSAMPLE invidiv
  if(mSet$dataSet$exp.type == "t1f"){
    considerations = list(indiv.present.multiple, indiv.no.cat.change)
    times = unique(mSet$dataSet$covars[, ..time_var][[1]])
  }else if(mSet$dataSet$exp.type == "t"){
    considerations = list(indiv.present.multiple)
    times = unique(mSet$dataSet$covars[, ..time_var][[1]])
  }else{
    considerations = list(indiv.present.multiple, indiv.cat.change)
  }
  
  # C. which individuals do we keep?
  keep.indiv <- Reduce(f = intersect, considerations)
  
  if(mSet$dataSet$exp.type %in% c("t","t1f")){
    keep.indiv.final = unique(overview.tbl[sample %in% keep.samp, c("individual", "variable")])
    keep.indiv = caret::downSample(keep.indiv.final$individual, as.factor(keep.indiv.final$variable))$x
  }
  
  req.groups <- unique(mSet$dataSet$covars[, ..stats_var][[1]])
  
  keep.samp <- unlist(pbapply::pbsapply(keep.indiv, function(indiv){
    # get samples matching
    # if multiple for each, only pick one
    overv = overview.tbl[individual == indiv]
    if(mSet$dataSet$exp.type %in% c("t", "t1f")){
      # all time points need a sample, otherwise nvm (TODO: add a check for if only one sample has timepoint 4 and the rest has 3, majority vote...)
      time_var = mSet$dataSet$time.var
      ind.times = mSet$dataSet$covars[individual == indiv, ..time_var][[1]]
      if(length(intersect(ind.times, times)) == length(times)){
        indiv.samp[individual == indiv]$sample
      }else{
        c()
      }
    }else{
      has.all.groups = all(req.groups %in% overv$variable)
      samps = if(!has.all.groups) c() else{
        spl.overv <- split(overv, overv$variable)
        minsamp = min(unlist(lapply(spl.overv, nrow)))
        Reduce("c",lapply(req.groups, function(grp){
          rows = spl.overv[[grp]]
          rows$sample[sample(1:nrow(rows), minsamp)]
        }))
      }
    } 
  }))
  
  if(length(keep.samp)> 2){
    mSet <- MetaboShiny::subset.mSet(mSet, 
                                     subset_var = "sample", 
                                     subset_group = keep.samp)
    mSet$dataSet$paired <- TRUE  
  }else{
    MetaboShiny::metshiAlert("Not enough samples for paired analysis!")
    return(NULL)
  }  
  mSet
}