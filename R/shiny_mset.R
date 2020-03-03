is.ordered.mSet <- function(mSet) {
  covarsMatch = all(mSet$dataSet$covars$sample == rownames(mSet$dataSet$norm))
  mSet$dataSet$exp.type <-
    gsub("^1f.",  "1f", mSet$dataSet$exp.type)
  expVarsMatch <- switch(
    mSet$dataSet$exp.type,
    "1f" = {
      all(mSet$dataSet$cls == mSet$dataSet$covars[, mSet$dataSet$exp.fac, with =
                                                    F][[1]])
    },
    "2f" = {
      all(
        mSet$dataSet$facA == mSet$dataSet$covars[, mSet$dataSet$facA.lbl, with =
                                                   F][[1]] &&
          mSet$dataSet$facB == mSet$dataSet$covars[, mSet$dataSet$facB.lbl, with =
                                                     F][[1]]
      )
    },
    "t" = {
      all(
        mSet$dataSet$exp.fac == mSet$dataSet$covars$individual &&
        mSet$dataSet$time.fac == mSet$dataSet$covars[, mSet$dataSet$facA.lbl.orig, with = 
                                                       F][[1]] && 
        mSet$dataSet$facA == mSet$dataSet$time.fac
      )
    },
    "t1f" = {
      all(
        mSet$dataSet$exp.fac == mSet$dataSet$covars[, mSet$dataSet$facA.lbl, with =
                                                      F][[1]] &&
          mSet$dataSet$time.fac == mSet$dataSet$covars[, mSet$dataSet$facB.lbl, with =
                                                         F][[1]]
        &&
          mSet$dataSet$facB == mSet$dataSet$time.fac &&
          mSet$dataSet$facA == mSet$dataSet$covars[, mSet$dataSet$facA.lbl, with =
                                                     F][[1]]
        &&
          mSet$dataSet$facB == mSet$dataSet$covars[, mSet$dataSet$facB.lbl, with =
                                                     F][[1]]
      )
    }
  )
  isOK <- covarsMatch & expVarsMatch
  return(isOK)
}

name.mSet <- function(mSet) {
  info_vec = c()
  if (mSet$dataSet$exp.type %in% c("t1f", "t")) {
    info_vec = c("timeseries")
  } else if (mSet$dataSet$paired) {
    info_vec = c(info_vec, "paired")
  }
  if (length(mSet$dataSet$subset) > 0) {
    subsetgroups = sapply(1:length(mSet$dataSet$subset), function(i) {
      if (names(mSet$dataSet$subset)[i] == "sample") {
        NULL
      } else{
        paste0(
          names(mSet$dataSet$subset)[i],
          "=",
          paste0(mSet$dataSet$subset[[i]], collapse = "+")
        )
      }
    })
  } else{
    subsetgroups = NULL
  }
  
  subsetgroups <- unlist(subsetgroups)
  
  if (length(info_vec) > 0) {
    extra_info <- paste0("(", paste0(info_vec, collapse = " & "), ")")
  } else{
    extra_info <- ""
  }
  
  if (grepl(mSet$dataSet$exp.type, pattern = "^1f")) {
    mSet$dataSet$exp.type <- "1f"
  }
  
  prefix.name <- switch(
    mSet$dataSet$exp.type,
    "1f" = {
      change_var <-
        if (length(mSet$dataSet$exp.var) > 1)
          mSet$dataSet$exp.var[1]
      else
        mSet$dataSet$exp.var
      change_var
    },
    "2f" = {
      # facB should be time if it is there...
      time.check = grepl("time", mSet$dataSet$exp.var)
      is.time = any(time.check)
      if (is.time) {
        print("time series potential")
        idx2 = which(time.check)
        idx1 = setdiff(c(1, 2), idx2)
      } else{
        idx1 = 1
        idx2 = 2
      }
      paste0(mSet$dataSet$exp.var[c(idx1, idx2)], collapse =
               "+")
    },
    "t" = {
      "time"
    },
    "t1f" = {
      change_var <-
        if (length(mSet$dataSet$exp.var) > 1)
          mSet$dataSet$exp.var[1]
      else
        mSet$dataSet$exp.var
      change_var
    }
  )
  mset_name = paste0(prefix.name, extra_info,
                     if (!is.null(subsetgroups))
                       tolower(paste0(":", subsetgroups, collapse = ","))
                     else
                       "")
  mset_name
}

reset.mSet <- function(mSet, fn) {
  origItem <-
    readRDS(file = fn)
  mSet$dataSet <- origItem$data
  mSet$analSet <- origItem$analysis
  mSet$settings <- origItem$settings
  return(mSet)
}

load.mSet <- function(mSet, name = mSet$dataSet$cls.name) {
  mSet$analSet <- mSet$storage[[name]]$analysis
  mSet$settings <- mSet$storage[[name]]$settings
  return(mSet)
}

store.mSet <- function(mSet, name = mSet$dataSet$cls.name) {
  mSet$storage[[name]]$analysis <- mSet$analSet
  mSet$storage[[name]]$settings <- mSet$settings
  return(mSet)
}

change.mSet <-
  function(mSet,
           stats_type,
           stats_var = NULL,
           time_var = NULL) {
    mSet$dataSet$exp.type <- stats_type
    mSet$dataSet$exp.var <- stats_var
    mSet$dataSet$exp.fac <- stats_var
    mSet$dataSet$time.var <- time_var
    if (grepl(mSet$dataSet$exp.type, pattern = "^1f")) {
      mSet$dataSet$exp.type <- "1f"
    }
    
    mSet <- switch(
      mSet$dataSet$exp.type,
      "1f" = {
        change_var <- if (length(stats_var) > 1)
          stats_var[1]
        else
          stats_var
        mSet$dataSet$exp.lbl <- change_var
        # change current variable of interest to user pick from covars table
        mSet$dataSet$cls <-
          as.factor(mSet$dataSet$covars[, ..change_var, with = F][[1]])
        # adjust bivariate/multivariate (2, >2)...
        mSet$dataSet$cls.num <-
          length(levels(mSet$dataSet$cls))
        # - - -
        mSet
      },
      "2f" = {
        time.check = grepl("time", mSet$dataSet$exp.var)
        is.time = any(time.check)
        if (is.time) {
          shiny::showNotification("One variable may be time-related. Using for visualisation...")
          idx2 = which(time.check)
          idx1 = setdiff(c(1, 2), idx2)
        } else{
          idx1 = 1
          idx2 = 2
        }
        mSet$dataSet$facA <-
          as.factor(mSet$dataSet$covars[, stats_var, with = F][[idx1]])
        mSet$dataSet$facB <-
          as.factor(mSet$dataSet$covars[, stats_var, with = F][[idx2]])
        mSet$dataSet$facA.lbl <- stats_var[idx1]
        mSet$dataSet$facB.lbl <- stats_var[idx2]
        mSet$dataSet$exp.type <- "2f"
        mSet$dataSet$exp.lbl <- stats_var
        # - - -
        mSet
      },
      "t" = {
        mSet <- MetaboAnalystR::SetDesignType(mSet, "time0")
        mSet$dataSet$exp.fac <-
          as.factor(mSet$dataSet$covars$individual)
        if (!any(duplicated(mSet$dataSet$exp.fac))) {
          MetaboShiny::metshiAlert(
            "This analysis needs multiple of the same sample in the 'individual' metadata column!"
          )
          return(NULL)
        }
        mSet$dataSet$exp.lbl <- "sample"
        mSet$dataSet$time.fac <- as.factor(mSet$dataSet$covars[, time_var, with = F][[1]])
        mSet$dataSet$exp.type <- "t"
        mSet$dataSet$facA <- mSet$dataSet$time.fac
        mSet$dataSet$facA.lbl <- "Time"
        mSet$dataSet$facA.lbl.orig <- time_var
        mSet$dataSet$paired <- TRUE
        # - - -
        mSet
      },
      "t1f" = {
        print("time series 1 factor")
        mSet <-
          MetaboAnalystR::SetDesignType(mSet, "time")
        if (!any(duplicated(as.factor(mSet$dataSet$covars$individual)))) {
          MetaboShiny::metshiAlert(
            "This analysis needs multiple of the same sample in the 'individual' metadata column!"
          )
          return(NULL)
        }
        change_var <-
          if (length(stats_var) > 1)
            stats_var[1]
        else
          stats_var
        print(time_var)
        
        mSet$dataSet$facA <-
          as.factor(mSet$dataSet$covars[, ..change_var, with = F][[1]])
        mSet$dataSet$facB <-
          as.factor(mSet$dataSet$covars[, ..time_var, with = F][[1]])
        mSet$dataSet$facA.lbl <- change_var
        mSet$dataSet$facB.lbl <- time_var
        mSet$dataSet$exp.fac <- mSet$dataSet$facA
        mSet$dataSet$time.fac <- mSet$dataSet$facB
        mSet$dataSet$exp.type <- "t1f"
        mSet$dataSet$exp.lbl <- change_var
        mSet$dataSet$paired <- TRUE
        
        # - - -
        mSet
      }
    )
    mSet$settings$exp.type = mSet$dataSet$exp.type
    mSet$settings$exp.var = mSet$dataSet$exp.var
    mSet$settings$exp.fac = mSet$dataSet$exp.fac
    mSet$settings$time.var = mSet$dataSet$time.var
    mSet$settings$paired = mSet$dataSet$paired
    # - - - - - - - - - -
    mSet$analSet <- NULL
    return(mSet)
  }

subset.mSet <- function(mSet, subset_var, subset_group) {
  if (!is.null(subset_var)) {
    keep.i <- which(mSet$dataSet$covars[[subset_var]] %in% subset_group)
    keep.samples <- mSet$dataSet$covars$sample[keep.i]
    mSet$dataSet$covars <-
      mSet$dataSet$covars[sample %in% keep.samples, ]
    
    tables = c("norm","prenorm")
    clss = c("cls", "prenorm.cls")
    combi.tbl = data.table::data.table(tbl = tables,
                                       cls = clss)
    
    if (!("subset" %in% names(mSet$settings))) {
      mSet$settings$subset <- list()
    }
    
    for (i in 1:nrow(combi.tbl)) {
      try({
        tbl = combi.tbl$tbl[i]
        cls = combi.tbl$cls[i]
        keep = which(rownames(mSet$dataSet[[tbl]]) %in% keep.samples)
        mSet$dataSet[[tbl]] <- mSet$dataSet[[tbl]][keep, ]
        sampOrder = match(rownames(mSet$dataSet[[tbl]]), mSet$dataSet$covars$sample)
        mSet$dataSet[[cls]] <-
          as.factor(mSet$dataSet$covars[sampOrder, mSet$dataSet$exp.lbl, with = F][[1]])
        if (cls == "cls") {
          mSet$dataSet$cls.num <- length(levels(mSet$dataSet[[cls]]))
          if ("facA" %in% names(mSet$dataSet)) {
            mSet$dataSet$facA <-
              as.factor(mSet$dataSet$covars[sampOrder, mSet$dataSet$facA.lbl, with =
                                              F][[1]])
            mSet$dataSet$facB <-
              as.factor(mSet$dataSet$covars[sampOrder, mSet$dataSet$facB.lbl, with =
                                              F][[1]])
          }
          if ("time.fac" %in% names(mSet$dataSet)) {
            mSet$dataSet$time.fac <-
              as.factor(mSet$dataSet$covars[sampOrder, mSet$dataSet$time.var, with =
                                              F][[1]])
            mSet$dataSet$exp.fac <-
              as.factor(mSet$dataSet$covars[sampOrder, mSet$dataSet$exp.var, with = F][[1]])
          }
        }
      }, silent = T)
    }
    mSet$dataSet$subset[[subset_var]] <- subset_group
    mSet$settings$subset[[subset_var]] <- subset_group
  }
  mSet
}

pair.mSet <- function(mSet) {
  stats_var = mSet$dataSet$exp.var
  time_var = mSet$dataSet$time.var
  
  overview.tbl <-
    data.table::as.data.table(
      cbind(
        sample = mSet$dataSet$covars$sample,
        individual = mSet$dataSet$covars$individual,
        variable = as.character(mSet$dataSet$covars[, ..stats_var][[1]])
      ) 
    )
  
  # which samples should we keep?
  # A. which individuals are present multiple times in the dataset?
  indiv.samp <- unique(overview.tbl[, c("individual", "sample")])
  indiv.var <- unique(overview.tbl[, c("individual", "variable")])
  
  indiv.present.multiple <-
    unique(indiv.samp$individual[c(which(duplicated(indiv.samp$individual)))])
  
  # A2: mandatory for time series 1 factor. which individuals remain in the same category?
  indiv.cat.change = indiv.var$individual[c(which(duplicated(indiv.var$individual)))]
  indiv.no.cat.change <-
    setdiff(indiv.present.multiple, indiv.cat.change)
  
  # B: we want to balance the classes...
  # DOWNSAMPLE invidiv
  if (mSet$dataSet$exp.type == "t1f") {
    indiv.times = unique(data.table::as.data.table(
      cbind(
        individual = mSet$dataSet$covars$individual,
        time = mSet$dataSet$covars[, ..time_var][[1]]
      )
    ))
    indiv.all.times <-
      indiv.times$individual[c(which(duplicated(indiv.times$individual)))]
    # which samples are present
    considerations = list(indiv.present.multiple,
                          indiv.no.cat.change,
                          indiv.all.times)
    times = unique(mSet$dataSet$covars[, ..time_var][[1]])
  } else if (mSet$dataSet$exp.type == "t") {
    considerations = list(indiv.present.multiple)
    times = unique(mSet$dataSet$covars[, ..time_var][[1]])
  } else{
    considerations = list(indiv.present.multiple, indiv.cat.change)
  }
  
  # C. which individuals do we keep?
  keep.indiv <- Reduce(f = intersect, considerations)
  
  req.groups <- unique(mSet$dataSet$covars[, ..stats_var][[1]])
  time_var = mSet$dataSet$time.var
  
  keep.samp <- pbapply::pbsapply(keep.indiv, function(indiv) {
    # get samples matching
    # if multiple for each, only pick one
    overv = overview.tbl[individual == indiv]
    if (mSet$dataSet$exp.type %in% c("t", "t1f")) {
      # all time points need a sample, otherwise nvm (TODO: add a check for if only one sample has timepoint 4 and the rest has 3, majority vote...)
      ind.times = unique(mSet$dataSet$covars[individual == indiv, ..time_var][[1]])
      if (length(ind.times) == length(times)) {
        samps = indiv.samp[individual == indiv]$sample
      } else{
        samps = c()
      }
    } else{
      has.all.groups = all(req.groups %in% overv$variable)
      samps = if (!has.all.groups)
        c()
      else{
        spl.overv <- split(overv, overv$variable)
        minsamp = min(unlist(lapply(spl.overv, nrow)))
        samps = Reduce("c", lapply(req.groups, function(grp) {
          rows = spl.overv[[as.character(grp)]]
          rows$sample[sample(1:nrow(rows), minsamp)]
        }))
      }
    }
    list(samps)
  })
  
  names(keep.samp) <- keep.indiv
  keep.samp = keep.samp[sapply(keep.samp, function(x)
    ! is.null(x))]
  
  if (mSet$dataSet$exp.type == c("t1f")) {
    indiv.final = unique(overview.tbl[sample %in% unlist(keep.samp), c("individual", "variable")])
    downsampled = caret::downSample(indiv.final$individual, as.factor(indiv.final$variable))
    keep.indiv = downsampled$x
    keep.samp <-
      as.character(unlist(keep.samp[which(names(keep.samp) %in% keep.indiv)]))
  }else{
    keep.samp <- as.character(unlist(keep.samp))
  }
  
  if (length(keep.samp) > 2) {
    mSet <- MetaboShiny::subset.mSet(mSet,
                                     subset_var = "sample",
                                     subset_group = keep.samp)
    mSet$settings$paired <- TRUE
    mSet$dataSet$paired <- TRUE
  } else{
    MetaboShiny::metshiAlert("Not enough samples for paired analysis!")
    return(NULL)
  }
  mSet
}