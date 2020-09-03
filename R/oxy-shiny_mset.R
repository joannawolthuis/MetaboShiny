#' @title Check if mSet sample order in peaktable data and metadata are the same
#' @param mSet mSet object
#' @return TRUE/FALSE
#' @rdname is.ordered.mSet
#' @export 
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

#' @title Generate mSet name
#' @description For MetShi display purposes. Takes subset and selected variable(s) to generate a name for this subexperiment.
#' @param mSet mSet object
#' @return character object
#' @rdname name.mSet
#' @export 
name.mSet <- function(mSet) {
  info_vec = c()
  if (mSet$dataSet$exp.type %in% c("t1f", "t")) {
    info_vec = c("timeseries")
  } else if (mSet$settings$paired) {
    info_vec = c(info_vec, "paired")
  }
  
  if("mz" %in% names(mSet$settings$subset)){
    info_vec = c(info_vec, "prematched m/z only")
  }
  
  if (length(mSet$settings$subset) > 0) {
    subsetgroups = sapply(1:length(mSet$settings$subset), function(i) {
      if (names(mSet$settings$subset)[i] %in% c("sample", "mz")) {
        NULL
      } else{
        paste0(
          names(mSet$settings$subset)[i],
          "=",
          paste0(mSet$settings$subset[[i]], collapse = "+")
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
  print(mset_name)
  mset_name
}

#' @title Reset mSet to original post-normalization format
#' @description Returns mSet to right after normalization, without subsetting or changing the default variable.
#' @param mSet mSet object
#' @param fn filename to load in
#' @return mSet object
#' @rdname reset.mSet
#' @export 
reset.mSet <- function(mSet, fn) {
  origItem <-
    readRDS(file = fn)
  mSet$dataSet <- origItem$data
  mSet$analSet <- origItem$analysis
  mSet$settings <- origItem$settings
  mSet$report <- list(mzStarred = data.table::data.table(mz = colnames(mSet$dataSet$norm),
                                                         star = c(FALSE)))  
  data.table::setkey(mSet$report$mzStarred, mz)
  return(mSet)
}

#' @title Load mSet from mSet internal storage
#' @description MetShi mSets store previous dataset results and settings in the mSet storage. This loads one of those datasets.
#' @param mSet mSet object
#' @param name Name of subexperiemnt to load, Default: mSet$dataSet$cls.name
#' @return mSet object
#' @rdname load.mSet
#' @export 
load.mSet <- function(mSet, name = mSet$dataSet$cls.name) {
  mSet$analSet <- mSet$storage[[name]]$analysis
  mSet$settings <- mSet$storage[[name]]$settings
  mSet$report <- mSet$storage[[name]]$report
  return(mSet)
}

#' @title Store mSet analysis and settings in mSet internal storage
#' @description MetShi mSets store previous dataset results and settings in the mSet storage. This saves the current mSet in there.
#' @param mSet mSet object
#' @param name PARAM_DESCRIPTION, Default: mSet$dataSet$cls.name
#' @return Name of current subexperiment
#' @rdname store.mSet
#' @export 
store.mSet <- function(mSet, name = mSet$dataSet$cls.name) {
  mSet$storage[[name]]$analysis <- mSet$analSet
  mSet$storage[[name]]$settings <- mSet$settings
  mSet$storage[[name]]$report <- mSet$report
  return(mSet)
}

#' @title Change mSet
#' @description Wrapper function to change statistics variable(s) and mode (one/two-factor, time series)
#' @param mSet mSet object
#' @param stats_type Statistics category to change to
#' @param stats_var Metadata variable(s) to do statistics on, Default: NULL
#' @param time_var Metadata variable representing time, Default: NULL
#' @return mSet object
#' @seealso 
#'  \code{\link[shiny]{showNotification}}
#'  \code{\link[MetaboAnalystR]{SetDesignType}}
#' @rdname change.mSet
#' @export 
#' @importFrom shiny showNotification
#' @importFrom MetaboAnalystR SetDesignType
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
          metshiAlert(
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
          metshiAlert(
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


# subset by mz
subset_mSet_mz <- function(mSet, keep.mzs) {
  if (length(keep.mzs) > 0) {
    tables = c("norm", "proc")
    clss = c("cls", "proc.cls")
    combi.tbl = data.table::data.table(tbl = tables)
    
    for (i in 1:nrow(combi.tbl)){
      try({
        tbl = combi.tbl$tbl[i]
        # convert back
        if(any(grepl(colnames(mSet$dataSet[[tbl]]), pattern="/"))){
          eachPPM = T
          ppm = ceiling(as.numeric(gsub(colnames(tbl), pattern="^.*/", replacement="")))
          fixedCols = gsub(colnames(mSet$dataSet[[tbl]]), pattern="/.*$", replacement="")
        }else{
          eachPPM = F
          ppm = mSet$ppm
          fixedCols = gsub(colnames(mSet$dataSet[[tbl]]), pattern="RT.*$", replacement="")
        }
        keep = which(fixedCols %in% keep.mzs)
        mSet$dataSet[[tbl]] <- mSet$dataSet[[tbl]][, keep]
        mSet$settings$subset$mz <- keep.mzs
      }, silent = T)
    }
  }
  mSet
}

#' @title Subset mSet
#' @description Subset the current dataset based on metadata variables
#' @param mSet mSet object
#' @param subset_var Metadata variable to subset on
#' @param subset_group Within that variable, the categories to include in subset
#' @return mSet object
#' @rdname subset_mSet
#' @export 
#' @importFrom data.table data.table
subset_mSet <- function(mSet, subset_var, subset_group) {
  if (!is.null(subset_var)) {
    keep.i <- which(mSet$dataSet$covars[[subset_var]] %in% subset_group)
    keep.samples <- mSet$dataSet$covars$sample[keep.i]
    mSet$dataSet$covars <-
      mSet$dataSet$covars[sample %in% keep.samples, ]
    
    tables = c("norm", "proc")
    clss = c("cls", "proc.cls")
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

#' @title Pair mSet samples
#' @description Some analyses require paired samples. This function takes care of downsampling to achieve that.
#' @param mSet mSet object
#' @return mSet object
#' @seealso 
#'  \code{\link[data.table]{as.data.table}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[caret]{downSample}}
#' @rdname pair.mSet
#' @export 
#' @importFrom data.table as.data.table
#' @importFrom pbapply pbsapply
#' @importFrom caret downSample
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
    mSet <- subset_mSet(mSet,
                        subset_var = "sample",
                        subset_group = keep.samp)
    mSet$settings$paired <- TRUE
    mSet$dataSet$paired <- TRUE
  } else{
    metshiAlert("Not enough samples for paired analysis!")
    return(NULL)
  }
  mSet
}
