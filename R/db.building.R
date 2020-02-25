#' @export
import.pat.csvs <- function(db.name,
                            metapath,
                            pospath,
                            negpath,
                            overwrite=FALSE,
                            rtree=TRUE,
                            make.full = TRUE,
                            ppm=2,
                            inshiny=F,
                            csvpath=lcl$paths$csv_loc,
                            wipe.regex = ".*_(?>POS|NEG)_[0+]*",
                            missperc=99){
  ppm = as.numeric(ppm)

  # METADATA
  metadata <- data.table::fread(metapath)
  metadata <- MetaboShiny::reformat.metadata(metadata)
  keep.cols = sapply(colnames(metadata), function(x) length(unique(metadata[,..x][[1]])) > 1) 
  metadata$sample <- gsub(metadata$sample, pattern = wipe.regex, replacement = "", perl=T)
  metadata.filt = unique(metadata[, ..keep.cols])
  if(!("individual" %in% colnames(metadata.filt))) metadata.filt$individual <- metadata.filt$sample
  #/METADATA
  
  print(pospath)
  print(negpath)
  print(metapath)
  
  poslist <- data.table::fread(pospath, header=T)
  neglist <- data.table::fread(negpath, header=T)
  
  # PIVOT IF WRONG SIDE AROUND - METABOLIGHTS DATA
  if(any(grepl("mzmed|mass_to_charge", colnames(poslist)))){
    setnames(poslist, "mass_to_charge", "mzmed", skip_absent = T)
    setnames(neglist, "mass_to_charge", "mzmed", skip_absent = T)
    poslist_melty <- data.table::melt(poslist_melty, id.vars=c("sample"), variable.name = "mzmed", value.name="into")  
    poslist <- data.table::dcast(poslist, mzmed ~ sample, value.var = "into")
    neglist_melty <- data.table::melt(neglist_melty, id.vars=c("sample"), variable.name = "mzmed", value.name="into")  
    neglist <- data.table::dcast(neglist, mzmed ~ sample, value.var = "into")
    # REQUIREMENTS: "sample" - mz1 mz2 mz3 mz3 ... 
    # TODO: check for metabolights samples!
  }
  
  if("label" %in% colnames(poslist)[1:5]){
    poslist[,label:=NULL]
    neglist[,label:=NULL]
  }

  # replace 0 with NA  
  if(!any(is.na(unlist(poslist[1:5,10:20])))){
    poslist[,(2:ncol(poslist)) := lapply(.SD,function(x){ ifelse(x == 0, NA, x)}), .SDcols = 2:ncol(poslist)]
    neglist[,(2:ncol(neglist)) := lapply(.SD,function(x){ ifelse(x == 0, NA, x)}), .SDcols = 2:ncol(neglist)]
  }
  
  # miss values
  miss_threshold = ceiling(missperc/100 * nrow(poslist))

  poslist <- poslist[,which(colSums(is.na(poslist)) <= miss_threshold),with=F] # at least one sample with non-na
  neglist <- neglist[,which(colSums(is.na(neglist)) <= miss_threshold),with=F] # at least one sample with non-na
  
  print(paste0("Remaining m/z values for positive mode:", ncol(poslist)))
  print(paste0("Remaining m/z values for negative mode:", ncol(neglist)))
  
  # REGEX SAMPLE NAMES if applicable
  if(wipe.regex != ""){
    poslist[, sample := gsub(wipe.regex, replacement = "", sample)]
    neglist[, sample := gsub(wipe.regex, replacement = "", sample)]
  } 
  
  # REMOVE EXTRA DECIMALS THAT ARE NOT WITHIN PPM MARGIN
  zeros_after_period <- function(x) {
    if (isTRUE(all.equal(round(x),x))) return (0) # y would be -Inf for integer values
    y <- log10(abs(x)-floor(abs(x)))   
    ifelse(isTRUE(all.equal(round(y),y)), -y-1, -ceiling(y))}
  
  
  poslist$sample <- as.character(poslist$sample)
  neglist$sample <- as.character(neglist$sample)
  metadata.filt$sample <- as.character(metadata.filt$sample)
  
  unique.samples = unique(metadata.filt$sample)
  
  # CHECK SAMPLES THAT ARE IN METADATA
  keep.samples.pos <- intersect(poslist$sample, unique.samples)
  keep.samples.neg <- intersect(neglist$sample, unique.samples)
  missing.samples = setdiff(unique.samples, unique(c(keep.samples.pos, keep.samples.neg)))
  if(length(missing.samples) > 0 ){
    if(inshiny){
      shiny::showNotification(paste0("Missing samples: ", paste0(missing.samples, collapse = ","), ". Please check your peak table < > metadata for naming mismatches!"))
    }else{
      print(paste0("Missing samples: ", paste0(missing.samples, collapse = ","), ". Please check your peak table < > metadata for naming mismatches!"))
    }
  }
  
  # CHECK IF CUSTOM PPMVALUE
  hasPPM = any(grepl(colnames(poslist), pattern = "\\|"))

  # ROUND MZ VALUES
  ismz <- suppressWarnings(which(!is.na(as.numeric(gsub(colnames(poslist), pattern="\\|.*$", replacement="")))))
  colnames(poslist)[ismz] <- pbapply::pbsapply(colnames(poslist)[ismz], function(mz){
      if(hasPPM){
        split.mz <- stringr::str_split(mz, "/")[[1]]
        mz = split.mz[1]
        ppm = split.mz[2]
      }
      ppmRange <- as.numeric(mz)/1e6 * as.numeric(ppm)
      zeros = sapply(ppmRange,zeros_after_period)
      decSpots = zeros + 1 # todo: verify this formula?
      roundedMz <- formatC(as.numeric(mz), digits = decSpots, format = "f")
      roundedPpm <-formatC(as.numeric(ppm), digits = 6, format = "f")
      newName = paste0(roundedMz, "+", if(hasPPM) paste0("/", roundedPpm) else "") 
      return(newName)
    })
  ismz <- suppressWarnings(which(!is.na(as.numeric(gsub(colnames(neglist), pattern="\\|.*$", replacement="")))))
  colnames(neglist)[ismz] <- pbapply::pbsapply(colnames(neglist)[ismz], function(mz){
    if(hasPPM){
      split.mz <- stringr::str_split(mz, "/")[[1]]
      mz = split.mz[1]
      ppm = split.mz[2]
    }
    ppmRange <- as.numeric(mz)/1e6 * as.numeric(ppm)
    zeros = sapply(ppmRange,zeros_after_period)
    decSpots = zeros + 1 # todo: verify this formula?
    roundedMz <- formatC(as.numeric(mz), digits = decSpots, format = "f")
    roundedPpm <-formatC(as.numeric(ppm), digits = 6, format = "f")
    newName = paste0(roundedMz, "-", if(hasPPM) paste0("/", roundedPpm) else "") 
    return(newName)
    })
  
  poslist <- poslist[sample %in% keep.samples.pos,]
  neglist <- neglist[sample %in% keep.samples.neg,]
  
  # replace commas with dots
  if(any(grepl(unlist(poslist[1:5,10:20]), pattern = ","))){
    for(j in seq_along(poslist)){
      set(poslist, i=which(grepl(pattern=",",poslist[[j]])), j=j, value=gsub(",",".",poslist[[j]]))
    }
    for(j in seq_along(neglist)){
      set(neglist, i=which(grepl(pattern=",",neglist[[j]])), j=j, value=gsub(",",".",neglist[[j]]))
    }
  }
    
  # merge metadata with peaks
  setkey(poslist, sample)
  setkey(neglist, sample)
  setkey(metadata, sample)
  
  mzlist <- merge(poslist, neglist, by = "sample")
  mzlist <- merge(metadata, mzlist, by = "sample")

  mz.meta <- getColDistribution(mzlist)
  exp.vars = mz.meta$meta
  mzlist[,(exp.vars) := lapply(.SD,function(x){ ifelse(x == "" | is.na(x) | x == "Unknown", "unknown", x)}), .SDcols = exp.vars]
  
  colnames(mzlist)[which(colnames(mzlist) == "time")] <- "Time"
  mzlist$sample <- gsub("[^[:alnum:]./_-]", "", mzlist$sample)
  # - - - - - - - - - - - - - - - - - - 
  data.table::fwrite(mzlist, csvpath)
  params = data.table(ppm=ppm,
                      ppmpermz = if(hasPPM) "yes" else "no")
  data.table::fwrite(params, gsub(csvpath, pattern="\\.csv", replacement="_params.csv"))
}

#' @export
build.pat.db <- function(db.name,
                         metapath,
                         pospath,
                         negpath,
                         overwrite=FALSE,
                         rtree=TRUE,
                         make.full = TRUE,
                         ppm=2,
                         inshiny=F,
                         wipe.regex = ".*_(?>POS|NEG)_[0+]*"){

  ppm = as.numeric(ppm)
  
  # =================
  
  if(overwrite==TRUE & file.exists(db.name)) file.remove(db.name)
  
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db.name)
  
  metadata <- data.table::fread(metapath)
  
  metadata <- MetaboShiny::reformat.metadata(metadata)

  keep.cols = sapply(colnames(metadata), function(x) length(unique(metadata[,..x][[1]])) > 1)  
  
  metadata$sample <- gsub(metadata$sample, pattern = wipe.regex, replacement = "", perl=T)
  
  metadata.filt = unique(metadata[, ..keep.cols])
  
  if(!("individual" %in% colnames(metadata.filt))) metadata.filt$individual <- metadata.filt$sample
  
  RSQLite::dbWriteTable(conn, "individual_data", metadata.filt, overwrite=TRUE) # insert into
  
  if("sample" %in% colnames(poslist)){
    # reshape
    poslist_melty <- data.table::melt(poslist[,-"label"], id.vars=c("sample"), variable.name = "mzmed", value.name="into")  
    neglist_melty <- data.table::melt(neglist[,-"label"], id.vars=c("sample"),variable.name = "mzmed", value.name="into")
    poslist <- data.table::dcast(poslist_melty, mzmed ~ sample, value.var = "into")
    neglist <- data.table::dcast(neglist_melty, mzmed ~ sample, value.var = "into")
  }
  
  adjust = which(!(colnames(poslist) %in% c("mzmed", "mass_to_charge")))
  colnames(poslist)[adjust] <- gsub(colnames(poslist)[adjust], pattern = wipe.regex, replacement = "", perl=T)
  colnames(neglist)[adjust] <- gsub(colnames(neglist)[adjust], pattern = wipe.regex, replacement = "", perl=T)
  
  zeros_after_period <- function(x) {
    if (isTRUE(all.equal(round(x),x))) return (0) # y would be -Inf for integer values
    y <- log10(abs(x)-floor(abs(x)))   
    ifelse(isTRUE(all.equal(round(y),y)), -y-1, -ceiling(y))} # corrects case ending with ..01
  
  unique.samples = unique(metadata.filt$sample)
  
  keep.samples.pos <- intersect(colnames(poslist), unique.samples)
  keep.samples.neg <- intersect(colnames(neglist), unique.samples)
  
  missing.samples = setdiff(unique.samples, unique(c(keep.samples.pos, keep.samples.neg)))
  if(length(missing.samples) > 0 ){
    if(inshiny){
      shiny::showNotification(paste0("Missing samples: ", paste0(missing.samples, collapse = ","), ". Please check your peak table < > metadata for naming mismatches!"))
    }else{
      print(paste0("Missing samples: ", paste0(missing.samples, collapse = ","), ". Please check your peak table < > metadata for naming mismatches!"))
    }
  }
  
  setnames(poslist, "mass_to_charge", "mzmed", skip_absent = T)
  setnames(neglist, "mass_to_charge", "mzmed", skip_absent = T)
  
  # - - fix errors in rounding m/z - -
  poslist$mzmed <- pbapply::pbsapply(poslist$mzmed, function(mz){
    ppmRange <- as.numeric(mz)/1e6 * as.numeric(ppm)
    zeros = sapply(ppmRange,zeros_after_period)
    decSpots = zeros + 1 # todo: verify this formula?
    roundedMz <- formatC(as.numeric(mz), digits = decSpots, format = "f")
    return(paste0(roundedMz, "+"))
  })
  
  neglist$mzmed <- pbapply::pbsapply(neglist$mzmed, function(mz){
    ppmRange <- as.numeric(mz)/1e6 * as.numeric(ppm)
    zeros = sapply(ppmRange,zeros_after_period)
    decSpots = zeros + 1 # todo: verify this formula?
    roundedMz <- formatC(as.numeric(mz), digits = decSpots, format = "f")
    return(paste0(roundedMz, "-"))
  })
  
  keepcols.pos <- c("mzmed", keep.samples.pos)
  keepcols.neg <- c("mzmed", keep.samples.neg)
  
  poslist <- poslist[, ..keepcols.pos]
  neglist <- neglist[, ..keepcols.neg]
  
  # replace commas with dots
  poslist <- data.table::as.data.table(sapply(poslist, gsub, pattern = ",", replacement= "."))
  neglist <- data.table::as.data.table(sapply(neglist, gsub, pattern = ",", replacement= "."))
  
 if(inshiny) setProgress(.20)
  
  mzvals <- data.table::data.table(mzmed = c(poslist$mzmed, neglist$mzmed),
                                   foundinmode = c(rep("positive", nrow(poslist)), rep("negative", nrow(neglist))))
  
  if(inshiny) setProgress(.30)
  
  poslist <- data.table::melt(poslist,#[,(rmv.cols) := NULL],
                              id.vars="mzmed",
                              variable.name="filename",
                              value.name="intensity",
                              variable.factor=TRUE
  )
  neglist <- data.table::melt(neglist,#[,(rmv.cols) := NULL],
                              id.vars="mzmed",
                              variable.name="filename",
                              value.name="intensity",
                              variable.factor=TRUE
  )
  if(inshiny) setProgress(.40)
  
  sql.make.int <- strwrap("CREATE TABLE mzintensities(
                          ID INTEGER PRIMARY KEY AUTOINCREMENT,
                          mzmed decimal(30,13),
                          filename text,
                          intensity float)", width=10000, simplify=TRUE)
  
  RSQLite::dbExecute(conn, sql.make.int)
  
  if(inshiny) setProgress(.60)
  
  # --- write intensities to table and index ---
  RSQLite::dbWriteTable(conn, "mzintensities", poslist, append=TRUE) # insert into
  RSQLite::dbWriteTable(conn, "mzintensities", neglist, append=TRUE) # insert into
  
  RSQLite::dbExecute(conn, "CREATE INDEX intindex ON mzintensities(filename,'mzmed',intensity)")
  RSQLite::dbExecute(conn, "CREATE INDEX intindex2 ON mzintensities('mzmed')")
  
  # ------------------------
  
  sql.make.meta <- strwrap("CREATE TABLE mzvals(
                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
                           mzmed decimal(30,13),
                           foundinmode text)", width=10000, simplify=TRUE)
  RSQLite::dbExecute(conn, sql.make.meta)
  RSQLite::dbExecute(conn, "create index mzfind on mzvals(mzmed, foundinmode);")
  
  if(inshiny) setProgress(.70)
  
  RSQLite::dbWriteTable(conn, "params", 
                        data.table::data.table(ppm=ppm), 
                        overwrite=T)
  
  # --- write vals to table ---
  RSQLite::dbWriteTable(conn, "mzvals", mzvals, append=TRUE) # insert into
  
  # --- make range table (choose if R*tree or not) ---
 
  if(inshiny) setProgress(.80)
  
  # --- cleanup ---
  RSQLite::dbExecute(conn, "VACUUM")
  
  if(inshiny) setProgress(.90)
  
  # ----------------
  RSQLite::dbDisconnect(conn)
}

