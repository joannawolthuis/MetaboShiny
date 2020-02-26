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
  keep.cols = colSums(is.na(metadata)) < nrow(metadata) & sapply(colnames(metadata), function(x) length(unique(metadata[,..x][[1]])) > 1) 
  metadata$sample <- gsub(metadata$sample, pattern = wipe.regex, replacement = "", perl=T)
  metadata.filt = unique(metadata[, ..keep.cols])
  if(!("individual" %in% colnames(metadata.filt))) metadata.filt$individual <- metadata.filt$sample
  #/METADATA
  
  print(pospath)
  print(negpath)
  print(metapath)
  
  poslist <- data.table::fread(pospath, header=T)
  neglist <- data.table::fread(negpath, header=T)
  
  hasRT=F
  
  # PIVOT IF WRONG SIDE AROUND - METABOLIGHTS DATA
  if(any(grepl("mass_to_charge", colnames(poslist)))){
    setnames(poslist, "mass_to_charge", "mzmed", skip_absent = T)
    setnames(neglist, "mass_to_charge", "mzmed", skip_absent = T)
    if(!is.na(poslist$retention_time[1])){
      hasRT = TRUE
      poslist$mzmed <- paste0(poslist$mzmed,"RT", poslist$retention_time)
      neglist$mzmed <- paste0(neglist$mzmed,"RT", neglist$retention_time)
    }
    rmcols = c("database_identifier", "chemical_formula", "smiles", "inchi", 
               "metabolite_identification", "fragmentation", "modifications", 
               "charge", "retention_time", "taxid", "species", "database", "database_version", 
               "reliability", "uri", "search_engine", "search_engine_score", 
               "smallmolecule_abundance_sub", "smallmolecule_abundance_stdev_sub", 
               "smallmolecule_abundance_std_error_sub")
    poslist_melty <- data.table::melt(poslist[,-..rmcols], id.vars=c("mzmed"), variable.name = "sample", value.name="into")  
    poslist <- data.table::dcast(poslist_melty, sample ~ mzmed, value.var = "into")
    neglist_melty <- data.table::melt(neglist[,-..rmcols], id.vars=c("mzmed"), variable.name = "sample", value.name="into")  
    neglist <- data.table::dcast(neglist_melty, sample ~ mzmed, value.var = "into")
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
    poslist[, sample := gsub(wipe.regex, replacement = "", sample, perl=T)]
    neglist[, sample := gsub(wipe.regex, replacement = "", sample, perl=T)]
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
  hasPPM = any(grepl(colnames(poslist), pattern = "/"))

  # ROUND MZ VALUES
  ismz <- suppressWarnings(which(!is.na(as.numeric(gsub(colnames(poslist), pattern="/.*$|RT.*$", replacement="")))))
  colnames(poslist)[ismz] <- pbapply::pbsapply(colnames(poslist)[ismz], function(mz){
      if(hasPPM | hasRT){
        split.mz <- stringr::str_split(mz, "/|RT")[[1]]
        mz = split.mz[1]
        ppm = split.mz[2]
      }
      ppmRange <- as.numeric(mz)/1e6 * as.numeric(ppm)
      zeros = sapply(ppmRange,zeros_after_period)
      decSpots = zeros + 1 # todo: verify this formula?
      roundedMz <- formatC(as.numeric(mz), digits = decSpots, format = "f")
      roundedPpm <-formatC(as.numeric(ppm), digits = 6, format = "f")
      newName = paste0(roundedMz, "+", if(hasPPM|hasRT) paste0(if(hasPPM)"/"else"RT", roundedPpm) else "") 
      return(newName)
    })
  ismz <- suppressWarnings(which(!is.na(as.numeric(gsub(colnames(neglist), pattern="/.*$|RT.*$", replacement="")))))
  colnames(neglist)[ismz] <- pbapply::pbsapply(colnames(neglist)[ismz], function(mz){
    if(hasPPM|hasRT){
      split.mz <- stringr::str_split(mz, "/|RT")[[1]]
      mz = split.mz[1]
      ppm = split.mz[2]
    }
    ppmRange <- as.numeric(mz)/1e6 * as.numeric(ppm)
    zeros = sapply(ppmRange,zeros_after_period)
    decSpots = zeros + 1 # todo: verify this formula?
    roundedMz <- formatC(as.numeric(mz), digits = decSpots, format = "f")
    roundedPpm <-formatC(as.numeric(ppm), digits = 6, format = "f")
    newName = paste0(roundedMz, "-",  if(hasPPM|hasRT) paste0(if(hasPPM)"/"else"RT", roundedPpm) else "") 
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
  mzlist <- merge(metadata.filt, mzlist, by = "sample")

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