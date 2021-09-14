getMissing <- function(peakpath, nrow=NULL){
  
  bigFile = utils:::format.object_size(file.info(peakpath)$size, "GB")
  reallyBig = if(bigFile == "0 Gb") F else T
  
  nrow = length(vroom::vroom_lines(peakpath)) - 1L
  
  print("Checking missing values...")
  
  skipCols=c("sample","label")
  
  con <- file(peakpath, "r")#, blocking = FALSE)
  header = readLines(con, n = 1) # empty
  cols = stringi::stri_split(header, fixed=",")[[1]]

  if(any(grepl("mass_to_charge", cols))){
    peaklist = data.table::fread(peakpath,fill=TRUE)
    data.table::setnames(peaklist, "mass_to_charge", "mzmed", skip_absent = T)
    if(!is.na(peaklist$retention_time[1])){
      hasRT = TRUE
      peaklist$mzmed <- paste0(peaklist$mzmed,"RT", peaklist$retention_time)
    }
    rmcols = c("database_identifier", "chemical_formula", "smiles", "inchi", 
               "metabolite_identification", "fragmentation", "modifications", 
               "charge", "retention_time", "taxid", "species", "database", "database_version", 
               "reliability", "uri", "search_engine", "search_engine_score", 
               "smallmolecule_abundance_sub", "smallmolecule_abundance_stdev_sub", 
               "smallmolecule_abundance_std_error_sub")
    peaklist_melty <- data.table::melt(peaklist[,-..rmcols], id.vars=c("mzmed"), variable.name = "sample", value.name="into")  
    peaklist <- data.table::dcast(peaklist_melty, sample ~ mzmed, value.var = "into")
    considerMe=which(!(tolower(colnames(peaklist)) %in% skipCols))
    peaklist = peaklist[, ..considerMe]
    totalMissing = colSums(peaklist == "0" | peaklist == 0 | peaklist == "" | is.na(peaklist))
  }else{
    if(!reallyBig){
      peaklist <- data.table::fread(peakpath,
                                    header=T,
                                    fill=TRUE)
      considerMe=which(!(tolower(colnames(peaklist)) %in% skipCols))
      peaklist = peaklist[, ..considerMe]
      totalMissing = colSums(peaklist == "0" | peaklist == 0 | peaklist == "" | is.na(peaklist))
    }else{
      considerMe=which(!(tolower(cols) %in% skipCols))
      mzs = cols[3:length(cols)]
      totalMissing <- rep(0, length(mzs))
      names(totalMissing) = mzs
      pbapply::pbsapply(2:nrow, function(i){
        line = readLines(con, n = 1) # empty
        splRow = stringr::str_split(line, pattern=",")[[1]]
        sampName = splRow[1]
        splRow = splRow[3:length(splRow)]
        isMissing = splRow == "0" | splRow == 0 | splRow == "" | is.na(splRow)
        totalMissing[isMissing] <<- totalMissing[isMissing] + 1
        NULL
      })
    }
  }
  list(missPerc = totalMissing, isMz = considerMe, nrows=nrow)
}

#' @title Merge metadata and peak tables into SQLITE database.
#' @description Wrapper function to integrate user data into the format used in MetaboShiny.
#' @param metapath Path to metadata csv
#' @param pospath Path to positive peaklist csv
#' @param negpath Path to negative peaklist csv
#' @param overwrite Overwrite existing database?, Default: FALSE
#' @param rtree Use RTree structure? Generally enables faster matching., Default: TRUE
#' @param ppm Parts per million error allowed., Default: 2
#' @param inshiny Running in shiny?, Default: F
#' @param csvpath Path to write resulting, Default: lcl$paths$csv_loc
#' @param wipe.regex Regex to apply to metadata sample names. Matches are removed from name, Default: '.*_(?>POS|NEG)_[0+]*'
#' @param missperc.mz Allowed percentage of samples missing for a  m/z, removed otherwise., Default: 99
#' @param missperc.samp Allowed percentage of m/z missing for a sample, removed otherwise., Default: 100
#' @seealso 
#'  \code{\link[data.table]{fread}},\code{\link[data.table]{melt.data.table}},\code{\link[data.table]{dcast.data.table}},\code{\link[data.table]{fwrite}}
#'  \code{\link[shiny]{showNotification}}
#'  \code{\link[stringr]{str_split}}
#' @rdname import.pat.csvs
#' @export 
#' @importFrom data.table fread melt dcast fwrite
#' @importFrom shiny showNotification
#' @importFrom stringr str_split
import.pat.csvs <- function(metapath,
                            pospath,
                            negpath,
                            overwrite = FALSE,
                            rtree = TRUE,
                            ppm = 2,
                            inshiny = F,
                            csvpath = lcl$paths$csv_loc,
                            wipe.regex = ".*_",
                            missperc.mz = 99,
                            missperc.samp = 100,
                            missList = c(pos=c(),neg=c()),
                            roundMz = T){
  ppm = as.numeric(ppm)

  metadata = NULL
  try({
    metadata <- data.table::fread(metapath,fill=TRUE)
    metadata <- reformat.metadata(metadata)
    keep.cols = colSums(is.na(metadata)) < nrow(metadata) & sapply(colnames(metadata), function(x) length(unique(metadata[,..x][[1]])) > 1) 
    metadata$sample <- gsub(metadata$sample, 
                            pattern = wipe.regex, 
                            replacement = "", 
                            perl=T)
    metadata = unique(metadata[, ..keep.cols])
    if(!("individual" %in% colnames(metadata))) metadata$individual <- metadata$sample
    meta_col_order = colnames(metadata)
    #meta_col_order = meta_col_order[meta_col_order != "sample"]
  },silent = T)
  
  if(is.null(metadata)){
    try({
      peaklist = data.table::fread(pospath,fill=TRUE)
    })
    try({
      peaklist = data.table::fread(negpath,fill=TRUE)
    })
    metadata = data.table::data.table(sample = peaklist$Sample,
                                      individual = peaklist$Sample,
                                      label = peaklist$Label)
    metadata$sample <- metadata$individual <- gsub(metadata$sample, 
                                                   pattern = wipe.regex, 
                                                   replacement = "", 
                                                   perl=T)
    meta_col_order <- c("sample", "individual", "label")
  }
  
  # check for all_unique or all_same columns and remove
  # check_meta_cols = setdiff(meta_col_order, c("sample", "individual"))
  # keep_meta_cols_bool = sapply(check_meta_cols, function(var){
  #   values = metadata[[var]]
  #   uniq_groups = table(values)
  #   is.uniq = length(uniq_groups) == nrow(metadata)
  #   is.same = length(uniq_groups) == 1
  #   !is.uniq & !is.same
  # })
  # 
  # keep_meta_cols = c("sample", "individual", names(which(keep_meta_cols_bool)))
  # metadata = metadata[, ..keep_meta_cols]
  
  samplesIn <- metadata$sample
  if(any(duplicated(samplesIn))){
    warning("'sample' column must be unique! Last instance will be used. Please correct to avoid unpredictable results.")
    metadata = metadata[!duplicated(sample)]
  }

  zeros_after_period <- function(x) {
    if (isTRUE(all.equal(round(x),x))) return (0) # y would be -Inf for integer values
    y <- log10(abs(x)-floor(abs(x)))   
    ifelse(isTRUE(all.equal(round(y),y)), -y-1, -ceiling(y))}
  
  peaklists = lapply(c("pos", "neg"), function(ionMode){
    
    peakpath = switch(ionMode,
                      "pos" = pospath,
                      "neg" = negpath)
    
    if(length(peakpath) == 0){
      return( list(ionMode = ionMode, 
                   peaktbl = data.table::data.table()))
    }
    
    print(paste("Importing", ionMode, "mode peaks!"))
    
    if(length(missList[[ionMode]]) > 0){
      missCounts = missList[[ionMode]]
    }else{
      missCounts = getMissing(peakpath)
    }
    
    bigFile = utils:::format.object_size(file.info(peakpath)$size, "GB")
    reallyBig = if(bigFile == "0 Gb") F else T

    nrows = length(vroom::vroom_lines(peakpath, altrep = TRUE, progress = TRUE)) - 1L
    miss_threshold_mz = ceiling(nrows * (missperc.mz/100))
    qualifies = which(missCounts$missPerc <= miss_threshold_mz)
    # missCounts = NULL
    # gc()
    
    # get actual data
    print("Subsetting table...")
    write_loc = tempfile()
    if(file.exists(write_loc)) file.remove(write_loc)
    con_read = base::file(peakpath, "r")
    cols = stringr::str_split(readLines(con_read, n = 1), ",|\\t")[[1]]
    
    # PIVOT IF WRONG SIDE AROUND - METABOLIGHTS DATA
    if(any(grepl("mass_to_charge", cols))){
      peaklist = data.table::fread(peakpath,fill=TRUE)
      data.table::setnames(peaklist, "mass_to_charge", "mzmed", skip_absent = T)
      if(!is.na(peaklist$retention_time[1])){
        hasRT = TRUE
        peaklist$mzmed <- paste0(peaklist$mzmed,"RT", peaklist$retention_time)
      }
      rmcols = c("database_identifier", "chemical_formula", "smiles", "inchi", 
                 "metabolite_identification", "fragmentation", "modifications", 
                 "charge", "retention_time", "taxid", "species", "database", "database_version", 
                 "reliability", "uri", "search_engine", "search_engine_score", 
                 "smallmolecule_abundance_sub", "smallmolecule_abundance_stdev_sub", 
                 "smallmolecule_abundance_std_error_sub")
      peaklist_melty <- data.table::melt(peaklist[,-..rmcols], id.vars=c("mzmed"), variable.name = "sample", value.name="into")  
      peaklist <- data.table::dcast(peaklist_melty, sample ~ mzmed, value.var = "into")
      peaklist$sample <- gsub(pattern = wipe.regex, "", peaklist$sample,perl = T)
      meta_row_order = match(metadata$sample, peaklist$sample)
      peaklist <- cbind(metadata[meta_row_order,..meta_col_order], 
                        peaklist[,-c("sample", "individual")])
    }else{
      row = c(tolower(c(colnames(metadata[,..meta_col_order]))),
                        cols[missCounts$isMz][qualifies])
      
      write(paste0(row, collapse=","),
            file = write_loc,
            append=TRUE)
      
      pbapply::pbsapply(1:nrows, function(i, con, qualifies){
        line = readLines(con_read, n = 1)
        splRow = stringr::str_split(line, ",")[[1]]
        sampName = splRow[1]
        sampName = gsub(sampName, 
                         pattern = wipe.regex, 
                         replacement = "", perl=T)
        label = splRow[2]
        splRow = splRow[missCounts$isMz]
        splRow[splRow == "0" | splRow == 0 | splRow == ""] <- NA
        try({
          row = c(metadata[sample == as.character(sampName), ..meta_col_order],
                  splRow[qualifies])
          write(paste0(row, collapse=","),
                file = write_loc,append=TRUE)
        })
      }, con = con_read, qualifies = qualifies)
      
      close.connection(con_read)
      
      peaklist = data.table::fread(file = write_loc, 
                                   sep = ",",
                                   #key = "sample",
                                   header = T)
    }
    
    if(missperc.samp < 100){
      miss_threshold_samp = ceiling((ncol(peaklist)-1) * (missperc.samp/100))
      rs = rowSums(is.na(peaklist))
      keep.samps.peak = which(rs <= miss_threshold_samp)
      if(length(keep.samps.peak) == 0){
        lowerMe = seq(missperc.samp, 100, 5)
        for(missPerc in lowerMe){
          miss_threshold_samp = ceiling(ncol(peaklist) * (missPerc/100))
          keep.samps.peak = which(rs <= miss_threshold_samp)
          if(length(keep.samps.peak) > 1){
            try({
              shiny::showNotification(paste0("No samples qualify with original missing m/z threshold! Changed to ", missPerc, " percent"))
            }, silent = T)
            break
          }
        }
      }
      peaklist <- peaklist[keep.samps.peak, ] # at least one sample with non-na
      gc()
      try({
        shiny::showNotification(paste0("Remaining samples:", nrow(peaklist)))
      }, silent=T)
    }
    
    if(missperc.mz < 100 & !reallyBig){
      miss_threshold_mz = ceiling(nrow(peaklist)*(missperc.mz/100))
      peaklist <- peaklist[,which(colSums(is.na(peaklist)) <= miss_threshold_mz), with=F] # at least one sample with non-na
    }
    
    try({
      msg=paste0("Remaining m/z:", ncol(peaklist))
      print(msg)
      shiny::showNotification(msg)
    }, silent=T)
    
    # CHECK SAMPLES THAT ARE IN METADATA
    keep.samples.peak <- if(typeof(metadata) == "list") intersect(peaklist$sample, samplesIn) else peaklist$sample
    missing.samples = if(typeof(metadata) == "list") samplesIn[which(!(samplesIn %in% keep.samples.peak))] else c()

    if(length(missing.samples) > 0 ){
      if(inshiny){
        shiny::showNotification(paste0("Missing samples: ", paste0(missing.samples, collapse = ","), ". Please check your peak table < > metadata for naming mismatches!"))
      }else{
        print(paste0("Missing samples: ", paste0(missing.samples, collapse = ","), ". Please check your peak table < > metadata for naming mismatches!"))
      }
    }
    
    hasRT = any(grepl(colnames(peaklist), pattern = "\\d+RT\\d+"))
    # CHECK IF CUSTOM PPMVALUE
    hasPPM = any(grepl(colnames(peaklist), pattern = "/"))
    
    # ROUND MZ VALUES
    
    subbed = gsub(x = colnames(peaklist), pattern="/.*$|RT.*$", replacement="")
    ismz <- suppressWarnings(which(!is.na(as.numeric(subbed))))
    
    if(roundMz){
      print("Rounding m/z to appropriate decimal for given ppm value...")
      colnames(peaklist)[ismz] <- pbapply::pbsapply(colnames(peaklist)[ismz], function(mz){
        if(hasPPM | hasRT){
          split.mz <- stringr::str_split(mz, "/|RT")[[1]]
          mz = split.mz[1]
          ppm = split.mz[2]
        }
        ppmRange <- as.numeric(mz)/1e6 * as.numeric(ppm)
        zeros = sapply(ppmRange,zeros_after_period)
        decSpots = zeros + 2 # todo: verify this formula?
        roundedMz <- formatC(as.numeric(mz), digits = decSpots, format = "f")
        roundedPpm <-formatC(as.numeric(ppm), digits = 6, format = "f")
        newName = paste0(roundedMz, 
                         if(ionMode=="pos") "+" else "-",
                         if(hasPPM|hasRT) paste0(if(hasPPM) "/" else "RT", roundedPpm) else "") 
        return(newName)
      })
    }else{
      colnames(peaklist)[ismz] <- paste0(colnames(peaklist)[ismz], ifelse(ionMode == "pos", "+", "-"))
    }
    
    peaklist <- peaklist[sample %in% keep.samples.peak,]
    print(colnames(peaklist)[1:30])
    gc()
    
    # replace commas with dots
    # if(any(grepl(unlist(peaklist[1:5,10:20]), pattern = ","))){
    #   for(j in seq_along(peaklist)){
    #     data.table::set(peaklist, i=which(grepl(pattern=",",peaklist[[j]])), j=j, value=gsub(",",".",peaklist[[j]]))
    #   }
    # }
    return(list(ionMode = ionMode, peaktbl = peaklist))
  })
  
  isNonEmpty = unlist(sapply(peaklists, function(l) nrow(l$peaktbl) > 0))

  if(all(isNonEmpty)){
    tbl1 = peaklists[[1]]$peaktbl
    tbl2 = peaklists[[2]]$peaktbl
    tbl.meta <- getColDistribution(tbl1)
    tbl.exp.vars = tbl.meta$meta
    mzlist <- merge(tbl1, tbl2, by = colnames(tbl1)[tbl.exp.vars])
  }else{
    mzlist = peaklists[[which(isNonEmpty)]]$peaktbl
  }
  
  hasRT=any(grepl(colnames(mzlist), pattern = "RT\\d+"))

  mz.meta <- getColDistribution(mzlist)
  exp.vars = mz.meta$meta
  mzlist[,(exp.vars) := lapply(.SD,function(x){ ifelse(x == "" | is.na(x) | x == "Unknown",
                                                       "unknown",
                                                       x)}), .SDcols = exp.vars]
  
  colnames(mzlist)[exp.vars] <- gsub("\\.x$", "", colnames(mzlist)[exp.vars])
  colnames(mzlist)[which(colnames(mzlist) == "time")] <- "Time"
  mzlist$sample <- gsub("[^[:alnum:]./_-]", "", mzlist$sample)
  
  print("Preview:")
  print(mzlist[1:2,1:40])
  # - - - - - - - - - - - - - - - - - - 
  data.table::fwrite(mzlist, csvpath)

  hasPPM = any(grepl(colnames(mzlist), pattern = "/"))
  params = data.table::data.table(ppm=ppm,
                      ppmpermz = if(hasPPM) "yes" else "no")
  data.table::fwrite(params, gsub(csvpath, pattern="\\.csv", replacement="_params.csv"))
}
