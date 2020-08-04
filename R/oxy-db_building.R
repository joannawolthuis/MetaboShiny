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
                            wipe.regex = ".*_(?>POS|NEG)_[0+]*",
                            missperc.mz = 99,
                            missperc.samp = 100,
                            roundMz = T){
  
  ppm = as.numeric(ppm)

  metadata = NULL
  try({
    metadata <- data.table::fread(metapath)
    metadata <- reformat.metadata(metadata)  
    keep.cols = colSums(is.na(metadata)) < nrow(metadata) & sapply(colnames(metadata), function(x) length(unique(metadata[,..x][[1]])) > 1) 
    metadata$sample <- gsub(metadata$sample, pattern = wipe.regex, replacement = "", perl=T)
    metadata = unique(metadata[, ..keep.cols])
    if(!("individual" %in% colnames(metadata))) metadata$individual <- metadata$sample
  },silent = T)
  
  samplesIn <- metadata$sample

  zeros_after_period <- function(x) {
    if (isTRUE(all.equal(round(x),x))) return (0) # y would be -Inf for integer values
    y <- log10(abs(x)-floor(abs(x)))   
    ifelse(isTRUE(all.equal(round(y),y)), -y-1, -ceiling(y))}
  
  peaklists <- lapply(c("pos", "neg"), function(ionMode){
    
    peakpath = switch(ionMode,
                      "pos" = pospath,
                      "neg" = negpath)
    
    if(length(peakpath) == 0) return(list(ionMode = ionMode, 
                                          peaktbl = data.table::data.table()))
    
    bigFile = utils:::format.object_size(file.info(peakpath)$size, "GB")
    reallyBig = if(bigFile == "0 Gb") F else T
    nrows = length(count.fields(peakpath, sep = ","))
    
    print(paste("Importing", ionMode, "mode peaks!"))
    if(!reallyBig){
      peaklist <- data.table::fread(peakpath,
                                   header=T)
    }else{
      # === POSITIVE MODE ===
      con = file(peakpath, "r")
      # get missing count
      missRes = lapply(1:nrows, function(i, con){
        line = readLines(con, n = 1)
        splRow = stringr::str_split(line, ",")[[1]]
        sampName = splRow[1]
        sampName =  gsub(sampName, pattern = wipe.regex, replacement = "")
        splRow = splRow[3:length(splRow)]
        if(i > 1 & sampName %in% samplesIn){
          as.list(splRow == "0" | splRow == 0 | splRow == "" | is.na(splRow))
        }else{
          NULL
        }
      }, con = con)
      close.connection(con)
      missTable = data.table::rbindlist(missRes[!sapply(missRes, is.null)])
      miss_threshold_mz = ceiling(nrows * (missperc.mz/100))
      qualifies = which(colSums(missTable) <= miss_threshold_mz)
      missTable <- NULL
      missRes <- NULL
      gc()
      
      # get actual data
      con = file(peakpath, "r")
      cols= readLines(con, n = 1)
      filtRes = lapply(2:nrows, function(i, con, qualifies){
        line = readLines(con, n = 1)
        splRow = stringr::str_split(line, ",")[[1]]
        sampName = splRow[1]
        sampName =  gsub(sampName, pattern = wipe.regex, replacement = "", perl=T)
        label =splRow[2]
        splRow = splRow[3:length(splRow)]
        if(sampName %in% samplesIn){
          as.list(c(sampName, label, splRow[qualifies]))
        }else{
          NULL
        }
      }, con = con, qualifies = qualifies)
      close.connection(con)
      peaklist = data.table::rbindlist(filtRes[!sapply(filtRes, is.null)])
      colnames(peaklist) <- stringr::str_split(cols, ",")[[1]][c(1,2,qualifies)]
      filtRes <- NULL
      gc()
    }
    colnames(peaklist)[1:2] <- tolower(colnames(peaklist)[1:2])
    
    hasRT=F
    
    if(is.null(metadata)){
      if("label" %in% tolower(colnames(peaklists[[1]])[1:5])){
        premeta = peaklist[,1:30]
        colnames(premeta) <- tolower(colnames(premeta))
        metadata = premeta[,c("sample", "label")]
      }else{
        stop("Requires either metadata uploaded OR a label/Label column!")
      }
    }
    
    # PIVOT IF WRONG SIDE AROUND - METABOLIGHTS DATA
    if(any(grepl("mass_to_charge", colnames(peaklist)))){
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
    }else{
      hasRT=any(grepl(colnames(peaklist), pattern = "RT\\d+"))
    }
    
    if("label" %in% colnames(peaklist)[1:5]){
      peaklist[,label:=NULL]
    }
    
    # miss values
    if(!any(is.na(unlist(peaklist[1:5,10:20])))){
      peaklist[,(2:ncol(peaklist)) := lapply(.SD,function(x){ ifelse(x == 0, NA, x)}), .SDcols = 2:ncol(peaklist)]
    }
    
    if(missperc.samp < 100){
      miss_threshold_samp = ceiling(ncol(peaklist) * (missperc.samp/100))
      keep.samps.peak = which(rowSums(is.na(peaklist)) <= miss_threshold_samp)
      peaklist <- peaklist[,keep.samps.peak, with=F] # at least one sample with non-na
      try({
        shiny::showNotification(paste0("Remaining samples:", nrow(peaklist)))
      }, silent=T)
    }
    
    if(missperc.mz < 100 & !reallyBig){
      miss_threshold_mz = ceiling(nrow(peaklist)*(missperc.mz/100))
      peaklist <- peaklist[,which(colSums(is.na(peaklist)) <= miss_threshold_mz), with=F] # at least one sample with non-na
    }
    
    try({
      shiny::showNotification(paste0("Remaining m/z:", ncol(peaklist)))
    }, silent=T)
    
    # REGEX SAMPLE NAMES if applicable
    if(wipe.regex != ""){
      peaklist[, sample := gsub(wipe.regex, replacement = "", sample, perl=T)]
    } 
    
    peaklist$sample <- as.character(peaklist$sample)
    metadata$sample <- as.character(metadata$sample)
    unique.samples = unique(metadata$sample)
    
    # CHECK SAMPLES THAT ARE IN METADATA
    keep.samples.peak <- intersect(peaklist$sample, unique.samples)

    missing.samples = unique.samples[which(!(unique.samples %in% keep.samples.peak))]
    
    if(length(missing.samples) > 0 ){
      if(inshiny){
        shiny::showNotification(paste0("Missing samples: ", paste0(missing.samples, collapse = ","), ". Please check your peak table < > metadata for naming mismatches!"))
      }else{
        print(paste0("Missing samples: ", paste0(missing.samples, collapse = ","), ". Please check your peak table < > metadata for naming mismatches!"))
      }
    }
    
    # CHECK IF CUSTOM PPMVALUE
    hasPPM = any(grepl(colnames(peaklist), pattern = "/"))
    
    # ROUND MZ VALUES
    
    subbed = gsub(x = colnames(peaklist), pattern="/.*$|RT.*$", replacement="")
    ismz <- suppressWarnings(which(!is.na(as.numeric(subbed))))
    
    if(roundMz){
      colnames(peaklist)[ismz] <- sapply(colnames(peaklist)[ismz], function(mz){
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
    }
    
    peaklist <- peaklist[sample %in% keep.samples.peak,]
    
    # replace commas with dots
    if(any(grepl(unlist(peaklist[1:5,10:20]), pattern = ","))){
      for(j in seq_along(peaklist)){
        data.table::set(peaklist, i=which(grepl(pattern=",",peaklist[[j]])), j=j, value=gsub(",",".",peaklist[[j]]))
      }
    }
    
    return(list(ionMode = ionMode, peaktbl = peaklist))
  })
  
  isNonEmpty = unlist(sapply(peaklists, function(l) nrow(l$peaktbl) > 0))

  if(all(isNonEmpty)){
    tbl1 = peaklists[[1]]$peaktbl
    data.table::setkey(tbl1, sample)
    tbl2 = peaklists[[2]]$peaktbl
    data.table::setkey(tbl2, sample)  
    mzlist <- merge(tbl1, tbl2, by = "sample")
  }else{
    mzlist = peaklists[[which(isNonEmpty)]]$peaktbl
  }
  
  data.table::setkey(metadata, sample)
  
  mzlist <- merge(metadata, mzlist, by = "sample")

  mz.meta <- getColDistribution(mzlist)
  exp.vars = mz.meta$meta
  mzlist[,(exp.vars) := lapply(.SD,function(x){ ifelse(x == "" | is.na(x) | x == "Unknown",
                                                       "unknown",
                                                       x)}), .SDcols = exp.vars]
  
  colnames(mzlist)[which(colnames(mzlist) == "time")] <- "Time"
  mzlist$sample <- gsub("[^[:alnum:]./_-]", "", mzlist$sample)
  # - - - - - - - - - - - - - - - - - - 
  data.table::fwrite(mzlist, csvpath)

  hasPPM = any(grepl(colnames(mzlist), pattern = "/"))
  params = data.table::data.table(ppm=ppm,
                      ppmpermz = if(hasPPM) "yes" else "no")
  data.table::fwrite(params, gsub(csvpath, pattern="\\.csv", replacement="_params.csv"))
}
