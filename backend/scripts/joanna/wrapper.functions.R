#' @export
wrapper.makedb <- function(source.db){
  
}

#' @export
wrapper.identify <- function(db.loc,
                             outlist.loc,
                             sources=c("internal", "hmdb", "chebi"),
                             isofilt=TRUE,
                             excl.adducts=c("PLACEHOLDER")){
  print(paste("Matching outlist peaks with", paste(sources, collapse=","), "databases...", sep=" "))
  # --- init ---
  identified.list <- pblapply(sources, FUN=function(source){
    results <- iden.code.binned(outlist.loc, 
                     file.path(db.loc, paste(source, ".full.db", sep="")), 
                     isofilt=isofilt,
                     excl.adducts=excl.adducts)
    results
  })
  # I would like a summed table too with all results (for later csv creation)
  # --- get results ---
  identified <- as.data.table(rbindlist(identified.list))
  # --- return ---
  unique(identified)
}

#' @export
get.csv <- function(patdb, 
                    time.series = F, 
                    group.by.adduct = F, 
                    exp.condition = "diet",
                    chosen.display = "mz",
                    relative.time = F,
                    match.table = "hmdb",
                    max.vals = -1,
                    outpath = "results.csv"){
  # --- announce some stuff ---
  max.cols <- if(max.vals == -1) "unlimited" else max.vals + 3
  cat(fn$paste("Creating csv for metabolomics analysis with max $max.cols columns.
 Will find $chosen.display values from matches found in the '$match.table' database.
 Chosen experimental conditon is '$exp.condition'.\n"))
  if(time.series) cat("You indicated that the dataset is a time series.\n")
  # --- road split ---
  if(group.by.adduct){
    cat("Will sum adduct peak intensities for each molecular formula.\n")
    group.line <- "group by i.filename, m.baseformula"
  } else{
    cat("Will not group adducts and display for each seperate peak.\n")
    group.line = "group by i.filename, i.mz"
  }
  if(chosen.display == "mz"){
    tab <- "i."
    filter.line = ""
    match.line = ""
  } else{
    tab <- "m."
    filter.line = "and not m.baseformula ISNULL"
    match.line = fn$paste("join matches_$match.table m on m.[mzmed.pgrp] = i.mz")
  }
  id.line <- fn$paste("group_concat(distinct $tab$chosen.display) as identifier")
  ## TODO TOMORROW : for per mz concatenate name /baseformula and adduct (ie. A[M+H], B[M-H]), causing the 'aggregate fun not found' error#
  # --- connect to patient db containing match values ---
  conn <- dbConnect(RSQLite::SQLite(), patdb)
  #dbExecute(conn, "drop table if exists avg_intensities")
  # --- create table with averaged intensities for the triplicates ---
  make.query <- strwrap(fn$paste("create table if not exists avg_intensities as
                                 select [mzmed.pgrp] as mz, filename, avg(intensity) as intensity 
                                 from mzintensities 
                                 group by filename, [mzmed.pgrp]"), width=10000, simplify=TRUE)
  dbExecute(conn, make.query)
  index.query <- "create index if not exists avg_index on avg_intensities(filename, mz)"
  # --- build result fetching query ---
  get.query = strwrap(fn$paste("select distinct i.filename, d.animal_internal_id, s.$exp.condition as label, d.sampling_date, $id.line, sum(i.intensity) as intensity
                               from avg_intensities i
                               join individual_data d
                               on i.filename = d.card_id
                               join setup s
                               on d.[group] = s.[group]
                               $match.line $filter.line $group.line"), width=10000, simplify=TRUE)
  print(get.query)
  
  # --- fetch results ---
  z = dbGetQuery(conn, get.query)
  z.dt <- as.data.table(z)
  # --- cast to right format ---
  cast.dt <- dcast.data.table(z.dt, animal_internal_id + sampling_date + label ~ identifier, fun.aggregate = sum, value.var = "intensity") # what to do w/ duplicates? 
  max.vals.final <- if(max.vals == -1) ncol(cast.dt) else max.vals + 3
  small.set <- cast.dt[,1:(max.vals.final)]
  # --- name for metaboanalyst ---
  names(small.set)[1:3] <- c("Sample", "Time", "Label")
  # --- make time series if necessary (this factorizes sampling date) ---
  small.set <- if(!time.series) small.set[,-"Time", with=F] else small.set
  if(time.series){
    small.set$Time <- as.numeric(as.factor(as.Date(small.set$Time)))
    small.set$Sample <-  paste(gsub(small.set$Sample,
                                    pattern="\\.\\d",
                                    replacement=""), 
                               small.set$Time, sep="T") # give each time an unique name
  }
  # --- measure file size ---
  size <- object.size(small.set)
  cat("Resulting file will be approximately ")
  print(size, units = "Mb")
  # ------- return -------
  return(small.set)
}