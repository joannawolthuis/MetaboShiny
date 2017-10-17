#' @export
get.csv <- function(patdb, 
                    time.series = F, 
                    exp.condition = "diet",
                    max.vals = -1,
                    group_adducts = F,
                    which_dbs = NA,
                    which_adducts = NA){
  library(reshape2)
  library(DBI)
  library(RSQLite)
  library(data.table)
  # --- announce some stuff ---
  max.cols <- if(max.vals == -1) "unlimited" else max.vals + 3
  cat(fn$paste("Creating csv for metabolomics analysis with max $max.cols columns.
 Chosen experimental conditon is '$exp.condition'.\n"))
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
  z = if(group_adducts){
    get_all_matches(exp.condition, 
                    pat.conn = conn,
                    which_dbs,
                    which_adducts)
  }else{
    dbGetQuery(conn, strwrap(fn$paste("select distinct i.filename, d.animal_internal_id, s.[$exp.condition] as label, d.sampling_date, i.mz as identifier, sum(i.intensity) as intensity
                      from avg_intensities i
                      join individual_data d
                      on i.filename = d.card_id
                      join setup s
                      on d.[group] = s.[group]
                      group by d.animal_internal_id, d.sampling_date, i.mz"), 
            width=10000, 
            simplify=TRUE)
            )
      }
  # --- fetch results ---
  z.dt <- as.data.table(z)
  # --- cast to right format ---
  cast.dt <- dcast.data.table(z.dt, animal_internal_id + sampling_date + label ~ identifier, fun.aggregate = sum, value.var = "intensity") # what to do w/ duplicates? 
  print(cast.dt)
  max.vals.final <- if(max.vals == -1) ncol(cast.dt) else max.vals + 3
  small.set <- cast.dt[,1:(max.vals.final)]
  # --- name for metaboanalyst ---
  names(small.set)[1:3] <- c("Sample", "Time", "Label")
  # --- make time series if necessary (this factorizes sampling date) ---
  small.set <- if(!time.series) small.set[,-"Time", with=F] else small.set
  if(time.series){
    small.set$Time <- as.numeric(as.factor(as.Date(small.set$Time)))
  }
  # --- measure file size ---
  size <- object.size(small.set)
  cat("Resulting file will be approximately ")
  print(size, units = "MB")
  # ------- return -------
  return(small.set)
}
