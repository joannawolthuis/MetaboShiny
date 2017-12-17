#' @export
get.csv <- function(patdb, 
                    time.series = T, 
                    exp.condition = "diet",
                    max.vals = -1,
                    group_adducts = T,
                    which_dbs = file.path(options$db_dir, "kegg.full.db"),
                    which_adducts = c("M+H", "M-H", "M"),
                    group_by = "mz",
                    var_table = "setup"){
  library(data.table)
  # --- announce some stuff ---
  adducts <- paste(which_adducts, collapse = ", ")
  groupfac <- group_by
  max.cols <- if(max.vals == -1) "unlimited" else max.vals + 3
  cat(fn$paste("Creating csv for metabolomics analysis with max $max.cols columns.
                - Chosen experimental conditon is '$exp.condition'.
                - Grouping by $groupfac
                - Using adducts $adducts"))
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)
  #dbExecute(conn, "drop table if exists avg_intensities")
  # --- create table with averaged intensities for the triplicates ---
  # make.query <- strwrap(fn$paste("create table if not exists avg_intensities as
  #                                select mzmed as mz, filename, avg(intensity) as intensity 
  #                                from mzintensities 
  #                                group by filename, mzmed"), width=10000, simplify=TRUE)
  # dbExecute(conn, make.query)
  # index.query <- "create index if not exists avg_index on avg_intensities(filename, mz)"
  # --- build result fetching query ---
  if(group_adducts){
    # if(is.null(which_adducts)) which_adducts <- c("M+H", "M-H", "M-Cl", "M+K", "M")
    # if(is.null(which_adducts)) which_dbs <- c(file.path(options$db_dir, "chebi.full.db"))
    # group_by = "baseformula"
    z = get_all_matches(exp.condition, 
                        pat.conn = conn,
                        which_dbs,
                        which_adducts,
                        group_by,
                        var_table)
    z.dt <- as.data.table(z)
    cast.dt <- dcast.data.table(z.dt, 
                                card_id + sampling_date + label ~ identifier, 
                                fun.aggregate = sum, 
                                value.var = "intensity") # what to do w/ duplicates? 
  }else{
    qa <- switch(var_table,
                 individual_data = fn$paste("d.[$exp.condition] as label"),
                 setup = fn$paste("s.[$exp.condition] as label")
                 )

    query <- strwrap(fn$paste("select distinct d.card_id, 
                                           $qa, 
                              d.sampling_date, 
                              i.mzmed as identifier,
                              i.intensity
                              from mzintensities i
                              join individual_data d
                              on i.filename = d.card_id
                              join setup s on d.[Group] = s.[Group]
                              group by d.card_id, 
                              d.sampling_date, 
                              i.mzmed"),
                     width=10000,
                     simplify=TRUE)
    print(query)
    z = RSQLite::dbGetQuery(conn, query)
    z.dt <- as.data.table(z)
    cast.dt <- dcast.data.table(z.dt, 
                                card_id + sampling_date + label ~ identifier, 
                                fun.aggregate = sum, 
                                value.var = "intensity") # what to do w/ duplicates? 
  }
  # --- cast to right format ---
  max.vals.final <- if(max.vals == -1) ncol(cast.dt) else max.vals + 3
  small.set <- cast.dt[,1:(max.vals.final)]
  # --- name for metaboanalyst ---
  names(small.set)[1:3] <- c("Sample", "Time", "Label")
  # --- make time series if necessary (this factorizes sampling date) ---
  if(time.series){
      small.set$Time <- as.numeric(as.factor(as.Date(small.set$Time)))
      small.set$Sample <- paste(small.set$Sample, as.character(small.set$Time), sep="_T")
  }
  # --- measure file size ---
  size <- object.size(small.set)
  cat("Resulting file will be approximately ")
  print(size, units = "MB")
  # ------- return -------
  return(small.set)
}
