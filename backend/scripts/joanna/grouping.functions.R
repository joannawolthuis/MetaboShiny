#' @export
group.patdata <- function(outlist.path, matchy="matches_hmdb", concatty="names", groupy="[mzmed.pgrp]", reshape=TRUE, group=FALSE,return.table=FALSE){
  library(data.table)
  library(pbapply)
  library(RSQLite)
  library(DBI)
  library(sqldf)
  library(gsubfn)
  library(reshape2)
  # ---- connect to db ----
  conn <- dbConnect(RSQLite::SQLite(), outlist.path)
  # -------------------------
  dbExecute(conn, fn$paste("drop table if exists csv_ingredients"))
  # what to include in name column? :-)
  which.sql <- strwrap(fn$paste("CREATE TEMP TABLE csv_ingredients AS
                                select distinct [mzmed.pgrp] as mz, filename, intensity
                                from mzintensities i indexed by intindex"))

  print("Finding matches per sample...")
  result <- as.data.table(dbExecute(conn, which.sql))
  # --- reshape ---
  NULL
  # --- return ---
  return(as.data.table(result))
}
# test_table <- group.patdata(patdb)
