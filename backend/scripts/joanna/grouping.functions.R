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
  dbExecute(conn, fn$paste("drop table if exists results_with_$concatty"))
  # what to include in name column? :-)
  which.sql <- switch(concatty,
                       names=strwrap(fn$paste("create table results_with_$concatty as
                                select distinct m.[mzmed.pgrp] as mz, group_concat(distinct m.compoundname) as compoundname, i.filename as filename, i.intensity as intensity
                                from $matchy m indexed by $matchy_idx1
                                join mzintensities i indexed by intindex
                                ON i.[mzmed.pgrp] = m.[mzmed.pgrp]
			                          GROUP BY m.$groupy"),width=10000, simplify=TRUE),
                       formulas=strwrap(fn$paste("create table results_with_$concatty as
select distinct m.[mzmed.pgrp] as mz, group_concat(distinct printf('%s(%d)', m.baseformula, m.basecharge)) as baseformula, i.filename as filename, i.intensity as intensity
                                from $matchy m indexed by $matchy_idx1
                                join mzintensities i indexed by intindex
                                ON i.[mzmed.pgrp] = m.[mzmed.pgrp]
			                          GROUP BY m.$groupy"),width=10000, simplify=TRUE),
                       identifiers=strwrap(fn$paste("create table results_with_$concatty as
select distinct m.[mzmed.pgrp] as mz, group_concat(distinct m.identifier) as identifier, i.filename as filename, i.intensity as intensity
                                from $matchy m
                                join mzintensities i indexed by intindex
                                ON i.[mzmed.pgrp] = m.[mzmed.pgrp]
			                          GROUP BY m.$groupy"),width=10000, simplify=TRUE),
                       all=strwrap(fn$paste("create table results_with_$concatty as
select distinct m.[mzmed.pgrp] as mz, group_concat(distinct printf('%s(%d):%s(%s)', m.baseformula, m.basecharge, m.compoundname, m.identifier)) as allinfo, i.filename as filename, i.intensity as intensity
                                from $matchy m indexed by $matchy_idx1
                                join mzintensities i indexed by intindex
                                ON i.[mzmed.pgrp] = m.[mzmed.pgrp]
			                          GROUP BY m.$groupy"),width=10000, simplify=TRUE))
  
  print(which.sql)
  print("Finding matches per sample...")
  result <- as.data.table(dbExecute(conn, which.sql))
  # --- reshape ---
  NULL
  # --- return ---
  return(as.data.table(result))
}
# test_table <- group.patdata(patdb)

#reformat.patdata(patdb)
#name.patdata(patdb, file.path(dbDir, "hmdb.base.db"))

#reshape with sql...

# -- join mzval and name columns ---

# --- small test set ---

"SELECT * FROM results_with_names WHERE filename = 'BSP20160322-99' OR filename = 'BSP20160322-98'"

# --- reshape into long ---

