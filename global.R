# === LOAD LIBRARIES ===

library(ggplot2)
library(DT)
library(DBI)
library(RSQLite)
library(gsubfn)
library(data.table)



# === SOURCE OWN CODE ===

# === SOURCE METABOANALST CODE ===

# === LATER FUNCTIONS ===

get_matches <- function(mz, table){
  # --- connect to db ---
  conn <- dbConnect(RSQLite::SQLite(), file.path("backend","data","time.trials.db"))
  # --- attach patient outlist and get mzmed pgrp values ---
  dbGetQuery(conn, fn$paste("select distinct compoundname as Compound, identifier as Identifier, Adduct as Adduct from $table where [mzmed.pgrp] like $mz"))
}