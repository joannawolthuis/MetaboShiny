# === GENERAL OPTIONS ===

options(stringsAsFactors = FALSE)
if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10000*1024^2)

# === LOAD LIBRARIES ===

library(ggplot2)
library(DT)
library(DBI)
library(RSQLite)
library(gsubfn)
library(data.table)
library(parallel)
library(pbapply)
library(enviPat)
library(plotly)
library(jsonlite)
library(shinyFiles)
library(stringr)

# === GET OPTIONS ===

setwd("/Users/jwolthuis/Google Drive/MetaboShiny")
# --- laod adduct table for general use ---
load("./backend/db/NeededFiles/AdductTable/AdductTableWKZ.RData")
sourceDir("./backend/scripts/joanna")
data(isotopes)
dbDir <<- file.path("./backend/db")
# --- beta stuff ---
session_cl <<- NA
mode <- "time"
opt_conn <- file(".conf")
options_raw <<- readLines(opt_conn)
exp_dir <<- str_match(options_raw[[1]], "(?<=')(.*)(?=')")[1,1]
proj_name <<- str_match(options_raw[[2]], "(?<=')(.*)(?=')")[1,1]
close(opt_conn)

# === SOURCE OWN CODE ===

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceAll <- function(where, 
                      which=c("general", "stats", "time", "enrich_path", "power_roc", "utils")){
  library(compiler)
  print("sourcing R code ... ");
  for(i in 1:length(which)){
    script.loc <- file.path(where, which[i])
    print(script.loc)
    files <- list.files(file.path(where, which[i]),full.names=TRUE, pattern=".R$");
    print(files)
    for(f in files){
      print(f)
      source(f)
    }
  }
  return("TRUE");
}

############ BUILD DATABASES ###############

# === SOURCE METABOANALST CODE ===

# === LATER FUNCTIONS ===

get_exp_vars <- function(){
  conn <- dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  dbGetQuery(conn, "PRAGMA table_info(setup)")$name
}

get_matches <- function(mz, table){
  # --- connect to db ---
  req("patdb")
  conn <- dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  # --- attach patient outlist and get mzmed pgrp values ---
  dbGetQuery(conn, fn$paste("select distinct compoundname as Compound, identifier as Identifier, Adduct as Adduct from $table where [mzmed.pgrp] like $mz"))
}