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

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}


test.matr <- PlotMBTimeProfile(cmpdNm = 125.986182503113)
test.matr

ggplot(data=test.matr) +
  geom_line(size=0.2, aes(x=time, y=abundance, group=sample, color=group)) +
  stat_summary(fun.y="mean", size=2, geom="line", aes(x=time, y=abundance, color=group, group=group))
  

# - make ggplot? -



# === GET OPTIONS ===

wd <- "/Users/jwolthuis/Google Drive/MetaboShiny"

# --- laod adduct table for general use ---

load(file.path(wd, "backend/umcfiles/adducts/AdductTableWKZ.RData"))
sourceDir(file.path(wd, "backend/scripts/joanna"))
data(isotopes)

# --- beta stuff ---
session_cl <<- NA
mode <- "time"

# --- check options ---

opt_conn <- file(".conf")
options_raw <<- readLines(opt_conn)
print(options_raw)
dbDir <<- str_match(options_raw[[1]], "(?<=')(.*)(?=')")[1,1]
exp_dir <<- str_match(options_raw[[2]], "(?<=')(.*)(?=')")[1,1]
proj_name <<- str_match(options_raw[[3]], "(?<=')(.*)(?=')")[1,1]

patdb <<- file.path(exp_dir, paste0(proj_name, ".db"))

close(opt_conn)

# === SOURCE OWN CODE ===

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

browse_db <- function(chosen.db){
  conn <- dbConnect(RSQLite::SQLite(), chosen.db) # change this to proper var later
  # --- browse ---
  result <- dbGetQuery(conn, "SELECT DISTINCT compoundname as Compound, baseformula as Formula, description as Description, charge as Charge FROM base")
  # --- result ---
  result
}

get_matches <- function(mz, chosen.db){
  # --- connect to db ---
  req("patdb")
  conn <- dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  # 0. Attach db
  query.zero <- fn$paste("ATTACH '$chosen.db' AS db")

  dbExecute(conn, query.zero)
  query.one <- fn$paste(strwrap(
    "CREATE TEMP TABLE unfiltered AS
    SELECT cpd.baseformula, cpd.adduct
    FROM mzvals mz
    JOIN mzranges rng ON rng.ID = mz.ID
    JOIN db.extended cpd indexed by e_idx2
    ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
    AND mz.foundinmode = cpd.foundinmode
    WHERE ABS(mz.[mzmed.pgrp] - $mz) < 0.000000000001",width=10000, simplify=TRUE))
  # 1. Find matches in range (reasonably fast <3)
  dbExecute(conn, query.one)
  #  2. get isotopes for these matchies (reverse search)
  query.two <- fn$paste(strwrap(
    "CREATE TEMP TABLE isotopes AS
    SELECT cpd.baseformula, cpd.adduct, cpd.isoprevalence, cpd.basecharge 
    FROM db.extended cpd indexed by e_idx1
    JOIN unfiltered u
    ON u.baseformula = cpd.baseformula
    AND u.adduct = cpd.adduct
    JOIN mzranges rng
    ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax"
    , width=10000, simplify=TRUE))
  dbExecute(conn, query.two)
  query.three <-  strwrap(
    "SELECT DISTINCT base.compoundname as Compound, base.identifier as Identifier, iso.adduct as Adduct, base.description as Description 
    FROM isotopes iso
    JOIN db.base base indexed by b_idx1
    ON base.baseformula = iso.baseformula AND
    base.charge = iso.basecharge
    GROUP BY iso.baseformula, iso.adduct
    HAVING COUNT(iso.isoprevalence > 99.99999999999) > 0"
    , width=10000, simplify=TRUE)
  # 3. get the info you want
  results <- dbGetQuery(conn,query.three)
  # --- results ---
  results
  }


get_mzs <- function(baseformula, charge, chosen.db){
  # --- connect to db ---
  req("patdb")
  conn <- dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  print(baseformula)
  print(charge)
  query.zero <- fn$paste("ATTACH '$chosen.db' AS db")
  print(query.zero)
  dbExecute(conn, query.zero)
  # search combo of baseformula and charge matching your choice and find all possible mzvals and adducts
  query.one <-  fn$paste(strwrap(
    "CREATE TEMP TABLE possible_options AS
    SELECT DISTINCT e.fullmz, e.adduct, e.isoprevalence
    FROM db.extended e
    WHERE e.baseformula = '$baseformula' 
    AND e.basecharge = $charge"
    , width=10000, simplify=TRUE))
  print(query.one)
  
  dbExecute(conn, query.one)
  
  # join with patdb
  query.two <- fn$paste(strwrap(
    "CREATE TEMP TABLE isotopes AS
    SELECT DISTINCT mz.[mzmed.pgrp], o.*
    FROM possible_options o
    JOIN mzranges rng
    ON o.fullmz BETWEEN rng.mzmin AND rng.mzmax
    JOIN mzvals mz
    ON rng.ID = mz.ID"))
  print(query.two)
  
  dbExecute(conn, query.two)
  
  # isofilter and only in
  
  query.three <-  strwrap(
    "SELECT DISTINCT iso.[mzmed.pgrp], adduct
    FROM isotopes iso
    GROUP BY iso.adduct
    HAVING COUNT(iso.isoprevalence > 99.99999999999) > 0"
    , width=10000, simplify=TRUE)
  print(query.three)
  
  results <- dbGetQuery(conn,query.three)
  print(results)
}
