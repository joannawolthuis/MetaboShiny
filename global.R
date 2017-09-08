# === GENERAL OPTIONS ===

options(stringsAsFactors = FALSE)
if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10000*1024^2)

# === LOAD LIBRARIES ===

library(shiny)
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
library(colourpicker)

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# - make ggplot? -

# === GET OPTIONS ===

wd <- "/Users/jwolthuis/Google Drive/MetaboShiny"
setwd(wd)
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
ppm <<- as.numeric(str_match(options_raw[[4]], "(?<=')(.*)(?=')")[1,1])

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
