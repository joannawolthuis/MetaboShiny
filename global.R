# === GENERAL OPTIONS ===

options(stringsAsFactors = FALSE)
if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10000*1024^2)

# === LOAD LIBRARIES ===

library(pacman)
library(shiny)
library(DT)
library(data.table)
library(shinyFiles)
# library(ggplot2)
# library(DBI)
# library(RSQLite)
# library(gsubfn)
# library(data.table)
# library(pbapply)
# library(enviPat)
# library(jsonlite)


# ------------------------

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

getOptions <- function(file.loc){
  opt_conn <- file(file.loc)
  # ----------------
  options_raw <<- readLines(opt_conn)
  close(opt_conn)
  # --- list-ify ---
  options <- list()
  for(line in options_raw){
    split  <- (strsplit(line, ' = '))[[1]]
    options[[split[[1]]]] = split[[2]]
  }
  # --- return ---
  options 
}

setOption <- function(file.loc, key, value){
  opt_conn <- file(file.loc)
  # -------------------------
  options <- getOptions(file.loc)
  # --- add new or change ---
  options[[key]] = value
  # --- list-ify ---
  new_options <- lapply(seq_along(options), FUN=function(i){
    line <- paste(names(options)[i], options[i], sep=" = ")
    line
  })
  print(new_options)
  writeLines(opt_conn, text = unlist(new_options))
  close(opt_conn)
}

# --- beta stuff ---

mode <- "time"

# === SOURCE OWN CODE ===

sourceAll <- function(where, 
                      which=c("general", "stats", "time", "enrich_path", "power_roc", "utils")){
  p_load(compiler)
  # ----------------------------
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

packages <<- c("data.table", "DBI", "RSQLite", "ggplot2", "minval", "enviPat",
               "plotly", "parallel", "shinyFiles", "curl", "httr", "pbapply", "sqldf", "plyr", "ChemmineR", "gsubfn", 
               "stringr", "plotly", "reshape2", "XML", "xlsx", "colourpicker", "DT","Rserve", "ellipse", 
               "scatterplot3d","pls", "caret", "lattice",
               "Cairo", "randomForest", "e1071","gplots", "som", "xtable",
               "RColorBrewer", "xcms","impute", "pcaMethods","siggenes",
               "globaltest", "GlobalAncova", "Rgraphviz","KEGGgraph",
               "preprocessCore", "genefilter", "pheatmap", "igraph",
               "RJSONIO", "SSPA", "caTools", "ROCR", "pROC", "sva")

# --------------------------
wdir <<- "/Users/jwolthuis/Google Drive/MetaboShiny"
setwd(wdir)
options <- getOptions(".conf")
dbDir <<- options$db_dir
exp_dir <<- options$work_dir
proj_name <<- options$proj_name
ppm <<- options$ppm
packages_installed <<- options$packages_installed
mz <- NULL

if(packages_installed == "Y"){
  p_load(char = packages, character.only = T)
  load(file.path(wdir, "backend/umcfiles/adducts/AdductTableWKZ.RData"))
  sourceDir(file.path(wdir, "backend/scripts/joanna"))
  data(isotopes, package = "enviPat")
}


print("loaded global settings")
